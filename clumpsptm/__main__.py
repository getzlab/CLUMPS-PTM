import pandas as pd
import argparse
import numpy as np
import os
from .clumps import transform_dx, init_alg, clumps
from .pdbstore import PdbStore
from .samplers import PTMSampler
from .utils import generate_clumpsptm_output
from agutil.parallel import parallelize2

def main():
    parser = argparse.ArgumentParser(description='Run CLUMPS-PTM.')
    parser.add_argument('-i', '--input', required=True, help='<Required> Input file.')
    parser.add_argument('-m','--maps', required=True, help='<Required> Mapping with index as indices that overlap input.')
    parser.add_argument('-w', '--weight', required=True, help='<Required> Weighting for CLUMPS-PTM (ex. logFC).')
    parser.add_argument('-s','--pdbstore',required=True, help='<Required> path to PDBStore directory.')
    parser.add_argument('-o','--output_dir', default=".", help='Output directory.')
    parser.add_argument('-x', '--xpo', default=[3, 4.5, 6, 8, 10], type=list,
        help='Soft threshold parameter for truncated Gaussian.')
    parser.add_argument('--threads', type=int, default=1, help='Number of threads for sampling.')
    parser.add_argument('-v', '--verbose', action='store_true', default=False, help='Verbosity.')
    parser.add_argument('-f', '--features', nargs="*", default=None, help='Assays to subset for.')
    parser.add_argument('-g', '--grouping', default=None, help='DE group to use.')
    parser.add_argument('-q', '--use_only_significant_sites', action='store_true', help='Only use significant sites for CLUMPS-PTM.')
    parser.add_argument('--min_sites', default=3, help='Minimum number of sites.')
    parser.add_argument('--subset', default=None, help='Subset sites.', choices=('positive','negative'))
    args = parser.parse_args()

    print("---------------------------------------------------------")
    print("------------------ C L U M P S - P T M ------------------")
    print("---------------------------------------------------------")

    # ----------------------------------
    # Load files
    # ----------------------------------
    mapped_sites_df = pd.read_csv(args.maps, sep='\t', index_col=0)
    os.makedirs(args.output_dir, exist_ok=True)
    os.makedirs(os.path.join(args.output_dir, "clumpsptm"), exist_ok=True)

    if args.verbose:
        print("   * creating output directory: {}".format(args.output_dir))

    if args.verbose:
        print("   * loading {} mapped sites".format(mapped_sites_df.shape[0]))

    # Load PDB Store
    pdbstore = PdbStore(args.pdbstore)

    # Load input file
    if args.input.endswith(".parquet"):
        de_df = pd.read_parquet(args.input)
        try:
            de_df = de_df.set_index("index")
        except:
            pass
    else:
        de_df = pd.read_csv(args.input, sep='\t')

        if 'index' in list(de_df):
            de_df = de_df.set_index("index")
        else:
            de_df = de_df.set_index(list(de_df)[0])

    de_df['id'] = de_df['id'].astype(str)

    if args.grouping is not None:
        de_df = de_df[de_df['id']==str(args.grouping)]
    if args.features is not None:
        de_df = de_df[de_df['feature'].isin(args.features)]

    # Combined Input with Mapping File
    de_df = mapped_sites_df.loc[np.intersect1d(mapped_sites_df.index, de_df.index)].join(
        de_df.drop(columns=np.intersect1d(de_df.columns, mapped_sites_df.columns))
    )

    # Subset for pos/neg sites
    if args.subset is not None:
        if args.verbose:
            print("   * subsetting for {} sites".format(args.subset))

        if args.subset == "positive":
            de_df = de_df[de_df[args.weight]>0]
        elif args.subset == "negative":
            de_df = de_df[de_df[args.weight]<0]
            de_df[args.weight] = -1 * de_df[args.weight]

    # Subset for significant sites
    # TODO: allow to select between p-value, fdr, and select threshold
    if args.use_only_significant_sites:
        print("   * filtering for significant sites below < 0.1 nominal-pvalue")
        de_df = de_df[de_df['P.Value']<0.1]

    # Subset for min sites
    gb = de_df.groupby("accession_number").size()
    de_df = de_df[de_df["accession_number"].isin(gb[gb>=args.min_sites].index)]
    de_df['pdb_res_i'] = de_df['pdb_res_i'].astype(int)

    if args.verbose:
        print("   * {} input sites mapped".format(de_df.shape[0]))
        print("   * {} proteins with > {} sites".format(np.unique(de_df["accession_number"]).shape[0], args.min_sites))

    # ----------------------------------
    # Run Clumps
    # ----------------------------------
    @parallelize2(maximum=args.threads)
    def run_clumps(acc):
        """
        Run CLUMPS PTM.
        """
        inputs_df = de_df[de_df['accession_number']==acc].copy()
        inputs_df = inputs_df.sort_values('pdb_res_i')
        pdb, chain = inputs_df.iloc[0][['pdb','chain']]

        # Rule for multi-site modifications
        inputs_df = inputs_df.sort_values(['pdb_res_i', args.weight], ascending=[True, False]).drop_duplicates('pdb_res_i')

        # Get PDB Information
        D,x,pdb_resnames = pdbstore.load_dm(pdb, chain)
        DDt = transform_dx(D, args.xpo)

        # Get matched numerical site
        pdb_struct_i_to_d_map = {x:i for i,x in enumerate(np.sort(list(x.keys())))}
        inputs_df['d'] = inputs_df['pdb_res_i'].apply(lambda x: pdb_struct_i_to_d_map[int(x)])
        inputs_df['value'] = inputs_df[args.weight]

        def _run(feature):
            """Run."""
            if feature == "acetylome":
                _inputs_df = inputs_df[inputs_df["feature"]==feature].copy()
                _res = ["LYS"]
            elif feature == "phosphoproteome":
                _inputs_df = inputs_df[inputs_df["feature"]==feature].copy()
                _res = ["SER","THR","TYR"]
            elif feature == "ptm":
                _inputs_df = inputs_df.copy()
                _res = ["LYS","SER","THR","TYR"]

            if _inputs_df.shape[0] ==0:
                # print("     > {} - {} | p-value: {}".format(acc, feature, "no sites."))
                return

            # Initialize
            init_dict = init_alg(_inputs_df, DDt)
            sam = PTMSampler(pdb_struct_i_to_d_map, pdb_resnames, res=_res)

            # Run CLUMPS
            pval, wap_scr, niter = clumps(init_dict, sam, DDt, args.xpo, max_rand=1e4, use_booster=False)

            # Output
            # print("     > {} - {} | p-value: {}".format(acc, feature, sum(pval) / (5*niter)))
            _inputs_df.loc[:,"clumpsptm_niter"] = niter
            _inputs_df.loc[:,"clumpsptm_pval"] = sum(pval) / (5*niter)

            res_input = [str(x) for x in _inputs_df['pdb_res_i']]
            res_sample = [str(x) for x in sam.ptm]

            _inputs_df.loc[:,"clumpsptm_input"] = "+".join(res_input)
            _inputs_df.loc[:,"clumpsptm_sample"] = "+".join(res_sample)
            _inputs_df.loc[:,"clumpsptm_input_n"] = len(res_input)
            _inputs_df.loc[:,"clumpsptm_sample_n"] = len(res_sample)
            _inputs_df.to_parquet(os.path.join(args.output_dir, "clumpsptm", "{}_{}_clumpsptm.parquet".format(acc, feature)))

        if "acetylome" in args.features:
            try:
                _run("acetylome")
            except:
                print("     > ERROR for {} - {}".format(acc, "acetylome"))
        if "phosphoproteome" in args.features:
            try:
                _run("phosphoproteome")
            except:
               print("     > ERROR for {} - {}".format(acc, "phosphoproteome"))
        if "acetylome" in args.features and "phosphoproteome" in args.features:
            try:
                _run("ptm")
            except:
               print("     > ERROR for {} - {}".format(acc, "ptm"))

    print("   * running using {} threads".format(args.threads))

    accession_nos = np.unique(de_df['accession_number'])

    tmp = [run_clumps(acc) for acc in accession_nos]
    results = [callback() for callback in tmp]

    # ----------------------------------
    # Generate Output File
    # ----------------------------------
    results_df = generate_clumpsptm_output(os.path.join(args.output_dir, "clumpsptm"))
    results_df.to_csv(os.path.join(args.output_dir, "clumpsptm_output.tsv"), sep="\t")

if __name__ == "__main__":
    main()
