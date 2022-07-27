import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

def dotplot(
    results_df,
    y="y",
    x="clumpsptm_fdr",
    thresh=0.1,
    sort_by="ptm",
    ax=None,
    shape_dict = {'ptm':'o', 'phosphoproteome':',', 'acetylome':'^'},
    name_dict = {'ptm':'Combined', 'phosphoproteome':'Phosphorylation', 'acetylome':'Acetylation'},
    color_dict = {'fdr':'red', 'pval':'orange', 'nsig':'lightgrey'},
    title = None,
    show_gene_only = True,
    ):
    """
    Plot summaryplot for CLUMPS-PTM.
    """
    import matplotlib.patches as mpatches

    def _assign_color(row):
        if row['clumpsptm_fdr'] <= thresh:
            return color_dict['fdr']
        elif row['clumpsptm_fdr'] <= 0.25:
            return color_dict['pval']
        else:
            return color_dict['nsig']

    # Select Structures to Plot (y-axis)
    y_to_plot = np.unique(results_df[
        (results_df['clumpsptm_pval']<=thresh) &
        (results_df['clumpsptm_sampler']==sort_by)
        ][y])

    if y_to_plot.shape[0] == 0:
        print("No Nominal P-Values < 0.1")
        return

    # Subset for Data
    data_to_plot = results_df[results_df[y].isin(y_to_plot)].copy()
    assays = set(data_to_plot['clumpsptm_sampler']) - {sort_by}
    data_to_plot['category'] = pd.Categorical(data_to_plot['clumpsptm_sampler'], [sort_by]+list(assays))

    # Log FDR & color
    data_to_plot.loc[data_to_plot[data_to_plot[x]==0].index,x] = np.min(
        data_to_plot[x][data_to_plot[x]>0]*0.1)

    data_to_plot['x'] = -np.log10(data_to_plot[x])
    data_to_plot['c'] = data_to_plot.apply(_assign_color, 1)

    # Sort values
    data_to_plot = data_to_plot.sort_values(["category","x"], ascending=[True, True])

    # Create figure
    if ax is None:
        fig, ax = plt.subplots(figsize=(4, 2 + int(data_to_plot.shape[0]/20)))

    # Plot
    for idx, cat in enumerate(data_to_plot['category'].cat.categories):
        ax.scatter(
            y=y, x="x", c="c",
            data=data_to_plot[data_to_plot['clumpsptm_sampler']==cat],
            linewidth=1, edgecolor="black", zorder=10-idx,
            marker=shape_dict[cat], label=name_dict[cat], s=100, alpha=0.6
        )

    for j in range(y_to_plot.shape[0]):
        ax.axhline(j, c='lightgrey', alpha=0.4)

    ax.set_ylim(-1, y_to_plot.shape[0]+.5)
    ax.axvline(-np.log10(thresh), linewidth=0.5, c='lightgrey')
    ax.tick_params(axis='y', which='major', labelsize=7)

    if 'pval' in x:
        ax.set_xlabel("-log_10 P-Val")
    elif 'fdr' in x:
        ax.set_xlabel("-log10 FDR")

    handles, labels = ax.get_legend_handles_labels()

    # Patches
    patch = mpatches.Patch(
        facecolor=color_dict['fdr'],
        label="< {} {}".format(0.1, "FDR"),
        edgecolor='k',
        linewidth=1,
    )

    handles.append(patch)

    patch = mpatches.Patch(
        facecolor=color_dict['pval'],
        label="< {} {}".format(0.25, "FDR"),
        edgecolor='k',
        linewidth=1,
    )

    handles.append(patch)

    ax.legend(handles=handles, loc='lower left', bbox_to_anchor=(1, 0))

    if title is None:
        ax.set_title(results_df['id'][0])

    if show_gene_only and y=='y':
        plt.draw()
        labels = [item.get_text().split(" | ")[1] for item in ax.get_yticklabels()]
        ax.set_yticklabels(labels)
