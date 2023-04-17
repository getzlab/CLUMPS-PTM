from setuptools import setup
from setuptools import find_namespace_packages

# Load the README file.
with open(file="README.md", mode="r") as readme_handle:
    long_description = readme_handle.read()

setup(
    name='clumps-ptm',
    author='Shankara Anand',
    version="0.0.6",
    author_email='sanand@broadinstitute.org',
    description='CLUMPS-PTM driver gene discovery using 3D protein structure (Getz Lab).',
    long_description=long_description,
    long_description_content_type="text/markdown",
    url='https://github.com/getzlab/CLUMPS-PTM',
    keywords='cancer, bioinformatics, genomics, proteomics, proteins, alphafold, post-translational modifications, phosphorylation, acetylation',
    python_requires='>=3.6',
    packages=[
        'clumpsptm',
        'clumpsptm.clumps',
        'clumpsptm.mapping',
        'clumpsptm.samplers',
        'clumpsptm.vis'
    ],
    install_requires=[
        "matplotlib>=3.2.1",
        "twobitreader>=3.1.7",
        "statsmodels>=0.9.0",
        "scipy>=1.3.0",
        "pyopenssl>=19.0.0",
        "prody>=1.10.10",
        "lxml>=4.3.4",
        "biopython>=1.73",
        "tqdm>=4.32.1",
        "agutil"
    ],
    entry_points={
        'console_scripts':[
            'clumpsptm = clumpsptm.__main__:main',
        ]
    },
    classifiers = [
        "Development Status :: 3 - Alpha",
        "Programming Language :: Python :: 3",
        "Intended Audience :: Science/Research",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
        "Typing :: Typed",
        "License :: OSI Approved :: MIT License",
    ],
)
