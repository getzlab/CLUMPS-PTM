from setuptools import setup

setup(
    name='clumpsptm',
    author='Shankara Anand',
    author_email='sanand@broadinstitute.org',
    description='CLUMPS-PTM driver gene discovery using 3D protein structure (Getz Lab).',
    install_requires=[
        "twobitreader>=3.1.7",
        "statsmodels>=0.9.0",
        "scipy>=1.3.0",
        "pyopenssl>=19.0.0",
        "prody>=1.10.10",
        "lxml>=4.3.4",
        "canine>=0.1.3",
        "biopython>=1.73",
        "tqdm>=4.32.1"
    ],
)
