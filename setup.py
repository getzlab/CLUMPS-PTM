from setuptools import setup
import sys

if sys.version_info < (3,6):
    import warnings
    print("This python version is not supported:")
    print(sys.version)
    print("CLUMPS-PTM requires python 3.6 or greater")
    sys.exit(1)

from clumpsptm import __version__ as version

setup(
    name='clumps-ptm',
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
        "biopython>=1.73",
        "tqdm>=4.32.1",
        "agutil"
    ],
    entry_points={
        'console_scripts':[
            'clumpsptm = clumpsptm.__main__:main',
        ]
    }
)
