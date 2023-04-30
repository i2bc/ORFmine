# -*- coding: utf-8 -*-
"""
Created on Fri Jul 24 15:08:04 2020

@author: nicolas
"""

import setuptools as st

MIN_PY_VER = "3.7"

REQUIRES = [
    'biopython==1.78',
    'cutadapt==4.1',
    'DendroPy==4.5.2',
    'gffutils==0.11.1',
    'scipy==1.9.3',
    'seaborn==0.11.0',
    'matplotlib==3.3.3',
    'pandas',
    'lightgbm',
    'joblib',
    'CairoSVG',
    'pysam==0.19.1',
    'requests==2.28.1',
    "urllib3==1.26.13",
    'six==1.12',
    'snakemake==7.16.0',
    'ete3==3.1.1',
    'pyHCA @ git+https://github.com/T-B-F/pyHCA.git',
]

PACKAGES = [
    'orfmine',
    'orfmine.orftrack',
    'orfmine.orftrack.lib',
    'orfmine.orftrack.scripts',
    'orfmine.orfold',
    'orfmine.orfold.lib',
    'orfmine.orfold.scripts',
    'orfmine.orfdate',
    'orfmine.orfribo',
    'orfmine.orfribo.lib',
    'orfmine.orfribo.scripts',
    'orfmine.utilities',
    'orfmine.utilities.lib',
]

st.setup(
    name='orfmine',
    python_requires=">={}".format(MIN_PY_VER),
    packages=PACKAGES,
    include_package_data=True,
    install_requires=REQUIRES,
    entry_points={
        'console_scripts': [
            'orftrack=orfmine.orftrack.main:main',
            'orfget=orfmine.orftrack.scripts.ORFget:main',
            'gff2prot=orfmine.orftrack.scripts.gff2prot:main',
            'orfold=orfmine.orfold.orfold:main',
            'orfplot=orfmine.orfold.scripts.plot_orfold:main',
            'orfdate=orfmine.orfdate.ORFdate:main',
            'orfribo=orfmine.orfribo.orfribo:main',
            'bam2reads=orfmine.orfribo.scripts.BAM2Reads:main',
            'orfstats=orfmine.orfribo.scripts.ORFstats:main',
            'merge_read_tables=orfmine.orfribo.scripts.concatenate:main',            
        ]
    }
)