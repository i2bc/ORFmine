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
    'requests==2.18',
    'six==1.12',
    'snakemake==7.16.0',
    'ete3==3.1.1',
    'pyHCA @ git+https://github.com/T-B-F/pyHCA.git',
]

PACKAGES = [
    'packages',
    'packages.orftrack',
    'packages.orftrack.lib',
    'packages.orftrack.scripts',
    'packages.orfold',
    'packages.orfold.lib',
    'packages.orfold.scripts',
    'packages.orfdate',
    'packages.orfribo',
    'packages.orfribo.lib',
    'packages.orfribo.scripts',
]

st.setup(
    name='orfmine',
    python_requires=">={}".format(MIN_PY_VER),
    packages=PACKAGES,
    include_package_data=True,
    install_requires=REQUIRES,
    entry_points={
        'console_scripts': [
            'orftrack=packages.orftrack.main:main',
            'orfget=packages.orftrack.scripts.ORFget:main',
            'gff2prot=packages.orftrack.scripts.gff2prot:main',
            'orfold=packages.orfold.orfold:main',
            'orfplot=packages.orfold.scripts.plot_orfold:main',
            'orfdate=packages.orfdate.ORFdate:main',
            'bam2reads=packages.orfribo.scripts.BAM2Reads:main',
            'orfstats=packages.orfribo.scripts.ORFstats:main',
            'merge_read_tables=packages.orfribo.scripts.concatenate:main',            
        ]
    }
)