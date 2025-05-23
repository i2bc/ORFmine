# -*- coding: utf-8 -*-
"""
Created on Fri Jul 24 15:08:04 2020

@author: nicolas
"""

import setuptools as st

MIN_PY_VER = "3.9"

REQUIRES = [
      'biopython==1.79',
      'cairocffi==1.7.1',
      'cairosvg==2.7.1',
      'contourpy==1.3.0',
      'cssselect2==0.7.0',
      'cycler==0.12.1',
      'defusedxml==0.7.1',
      'ete3==3.1.3',
      'fonttools==4.54.1',
      'freetype-py==2.5.1',
      'kiwisolver==1.4.7',
      'lightgbm==4.5.0',
      'matplotlib==3.9.2',
      'pillow==10.4.0',
      'scipy==1.13.1',
      'seaborn==0.13.2',
      'tinycss2==1.3.0',
      'webencodings==0.5.1',
      'snakemake==7.16.0',
      'multiqc==1.25.1',
      'DendroPy==4.5.2',
      'PyQt5==5.15.11',
      'pulp==2.7.0',
      #'pyHCA @ git+https://github.com/T-B-F/pyHCA.git',
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
    'orfmine.orfdate.lib',
    'orfmine.orfribo',
    'orfmine.orfribo.lib',
    'orfmine.orfribo.scripts',
    'orfmine.utilities',
    'orfmine.utilities.lib',
]

st.setup(
    name='orfmine',
    version='3.0.17',
    python_requires=">={}".format(MIN_PY_VER),
    packages=PACKAGES,
    package_data={
        "orfmine.orfribo": ["config.yaml", "Snakefile", "Rscripts/*"],
        "orfmine.orfold": ["data/*.tab"],
        "orfmine.utilities": ["data/*"],
    },
    include_package_data=True,
    install_requires=REQUIRES,
    entry_points={
        'console_scripts': [
            'orftrack=orfmine.orftrack.orftrack:main',
            'orfget=orfmine.orftrack.scripts.ORFget:main',
            'gff2prot=orfmine.orftrack.scripts.gff2prot:main',
            'orfold=orfmine.orfold.orfold:main',
            'orfplot=orfmine.orfold.scripts.orfplot:main',
            'orfdate=orfmine.orfdate.orfdate:main',
            'orfribo=orfmine.orfribo.orfribo:main',
            'bam2reads=orfmine.orfribo.scripts.BAM2Reads:main',
            'orfstats=orfmine.orfribo.scripts.ORFstats:main',
            'selected_length=orfmine.orfribo.scripts.selected_length:main',
            'report=orfmine.orfribo.scripts.report:main',
            'concatenate=orfmine.orfribo.scripts.concatenate:main', 
            'orfmine=orfmine.scripts:main',
        ]
    }
)
