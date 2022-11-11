# -*- coding: utf-8 -*-
"""
Created on Fri Jul 24 15:08:04 2020

@author: nicolas
"""

import setuptools as st

MIN_PY_VER = "3.7"

REQUIRES = [
    'biopython==1.78',
    'DendroPy==4.5.2',
    'scipy==1.9.3',
    'seaborn==0.11.0',
    'matplotlib==3.3.3',
    'scikit-learn>=0.19',
    'pandas',
    'lightgbm',
    'joblib',
    'CairoSVG',
    'requests>=2.18',
    'six>=1.11',
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
]

st.setup(
    name='orfmine',
    python_requires=">={}".format(MIN_PY_VER),
    packages=PACKAGES,
    install_requires=REQUIRES,
    entry_points={
        'console_scripts': [
            'orftrack=packages.orftrack.main:main',
            'orfget=packages.orftrack.scripts.ORFget:main',
            'orfold=packages.orfold.orfold:main',
            'orfplot=packages.orfold.scripts.plot_orfold:main',
            'orfdate=packages.orfdate.ORFdate:main'
        ]
    }
)