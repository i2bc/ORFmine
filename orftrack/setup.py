# -*- coding: utf-8 -*-
"""
Created on Fri Jul 24 15:08:04 2020

@author: nicolas
"""

import setuptools as st

packages = ['orftrack', 'orftrack.lib', 'orftrack.scripts']

st.setup(
    name='orftrack',
    python_requires='>=3',
    packages=packages,
    install_requires=['biopython>=1.78'],
    entry_points={
        'console_scripts': [
            'orftrack=orftrack.main:main',
            'orfget=orftrack.scripts.ORFget:main'
        ]
    }
)
