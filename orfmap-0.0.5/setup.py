# -*- coding: utf-8 -*-
"""
Created on Fri Jul 24 15:08:04 2020

@author: nicolas
"""

import setuptools as st

packages = ['orfmap', 'orfmap.lib', 'orfmap.scripts']

st.setup(name='orfmap',
         python_requires='>=3',
         version='0.0',
         packages=packages,
         entry_points={'console_scripts': ['orfmap=orfmap.main:main',
                                           'gff2prot=orfmap.scripts.gff2prot:main']}
         )
