# -*- coding: utf-8 -*-
"""
Created on Fri Jul 24 15:08:04 2020

@author: nicolas
"""
import setuptools as st
#import os

#os.system('python3 ./extract_iupred_functions.py')

packages = ['orfold', 'orfold.lib', 'orfold.scripts','orfold.softwares']

st.setup(name='orfold',
         python_requires='>=3.6',
         packages=packages,
         install_requires=['biopython>=1.78',
         				   'scipy>=1.5.4',
         				   'seaborn>=0.11.0',
         				   'matplotlib>=3.3.3'],
         entry_points={'console_scripts': ['orfold=orfold.orfold:main',
                                           #'orfget=orfold.scripts.ORFget:main',
                                           'orfplot=orfold.scripts.plot_orfold:main']}
         )
