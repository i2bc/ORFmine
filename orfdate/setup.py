# -*- coding: utf-8 -*-
import setuptools as st

packages = ['orfdate']

st.setup(name='orfdate',
    python_requires='>=3.6',
    packages=packages,
    entry_points={
        'console_scripts': [
            'orfdate=orfdate.ORFdate:run'
        ]
    }
)
