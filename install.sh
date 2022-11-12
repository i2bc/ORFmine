#!/bin/bash

# ORFTRACK
# --------
# install ORFtrack
cd /ORFmine/orftrack/
/bin/miniconda3/envs/$1/bin/pip install -e .


# ORFOLD
# --------
cd /ORFmine/orfold_v1/

# ensure softwares mounting point exists (for iupred & tango)
if [ ! -d "orfold/softwares" ];
then
    mkdir -p orfold/softwares
fi

# Install pyHCA
/bin/miniconda3/envs/$1/bin/pip install -e /ORFmine/orfold_v1/orfold/pyHCA/

# install orfold
/bin/miniconda3/envs/$1/bin/pip install -e .


# Then we extract the iupred functions from the iupred code
# if [ -f ./orfold/softwares/iupred2a.py ]
# then
#     /bin/miniconda3/envs/$1/bin/python3 /ORFmine/orfold_v1/extract_iupred_functions.py
# fi


# ORFDATE
# --------
# install ORFdate
cd /ORFmine/orfdate/
/bin/miniconda3/envs/$1/bin/pip install -e .


# Back to root directory
cd /
