#!/bin/bash

# We enter in the ORFtrack directory
cd /ORFmine/orftrack/

# # We uninstall the previous ORFtrack
# /bin/miniconda3/envs/ORFmine_env/bin/pip uninstall orftrack

# We install ORFtrack
/bin/miniconda3/envs/$1/bin/pip install -e .


# We leave ORFtrack directory and enter in the ORFold directory
cd /ORFmine/orfold_v1/

# We prepare the folder for IUPred and Tango
if [ ! -d "orfold/softwares" ];
then
    mkdir -p orfold/softwares
fi

# # We uninstall the previous orfold
# /bin/miniconda3/envs/ORFmine_env/bin/pip uninstall orfold

# Then we extract the iupred functions from the iupred code
if [ -f ./orfold/softwares/iupred2a.py ]
then
    /bin/miniconda3/envs/$1/bin/python3 /ORFmine/orfold_v1/extract_iupred_functions.py
fi

# Install pyHCA for ORFold
/bin/miniconda3/envs/$1/bin/pip install -e /ORFmine/orfold_v1/orfold/pyHCA/

# And reinstall orfold
/bin/miniconda3/envs/$1/bin/pip install -e .

# We enter in the ORFdate directory
cd /ORFmine/orfdate/

# # We uninstall the previous ORFdate
# /bin/miniconda3/envs/ORFmine_env/bin/pip uninstall orfdate

# We install ORFdate
/bin/miniconda3/envs/$1/bin/pip install -e .

# Back to root directory
cd /
