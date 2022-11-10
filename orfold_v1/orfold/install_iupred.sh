#!/bin/bash


# We leave ORFtrack directory and enter in the ORFold directory
cd /ORFmine/orfold_v1/

# # create softwares directory for orfold
# mkdir ./orfold/softwares/

# Then we extract the iupred functions from the iupred code
if [ -f ./orfold/softwares/iupred2a.py ]
then
    /bin/miniconda3/envs/ORFmine_env/bin/python3 /ORFmine/orfold_v1/extract_iupred_functions.py
fi

# We uninstall the previous orfold
/bin/miniconda3/envs/$1/bin/pip uninstall orfold

# And reinstall the new one
/bin/miniconda3/envs/$1/bin/pip install -e .

# Back to root directory
cd /
