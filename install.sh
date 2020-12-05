#!/bin/bash

# We enter in the ORFmap directory
cd orfmap-0.0.5

# We uninstall the previous ORFmap
pip uninstall orfmap

# We install ORfmap
python setup.py install

# We leave ORFmap directory
cd ..

#We enter in the ORFold directory
cd orfold_v1

# First we change the path of the softwares in the orfold script
courent_path=$(pwd)
sed -i '' "s|.*softwares_path=.*|softwares_path=\"${courent_path}\"|g" ./orfold/orfold.py

# Then we extract the iupred functions from the iupred code
if [ -f ./orfold/softwares/iupred2a.py ]
then
python extract_iupred_functions.py
fi

# We uninstall the previous orfold
pip uninstall orfold

# And reinstall the new one
pip install .


