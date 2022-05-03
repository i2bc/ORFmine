#!/bin/bash

# We enter in the ORFtrack directory
cd orftrack

# We uninstall the previous ORFtrack
pip uninstall orftrack

# We install ORFtrack
python setup.py install

# We leave ORFtrack directory
cd ..

#We enter in the ORFold directory
cd orfold_v1

# create softwares directory for orfold 
mkdir orfold/softwares 

# First we change the path of the softwares in the orfold script
courent_path=$(pwd)
if [ $(uname -s) == "Darwin" ]
  then
  sed -i '' "s|.*softwares_path=.*|softwares_path=\"${courent_path}\"|g" ./orfold/orfold.py
elif [ $(uname -s) == "Linux" ]
  then
  sed -i "s|.*softwares_path=.*|softwares_path=\"${courent_path}\"|g" ./orfold/orfold.py
 fi
 
 
# Second we change the path of the data in the plot_orfold script
sed -i '' "s|.*glo_ref_path=.*|glo_ref_path=\"${courent_path}\"|g" ./orfold/scripts/plot_orfold.py

# Then we extract the iupred functions from the iupred code
if [ -f ./orfold/softwares/iupred2a.py ]
then
python extract_iupred_functions.py
fi

# We uninstall the previous orfold
pip uninstall orfold

# And reinstall the new one
pip install .


