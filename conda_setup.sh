#!/bin/bash -e

#make sure you have curl installed
type curl 2>/dev/null || { echo "curl is not installed!!!" ; exit 1 ; }

CURRENTPATH=`pwd`

#this repository holds the minimal installers for Conda and Mamba specific to conda-forge, 
curl -L -O "https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-$(uname)-$(uname -m).sh"
bash Miniforge3-$(uname)-$(uname -m).sh -b -p "${CURRENTPATH}/conda"
rm -f Miniforge3-$(uname)-$(uname -m).sh
source "${CURRENTPATH}/conda/etc/profile.d/conda.sh"
conda activate base
conda update -n base conda -y

#update base environment to run RCSB_cmap.py
conda install -c conda-forge biopython -y
conda install -c conda-forge matplotlib -y
conda install pandas -y
conda install scipy -y

echo "###________________________________________________________________________###"
echo "The python environment for RCSB_cmap.py has been prepared"
echo " "
echo "Activate this environment with:  source ${CURRENTPATH}/conda/bin/activate"
echo " "
echo "Run the script python RCSB_cmap.py -pdb ####"
echo "###________________________________________________________________________###"

