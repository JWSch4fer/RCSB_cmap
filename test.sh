#!/bin/bash -e

#expecting this to be run from the repositories home directory
CURRENTPATH=`pwd`
source "${CURRENTPATH}/conda/bin/activate"

echo "###________________________________________________________________________###"
echo " "
echo "Running tests in examples/"
echo " "
echo "###________________________________________________________________________###"

#move into the examples directory
cd examples/

    #this one should fail 
    python ../RCSB_cmap.py > ${CURRENTPATH}/examples/test_output.txt;
    #download a monomer
    python ../RCSB_cmap.py -pdb 1t4y>> ${CURRENTPATH}/examples/test_output.txt;
    #local monomer
    python ../RCSB_cmap.py -pdb 1t4y.pdb>> ${CURRENTPATH}/examples/test_output.txt;
    #download an oligomer
    python ../RCSB_cmap.py -pdb 3egm>> ${CURRENTPATH}/examples/test_output.txt;
    #compare two conformations side by side (dual-fold contact map)
    python ../RCSB_cmap.py -pdb 6c6s -chain D -pdb2 5ond -chain2 A>> ${CURRENTPATH}/examples/test_output.txt;

cd ../

echo "###________________________________________________________________________###"
echo " "
echo "Tests complete:"
echo "6 png files should be in examples/"
echo "output was captured in examples/test_output.txt"
echo " "
echo "###________________________________________________________________________###"

