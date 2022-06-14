#!/usr/bin/env bash

# H2, 1angstrom, STO-3G
abspath="$(cd $(dirname $0); pwd)" # Get the path of this script
cd "$abspath" || return 1 # cd to the path of this script
binarypath="../../bin" # Binary path is specified by relative path

# Object files path
r4dcasci="$binarypath/r4dcascicoexe"
r4dcaspt2="$binarypath/r4dcaspt2ocoexe"

# Object files check
echo "Object file check"
ls $binarypath
echo "$r4dcasci, $r4dcaspt2"

# Output path
outpath="$abspath/H2.caspt2.out"

# Ouput file path check

echo "Ouput file path check : $outpath"

echo "h2.sh abspath $abspath $binarypath $r4dcasci $r4dcaspt2"
# OMP_NUM_THREADS=4
# Execute CASCI/CASPT2 programs
$r4dcasci &> $outpath
$r4dcaspt2 &>> $outpath

ls $abspath
cat $outpath
# Remove unnecessary files
rm [A-H]*int* MDCINTNEW NEWCICOEFF CIMAT* e0after EPS TRANSFOCK
