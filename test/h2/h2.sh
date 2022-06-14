#!/usr/bin/env bash

# H2, 1angstrom, STO-3G
abspath="$(dirname $0)"
binarypath="$abspath/../../bin"
r4dcasci="$binarypath/r4dcascicoexe"
r4dcaspt2="$binarypath/r4dcaspt2ocoexe"
echo "h2.sh abspath $abspath $binarypath $r4dcasci $r4dcaspt2"
# OMP_NUM_THREADS=4
# Execute CASCI/CASPT2 programs
$r4dcasci &> H2.caspt2.out
$r4dcaspt2 &>> H2.caspt2.out

# Remove unnecessary files
rm [A-H]*int* MDCINTNEW NEWCICOEFF CIMAT* e0after EPS TRANSFOCK
