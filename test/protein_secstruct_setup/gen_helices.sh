#!/bin/bash

export AMBERHOME=/opt/amber

cat > leap.in << EOF
source leaprc.ff14SB

mola = sequence { ALA ALA ALA ALA ALA ALA ALA ALA ALA ALA ALA ALA ALA ALA ALA ALA ALA ALA ALA ALA ALA ALA ALA ALA ALA ALA ALA ALA ALA ALA ALA ALA }

molaa = sequence { ALA ALA ALA ALA ALA ALA ALA ALA ALA ALA ALA ALA ALA ALA ALA ALA ALA ALA ALA ALA ALA ALA ALA ALA ALA ALA ALA ALA ALA ALA ALA ALA }

molaaa = sequence { ALA ALA ALA ALA ALA ALA ALA ALA ALA ALA ALA ALA ALA ALA ALA ALA ALA ALA ALA ALA ALA ALA ALA ALA ALA ALA ALA ALA ALA ALA ALA ALA }

molt = sequence { ALA ALA ALA ALA ALA ALA ALA ALA ALA ALA ALA ALA ALA ALA ALA ALA ALA ALA ALA ALA ALA ALA ALA ALA ALA ALA ALA ALA ALA ALA ALA ALA }

molp = sequence { ALA ALA ALA ALA ALA ALA ALA ALA ALA ALA ALA ALA ALA ALA ALA ALA ALA ALA ALA ALA ALA ALA ALA ALA ALA ALA ALA ALA ALA ALA ALA ALA }


impose mola { 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30 31 32  }  { { "N" "CA" "C" "N" -40.0 } { "C" "N" "CA" "C" -60.0 } }

impose molaa { 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30 31 32 }  { { "N" "CA" "C" "N" -45.0 } { "C" "N" "CA" "C" -60.0 } }

impose molaaa { 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30 31 32}  { { "N" "CA" "C" "N" -47.0 } { "C" "N" "CA" "C" -57.0 } }

impose molt { 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30 31 32  }  { { "N" "CA" "C" "N" -27.0 } { "C" "N" "CA" "C" -49.0 } }

impose molp { 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30 31 32  }  { { "N" "CA" "C" "N" -70.0 } { "C" "N" "CA" "C" -57.0 } }


savePDB mola alpha_-60_-40.pdb

savePDB molaa alpha_-60_-45.pdb

savePDB molaaa alpha_-57_-47.pdb

savePDB molt threeTen_-49_-27.pdb

savePDB molp pi_-57_-70.pdb

quit
EOF
$AMBERHOME/bin/tleap -f leap.in
