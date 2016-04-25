#!/bin/bash

gfortran -c liblawson.f90
ar rcv liblawson.a liblawson.o
ranlib liblawson.a
rm liblawson.o

mv liblawson.a ../lib/

#Compile the subroutines first...
gfortran -g -fbacktrace -fbounds-check -Wall -Wextra -Wconversion -pedantic -fcheck=all -c Kerr.f90\
 -o ../lib/Kerr.o
gfortran -g -fbacktrace -fbounds-check -Wall -Wextra -Wconversion -pedantic -fcheck=all -c Kopen.f90\
 -o ../lib/Kopen.o
gfortran -g -fbacktrace -fbounds-check -Wall -Wextra -Wconversion -pedantic -fcheck=all -c Kpitchang.f90\
 -o ../lib/Kpitchang.o
gfortran -g -fbacktrace -fbounds-check -Wall -Wextra -Wconversion -pedantic -fcheck=all -c Kprinter.f90\
 -o ../lib/Kprinter.o
gfortran -g -fbacktrace -fbounds-check -Wall -Wextra -Wconversion -pedantic -fcheck=all -c Kprinter_optHelix.f90\
 -o ../lib/Kprinter_optHelix.o
gfortran -g -fbacktrace -fbounds-check -Wall -Wextra -Wconversion -pedantic -fcheck=all -c Kread_control_file.f90\
 -o ../lib/Kread_control_file.o
gfortran -g -fbacktrace -fbounds-check -Wall -Wextra -Wconversion -pedantic -fcheck=all -c Kreadpdb.f90\
 -o ../lib/Kreadpdb.o
gfortran -g -fbacktrace -fbounds-check -Wall -Wextra -Wconversion -pedantic -fcheck=all -c Kreadtraj.f90\
 -o ../lib/Kreadtraj.o
gfortran -g -fbacktrace -fbounds-check -Wall -Wextra -Wconversion -pedantic -fcheck=all -c Krotate.f90\
 -o ../lib/Krotate.o

#Compile in debug mode...
gfortran -g -fbacktrace -fbounds-check -Wall -Wextra -Wconversion -pedantic -fcheck=all kHelixv2.f90 \
../lib/Kerr.o \
../lib/Kopen.o \
../lib/Kpitchang.o \
../lib/Kprinter.o \
../lib/Kprinter_optHelix.o \
../lib/Kread_control_file.o \
../lib/Kreadpdb.o \
../lib/Kreadtraj.o \
../lib/Krotate.o \
-L../lib -llawson -o kHelixv2.o

#-----------------------------------
cat > Ktest.input << EOF
ktest 1
EOF
$KHELIXHOME/kHelixv2.o Ktest.input
#-----------------------------------
