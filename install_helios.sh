#!/bin/bash

#make directories only if they do not already exist
mkdir -p bin
mkdir -p src
mkdir -p lib

cd ./src


#now compile the library
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
gfortran -g -fbacktrace -fbounds-check -Wall -Wextra -Wconversion -pedantic -fcheck=all helios.f90 \
../lib/Kerr.o \
../lib/Kopen.o \
../lib/Kpitchang.o \
../lib/Kprinter.o \
../lib/Kprinter_optHelix.o \
../lib/Kread_control_file.o \
../lib/Kreadpdb.o \
../lib/Kreadtraj.o \
../lib/Krotate.o \
-L../lib -llawson -o ../bin/helios.o

cd ../

<< 'commentOUT1'
echo
echo
echo -
echo - - -
commentOUT1

#-----------------------------------
export KHELIXHOME=$(pwd)
export PATH=$PATH:$KHELIXHOME
#-----------------------------------

#<< 'commentOUT3'
echo Please add this to your bashrc '(or bash_profile)':
echo ----------------------------------------------
echo export KHELIXHOME=$(pwd)
echo 'export PATH=$PATH:$KHELIXHOME'
echo ----------------------------------------------
echo 
echo -
echo - - -
echo By executing this script, the above paths have already 
echo  been exported to the shell. You may have to open a 
echo   new shell, or "set", if you-re on a Mac '8^)'
echo
echo -
echo - - -
echo Now the script will run kHelix to Hello World:
echo 
echo -
echo - - - 
#commentOUT3

#-----------------------------------
cat > Ktest.input << EOF
ktest 1
EOF
$KHELIXHOME/bin/heliosv2.o Ktest.input
#-----------------------------------

#<< 'commentOUT2'
echo 
echo -
echo - - -
echo The examples in the test_1, test_2, test_3 directories
echo  have simple bash scripts that generate a kHelix input
echo   file and then execute based on that input file. The
echo    example directories have multiple coordinates files
echo     that will be looped over by the run script.
echo Please play around with the run scripts!
echo If you are having problems that Google can-t help you 
echo  with, please contact Kevin at 84hauser@gmail.com....
echo Rock on...
#commentOUT2
