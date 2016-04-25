#!/bin/bash

#make directories only if they do not already exist
echo "making two folders: ./bin ./lib"
mkdir -p bin
mkdir -p lib

echo "Moving into the src directory to start building"
cd ./src


#now compile the library
echo "Making Lawson-s SVD routine linkable"
gfortran -c liblawson.f90
ar rcv liblawson.a liblawson.o
ranlib liblawson.a
rm liblawson.o

echo "Moving the linkable library liblawson.a into the lib folder"
mv liblawson.a ../lib/

#Compile the subroutines first...
echo "Making Kerr.f90 -- this code handles errors..."
gfortran -g -fbacktrace -fbounds-check -Wall -Wextra -Wconversion -pedantic -fcheck=all -c Kerr.f90\
 -o ../lib/Kerr.o
echo "Making Kopen.f90 -- this code handles file io..."
gfortran -g -fbacktrace -fbounds-check -Wall -Wextra -Wconversion -pedantic -fcheck=all -c Kopen.f90\
 -o ../lib/Kopen.o
echo "Making Kpitchang.f90 -- this code calculates helical pitch and sweep..."
gfortran -g -fbacktrace -fbounds-check -Wall -Wextra -Wconversion -pedantic -fcheck=all -c Kpitchang.f90\
 -o ../lib/Kpitchang.o
echo "Making Kprinter.f90 -- this code handles how information is formatted when writing to disk..."
gfortran -g -fbacktrace -fbounds-check -Wall -Wextra -Wconversion -pedantic -fcheck=all -c Kprinter.f90\
 -o ../lib/Kprinter.o
echo "Making Kprint_optHelix.f90 -- this code handles formatted printing of the optimal helical axis, in PDB - RCSB.org format..."
gfortran -g -fbacktrace -fbounds-check -Wall -Wextra -Wconversion -pedantic -fcheck=all -c Kprinter_optHelix.f90\
 -o ../lib/Kprinter_optHelix.o
echo "Making Kread_control_file.f90 -- this code handles user input"
gfortran -g -fbacktrace -fbounds-check -Wall -Wextra -Wconversion -pedantic -fcheck=all -c Kread_control_file.f90\
 -o ../lib/Kread_control_file.o
echo "Making Kreadpdb.f90 -- this code handles PDB - RCSB.org input file reading..."
gfortran -g -fbacktrace -fbounds-check -Wall -Wextra -Wconversion -pedantic -fcheck=all -c Kreadpdb.f90\
 -o ../lib/Kreadpdb.o
echo "Making Kreadtraj.f90 -- this code handles Amber - AmberMD.org input file reading..."
gfortran -g -fbacktrace -fbounds-check -Wall -Wextra -Wconversion -pedantic -fcheck=all -c Kreadtraj.f90\
 -o ../lib/Kreadtraj.o
echo "Making Krotate.f90 -- this code applies Evangelos Coutsias- rotation algebra..."
gfortran -g -fbacktrace -fbounds-check -Wall -Wextra -Wconversion -pedantic -fcheck=all -c Krotate.f90\
 -o ../lib/Krotate.o

#Compile in debug mode...
echo "Making Helios... ... ..."
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

#<< 'commentOUT1'
#-----------------------------------
echo
echo
echo -
echo - - -
#-----------------------------------
#commentOUT1

<< 'commentOUT2'
#-----------------------------------
export KHELIXHOME=$(pwd)
export PATH=$PATH:$KHELIXHOME
#-----------------------------------
<< commentOUT2

#<< 'commentOUT3'
#-----------------------------------
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
#-----------------------------------
#commentOUT3

#<<'commentOUT4'
#-----------------------------------
cat > Ktest.input << EOF
ktest 1
EOF
$KHELIXHOME/bin/heliosv2.o Ktest.input
#-----------------------------------
#<< commentOUT4

#<< 'commentOUT5'
#-----------------------------------
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
#-----------------------------------
#commentOUT5
