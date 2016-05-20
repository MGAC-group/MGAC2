
================================================================================
INSTALL INSTRUCTIONS:

Requires:

-cmake 2.8.12 or greater
-mpich2 3.1.4 or greater
-gcc 4.9.2 or greater (C++11 support)
-QE 5.0.2 (other versions NOT tested)
-objcopy (should be standard on linux)

Intel compilers/MPI optional but not required
-Intel compilers v2015.1.133 or greater
-IMPI 5.0.1.035 or greater

On generic systems without Intel compilers:

mkdir mgac-build
cd mgac-build
cmake -D CMAKE_CXX_COMPILER=g++ -D CMAKE_C_COMPILER=gcc ../mgac-redux
make

MAKE SURE TO SOURCE MPICH2 VARS!

On CHPC systems with Intel compilers:

module load cmake intel impi gcc/4.9.2
mkdir mgac-build
cd mgac-build
cmake -D CMAKE_CXX_COMPILER=icpc -D CMAKE_C_COMPILER=icc ../mgac-redux
make

MODULE LOAD ORDER IS IMPORTANT! gcc/4.9.2 must be loaded last

To run make sure that the gcc/4.9.2 libraries are accessible. If using mpi 3.1.4 or greater mgac can use either mpich2 or impi to run. 

================================================================================
RUNNING NOTES:

To see the help:
mgac.x -h
mgac.x --help

A simple run is executed by:
mgac.x -i inputfile

A restart is the same, with a reference to the sql file:
mgac.x -i inputfile -r restart.sq3

A list of spacegroups can be listed by:
mgac.x -l 

To convert an output sqlite file to CIF format. This automatically sorts structures:
mgac.x -i output.sq3 -c file.cif

The number of best structures for CIF conversion defaults to 100. Use the -s flag to specify a different number:
mgac.x -i output.sq3 -c file.cif -s 200

A different table can be selected from the default "structs" table using -b:
mgac.x -i output.sq3 -c file.cif -s 100 -b precluster

An input template can be quickly generated from a CIF using -t:
mgac.x -i input.cif -t template.xml

The input template still requires all other parameters to be manually added! 

The plane for molecules must be specified using -p in comma delimited form using the atom labels specified in the CIF:
mgac.x -i input.cif -t template.xml -p C1,C2,N1

NOTE: the template generation does not include bond/angle/limitations or dihedrals. Those must be specified manually still.

================================================================================
INPUT FORMAT


================================================================================
DEVELOPMENT NOTES:

A hook for steady state methods exists in gasp2.cpp at line 1304, This is where it should likely be implemented, although further changes to the GASPcontrol serverprogram will be needed.

SQL stuff: 



Outlier filter testing: There is a function designed to test if a structure energy is way out of range. (GASP2pop.cpp, line 1496). It isn't working right. 

