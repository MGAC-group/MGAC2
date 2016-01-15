INSTALL INSTRUCTIONS:

Requires:

-cmake 2.8.12 or greater
-mpich2 3.1.4 or greater
-gcc 4.9.2 or greater (C++11 support)
-QE 5.0.2 (other versions NOT tested)

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

