cmake_minimum_required (VERSION 2.8.12)


set(MPI_ROOT "/uufs/chpc.utah.edu/sys/pkg/intel/ics/impi/5.0.1.035/intel64")
set(MPI_C_COMPILER "${MPI_ROOT}/bin/mpiicc")
set(MPI_C_INCLUDE_PATH "${MPI_ROOT}/include")
set(MPI_C_LIBRARIES "${MPI_ROOT}/lib/debug/libmpi.so")
set(MPI_CXX_INCLUDE_PATH "${MPI_ROOT}/include")
set(MPI_CXX_LIBRARIES "${MPI_ROOT}/lib/debug/libmpi.so")
set(MPI_Fortran_COMPILER "${MPI_ROOT}/bin/mpiifort")
set(MPI_Fortran_INCLUDE_PATH "${MPI_ROOT}/include")
set(MPI_Fortran_LIBRARIES "${MPI_ROOT}/lib/libmpifort.so" "${MPI_ROOT}/lib/libmpi.so")
set(MPIEXEC "mpirun")

project (mgac2)

find_package(MPI REQUIRED)

add_compile_options(-std=c++11)

include_directories(${MPI_INCLUDE_PATH})
include_directories("${PROJECT_SOURCE_DIR}/include")
include_directories("${PROJECT_SOURCE_DIR}/vendor/svl-1.5/include")
include_directories("${PROJECT_SOURCE_DIR}/vendor/svl-1.5/src")

#compile svl
SET( SVL_SRC ${PROJECT_SOURCE_DIR}/vendor/svl-1.5/src )
add_library(
svl
	${SVL_SRC}/LibSVL.cpp 
	${SVL_SRC}/Basics.cpp
)

#compile spacegroups.xml
add_custom_target(spacegroups.o
COMMAND cp ${PROJECT_SOURCE_DIR}/src/spacegroups.xml ${CMAKE_BINARY_DIR}
COMMAND objcopy --input binary --output elf64-x86-64 --binary-architecture i386 spacegroups.xml spacegroups.o )

#compile gasp2 library
add_library( 
gasp2
	src/gasp2.cpp 
	src/gasp2common.cpp
	src/gasp2param.cpp 
	src/gasp2pop.cpp 
	src/gasp2qe.cpp 
	src/gasp2struct.cpp
	src/tinyxml2.cpp
)
target_link_libraries( gasp2 ${MPI_LIBRARIES} )
target_link_libraries( gasp2 svl )
target_link_libraries( gasp2 ${CMAKE_BINARY_DIR}/spacegroups.o )
add_dependencies( gasp2 spacegroups.o )
add_dependencies( gasp2 svl )

#compile mgac.x program
add_executable(
mgac.x 
	src/mgac.cpp
)
target_link_libraries(
mgac.x
	gasp2 uuid
)
add_dependencies( mgac.x gasp2 )
