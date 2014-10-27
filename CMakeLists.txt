#Specify the version being used aswell as the language
cmake_minimum_required(VERSION 2.8)
#Name your project here
project(pigasus)

#set(PIGASUS_BUILD_DIR $ENV{PIGASUS_BUILD_DIR} )
set(PIGASUS_BUILD_DIR $ENV{PWD} )

file(MAKE_DIRECTORY ${PIGASUS_BUILD_DIR}/lib)

set(PIGASUS_DIRECTORY_TOOLS          ${PIGASUS_BUILD_DIR}/fortran/tools)
set(PIGASUS_DIRECTORY_TRACELOG       ${PIGASUS_BUILD_DIR}/fortran/tracelog)
set(PIGASUS_DIRECTORY_GEOMETRIES     ${PIGASUS_BUILD_DIR}/fortran/geometries)
set(PIGASUS_DIRECTORY_BLACKBOXES     ${PIGASUS_BUILD_DIR}/fortran/blackboxes)
set(PIGASUS_DIRECTORY_CONNECTIVITIES ${PIGASUS_BUILD_DIR}/fortran/connectivities)
set(PIGASUS_DIRECTORY_GRIDS          ${PIGASUS_BUILD_DIR}/fortran/grids)
set(PIGASUS_DIRECTORY_SPM            ${PIGASUS_BUILD_DIR}/fortran/spm)
set(PIGASUS_DIRECTORY_FEM            ${PIGASUS_BUILD_DIR}/fortran/fem)
set(PIGASUS_DIRECTORY_ASSEMBLY       ${PIGASUS_BUILD_DIR}/fortran/assembly)
set(PIGASUS_DIRECTORY_CORE           ${PIGASUS_BUILD_DIR}/python/core)

#file(MAKE_DIRECTORY ${PIGASUS_DIRECTORY_CORE})

include_directories(
	${PIGASUS_DIRECTORY_TOOLS} 
	${PIGASUS_DIRECTORY_TRACELOG}
	${PIGASUS_DIRECTORY_GEOMETRIES}
	${PIGASUS_DIRECTORY_BLACKBOXES}
	${PIGASUS_DIRECTORY_CONNECTIVITIES}
	${PIGASUS_DIRECTORY_GRIDS}
	${PIGASUS_DIRECTORY_SPM}
	${PIGASUS_DIRECTORY_FEM}
	${PIGASUS_DIRECTORY_ASSEMBLY})

set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${PIGASUS_BUILD_DIR}/lib)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${PIGASUS_BUILD_DIR}/lib)

##########################################################################
ADD_CUSTOM_TARGET( clear
   COMMAND find . -name *.mod -delete
   COMMAND find . -name *.a -delete
   COMMAND find . -name *.o -delete
   COMMAND find . -name *.tp -delete
   COMMAND find . -name *.cmake -delete
   COMMAND find . -name *.xmf -delete
   COMMAND find . -name *.h5 -delete
   COMMAND find . -name *.gnu -delete
   COMMAND find . -name *.dat -delete
   COMMAND find . -name *.mat -delete
   COMMAND find . -name *.mtx -delete
   COMMAND find . -name *.txt -delete
   COMMAND find . -name *.log -delete
   COMMAND rm -rf bin
   COMMAND rm -rf fortran 
   COMMAND rm -rf python 
   COMMAND rm -rf  build 

   WORKING_DIRECTORY ${CMAKE_BINARY_DIR}
   COMMENT "Delete all .a, .mod, .o, .tp, .cmake, .xmf, .h5, .gnu, .dat, .mat, .mtx, .txt, .log and folders /modules, /bin" VERBATIM)
##########################################################################

enable_language (Fortran)

##########################################################################
# FFLAGS depend on the compiler
##########################################################################
get_filename_component (Fortran_COMPILER_NAME ${CMAKE_Fortran_COMPILER} NAME)

if (Fortran_COMPILER_NAME MATCHES "gfortran.*")
  # gfortran
  set (CMAKE_Fortran_FLAGS_RELEASE "-O3 -fPIC")
  set (CMAKE_Fortran_FLAGS_DEBUG   "-O3 -fPIC -g")
elseif (Fortran_COMPILER_NAME MATCHES "ifort.*")
  # ifort (untested)
  set (CMAKE_Fortran_FLAGS_RELEASE "-f77rtl -O3 -fpic")
  set (CMAKE_Fortran_FLAGS_DEBUG   "-f77rtl -O0 -fpic -g")
elseif (Fortran_COMPILER_NAME MATCHES "g77")
  # g77
  set (CMAKE_Fortran_FLAGS_RELEASE "-funroll-all-loops -fno-f2c -O3 -m32")
  set (CMAKE_Fortran_FLAGS_DEBUG   "-fno-f2c -O0 -g -m32")
else (Fortran_COMPILER_NAME MATCHES "gfortran.*")
  message ("CMAKE_Fortran_COMPILER full path: " ${CMAKE_Fortran_COMPILER})
  message ("Fortran compiler: " ${Fortran_COMPILER_NAME})
  message ("No optimized Fortran compiler flags are known, we just try -O2...")
  set (CMAKE_Fortran_FLAGS_RELEASE "-O2")
  set (CMAKE_Fortran_FLAGS_DEBUG   "-O0 -g")
endif (Fortran_COMPILER_NAME MATCHES "gfortran.*")
##########################################################################



##########################################################################
add_subdirectory(fortran)
add_subdirectory(python)
##########################################################################
