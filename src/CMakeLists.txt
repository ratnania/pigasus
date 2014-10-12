#Specify the version being used aswell as the language
cmake_minimum_required(VERSION 2.6)
#Name your project here
project(pigasus)

#set(PIGASUS_BUILD_DIR $ENV{PIGASUS_BUILD_DIR} )
set(PIGASUS_BUILD_DIR $ENV{PWD} )

file(MAKE_DIRECTORY ${PIGASUS_BUILD_DIR}/lib)

set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${PIGASUS_BUILD_DIR}/lib)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${PIGASUS_BUILD_DIR}/lib)

enable_language (Fortran)

##########################################################################
# FFLAGS depend on the compiler
##########################################################################
get_filename_component (Fortran_COMPILER_NAME ${CMAKE_Fortran_COMPILER} NAME)

if (Fortran_COMPILER_NAME MATCHES "gfortran.*")
  # gfortran
  set (CMAKE_Fortran_FLAGS_RELEASE "-O3 -fPIC")
  set (CMAKE_Fortran_FLAGS_DEBUG   "-O3 -fPIC -g")
  set (CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS}  -fbounds-check")
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
# DEBUG DETAIL FLAG
##########################################################################
set (CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS}  -D_TRACE -D_DEBUG")
message("CMAKE_Fortran_FLAGS:   ${CMAKE_Fortran_FLAGS}")
##########################################################################

##########################################################################
add_subdirectory(tools)
add_subdirectory(tracelog)
add_subdirectory(geometries)
add_subdirectory(blackboxes)
add_subdirectory(connectivities)
add_subdirectory(grids)
add_subdirectory(spm)
add_subdirectory(fem)
add_subdirectory(assembly)
##########################################################################
