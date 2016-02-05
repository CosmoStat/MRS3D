cmake_minimum_required(VERSION 2.6)

set(CMAKE_Fortran_COMPILER_INIT gfortran)

project(MRS3D)
enable_language (Fortran)

# Checking environment variables
SET(IDL $ENV{IDL_DIR})
SET(CUSTOM_INSTALL $ENV{INSTALL_DIR})

# idl dependencies
if(IDL)
include_directories("$ENV{IDL_DIR}/external/include")
if(UNIX)
link_directories("$ENV{IDL_DIR}/bin/bin.linux.x86_64")
elseif(APPLE)
link_directories("$ENV{IDL_DIR}/bin/bin.darwin.x86_64")
endif(UNIX)
endif(IDL)

link_directories($ENV{CFITSIO_DIR})
set(CFITSIO_LIBRARIES "-lcfitsio")


########################### HEALPIX INFO  ########################### 

#
# HEALPIX cxx support
#
if (DEFINED ENV{HEALPIX})
# HEALPIX Fortran support
include_directories($ENV{HEALPIX}/include)
link_directories($ENV{HEALPIX}/lib)
include_directories($ENV{HEALPIX}/src/cxx/$ENV{HEALPIX_TARGET}/include)
# This makes sure cmake won't be looking for the libraries elsewhere
find_library(HEALPIX_CXX NAMES libhealpix_cxx.a PATHS $ENV{HEALPIX}/src/cxx/$ENV{HEALPIX_TARGET}/lib NO_DEFAULT_PATH)
find_library(CXX_SUPPORT NAMES libcxxsupport.a  PATHS $ENV{HEALPIX}/src/cxx/$ENV{HEALPIX_TARGET}/lib NO_DEFAULT_PATH)
find_library(PSHT NAMES libsharp.a PATHS $ENV{HEALPIX}/src/cxx/$ENV{HEALPIX_TARGET}/lib NO_DEFAULT_PATH)
find_library(FFTPACK NAMES libfftpack.a PATHS $ENV{HEALPIX}/src/cxx/$ENV{HEALPIX_TARGET}/lib NO_DEFAULT_PATH)
find_library(C_UTILS NAMES libc_utils.a PATHS $ENV{HEALPIX}/src/cxx/$ENV{HEALPIX_TARGET}/lib NO_DEFAULT_PATH)
set(HEALPIX_LIBRARIES ${HEALPIX_CXX} ${CXX_SUPPORT} ${PSHT} ${FFTPACK} ${C_UTILS})
else(DEFINED ENV{HEALPIX})
message("==================================================")
message(" ERROR : HEALPIX NOT FOUND.")
message("==================================================")
endif (DEFINED ENV{HEALPIX})



set(CMAKE_CXX_FLAGS "-O3 -fomit-frame-pointer -fno-common -fPIC -fopenmp")
set(CMAKE_Fortran_FLAGS "-O3 -DGFORTRAN -fno-second-underscore -fdefault-integer-8 -fPIC -frecursive -fopenmp")


#
# Build FastDSBT Library and IDL interface
#
#
add_library(fastDSBT STATIC src/cxx/fastDSBT.cpp src/cxx/sbtools.cpp src/f90/rfsfns.f src/f90/sbtools.f90)
target_link_libraries(fastDSBT ${HEALPIX_LIBRARIES} ${CFITSIO_LIBRARIES} "-lgfortran")
if(IDL)
add_library(fastDSBT_IDL SHARED src/cxx/fastDSBT_IDL.cpp)
target_link_libraries(fastDSBT_IDL fastDSBT)
endif(IDL)


# 
# Build MRS3D Library and IDL interface
#
add_library(mrs3d STATIC src/cxx/almn_fitsio.cpp src/cxx/almn_cln_tools.cpp src/cxx/cln_fitsio.cpp src/cxx/field_fitsio.cpp)
target_link_libraries(mrs3d ${HEALPIX_LIBRARIES} fastDSBT ${CFITSIO_LIBRARIES})
if(IDL)
add_library(mrs3d_IDL SHARED src/cxx/MRS3D_IDL.cpp)
target_link_libraries(mrs3d_IDL mrs3d)
endif(IDL)

#
# Build MRS3D Tools
#
set(CMAKE_C_COMPILER mpicc)
set(CMAKE_CXX_COMPILER mpicxx)

add_executable(gen_qln src/cxx/gen_qln.cpp)
target_link_libraries(gen_qln mrs3d)

add_executable(almn2wavelet src/cxx/almn2wavelet.cpp)
target_link_libraries(almn2wavelet mrs3d)

add_executable(waveletThresholding src/cxx/waveletThresholding.cpp)
target_link_libraries(waveletThresholding mrs3d)

add_executable(field2almn src/cxx/field2almn.cpp)
target_link_libraries(field2almn mrs3d)

add_executable(almn2field src/cxx/almn2field.cpp)
target_link_libraries(almn2field mrs3d)

add_executable(backprojection src/cxx/backprojection.cpp)
target_link_libraries(backprojection mrs3d)


# add a target to generate API documentation with Doxygen
find_package(Doxygen)
if(DOXYGEN_FOUND)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/Doxyfile.in ${CMAKE_CURRENT_BINARY_DIR}/Doxyfile @ONLY)
add_custom_target(doc ${DOXYGEN_EXECUTABLE} ${CMAKE_CURRENT_BINARY_DIR}/Doxyfile
WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
COMMENT "Generating API documentation with Doxygen" VERBATIM
)
endif(DOXYGEN_FOUND)

#
# Install (by default in the project directory)
#
set(CMAKE_INSTALL_PREFIX ${PROJECT_SOURCE_DIR})
if(CUSTOM_INSTALL)
set(CMAKE_INSTALL_PREFIX $ENV{INSTALL_DIR})
endif(CUSTOM_INSTALL)

install(TARGETS almn2wavelet waveletThresholding field2almn almn2field gen_qln DESTINATION bin)
install(TARGETS mrs3d fastDSBT DESTINATION lib)
if(IDL)
install(TARGETS mrs3d_IDL fastDSBT_IDL DESTINATION lib)
endif(IDL)

FILE(GLOB files "${PROJECT_SOURCE_DIR}/src/cxx/*.h")
INSTALL(FILES ${files} DESTINATION include)