cmake_minimum_required(VERSION 3.20)

project(summa_pattern LANGUAGES Fortran)

set(CMAKE_Fortran_BIN_DIRECTORY ${PROJECT_BINARY_DIR}/bin)
set(CMAKE_Fortran_LIB_DIRECTORY ${PROJECT_BINARY_DIR}/lib)
set(CMAKE_Fortran_MODULE_DIRECTORY ${PROJECT_BINARY_DIR}/modules)

set(CMAKE_BUILD_TYPE "Debug" CACHE STRING "" FORCE)

# Compiler flags
set(GCC_DEBUG
     -g                                           # Generate symbols
     -fbacktrace                                  # symbolic stack traceback
     -ffpe-trap=invalid,zero,overflow,underflow   # control over floating-point exception
     -finit-real=nan                              #  All real scalars are initialised to NaN
     -fcheck=all                                  # Enable all run-time test of -fcheck: array-temps, bits, bounds, do, mem, pointer, recursion
    )
           
set(GCC_RELEASE
    -O3
   )

# Interface to control scope of compiler flags and options
add_library(common_flags INTERFACE)
# Build-specific flags
target_compile_options(common_flags INTERFACE 
   $<$<AND:$<CONFIG:Debug>,$<Fortran_COMPILER_ID:GNU>>:${GCC_DEBUG}>
   $<$<AND:$<CONFIG:Release>,$<Fortran_COMPILER_ID:GNU>>:${GCC_RELEASE}>
   )

# External libraries
find_package(MPI 3 COMPONENTS Fortran)
find_package(BLAS REQUIRED)
find_package(LAPACK REQUIRED)

add_executable(summa)
target_compile_options(summa PRIVATE)
target_link_libraries(summa PRIVATE 
   MPI::MPI_Fortran 
   BLAS::BLAS 
   LAPACK::LAPACK
   common_flags)

add_executable(summa_overlapping)
target_compile_options(summa_overlapping PRIVATE)
target_link_libraries(summa_overlapping PRIVATE 
   MPI::MPI_Fortran 
   BLAS::BLAS 
   LAPACK::LAPACK
   common_flags)

add_subdirectory(src)
