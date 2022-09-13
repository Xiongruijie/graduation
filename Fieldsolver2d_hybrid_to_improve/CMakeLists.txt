##
#  CMake script
##

CMAKE_MINIMUM_REQUIRED(VERSION 3.0)

# Set the name of the project and target:
SET(TARGET "fieldsolver2d")
SET(CHOLMOD_HEADER_PATH "/usr/include/suitesparse")
SET(CHOLMOD_LIBRARY_PATH "/usr/lib/x86_64-linux-gnu")
SET(EIGEN_HEADER_PATH "/usr/include/eigen3")
#SET(COMMON_HEADER_PATH "../COMMON")
SET(SRC_HEADER_PATH "./src")
SET(FRW_HEADER_PATH "./src/frw")

# Project name
project(${TARGET})

# Compiler options
add_definitions("-Wall -g")

# Declare all source files the target consists of:
SET(TARGET_SRC
  src/${TARGET}.cc
  # You can specify additional files here!
  )

# Add header path
#include_directories(${CHOLMOD_HEADER_PATH} ${COMMON_HEADER_PATH} ${SRC_HEADER_PATH})
include_directories(${CHOLMOD_HEADER_PATH} ${EIGEN_HEADER_PATH} ${SRC_HEADER_PATH} ${FRW_HEADER_PATH})

# Add library path
link_directories(${CHOLMOD_LIBRARY_PATH})

# Add link libraries
link_libraries(cholmod amd colamd suitesparseconfig pthread)

# Generate executable file
add_executable(${TARGET} ${TARGET_SRC})
