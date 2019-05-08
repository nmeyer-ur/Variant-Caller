# cmake version
cmake_minimum_required(VERSION 2.8)

# project
project(ec)

# sources
set(SOURCES src/ec.cpp)
add_executable(ec ${SOURCES})

# compiler
set(CMAKE_CXX_COMPILER "g++")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++14 -Wall -Wextra -pedantic -funsafe-loop-optimizations -Ofast -march=native -fmax-errors=5")
#set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++14 -Wall -Wextra -pedantic -funsafe-loop-optimizations -Ofast -march=native -fmax-errors=5 -g")

# zlib
find_package(ZLIB)
if (ZLIB_FOUND)
  include_directories(${ZLIB_INCLUDE_DIRS})
  target_link_libraries(ec ${ZLIB_LIBRARIES})
endif()

# headers and libraries directories
set(HTSLIB_INC "" CACHE FILEPATH "Path to htslib includes")
set(HTSLIB_LIB "" CACHE FILEPATH "Path to htslib library")
set(EIGEN_INC "" CACHE FILEPATH "Path to Eigen includes")

# include directories
include_directories(
"include/"
${HTSLIB_INC}
${EIGEN_INC}
)

# library directory
link_directories(
"/home/nils/genomics/Fraunhofer/error_correction/c/htslib"
)

# search htslib in default path
find_library(
HTS
NAMES "hts" "libhts" "libhts.so"
)

# search htslib in alternative path
find_library(
HTS
NAMES "hts" "htslib" "htslib.so"
PATHS ${HTSLIB_LIB}
NO_DEFAULT_PATH
)

if(NOT HTS)
  message(FATAL_ERROR "-- htslib not found")
else()
  message("-- Found htslib: ${HTS}")
  target_link_libraries(ec ${HTS})
endif()

# OpenMP support
find_package(OpenMP)

if (OPENMP_FOUND)
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
else()
  message("-- Threading will be disabled")
endif()
