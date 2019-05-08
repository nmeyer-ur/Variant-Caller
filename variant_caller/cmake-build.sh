#!/bin/bash

# compiler -- gcc recommended
CXX=g++-7
#CXX=clang-5.0

# directories
BUILD_DIR="build"
HTSLIB_INC=""
HTSLIB_LIB=""
EIGEN_INC=""

# build
mkdir -p ${BUILD_DIR}
cd ${BUILD_DIR}

cmake .. -DHTSLIB_INC=${HTSLIB_INC} -DHTSLIB_LIB=${HTSLIB_LIB} -DEIGEN_INC=${EIGEN_INC} -DCMAKE_CXX_COMPILER=${CXX}
make

cd ..

cp ${BUILD_DIR}/ec .

#./cmake-clean.sh
