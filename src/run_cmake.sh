#!/bin/sh

# Change this to change build type
buildtype=Debug
#buildtype=Release
#buildtype=RelWithDebInfo
#buildtype=MinSizeRel
#buildtype=gcov
#buildtype=Callgrind

# Always build CShore
echo "Building all versions of the CShore library"
echo ""
cd cshore
./make_cshore_lib.sh
cd ..

# Now run CMake for CoastalME
echo ""
echo "================================================================="
echo ""
echo "CoastalME: starting CMake for Linux (using gcc, $buildtype build)"

rm -f CMakeCache.txt
#cmake -G "Unix Makefiles" -DCMAKE_BUILD_TYPE=$buildtype .
cmake -DCMAKE_BUILD_TYPE=$buildtype . -G"CodeBlocks - Unix Makefiles"

echo ""
echo "================================================================="
