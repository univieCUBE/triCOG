#!/bin/bash

# Make sure we are in the correct directory
cd "$(dirname "$0")"

# Remove old build files
echo "Cleaning old build files..."
rm -rf build
mkdir build

# Configure the build system
echo "Configure the build system..."
cmake -S . -B build

# Compile build in parallel
echo "Compile the project..."
cmake --build build -j $(nproc 2>/dev/null || sysctl -n hw.ncpu 2>/dev/null || echo 1)

echo "Finished compilation!"