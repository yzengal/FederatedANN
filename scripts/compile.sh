#!/bin/bash

DIR_PATH="../src/build"

if [ -d "$DIR_PATH" ]; then
  echo "Directory $DIR_PATH already exists. Removing it..."
  rm -rf "$DIR_PATH"
fi

mkdir -p ../src/build
cd ../src/build
cmake .. -DUSE_SGX=OFF -DSGX_HW=OFF -DLOCAL_DEBUG=OFF
make -j 4

cd ../../scripts

