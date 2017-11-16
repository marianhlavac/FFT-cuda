#!/bin/sh

if [ "$ENV_STAR" = "1" ]
then
  echo "Compiling binaries for star cluster..."
  pgc++ -acc -ta=nvidia,time -Minfo=accel -o bin/fft-cuda src/fft-cuda.cpp
fi

echo "Compiling binaries..."
g++ -o bin/fft-cuda src/fft-cuda.cpp 
g++ -o bin/create-data src/create-data.cpp 

echo "Compiled."