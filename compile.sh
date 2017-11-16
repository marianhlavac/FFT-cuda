#!/bin/sh

if [ "$ENV_STAR" = "1" ]
then
  echo "Compiling binaries for star cluster..."
  pgc++ -acc -ta=nvidia,time -Minfo=accel -o bin/fft-cuda-acc src/fft-cuda.cpp
else
  echo "Compiling tool for creating data..."
  g++ -o bin/create-data src/create-data.cpp 
fi

echo "Compiling FFT-cuda binaries..."
g++ -o bin/fft-cuda src/fft-cuda.cpp 


echo "Compiled."