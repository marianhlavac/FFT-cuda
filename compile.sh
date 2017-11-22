#!/bin/sh

if [ "$ENV_STAR" = "1" ]
then
  echo "Compiling binaries for star cluster..."
  pgc++ -acc -fast -Minfo=accel -o bin/fft-cuda-acc -std=c++11 src/fft-cuda.cpp
else
  echo "Compiling tool for creating data..."
  g++ -o bin/create-data -std=c++11 src/create-data.cpp 
fi

if [ "$ENV_LOCAL" = "1" ]
then
  echo "Compiling OpenACC binaries on local machine..."
  /opt/pgi/osx86-64/2017/bin/pgc++ -ta=multicore -fast -Minfo=all,ccff -o bin/fft-cuda-acc src/fft-cuda.cpp
else
  echo "Compiling FFT-cuda binaries..."
  g++ -o bin/fft-cuda -std=c++11 src/fft-cuda.cpp -g
fi

echo "Compiled."
