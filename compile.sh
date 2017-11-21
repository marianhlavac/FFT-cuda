#!/bin/sh

if [ "$ENV_STAR" = "1" ]
then
  echo "Compiling binaries for star cluster..."
  cp src/fft-cuda.cpp src/fft-cuda-acc.cpp
  echo -e "#define STAR\n$(cat src/fft-cuda-acc.cpp)" > src/fft-cuda-acc.cpp
  pgc++ -acc -ta=nvidia:managed,time -Minfo=accel -o bin/fft-cuda-acc src/fft-cuda-acc.cpp
else
  echo "Compiling tool for creating data..."
  g++ -o bin/create-data src/create-data.cpp 
fi

if [ "$ENV_LOCAL" = "1" ]
then
  echo "Compiling OpenACC binaries on local machine..."
  /opt/pgi/osx86-64/2017/bin/pgc++ -ta=multicore -fast -Minfo=all,ccff -o bin/fft-cuda-acc src/fft-cuda.cpp
else
  echo "Compiling FFT-cuda binaries..."
  g++ -o bin/fft-cuda src/fft-cuda.cpp -g
fi

echo "Compiled."