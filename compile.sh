#!/bin/sh

echo "Compiling src/fft-cuda.cu"
nvcc src/fft-cuda.cu -std=c++11 -o bin/fft-cuda

echo "Compiled."
