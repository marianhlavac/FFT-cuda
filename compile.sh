#!/bin/sh

echo "Compiling src/fft-cuda.cu"
nvcc src/fft-cuda.cu -o bin/fft-cuda

echo "Compiled."
