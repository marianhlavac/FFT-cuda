#!/bin/sh

g++ src/fft-cuda.cpp -o bin/fft-cuda -g
g++ src/create-data.cpp -o bin/create-data