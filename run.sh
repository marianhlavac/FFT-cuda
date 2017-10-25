#!/bin/sh

mkdir -p bin
./compile.sh
./bin/fft-cuda $1 $2 > bin/outp.dat
python3 display-chart.py $1 bin/outp.dat
