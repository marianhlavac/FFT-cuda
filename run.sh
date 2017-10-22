#!/bin/sh

./compile.sh
./bin/fft-cuda $1 $2 > bin/outp.dat
python3 display-chart $1 bin/outp.dat
