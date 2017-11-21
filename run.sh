mkdir -p bin
./compile.sh
./bin/fft-cuda $1 bin/outp.dat $2 
python3 display-chart.py $1 bin/outp.dat