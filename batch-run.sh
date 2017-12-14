echo "Compiling FFT-cuda to run the latest version..."
./compile.sh

echo "Computing folder 50Hz..."
./bin/fft-cuda data/50Hz > outputs.csv

echo "Computing folder 50Hz+500Hz..."
./bin/fft-cuda data/50Hz+500Hz | tail -n +2 >> outputs.csv

echo "Computing folder 50Hz+505Hz+12000Hz..."
./bin/fft-cuda data/50Hz+505Hz+12000Hz | tail -n +2 >> outputs.csv

echo "Done."