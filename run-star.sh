mkdir -p bin
#$ -S /bin/sh
#$ -cwd
#$ -e .
#$ -o .
#$ -pe ompi 1
#$ -q gpu_02.q
./bin/fft-cuda-acc data/50Hz/
./bin/fft-cuda-acc data/50Hz+500Hz/
./bin/fft-cuda-acc data/500Hz+505Hz+12000Hz
