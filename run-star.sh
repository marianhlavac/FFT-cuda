mkdir -p bin
#$ -S /bin/sh
#$ -cwd
#$ -e .
#$ -o .
#$ -pe ompi 1
#$ -q gpu_02.q
./bin/fft-cuda-acc $1 bin/outp.dat $2