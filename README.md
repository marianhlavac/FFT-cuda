# FFT-cuda
Fast Fourier Transform implementation, computable on CUDA computing platform

> 😵 Currenty in WIP state, computing just sequentially. CUDA support coming soon. Ain't nobody got time for dat.

## How to run
Requires **g++** (duh) and **python 3** to be installed. Also some Python dependencies specified in `requirements.txt`.

I probably hate Makefile. Also too lazy to set-up CMake. You can hate me too.

Install Python dependencies

```
 $ pip3 install -r requirements.txt
```

Then the run script does everything.

```
 $ run.sh [input_data_file] [sampling_rate]
```