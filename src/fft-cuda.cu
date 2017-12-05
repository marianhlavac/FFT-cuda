#include <iostream>
#include <fstream>
#include <sstream>
#include <chrono>
#include <vector>
#include <cmath>
#include <dirent.h>
#include <cstring>

using namespace std;

#define CUDA_N_THR 128

// Complex numbers data type
typedef float2 Cplx;

static __device__ __host__ inline Cplx CplxAdd(Cplx a, Cplx b) {
  Cplx c; c.x = a.x + b.x; c.y = a.y + b.y; return c;
}

static __device__ __host__ inline Cplx CplxInv(Cplx a) {
  Cplx c; c.x = -a.x; c.y = -a.y; return c;
}

static __device__ __host__ inline Cplx CplxMul(Cplx a, Cplx b) {
  Cplx c; c.x = a.x * b.x - a.y + b.y; c.y = a.x * b.y + a.y * b.x; return c;
}

static __device__ __host__ inline Cplx CplxMul(float a, Cplx b) {
  Cplx c; c.x = a * b.x; c.y = a * b.y; return c;
}

__global__ void bitrev_reorder(Cplx* __restrict__ r, Cplx* __restrict__ d, int s) {
  int id = blockIdx.x * CUDA_N_THR + threadIdx.x;
  int j = id, k = 0;
  
  // ASK: Nejde reverse bitu GPUckem?
  for (int l = 1; l <= s; l++) {
    k = k * 2 + (j & 1);
    j >>= 1;
  }
  
  r[k] = d[id];
}

__global__ void inplace_fft(Cplx* __restrict__ r, int m, int n) {
  int j = (blockIdx.x * CUDA_N_THR + threadIdx.x) * m;
  
  for (int k = 0; k < m / 2; k++) {
    if (j + k + m / 2 < n) { 
      Cplx t, u;
      
      t.x = __cosf((2.0 * M_PI * k) / (1.0 * m));
      t.y = -__sinf((2.0 * M_PI * k) / (1.0 * m));
      
      u = r[j + k];
      t = CplxMul(t, r[j + k + m / 2]);

      r[j + k] = CplxAdd(u, t);
      r[j + k + m / 2] = CplxAdd(u, CplxInv(t));
    }
  }
}

/**
 * Runs in-place Iterative Fast Fourier Transformation.
 */
void fft(Cplx* __restrict__ d, size_t n) {
  Cplx *r, *dn;
  size_t data_size = n * sizeof(Cplx);
  
  // Copy data to GPU
  cudaMalloc((void**)&r, data_size);
  cudaMalloc((void**)&dn, data_size);
  cudaMemcpy(dn, d, data_size, cudaMemcpyHostToDevice);
  
  // Bit-reversal reordering
  int s = log2(n);
  bitrev_reorder<<<ceil(n / CUDA_N_THR), CUDA_N_THR>>>(r, dn, s);
  
  // Synchronize
  cudaDeviceSynchronize();
  
  // Iterative FFT
  for (int i = 1; i <= s; i++) {
    int m = 1 << i;
    inplace_fft<<<ceil(n / m / CUDA_N_THR), CUDA_N_THR / m>>>(r, m, n);
  }
  
  // Copy data from GPU
  Cplx* result;
  result = (Cplx*)malloc(data_size / 2);
  cudaMemcpy(result, r, data_size / 2, cudaMemcpyDeviceToHost);
}

/**
 * Reads numeric data from a file.
 */
void read_file(const char* filename, vector<Cplx>& out) {
  ifstream file;
  file.open(filename);
  
  if (file.is_open()) {
    while (!file.eof()) {
      Cplx val;
      if (file >> val.x) {
        val.y = 0;
        out.push_back(val);
      }
    }
  } else {
    cerr << "Can't open file " << filename << " for reading." << endl;
  }
  
  file.close();
}

void save_results(const char* filename, Cplx* result, size_t count, int sample_rate) {
  char* outfilename = new char[512];
  strcpy(outfilename, filename);
  strcat(outfilename, ".out");
  ofstream outfile;
  outfile.open (outfilename);
  outfile.precision(4);
  outfile << "frequency, value" << endl;
  for (int i = 0; i < count / 2; i++) {
      outfile << i * ((float)sample_rate/count) << "," << result[i].x << endl;
  }
  outfile.close();
}

void compute(Cplx* buffer, size_t count, int sample_rate, const char* filename) {
  // Start the stopwatch
  // auto start = chrono::high_resolution_clock::now();
  
  // Run FFT algorithm with loaded data
  fft(buffer, count);
  
  // Log the elapsed time
  // auto finish = chrono::high_resolution_clock::now();
  // auto microseconds = chrono::duration_cast<std::chrono::microseconds>(finish-start);
  // 
  // cout << "elapsed" << microseconds.count() << "ms" << endl;
  
  // Save the computed data
  save_results(filename, buffer, count, sample_rate);
}

void compute_file(const char* folder, const char* filename, 
                  const char* sample_count, int sample_rate) {
  vector<Cplx> buffer;
  
  // Read data file
  read_file(filename, buffer);
  int count = buffer.size();
  
  // Is power of 2?
  if (count & (count-1)) {
    cerr << "Input data sample count have to be power of two." << endl;
  } else {
    compute(&buffer[0], count, sample_rate, filename);
  }
}

int main(int argc, char** argv) {
  srand (time(NULL));
  
  // Deal with program arguments
  if (argc < 2) {
    cerr << "Usage: " << argv[0] << " [input_folder]"; return 2;
  }
  
  compute_file("data/50Hz", "data/50Hz/128smp@44100.dat", "128", 44100);

  // Compute all files in folder
  // DIR* dirp = opendir(argv[1]);
  // struct dirent *epdf;
  // while ((epdf = readdir(dirp)) != NULL) {
  //   size_t len = strlen(epdf->d_name);
  //   if (strcmp(epdf->d_name,".") != 0 && strcmp(epdf->d_name,"..") != 0
  //       && strcmp(&epdf->d_name[len-3], "dat") == 0) {
  //     stringstream fname(epdf->d_name);
  //     string samples, sr;
  // 
  //     getline(fname, samples, '@');
  //     getline(fname, sr, '.');
  // 
  //     char* fold = new char[512];
  //     strcpy(fold, argv[1]);
  // 
  //     compute_file(argv[1], strcat(strcat(fold, "/"), epdf->d_name), samples.c_str(), atoi(sr.c_str()));
  //   }
  // }
  // closedir(dirp);
  
  return 0;
}