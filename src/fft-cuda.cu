#include <iostream>
#include <fstream>
#include <sstream>
#include <chrono>
#include <vector>
#include <cmath>
#include <dirent.h>
#include <cstring>

using namespace std;

#define N_REPEAT    3

// Complex numbers data type
typedef float2 Cplx;

// Complex numbers operations
static __device__ __host__ inline Cplx CplxAdd(Cplx a, Cplx b) {
  Cplx c; c.x = a.x + b.x; c.y = a.y + b.y; return c;
}

static __device__ __host__ inline Cplx CplxInv(Cplx a) {
  Cplx c; c.x = -a.x; c.y = -a.y; return c;
}

static __device__ __host__ inline Cplx CplxMul(Cplx a, Cplx b) {
  Cplx c; c.x = a.x * b.x - a.y + b.y; c.y = a.x * b.y + a.y * b.x; return c;
}

/**
 * Reorders array by bit-reversing the indexes.
 */
__global__ void bitrev_reorder(Cplx* __restrict__ r, Cplx* __restrict__ d, int s, size_t nthr) {
  int id = blockIdx.x * nthr + threadIdx.x;
  r[__brev(id) >> (32 - s)] = d[id];
}

/**
 * Inner part of FFT loop. Contains the procedure itself.
 */
__device__ void inplace_fft_inner(Cplx* __restrict__ r, int j, int k, int m, int n) {
  if (j + k + m / 2 < n && k < m / 2) { 
    Cplx t, u;
    
    t.x = __cosf((2.0 * M_PI * k) / (1.0 * m));
    t.y = -__sinf((2.0 * M_PI * k) / (1.0 * m));
    
    u = r[j + k];
    t = CplxMul(t, r[j + k + m / 2]);

    r[j + k] = CplxAdd(u, t);
    r[j + k + m / 2] = CplxAdd(u, CplxInv(t));
  }
}

/**
 * Middle part of FFT for small scope paralelism.
 */
__global__ void inplace_fft(Cplx* __restrict__ r, int j, int m, int n, size_t nthr) {
  int k = blockIdx.x * nthr + threadIdx.x;
  inplace_fft_inner(r, j, k, m, n);
}

/**
 * Outer part of FFT for large scope paralelism.
 */
__global__ void inplace_fft_outer(Cplx* __restrict__ r, int m, int n, size_t nthr) {
  int j = (blockIdx.x * nthr + threadIdx.x) * m;
  
  for (int k = 0; k < m / 2; k++) {
    inplace_fft_inner(r, j, k, m, n);
  }
}

/**
 * Runs in-place Iterative Fast Fourier Transformation.
 */
void fft(Cplx* __restrict__ d, size_t n, size_t threads, int balance) {
  size_t data_size = n * sizeof(Cplx);
  Cplx *r, *dn;
  
  // Copy data to GPU
  cudaMalloc((void**)&r, data_size);
  cudaMalloc((void**)&dn, data_size);
  cudaMemcpy(dn, d, data_size, cudaMemcpyHostToDevice);
  
  // Bit-reversal reordering
  int s = log2(n);
  bitrev_reorder<<<ceil(n / threads), threads>>>(r, dn, s, threads);
  
  // Synchronize
  cudaDeviceSynchronize();
  
  // Iterative FFT (with loop paralelism balancing)
  for (int i = 1; i <= s; i++) {
    int m = 1 << i;
    if (n/m > balance) {
      inplace_fft_outer<<<ceil((float)n / m / threads), threads>>>(r, m, n, threads);
    } else {
      for (int j = 0; j < n; j += m) {
        float repeats = m / 2;
        inplace_fft<<<ceil(repeats / threads), threads>>>(r, j, m, n, threads);
      }
    }
  }
  
  // Copy data from GPU & free the memory blocks
  Cplx* result;
  result = (Cplx*)malloc(data_size / 2);
  cudaMemcpy(result, r, data_size / 2, cudaMemcpyDeviceToHost);
  cudaFree(r);
  cudaFree(dn);
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

/**
 * Saves the result data to an output file.
 */
void save_results(const char* filename, Cplx* result, size_t count, int sample_rate) {
  char* outfilename = new char[512];
  
  // Compose the output filename
  strcpy(outfilename, filename);
  strcat(outfilename, ".out");
  
  // Create the file
  ofstream outfile;
  outfile.open (outfilename);
  outfile.precision(4);
  
  // Save the data
  outfile << "frequency, value" << endl;
  for (int i = 0; i < count / 2; i++) {
      outfile << i * ((float)sample_rate/count) << "," << result[i].x << endl;
  }
  
  outfile.close();
}

void compute_file(const char* filename, int sample_rate, size_t threads, int balance) {
  vector<Cplx> buffer;
  
  // Read the file
  read_file(filename, buffer);
  int count = buffer.size();
  
  // Display active computation
  cout << filename << "," << count << "," << sample_rate << "," << threads << 
          "," << balance;
  cout.flush();
  
  // Start the stopwatch
  auto start = chrono::high_resolution_clock::now();
  
  // Run FFT algorithm with loaded data
  fft(&buffer[0], count, threads, balance);
  
  // Log the elapsed time
  auto finish = chrono::high_resolution_clock::now();
  auto microseconds = chrono::duration_cast<std::chrono::microseconds>(finish-start);
  
  cout << "," << microseconds.count() << endl;
  
  // Save the computed data
  save_results(filename, &buffer[0], count, sample_rate);
}

int main(int argc, char** argv) {
  srand (time(NULL));
  
  // Deal with program arguments
  if (argc < 2) {
    cerr << "Usage: " << argv[0] << " [input_folder]"; return 2;
  }

  // Initialize CUDA
  cudaFree(0);
  
  // Print out the CSV header
  cout << "file,samples,sample_rate,threads,balance,elapsed_us" << endl;
  
  // Read the folder
  DIR* dirp = opendir(argv[1]);
  struct dirent *epdf;
  
  // Compute all files in the folder
  while ((epdf = readdir(dirp)) != NULL) {
    size_t len = strlen(epdf->d_name);
    
    // Pick only .dat files
    if (strcmp(epdf->d_name,".") != 0 && strcmp(epdf->d_name,"..") != 0
        && strcmp(&epdf->d_name[len-3], "dat") == 0) {
      stringstream fname(epdf->d_name);
      string samples, sr;
  
      // Read file properties
      getline(fname, samples, '@');
      getline(fname, sr, '.');
  
      char* fold = new char[512];
      strcpy(fold, argv[1]);
      int smp = atoi(samples.c_str());
  
      // Compute for all set parameters
      for (int th = 0; th <= 1024; th <<= 1) {
        if (th == 0) th = 1;
        for (int bal = 2; bal <= smp / 2; bal <<= 1)
        for (int r = 0; r < N_REPEAT; r++) {
          char fname[512];
          strcpy(fname, fold);
          strcat(strcat(fname, "/"), epdf->d_name);
          compute_file(fname, atoi(sr.c_str()), th, bal);
        }
      }
    }
  }
  closedir(dirp);
  
  return 0;
}