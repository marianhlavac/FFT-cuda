#include <iostream>
#include <fstream>
#include <sstream>
#include <chrono>
#include <vector>
#include <cmath>
#include <dirent.h>

using namespace std;

#define NUM_GANGS_START 2
#define NUM_GANGS_END 512

/**
 * In-place Iterative Fast Fourier Transformation with OpenACC support.
 */
template<int T>
float* fft(float* __restrict__ x, size_t n) {
  size_t n2 = n * 2;
  int s = log2(n);
  float* r = new float[n2];
  
  #pragma acc data copy(x[0:n]) create(r[0:n2]) copyout(r[0:n])
  {
    // Bit-reversal reordering
    #pragma acc parallel loop
    for (int i = 0; i <= n-1; i++) {
      int j = i, k = 0;
      
      #pragma acc loop seq
      for (int l = 1; l <= s; l++) {
        k = k * 2 + (j & 1);
        j >>= 1;
      }
      
      r[k*2]= x[i];
      r[k*2 + 1] = 0;
    }
    
    #pragma acc loop seq
  	for (int i = 1; i <= s; i++) {
  		int m = 1 << i;
    
      #pragma acc parallel loop num_gangs(T)
  		for (int j = 0; j < n; j += m) {
        
        #pragma acc loop independent
        for(int k = 0; k < m/2; k++) {
          if (j + k + m / 2 < n) { 
            float t[2], u[2], tr[2];
            
            t[0] = cos((2.0*M_PI*k)/(1.0*m));
            t[1] = -sin((2.0*M_PI*k)/(1.0*m));
            
            size_t ridx = (j + k) * 2;
      			u[0] = r[ridx];
            u[1] = r[ridx + 1];
            
            size_t ridx2 = ridx + m;
            tr[0] = t[0] * r[ridx2] - t[1] * r[ridx2 + 1];
      			tr[1] = t[0] * r[ridx2] + t[1] * r[ridx2];
      
            t[0] = tr[0];
            t[1] = tr[1];
      
            r[ridx] = u[0] + t[0]; 
      			r[ridx + 1] = u[1] + t[1]; 
      
            r[ridx2] = u[0] - t[0];
      			r[ridx2 + 1] = u[1] - t[1];
            
            // Whew. That was ugly. Eww.
          }
    		}
        
      }
      
  	}
  }

  return r;
}

/**
 * Reads numeric data from a file.
 */
void read_file(const char* filename, vector<float>& out) {
  ifstream file;
  file.open(filename);
  
  if (file.is_open()) {
    while (!file.eof()) {
      float val;
      if (file >> val) {
        out.push_back(val);
      }
    }
  } else {
    cerr << "Can't open file " << filename << " for reading." << endl;
  }
  
  file.close();
}

void compute(int num_gangs, float* buffer, size_t count, int sample_rate, const char* filename) {
  float time = 0;
  
  float* result = new float[1];
  
  for (int i = 0; i < 4; i++) {
    // Start the stopwatch
    auto start = chrono::high_resolution_clock::now();
    
    // Run FFT algorithm with loaded data
    switch (num_gangs) {
      case 2: result = fft<2>(buffer, count); break;
      case 4: result = fft<4>(buffer, count); break;
      case 8: result = fft<8>(buffer, count); break;
      case 16: result = fft<16>(buffer, count); break;
      case 32: result = fft<32>(buffer, count); break;
      case 64: result = fft<64>(buffer, count); break;
      case 128: result = fft<128>(buffer, count); break;
      case 256: result = fft<256>(buffer, count); break;
      case 512: result = fft<512>(buffer, count); break;
    }
    // Log the elapsed time
    auto finish = chrono::high_resolution_clock::now();
    auto microseconds = chrono::duration_cast<std::chrono::microseconds>(finish-start);
    
    time += microseconds.count();
  }
  
  cout << count << "," << sample_rate << "," << num_gangs << "," << time/4 << endl;
  
  // Save the computed data
  char* outfilename = new char[512];
  strcpy(outfilename, filename);
  strcat(outfilename, ".out");
  ofstream outfile;
  outfile.open (outfilename);
  outfile.precision(4);
  outfile << "frequency, value" << endl;
  for (int i = 0; i < count / 2; i++) {
      outfile << i * ((float)sample_rate/count) << "," << result[i] << endl;
  }
  outfile.close();
}

void compute_file(const char* folder, const char* filename, 
                  const char* sample_count, int sample_rate) {
  vector<float> buffer;
  
  // Read data file
  read_file(filename, buffer);
  int count = buffer.size();
  
  // Is power of 2?
  if (count & (count-1)) {
    cerr << "Input data sample count have to be power of two." << endl;
  } else {
    // Go through initialization
    fft<1>(&buffer[0], count);
    
    // Run compute for various num_gangs
    for (int ng = NUM_GANGS_START; ng <= NUM_GANGS_END; ng <<= 1) {
      cout << folder << ",";
      compute(ng, &buffer[0], count, sample_rate, filename);
    }
  }
}

int main(int argc, char** argv) {
  srand (time(NULL));
  
  // Deal with program arguments
  if (argc < 2) {
    cerr << "Usage: " << argv[0] << " [input_folder]"; return 2;
  }

  // Compute all files in folder
  DIR* dirp = opendir(argv[1]);
  struct dirent *epdf;
  while ((epdf = readdir(dirp)) != NULL) {
    size_t len = strlen(epdf->d_name);
    if (strcmp(epdf->d_name,".") != 0 && strcmp(epdf->d_name,"..") != 0
        && strcmp(&epdf->d_name[len-3], "dat") == 0) {
      stringstream fname(epdf->d_name);
      string samples, sr;
      
      getline(fname, samples, '@');
      getline(fname, sr, '.');
      
      char* fold = new char[512];
      strcpy(fold, argv[1]);
      
      compute_file(argv[1], strcat(strcat(fold, "/"), epdf->d_name), samples.c_str(), atoi(sr.c_str()));
    }
  }
    
  closedir(dirp);
  
  return 0;
}