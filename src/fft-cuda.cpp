#include <complex>
#include <iostream>
#include <fstream>
#include <valarray>
#include <chrono>
#include <vector>
 
using namespace std;

/**
 * In-place Iterative Fast Fourier Transformation.
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
          }
    		}
        
      }
      
  	}
  }

  return r;
}

/**
 * Converts array of real numbers to complex numbers array.
 */
// double** real_to_complex(vector<double> real_array) {
//   double** ary = new double*[real_array.size()];
//   for (int i = 0; i < real_array.size(); i++) {
//       ary[i] = new double[2];
//       ary[i][0] = real_array[i];
//       ary[i][1] = 0.0;
//   }
// 
//   return ary;
// }

/**
 * Reads numeric data from a file.
 */
void read_file(char* filename, vector<float>& out) {
  ifstream file;
  file.open(filename);
  
  if (file.is_open()) {
    while (!file.eof()) {
      float val;
      if (file >> val) {
        out.push_back(val);
      }
    }
  }
  
  file.close();
}

int main(int argc, char** argv) {
  srand (time(NULL));
  
  // Deal with program arguments
  if (argc < 3) {
    cerr << "Usage: " << argv[0] << " [input_file] [output_file] [sample_rate]"; return 2;
  }
  
  vector<float> buffer;
  int sample_rate = atoi(argv[3]); 
  
  // Read data file
  read_file(argv[1], buffer);
  int count = buffer.size();
  
  // Is power of 2?
  if (count & (count-1)) {
    cerr << "Input data sample count have to be power of two." << endl;
    return 2;
  }
  
  // Print out stats
  cout << "Data info:" << endl << "Sample count: " << count << endl << "Sample rate: " << sample_rate << endl;
  
  // Go through initialization
  fft<1>(&buffer[0], count);
  
  float time = 0;
  
  for (int i = 0; i < 4; i++) {
    // Start the stopwatch
    auto start = chrono::high_resolution_clock::now();
    
    // Run FFT algorithm with loaded data
    float* result = fft<32>(&buffer[0], count);
    
    // Log the elapsed time
    auto finish = chrono::high_resolution_clock::now();
    auto microseconds = chrono::duration_cast<std::chrono::microseconds>(finish-start);
    
    time += microseconds.count();
  }
  
  cout << "For num_gangs(32) program took " << time/4 << " us to finish." << endl;
  
  // // Save the computed data
  // ofstream outfile;
  // outfile.open (argv[2]);
  // outfile.precision(4);
  // outfile << "frequency, value" << endl;
  // for (int i = 0; i < count / 2; i++) {
  //     outfile << i * ((float)sample_rate/count) << "," << result[i] << endl;
  // }
  // outfile.close();
  
  return 0;
}