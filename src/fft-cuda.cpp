#include <complex>
#include <iostream>
#include <fstream>
#include <valarray>
#include <time.h>
#include <vector>
 
using namespace std;
 
typedef complex<double> Complex;
typedef valarray<Complex> CplxArray;

/**
 * In-place Iterative Fast Fourier Transformation.
 */
CplxArray fft(CplxArray& x) {
  #pragma acc kernels 
  {
    const size_t n = x.size();
    CplxArray r(n);
  	int s = log2(n);
  	
    // Bit-reversal reordering
  	for (int i = 0; i <= n-1; i++) {
  		int j = i, k = 0;
      
  		for (int l = 1; l <= s; l++) {
  			k = k * 2 + (j & 1);
        j >>= 1;
      }

  		r[k] = x[i];
  	}
    
  	for (int i = 1; i <= s; i++) {
  		int m = 1 << i;
      
  		for (int j = 0; j < n; j += m) for(int k = 0; k < m/2; k++) {
        Complex t = polar(1.0, -2 * M_PI * k / m);
  			Complex u = r[j + k];
  			t *= r[j + k + m / 2];

  			r[j + k] = u + t; 
  			r[j + k + m / 2] = u - t;
  		}
  	}
    
  	return r;
  }
}

/**
 * Converts array of real numbers to complex numbers array.
 */
CplxArray* real_to_complex(vector<double> real_array) {
  CplxArray* ary = new CplxArray(real_array.size());
  
  for (int i = 0; i < real_array.size(); i++) {
    (*ary)[i] = complex<double>(real_array[i], 0.0);
  }
  
  return ary;
}

/**
 * Reads numeric data from a file.
 */
void read_file(char* filename, vector<double>& out) {
  ifstream file;
  file.open(filename);
  
  if (file.is_open()) {
    while (!file.eof()) {
      double val;
      if (file >> val) {
        out.push_back(val);
      }
    }
  }
  
  file.close();
}

int main(int argc, char** argv) {
  // Deal with program arguments
  if (argc < 3) {
    cerr << "Usage: " << argv[0] << " [input_file] [output_file] [sample_rate]"; return 2;
  }
  
  vector<double> buffer;
  int sample_rate = atoi(argv[3]); 
  
  // Read data file
  read_file(argv[1], buffer);
  int count = buffer.size();
  CplxArray data = *real_to_complex(buffer);
  
  // Is power of 2?
  if (count & (count-1)) {
    cerr << "Input data sample count have to be power of two." << endl;
    return 2;
  }
  
  // Print out stats
  cout << "Data info:" << endl << "Sample count: " << count << endl << "Sample rate: " << sample_rate << endl;
  
  // Start the stopwatch
  clock_t t;
	t = clock();
  
  // Run FFT algorithm with loaded data
  CplxArray result = fft(data);
  
  // Print out the total elapsed time
  t = clock() - t;
  cout << "Program took " << t*1000.0/CLOCKS_PER_SEC << "ms to finish." << endl;
  
  // Save the computed data
  ofstream outfile;
  outfile.open (argv[2]);
  outfile.precision(4);
  outfile << "frequency, value" << endl;
  for (int i = 0; i < count / 2; i++) {
      outfile << i * ((double)sample_rate/count) << "," << abs(result[i]) << endl;
  }
  outfile.close();
  
  return 0;
}