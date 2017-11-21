#include <complex>
#include <iostream>
#include <fstream>
#include <valarray>
#include <ctime>
#include <vector>

#ifdef STAR
  #include <omp.h> 
#endif
 
using namespace std;
 
typedef complex<double> Complex;
typedef valarray<Complex> CplxArray;

/**
 * In-place Iterative Fast Fourier Transformation.
 */
CplxArray fft(CplxArray& x) {
  const size_t n = x.size();
  CplxArray r(n);

  #pragma acc parallel 
  #pragma acc data copy(r[0:n]) 
  {
  	int s = log2(n);
  	
    // Bit-reversal reordering
    #pragma acc loop
  	for (int i = 0; i <= n-1; i++) {
  		int j = i, k = 0;
      
      #pragma acc loop
  		for (int l = 1; l <= s; l++) {
  			k = k * 2 + (j & 1);
        j >>= 1;
      }

  		r[k] = x[i];
  	}
    
    #pragma acc loop
  	for (int i = 1; i <= s; i++) {
  		int m = 1 << i;
      
      #pragma acc loop
  		for (int j = 0; j < n; j += m) {
        #pragma acc loop
        for(int k = 0; k < m/2; k++) {
          Complex t = polar(1.0, -2 * M_PI * k / m);
    			Complex u = r[j + k];
    			t *= r[j + k + m / 2];

    			r[j + k] = u + t; 
    			r[j + k + m / 2] = u - t;
    		}
      }
  	}
  }
    
  return r;
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
  srand (time(NULL));
  
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
  #ifdef STAR
    double start = omp_get_wtime();
  #else
    clock_t t;
  	t = clock();
  #endif
  
  
  // Run FFT algorithm with loaded data
  CplxArray result = fft(data);
  
  // Log the elapsed time
  ofstream logfile;
  logfile.open("log.txt");
  time_t now = time(0);
  logfile << ctime(&now) << argv[1] << endl;
  #ifdef STAR
    const double seconds = omp_get_wtime() - start;
  #else
  	t = clock() - t;
    const double seconds = (double)t / CLOCKS_PER_SEC;
  #endif
  
  logfile << "Program took " << seconds << " secs to finish." << endl;
  cout << "Program took " << seconds << " secs to finish." << endl;
  
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