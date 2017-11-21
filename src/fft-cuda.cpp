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

double** twodim_double_ary(size_t n) {
  double** ary = new double*[n];
  for (int i = 0; i < n; i++) {
      ary[i] = new double[2];
  }
  return ary;
}

/**
 * In-place Iterative Fast Fourier Transformation.
 */
double** fft(double** x, size_t n) {
  double** r = twodim_double_ary(n);

  #pragma acc data copy(r[0:n][0:2]) 
  {
  	int s = log2(n);
  	
    // Bit-reversal reordering
    #pragma acc parallel loop
  	for (int i = 0; i <= n-1; i++) {
  		int j = i, k = 0;
      
      #pragma acc loop
  		for (int l = 1; l <= s; l++) {
  			k = k * 2 + (j & 1);
        j >>= 1;
      }
      
  		r[k][0] = x[i][0];
  		r[k][1] = x[i][1];
  	}
    
    #pragma acc parallel loop
  	for (int i = 1; i <= s; i++) {
  		int m = 1 << i;
      
      #pragma acc loop
  		for (int j = 0; j < n; j += m) {
        #pragma acc loop
        for(int k = 0; k < m/2; k++) {
          // Is PGP dumb, or what? This is ugly. TODO: fix
          double t[2], u[2], tr[2];
          t[0] = cos((2.0*M_PI*k)/(1.0*m));
          t[1] = -sin((2.0*M_PI*k)/(1.0*m));
    			u[0] = r[j + k][0];
          u[1] = r[j + k][1];
          
          double* rz = r[j + k + m / 2];
          
          tr[0] = t[0] * rz[0] - t[1] * rz[1];
    			tr[1] = t[0] * rz[1] + t[1] * rz[0];
          
          t[0] = tr[0];
          t[1] = tr[1];

          r[j + k][0] = u[0] + t[0]; 
    			r[j + k][1] = u[1] + t[1]; 

          rz[0] = u[0] - t[0];
    			rz[1] = u[1] - t[1];
    		}
      }
  	}
  }

  return r;
}

/**
 * Converts array of real numbers to complex numbers array.
 */
double** real_to_complex(vector<double> real_array) {
  double** ary = new double*[real_array.size()];
  for (int i = 0; i < real_array.size(); i++) {
      ary[i] = new double[2];
      ary[i][0] = real_array[i];
      ary[i][1] = 0.0;
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
  double** data = real_to_complex(buffer);
  
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
  double** result = fft(data, count);
  
  // Log the elapsed time
  #ifdef STAR
    const double seconds = omp_get_wtime() - start;
  #else
  	t = clock() - t;
    const double seconds = (double)t / CLOCKS_PER_SEC;
  #endif
  cout << "Program took " << seconds << " secs to finish." << endl;
  
  // Save the computed data
  ofstream outfile;
  outfile.open (argv[2]);
  outfile.precision(4);
  outfile << "frequency, value" << endl;
  for (int i = 0; i < count / 2; i++) {
      outfile << i * ((double)sample_rate/count) << "," << abs(result[i][0]) << endl;
  }
  outfile.close();
  
  return 0;
}