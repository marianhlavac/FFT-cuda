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
 * Recursive, in-place Fast-Fourier Transformation.
 * https://rosettacode.org/wiki/Fast_Fourier_transform#C.2B.2B
 */
void fft(CplxArray& x) {
    const size_t N = x.size();
    if (N <= 1) return;
 
    CplxArray even = x[slice(0, N/2, 2)];
    CplxArray odd = x[slice(1, N/2, 2)];
 
    fft(even); fft(odd);
 
    for (size_t k = 0; k < N/2; ++k) {
        Complex t = polar(1.0, -2 * M_PI * k / N) * odd[k];
        x[k] = even[k] + t;
        x[k+N/2] = even[k] - t;
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
      file >> val;
      out.push_back(val);
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
  
  // Start the stopwatch
  clock_t t;
	t = clock();
  
  // Run FFT algorithm with loaded data
  fft(data);
  
  // Print out the total elapsed time
  t = clock() - t;
  cout << "Program took " << t*1000.0/CLOCKS_PER_SEC << "ms to finish." << endl;
  
  // Save the computed data
  ofstream outfile;
  outfile.open (argv[2]);
  outfile.precision(4);
  outfile << "frequency, value" << endl;
  for (int i = 0; i < count / 2; i++) {
      outfile << i * ((double)sample_rate/count) << "," << abs(data[i]) << endl;
  }
  outfile.close();
  
  return 0;
}