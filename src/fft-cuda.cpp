#include <complex>
#include <iostream>
#include <fstream>
#include <valarray>
 
using namespace std;
 
typedef complex<double> Complex;
typedef valarray<Complex> CplxArray;

/**
 * Recursive, in-place Fast-Fourier Transformation.
 */
// https://rosettacode.org/wiki/Fast_Fourier_transform#C.2B.2B
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
CplxArray* real_to_complex(double* real_array, int length) {
  CplxArray* ary = new CplxArray(length);
  
  for (int i = 0; i < length; i++) {
    (*ary)[i] = complex<double>(real_array[i], 0.0);
  }
  
  return ary;
}

/**
 * Reads numeric data from a file.
 */
int read_file(char* filename, double* out_ary) {
  ifstream file;
  int pos = 0;
  file.open(filename);
  
  if (file.is_open()) {
    while (!file.eof()) {
      file >> out_ary[pos++];
    }
  }
  
  file.close();
  return pos;
}

int main(int argc, char** argv) {
  // Deal with program arguments
  if (argc < 3) {
    cerr << "Usage: " << argv[0] << " [input_file] [sample_rate]"; return 2;
  }
  int sample_rate = atoi(argv[2]);
  int count = read_file(argv[1], buffer); 
  
  // Read data file
  double buffer[1024];
  CplxArray data = *real_to_complex(buffer, count);
  
  // Run FFT algorithm with loaded data
  fft(data);
  
  // Print out the computed data
  cout.precision(4);
  for (int i = 0; i < count / 2; i++) {
      cout << i << "\t" << i * ((double)sample_rate/count) << "-" << (i+1) * ((double)sample_rate/count) << "Hz\t" << abs(data[i]) << endl;
  }
  
  return 0;
}