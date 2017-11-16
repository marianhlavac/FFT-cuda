#include <iostream>
#include <cmath>
#include <fstream>

using namespace std;

/**
 * Generates sinusoid oscillation.
 */
double gen_freq(double time, double freq, double phase, double amplitude) {
  return sin(2.0f * (time * M_PI + phase) * freq ) * amplitude;
}

/**
 * Generator function used for generating the data in this program.
 * Is intended to be changed here, to generate various patterns.
 */
double gen_func(double time) {
  // -!- Change the generator function here:
  return gen_freq(time, 12000.0f, 0, 0.05f) + gen_freq(time, 55.0f, 0, 0.375f) 
    + gen_freq(time, 56.0f, 0, 0.375f) + gen_freq(time, 13420.0f, 0, 0.012f);
}

int main(int argc, char** argv) {
  // Deal with program arguments
  if (argc < 4) {
    cerr << "Usage: " << argv[0] << " [output_file] [sample_rate] [samples]" << endl; 
    return 2;
  }
  int sample_rate = atoi(argv[2]);
  int samples_count = atoi(argv[3]);
  
  // Open file for writing
  ofstream out_file;
  out_file.open(argv[1]);
  
  // Calculate output parameters
  double sample_length = 1.0f / sample_rate;
  double length = sample_length * samples_count;
  
  // Generate the data, saving into the file
  if (out_file.is_open()) {
    for (double t = 0; t < length; t += sample_length) {
        out_file << gen_func(t) << endl;
    }
  } else {
    cerr << "Error: File " << argv[1] << " could not be open." << endl;
    return 1;
  }
  
  cout << "Done." << endl;
  out_file.close();
  return 0;
}