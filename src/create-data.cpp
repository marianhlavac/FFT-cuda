#include <iostream>
#include <cmath>
#include <fstream>
#include <sstream>

#define GENERATOR (gen_freq(time, 500.0f, 0, 0.1f) + gen_freq(time, 505.0f, 0.5f, 0.1f) + gen_freq(time, 12000.0f, 0.5f, 0.1f))
#define GENERATOR_NAME "500Hz+505Hz+12000Hz"
#define SAMPLE_COUNT_START 128
#define SAMPLE_COUNT_END 65536

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
  return GENERATOR;
}

void generate_file(int sample_rate, int sample_count) {
  // Create a name
  stringstream ss;
  ss << "data/" << GENERATOR_NAME << "/" << sample_count << 
        "smp@" << sample_rate << ".dat";
  
  // Open file for writing
  ofstream out_file;
  out_file.open(ss.str());
  
  // Calculate output parameters
  double sample_length = 1.0f / sample_rate;
  double length = sample_length * sample_count;
  
  // Generate the data, saving into the file
  if (out_file.is_open()) {
    for (double t = 0; t < length; t += sample_length) {
        out_file << gen_func(t) << endl;
    }
    cout << "Done." << endl;
  } else {
    cerr << "Error: File " << ss.str() << " could not be open." << endl;
  }
  
  out_file.close();
}

int main(int argc, char** argv) {
  // Deal with program arguments
  if (argc < 2) {
    cerr << "Usage: " << argv[0] << " [sample_rate]" << endl; 
    return 2;
  }
  int sample_rate = atoi(argv[1]);
  
  for (int sc = SAMPLE_COUNT_START; sc <= SAMPLE_COUNT_END; sc <<= 1) {
    generate_file(sample_rate, sc);
  }
  
  return 0;
}