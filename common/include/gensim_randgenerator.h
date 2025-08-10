#ifndef __gensim_randgenerator_h
#define __gensim_randgenerator_h

#include <cmath>
#include <ctime>
#include <cstdlib>

using namespace std;

class GenSim_RandGenerator {

 public:

  GenSim_RandGenerator();
  ~GenSim_RandGenerator();

  double Uniform_Dist(double min, double max);
  double Normal_Dist(double mean, double sigma);

};

#endif
