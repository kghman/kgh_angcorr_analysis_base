#include "gensim_randgenerator.h"

GenSim_RandGenerator::GenSim_RandGenerator() {srand(time(0));}

GenSim_RandGenerator::~GenSim_RandGenerator() {}

double GenSim_RandGenerator::Uniform_Dist(double min, double max) {

  if (min > max) {
    double temp = min;
    min = max;
    max = temp;
  }

  double interval = max - min;

  double prob;
  if (interval != 0)
    prob = 1/interval;
  else prob = 0;

  double randnum = (rand()%1000)/1000.;

  return min + interval*randnum;

}

double GenSim_RandGenerator::Normal_Dist(double mean, double sigma) {

  if (sigma == 0) sigma = 0.001;

  double range = 10*sigma;

  double dist[10000];
  double x = mean - range/2;
  for (int i=0; i<10000; i++) {
    dist[i] = exp(-1*(pow((x-mean),2)/(2*pow(sigma,2))));
    dist[i] /= sqrt(2*M_PI*sigma*sigma);
    x+=range/10000.;
  }

  double randnum = GenSim_RandGenerator::Uniform_Dist(-range/2,range/2);

  double currentval = -range/2;
  double returnval = mean;
  for (int i=0; i<10000; i++) {
    if (currentval >= randnum) {
      returnval = currentval;
      break;
    }
    currentval += dist[i];
  }

  return returnval;

}
