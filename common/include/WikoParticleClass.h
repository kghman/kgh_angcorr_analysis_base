#ifndef __wikoparticleclass_h
#define __wikoparticleclass_h

#include <complex>

class WikoParticleClass {

 public:

  WikoParticleClass(double iIi, double iGam=1.0); //presuming 0+ final state of heavy decay (fewer params)
  ~WikoParticleClass();

  double Get_IiPi();
  double Get_BR();
  int Get_Max_Rank();
  int Get_L();

  void Set_Ii(double iIi);
  void Set_Gam(double iGam);

  double Decay_Coefficient(int k);
  double Spherical_Harmonic(int k, int q, double theta, double phi); //angle inputs in radians
  double Spherical_Harmonic_Cosine(int k, int q, double cos_theta, double phi);
  std::complex<double> Complex_Spherical_Harmonic(int k, int q, double theta, double phi); //angle inputs in radians
  std::complex<double> Complex_Spherical_Harmonic_Cosine(int k, int q, double cos_theta, double phi);

 private:

  double Ii, If; //initial and final state spins
  int Pi, Pf; //initial and final state parity (+/-1)
  bool odd_L;
  int max_rank;
  double j,s; //(relative) momenta and spin
  int l; //same but orbital
  double Gam; //branching ratio (input)
  
  double *coeffs;

  bool Check_kq(int k, int q=0);
  double Round_Spin(double ispin);
  void Generate_Coeffs();

};

#endif
