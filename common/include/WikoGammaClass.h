#ifndef __wikogammaclass_h
#define __wikogammaclass_h

#include <complex>

class WikoGammaClass {

 public:

  WikoGammaClass(double iIi, double iIf, double idelta21=0, double idelta31=0);
  ~WikoGammaClass();

  int Get_Max_Rank();
  int L();
  int LP();
  int LPP();
  int Pi();
  int Pf();
  double Ii();
  double If();
  double Delta21();
  double Delta31();

  void Set_Ii(double iIi);
  void Set_If(double iIf);
  void Set_Spins(double iIi, double iIf);
  void Set_Delta21(double idelta21);
  void Set_Delta31(double idelta31);
  void Set_Deltas(double idelta21, double idelta31);

  double Decay_Coefficient(int k);
  double Spherical_Harmonic(int k, int q, double theta, double phi); //angle inputs in radians
  double Spherical_Harmonic_Cosine(int k, int q, double cos_theta, double phi);
  std::complex<double> Complex_Spherical_Harmonic(int k, int q, double theta, double phi); //angle inputs in radians
  std::complex<double> Complex_Spherical_Harmonic_Cosine(int k, int q, double cos_theta, double phi);

 private:

  double Ii_v, If_v; //initial and final state spins
  int Pi_v, Pf_v; //initial and final state parity (+/-1)
  bool odd_L;
  int max_rank;
  int L_v, LP_v, LPP_v; //L = minimum allowed decay (1 = dipole, 2 = quadrupole, etc.); LP = L+1, etc.
  double delta21, delta31; //mixing ratios for LP/L and LPP/L respectively

  double *coeffs;

  bool Check_kq(int k, int q=0);
  double Round_Spin(double ispin);
  void Set_Momenta();
  void Set_Max_Rank();
  void Generate_Coeffs();

};

#endif
