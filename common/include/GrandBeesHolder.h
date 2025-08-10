#ifndef __grandbeesholder_h
#define __grandbeesholder_h

#include <complex>
#include <fstream>

#include "WikoParticleClass.h"
#include "BeesHolder.h"

/*

  It is PRESUMED that the reaction which populated each residual state (target, beam species,
  energy, etc.) is the same for them all

  It is also presumed that each state is unbound -- necessary since this formalism
  is only for states which physically overlap in their widths, and therefore interfere
  on this level

  The idea (based on R. Huby and Q. K. K. Liu, NPA122 1968) is that the reaction amplitudes
  for each state are generated assuming they are bound, then the resonance energies and widths
  are inputted separately as "weights" later

  Assuming decay down to a 0+ final state

*/

class GrandBeesHolder {

 public:

  //inputs: number of states/fort.37 files (required)
  //        array of fort.37 file names (required)
  //        array of resonance energies Er (for energy averaging)
  //        array of resonance widths Gamma (same)
  //        array of branching ratios per decay component
  //        charge numbers and masses (MeV) of decay products for Coulomb interference
  //        euler angles (in radians; optional)
  GrandBeesHolder(int inum_states,
		  char **ifile_paths,
                  double *iEr,
                  double *iGamma,
		  double *iBR,
		  int *iparity,
		  int iZ1,
		  int iZ2,
		  double iM1,
		  double iM2,
		  double *iE12,
		  double ieu_a=0,
		  double ieu_b=0,
		  double ieu_g=0);
  ~GrandBeesHolder();

  double Euler_Alpha();
  double Euler_Beta();
  double Euler_Gamma();

  int Num_States();
  int Max_Rank();
  int Num_Angles();
  int Num_Tot_MSubs();

  double Get_Er(int state);
  double Get_Gamma(int state);
  double Get_BR(int state);
  double J(int state);
  int    Parity(int state);
  double Jd(int state);
  double Ld(int state);
  double Beam_E(); //MeV (lab)
  double Total_XSec(); //mb
  double DTheta(); //radians

  std::complex<double> Get_Eta(int state1, int state2); //resonance energy averaging factor
  
  bool Files_Loaded();

  double Angle(int angindex); //radians
  double XSec(int angindex); //mb/sr
  double DensMat_Real(int state1, int state2, double m1, double m2, int angindex);
  double DensMat_Imag(int state1, int state2, double m1, double m2, int angindex);
  double Bkq_Real(int state1, int state2, int k, int q, int angindex);
  double Bkq_Imag(int state1, int state2, int k, int q, int angindex);
  double akq_Real(int k, int q, int angindex);
  double akq_Imag(int k, int q, int angindex);
  double Ak(int state, int k);

  double* Angle();
  double* XSec();
  double* DensMat_Real(int state1, int state2, double m1, double m2);
  double* DensMat_Imag(int state1, int state2, double m1, double m2);
  double* Bkq_Real(int state1, int state2, int k, int q);
  double* Bkq_Imag(int state1, int state2, int k, int q);
  double* akq_Real(int k, int q);
  double* akq_Imag(int k, int q);

  std::complex<double> Bkq_Complex(int state1, int state2, int k, int q, int angindex);
  std::complex<double> akq_Complex(int k, int q, int angindex);

  //borrowed from the wiko class
  double Spherical_Harmonic(int k, int q, double theta, double phi); //angle inputs in radians
  double Spherical_Harmonic_Cosine(int k, int q, double cos_theta, double phi);
  std::complex<double> Complex_Spherical_Harmonic(int k, int q, double theta, double phi); //angle inputs in radians
  std::complex<double> Complex_Spherical_Harmonic_Cosine(int k, int q, double cos_theta, double phi);
  
  void Set_Euler_Angles(double ieu_a, double ieu_b, double ieu_g); //expects radians
  bool Set_File_Paths(char **ifile_paths);

  bool Generate_Tensors();
  bool Generate_Tensors(char **ifile_paths);

  int Factorial(int n);

  double RotMat(double j, double m, double n);

  std::complex<double> Coulomb_phase_shift(int l, double Ruthparam);
  
  //have individual BeesHolders and decay coefficients at public access
  WikoParticleClass **indiv_wiko;
  BeesHolder        **indiv_bees;
  
 private:

  const double THRESHOLD;

  double JT,JP,JE; //target, projectile, ejectile
  double beamE; //beam energy in lab
  double total_xsec; //total cross-section, integrated (mb)

  //Euler angles (radians)
  double
    euler_alpha,
    euler_beta,
    euler_gamma;

  int num_states;

  //each property of the residual states now depends on the number of states
  double *Er; //resonance energy
  double *Gamma; //resonance widths (parameters)
  double *BR; //branching ratios (decay probabilities)
  std::complex<double> **eta; //weight factor due to energy integral over finite width
  std::complex<double> *Cphase; //Coulomb phases for fixed LP (one per state)
  double *JR, *JD; //residual spin and decay spin
  int    *PAR, *LD; //residual state parity and decay orbital ang mom
  double sp; //spin of proton
  double Jf; //final state spin of heavy decay residual (assumed 0+ for now)
  int *num_msubs;

  //kinematic specifics for the Coulomb phase shift calculation
  double *Ruthparam;
  int Z1, Z2; //charges of two decay particles (should be same for all states)
  double M1, M2, MU; //masses in MeV as well as reduced mass to be calculated
  double *E12, *k; //relative energy of two decay particles (input), and corresponding wave number
  
  bool unbound;

  int max_rank;
  int num_ranks;
  int num_tot_msubs;
  int num_angs;

  char **file_paths;
  std::ifstream *fort37_files;
  bool files_loaded;

  double *angle;
  double *xsec; //in units of mb/sr

  double *return_arr; //return array for certain functions

  std::complex<double> *****densmat;
  std::complex<double> *****Bkq;
  std::complex<double> ***akq;
  
  std::complex<double> eye;

  void Initialize_Tensors();
  void Clear_Tensors();

  bool Check_MSubs_Input(int state1, int state2, double m1, double m2, int angindex=0);
  bool Check_ijkq_Input(int state1, int state2, int k, int q=0, int angindex=0);
  
};

#endif
