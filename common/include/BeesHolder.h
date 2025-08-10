#ifndef __beesholder_h
#define __beesholder_h

#include <complex>
#include <fstream>

class BeesHolder {

 public:

  //inputs: path of fort.37 file (required), euler angles (radians, optional)
  BeesHolder(char ifile_path[],
	     double ieu_a=0,
	     double ieu_b=0,
	     double ieu_g=0);
  ~BeesHolder();

  double Euler_Alpha();
  double Euler_Beta();
  double Euler_Gamma();

  int Max_Rank();
  int Num_Angles();
  int Num_MSubs();

  double J();
  double Beam_E(); //MeV (lab)
  double Total_XSec(); //mb
  double DTheta(); //radians

  bool File_Loaded();

  double Angle(int angindex); //radians
  double XSec(int angindex); //mb/sr
  double MSubs(double m, int angindex);
  double MSubs_Frac(double m, int angindex);
  double MSubs_Integ(double m);
  double MSubs_Integ_Frac(double m);
  double DensMat_Real(double m1, double m2, int angindex);
  double DensMat_Imag(double m1, double m2, int angindex);
  double Bkq_Real(int k, int q, int angindex);
  double Bkq_Imag(int k, int q, int angindex);

  double* Angle();
  double* XSec();
  double* MSubs(double m);
  double* MSubs_Frac(double m);
  double* MSubs_Integ();
  double* MSubs_Integ_Frac();
  double* DensMat_Real(double m1, double m2);
  double* DensMat_Imag(double m1, double m2);
  double* Bkq_Real(int k, int q);
  double* Bkq_Imag(int k, int q);

  std::complex<double> Bkq_Complex(int k, int q, int angindex);

  double** MSubs();
  double** MSubs_Frac();

  void Set_Euler_Angles(double ieu_a, double ieu_b, double ieu_g); //expects radians
  bool Set_File_Path(char ifile_path[]);

  bool Generate_Tensors();
  bool Generate_Tensors(char ifile_path[]);

  void Toggle_Norm_to_B00();
  bool Norm_to_B00();

  int Factorial(int n);

  double RotMat(double j, double m, double n);
  
 private:

  const double THRESHOLD;

  double JT,JP,JE,JR; //target, projectile, ejectile, residual
  double beamE; //beam energy in lab
  double total_xsec; //total cross-section, integrated (mb)

  //Euler angles (radians)
  double
    euler_alpha,
    euler_beta,
    euler_gamma;

  int max_rank;
  int num_ranks;
  int num_msubs;
  int num_angs;

  char file_path[100];
  std::ifstream fort37_file;
  bool file_loaded;

  double dtheta;

  double *angle;
  double *xsec; //in units of mb/sr
  double **msubs; //in units of xsec (mb/sr)
  double **msubs_frac; //in relative units (fractional)
  double *msubs_integ; //same, but integrated over all angles
  double *msubs_integ_frac;

  double *return_arr; //return array for certain functions

  std::complex<double> ***densmat;
  std::complex<double> ***Bkq;

  bool norm_to_B00;

  void Initialize_Tensors();
  void Clear_Tensors();

  bool Check_MSubs_Input(double m1, double m2=0, int angindex=0);
  bool Check_kq_Input(int k, int q=0, int angindex=0);

};

#endif
