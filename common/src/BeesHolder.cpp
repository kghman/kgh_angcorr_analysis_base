#include <iostream>
#include <fstream>
#include <cstring>
#include <cstdlib>

#include "wiko_wiko5.hpp"
#include "BeesHolder.h"

using namespace std;

BeesHolder::BeesHolder(char ifile_path[],
		       double ieu_a,
		       double ieu_b,
		       double ieu_g) : THRESHOLD(1E-10) {

  JT=0;
  JP=0;
  JE=0;
  JR=0;
  beamE=0;
  max_rank=0;
  num_ranks=1;
  num_msubs=1;
  num_angs=1;
  total_xsec=0;
  dtheta=0;
  angle = new double[num_angs];
  xsec  = new double[num_angs];
  return_arr = new double[num_angs];
  for (int ang=0; ang<num_angs; ang++) {
    angle[ang]=0;
    xsec[ang]=0;
    return_arr=0;
  }
  msubs = new double*[num_msubs];
  msubs_frac = new double*[num_msubs];
  msubs_integ = new double[num_msubs];
  msubs_integ_frac = new double[num_msubs];
  for (int m=0; m<num_msubs; m++) {
    msubs[m] = new double[num_angs];
    msubs_frac[m] = new double[num_angs];
    msubs_integ[m] = 0;
    msubs_integ_frac[m] = 0;
    for (int ang=0; ang<num_angs; ang++) {
      msubs[m][ang]=0;
      msubs_frac[m][ang]=0;
    }
  }
  Initialize_Tensors();

  norm_to_B00 = true;

  euler_alpha = ieu_a;
  euler_beta  = ieu_b;
  euler_gamma = ieu_g;

  strcpy(file_path,ifile_path);
  fort37_file.open(file_path);
  if (!fort37_file) file_loaded = false;
  else              file_loaded = true;

}

BeesHolder::~BeesHolder() {
  delete [] angle;
  delete [] xsec;
  delete [] return_arr;
  for (int m=0; m<num_msubs; m++) {
    delete [] msubs[m];
    delete [] msubs_frac[m];
  }
  delete [] msubs;
  delete [] msubs_frac;
  delete [] msubs_integ;
  delete [] msubs_integ_frac;

  Clear_Tensors();
  fort37_file.close();
}

double BeesHolder::Euler_Alpha() {return euler_alpha;}
double BeesHolder::Euler_Beta()  {return euler_beta;}
double BeesHolder::Euler_Gamma() {return euler_gamma;}

int BeesHolder::Max_Rank()   {return max_rank;}
int BeesHolder::Num_Angles() {return num_angs;}
int BeesHolder::Num_MSubs()  {return num_msubs;}

double BeesHolder::J()          {return JR;}
double BeesHolder::Beam_E()     {return beamE;}
double BeesHolder::Total_XSec() {return total_xsec;}
double BeesHolder::DTheta()     {return dtheta;}

bool BeesHolder::File_Loaded() {return file_loaded;}

double BeesHolder::Angle(int angindex)
{return Check_kq_Input(0,0,angindex) ? angle[angindex] : 0.0;}
double BeesHolder::XSec(int angindex)
{return Check_kq_Input(0,0,angindex) ? xsec[angindex] : 0.0;}
double BeesHolder::MSubs(double m, int angindex)
{return Check_MSubs_Input(m,m,angindex) ? msubs[(int)(JR+m)][angindex] : 0.0;}
double BeesHolder::MSubs_Frac(double m, int angindex)
{return Check_MSubs_Input(m,m,angindex) ? msubs_frac[(int)(JR+m)][angindex] : 0.0;}
double BeesHolder::MSubs_Integ(double m)
{return Check_MSubs_Input(m,m) ? msubs_integ[(int)(JR+m)] : 0.0;}
double BeesHolder::MSubs_Integ_Frac(double m)
{return Check_MSubs_Input(m,m) ? msubs_integ_frac[(int)(JR+m)] : 0.0;}
double BeesHolder::DensMat_Real(double m1, double m2, int angindex)
{return Check_MSubs_Input(m1,m2,angindex) ? real(densmat[(int)(JR+m1)][(int)(JR+m2)][angindex]) : 0.0;}
double BeesHolder::DensMat_Imag(double m1, double m2, int angindex)
{return Check_MSubs_Input(m1,m2,angindex) ? imag(densmat[(int)(JR+m1)][(int)(JR+m2)][angindex]) : 0.0;}
double BeesHolder::Bkq_Real(int k, int q, int angindex)
{return Check_kq_Input(k,q,angindex) ? real(Bkq[k][k+q][angindex]) : 0.0;}
double BeesHolder::Bkq_Imag(int k, int q, int angindex)
{return Check_kq_Input(k,q,angindex) ? imag(Bkq[k][k+q][angindex]) : 0.0;}

double* BeesHolder::Angle() {return angle;}
double* BeesHolder::XSec()  {return xsec;}
double* BeesHolder::MSubs(double m) {return Check_MSubs_Input(m,m) ? msubs[(int)(JR+m)] : 0;}
double* BeesHolder::MSubs_Frac(double m) {return Check_MSubs_Input(m,m) ? msubs_frac[(int)(JR+m)] : 0;}
double* BeesHolder::MSubs_Integ()      {return msubs_integ;}
double* BeesHolder::MSubs_Integ_Frac() {return msubs_integ_frac;}

double* BeesHolder::DensMat_Real(double m1, double m2) {
  if (!(Check_MSubs_Input(m1,m2))) return 0;
  for (int ang=0; ang<num_angs; ang++)
    return_arr[ang] = real(densmat[(int)(JR+m1)][(int)(JR+m2)][ang]);
  return return_arr;
}

double* BeesHolder::DensMat_Imag(double m1, double m2) {
  if (!(Check_MSubs_Input(m1,m2))) return 0;
  for (int ang=0; ang<num_angs; ang++)
    return_arr[ang] = imag(densmat[(int)(JR+m1)][(int)(JR+m2)][ang]);
  return return_arr;
}

double* BeesHolder::Bkq_Real(int k, int q) {
  if (!(Check_kq_Input(k,q))) return 0;
  for (int ang=0; ang<num_angs; ang++)
    return_arr[ang] = real(Bkq[k][k+q][ang]);
  return return_arr;
}

double* BeesHolder::Bkq_Imag(int k, int q) {
  if (!(Check_kq_Input(k,q))) return 0;
  for (int ang=0; ang<num_angs; ang++)
    return_arr[ang] = imag(Bkq[k][k+q][ang]);
  return return_arr;
}

complex<double> BeesHolder::Bkq_Complex(int k, int q, int angindex) {
  return (Check_kq_Input(k,q,angindex) ? Bkq[k][k+q][angindex] : (complex<double>)(0.,0.));
}

double** BeesHolder::MSubs()      {return msubs;}
double** BeesHolder::MSubs_Frac() {return msubs_frac;}

void BeesHolder::Set_Euler_Angles(double ieu_a, double ieu_b, double ieu_g) {

  euler_alpha = ieu_a;
  euler_beta  = ieu_b;
  euler_gamma = ieu_g;

}

bool BeesHolder::Set_File_Path(char ifile_path[]) {
  strcpy(file_path,ifile_path);
  fort37_file.close();
  fort37_file.open(file_path);
  if (!fort37_file) file_loaded = false;
  else              file_loaded = true;
  return file_loaded;
}

bool BeesHolder::Generate_Tensors() {

  delete [] angle;
  delete [] xsec;
  delete [] return_arr;
  for (int m=0; m<num_msubs; m++) {
    delete [] msubs[m];
    delete [] msubs_frac[m];
  }
  delete [] msubs;
  delete [] msubs_frac;
  delete [] msubs_integ;
  delete [] msubs_integ_frac;

  Clear_Tensors();

  double throwaway;
  //             ja    jb    jap   jbp    (Tostevin notation)
  fort37_file >> JP >> JT >> JE >> JR >> num_angs >> throwaway >> beamE;
  int num_P_msubs = (int)(2*JP+1);
  int num_T_msubs = (int)(2*JT+1);
  int num_E_msubs = (int)(2*JE+1);
  num_msubs = (int)(2*JR+1);
  
  max_rank = (int)(2*JR);
  num_ranks = max_rank+1;

  total_xsec = 0;
  dtheta = 0;

  angle = new double[num_angs];
  xsec  = new double[num_angs];
  return_arr = new double[num_angs];
  for (int ang=0; ang<num_angs; ang++) {
    angle[ang]=0;
    xsec[ang]=0;
    return_arr[ang]=0;
  }
  msubs = new double*[num_msubs];
  msubs_frac = new double*[num_msubs];
  msubs_integ = new double[num_msubs];
  msubs_integ_frac = new double[num_msubs];
  for (int m=0; m<num_msubs; m++) {
    msubs[m] = new double[num_angs];
    msubs_frac[m] = new double[num_angs];
    msubs_integ[m] = 0;
    msubs_integ_frac[m] = 0;
    for (int ang=0; ang<num_angs; ang++) {
      msubs[m][ang]=0;
      msubs_frac[m][ang]=0;
    }
  }

  Initialize_Tensors();

  double real_part=0, imag_part=0;
  complex<double> f[num_P_msubs][num_T_msubs][num_E_msubs][num_msubs];
  complex<double> f_rot[num_P_msubs][num_T_msubs][num_E_msubs][num_msubs];

  double av_fac = 1./((2*JT+1)*(2*JP+1)); //average over initial states
  double fmsq_to_mb = 10.;

  for (int ang=0; ang<num_angs; ang++) {

    for (int mP=0; mP<num_P_msubs; mP++) {
      for (int mT=0; mT<num_T_msubs; mT++) {
	for (int mE=0; mE<num_E_msubs; mE++) {
	  for (int mR=0; mR<num_msubs; mR++) {
	    f[mP][mT][mE][mR] = complex<double>(0.,0.);
	    f_rot[mP][mT][mE][mR] = complex<double>(0.,0.);
	  }
	}
      }
    }

    fort37_file >> angle[ang];
    angle[ang] *= M_PI/180.;

    for (int mP=0; mP<num_P_msubs; mP++) {
      for (int mT=0; mT<num_T_msubs; mT++) {
	for (int mE=0; mE<num_E_msubs; mE++) {
	  for (int mR=0; mR<num_msubs; mR++) {
	    fort37_file >> real_part >> imag_part;
	    f[mP][mT][mE][mR] = complex<double>(real_part,imag_part);
	  }
	}
      }
    }

    for (int mP=0; mP<num_P_msubs; mP++) {
      for (int mT=0; mT<num_T_msubs; mT++) {
	for (int mE=0; mE<num_E_msubs; mE++) {
	  for (int mR1=0; mR1<num_msubs; mR1++) {
	    double mR1_val = mR1 - JR;
	    for (int mR2=0; mR2<num_msubs; mR2++) {
	      double mR2_val = mR2 - JR;
	      //here is defined the complex conjugate of the Brink & Satchler formulas,
	      //since we're rotating a FINAL state
	      double phase_ang = euler_alpha*mR1_val + euler_gamma*mR2_val;
	      complex<double> phase(cos(phase_ang),sin(phase_ang));
	      f_rot[mP][mT][mE][mR2] += f[mP][mT][mE][mR1]*phase*RotMat(JR,mR1_val,mR2_val);
	    }
	  }
	}
      }
    }

    //construct rho...
    for (int mR1=0; mR1<num_msubs; mR1++) {
      for (int mR2=0; mR2<num_msubs; mR2++) {
	for (int mP=0; mP<num_P_msubs; mP++) {
	  for (int mT=0; mT<num_T_msubs; mT++) {
	    for (int mE=0; mE<num_E_msubs; mE++) {
	      densmat[mR1][mR2][ang] += f_rot[mP][mT][mE][mR1]*conj(f_rot[mP][mT][mE][mR2])*av_fac;
	    }
	  }
	}
      }
    }

    //check for small numbers and set to zero...
    for (int m1=0; m1<num_msubs; m1++) {
      for (int m2=0; m2<num_msubs; m2++) {
	double realpart=real(densmat[m1][m2][ang]);
	double imagpart=imag(densmat[m1][m2][ang]);
	if (abs(realpart) < THRESHOLD) realpart=0;
	if (abs(imagpart) < THRESHOLD) imagpart=0;
	densmat[m1][m2][ang] = complex<double>(realpart,imagpart);
      }
    }

    //construct msubs and xsec...NOT in rotated frame (apparently the only correct way via Jeff T.)
    for (int mP=0; mP<num_P_msubs; mP++) {
      for (int mT=0; mT<num_T_msubs; mT++) {
	for (int mE=0; mE<num_E_msubs; mE++) {
	  for (int mR=0; mR<num_msubs; mR++) {
	    msubs[mR][ang] += norm(f[mP][mT][mE][mR])*fmsq_to_mb*av_fac;
	    xsec[ang]      += norm(f[mP][mT][mE][mR])*fmsq_to_mb*av_fac;
	  }
	}
      }
    }

    //construct msubs_frac...
    for (int m=0; m<num_msubs; m++)
      msubs_frac[m][ang] = msubs[m][ang]/xsec[ang];

  } //end ang loop
  
  //construct integrated msubs...
  dtheta = angle[1] - angle[0];
  for (int ang=0; ang<num_angs; ang++) {
    for (int m=0; m<num_msubs; m++) {
      msubs_integ[m] += msubs[m][ang]*2*M_PI*sin(angle[ang])*dtheta;
    }
    total_xsec += xsec[ang]*2*M_PI*sin(angle[ang])*dtheta;
  }

  //construct integrated msubs_frac...
  for (int m=0; m<num_msubs; m++)
    msubs_integ_frac[m] = msubs_integ[m]/total_xsec;

  //construct Bkq...
  for (int k=0; k<num_ranks; k++) {
    for (int q=0; q<2*k+1; q++) {
      for (int ang=0; ang<num_angs; ang++) {
	for (int m1=0; m1<num_msubs; m1++) {
	  for (int m2=0; m2<num_msubs; m2++) {
	    Bkq[k][q][ang] += 1*
	                      sqrt(2*k+1)*
                              densmat[m1][m2][ang]*
                              clebsch(JR,k,JR,m1-JR,q-k,m2-JR);
	  }
	}
      }
    }
  }
  
  //check for small numbers and set to zero...  
  //and normalize to B00 (real part) if toggled
  for (int k=max_rank; k>=0; k--) {
    for (int q=0; q<2*k+1; q++) {
      for (int ang=0; ang<num_angs; ang++) {
	double realpart, imagpart;
	if (norm_to_B00 && abs(real(Bkq[0][0][ang])) > THRESHOLD) {
	  realpart = real(Bkq[k][q][ang])/real(Bkq[0][0][ang]);
	  imagpart = imag(Bkq[k][q][ang])/real(Bkq[0][0][ang]);
	}
	else {
	  realpart = real(Bkq[k][q][ang]);
	  imagpart = imag(Bkq[k][q][ang]);
	}
	if (abs(realpart) < THRESHOLD) realpart=0;
	if (abs(imagpart) < THRESHOLD) imagpart=0;
	Bkq[k][q][ang]  = complex<double>(realpart,imagpart);
      }
    }
  }
  
  return true;

}

bool BeesHolder::Generate_Tensors(char ifile_path[]) {

  strcpy(file_path,ifile_path);
  return Generate_Tensors();

}

void BeesHolder::Toggle_Norm_to_B00() {norm_to_B00 = !norm_to_B00;}

bool BeesHolder::Norm_to_B00() {return norm_to_B00;}

void BeesHolder::Initialize_Tensors() {

  densmat = new complex<double>**[num_msubs];
  for (int m1=0; m1<num_msubs; m1++) {
    densmat[m1] = new complex<double>*[num_msubs];
    for (int m2=0; m2<num_msubs; m2++) {
      densmat[m1][m2] = new complex<double>[num_angs];
      for (int ang=0; ang<num_angs; ang++) {
	densmat[m1][m2][ang] = complex<double>(0.,0.);
      }
    }
  }

  Bkq = new complex<double>**[num_ranks];
  for (int k=0; k<num_ranks; k++) {
    Bkq[k] = new complex<double>*[2*k+1];
    for (int q=0; q<2*k+1; q++) {
      Bkq[k][q] = new complex<double>[num_angs];
      for (int ang=0; ang<num_angs; ang++) {
	Bkq[k][q][ang] = complex<double>(0.,0.);
      }
    }
  }

}

void BeesHolder::Clear_Tensors() {

  for (int m1=0; m1<num_msubs; m1++) {
    for (int m2=0; m2<num_msubs; m2++) {
      delete [] densmat[m1][m2];
    }
    delete [] densmat[m1];
  }
  delete [] densmat;

  for (int k=0; k<num_ranks; k++) {
    for (int q=0; q<2*k+1; q++) {
      delete [] Bkq[k][q];
    }
    delete [] Bkq[k];
  }
  delete [] Bkq;

}

bool BeesHolder::Check_MSubs_Input(double m1, double m2, int angindex) {
  return (angindex >= 0 && angindex < num_angs &&
	  m1 >= -JR && m1 <= JR &&
	  m2 >= -JR && m2 <= JR) ? true : false;
}

bool BeesHolder::Check_kq_Input(int k, int q, int angindex) {
  return (angindex >= 0 && angindex < num_angs &&
	  k >= 0 && k <= max_rank &&
	  q >= -k && q <= k) ? true : false;
}

double BeesHolder::RotMat(double j, double m, double n) {

  if (j < 0) return 0;
  if (abs(m) > j || abs(n) > j) return 0;

  //Brink & Satchler "Ang. Mom." (ed 3) p.22
  double return_val = 0;

  for (int t=0; ((int)(j+m)-t >= 0 && (int)(j-n)-t >= 0); t++) {
    if ((int)(n-m)+t < 0) continue;
    double sum_val = pow(-1,t);
    sum_val *= sqrt(Factorial((int)(j+m)));
    sum_val *= sqrt(Factorial((int)(j-m)));
    sum_val *= sqrt(Factorial((int)(j+n)));
    sum_val *= sqrt(Factorial((int)(j-n)));
    sum_val /= Factorial((int)(j+m)-t);
    sum_val /= Factorial((int)(j-n)-t);
    sum_val /= Factorial((int)(n-m)+t);
    sum_val /= Factorial(t);
    sum_val *= pow(cos(euler_beta/2),(int)(2*j+m-n)-2*t);
    sum_val *= pow(sin(euler_beta/2),(int)(n-m)+2*t);
    
    return_val += sum_val;
  }

  return return_val;

}

int BeesHolder::Factorial(int n) {

  if (n<0) return 0;

  switch (n) {
  case 0:
  case 1: return 1;
  default: return n*Factorial(n-1);
  }

}
