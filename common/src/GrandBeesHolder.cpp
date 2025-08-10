#include <iostream>
#include <fstream>
#include <cstring>
#include <cstdlib>

#include "wiko_wiko5.hpp"
#include "GrandBeesHolder.h"
#include "constants.h"

#include <gsl_sf.h>

using namespace std;

GrandBeesHolder::GrandBeesHolder(int inum_states,
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
				 double ieu_a,
				 double ieu_b,
				 double ieu_g) : THRESHOLD(1E-10), eye(0.,1.) {

  num_states = inum_states;

  Er = new double[num_states];
  Gamma = new double[num_states];
  BR = new double[num_states];
  eta = new std::complex<double>*[num_states];
  Cphase = new std::complex<double>[num_states];
  JR  = new double[num_states];
  PAR = new int[num_states];
  JD = new double[num_states];
  LD = new int[num_states];
  num_msubs = new int[num_states];
  indiv_wiko = new WikoParticleClass*[num_states];
  indiv_bees = new BeesHolder*[num_states];
  file_paths = new char*[num_states];
  fort37_files = new ifstream[num_states];

  sp = 0.5;
  Jf = 0;

  Ruthparam = new double[num_states];
  E12 = new double[num_states];
  k = new double[num_states];
  Z1 = iZ1;
  Z2 = iZ2;
  M1 = iM1;
  M2 = iM2;
  MU = M1*M2/(M1+M2);
  
  unbound = true;

  for (int i=0; i<num_states; i++) {
    Er[i] = iEr[i];
    Gamma[i] = iGamma[i];
    BR[i] = iBR[i];
    eta[i] = new std::complex<double>[num_states];    
    JR[i]=0;
    PAR[i]=iparity[i]; //either +1 or -1
    JD[i]=0;
    LD[i]=0;

    E12[i] = iE12[i];
    k[i] = sqrt(2*MU*E12[i])/HBAR_C;
    Ruthparam[i] = Z1*Z2*MU*FINE_STRUCTURE/(HBAR_C*k[i]);
    
    num_msubs[i]=0;
    indiv_bees[i] = new BeesHolder(ifile_paths[i],ieu_a,ieu_b,ieu_g);
    file_paths[i] = new char[100];
  }

  for (int i=0; i<num_states; i++)
    if (Gamma[i]==0) unbound = false;
  
  for (int i=0; i<num_states; i++) {
    for (int j=0; j<num_states; j++) {
      eta[i][j] = (unbound?(sqrt(Gamma[i]*Gamma[j])*(0.5*(Gamma[i]+Gamma[j])-eye*(Er[i]-Er[j])))/(0.25*pow(Gamma[i]+Gamma[j],2)+pow(Er[i]-Er[j],2)):std::complex<double>(1.,0.));
    }
  }

  JT=0;
  JP=0;
  JE=0;
  beamE=0;
  max_rank=0;
  num_ranks=1;
  num_tot_msubs=1;
  num_angs=1;
  total_xsec=0;
  angle = new double[num_angs];
  xsec  = new double[num_angs];
  return_arr = new double[num_angs];
  for (int ang=0; ang<num_angs; ang++) {
    angle[ang]=0;
    xsec[ang]=0;
    return_arr=0;
  }
  Initialize_Tensors();

  euler_alpha = ieu_a;
  euler_beta  = ieu_b;
  euler_gamma = ieu_g;

  for (int i=0; i<num_states; i++) {
    strcpy(file_paths[i],ifile_paths[i]);
    fort37_files[i].open(file_paths[i]);
    if (!fort37_files[i]) files_loaded = false;
    else                  files_loaded = true;
  }

}

GrandBeesHolder::~GrandBeesHolder() {

  Clear_Tensors();

  delete [] Er;
  delete [] Gamma;
  delete [] BR;
  delete [] JR;
  delete [] PAR;
  delete [] JD;
  delete [] LD;
  delete [] num_msubs;
  delete [] angle;
  delete [] xsec;
  delete [] return_arr;
  delete [] E12;
  delete [] k;
  delete [] Ruthparam;
  delete [] Cphase;

  delete [] indiv_wiko;
  delete [] indiv_bees;
  for (int i=0; i<num_states; i++) {
    delete [] eta[i];
    delete [] file_paths[i];
    fort37_files[i].close();
  }
  delete [] eta;
  delete [] file_paths;
  delete [] fort37_files;

}

double GrandBeesHolder::Euler_Alpha() {return euler_alpha;}
double GrandBeesHolder::Euler_Beta()  {return euler_beta;}
double GrandBeesHolder::Euler_Gamma() {return euler_gamma;}

int GrandBeesHolder::Num_States()    {return num_states;}
int GrandBeesHolder::Max_Rank()      {return max_rank;}
int GrandBeesHolder::Num_Angles()    {return num_angs;}
int GrandBeesHolder::Num_Tot_MSubs() {return num_tot_msubs;}

double GrandBeesHolder::Get_Er(int state)    {return (state>=0&&state<num_states?Er[state]:0);}
double GrandBeesHolder::Get_Gamma(int state) {return (state>=0&&state<num_states?Gamma[state]:0);}
double GrandBeesHolder::Get_BR(int state)    {return (state>=0&&state<num_states?BR[state]:0);}
double GrandBeesHolder::J(int state)         {return (state>=0&&state<num_states?JR[state]:0);}
int    GrandBeesHolder::Parity(int state)    {return (state>=0&&state<num_states?PAR[state]:0);}
double GrandBeesHolder::Jd(int state)        {return (state>=0&&state<num_states?JD[state]:0);}
double GrandBeesHolder::Ld(int state)        {return (state>=0&&state<num_states?LD[state]:0);}
double GrandBeesHolder::Beam_E()             {return beamE;}
double GrandBeesHolder::Total_XSec()         {return total_xsec;}
double GrandBeesHolder::DTheta()             {return indiv_bees[0]->DTheta();}

std::complex<double> GrandBeesHolder::Get_Eta(int state1, int state2)
{return Check_ijkq_Input(state1,state2,0,0,0) ? eta[state1][state2] : std::complex<double>(0.,0.);}

bool GrandBeesHolder::Files_Loaded() {return files_loaded;}

double GrandBeesHolder::Angle(int angindex)
{return Check_ijkq_Input(0,0,0,0,angindex) ? angle[angindex] : 0.0;}
double GrandBeesHolder::XSec(int angindex)
{return Check_ijkq_Input(0,0,0,0,angindex) ? xsec[angindex] : 0.0;}
double GrandBeesHolder::DensMat_Real(int state1, int state2, double m1, double m2, int angindex)
{return Check_MSubs_Input(state1,state2,m1,m2,angindex) ? real(densmat[state1][state2][(int)(JR[state1]+m1)][(int)(JR[state2]+m2)][angindex]) : 0.0;}
double GrandBeesHolder::DensMat_Imag(int state1, int state2, double m1, double m2, int angindex)
{return Check_MSubs_Input(state1,state2,m1,m2,angindex) ? imag(densmat[state1][state2][(int)(JR[state1]+m1)][(int)(JR[state2]+m2)][angindex]) : 0.0;}
double GrandBeesHolder::Bkq_Real(int state1, int state2, int k, int q, int angindex)
{return Check_ijkq_Input(state1,state2,k,q,angindex) ? real(Bkq[state1][state2][k][k+q][angindex]) : 0.0;}
double GrandBeesHolder::Bkq_Imag(int state1, int state2, int k, int q, int angindex)
{return Check_ijkq_Input(state1,state2,k,q,angindex) ? imag(Bkq[state1][state2][k][k+q][angindex]) : 0.0;}
double GrandBeesHolder::akq_Real(int k, int q, int angindex)
{return Check_ijkq_Input(0,0,k,q,angindex) ? real(akq[k][k+q][angindex]) : 0.0;}
double GrandBeesHolder::akq_Imag(int k, int q, int angindex)
{return Check_ijkq_Input(0,0,k,q,angindex) ? imag(akq[k][k+q][angindex]) : 0.0;}
double GrandBeesHolder::Ak(int state, int k)
{return Check_ijkq_Input(state,0,k,0,0) ? indiv_wiko[state]->Decay_Coefficient(k) : 0.0;}

double* GrandBeesHolder::Angle() {return angle;}
double* GrandBeesHolder::XSec()  {return xsec;}

double* GrandBeesHolder::DensMat_Real(int state1, int state2, double m1, double m2) {
  if (!(Check_MSubs_Input(state1,state2,m1,m2))) return 0;
  for (int ang=0; ang<num_angs; ang++)
    return_arr[ang] = real(densmat[state1][state2][(int)(JR[state1]+m1)][(int)(JR[state2]+m2)][ang]);
  return return_arr;
}

double* GrandBeesHolder::DensMat_Imag(int state1, int state2, double m1, double m2) {
  if (!(Check_MSubs_Input(state1,state2,m1,m2))) return 0;
  for (int ang=0; ang<num_angs; ang++)
    return_arr[ang] = imag(densmat[state1][state2][(int)(JR[state1]+m1)][(int)(JR[state2]+m2)][ang]);
  return return_arr;
}

double* GrandBeesHolder::Bkq_Real(int state1, int state2, int k, int q) {
  if (!(Check_ijkq_Input(state1,state2,k,q))) return 0;
  for (int ang=0; ang<num_angs; ang++)
    return_arr[ang] = real(Bkq[state1][state2][k][k+q][ang]);
  return return_arr;
}

double* GrandBeesHolder::Bkq_Imag(int state1, int state2, int k, int q) {
  if (!(Check_ijkq_Input(state1,state2,k,q))) return 0;
  for (int ang=0; ang<num_angs; ang++)
    return_arr[ang] = imag(Bkq[state1][state2][k][k+q][ang]);
  return return_arr;
}

double* GrandBeesHolder::akq_Real(int k, int q) {
  if (!(Check_ijkq_Input(0,0,k,q))) return 0;
  for (int ang=0; ang<num_angs; ang++)
    return_arr[ang] = real(akq[k][k+q][ang]);
  return return_arr;
}

double* GrandBeesHolder::akq_Imag(int k, int q) {
  if (!(Check_ijkq_Input(0,0,k,q))) return 0;
  for (int ang=0; ang<num_angs; ang++)
    return_arr[ang] = imag(akq[k][k+q][ang]);
  return return_arr;
}

complex<double> GrandBeesHolder::Bkq_Complex(int state1, int state2, int k, int q, int angindex) {
  return (Check_ijkq_Input(state1,state2,k,q,angindex) ? Bkq[state1][state2][k][k+q][angindex] : (complex<double>)(0.,0.));
}

complex<double> GrandBeesHolder::akq_Complex(int k, int q, int angindex) {
  return (Check_ijkq_Input(0,0,k,q,angindex) ? akq[k][k+q][angindex] : (complex<double>)(0.,0.));
}

double GrandBeesHolder::Spherical_Harmonic(int k, int q, double theta, double phi)
{return indiv_wiko[0]->Spherical_Harmonic(k,q,theta,phi);}
double GrandBeesHolder::Spherical_Harmonic_Cosine(int k, int q, double cos_theta, double phi)
{return indiv_wiko[0]->Spherical_Harmonic_Cosine(k,q,cos_theta,phi);}
std::complex<double> GrandBeesHolder::Complex_Spherical_Harmonic(int k, int q, double theta, double phi)
{return indiv_wiko[0]->Complex_Spherical_Harmonic(k,q,theta,phi);}
std::complex<double> GrandBeesHolder::Complex_Spherical_Harmonic_Cosine(int k, int q, double cos_theta, double phi)
{return indiv_wiko[0]->Complex_Spherical_Harmonic_Cosine(k,q,cos_theta,phi);}

void GrandBeesHolder::Set_Euler_Angles(double ieu_a, double ieu_b, double ieu_g) {

  euler_alpha = ieu_a;
  euler_beta  = ieu_b;
  euler_gamma = ieu_g;

}

bool GrandBeesHolder::Set_File_Paths(char **ifile_paths) {
  for (int i=0; i<num_states; i++) {
    strcpy(file_paths[i],ifile_paths[i]);
    fort37_files[i].close();
    fort37_files[i].open(file_paths[i]);
    if (!fort37_files[i]) files_loaded = false;
    else                  files_loaded = true;
  }
  return files_loaded;
}

bool GrandBeesHolder::Generate_Tensors() {

  for (int i=0; i<num_states; i++)
    indiv_bees[i]->Generate_Tensors();

  delete [] angle;
  delete [] xsec;
  delete [] return_arr;

  Clear_Tensors();

  int num_P_msubs = 0;
  int num_T_msubs = 0;
  int num_E_msubs = 0;

  num_tot_msubs = 0;

  double throwaway;
  for (int i=0; i<num_states; i++) {
    //                 ja    jb    jap   jbp    (Tostevin notation)
    fort37_files[i] >> JP >> JT >> JE >> JR[i] >> num_angs >> throwaway >> beamE;
    //these will overwrite, hence we're assuming they're the same for all files
    num_P_msubs = (int)(2*JP+1);
    num_T_msubs = (int)(2*JT+1);
    num_E_msubs = (int)(2*JE+1);
    //these, however, we want to compile all together ~ dimension size
    num_msubs[i] = (int)(2*JR[i]+1);
    num_tot_msubs += num_msubs[i];

    indiv_wiko[i] = new WikoParticleClass(PAR[i]*JR[i],Gamma[i]);
    
    JD[i] = JR[i]; //assumption of decaying to 0+ state
    if (PAR[i]==+1) //no change
      LD[i] = ((int)(JD[i]-sp)%2==0?(int)(JD[i]-sp):(int)(JD[i]+sp));
    else
      LD[i] = ((int)(JD[i]-sp)%2==1?(int)(JD[i]-sp):(int)(JD[i]+sp));

    Cphase[i] = std::complex<double>(exp(-eye*std::complex<double>(LD[i]*M_PI/2,0) + eye*Coulomb_phase_shift(LD[i],Ruthparam[i])));
    
  }
  
  max_rank=0;
  for (int i=0; i<num_states; i++) {
    for (int j=0; j<num_states; j++) {
      if ((int)(JR[i]+JR[j]) > max_rank)
	max_rank = (int)(JR[i]+JR[j]);
    }
  }
  num_ranks = max_rank+1;
  
  total_xsec = 0; //as in angle integrated AND fully coupled

  angle = new double[num_angs];
  xsec  = new double[num_angs]; //fully coupled angle-dependent
  return_arr = new double[num_angs];
  for (int ang=0; ang<num_angs; ang++) {
    angle[ang]=0;
    xsec[ang]=0;
    return_arr[ang]=0;
  }

  Initialize_Tensors();

  double real_part=0, imag_part=0;
  complex<double> ******f_bax = new complex<double>*****[num_states]; //beam axis quant.
  complex<double> ******f_rot = new complex<double>*****[num_states]; //rotated by euler angs
  for (int i=0; i<num_states; i++) {
    f_bax[i] = new complex<double>****[num_angs];
    f_rot[i] = new complex<double>****[num_angs];
    for (int a=0; a<num_angs; a++) {
      f_bax[i][a] = new complex<double>***[num_P_msubs];
      f_rot[i][a] = new complex<double>***[num_P_msubs];
      for (int mP=0; mP<num_P_msubs; mP++) {
	f_bax[i][a][mP] = new complex<double>**[num_T_msubs];
	f_rot[i][a][mP] = new complex<double>**[num_T_msubs];
	for (int mT=0; mT<num_T_msubs; mT++) {
	  f_bax[i][a][mP][mT] = new complex<double>*[num_E_msubs];
	  f_rot[i][a][mP][mT] = new complex<double>*[num_E_msubs];
	  for (int mE=0; mE<num_E_msubs; mE++) {
	    f_bax[i][a][mP][mT][mE] = new complex<double>[num_msubs[i]];
	    f_rot[i][a][mP][mT][mE] = new complex<double>[num_msubs[i]];
	    for (int mR=0; mR<num_msubs[i]; mR++) {
	      f_bax[i][a][mP][mT][mE][mR] = complex<double>(0.,0.);
	      f_rot[i][a][mP][mT][mE][mR] = complex<double>(0.,0.);
	    }
	  }
	}
      }
    }
  }

  double av_fac = 1./((2*JT+1)*(2*JP+1)); //average over initial states
  double fmsq_to_mb = 10.;

  for (int i=0; i<num_states; i++) {

    for (int ang=0; ang<num_angs; ang++) {

      fort37_files[i] >> angle[ang]; //it's presumed that each state has the same angle binning
      angle[ang] *= M_PI/180.;

      for (int mP=0; mP<num_P_msubs; mP++) {
	for (int mT=0; mT<num_T_msubs; mT++) {
	  for (int mE=0; mE<num_E_msubs; mE++) {
	    for (int mR=0; mR<num_msubs[i]; mR++) {
	      fort37_files[i] >> real_part >> imag_part;
	      f_bax[i][ang][mP][mT][mE][mR] = complex<double>(real_part,imag_part);
	    }
	  }
	}
      }

      for (int mP=0; mP<num_P_msubs; mP++) {
	for (int mT=0; mT<num_T_msubs; mT++) {
	  for (int mE=0; mE<num_E_msubs; mE++) {
	    for (int mR1=0; mR1<num_msubs[i]; mR1++) {
	      double mR1_val = mR1 - JR[i];
	      for (int mR2=0; mR2<num_msubs[i]; mR2++) {
		double mR2_val = mR2 - JR[i];
		//here is defined the complex conjugate of the Brink & Satchler formulas,
		//since we're rotating a FINAL state
		double phase_ang = euler_alpha*mR1_val + euler_gamma*mR2_val;
		complex<double> phase(cos(phase_ang),sin(phase_ang));
		f_rot[i][ang][mP][mT][mE][mR2] += f_bax[i][ang][mP][mT][mE][mR1]*phase*RotMat(JR[i],mR1_val,mR2_val);
	      }
	    }
	  }
	}
      }

    }
    
  }
  
  //now construct E-INTEGRATED rho...
  for (int i=0; i<num_states; i++) {
    for (int j=0; j<num_states; j++) {
      for (int ang=0; ang<num_angs; ang++) {
	for (int mR1=0; mR1<num_msubs[i]; mR1++) {
	  for (int mR2=0; mR2<num_msubs[j]; mR2++) {
	    for (int mP=0; mP<num_P_msubs; mP++) {
	      for (int mT=0; mT<num_T_msubs; mT++) {
		for (int mE=0; mE<num_E_msubs; mE++) {
		  densmat[i][j][mR1][mR2][ang] += f_rot[i][ang][mP][mT][mE][mR1]*conj(f_rot[j][ang][mP][mT][mE][mR2])*eta[i][j]*av_fac;
		}
	      }
	    }
	  }
	}
      }
    }
  }

  //check for small numbers and set to zero...
  for (int i=0; i<num_states; i++) {
    for (int j=0; j<num_states; j++) {
      for (int ang=0; ang<num_angs; ang++) {
	for (int m1=0; m1<num_msubs[i]; m1++) {
	  for (int m2=0; m2<num_msubs[j]; m2++) {
	    double realpart=real(densmat[i][j][m1][m2][ang]);
	    double imagpart=imag(densmat[i][j][m1][m2][ang]);
	    if (abs(realpart) < THRESHOLD) realpart=0;
	    if (abs(imagpart) < THRESHOLD) imagpart=0;
	    densmat[i][j][m1][m2][ang] = complex<double>(realpart,imagpart);
	  }
	}
      }
    }
  }

  //construct Bkq...
  for (int i=0; i<num_states; i++) {
    for (int j=0; j<num_states; j++) {
      for (int k=0; k<num_ranks; k++) {
	for (int q=0; q<2*k+1; q++) {
	  for (int ang=0; ang<num_angs; ang++) {
	    for (int m1=0; m1<num_msubs[i]; m1++) {
	      for (int m2=0; m2<num_msubs[j]; m2++) {
		Bkq[i][j][k][q][ang] += 1*
		                        sqrt(2*k+1)*
		                        densmat[i][j][m1][m2][ang]*
		                        clebsch(JR[i],k,JR[j],m1-JR[i],q-k,m2-JR[j]);
	      }
	    }
	  }
	}
      }
    }
  }
  
  //check for small numbers and set to zero...
  for (int i=0; i<num_states; i++) {
    for (int j=0; j<num_states; j++) {
      for (int k=max_rank; k>=0; k--) {
	for (int q=0; q<2*k+1; q++) {
	  for (int ang=0; ang<num_angs; ang++) {
	    double realpart, imagpart;
	    realpart = real(Bkq[i][j][k][q][ang]);
	    imagpart = imag(Bkq[i][j][k][q][ang]);
	    if (abs(realpart) < THRESHOLD) realpart=0;
	    if (abs(imagpart) < THRESHOLD) imagpart=0;
	    Bkq[i][j][k][q][ang] = complex<double>(realpart,imagpart);
	  }
	}
      }
    }
  }
  
  //form xsec out of B00 terms
  double dtheta = angle[1]-angle[0];
  for (int i=0; i<num_states; i++) {
    for (int j=0; j<num_states; j++) {
      for (int ang=0; ang<num_angs; ang++) {
	xsec[ang]  += real(Bkq[i][j][0][0][ang])*fmsq_to_mb;
	total_xsec += real(Bkq[i][j][0][0][ang])*fmsq_to_mb*2*M_PI*sin(angle[ang])*dtheta;
      }
    }
  }
  
  //then finally construct the coupled akq which form W(theta,phi) = sum_{kq}a_{kq}Y_{kq}(theta,phi)
  for (int i=0; i<num_states; i++) {
    for (int j=0; j<num_states; j++) {
      for (int k=0; k<num_ranks; k++) {
	for (int q=0; q<2*k+1; q++) {
	  for (int ang=0; ang<num_angs; ang++) {
	    akq[k][q][ang] += 1/sqrt(4*M_PI*(2*k+1))*
	      pow(-1,(int)(k+Jf+JR[i]+3*sp))* //Brink & Satchler phases for Wigner-6j conversion
	      sqrt(BR[i]*BR[j])*
	      Cphase[i]*conj(Cphase[j])*
	      sqrt((2*JR[i]+1)*(2*JD[i]+1)*(2*JD[j]+1)*(2*LD[i]+1)*(2*LD[j]+1))*
	      Bkq[i][j][k][q][ang]*
	      wign6j(JD[j],JD[i],(double)k,JR[i],JR[j],Jf)*
              wign6j(LD[j],LD[i],(double)k,JD[i],JD[j],sp)*
	      clebsch(LD[j],LD[i],k,0,0,0);
	  }
	}
      }
    }
  }
  
  for (int i=0; i<num_states; i++) {
    for (int a=0; a<num_angs; a++) {
      for (int mP=0; mP<num_P_msubs; mP++) {
	for (int mT=0; mT<num_T_msubs; mT++) {
	  for (int mE=0; mE<num_E_msubs; mE++) {
	    delete [] f_bax[i][a][mP][mT][mE];
	    delete [] f_rot[i][a][mP][mT][mE];
	  }
	  delete [] f_bax[i][a][mP][mT];
	  delete [] f_rot[i][a][mP][mT];       
	}
	delete [] f_bax[i][a][mP];
	delete [] f_rot[i][a][mP];
      }
      delete [] f_bax[i][a];
      delete [] f_rot[i][a];
    }
    delete [] f_bax[i];
    delete [] f_rot[i];
  }
  delete [] f_bax;
  delete [] f_rot;

  return true;

}

bool GrandBeesHolder::Generate_Tensors(char **ifile_paths) {

  for (int i=0; i<num_states; i++)
    strcpy(file_paths[i],ifile_paths[i]);
  return Generate_Tensors();

}

void GrandBeesHolder::Initialize_Tensors() {

  densmat = new complex<double>****[num_states];
  for (int i=0; i<num_states; i++) {
    densmat[i] = new complex<double>***[num_states];
    for (int j=0; j<num_states; j++) {
      densmat[i][j] = new complex<double>**[num_msubs[i]];
      for (int m1=0; m1<num_msubs[i]; m1++) {
	densmat[i][j][m1] = new complex<double>*[num_msubs[j]];
	for (int m2=0; m2<num_msubs[j]; m2++) {
	  densmat[i][j][m1][m2] = new complex<double>[num_angs];
	  for (int ang=0; ang<num_angs; ang++) {
	    densmat[i][j][m1][m2][ang] = complex<double>(0.,0.);
	  }
	}
      }
    }
  }

  Bkq = new complex<double>****[num_states];
  for (int i=0; i<num_states; i++) {
    Bkq[i] = new complex<double>***[num_states];
    for (int j=0; j<num_states; j++) {
      Bkq[i][j] = new complex<double>**[num_ranks];
      for (int k=0; k<num_ranks; k++) {
	Bkq[i][j][k] = new complex<double>*[2*k+1];
	for (int q=0; q<2*k+1; q++) {
	  Bkq[i][j][k][q] = new complex<double>[num_angs];
	  for (int ang=0; ang<num_angs; ang++) {
	    Bkq[i][j][k][q][ang] = complex<double>(0.,0.);
	  }
	}
      }
    }
  }

  akq = new complex<double>**[num_ranks];
  for (int k=0; k<num_ranks; k++) {
    akq[k] = new complex<double>*[2*k+1];
    for (int q=0; q<2*k+1; q++) {
      akq[k][q] = new complex<double>[num_angs];
      for (int ang=0; ang<num_angs; ang++) {
	akq[k][q][ang] = complex<double>(0.,0.);
      }
    }
  }

}

void GrandBeesHolder::Clear_Tensors() {

  for (int i=0; i<num_states; i++) {
    for (int j=0; j<num_states; j++) {
      for (int m1=0; m1<num_msubs[i]; m1++) {
	for (int m2=0; m2<num_msubs[j]; m2++) {
	  delete [] densmat[i][j][m1][m2];
	}
	delete [] densmat[i][j][m1];
      }
      delete [] densmat[i][j];
    }
    delete [] densmat[i];
  }
  delete [] densmat;

  for (int i=0; i<num_states; i++) {
    for (int j=0; j<num_states; j++) {
      for (int k=0; k<num_ranks; k++) {
	for (int q=0; q<2*k+1; q++) {
	  delete [] Bkq[i][j][k][q];
	}
	delete [] Bkq[i][j][k];
      }
      delete [] Bkq[i][j];
    }
    delete [] Bkq[i];
  }
  delete [] Bkq;

  for (int k=0; k<num_ranks; k++) {
    for (int q=0; q<2*k+1; q++) {
      delete [] akq[k][q];
    }
    delete [] akq[k];
  }
  delete [] akq;

}

bool GrandBeesHolder::Check_MSubs_Input(int state1, int state2, double m1, double m2, int angindex) {
  if (state1 < 0 || state2 < 0) return false;
  return (angindex >= 0 && angindex < num_angs &&
	  abs(m1) <= JR[state1] &&
	  abs(m2) <= JR[state2] &&
	  state1 < num_states &&
	  state2 < num_states) ? true : false;
}

bool GrandBeesHolder::Check_ijkq_Input(int state1, int state2, int k, int q, int angindex) {
  return (angindex >= 0 && angindex < num_angs &&
	  k >= 0 && k <= max_rank &&
	  q >= -k && q <= k &&
	  state1 >= 0 && state1 < num_states &&
	  state2 >= 0 && state2 < num_states) ? true : false;
}

double GrandBeesHolder::RotMat(double j, double m, double n) {

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

int GrandBeesHolder::Factorial(int n) {

  if (n<0) return 0;

  switch (n) {
  case 0:
  case 1: return 1;
  default: return n*Factorial(n-1);
  }

}

std::complex<double> GrandBeesHolder::Coulomb_phase_shift(int l, double Ruthparam) {

  std::complex<double> eye(0.,1.);
  std::complex<double> neyeovertwo(0.,-0.5);

  gsl_sf_result gam_top_lnr;
  gsl_sf_result gam_top_arg;
  gsl_sf_result gam_bot_lnr;
  gsl_sf_result gam_bot_arg;

  gsl_sf_lngamma_complex_e(l+1, Ruthparam,&gam_top_lnr,&gam_top_arg);
  gsl_sf_lngamma_complex_e(l+1,-Ruthparam,&gam_bot_lnr,&gam_bot_arg);
  return neyeovertwo*(gam_top_lnr.val - gam_bot_lnr.val + eye*(gam_top_arg.val - gam_bot_arg.val));
  
}
