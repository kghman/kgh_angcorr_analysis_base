#include <cmath>
#include <complex>
#include <cstdlib>
#include <iostream>

#include <gsl_sf.h>

#include "wiko_wiko5.hpp"
#include "WikoClass.h"

using namespace std;

WikoClass::WikoClass(int iparticle_type, double iIi, double iIf, double idelta21, double idelta31) {

  particle_type = iparticle_type;

  if (iIi < 0) Pi_v =-1; //negative parity
  else         Pi_v =+1;

  if (iIf < 0) Pf_v =-1; //negative parity
  else         Pf_v =+1;
  
  if (Pi_v != Pf_v) odd_L = true;
  else              odd_L = false;

  Ii_v = Round_Spin(iIi);
  If_v = Round_Spin(iIf);
  Set_Momenta();
  Set_Max_Rank();

  delta21 = idelta21;
  delta31 = idelta31;

  coeffs = new double[max_rank+1];

  Generate_Coeffs();
}

WikoClass::~WikoClass() {delete [] coeffs;}

int WikoClass::Particle_Type() {return particle_type;}
int WikoClass::Max_Rank()      {return max_rank;}
int WikoClass::L()             {return L_v;}
int WikoClass::LP()            {return LP_v;}
int WikoClass::LPP()           {return LPP_v;}
int WikoClass::Pi()            {return Pi_v;}
int WikoClass::Pf()            {return Pf_v;}
double WikoClass::Ii()         {return Ii_v;}
double WikoClass::If()         {return If_v;}
double WikoClass::Delta21()    {return delta21;}
double WikoClass::Delta31()    {return delta31;}

void WikoClass::Set_Ii(double iIi) {
  if (iIi < 0) Pi_v =-1;
  else         Pi_v =+1;
  if (Pi_v != Pf_v) odd_L = true;
  else              odd_L = false;
  Ii_v = Round_Spin(iIi);
  Set_Momenta();
  Set_Max_Rank();
  Generate_Coeffs();
}
void WikoClass::Set_If(double iIf) {
  if (iIf < 0) Pf_v =-1;
  else         Pf_v =+1;
  if (Pi_v != Pf_v) odd_L = true;
  else              odd_L = false;
  If_v = Round_Spin(iIf);
  Set_Momenta();
  Set_Max_Rank();
  Generate_Coeffs();
}
void WikoClass::Set_Spins(double iIi, double iIf) {
  if (iIi < 0) Pi_v =-1;
  else         Pi_v =+1;
  if (iIf < 0) Pf_v =-1;
  else         Pf_v =+1;
  if (Pi_v != Pf_v) odd_L = true;
  else              odd_L = false;
  Ii_v = Round_Spin(iIi);
  If_v = Round_Spin(iIf);
  Set_Momenta();
  Set_Max_Rank();
  Generate_Coeffs();
}
void WikoClass::Set_Delta21(double idelta21) {
  delta21 = idelta21;
  Generate_Coeffs();
}
void WikoClass::Set_Delta31(double idelta31) {
  delta31 = idelta31;
  Generate_Coeffs();
}
void WikoClass::Set_Deltas(double idelta21, double idelta31) {
  delta21 = idelta21;
  delta31 = idelta31;
  Generate_Coeffs();
}

double WikoClass::Decay_Coefficient(int k) {
  return (Check_kq(k) ? coeffs[k] : 0);
}

double WikoClass::Spherical_Harmonic(int k, int q, double theta, double phi) {
  return (abs(q)<=k ? spherical_harmonic(k,q,theta,phi) : 0);
}

double WikoClass::Spherical_Harmonic_Cosine(int k, int q, double cos_theta, double phi) {
  return (abs(q)<=k ? spherical_harmonic_cosine(k,q,cos_theta,phi) : 0);
}

complex<double> WikoClass::Complex_Spherical_Harmonic(int k, int q, double theta, double phi) {
  return (abs(q)<=k ? complex_spherical_harmonic(k,q,theta,phi) : (complex<double>)(0.,0.));
}

complex<double> WikoClass::Complex_Spherical_Harmonic_Cosine(int k, int q, double cos_theta, double phi) {
  return (abs(q)<=k ? complex_spherical_harmonic_cosine(k,q,cos_theta,phi) : (complex<double>)(0.,0.));
}

bool WikoClass::Check_kq(int k, int q) {

  if (k > max_rank) return false;
  if (abs(q) > k) return false;
  return true;

}

double WikoClass::Round_Spin(double ispin) {

  ispin = abs(ispin);
  ispin = (int)(2*ispin);
  return ispin/2.;

}

void WikoClass::Set_Momenta() {

  if (particle_type != 0) L_v = max(1,(int)(Ii_v-If_v)); //gamma
  else {
    double Idiff = abs(Ii_v-If_v); //ALWAYS odd-half-integer for prot/neut decay
    L_v = (int)(Idiff-0.5+odd_L);
  }

  //only for gamma case (for now), but calculate for both
  LP_v = L_v+1;
  LPP_v = LP_v+1;

}

void WikoClass::Set_Max_Rank() {
  if (particle_type != 0) max_rank = 2*LPP_v; //gamma
  else                    max_rank = (int)(2*(Ii_v+If_v)); //spin-1/2 particle
}

void WikoClass::Generate_Coeffs() {

  delete [] coeffs;
  coeffs = new double[max_rank+1];

  if (particle_type != 0) { //gamma
    for (int k=0; k<=max_rank; k++) {
      coeffs[k] = (wiko_f3(L_v,L_v,Ii_v,If_v,k,k,0) +
		  2*delta21*wiko_f3(L_v,LP_v,Ii_v,If_v,k,k,0) +
		  delta21*delta21*wiko_f3(LP_v,LP_v,Ii_v,If_v,k,k,0) +
		  2*delta31*wiko_f3(L_v,LPP_v,Ii_v,If_v,k,k,0) +
		  delta31*delta31*wiko_f3(LPP_v,LPP_v,Ii_v,If_v,k,k,0) +
		  2*delta21*delta31*wiko_f3(LP_v,LPP_v,Ii_v,If_v,k,k,0));
      coeffs[k] /= (1 + delta21*delta21 + delta31*delta31);
      coeffs[k] /= sqrt(2*k+1);
    }
  }
  else { //spin-1/2 particle -- FOR NOW with a SINGULAR FIXED {j,l}
    double s=0.5; //spin
    if (If_v==0) {
      for (int k=0; k<=max_rank; k++) {
	double j = abs(Ii_v-If_v); //fixed for now, and only include L_v
	double gam = 1; //branching ratio, will divide out in this case but keep track of it
	coeffs[k]  = wign6j(j,j,(double)k,Ii_v,Ii_v,If_v); //W --> 6j using conv of Brink & Satchler
	coeffs[k] *= wign6j((double)L_v,(double)L_v,k,j,j,s);
	coeffs[k] *= clebsch((double)L_v,(double)L_v,(double)k,0,0,0);
	coeffs[k] *= gam*sqrt(2*Ii_v+1)*(2*j+1)*(2*L_v+1);
	coeffs[k] /= sqrt(4*M_PI*(2*k+1));
	coeffs[k] *= pow(-1,(int)(k+If_v+Ii_v+3*s));
      }
    }
    else { //must use the general case of mixed L's -- ONLY IMPLEMENTED IN GDM CLASS
      double j_max = Ii_v + If_v;
      double j_min = abs(Ii_v - If_v);
      for (int k=0; k<=max_rank; k++) {
	coeffs[k] = 0;
	for (double j=j_min; j<=j_max; j++) {
	  for (double jp=j_min; jp<=j_max; jp++) {
	    for (double l=abs(j-s); l<=(j+s); l++) {
	      if (((int)l)%2==0 && odd_L) continue;
	      if (((int)l)%2==1 && !odd_L) continue;
	      for (double lp=abs(jp-s); lp<=(jp+s); lp++) {
		if (((int)lp)%2==0 && odd_L) continue;
		if (((int)lp)%2==1 && !odd_L) continue;
		double GL  = 1; //L partial width
		double GLP = 1; //Lprime partial width
		double CL  = 1; //L phase factor
		double CLP = 1; //Lprime phase factor
		double val = 0;
		val  = wign6j(jp,j,(double)k,Ii_v,Ii_v,If_v);
		val *= wign6j(lp,l,(double)k,j,jp,s);
		val *= clebsch(l,lp,(double)k,0,0,0);
		val *= GL*GLP*CL*CLP;
		val *= sqrt(2*Ii_v+1)*(2*j+1)*(2*l+1)*(2*jp+1)*(2*lp+1);
		val /= sqrt(4*M_PI*(2*k+1));
		val *= pow(-1,(int)(l+lp+2*j+2*jp+If_v+3*Ii_v+s));
		coeffs[k] += val;
		if (val != 0) cout << "k="  << k  << " | "
				   << "j="  << j  << " | "
				   << "jp=" << jp << " | "
				   << "l="  << l  << " | "
				   << "lp=" << lp << endl;		
	      }
	    }
	  }
	}
      }
    }
  }
 
}
