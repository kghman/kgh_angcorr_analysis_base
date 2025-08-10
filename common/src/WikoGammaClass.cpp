#include <cmath>
#include <complex>
#include <cstdlib>
#include <iostream>

#include <gsl_sf.h>

#include "wiko_wiko5.hpp"
#include "WikoGammaClass.h"

using namespace std;

WikoGammaClass::WikoGammaClass(double iIi, double iIf, double idelta21, double idelta31) {

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

WikoGammaClass::~WikoGammaClass() {delete [] coeffs;}

int WikoGammaClass::Get_Max_Rank()  {return max_rank;}
int WikoGammaClass::L()             {return L_v;}
int WikoGammaClass::LP()            {return LP_v;}
int WikoGammaClass::LPP()           {return LPP_v;}
int WikoGammaClass::Pi()            {return Pi_v;}
int WikoGammaClass::Pf()            {return Pf_v;}
double WikoGammaClass::Ii()         {return Ii_v;}
double WikoGammaClass::If()         {return If_v;}
double WikoGammaClass::Delta21()    {return delta21;}
double WikoGammaClass::Delta31()    {return delta31;}

void WikoGammaClass::Set_Ii(double iIi) {
  if (iIi < 0) Pi_v =-1;
  else         Pi_v =+1;
  if (Pi_v != Pf_v) odd_L = true;
  else              odd_L = false;
  Ii_v = Round_Spin(iIi);
  Set_Momenta();
  Set_Max_Rank();
  Generate_Coeffs();
}
void WikoGammaClass::Set_If(double iIf) {
  if (iIf < 0) Pf_v =-1;
  else         Pf_v =+1;
  if (Pi_v != Pf_v) odd_L = true;
  else              odd_L = false;
  If_v = Round_Spin(iIf);
  Set_Momenta();
  Set_Max_Rank();
  Generate_Coeffs();
}
void WikoGammaClass::Set_Spins(double iIi, double iIf) {
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
void WikoGammaClass::Set_Delta21(double idelta21) {
  delta21 = idelta21;
  Generate_Coeffs();
}
void WikoGammaClass::Set_Delta31(double idelta31) {
  delta31 = idelta31;
  Generate_Coeffs();
}
void WikoGammaClass::Set_Deltas(double idelta21, double idelta31) {
  delta21 = idelta21;
  delta31 = idelta31;
  Generate_Coeffs();
}

double WikoGammaClass::Decay_Coefficient(int k) {
  return (Check_kq(k) ? coeffs[k] : 0);
}

double WikoGammaClass::Spherical_Harmonic(int k, int q, double theta, double phi) {
  return (abs(q)<=k ? spherical_harmonic(k,q,theta,phi) : 0);
}

double WikoGammaClass::Spherical_Harmonic_Cosine(int k, int q, double cos_theta, double phi) {
  return (abs(q)<=k ? spherical_harmonic_cosine(k,q,cos_theta,phi) : 0);
}

complex<double> WikoGammaClass::Complex_Spherical_Harmonic(int k, int q, double theta, double phi) {
  return (abs(q)<=k ? complex_spherical_harmonic(k,q,theta,phi) : (complex<double>)(0.,0.));
}

complex<double> WikoGammaClass::Complex_Spherical_Harmonic_Cosine(int k, int q, double cos_theta, double phi) {
  return (abs(q)<=k ? complex_spherical_harmonic_cosine(k,q,cos_theta,phi) : (complex<double>)(0.,0.));
}

bool WikoGammaClass::Check_kq(int k, int q) {

  if (k > max_rank) return false;
  if (abs(q) > k) return false;
  return true;

}

double WikoGammaClass::Round_Spin(double ispin) {

  ispin = abs(ispin);
  ispin = (int)(2*ispin);
  return ispin/2.;

}

void WikoGammaClass::Set_Momenta() {

  L_v = max(1,(int)(Ii_v-If_v));
  LP_v = L_v+1;
  LPP_v = LP_v+1;

}

void WikoGammaClass::Set_Max_Rank() {
  max_rank = 2*LPP_v; //gamma
}

void WikoGammaClass::Generate_Coeffs() {

  delete [] coeffs;
  coeffs = new double[max_rank+1];

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
