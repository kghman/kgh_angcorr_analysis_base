#include <cmath>
#include <complex>
#include <cstdlib>
#include <iostream>

#include <gsl_sf.h>

#include "wiko_wiko5.hpp"
#include "WikoParticleClass.h"

using namespace std;

WikoParticleClass::WikoParticleClass(double iIi, double iGam) {

  if (iIi < 0) Pi =-1; //negative parity
  else         Pi =+1;

  Gam = (iGam!=0?iGam:1);
  
  If = 0.0; //MAIN ASSUMPTION
  Pf = +1;
    
  if (Pi != Pf) odd_L = true;
  else          odd_L = false;

  Ii = Round_Spin(iIi);
  //the following automatically accounts for If_v=0+
  max_rank = (int)(2*Ii);
  j = Ii;
  s = 0.5;
  int plsval = (int)(j+s);
  int mnsval = (int)(j-s);
  if (odd_L) l = plsval%2==1?plsval:mnsval;
  else       l = plsval%2==0?plsval:mnsval;

  coeffs = new double[max_rank+1];

  Generate_Coeffs();
  
}

WikoParticleClass::~WikoParticleClass() {delete [] coeffs;}

double WikoParticleClass::Get_IiPi()  {return Ii*Pi;}
double WikoParticleClass::Get_BR()    {return Gam;}
int WikoParticleClass::Get_Max_Rank() {return max_rank;}
int WikoParticleClass::Get_L()        {return l;}

void WikoParticleClass::Set_Ii(double iIi) {
  if (iIi < 0) Pi =-1;
  else         Pi =+1;
  if (Pi != Pf) odd_L = true;
  else          odd_L = false;
  Ii = Round_Spin(iIi);
  max_rank = (int)(2*Ii);
  Generate_Coeffs();
}

void WikoParticleClass::Set_Gam(double iGam) {Gam = iGam;}

double WikoParticleClass::Decay_Coefficient(int k) {
  return (Check_kq(k) ? coeffs[k] : 0);
}

double WikoParticleClass::Spherical_Harmonic(int k, int q, double theta, double phi) {
  return (abs(q)<=k ? spherical_harmonic(k,q,theta,phi) : 0);
}

double WikoParticleClass::Spherical_Harmonic_Cosine(int k, int q, double cos_theta, double phi) {
  return (abs(q)<=k ? spherical_harmonic_cosine(k,q,cos_theta,phi) : 0);
}

complex<double> WikoParticleClass::Complex_Spherical_Harmonic(int k, int q, double theta, double phi) {
  return (abs(q)<=k ? complex_spherical_harmonic(k,q,theta,phi) : (complex<double>)(0.,0.));
}

complex<double> WikoParticleClass::Complex_Spherical_Harmonic_Cosine(int k, int q, double cos_theta, double phi) {
  return (abs(q)<=k ? complex_spherical_harmonic_cosine(k,q,cos_theta,phi) : (complex<double>)(0.,0.));
}

bool WikoParticleClass::Check_kq(int k, int q) {

  if (k > max_rank) return false;
  if (abs(q) > k) return false;
  return true;

}

double WikoParticleClass::Round_Spin(double ispin) {

  ispin = abs(ispin);
  ispin = (int)(2*ispin);
  return ispin/2.;

}

void WikoParticleClass::Generate_Coeffs() {

  delete [] coeffs;
  coeffs = new double[max_rank+1];

  for (int k=0; k<=max_rank; k++) {
    coeffs[k]  = wign6j(j,j,(double)k,Ii,Ii,If); //W --> 6j using conv of Brink & Satchler
    coeffs[k] *= wign6j((double)l,(double)l,(double)k,j,j,s);
    coeffs[k] *= clebsch((double)l,(double)l,(double)k,0,0,0);
    coeffs[k] *= Gam*sqrt(2*Ii+1)*(2*j+1)*(2*l+1);
    coeffs[k] /= sqrt(4*M_PI*(2*k+1));
    coeffs[k] *= pow(-1,(int)(k+If+Ii+3*s));
  }
 
}
