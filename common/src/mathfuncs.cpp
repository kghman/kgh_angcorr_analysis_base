#include "mathfuncs.h"

#include <cmath>

using namespace std;

double Gaussian(double E, double M, double FWHM, double scale) {

  return scale*sqrt(4*log(2)/M_PI)/FWHM*exp(-4*log(2)*(E-M)*(E-M)/(FWHM*FWHM));

}

double Breit_Wigner(double E, double M, double Gamma, double scale) {

  return scale*(Gamma/(2*M_PI))/((E-M)*(E-M) + Gamma*Gamma/4);

}

double Breit_Wigner_Rel(double E, double M, double Gamma, double scale) {

  return scale*k_fac(M,Gamma)/((E*E-M*M)*(E*E-M*M) + M*M*Gamma*Gamma);

}

double Gaussian2(double *var_arr, double *param_arr) {

  double E = var_arr[0];
  double fac = param_arr[0];
  double M = param_arr[1];
  double G = param_arr[2];

  return Gaussian(E,M,G,fac);
  
}

double Breit_Wigner2(double *var_arr, double *param_arr) {

  double E = var_arr[0];
  double fac = param_arr[0];
  double M = param_arr[1];
  double G = param_arr[2];

  return Breit_Wigner(E,M,G,fac);
  
}

double Breit_Wigner_Rel2(double *var_arr, double *param_arr) {

  double E = var_arr[0];
  double fac = param_arr[0];
  double M = param_arr[1];
  double G = param_arr[2];

  return Breit_Wigner_Rel(E,M,G,fac);
  
}

double k_fac(double M, double Gamma) {
  
  double k = 2*sqrt(2)*M*M*Gamma*sqrt(M*M+Gamma*Gamma);
  k /= M_PI*sqrt(M*M+M*sqrt(M*M+Gamma*Gamma));
  
  return k;

}
