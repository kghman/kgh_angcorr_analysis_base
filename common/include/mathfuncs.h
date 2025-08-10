#ifndef __MATHFUNCS_H
#define __MATHFUNCS_H

double Gaussian(double E, double M, double FWHM, double scale);
double Breit_Wigner(double E, double M, double Gamma, double scale);
double Breit_Wigner_Rel(double E, double M, double Gamma, double scale);
double Gaussian2(double *var_arr, double *param_arr);
double Breit_Wigner2(double *var_arr, double *param_arr);
double Breit_Wigner_Rel2(double *var_arr, double *param_arr);

double k_fac(double M, double Gamma);

#endif
