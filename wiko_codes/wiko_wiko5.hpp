#include <cmath>
#include <complex>
#include <cstdlib>
#include <cstdio>

#include <gsl_sf.h>

#include "wiko_types5.h"
#define MAX_LAMBDA_INDEX 30

using namespace std;

#if 1
double fakul (int a);
double minuseinshoch(int a);
complex<double> ihoch(int a);

double wiko_angular_function(double theta1,double theta2,double phi,
			     int lam1,int lam,int lam2);

complex<double> complex_spherical_harmonic(int lam,int q,double theta,double phi);

double spherical_harmonic(int lam,int q,double theta,double phi);

double wiko_triple_angular_function(double theta1,double theta2,double theta3,
				    double phi1,double phi2,
				    int lam1,int lam2,int lam3,
				    int lam4,int lam5);

double wiko_orient_gaussian(int new_spins,double i1,double sigma,
			    int k1);
double wiko_dumb_orient_gaussian(int new_spins,double i1,double sigma,
			    int k1);
double wiko_orient_single_m(int new_spins,double i1,double sigma,
			    int k1);
double wiko_not_oriented(int new_spins,double i1,double sigma,
			    int k1);

double legen(int k, int q, double tet);
void daslgf_(int *,double *,int *,int *,double *);

double wiko_f3(int l1, int l2, double i2,double i1,
	       int lam1, int lam, int lam2);

double wiko_alpha_f3(int l1, int l2, double i2,double i1,
	       int lam1, int lam, int lam2);

double wiko_prot_f1(double l1,double l2,double i1,double i2,int lam);

double wiko_corr_perturbation(int l1,int l,int l2,
			      double tau_sig,double tau,
			      double k2,double k4,double k6);

double wiko_double_corr(wiko_setup *setup,
			int new_spins,
			double i1,double i2,double i3,double i4,
			double sigma,
			double delta1,double delta2,
			int k1,int k,int k2);

double wiko_triple_corr(wiko_setup *setup,
			int new_spins,
			double i1,double i2,
			double i3,double i4,
			double i5,double i6,
			double sigma,
			double delta1,
			double delta2,
			double delta3,
			int l1,int l2,int l3,int l4,int l5);


double wiko_dis_from_oriented(int new_spins,
			      double i1,double i2,
			      double sigma,
			      double delta,
			      int k);

double wiko_deo_fakt(double i1,double i2,int l,int lambda);

double wiko_deorient(int new_spins,
		     double i1,double i2,
		     double delta,
		     int lambda);

int wiko_init_double_setup(wiko_setup *setup,
			   double (*orient)(int new_spins,double i1,double sigma,int lam),
			   double (*f_1)(int l1,int l2,double i1,double i2,
					 int lam1,int lam,int lam2),
			   double (*f_2)(int l1,int l2,double i1,double i2,
					 int lam1,int lam,int lam2)
			   );
			   
int wiko_init_triple_setup(wiko_setup *setup,
			   double (*orient)(int new_spins,double i1,double sigma,int lam),
			   double (*f_1)(int l1,int l2,double i1,double i2,
					 int lam1,int lam,int lam2),
			   double (*f_2)(int l1,int l2,double i1,double i2,
					 int lam1,int lam,int lam2),
			   double (*f_3)(int l1,int l2,double i1,double i2,
					 int lam1,int lam,int lam2)
			   );

double gg_partial_angular_function(double theta1,
				   double phi1,
				   double theta2,
				   double phi2,
				   int lam0,int lam1,int lam2,
				   int q0);
double partial_g_corr_old(int new_spins,
			  int l0,int q0,
			  double i2,double i3,
			  double theta1,double phi1);
double partial_g_corr(int new_spins,
		      int l0,int q0,
		      double i2,double i3,
		      int l1, int l2,double delta,
		      double theta1,double phi1);
double partial_g_corr_3mult(int new_spins,
			    int l0,int q0,
			    double i2,double i3,
			    int l1, int l2, int l3,
			    double delta21, double delta31,
			    double theta1,double phi1);
double partial_p_corr(int new_spins,
		      int l0,int q0,
		      double i2,double i3,
		      double theta1,double phi1);
double partial_gg_corr(int new_spins,
		       int l0,int q0,
		       double i2,double i3,double i4,
		       double theta1,double phi1,
		       double theta2,double phi2);
double partial_alpha_alpha_corr(int new_spins,
				double *imaginary,
				double i1,double i2,
				int lmin,
				double delta,
				double theta,
				double phi,
				int k,int q);
double partial_alpha_alpha_corr_lower(int new_spins,
				      double *imaginary,
				      double i2,double i3,
				      int lmin,
				      double delta,
				      double theta,
				      double phi,
				      int k,int q);
#else
double wiko_angular_function();
double spherical_harmonic();
double wiko_triple_angular_function();
double wiko_orient_function();
double wiko_dumb_orient_function();
double legen();
void daslgf_();
double wiko_f3();
double wiko_alpha_f3();
double wiko_f1();
double wiko_prot_f1();
double wiko_corr_perturbation();
double wiko_corr_from_oriented();
double wiko_alphagamma_corr_from_oriented();
double wiko_corr_from_oriented_octu();
double wiko_dis_from_oriented();
double wiko_deo_fakt();
double wiko_deorient();
double gg_partial_angular_function();
double partial_g_corr();
double partial_gg_corr();
double partial_alpha_alpha_corr();
#endif



/* clebsch-gordan routinen : */	 
#if 1
double clebsch(double i1,double i,double i2,double m1,double m,double m2);
double dclebg_(double *,double *,double *,double *,double *,double *);

double wign3j(double a ,double b,double c,
	      double xx,double yy,double zz);
double dwig3j_(double *,double *,double *,double *,double *,double *);

double wign6j(double a,double b,double c,
	      double xx,double yy,double zz);
double dwig6j_(double *,double *,double *,double *,double *,double *);

double wign9j(double a,double b,double c,
	      double xx,double yy,double zz,double gg,double hh,double pp);
double dwig9j_(double *,double *,double *,\
	       double *,double *,double *,\
	       double *,double *,double *);

#else
double clebsch();
double dclebg_();

double wign3j();
double dwig3j_();

double wign6j();
double dwig6j_();

double wign9j();
double dwig9j_();
#endif

#define max(a,b) ((a)>(b)?(a):(b))
#define min(a,b) ((a)<(b)?(a):(b))

#define loop_over_lambdas(max_index,lam1,lam,lam2,i)			\
   for (lam=0,i=0;lam<=max_index;lam+=2)				\
     for (lam2=0;lam2<=max_index;lam2+=2)				\
       for (lam1=abs(lam-lam2);lam1<=min(max_index,lam+lam2);lam1+=2,i++)

#define triple_loop_over_lambdas(max_index,lam1,lam2,lam3,lam4,lam5,i)	    \
   for (lam1=0,i=0;lam1<=max_index;lam1+=2)				    \
    for (lam2=0;lam2<=max_index;lam2+=2)				    \
     for (lam3=abs(lam1-lam2);lam3<=min(max_index,lam1+lam2);lam3+=2)	    \
      for (lam4=0;lam4<=max_index;lam4+=2)				    \
       for (lam5=abs(lam3-lam4);lam5<=min(max_index,lam3+lam4);lam5+=2,i++) 

