#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <complex.h>

#include <gsl/gsl_sf.h>
#include "wiko_wiko5.h"

/* angular correlation subroutines in C (finally) */

#define min(a,b) ((a)<(b)?(a):(b))
#define d(a,b) ((a)==(b)?(1.):(0.))

/***********************/
/* auxillary functions */
/***********************/

double fakul(a)
int a;
{
  static int first_time=1;
  static double faks[50];
  int i;
  double fak; 
  if (first_time){
    first_time=0;
    faks[0]=1;
    for (fak=1.,i=1;i<50;i++){
      fak*=(double )i;
      faks[i]=fak;
    }
  }
  if (a>=50){
    fprintf(stderr,"factorial coefficient %d not supported\n",a);
  }

  return(faks[a]);
}

double minuseinshoch(a)
int a;
{
  if (abs(a)%2)
    return(-1.);
  else
    return(1.);
}

double complex ihoch(int a)
{
  switch(abs(a)%4)
    {
    case 0: return(1);
    case 1: return(I);
    case 2: return(-1);
    case 3: return(-I);
    }
}

/***********************/
/* geometry functions **/
/***********************/

/* angular function for tensor koefficients lam1,lam,lam2 */
/* see KSW: formula 24 */
double wiko_angular_function(theta1,theta2,phi,lam1,lam,lam2)
     double theta1,theta2,phi;
     int lam1,lam,lam2;
{
  double q;
  int lamstr;
  double sum=0.;
  lamstr=min(lam,lam2);
  for (q=0.;q<=lamstr;q+=1.){
    sum+= (2.-d(q,0.)) * clebsch((double)lam1,(double)lam,(double)lam2,
				 0.,(double)q,(double)q) * 
      sqrt((2.*lam+1)/(2.*lam2+1) * 
	   (fakul((int)(lam-q))*fakul((int)(lam2-q))) 
	   / (fakul((int)(lam+q))*fakul((int)(lam2+q)))
	   )* 
      legen(lam,(int)q,theta1) * legen(lam2,(int)q,theta2) * 
      cos (q*phi*M_PI/180.);
  }
  return(sum);
}


double spherical_harmonic(lam,q,theta,phi)
int lam,q;
double theta,phi;
{
  double value;
  value=1./sqrt(4.*M_PI);
  if (q>=0){
    value *= (minuseinshoch(q) * cos(q*phi*M_PI/180.));
  }
  else {
    value *= cos(q*phi*M_PI/180.);
  }
  value *= sqrt((2*lam+1)*fakul(lam-abs(q))/fakul(lam+abs(q)));
  return(value*legen(lam,abs(q),theta));
}

double complex complex_spherical_harmonic(lam,q,theta,phi)
int lam,q;
double theta,phi;
{
  double complex cvalue;

  if (abs(q)>lam)
    return(0);
  cvalue = cexp(I*phi*q*M_PI/180.)/sqrt(4.*M_PI);
  if (q>=0){
    cvalue *= minuseinshoch(q);
  }
  cvalue *= sqrt((2*lam+1)*fakul(lam-abs(q))/fakul(lam+abs(q)));
  cvalue*= gsl_sf_legendre_Plm (lam,abs(q),cos(theta*M_PI/180.));
  return(cvalue);
}

/* angular function for triple-DCO, still experimental */

double wiko_triple_angular_function(theta1,theta2,theta3,phi1,phi2,
				    lam1,lam2,lam3,lam4,lam5)
double theta1,theta2,theta3,phi1,phi2;
int lam1,lam2,lam3,lam4,lam5;
{
  int q1,q2,q3;
  double sum=0.,fak=0.;
  for (q1=0;q1<=min(lam2,lam3);q1++){
    for (q2=-lam4;q2<=lam4;q2++){
      q3=-q1-q2;
      if (abs(q3)<=lam5){
	fak =1;
	fak *= minuseinshoch(q1);
	fak /= sqrt(2.*lam5+1);
	fak *=(2.-d(q1,0.));
	fak *=wign3j((double)lam3,(double)lam2,(double)lam1,
		     (double)-q1,(double)q1,0.);
	fak *=wign3j((double)lam5,(double)lam4,(double)lam3,
		     (double)q3,(double)q2,(double)q1);
	fak *=spherical_harmonic(lam2,q1,theta1,0.);
	fak *=spherical_harmonic(lam4,q2,theta2,phi1);
	fak *=spherical_harmonic(lam5,q3,theta3,phi2);
	fak *=sqrt(64*M_PI*M_PI*M_PI);
	sum +=fak;
      }
    }
  }	
  return(sum);
}

/***********************/
/*orientation functions*/
/***********************/


/* calculates the orientation tensor B_lamda of
   state i1 gaussian m-population with width sigma */

double wiko_orient_gaussian(new_spins,i1,sigma,lam)
     int new_spins;double i1;double sigma;
     int lam;
{
  double m;
  int l;
  double sum=0.,norm=0.;
  static double *c[MAX_LAMBDA_INDEX/2+1];
  static double max_spin=-1.;
  static double last_sigma,b[MAX_LAMBDA_INDEX/2+1];
  
  if (new_spins){
    if (i1>max_spin){
      for (l=0;l<=MAX_LAMBDA_INDEX;l+=2){
	if (c[l/2])
	  free(c[l/2]);
	c[l/2]=(double *)malloc(((int)i1+1)*sizeof(double));
      }
      max_spin=i1;
    }
    for (l=0;l<=MAX_LAMBDA_INDEX;l+=2){
      for (m=i1;m>=0;m-=1.){
	c[l/2][(int)m]=
	  clebsch(i1,i1,(double)l,(double)-m,(double)m,0.);
      }
    }
  }
  if ((new_spins)||(last_sigma!=sigma)){
    for (l=0;l<=MAX_LAMBDA_INDEX;l+=2){
      sum=0;norm=0;
      for (m=i1;m>=0;m-=1.){
	sum+= (2-d(m,0))*minuseinshoch((int)(i1+m))*
	  c[l/2][(int)m]*
	  exp(-(m*m)/(2*sigma*sigma));
	norm+=(2-d(m,0))*exp(-(m*m)/(2*sigma*sigma));
      }
      b[l/2]=sqrt(2.*i1+1)*(sum/norm);    
    }
    last_sigma=sigma;
  }
  return(b[lam/2]);
}

double wiko_orient_single_m(new_spins,i1,them,lam)
     int new_spins;double i1;double them;
     int lam;
{
  int l;
  double m;
  static double *c[MAX_LAMBDA_INDEX/2+1];
  static double max_spin=-1.;
  static double last_m,b[MAX_LAMBDA_INDEX/2+1];
  
  if (new_spins){
    if (i1>max_spin){
      for (l=0;l<=MAX_LAMBDA_INDEX;l+=2){
	if (c[l/2])
	  free(c[l/2]);
	c[l/2]=(double *)malloc(((int)i1+1)*sizeof(double));
      }
      max_spin=i1;
    }
    for (l=0;l<=MAX_LAMBDA_INDEX;l+=2){
      for (m=i1;m>=0.;m-=1.){
	c[l/2][(int)m]=
	  wign3j(i1,i1,(double)l,(double)-m,(double)m,0.)*sqrt(2*l+1);
      }
    }
  }
  if ((new_spins)||(last_m!=them)){
    for (l=0;l<=MAX_LAMBDA_INDEX;l+=2){
      b[l/2]=sqrt((2.*i1+1.))
	*minuseinshoch(i1+them)*c[l/2][(int)them];
    }
    last_m=them;
  }
  return(b[lam/2]);
}

double wiko_dumb_orient_gaussian(new_spins,i1,sigma,lam)
     int new_spins;double i1;double sigma;
     int lam;
{
  double m;
  double sum=0.,norm=0.;
  
  sum=0;norm=0;
  for (m=i1;m>=-i1;m-=1.){
    sum+= minuseinshoch((int)(i1+m))*
      clebsch(i1,i1,(double)lam,-m,m,0.)*
      exp(-(m*m)/(2*sigma*sigma));
    norm+=exp(-(m*m)/(2*sigma*sigma));
  }
  return(sqrt(2.*i1+1.)*(sum/norm));
}

double wiko_not_oriented(new_spins,i1,sigma,lam)
     int new_spins;double i1;double sigma;
     int lam;
{
  if (lam==0)
    return(1.);
  else
    return(0.);
}

/*****************************************/
/*f-coefficients for different particles */
/******************************************/

/* generalized F-Coefficient for gamma transition */
/* see KSW: formula 46 l1=l l2=l'*/
double wiko_f3(l1,l2,i1,i2,lam1,lam,lam2)
     int l1; int l2; double i1;double i2;
     int lam1; int lam; int lam2;
{
  if (lam2!=0){
    /* real generalized f-coefficient */
    return(sqrt((2.*i1+1)*(2.*i2+1)*(2.*l1+1)*(2.*l2+1)*
		(2.*lam+1)*(2*lam1+1)*(2*lam2+1))*
	   minuseinshoch(l2+lam1+lam2+1)*
	   wign3j((double)l1,(double)l2,(double)lam,1.,-1.,0.) *
	   wign9j(i2,(double)l1,i1,
		  i2,(double)l2,i1,
		  (double)lam2,(double)lam,(double)lam1));
  }
  else{
    /* reduces to ordinary f-coefficient, speed */
    if (lam1==lam){
      return(sqrt((2.*i1+1)*(2.*l1+1)*(2.*l2+1)*(2.*lam+1))*
	     minuseinshoch((int)(i1+i2)+1)*
	     wign3j((double)l1,(double)l2,(double)lam,1.,-1.,0.) *
	     wign6j((double)l1,(double)l2,(double)lam,i1,i1,i2));
    }
    else{
      return(0.);
    }
  }
}

/* calculates the correlation tensor for alpha particles */
/* b_lambda*F_lambda: see Hamilton: formula 12.338-12.341 */
double wiko_alpha_f3(l1,l2,i1,i2,lam1,lam,lam2)
     int l1; int l2; double i1;double i2;
     int lam1; int lam; int lam2;
{
  double r1,r2,r3,r4;

  /* note that there is a factor -1 in front of the alpha-particle coefficient */
  /* see Hamilton 12.340 */
  if (lam2!=0){
    r1=-1.*sqrt((2.*i1+1)*(2.*i2+1)*(2.*l1+1)*(2.*l2+1)*
	    (2.*lam+1)*(2*lam1+1)*(2*lam2+1));
    r2=minuseinshoch(l2+lam1+lam2+1);
    r3=wign3j((double)l1,(double)l2,(double)lam,0.,0.,0.);
    r4=wign9j(i2,(double)l1,i1,
	      i2,(double)l2,i1,
	      (double)lam2,(double)lam,(double)lam1);
  }else{
    if (lam1==lam){
      r1=-1.*sqrt((2.*i1+1)*(2.*l1+1)*(2.*l2+1)*(2.*lam+1));
      r2=minuseinshoch((int)(i1+i2)-1);
      r3=wign3j((double)l1,(double)l2,(double)lam,0.,0.,0.);
      r4=wign6j((double)l1,(double)l2,(double)lam,i1,i1,i2);
    } 
    else{
      return(0);
    }
  }
  return(r1*r2*r3*r4);

}

/* hi, MSU! */
double wiko_prot_f1(j1,j2,i1,i2,lam)
double j1; double j2; double i1; double i2;int lam;
{  
  return(sqrt((2.*i1+1)*(2.*j1+1)*(2.*j2+1)*(2.*lam+1))*
	 minuseinshoch((int)(i1+i2)+1)*
	 wign3j(j1,j2,(double)lam,0.5,-0.5,0.) *
	 wign6j(j1,j2,(double)lam,i1,i1,i2));
}


double wiko_corr_perturbation(l1,l,l2,tau_sig,tau,k2,k4,k6)
     int l1;int l;int l2;
     double tau_sig;double tau;
     double k2;double k4;double k6;
{
  double fak=1.;
  
  if (l1==2){
    fak *= 1./(1.+k2*tau_sig);
  }
  if (l1==4){
    fak *= 1./(1.+k4*tau_sig);
  }
  if (l2==2){
    fak *= 1./(1.+k2*tau);
  }
  if (l2==4){
    fak *= 1./(1.+k4*tau);
  }
  return(fak);
}

/* calculates the correlation coefficient for two-fold correlation from 
   oriented state */

int wiko_init_double_setup(wiko_setup *setup,
			   double (*orient)(),
			   double (*f_1)(),
			   double (*f_2)()
			   )
{
  setup->orient_function=orient;
  setup->fold=2;
  (setup->f_function1)=f_1;
  (setup->f_function2)=f_2;
  return(0);
}

int wiko_init_triple_setup(wiko_setup *setup,
			   double (*orient)(),
			   double (*f_1)(),
			   double (*f_2)(),
			   double (*f_3)()
			   )
{
  setup->orient_function=orient;
  setup->fold=3;
  setup->f_function1=f_1;
  setup->f_function2=f_2;
  setup->f_function3=f_3;
  return(0);
}


double wiko_double_corr(setup,new_spins,i1,i2,i3,i4,sigma,delta1,delta2,
			k1,k,k2)
     wiko_setup *setup;
     int new_spins;
     double i1;double i2;double i3;double i4;
     double sigma;
     double delta1;double delta2;
     int k1;int k;int k2;
{
  static double f3[3][MAX_LAMBDA_INDEX/2+1][MAX_LAMBDA_INDEX/2+1][MAX_LAMBDA_INDEX/2+1];
  static double f1[3][MAX_LAMBDA_INDEX/2+1];
  double a1,a3;
  double b;
  
  int max_lambda;
  int lam1,lam,lam2,i;
  int l1,l,l2;
  int l1a,l1b,l2a,l2b;

  l1a=max(1,(int)(fabs(i2-i1)));
  l1b=l1a+1;
  l2a=max(1,(int)(fabs(i4-i3)));
  l2b=l2a+1;
  max_lambda = 2* max(i1,max(l1b,l2b));
  if (new_spins){
    /* calculate arrays of F-coefficients for the given spins and mults. */
    /* the variation of delta's and sigmas can then be calculated without */
    /* repeating this initialisation (set new_spins=0) */
    loop_over_lambdas(max_lambda,lam1,lam,lam2,i){
      l1=lam1>>1;
      l=lam>>1;
      l2=lam2>>1;

      f3[0][l1][l][l2]=(setup->f_function1)(l1a,l1a,i1,i2,lam1,lam,lam2);
      f3[1][l1][l][l2]=(setup->f_function1)(l1a,l1b,i1,i2,lam1,lam,lam2);
      f3[2][l1][l][l2]=(setup->f_function1)(l1b,l1b,i1,i2,lam1,lam,lam2);

      f1[0][l2]=setup->f_function2(l2a,l2a,i3,i4,lam2,lam2,0.);
      f1[1][l2]=setup->f_function2(l2a,l2b,i3,i4,lam2,lam2,0.);
      f1[2][l2]=setup->f_function2(l2b,l2b,i3,i4,lam2,lam2,0.);

    }
  }
  l1=k1>>1;
  l=k>>1;
  l2=k2>>1;
  b=setup->orient_function(new_spins,i1,sigma,k1);
  a3=(1./(1+delta1*delta1)*
      (f3[0][l1][l][l2]+
       f3[1][l1][l][l2]*2*delta1+
       f3[2][l1][l][l2]*delta1*delta1));
  a1=(1./(1+delta2*delta2)*
      (f1[0][l2]+
       f1[1][l2]*2*delta2+
       f1[2][l2]*delta2*delta2));
  return(b*a1*a3);
} 

double wiko_triple_corr(setup,
			new_spins,i1,i2,i3,i4,i5,i6,
			sigma,delta1,delta2,delta3,
			k1,k2,k3,k4,k5)
wiko_setup *setup;
int new_spins;
double i1,i2,i3,i4,i5,i6;
double sigma;
double delta1,delta2,delta3;
int k1,k2,k3,k4,k5;
{
  static double f3_1[3][MAX_LAMBDA_INDEX/2+1]
    [MAX_LAMBDA_INDEX/2+1]
    [MAX_LAMBDA_INDEX/2+1];
  static double f3_2[3][MAX_LAMBDA_INDEX/2+1]
    [MAX_LAMBDA_INDEX/2+1]
    [MAX_LAMBDA_INDEX/2+1];
  static double f1[3][MAX_LAMBDA_INDEX/2+1];
  double a1,a3_1,a3_2;
  double b;
  double result;
  
  int lam1,lam2,lam3,i;
  int l1,l2,l3,l4,l5;
  int l1a,l1b,l2a,l2b,l3a,l3b;

  l1a=max(1,(int)(fabs(i2-i1)));
  l1b=l1a+1;
  l2a=max(1,(int)(fabs(i4-i3)));
  l2b=l2a+1;
  l3a=max(1,(int)(fabs(i6-i5)));
  l3b=l2a+1;
  if (new_spins){
    /* calculate arrays of F-coefficients for the given spins and mults. */
    /* the variation of delta's and sigmas can then be calculated without */
    /* repeating this initialisation (set new_spins=0) */
    loop_over_lambdas(MAX_LAMBDA_INDEX,lam1,lam2,lam3,i){
      l1=lam1>>1;
      l2=lam2>>1;
      l3=lam3>>1;
      f3_1[0][l1][l2][l3]=0.;
      f3_1[1][l1][l2][l3]=0.;
      f3_1[2][l1][l2][l3]=0.;
      f3_2[0][l1][l2][l3]=0.;
      f3_2[1][l1][l2][l3]=0.;
      f3_2[2][l1][l2][l3]=0.;
      f1[0][l3]=0.;
      f1[1][l3]=0.;
      f1[2][l3]=0.;
    }
    loop_over_lambdas(l1b,lam1,lam2,lam3,i){
      l1=lam1>>1;
      l2=lam2>>1;
      l3=lam3>>1;
      f3_1[0][l1][l2][l3]=setup->f_function1(l1a,l1a,i1,i2,lam1,lam2,lam3);
      f3_1[1][l1][l2][l3]=setup->f_function1(l1a,l1b,i1,i2,lam1,lam2,lam3);
      f3_1[2][l1][l2][l3]=setup->f_function1(l1b,l1b,i1,i2,lam1,lam2,lam3);
    }
    loop_over_lambdas(l2b,lam1,lam2,lam3,i){
      l1=lam1>>1;
      l2=lam2>>1;
      l3=lam3>>1;
      f3_2[0][l1][l2][l3]=setup->f_function2(l2a,l2a,i3,i4,lam1,lam2,lam3);
      f3_2[1][l1][l2][l3]=setup->f_function2(l2a,l2b,i3,i4,lam1,lam2,lam3);
      f3_2[2][l1][l2][l3]=setup->f_function2(l2b,l2b,i3,i4,lam1,lam2,lam3);
    }
    for (lam3=0;lam3<=l3b;lam3+=2){
      l3=lam3>>1;
      f1[0][l3]=setup->f_function3(l3a,l3a,i5,i6,lam3,lam3,0);
      f1[1][l3]=setup->f_function3(l3a,l3b,i5,i6,lam3,lam3,0);
      f1[2][l3]=setup->f_function3(l3b,l3b,i5,i6,lam3,lam3,0);
    }
  }
  l1=k1>>1;
  l2=k2>>1;
  l3=k3>>1;
  l4=k4>>1;
  l5=k5>>1;
  b=setup->orient_function(new_spins,i1,sigma,k1);
  a3_1=(1./(1+delta1*delta1)*
	(f3_1[0][l1][l2][l3]+
	 f3_1[1][l1][l2][l3]*2*delta1+
	 f3_1[2][l1][l2][l3]*delta1*delta1));
  a3_2=(1./(1+delta2*delta2)*
	(f3_2[0][l3][l4][l5]+
	 f3_2[1][l3][l4][l5]*2*delta2+
	 f3_2[2][l3][l4][l5]*delta2*delta2));
  a1=(1./(1+delta3*delta3)*
      (f1[0][l5]+
       f1[1][l5]*2*delta3+
       f1[2][l5]*delta3*delta3));
  result=b*a1*a3_1*a3_2;
  return(result);
}
				      
double wiko_dis_from_oriented(new_spins,i1,i2,sigma,delta,k)
     int new_spins;
     double i1;double i2;
     double sigma;
     double delta;
     int k;
{
  double a,b;
  int la,lb;

  la=max(1,(int)(fabs(i2-i1)));
  lb=la+1;
  b=wiko_orient_gaussian(1,i1,sigma,k);
  a=(1./(1+delta*delta))*
    (wiko_f3(la,la,i1,i2,k,k,0.)+2*delta*wiko_f3(la,lb,i1,i2,k,k,0.)+
     delta*delta*wiko_f3(lb,lb,i1,i2,k,0.,0.));
  return(b*a);
}

double wiko_deo_fakt(i1,i2,l,lambda)
     double i1;double i2;int l;int lambda;
{
  return(minuseinshoch((int)(i1+i2+l+lambda))*
	 sqrt((2*i1+1)*(2*i2+1))*wign6j(i2,i2,lambda,i1,i1,l));
}

double wiko_deorient(new_spins,i1,i2,delta,lambda)
     int new_spins;
     double i1;double i2;
     double delta;
     int lambda;
{
  static double u[2][MAX_LAMBDA_INDEX/2+1];

  int lam,l1a,l1b,l;

  l1a=max(1,(int)(fabs(i2-i1)));
  l1b=l1a+1;
  if (new_spins){
    for (lam=0;lam<=MAX_LAMBDA_INDEX;lam+=2){
      l=(int)(lam/2);
      u[0][l]=(wiko_deo_fakt(i1,i2,l1a,lam));
      u[1][l]=(wiko_deo_fakt(i1,i2,l1b,lam));
    }
  }
  l=(int)(lambda/2);
  return((u[0][l]+delta*delta*u[1][l])/(1+delta*delta));
}

/**********************************************************/
/****** routines for alpha-gamma-gamma correlations *******/
/**********************************************************/

double gg_partial_angular_function(theta1,phi1,theta2,phi2,lam0,lam1,lam2,q0)
     double theta1,phi1,theta2,phi2;
     int lam0,lam1,lam2,q0;
{
  int q1,q2;
  double sum=0.;
  
  for (q2= -lam2;q2<=lam2;q2++){
    q1=-(q0+q2);
    if (abs(q1)<=lam1){
      sum+= 
	wign3j((double)lam2,(double)lam1,(double)lam0,
	       (double)q2,(double)q1,(double)q0) *
	spherical_harmonic(lam1,q1,theta1,0.)*
	spherical_harmonic(lam2,q2,theta2,0.)*
	(cos(q1*phi1*(M_PI/180.))*cos(q2*phi2*(M_PI/180.))
	 +sin(q1*phi1*(M_PI/180.))*sin(q2*phi2*(M_PI/180.)));
    }
  }
  sum *= (1./sqrt(2*lam2+1.));
  return(sum);
}

#define dsin(a) sin((a)*M_PI/180.)
#define dcos(a) cos((a)*M_PI/180.)


double partial_g_corr_old(new_spins,
			  l0,q0,
			  i2,i3,
			  theta1,phi1)
int new_spins,l0,q0;
double i2,i3,theta1,phi1;
{
  
  int lam0;
  int lg1;
  double sum=0;
  static double f1[4];

  lg1=(int)fabs(i2-i3);
  if (lg1==0 && i2+i3>0) lg1++;
  if (new_spins){
    for (lam0=0;lam0<=2*lg1;lam0+=2){
      f1[lam0/2]=wiko_f3(lg1,lg1,i2,i3,lam0,lam0,0);
    }
  }
  lam0=l0;
  // This normalization was changed in accordance with the Hamilton
  // "experimentally" verified through a 0-1-0 alpha-gamma cascade  
  sum =
    (1./sqrt(2*lam0+1))*
    f1[lam0/2]*
    spherical_harmonic(lam0,q0,theta1,0.)*
    cos(q0*phi1*M_PI/180.);
  return(sum);
}

double partial_g_corr(new_spins,
		      l0,q0,
		      i2,i3,
		      l1,l2,delta,
		      theta1,phi1)
     int new_spins,l0,q0,l1,l2;
double i2,i3,theta1,phi1,delta;
{
  
  int lam0;
  double sum=0;
  static double f1[4];

  //delta is multipolarity ratio l2/l1
  if (new_spins){
    for (lam0=0;lam0<=2*l2;lam0+=2){
      f1[lam0/2]=(wiko_f3(l1,l1,i2,i3,lam0,lam0,0) +
		  2*delta*wiko_f3(l1,l2,i2,i3,lam0,lam0,0) +
		  delta*delta*wiko_f3(l2,l2,i2,i3,lam0,lam0,0));
      f1[lam0/2] /= (1 + delta*delta);
    }
  }
  lam0=l0;
  // This normalization was changed in accordance with the Hamilton
  // "experimentally" verified through a 0-1-0 alpha-gamma cascade  
  sum =
    (1./sqrt(2*lam0+1))*
    f1[lam0/2]*
    spherical_harmonic(lam0,q0,theta1,0.)*
    cos(q0*phi1*M_PI/180.);
  return(sum);
}

double partial_p_corr(new_spins,
		      l0,q0,
		      i2,i3,
		      theta1,phi1)
int new_spins,l0,q0;
double i2,i3,theta1,phi1;
{
  
  int lam0;
  double jp1;
  double sum=0;
  static double f1[10];
  int i=0;

  //for proton transfer/decay, have to sum all possible
  //TOTAL spins allowed between initial target state and
  //residual state (or vice versa for decay) -- kgh, 20190815

  if (new_spins){
    for (i=0; i<10; i++) f1[i]=0;
    for (jp1=fabs(i2-i3); jp1<=(i2+i3); jp1++) { //sum over possible spins of decay particle
      for (lam0=0;lam0<=2*jp1;lam0+=2){ //sum over possible multipoles for a given decay spin
	f1[lam0/2]+=wiko_prot_f1(jp1,jp1,i2,i3,lam0);
      }
    }
  }
  lam0=l0;
  sum =
    (1./sqrt(2*lam0+1))*
    f1[lam0/2]*
    spherical_harmonic(lam0,q0,theta1,0.)*
    cos(q0*phi1*M_PI/180.);
  return(sum);
}

double partial_gg_corr(new_spins,
		       l0,q0,
		       i2,i3,i4,
		       theta1,phi1,
		       theta2,phi2)
int new_spins,l0,q0;
double i2,i3,i4,theta1,phi1,theta2,phi2;
{
  
  int lam0,lam1,lam2;
  int lg1,lg2;
  double sum=0;
  static double f3[8][4][4],f1[4];

  lg1=(int)fabs(i2-i3);
  lg2=(int)fabs(i3-i4);
  if (new_spins){
    /* buffer the f-coefficients */
    for (lam0=0;lam0<16;lam0+=2){
      for (lam1=0;lam1<=2*lg1;lam1+=2){
	for (lam2=0;lam2<=2*lg2;lam2+=2){
	  f3[lam0/2][lam1/2][lam2/2]=wiko_f3(lg1,lg1,i2,i3,lam0,lam1,lam2);
	}
      }
    }
    for (lam2=0;lam2<=2*lg2;lam2+=2){
      f1[lam2/2]=wiko_f3(lg2,lg2,i3,i4,lam2,lam2,0);
    }
  }
  lam0=l0;
  for (lam1=0;lam1<=2*lg1;lam1+=2){
    for (lam2=0;lam2<=2*lg2;lam2+=2){
      if ((abs(lam0-lam1)<=lam2)&&(lam0+lam1>=lam2)){
	sum +=
	  f3[lam0/2][lam1/2][lam2/2]*
	  f1[lam2/2]*
	  gg_partial_angular_function
	  (theta1,phi1,theta2,phi2,lam0,lam1,lam2,q0);
      }
    }
  }
  return(sum);
}

double partial_alpha_alpha_corr_clebsch(int new_spins,
					double *imaginary,
					double i2,double i3,
					int lmin,
					double delta,
					double theta2,
					double phi,
					int k,int q)
{
  int mba,mbb;
  int lba,lbb;
  int lb1,lb2;
  int i;
  complex double sum=0.,value;
  double mat;

  lb1=lmin;
  lb2=lb1+2;
  for (i=0;i<1;i++){
    switch(i){
    case 0:
      lba=lb1;lbb=lb1;mat=1.;break;
    case 1:
      lba=lb1;lbb=lb2;mat=2*delta;break;
    case 2:
      lba=lb2;lbb=lb2;mat=delta*delta;break;
    }
    for (mba=-lba;mba<=lba;mba++){
      mbb=mba-q;
      if ((abs(mbb)<=min(lbb,i3))&&(abs(mba)<=min(lba,i3))){
	value=minuseinshoch((int)(lbb-mbb));
	value *= clebsch((double)i3,(double)lba,(double)i2,
			 (double)mba,(double)-mba,0.);
	value *= clebsch((double)i3,(double)lbb,(double)i2,
			 (double)mbb,(double)-mbb,0.);
	value *= clebsch((double)i3,(double)i3,(double)k,
			 (double)mba,(double)-mbb,(double)q);
	value *= spherical_harmonic(lba,-mba,theta2,phi);
	value *= conj(spherical_harmonic(lbb,-mbb,theta2,phi));

	sum += mat* value;
      }
    }
  }
  sum *= (4.*M_PI)/(1.+delta*delta);
  *imaginary=cimag(sum);
  return(creal(sum));
}

double partial_alpha_alpha_corr(int new_spins,
				double *imaginary,
				double i1,double i2,
				int lmin,
				double delta,
				double theta,
				double phi,
				int k,int q)
{
  int mba,mbb;
  int lba,lbb;
  int lb1,lb2;
  int i;
  complex double sum=0.,value,value_ad;
  double mat;

  lb1=lmin;
  lb2=lb1+2;
  // Here we only use stretched coupling (i==0)
  for (i=0;i<1;i++){
    switch(i){
    case 0:
      lba=lb1;lbb=lb1;mat=1.;break;
    case 1:
      lba=lb1;lbb=lb2;mat=2*delta;break;
    case 2:
      lba=lb2;lbb=lb2;mat=delta*delta;break;
    }
    for (mba=-lba;mba<=lba;mba++){
      mbb= mba-q;
      if ((abs(mba)<=i2)&&(abs(mbb)<=i2)){
	value=minuseinshoch((int)(i2-mbb))*(2.*i1+1)*sqrt(2.*k+1.);
	value *= wign3j((double)i2,(double)lba,(double)i1,
			(double)-mba,(double)mba,0.);
	value *= wign3j((double)i2,(double)lbb,(double)i1,
			(double)-mbb,(double)mbb,0.);
	value *= wign3j((double)i2,(double)i2,(double)k,
			(double)-mba,(double)mbb,(double)q);
	value_ad = complex_spherical_harmonic(lba,-mba,theta,phi)
  	  * conj(complex_spherical_harmonic(lbb,-mbb,theta,phi));
	
	sum += value*value_ad;
      }
    }	
    
  }
  sum *= (4.*M_PI)/(1.+delta*delta);
  *imaginary=cimag(sum);
  return(creal(sum));
}

double partial_alpha_alpha_corr_lower(int new_spins,
				      double *imaginary,
				      double i2,double i3,
				      int lmin,
				      double delta,
				      double theta,
				      double phi,
				      int k,int q)
{
  int mba,mbb;
  int lba,lbb;
  int lb1,lb2;
  int i;
  complex double sum=0.,value,value_ad;
  double mat;

  lb1=lmin;
  lb2=lb1+2;
  // Here we only use stretched coupling (i==0)
  for (i=0;i<1;i++){
    switch(i){
    case 0:
      lba=lb1;lbb=lb1;mat=1.;break;
    case 1:
      lba=lb1;lbb=lb2;mat=2*delta;break;
    case 2:
      lba=lb2;lbb=lb2;mat=delta*delta;break;
    }
    for (mba=-lba;mba<=lba;mba++){
      mbb=mba-q;
      if ((abs(mbb)<=min(lbb,i2))&&(abs(mba)<=min(lba,i2))){
	value=minuseinshoch((int)(i2+mba))*(2.*i3+1)*sqrt(2.*k+1.);
	value *= wign3j((double)i3,(double)lba,(double)i2,
			0.,(double)-mba,(double)mba);
	value *= wign3j((double)i3,(double)lbb,(double)i2,
			0.,(double)-mbb,(double)mbb);
	value *= wign3j((double)i2,(double)i2,(double)k,
			(double)-mba,(double)mbb,(double)q);
	value_ad = complex_spherical_harmonic(lba,mba,theta,phi)
	  * conj(complex_spherical_harmonic(lbb,mbb,theta,phi));
	
	sum += value*value_ad;
      }
    }	

  }
  sum *= (4.*M_PI)/(1.+delta*delta);
  *imaginary=cimag(sum);
  return(creal(sum));
}



/************************************/
/* routines as interface to gsl */
/************************************/

double legen(int k, int q, double tet)
{
  double costet = cos(tet*M_PI/180.);
  return gsl_sf_legendre_Plm (k,q,costet);
}

double clebsch(i1,i2,i,m1,m2,m)
     double i1;double i;double i2;double m1;double m;double m2;
{
  double r;
  r = minuseinshoch(m+i1-i2)*sqrt(2*i+1)
    *gsl_sf_coupling_3j((int)(2*i1),(int)(2*i2),(int)(2*i),
			(int)(2*m1),(int)(2*m2),(int)(-2*m));
#if 0
  if (fabs(r)<1e-20)
    r=0.;
#endif
  return(r);
}

double wign3j(a,b,c,xx,yy,zz)
     double a; double b; double c; double xx; double yy; double zz;
{
  double r;
  r=gsl_sf_coupling_3j((int)(2*a),(int)(2*b),(int)(2*c),
		       (int)(2*xx),(int)(2*yy),(int)(2*zz));
#if 0
  if (fabs(r)<1e-20)
    r=0.;
#endif
  return(r);
}

double wign6j(a,b,c,xx,yy,zz)
     double a, b, c, xx, yy, zz;
{
  double r;
  r=gsl_sf_coupling_6j((int)(2*a),(int)(2*b),(int)(2*c),
		       (int)(2*xx),(int)(2*yy),(int)(2*zz));
#if 0
  if (fabs(r)<1e-20)
    r=0.;
#endif
  return(r);
}

double wign9j(a, b, c, xx, yy, zz, gg, hh, pp)
     double a, b, c, xx, yy, zz, gg, hh, pp;
{
  double r;
  r=gsl_sf_coupling_9j((int)(2*a),(int)(2*b),(int)(2*c),
		       (int)(2*xx),(int)(2*yy),(int)(2*zz),
		       (int)(2*gg),(int)(2*hh),(int)(2*pp));
#if 0
  if (fabs(r)<1e-20)
    r=0.;
#endif
  return(r);
}



