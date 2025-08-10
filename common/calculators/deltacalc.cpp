//expects as inputs: L2 B(L2) L1 B(L1)
//L's are multipoles (integers)
//B's are reduced transition rates in W.u.
//kgh 20191202

#include <iostream>
#include <cstdlib>
#include <cmath>

using namespace std;

int double_factorial(int n);
double Weisskopf_E(int L, int A);
double Weisskopf_M(int L, int A);
double Rate_E(double B, int L, double Egamma);
double Rate_M(double B, int L, double Egamma);

int main(int argc, char* argv[]) {

  int A = 25;
  double Egamma = 0.974742; //MeV

  int E2 = 1; //0 if magnetic, 1 if electric
  int L2 = 2; //multipolarity
  double B2 = 0.91; //red. trans. rate in W.u.

  int E1 = 0;
  int L1 = 1;
  double B1 = 0.00113;

  double rate2 = 0;
  double rate1 = 0;

  if (E2) {
    B2 *= Weisskopf_E(L2,A);
    rate2 = Rate_E(B2,L2,Egamma);
  }
  else {
    B2 *= Weisskopf_M(L2,A);
    rate2 = Rate_M(B2,L2,Egamma);
  }

  if (E1) {
    B1 *= Weisskopf_E(L1,A);
    rate1 = Rate_E(B1,L1,Egamma);
  }
  else {
    B1 *= Weisskopf_M(L1,A);
    rate1 = Rate_M(B1,L1,Egamma);
  }

  cout << B2 << " " << B1 << endl;
    //       << rate2 << " " << rate1 << " " << sqrt(rate2/rate1) << endl;

  return 0;

}

int double_factorial(int n) {

  if (n!=0 && n!=-1) return n*double_factorial(n-2);
  return 1;

}

double Weisskopf_E(int L, int A) {

  return pow(1.2,2*L)/(4*M_PI)*pow((3./(L+3)),2)*pow(A,2.*L/3); //units of e^2*fm^(2L)

}

double Weisskopf_M(int L, int A) {

  return 10.*pow(1.2,2*L-2)/M_PI*pow((3./(L+3)),2)*pow(A,(2.*L-2)/3); //units of (nuclear magneton)^2*(fm)^(2L-2)

}

double Rate_E(double B, int L, double Egamma) {

  double hbar = 6.582119569E-22; //Mev*s
  double c = 2.99792458E23; //fm/s

  return B*8*M_PI*pow(Egamma/hbar,2*L+1)*(L+1)/(L*pow(c,2*L)*pow(double_factorial(2*L+1),2));

}

double Rate_M(double B, int L, double Egamma) {

  double hbar = 6.582119569E-22; //Mev*s
  double c = 2.99792458E23; //fm/s

  double
    e = 1, //fund charge
    hbarc = 197.327, //MeV*fm
    mp = 938.272; //MeV
 
  double nuc_mag = e*hbarc*c/(2*mp); //e*fm

  return B*8*M_PI*pow(Egamma/hbar,2*L+1)*(L+1)/(L*pow(c,2*L-2)*pow(nuc_mag,2)*pow(double_factorial(2*L+1),2));

}
