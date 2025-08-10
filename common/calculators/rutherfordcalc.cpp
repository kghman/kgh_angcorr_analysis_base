// idea is to fix masses, etc. in code and then input angle through command line

#include <iostream>
#include <cmath>
#include <cstdlib>

using namespace std;

double Rutherford_XSec(double eta, double k, double theta);

int main(int argc, char** argv) {

  const double HBARC             = 197.327; //MeV*fm
  const double FSC               = 1./(137.036); //fine structure const.
  const double UTOMEV            = 931.4940954; //MeV per u
  const double RESTMASS_ELECTRON = 0.000548579909; //amu

  //target
  const int ZT =  6;
  const int AT = 12;
  const double MT = (12. - ZT*RESTMASS_ELECTRON)*UTOMEV;

  //projectile (beam)
  const int ZP = 2;
  const int AP = 3;
  const double MP = (3.016029 - ZP*RESTMASS_ELECTRON)*UTOMEV;
  const double Ebeam = 18; //lab MeV

  const double ECM = Ebeam*AT/(AT+AP); //MeV

  const double MU  = MT*MP/(MT + MP);
  const double K   = sqrt(2*MU*ECM)/HBARC;
  const double ETA = ZT*ZP*MU*FSC/(HBARC*K);

  if (argc!=2) {
    cout << "Enter scattering angle (in CM degrees)." << endl;
    return 1;
  }

  double theta = atof(argv[1])*M_PI/180.;

  cout << "For (" << ZP << "," << AP << ") + (" << ZT << "," << AT << ") @ " << Ebeam << " MeV (lab) -- RXsec(" << theta*180./M_PI << ") = " << 10.*Rutherford_XSec(ETA,K,theta) << " mb/sr" << endl;

  return 0;

}

double Rutherford_XSec(double eta, double k, double theta) {return eta*eta/(4*k*k*pow(sin(theta/2),4));}
