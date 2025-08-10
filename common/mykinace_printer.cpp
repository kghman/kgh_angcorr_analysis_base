#include <iostream>

#include "Reaction.h"

using namespace std;

int main(int argc, char** argv) {

  int TZ = 6, TA = 12;
  int PZ = 2, PA = 3;
  int EZ = 2, EA = 3;

  double beamE = 24.0;

  double Texc = 0;
  double Pexc = 0;
  double Eexc = 0;
  double Rexc = 0;

  double minang = 5*DEGTORAD;
  double maxang = 90*DEGTORAD;
  double angincr = 5*DEGTORAD;
  
  Reaction rxn(1, TZ,TA, PZ,PA, EZ,EA);
  rxn.Run_Kinematics(beamE);

  double ang = minang;

  cout << "Reaction (" << TZ << "," << TA << ") + (" << PZ << "," << PA << ") --> (" << EZ << "," << EA << ") + ("
       << TZ+PZ-EZ << "," << TA+PA-EA << ") at " << beamE << " MeV beam energy..." << endl
       << "--> Excitation energies: ExcT = " << Texc << " | ExcP = " << Pexc << " | ExcE = " << Eexc << " | ExcR = " << Rexc << " || Q-value = " << rxn.Q_Value() <<  " || Ecm = " << rxn.Entrance_Channel_KEcm() << " MeV" << endl;

  cout << "Lab Ang (deg.)   Lab E (MeV)   CM Ang (deg.)    LabToCm   Rutherford (mb/sr)" << endl;

  while (ang <= maxang) {
    rxn.Run_Kinematics(beamE, ang, 0, Texc, Pexc, Eexc, Rexc);
    cout << rxn.Get_theta_lab(3)*RADTODEG << "   " << rxn.Get_KE_lab(3) << "   " << rxn.Get_theta_cm(3)*RADTODEG << "   " << rxn.LabToCM_Factor(rxn.Get_theta_cm(3),beamE,Eexc,Rexc) << "   " << rxn.Rutherford_XSec(rxn.Get_theta_cm(3)) << endl;
    ang += angincr;
  }

  return 0;

}
