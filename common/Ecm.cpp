#include <iostream>
#include <cstdlib>
#include <cmath>

#include "Reaction.h"
#include "constants.h"

using namespace std;

int main(int argc, char* argv[]) {

  int ZT=atoi(argv[1]), AT=atoi(argv[2]),
      ZP=atoi(argv[4]),  AP=atoi(argv[5]);

  double TT = atof(argv[3]);
  double TP = atof(argv[6]);

  Reaction rxn(1, ZT,AT, ZP,AP);

  double MP = rxn.Projectile_Restmass();
  double MT = rxn.Target_Restmass();

  double EP = TP + MP;
  double ET = TT + MT;
  double PP = sqrt(TP*TP + 2*TP*MP);
  double PT = sqrt(TT*TT + 2*TT*MT);

  double Ecm = sqrt(pow(EP+ET,2) - pow(PP+PT,2));
  double Tcm = Ecm - (MP + MT);

  cout << "For P(" << ZP << "," << AP << ", " << TP << " MeV) + "
       <<     "T(" << ZT << "," << AT << ", " << TT << " MeV) --> "
       << "Tcm = " << Tcm << " MeV" << endl;

  return 0;

}
