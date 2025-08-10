#include "FrameConverter.h"

#include <iostream>
#include <fstream>
#include <cmath>

#include "constants.h"

using namespace std;

FrameConverter::FrameConverter(int ZP,int AP,
			       int ZT,int AT,
			       int ZE,int AE,
			       int ZD1,int AD1) {

  int ZR = ZP + ZT - ZE;
  int AR = AP + AT - AE;

  proj_lab.Z   = ZP; proj_lab.A   = AP;
  proj_rxncm.Z = ZP; proj_rxncm.A = AP;
  proj_deccm.Z = ZP; proj_deccm.A = AP;
  proj_brkcm.Z = ZP; proj_brkcm.A = AP;

  targ_lab.Z   = ZT; targ_lab.A   = AT;
  targ_rxncm.Z = ZT; targ_rxncm.A = AT;
  targ_deccm.Z = ZT; targ_deccm.A = AT;
  targ_brkcm.Z = ZT; targ_brkcm.A = AT;

  resi_lab.Z   = ZR; resi_lab.A   = AR;
  resi_rxncm.Z = ZR; resi_rxncm.A = AR;
  resi_deccm.Z = ZR; resi_deccm.A = AR;
  resi_brkcm.Z = ZR; resi_brkcm.A = AR;

  ejec_lab.Z   = ZE; ejec_lab.A   = AE;
  ejec_rxncm.Z = ZE; ejec_rxncm.A = AE;
  ejec_deccm.Z = ZE; ejec_deccm.A = AE;
  ejec_brkcm.Z = ZE; ejec_brkcm.A = AE;

  dec1_lab.Z   = ZD1; dec1_lab.A   = AD1;
  dec1_rxncm.Z = ZD1; dec1_rxncm.A = AD1;
  dec1_deccm.Z = ZD1; dec1_deccm.A = AD1;
  dec1_brkcm.Z = ZD1; dec1_brkcm.A = AD1;

  dec2_lab.Z   = ZR-ZD1; dec2_lab.A   = AR-AD1;
  dec2_rxncm.Z = ZR-ZD1; dec2_rxncm.A = AR-AD1;
  dec2_deccm.Z = ZR-ZD1; dec2_deccm.A = AR-AD1;
  dec2_brkcm.Z = ZR-ZD1; dec2_brkcm.A = AR-AD1;

  if (!GetMasses()) cout << "***FrameConverter() warning: problem loading rest masses\n";

  Q_rxn = proj_lab.M + targ_lab.M - resi_lab.M - ejec_lab.M;
  Q_dec = resi_lab.M - dec1_lab.M - dec2_lab.M;
  Q_brk = proj_lab.M - ejec_lab.M - dec1_lab.M;

  zero_structs();

  quiet_mode = false;

}

FrameConverter::~FrameConverter() {}

bool FrameConverter::Calculate(double beamE,
			       double ejec_theta,
			       double dec1_theta,
			       double dec1_phi,
			       double proj_exc,
			       double targ_exc,
			       double resi_exc,
			       double ejec_exc,
			       double dec1_exc,
			       double dec2_exc) { //inputs are MeV and radians

  //input quantities are also all supposed to be in the LAB frame (target at rest)

  zero_structs();

  proj_lab.Eexc = proj_exc; proj_rxncm.Eexc = proj_exc; proj_deccm.Eexc = proj_exc;
  targ_lab.Eexc = targ_exc; targ_rxncm.Eexc = targ_exc; targ_deccm.Eexc = targ_exc;
  resi_lab.Eexc = resi_exc; resi_rxncm.Eexc = resi_exc; resi_deccm.Eexc = resi_exc;
  ejec_lab.Eexc = ejec_exc; ejec_rxncm.Eexc = ejec_exc; ejec_deccm.Eexc = ejec_exc;
  dec1_lab.Eexc = dec1_exc; dec1_rxncm.Eexc = dec1_exc; dec1_deccm.Eexc = dec1_exc;
  dec2_lab.Eexc = dec2_exc; dec2_rxncm.Eexc = dec2_exc; dec2_deccm.Eexc = dec2_exc;

  proj_lab.Mtot   = proj_lab.M + proj_lab.Eexc;
  proj_rxncm.Mtot = proj_lab.Mtot;
  proj_deccm.Mtot = proj_lab.Mtot;

  targ_lab.Mtot   = targ_lab.M + targ_lab.Eexc;
  targ_rxncm.Mtot = targ_lab.Mtot;
  targ_deccm.Mtot = targ_lab.Mtot;

  resi_lab.Mtot   = resi_lab.M + resi_lab.Eexc;
  resi_rxncm.Mtot = resi_lab.Mtot;
  resi_deccm.Mtot = resi_lab.Mtot;

  ejec_lab.Mtot   = ejec_lab.M + ejec_lab.Eexc;
  ejec_rxncm.Mtot = ejec_lab.Mtot;
  ejec_deccm.Mtot = ejec_lab.Mtot;

  dec1_lab.Mtot   = dec1_lab.M + dec1_lab.Eexc;
  dec1_rxncm.Mtot = dec1_lab.Mtot;
  dec1_deccm.Mtot = dec1_lab.Mtot;

  dec2_lab.Mtot   = dec2_lab.M + dec2_lab.Eexc;
  dec2_rxncm.Mtot = dec2_lab.Mtot;
  dec2_deccm.Mtot = dec2_lab.Mtot;

  proj_lab.KE   = beamE;
  proj_lab.Etot = proj_lab.KE + proj_lab.Mtot;
  proj_lab.ptot = sqrt(pow(proj_lab.Etot,2) - pow(proj_lab.Mtot,2));
  proj_lab.pz   = proj_lab.ptot;

  ejec_lab.theta = ejec_theta;
  ejec_lab.phi   = M_PI;
  dec1_lab.theta = dec1_theta;
  dec1_lab.phi   = dec1_phi;

  targ_lab.Etot = targ_lab.Mtot;

  Q_rxn = proj_lab.Mtot + targ_lab.Mtot - resi_lab.Mtot - ejec_lab.Mtot;
  Q_dec = resi_lab.Mtot - dec1_lab.Mtot - dec2_lab.Mtot;

  proj_vec.SetPxPyPzE(proj_lab.px,
		      proj_lab.py,
		      proj_lab.pz,
		      proj_lab.Etot);

  targ_vec.SetPxPyPzE(targ_lab.px,
		      targ_lab.py,
		      targ_lab.pz,
		      targ_lab.Etot);

  rxn_parent_vec = proj_vec + targ_vec;
  rxn_boost_vec = rxn_parent_vec.BoostVector();

  ejec_lab.ptot = pE_lab();
  ejec_lab.Etot = sqrt(ejec_lab.ptot*ejec_lab.ptot + ejec_lab.Mtot*ejec_lab.Mtot);
  ejec_lab.px   = ejec_lab.ptot*sin(ejec_lab.theta)*cos(ejec_lab.phi);
  ejec_lab.py   = ejec_lab.ptot*sin(ejec_lab.theta)*sin(ejec_lab.phi);
  ejec_lab.pz   = ejec_lab.ptot*cos(ejec_lab.theta);

  ejec_vec.SetPxPyPzE(ejec_lab.px,
		      ejec_lab.py,
		      ejec_lab.pz,
		      ejec_lab.Etot);

  resi_vec = rxn_parent_vec - ejec_vec;
  dec_boost_vec = resi_vec.BoostVector();

  resi_lab.theta = resi_vec.Theta();
  resi_lab.ptot  = sqrt(pow(resi_vec.E(),2) - pow(resi_vec.M(),2));

  dec1_lab.ptot = pD1_lab();
  dec1_lab.Etot = sqrt(dec1_lab.ptot*dec1_lab.ptot + dec1_lab.Mtot*dec1_lab.Mtot);
  dec1_lab.px   = dec1_lab.ptot*sin(dec1_lab.theta)*cos(dec1_lab.phi);
  dec1_lab.py   = dec1_lab.ptot*sin(dec1_lab.theta)*sin(dec1_lab.phi);
  dec1_lab.pz   = dec1_lab.ptot*cos(dec1_lab.theta);

  dec1_vec.SetPxPyPzE(dec1_lab.px,
		      dec1_lab.py,
		      dec1_lab.pz,
		      dec1_lab.Etot);

  dec2_vec = resi_vec - dec1_vec;

  //store all proper values
  proj_lab.px    = proj_vec.Px();
  proj_lab.py    = proj_vec.Py();
  proj_lab.pz    = proj_vec.Pz();
  proj_lab.Etot  = proj_vec.E();
  proj_lab.ptot  = sqrt(pow(proj_lab.Etot,2) - pow(proj_lab.Mtot,2));
  proj_lab.KE    = proj_lab.Etot - proj_lab.Mtot;
  proj_lab.theta = proj_vec.Theta();
  proj_lab.phi   = proj_vec.Phi();

  targ_lab.px    = targ_vec.Px();
  targ_lab.py    = targ_vec.Py();
  targ_lab.pz    = targ_vec.Pz();
  targ_lab.Etot  = targ_vec.E();
  targ_lab.ptot  = sqrt(pow(targ_lab.Etot,2) - pow(targ_lab.Mtot,2));
  targ_lab.KE    = targ_lab.Etot - targ_lab.Mtot;
  targ_lab.theta = targ_vec.Theta();
  targ_lab.phi   = targ_vec.Phi();

  resi_lab.px    = resi_vec.Px();
  resi_lab.py    = resi_vec.Py();
  resi_lab.pz    = resi_vec.Pz();
  resi_lab.Etot  = resi_vec.E();
  resi_lab.ptot  = sqrt(pow(resi_lab.Etot,2) - pow(resi_lab.Mtot,2));
  resi_lab.KE    = resi_lab.Etot - resi_lab.Mtot;
  resi_lab.theta = resi_vec.Theta();
  resi_lab.phi   = resi_vec.Phi();

  ejec_lab.px    = ejec_vec.Px();
  ejec_lab.py    = ejec_vec.Py();
  ejec_lab.pz    = ejec_vec.Pz();
  ejec_lab.Etot  = ejec_vec.E();
  ejec_lab.ptot  = sqrt(pow(ejec_lab.Etot,2) - pow(ejec_lab.Mtot,2));
  ejec_lab.KE    = ejec_lab.Etot - ejec_lab.Mtot;
  ejec_lab.theta = ejec_vec.Theta();
  ejec_lab.phi   = ejec_vec.Phi();

  dec1_lab.px    = dec1_vec.Px();
  dec1_lab.py    = dec1_vec.Py();
  dec1_lab.pz    = dec1_vec.Pz();
  dec1_lab.Etot  = dec1_vec.E();
  dec1_lab.ptot  = sqrt(pow(dec1_lab.Etot,2) - pow(dec1_lab.Mtot,2));
  dec1_lab.KE    = dec1_lab.Etot - dec1_lab.Mtot;
  dec1_lab.theta = dec1_vec.Theta();
  dec1_lab.phi   = dec1_vec.Phi();

  dec2_lab.px    = dec2_vec.Px();
  dec2_lab.py    = dec2_vec.Py();
  dec2_lab.pz    = dec2_vec.Pz();
  dec2_lab.Etot  = dec2_vec.E();
  dec2_lab.ptot  = sqrt(pow(dec2_lab.Etot,2) - pow(dec2_lab.Mtot,2));
  dec2_lab.KE    = dec2_lab.Etot - dec2_lab.Mtot;
  dec2_lab.theta = dec2_vec.Theta();
  dec2_lab.phi   = dec2_vec.Phi();

  proj_vec.Boost(-rxn_boost_vec);
  targ_vec.Boost(-rxn_boost_vec);
  resi_vec.Boost(-rxn_boost_vec);
  ejec_vec.Boost(-rxn_boost_vec);
  dec1_vec.Boost(-rxn_boost_vec);
  dec2_vec.Boost(-rxn_boost_vec);

  proj_rxncm.px    = proj_vec.Px();
  proj_rxncm.py    = proj_vec.Py();
  proj_rxncm.pz    = proj_vec.Pz();
  proj_rxncm.Etot  = proj_vec.E();
  proj_rxncm.ptot  = sqrt(pow(proj_rxncm.Etot,2) - pow(proj_rxncm.Mtot,2));
  proj_rxncm.KE    = proj_rxncm.Etot - proj_rxncm.Mtot;
  proj_rxncm.theta = proj_vec.Theta();
  proj_rxncm.phi   = proj_vec.Phi();

  targ_rxncm.px    = targ_vec.Px();
  targ_rxncm.py    = targ_vec.Py();
  targ_rxncm.pz    = targ_vec.Pz();
  targ_rxncm.Etot  = targ_vec.E();
  targ_rxncm.ptot  = sqrt(pow(targ_rxncm.Etot,2) - pow(targ_rxncm.Mtot,2));
  targ_rxncm.KE    = targ_rxncm.Etot - targ_rxncm.Mtot;
  targ_rxncm.theta = targ_vec.Theta();
  targ_rxncm.phi   = targ_vec.Phi();

  resi_rxncm.px    = resi_vec.Px();
  resi_rxncm.py    = resi_vec.Py();
  resi_rxncm.pz    = resi_vec.Pz();
  resi_rxncm.Etot  = resi_vec.E();
  resi_rxncm.ptot  = sqrt(pow(resi_rxncm.Etot,2) - pow(resi_rxncm.Mtot,2));
  resi_rxncm.KE    = resi_rxncm.Etot - resi_rxncm.Mtot;
  resi_rxncm.theta = resi_vec.Theta();
  resi_rxncm.phi   = resi_vec.Phi();

  ejec_rxncm.px    = ejec_vec.Px();
  ejec_rxncm.py    = ejec_vec.Py();
  ejec_rxncm.pz    = ejec_vec.Pz();
  ejec_rxncm.Etot  = ejec_vec.E();
  ejec_rxncm.ptot  = sqrt(pow(ejec_rxncm.Etot,2) - pow(ejec_rxncm.Mtot,2));
  ejec_rxncm.KE    = ejec_rxncm.Etot - ejec_rxncm.Mtot;
  ejec_rxncm.theta = ejec_vec.Theta();
  ejec_rxncm.phi   = ejec_vec.Phi();

  dec1_rxncm.px    = dec1_vec.Px();
  dec1_rxncm.py    = dec1_vec.Py();
  dec1_rxncm.pz    = dec1_vec.Pz();
  dec1_rxncm.Etot  = dec1_vec.E();
  dec1_rxncm.ptot  = sqrt(pow(dec1_rxncm.Etot,2) - pow(dec1_rxncm.Mtot,2));
  dec1_rxncm.KE    = dec1_rxncm.Etot - dec1_rxncm.Mtot;
  dec1_rxncm.theta = dec1_vec.Theta();
  dec1_rxncm.phi   = dec1_vec.Phi();

  dec2_rxncm.px    = dec2_vec.Px();
  dec2_rxncm.py    = dec2_vec.Py();
  dec2_rxncm.pz    = dec2_vec.Pz();
  dec2_rxncm.Etot  = dec2_vec.E();
  dec2_rxncm.ptot  = sqrt(pow(dec2_rxncm.Etot,2) - pow(dec2_rxncm.Mtot,2));
  dec2_rxncm.KE    = dec2_rxncm.Etot - dec2_rxncm.Mtot;
  dec2_rxncm.theta = dec2_vec.Theta();
  dec2_rxncm.phi   = dec2_vec.Phi();

  proj_vec.Boost(rxn_boost_vec);
  targ_vec.Boost(rxn_boost_vec);
  resi_vec.Boost(rxn_boost_vec);
  ejec_vec.Boost(rxn_boost_vec);
  dec1_vec.Boost(rxn_boost_vec);
  dec2_vec.Boost(rxn_boost_vec);

  proj_vec.Boost(-dec_boost_vec);
  targ_vec.Boost(-dec_boost_vec);
  resi_vec.Boost(-dec_boost_vec);
  ejec_vec.Boost(-dec_boost_vec);
  dec1_vec.Boost(-dec_boost_vec);
  dec2_vec.Boost(-dec_boost_vec);

  proj_deccm.px    = proj_vec.Px();
  proj_deccm.py    = proj_vec.Py();
  proj_deccm.pz    = proj_vec.Pz();
  proj_deccm.Etot  = proj_vec.E();
  proj_deccm.ptot  = sqrt(pow(proj_deccm.Etot,2) - pow(proj_deccm.Mtot,2));
  proj_deccm.KE    = proj_deccm.Etot - proj_deccm.Mtot;
  proj_deccm.theta = proj_vec.Theta();
  proj_deccm.phi   = proj_vec.Phi();

  targ_deccm.px    = targ_vec.Px();
  targ_deccm.py    = targ_vec.Py();
  targ_deccm.pz    = targ_vec.Pz();
  targ_deccm.Etot  = targ_vec.E();
  targ_deccm.ptot  = sqrt(pow(targ_deccm.Etot,2) - pow(targ_deccm.Mtot,2));
  targ_deccm.KE    = targ_deccm.Etot - targ_deccm.Mtot;
  targ_deccm.theta = targ_vec.Theta();
  targ_deccm.phi   = targ_vec.Phi();

  resi_deccm.px    = resi_vec.Px();
  resi_deccm.py    = resi_vec.Py();
  resi_deccm.pz    = resi_vec.Pz();
  resi_deccm.Etot  = resi_vec.E();
  resi_deccm.ptot  = sqrt(pow(resi_deccm.Etot,2) - pow(resi_deccm.Mtot,2));
  resi_deccm.KE    = resi_deccm.Etot - resi_deccm.Mtot;
  resi_deccm.theta = resi_vec.Theta();
  resi_deccm.phi   = resi_vec.Phi();

  ejec_deccm.px    = ejec_vec.Px();
  ejec_deccm.py    = ejec_vec.Py();
  ejec_deccm.pz    = ejec_vec.Pz();
  ejec_deccm.Etot  = ejec_vec.E();
  ejec_deccm.ptot  = sqrt(pow(ejec_deccm.Etot,2) - pow(ejec_deccm.Mtot,2));
  ejec_deccm.KE    = ejec_deccm.Etot - ejec_deccm.Mtot;
  ejec_deccm.theta = ejec_vec.Theta();
  ejec_deccm.phi   = ejec_vec.Phi();

  dec1_deccm.px    = dec1_vec.Px();
  dec1_deccm.py    = dec1_vec.Py();
  dec1_deccm.pz    = dec1_vec.Pz();
  dec1_deccm.Etot  = dec1_vec.E();
  dec1_deccm.ptot  = sqrt(pow(dec1_deccm.Etot,2) - pow(dec1_deccm.Mtot,2));
  dec1_deccm.KE    = dec1_deccm.Etot - dec1_deccm.Mtot;
  dec1_deccm.theta = dec1_vec.Theta();
  dec1_deccm.phi   = dec1_vec.Phi();

  dec2_deccm.px    = dec2_vec.Px();
  dec2_deccm.py    = dec2_vec.Py();
  dec2_deccm.pz    = dec2_vec.Pz();
  dec2_deccm.Etot  = dec2_vec.E();
  dec2_deccm.ptot  = sqrt(pow(dec2_deccm.Etot,2) - pow(dec2_deccm.Mtot,2));
  dec2_deccm.KE    = dec2_deccm.Etot - dec2_deccm.Mtot;
  dec2_deccm.theta = dec2_vec.Theta();
  dec2_deccm.phi   = dec2_vec.Phi();

  proj_vec.Boost(dec_boost_vec);
  targ_vec.Boost(dec_boost_vec);
  resi_vec.Boost(dec_boost_vec);
  ejec_vec.Boost(dec_boost_vec);
  dec1_vec.Boost(dec_boost_vec);
  dec2_vec.Boost(dec_boost_vec);

}

bool FrameConverter::Calculate_Breakup(double beamE,
				       double ejec_theta,
				       double proj_exc,
				       double ejec_exc,
				       double dec1_exc) {

  zero_structs();

  proj_lab.Eexc = proj_exc; proj_brkcm.Eexc = proj_exc;
  ejec_lab.Eexc = ejec_exc; ejec_brkcm.Eexc = ejec_exc;
  dec1_lab.Eexc = dec1_exc; dec1_brkcm.Eexc = dec1_exc;

  proj_lab.Mtot   = proj_lab.M + proj_lab.Eexc;
  proj_brkcm.Mtot = proj_lab.Mtot;

  ejec_lab.Mtot   = ejec_lab.M + ejec_lab.Eexc;
  ejec_brkcm.Mtot = ejec_lab.Mtot;

  dec1_lab.Mtot   = dec1_lab.M + dec1_lab.Eexc;
  dec1_brkcm.Mtot = dec1_lab.Mtot;

  proj_lab.KE   = beamE;
  proj_lab.Etot = proj_lab.KE + proj_lab.Mtot;
  proj_lab.ptot = sqrt(pow(proj_lab.Etot,2) - pow(proj_lab.Mtot,2));
  proj_lab.pz   = proj_lab.ptot;

  ejec_lab.theta = ejec_theta;
  ejec_lab.phi   = M_PI;

  Q_brk = proj_lab.Mtot - ejec_lab.Mtot - dec1_lab.Mtot;

  proj_vec.SetPxPyPzE(proj_lab.px,
		      proj_lab.py,
		      proj_lab.pz,
		      proj_lab.Etot);

  brk_boost_vec = proj_vec.BoostVector();

  ejec_lab.ptot = pE_breakup_lab();
  ejec_lab.Etot = sqrt(ejec_lab.ptot*ejec_lab.ptot + ejec_lab.Mtot*ejec_lab.Mtot);
  ejec_lab.px   = ejec_lab.ptot*sin(ejec_lab.theta)*cos(ejec_lab.phi);
  ejec_lab.py   = ejec_lab.ptot*sin(ejec_lab.theta)*sin(ejec_lab.phi);
  ejec_lab.pz   = ejec_lab.ptot*cos(ejec_lab.theta);

  ejec_vec.SetPxPyPzE(ejec_lab.px,
		      ejec_lab.py,
		      ejec_lab.pz,
		      ejec_lab.Etot);

  dec1_vec = proj_vec - ejec_vec;

  proj_lab.px    = proj_vec.Px();
  proj_lab.py    = proj_vec.Py();
  proj_lab.pz    = proj_vec.Pz();
  proj_lab.Etot  = proj_vec.E();
  proj_lab.ptot  = sqrt(pow(proj_lab.Etot,2) - pow(proj_lab.Mtot,2));
  proj_lab.KE    = proj_lab.Etot - proj_lab.Mtot;
  proj_lab.theta = proj_vec.Theta();
  proj_lab.phi   = proj_vec.Phi();

  ejec_lab.px    = ejec_vec.Px();
  ejec_lab.py    = ejec_vec.Py();
  ejec_lab.pz    = ejec_vec.Pz();
  ejec_lab.Etot  = ejec_vec.E();
  ejec_lab.ptot  = sqrt(pow(ejec_lab.Etot,2) - pow(ejec_lab.Mtot,2));
  ejec_lab.KE    = ejec_lab.Etot - ejec_lab.Mtot;
  ejec_lab.theta = ejec_vec.Theta();
  ejec_lab.phi   = ejec_vec.Phi();

  dec1_lab.px    = dec1_vec.Px();
  dec1_lab.py    = dec1_vec.Py();
  dec1_lab.pz    = dec1_vec.Pz();
  dec1_lab.Etot  = dec1_vec.E();
  dec1_lab.ptot  = sqrt(pow(dec1_lab.Etot,2) - pow(dec1_lab.Mtot,2));
  dec1_lab.KE    = dec1_lab.Etot - dec1_lab.Mtot;
  dec1_lab.theta = dec1_vec.Theta();
  dec1_lab.phi   = dec1_vec.Phi();

  proj_vec.Boost(-brk_boost_vec);
  ejec_vec.Boost(-brk_boost_vec);
  dec1_vec.Boost(-brk_boost_vec);

  proj_brkcm.px    = proj_vec.Px();
  proj_brkcm.py    = proj_vec.Py();
  proj_brkcm.pz    = proj_vec.Pz();
  proj_brkcm.Etot  = proj_vec.E();
  proj_brkcm.ptot  = sqrt(pow(proj_brkcm.Etot,2) - pow(proj_brkcm.Mtot,2));
  proj_brkcm.KE    = proj_brkcm.Etot - proj_brkcm.Mtot;
  proj_brkcm.theta = proj_vec.Theta();
  proj_brkcm.phi   = proj_vec.Phi();

  ejec_brkcm.px    = ejec_vec.Px();
  ejec_brkcm.py    = ejec_vec.Py();
  ejec_brkcm.pz    = ejec_vec.Pz();
  ejec_brkcm.Etot  = ejec_vec.E();
  ejec_brkcm.ptot  = sqrt(pow(ejec_brkcm.Etot,2) - pow(ejec_brkcm.Mtot,2));
  ejec_brkcm.KE    = ejec_brkcm.Etot - ejec_brkcm.Mtot;
  ejec_brkcm.theta = ejec_vec.Theta();
  ejec_brkcm.phi   = ejec_vec.Phi();

  dec1_brkcm.px    = dec1_vec.Px();
  dec1_brkcm.py    = dec1_vec.Py();
  dec1_brkcm.pz    = dec1_vec.Pz();
  dec1_brkcm.Etot  = dec1_vec.E();
  dec1_brkcm.ptot  = sqrt(pow(dec1_brkcm.Etot,2) - pow(dec1_brkcm.Mtot,2));
  dec1_brkcm.KE    = dec1_brkcm.Etot - dec1_brkcm.Mtot;
  dec1_brkcm.theta = dec1_vec.Theta();
  dec1_brkcm.phi   = dec1_vec.Phi();
  
  proj_vec.Boost(brk_boost_vec);
  ejec_vec.Boost(brk_boost_vec);
  dec1_vec.Boost(brk_boost_vec);
  
}

bool FrameConverter::Calculate_DecCM(double beamE,
				     double ejec_theta,
				     double dec1_theta_decCM,
				     double dec1_phi_decCM,
				     double proj_exc,
				     double targ_exc,
				     double resi_exc,
				     double ejec_exc,
				     double dec1_exc,
				     double dec2_exc) { //inputs are MeV and radians

  //input quantities are also all supposed to be in the LAB frame (target at rest)
  //EXCEPT dec1 theta/phi which are in decay (residual) rest frame

  zero_structs();

  proj_lab.Eexc = proj_exc; proj_rxncm.Eexc = proj_exc; proj_deccm.Eexc = proj_exc;
  targ_lab.Eexc = targ_exc; targ_rxncm.Eexc = targ_exc; targ_deccm.Eexc = targ_exc;
  resi_lab.Eexc = resi_exc; resi_rxncm.Eexc = resi_exc; resi_deccm.Eexc = resi_exc;
  ejec_lab.Eexc = ejec_exc; ejec_rxncm.Eexc = ejec_exc; ejec_deccm.Eexc = ejec_exc;
  dec1_lab.Eexc = dec1_exc; dec1_rxncm.Eexc = dec1_exc; dec1_deccm.Eexc = dec1_exc;
  dec2_lab.Eexc = dec2_exc; dec2_rxncm.Eexc = dec2_exc; dec2_deccm.Eexc = dec2_exc;

  proj_lab.Mtot   = proj_lab.M + proj_lab.Eexc;
  proj_rxncm.Mtot = proj_lab.Mtot;
  proj_deccm.Mtot = proj_lab.Mtot;

  targ_lab.Mtot   = targ_lab.M + targ_lab.Eexc;
  targ_rxncm.Mtot = targ_lab.Mtot;
  targ_deccm.Mtot = targ_lab.Mtot;

  resi_lab.Mtot   = resi_lab.M + resi_lab.Eexc;
  resi_rxncm.Mtot = resi_lab.Mtot;
  resi_deccm.Mtot = resi_lab.Mtot;

  ejec_lab.Mtot   = ejec_lab.M + ejec_lab.Eexc;
  ejec_rxncm.Mtot = ejec_lab.Mtot;
  ejec_deccm.Mtot = ejec_lab.Mtot;

  dec1_lab.Mtot   = dec1_lab.M + dec1_lab.Eexc;
  dec1_rxncm.Mtot = dec1_lab.Mtot;
  dec1_deccm.Mtot = dec1_lab.Mtot;

  dec2_lab.Mtot   = dec2_lab.M + dec2_lab.Eexc;
  dec2_rxncm.Mtot = dec2_lab.Mtot;
  dec2_deccm.Mtot = dec2_lab.Mtot;

  proj_lab.KE   = beamE;
  proj_lab.Etot = proj_lab.KE + proj_lab.Mtot;
  proj_lab.ptot = sqrt(pow(proj_lab.Etot,2) - pow(proj_lab.Mtot,2));
  proj_lab.pz   = proj_lab.ptot;

  ejec_lab.theta = ejec_theta;
  ejec_lab.phi   = M_PI;
  //***chief difference:
  dec1_deccm.theta = dec1_theta_decCM;
  dec1_deccm.phi   = dec1_phi_decCM;

  targ_lab.Etot   = targ_lab.Mtot;

  Q_rxn = proj_lab.Mtot + targ_lab.Mtot - resi_lab.Mtot - ejec_lab.Mtot;
  Q_dec = resi_lab.Mtot - dec1_lab.Mtot - dec2_lab.Mtot;

  proj_vec.SetPxPyPzE(proj_lab.px,
		      proj_lab.py,
		      proj_lab.pz,
		      proj_lab.Etot);

  targ_vec.SetPxPyPzE(targ_lab.px,
		      targ_lab.py,
		      targ_lab.pz,
		      targ_lab.Etot);

  rxn_parent_vec = proj_vec + targ_vec;
  rxn_boost_vec = rxn_parent_vec.BoostVector();

  ejec_lab.ptot = pE_lab();
  ejec_lab.Etot = sqrt(ejec_lab.ptot*ejec_lab.ptot + ejec_lab.Mtot*ejec_lab.Mtot);
  ejec_lab.px   = ejec_lab.ptot*sin(ejec_lab.theta)*cos(ejec_lab.phi);
  ejec_lab.py   = ejec_lab.ptot*sin(ejec_lab.theta)*sin(ejec_lab.phi);
  ejec_lab.pz   = ejec_lab.ptot*cos(ejec_lab.theta);

  ejec_vec.SetPxPyPzE(ejec_lab.px,
		      ejec_lab.py,
		      ejec_lab.pz,
		      ejec_lab.Etot);

  resi_vec = rxn_parent_vec - ejec_vec;
  dec_boost_vec = resi_vec.BoostVector();

  resi_lab.theta = resi_vec.Theta();
  resi_lab.ptot  = sqrt(pow(resi_vec.E(),2) - pow(resi_vec.M(),2));

  //now recall that our dec1 input angles were in the decCM frame
  //***(currently using a NON-RELATIVISTIC result for dec1_deccm.ptot)
  dec1_deccm.ptot = sqrt(2*Q_dec/(1/dec1_lab.Mtot + 1/dec2_lab.Mtot));
  dec1_deccm.Etot = sqrt(dec1_deccm.ptot*dec1_deccm.ptot + dec1_deccm.Mtot*dec1_deccm.Mtot);
  dec1_deccm.px   = dec1_deccm.ptot*sin(dec1_deccm.theta)*cos(dec1_deccm.phi);
  dec1_deccm.py   = dec1_deccm.ptot*sin(dec1_deccm.theta)*sin(dec1_deccm.phi);
  dec1_deccm.pz   = dec1_deccm.ptot*cos(dec1_deccm.theta);

  dec1_vec.SetPxPyPzE(dec1_deccm.px,
		      dec1_deccm.py,
		      dec1_deccm.pz,
		      dec1_deccm.Etot);

  //then boost back to lab frame before filling all values
  dec1_vec.Boost(dec_boost_vec);

  dec2_vec = resi_vec - dec1_vec;

  //store all proper values
  proj_lab.px    = proj_vec.Px();
  proj_lab.py    = proj_vec.Py();
  proj_lab.pz    = proj_vec.Pz();
  proj_lab.Etot  = proj_vec.E();
  proj_lab.ptot  = sqrt(pow(proj_lab.Etot,2) - pow(proj_lab.Mtot,2));
  proj_lab.KE    = proj_lab.Etot - proj_lab.Mtot;
  proj_lab.theta = proj_vec.Theta();
  proj_lab.phi   = proj_vec.Phi();

  targ_lab.px    = targ_vec.Px();
  targ_lab.py    = targ_vec.Py();
  targ_lab.pz    = targ_vec.Pz();
  targ_lab.Etot  = targ_vec.E();
  targ_lab.ptot  = sqrt(pow(targ_lab.Etot,2) - pow(targ_lab.Mtot,2));
  targ_lab.KE    = targ_lab.Etot - targ_lab.Mtot;
  targ_lab.theta = targ_vec.Theta();
  targ_lab.phi   = targ_vec.Phi();

  resi_lab.px    = resi_vec.Px();
  resi_lab.py    = resi_vec.Py();
  resi_lab.pz    = resi_vec.Pz();
  resi_lab.Etot  = resi_vec.E();
  resi_lab.ptot  = sqrt(pow(resi_lab.Etot,2) - pow(resi_lab.Mtot,2));
  resi_lab.KE    = resi_lab.Etot - resi_lab.Mtot;
  resi_lab.theta = resi_vec.Theta();
  resi_lab.phi   = resi_vec.Phi();

  ejec_lab.px    = ejec_vec.Px();
  ejec_lab.py    = ejec_vec.Py();
  ejec_lab.pz    = ejec_vec.Pz();
  ejec_lab.Etot  = ejec_vec.E();
  ejec_lab.ptot  = sqrt(pow(ejec_lab.Etot,2) - pow(ejec_lab.Mtot,2));
  ejec_lab.KE    = ejec_lab.Etot - ejec_lab.Mtot;
  ejec_lab.theta = ejec_vec.Theta();
  ejec_lab.phi   = ejec_vec.Phi();

  dec1_lab.px    = dec1_vec.Px();
  dec1_lab.py    = dec1_vec.Py();
  dec1_lab.pz    = dec1_vec.Pz();
  dec1_lab.Etot  = dec1_vec.E();
  dec1_lab.ptot  = sqrt(pow(dec1_lab.Etot,2) - pow(dec1_lab.Mtot,2));
  dec1_lab.KE    = dec1_lab.Etot - dec1_lab.Mtot;
  dec1_lab.theta = dec1_vec.Theta();
  dec1_lab.phi   = dec1_vec.Phi();

  dec2_lab.px    = dec2_vec.Px();
  dec2_lab.py    = dec2_vec.Py();
  dec2_lab.pz    = dec2_vec.Pz();
  dec2_lab.Etot  = dec2_vec.E();
  dec2_lab.ptot  = sqrt(pow(dec2_lab.Etot,2) - pow(dec2_lab.Mtot,2));
  dec2_lab.KE    = dec2_lab.Etot - dec2_lab.Mtot;
  dec2_lab.theta = dec2_vec.Theta();
  dec2_lab.phi   = dec2_vec.Phi();

  proj_vec.Boost(-rxn_boost_vec);
  targ_vec.Boost(-rxn_boost_vec);
  resi_vec.Boost(-rxn_boost_vec);
  ejec_vec.Boost(-rxn_boost_vec);
  dec1_vec.Boost(-rxn_boost_vec);
  dec2_vec.Boost(-rxn_boost_vec);

  proj_rxncm.px    = proj_vec.Px();
  proj_rxncm.py    = proj_vec.Py();
  proj_rxncm.pz    = proj_vec.Pz();
  proj_rxncm.Etot  = proj_vec.E();
  proj_rxncm.ptot  = sqrt(pow(proj_rxncm.Etot,2) - pow(proj_rxncm.Mtot,2));
  proj_rxncm.KE    = proj_rxncm.Etot - proj_rxncm.Mtot;
  proj_rxncm.theta = proj_vec.Theta();
  proj_rxncm.phi   = proj_vec.Phi();

  targ_rxncm.px    = targ_vec.Px();
  targ_rxncm.py    = targ_vec.Py();
  targ_rxncm.pz    = targ_vec.Pz();
  targ_rxncm.Etot  = targ_vec.E();
  targ_rxncm.ptot  = sqrt(pow(targ_rxncm.Etot,2) - pow(targ_rxncm.Mtot,2));
  targ_rxncm.KE    = targ_rxncm.Etot - targ_rxncm.Mtot;
  targ_rxncm.theta = targ_vec.Theta();
  targ_rxncm.phi   = targ_vec.Phi();

  resi_rxncm.px    = resi_vec.Px();
  resi_rxncm.py    = resi_vec.Py();
  resi_rxncm.pz    = resi_vec.Pz();
  resi_rxncm.Etot  = resi_vec.E();
  resi_rxncm.ptot  = sqrt(pow(resi_rxncm.Etot,2) - pow(resi_rxncm.Mtot,2));
  resi_rxncm.KE    = resi_rxncm.Etot - resi_rxncm.Mtot;
  resi_rxncm.theta = resi_vec.Theta();
  resi_rxncm.phi   = resi_vec.Phi();

  ejec_rxncm.px    = ejec_vec.Px();
  ejec_rxncm.py    = ejec_vec.Py();
  ejec_rxncm.pz    = ejec_vec.Pz();
  ejec_rxncm.Etot  = ejec_vec.E();
  ejec_rxncm.ptot  = sqrt(pow(ejec_rxncm.Etot,2) - pow(ejec_rxncm.Mtot,2));
  ejec_rxncm.KE    = ejec_rxncm.Etot - ejec_rxncm.Mtot;
  ejec_rxncm.theta = ejec_vec.Theta();
  ejec_rxncm.phi   = ejec_vec.Phi();

  dec1_rxncm.px    = dec1_vec.Px();
  dec1_rxncm.py    = dec1_vec.Py();
  dec1_rxncm.pz    = dec1_vec.Pz();
  dec1_rxncm.Etot  = dec1_vec.E();
  dec1_rxncm.ptot  = sqrt(pow(dec1_rxncm.Etot,2) - pow(dec1_rxncm.Mtot,2));
  dec1_rxncm.KE    = dec1_rxncm.Etot - dec1_rxncm.Mtot;
  dec1_rxncm.theta = dec1_vec.Theta();
  dec1_rxncm.phi   = dec1_vec.Phi();

  dec2_rxncm.px    = dec2_vec.Px();
  dec2_rxncm.py    = dec2_vec.Py();
  dec2_rxncm.pz    = dec2_vec.Pz();
  dec2_rxncm.Etot  = dec2_vec.E();
  dec2_rxncm.ptot  = sqrt(pow(dec2_rxncm.Etot,2) - pow(dec2_rxncm.Mtot,2));
  dec2_rxncm.KE    = dec2_rxncm.Etot - dec2_rxncm.Mtot;
  dec2_rxncm.theta = dec2_vec.Theta();
  dec2_rxncm.phi   = dec2_vec.Phi();

  proj_vec.Boost(rxn_boost_vec);
  targ_vec.Boost(rxn_boost_vec);
  resi_vec.Boost(rxn_boost_vec);
  ejec_vec.Boost(rxn_boost_vec);
  dec1_vec.Boost(rxn_boost_vec);
  dec2_vec.Boost(rxn_boost_vec);

  proj_vec.Boost(-dec_boost_vec);
  targ_vec.Boost(-dec_boost_vec);
  resi_vec.Boost(-dec_boost_vec);
  ejec_vec.Boost(-dec_boost_vec);
  dec1_vec.Boost(-dec_boost_vec);
  dec2_vec.Boost(-dec_boost_vec);

  proj_deccm.px    = proj_vec.Px();
  proj_deccm.py    = proj_vec.Py();
  proj_deccm.pz    = proj_vec.Pz();
  proj_deccm.Etot  = proj_vec.E();
  proj_deccm.ptot  = sqrt(pow(proj_deccm.Etot,2) - pow(proj_deccm.Mtot,2));
  proj_deccm.KE    = proj_deccm.Etot - proj_deccm.Mtot;
  proj_deccm.theta = proj_vec.Theta();
  proj_deccm.phi   = proj_vec.Phi();

  targ_deccm.px    = targ_vec.Px();
  targ_deccm.py    = targ_vec.Py();
  targ_deccm.pz    = targ_vec.Pz();
  targ_deccm.Etot  = targ_vec.E();
  targ_deccm.ptot  = sqrt(pow(targ_deccm.Etot,2) - pow(targ_deccm.Mtot,2));
  targ_deccm.KE    = targ_deccm.Etot - targ_deccm.Mtot;
  targ_deccm.theta = targ_vec.Theta();
  targ_deccm.phi   = targ_vec.Phi();

  resi_deccm.px    = resi_vec.Px();
  resi_deccm.py    = resi_vec.Py();
  resi_deccm.pz    = resi_vec.Pz();
  resi_deccm.Etot  = resi_vec.E();
  resi_deccm.ptot  = sqrt(pow(resi_deccm.Etot,2) - pow(resi_deccm.Mtot,2));
  resi_deccm.KE    = resi_deccm.Etot - resi_deccm.Mtot;
  resi_deccm.theta = resi_vec.Theta();
  resi_deccm.phi   = resi_vec.Phi();

  ejec_deccm.px    = ejec_vec.Px();
  ejec_deccm.py    = ejec_vec.Py();
  ejec_deccm.pz    = ejec_vec.Pz();
  ejec_deccm.Etot  = ejec_vec.E();
  ejec_deccm.ptot  = sqrt(pow(ejec_deccm.Etot,2) - pow(ejec_deccm.Mtot,2));
  ejec_deccm.KE    = ejec_deccm.Etot - ejec_deccm.Mtot;
  ejec_deccm.theta = ejec_vec.Theta();
  ejec_deccm.phi   = ejec_vec.Phi();

  dec1_deccm.px    = dec1_vec.Px();
  dec1_deccm.py    = dec1_vec.Py();
  dec1_deccm.pz    = dec1_vec.Pz();
  dec1_deccm.Etot  = dec1_vec.E();
  dec1_deccm.ptot  = sqrt(pow(dec1_deccm.Etot,2) - pow(dec1_deccm.Mtot,2));
  dec1_deccm.KE    = dec1_deccm.Etot - dec1_deccm.Mtot;
  dec1_deccm.theta = dec1_vec.Theta();
  dec1_deccm.phi   = dec1_vec.Phi();

  dec2_deccm.px    = dec2_vec.Px();
  dec2_deccm.py    = dec2_vec.Py();
  dec2_deccm.pz    = dec2_vec.Pz();
  dec2_deccm.Etot  = dec2_vec.E();
  dec2_deccm.ptot  = sqrt(pow(dec2_deccm.Etot,2) - pow(dec2_deccm.Mtot,2));
  dec2_deccm.KE    = dec2_deccm.Etot - dec2_deccm.Mtot;
  dec2_deccm.theta = dec2_vec.Theta();
  dec2_deccm.phi   = dec2_vec.Phi();

  proj_vec.Boost(dec_boost_vec);
  targ_vec.Boost(dec_boost_vec);
  resi_vec.Boost(dec_boost_vec);
  ejec_vec.Boost(dec_boost_vec);
  dec1_vec.Boost(dec_boost_vec);
  dec2_vec.Boost(dec_boost_vec);

}

double FrameConverter::LabToRxnCM_Ejec_Angle(double angle, double Beam_E, double Exc_Proj, double Exc_Targ, double Exc_Ejec, double Exc_Resi) {

  //variables for easier readibility
  double
    TRP = proj_lab.M + Exc_Proj,
    TRT = targ_lab.M + Exc_Proj,
    TRE = ejec_lab.M + Exc_Ejec,
    TRR = resi_lab.M + Exc_Resi,
    Q   = TRP + TRT - TRE - TRR,
    E   = Beam_E;

  //this one I really had to twist Illiadis around...he gives the formula for CM-->lab conversion,
  //which isn't so easily invertible. Trust me, this works, even if it isn't exactly obvious.
  double vcm = (TRP)/(TRP+TRT)*sqrt(2*E/TRP);
  double vbp = sqrt(2*TRR/(TRE*(TRE+TRR))*(Q + E*TRT/(TRP+TRT)));
  double velocity_ratio = vcm/vbp;
  velocity_ratio = 1/(cos(angle) + sqrt(1/(velocity_ratio*velocity_ratio) - sin(angle)*sin(angle)));

  return atan2(sin(angle),cos(angle) - velocity_ratio);

}

double FrameConverter::RxnCMToLab_Ejec_Angle(double angle, double Beam_E, double Exc_Proj, double Exc_Targ, double Exc_Ejec, double Exc_Resi) {

  double
    TRP = proj_lab.M + Exc_Proj,
    TRT = targ_lab.M + Exc_Targ,
    TRE = ejec_lab.M + Exc_Ejec,
    TRR = resi_lab.M + Exc_Resi,
    Q   = TRP + TRT - TRE - TRR,
    E   = Beam_E;

  double vcm = (TRP)/(TRP+TRT)*sqrt(2*E/TRP);
  double vbp = sqrt(2*TRR/(TRE*(TRE+TRR))*(Q + E*TRT/(TRP+TRT)));
  double velocity_ratio = vcm/vbp;

  return atan2(sin(angle),cos(angle) + velocity_ratio);
  
}

double FrameConverter::LabToRxnCM_Ejec_Angle_Illiadis(double angle, double Beam_E, double Exc_Proj, double Exc_Targ, double Exc_Ejec, double Exc_Resi) {

  //variables for easier readibility
  double
    TRP = proj_lab.M + Exc_Proj,
    TRT = targ_lab.M + Exc_Targ,
    TRE = ejec_lab.M + Exc_Ejec,
    TRR = resi_lab.M + Exc_Resi,
    Q   = TRP + TRT - TRE - TRR,
    E   = Beam_E;

  //this one I really had to twist Illiadis around...he gives the formula for CM-->lab conversion,
  //which isn't so easily invertible. Trust me, this works, even if it isn't exactly obvious.
  double velocity_ratio = sqrt(TRP*TRE*E/((TRR*(TRR + TRE))*Q + TRR*(TRR + TRE - TRP)*E));
  velocity_ratio = 1/(cos(angle) + sqrt(1/(velocity_ratio*velocity_ratio) - sin(angle)*sin(angle)));

  return atan2(sin(angle),cos(angle) - velocity_ratio);

}

double FrameConverter::RxnCMToLab_Ejec_Angle_Illiadis(double angle, double Beam_E, double Exc_Proj, double Exc_Targ, double Exc_Ejec, double Exc_Resi) {

  double
    TRP = proj_lab.M + Exc_Proj,
    TRT = targ_lab.M + Exc_Proj,
    TRE = ejec_lab.M + Exc_Ejec,
    TRR = resi_lab.M + Exc_Resi,
    Q   = TRP + TRT - TRE - TRR,
    E   = Beam_E;

  double velocity_ratio = sqrt(TRP*TRE*E/((TRR*(TRR + TRE))*Q + TRR*(TRR + TRE - TRP)*E));

  return atan2(sin(angle),cos(angle) + velocity_ratio);
  
}

double FrameConverter::LabToRxnCM_Ejec_Factor() {

  double
    RP  = proj_lab.Mtot,
    TRE = ejec_lab.Mtot,
    TRR = resi_lab.Mtot,
    Q   = Q_rxn,
    E   = proj_lab.KE;

  double velocity_ratio = sqrt(RP*TRE*E/((TRR*(TRR + TRE))*Q + TRR*(TRR + TRE - RP)*E));

  double numerator = 1 + velocity_ratio*cos(ejec_rxncm.theta);
  double denominator = pow(sin(ejec_rxncm.theta),2) + pow(cos(ejec_rxncm.theta) + velocity_ratio,2);
  denominator = pow(denominator,1.5);
  
  return numerator/denominator;

}

double FrameConverter::LabToRxnCM_Ejec_Factor_Incr() {

  double returnval = 0;

  double incr_val = 1E-3;

  double fx=0, fxp=0, hx=0, hxp=0;

  fx = cos(ejec_rxncm.theta);
  hx = cos(ejec_lab.theta);

  ejec_lab.theta += incr_val;

  Calculate(proj_lab.KE,
	    ejec_lab.theta,
	    dec1_lab.theta,
	    dec1_lab.phi,
	    proj_lab.Eexc,
	    targ_lab.Eexc,
	    resi_lab.Eexc,
	    ejec_lab.Eexc,
	    dec1_lab.Eexc,
	    dec2_lab.Eexc);

  fxp = cos(ejec_rxncm.theta);
  hxp = cos(ejec_lab.theta);

  ejec_lab.theta -= incr_val; //reset

  Calculate(proj_lab.KE,
	    ejec_lab.theta,
	    dec1_lab.theta,
	    dec1_lab.phi,
	    proj_lab.Eexc,
	    targ_lab.Eexc,
	    resi_lab.Eexc,
	    ejec_lab.Eexc,
	    dec1_lab.Eexc,
	    dec2_lab.Eexc);

  returnval += (hxp - hx)/(fxp - fx);

  return returnval;

}

double FrameConverter::LabToRxnCM_Resi_Factor_Incr() {

  double returnval = 0;

  double incr_val = 1E-3;

  double fx=0, fxp=0, hx=0, hxp=0;

  fx = cos(resi_rxncm.theta);
  hx = cos(resi_lab.theta);

  //since input is ejectile info, have to calculate corresponding angle
  resi_lab.theta += incr_val;
  resi_vec.SetTheta(resi_lab.theta);
  ejec_vec = rxn_parent_vec - resi_vec;
  ejec_lab.theta = ejec_vec.Theta();

  Calculate(proj_lab.KE,
	    ejec_lab.theta,
	    dec1_lab.theta,
	    dec1_lab.phi,
	    proj_lab.Eexc,
	    targ_lab.Eexc,
	    resi_lab.Eexc,
	    ejec_lab.Eexc,
	    dec1_lab.Eexc,
	    dec2_lab.Eexc);

  fxp = cos(resi_rxncm.theta);
  hxp = cos(resi_lab.theta);

  resi_lab.theta -= incr_val; //reset
  resi_vec.SetTheta(resi_lab.theta);
  ejec_vec = rxn_parent_vec - resi_vec;
  ejec_lab.theta = ejec_vec.Theta();

  Calculate(proj_lab.KE,
	    ejec_lab.theta,
	    dec1_lab.theta,
	    dec1_lab.phi,
	    proj_lab.Eexc,
	    targ_lab.Eexc,
	    resi_lab.Eexc,
	    ejec_lab.Eexc,
	    dec1_lab.Eexc,
	    dec2_lab.Eexc);

  returnval += (hxp - hx)/(fxp - fx);
    
  return returnval;

}

double FrameConverter::LabToDecCM_Resi_Factor_Incr() {

  double returnval = 0;

  double incr_val = 1E-3;

  double fx=0, fxp=0, hx=0, hxp=0;

  fx = cos(resi_deccm.theta);
  hx = cos(resi_lab.theta);

  resi_lab.theta += incr_val;
  resi_vec.SetTheta(resi_lab.theta);
  ejec_vec = rxn_parent_vec - resi_vec;
  ejec_lab.theta = ejec_vec.Theta();

  Calculate(proj_lab.KE,
	    ejec_lab.theta,
	    dec1_lab.theta,
	    dec1_lab.phi,
	    proj_lab.Eexc,
	    targ_lab.Eexc,
	    resi_lab.Eexc,
	    ejec_lab.Eexc,
	    dec1_lab.Eexc,
	    dec2_lab.Eexc);

  fxp = cos(resi_deccm.theta);
  hxp = cos(resi_lab.theta);

  resi_lab.theta -= incr_val; //reset
  resi_vec.SetTheta(resi_lab.theta);
  ejec_vec = rxn_parent_vec - resi_vec;
  ejec_lab.theta = ejec_vec.Theta();

  Calculate(proj_lab.KE,
	    ejec_lab.theta,
	    dec1_lab.theta,
	    dec1_lab.phi,
	    proj_lab.Eexc,
	    targ_lab.Eexc,
	    resi_lab.Eexc,
	    ejec_lab.Eexc,
	    dec1_lab.Eexc,
	    dec2_lab.Eexc);

  returnval += (hxp - hx)/(fxp - fx);
    
  return returnval;

}

double FrameConverter::LabToRxnCM_Dec1_Factor_Incr() {

  double returnval = 0;

  double incr_val = 1E-3;

  double fx=0, fxp=0, hx=0, hxp=0;

  fx = cos(dec1_rxncm.theta);
  hx = cos(dec1_lab.theta);

  dec1_lab.theta += incr_val;

  Calculate(proj_lab.KE,
	    ejec_lab.theta,
	    dec1_lab.theta,
	    dec1_lab.phi,
	    proj_lab.Eexc,
	    targ_lab.Eexc,
	    resi_lab.Eexc,
	    ejec_lab.Eexc,
	    dec1_lab.Eexc,
	    dec2_lab.Eexc);

  fxp = cos(dec1_rxncm.theta);
  hxp = cos(dec1_lab.theta);

  dec1_lab.theta -= incr_val; //reset

  Calculate(proj_lab.KE,
	    ejec_lab.theta,
	    dec1_lab.theta,
	    dec1_lab.phi,
	    proj_lab.Eexc,
	    targ_lab.Eexc,
	    resi_lab.Eexc,
	    ejec_lab.Eexc,
	    dec1_lab.Eexc,
	    dec2_lab.Eexc);

  returnval += (hxp - hx)/(fxp - fx);
    
  return returnval;

}

double FrameConverter::LabToDecCM_Dec1_Factor_Incr() {

  double returnval = 1;

  double incr_val = 1E-3;

  double fx=0, fxp=0, hx=0, hxp=0;

  fx = cos(dec1_deccm.theta);
  hx = cos(dec1_lab.theta);
  dec1_lab.theta += incr_val;
  Calculate(proj_lab.KE,
	    ejec_lab.theta,
	    dec1_lab.theta,
	    dec1_lab.phi,
	    proj_lab.Eexc,
	    targ_lab.Eexc,
	    resi_lab.Eexc,
	    ejec_lab.Eexc,
	    dec1_lab.Eexc,
	    dec2_lab.Eexc);
  fxp = cos(dec1_deccm.theta);
  hxp = cos(dec1_lab.theta);
  dec1_lab.theta -= incr_val; //reset
  Calculate(proj_lab.KE,
	    ejec_lab.theta,
	    dec1_lab.theta,
	    dec1_lab.phi,
	    proj_lab.Eexc,
	    targ_lab.Eexc,
	    resi_lab.Eexc,
	    ejec_lab.Eexc,
	    dec1_lab.Eexc,
	    dec2_lab.Eexc);
  returnval *= (fxp - fx)/(hxp - hx);
  //now for phi part

  fx = dec1_deccm.phi;
  hx = dec1_lab.phi;
  dec1_lab.phi += incr_val;
  Calculate(proj_lab.KE,
	    ejec_lab.theta,
	    dec1_lab.theta,
	    dec1_lab.phi,
	    proj_lab.Eexc,
	    targ_lab.Eexc,
	    resi_lab.Eexc,
	    ejec_lab.Eexc,
	    dec1_lab.Eexc,
	    dec2_lab.Eexc);
  fxp = dec1_deccm.phi;
  hxp = dec1_lab.phi;
  dec1_lab.phi -= incr_val; //reset
  Calculate(proj_lab.KE,
	    ejec_lab.theta,
	    dec1_lab.theta,
	    dec1_lab.phi,
	    proj_lab.Eexc,
	    targ_lab.Eexc,
	    resi_lab.Eexc,
	    ejec_lab.Eexc,
	    dec1_lab.Eexc,
	    dec2_lab.Eexc);
  returnval *= (fxp - fx)/(hxp - hx);
    
  //second term...
  double second_term = 1;

  fx = cos(dec1_deccm.theta);
  hx = dec1_lab.phi;
  dec1_lab.phi += incr_val;
  Calculate(proj_lab.KE,
	    ejec_lab.theta,
	    dec1_lab.theta,
	    dec1_lab.phi,
	    proj_lab.Eexc,
	    targ_lab.Eexc,
	    resi_lab.Eexc,
	    ejec_lab.Eexc,
	    dec1_lab.Eexc,
	    dec2_lab.Eexc);
  fxp = cos(dec1_deccm.theta);
  hxp = dec1_lab.phi;
  dec1_lab.phi -= incr_val; //reset
  Calculate(proj_lab.KE,
	    ejec_lab.theta,
	    dec1_lab.theta,
	    dec1_lab.phi,
	    proj_lab.Eexc,
	    targ_lab.Eexc,
	    resi_lab.Eexc,
	    ejec_lab.Eexc,
	    dec1_lab.Eexc,
	    dec2_lab.Eexc);
  second_term *= (fxp - fx)/(hxp - hx);
  //now for phi part

  fx = dec1_deccm.phi;
  hx = cos(dec1_lab.theta);
  dec1_lab.theta += incr_val;
  Calculate(proj_lab.KE,
	    ejec_lab.theta,
	    dec1_lab.theta,
	    dec1_lab.phi,
	    proj_lab.Eexc,
	    targ_lab.Eexc,
	    resi_lab.Eexc,
	    ejec_lab.Eexc,
	    dec1_lab.Eexc,
	    dec2_lab.Eexc);
  fxp = dec1_deccm.phi;
  hxp = cos(dec1_lab.theta);
  dec1_lab.theta -= incr_val; //reset
  Calculate(proj_lab.KE,
	    ejec_lab.theta,
	    dec1_lab.theta,
	    dec1_lab.phi,
	    proj_lab.Eexc,
	    targ_lab.Eexc,
	    resi_lab.Eexc,
	    ejec_lab.Eexc,
	    dec1_lab.Eexc,
	    dec2_lab.Eexc);
  second_term *= (fxp - fx)/(hxp - hx);

  returnval -= second_term;

  return 1/returnval;

}

double FrameConverter::LabToDecCM_Dec1_Factor_ThetaOnly_Incr() {

  double returnval = 0;

  double incr_val = 1E-3;

  double fx=0, fxp=0, hx=0, hxp=0;

  fx = cos(dec1_deccm.theta);
  hx = cos(dec1_lab.theta);

  dec1_lab.theta += incr_val;

  Calculate(proj_lab.KE,
	    ejec_lab.theta,
	    dec1_lab.theta,
	    dec1_lab.phi,
	    proj_lab.Eexc,
	    targ_lab.Eexc,
	    resi_lab.Eexc,
	    ejec_lab.Eexc,
	    dec1_lab.Eexc,
	    dec2_lab.Eexc);

  fxp = cos(dec1_deccm.theta);
  hxp = cos(dec1_lab.theta);

  dec1_lab.theta -= incr_val; //reset

  Calculate(proj_lab.KE,
	    ejec_lab.theta,
	    dec1_lab.theta,
	    dec1_lab.phi,
	    proj_lab.Eexc,
	    targ_lab.Eexc,
	    resi_lab.Eexc,
	    ejec_lab.Eexc,
	    dec1_lab.Eexc,
	    dec2_lab.Eexc);

  returnval += (hxp - hx)/(fxp - fx);
    
  return returnval;

}

double FrameConverter::pE_lab() {

  double MER = ejec_lab.Mtot/resi_lab.Mtot;
  double MRP = resi_lab.Mtot/proj_lab.Mtot;
  double MEP = ejec_lab.Mtot/proj_lab.Mtot;

  double sqrt_arg = -1*pow(MER*proj_lab.ptot*sin(ejec_lab.theta),2) +
                    proj_lab.ptot*proj_lab.ptot*MER*(MRP - 1 + MEP) +
                    2*Q_rxn*ejec_lab.Mtot*(1 + MER);

  if (sqrt_arg < 0) {
    if (!quiet_mode)
      cout << "***pE_lab() WARNING: sqrt argument is NEGATIVE; returning zero for pE..." << endl;
    return 0;
  }

  double plusval  = (MER*proj_lab.ptot*cos(ejec_lab.theta) + sqrt(sqrt_arg))/(1 + MER);
  double minusval = (MER*proj_lab.ptot*cos(ejec_lab.theta) - sqrt(sqrt_arg))/(1 + MER);

  if (plusval >= 0 && minusval >= 0)
    if (!quiet_mode)
      cout << "***pE_lab() WARNING: there are TWO positive sol'ns for pE; default return plusval..." << endl;

  double returnval;

  if (plusval >= 0) {
    returnval = plusval;
  }
  else if (minusval >= 0) {
    if (!quiet_mode)
      cout << "***pE_lab() WARNING: plusval <= 0; returning minusval for pE..." << endl;
    returnval = minusval;
  }
  else {
    if (!quiet_mode)
      cout << "***pE_lab() WARNING: plusval <= 0 AND minusval <= 0; returning zero for pE..." << endl;
    returnval = 0;
  }

  return returnval;

}

double FrameConverter::pE_breakup_lab() {


  double sqrt_arg = pow(proj_lab.ptot*cos(ejec_lab.theta),2) -
                    2*(1+dec1_lab.Mtot/ejec_lab.Mtot)*(-dec1_lab.Mtot*(proj_lab.KE+Q_brk)+proj_lab.Mtot*proj_lab.KE);

  if (sqrt_arg < 0) {
    if (!quiet_mode)
      cout << "***pE_breakup_lab() WARNING: sqrt argument is NEGATIVE; returning zero for pE_breakup..." << endl;
    return 0;
  }

  double plusval  = (proj_lab.ptot*cos(ejec_lab.theta) + sqrt(sqrt_arg))/(1 + dec1_lab.Mtot/ejec_lab.Mtot);
  double minusval = (proj_lab.ptot*cos(ejec_lab.theta) - sqrt(sqrt_arg))/(1 + dec1_lab.Mtot/ejec_lab.Mtot);

  if (plusval >= 0 && minusval >= 0)
    if (!quiet_mode)
      cout << "***pE_breakup_lab() WARNING: there are TWO positive sol'ns for pE_breakup; default return plusval..." << endl;

  double returnval;

  if (plusval >= 0) {
    returnval = plusval;
  }
  else if (minusval >= 0) {
    if (!quiet_mode)
      cout << "***pE_breakup_lab() WARNING: plusval <= 0; returning minusval for pE_breakup..." << endl;
    returnval = minusval;
  }
  else {
    if (!quiet_mode)
      cout << "***pE_breakup_lab() WARNING: plusval <= 0 AND minusval <= 0; returning zero for pE_breakup..." << endl;
    returnval = 0;
  }

  return returnval;

}

double FrameConverter::pD1_lab() {

  double MD12 = dec1_lab.Mtot/dec2_lab.Mtot;
  double M2R  = dec2_lab.Mtot/resi_lab.Mtot;
  double angpart = sin(dec1_lab.theta)*sin(resi_lab.theta)*cos(dec1_lab.phi) + cos(dec1_lab.theta)*cos(resi_lab.theta);

  double sqrt_arg = pow(MD12*resi_lab.ptot*angpart,2) + (1 + MD12)*(2*Q_dec*dec1_lab.Mtot + resi_lab.ptot*resi_lab.ptot*MD12*(M2R - 1));

  if (sqrt_arg < 0) {
    if (!quiet_mode)
      cout << "***pD1_lab() WARNING: sqrt argument is NEGATIVE; returning zero for pD1..." << endl;
    return 0;
  }

  double plusval  = (MD12*resi_lab.ptot*angpart + sqrt(sqrt_arg))/(1 + MD12);
  double minusval = (MD12*resi_lab.ptot*angpart - sqrt(sqrt_arg))/(1 + MD12);

  if (plusval >= 0 && minusval >= 0)
    if (!quiet_mode)
      cout << "***pD1_lab() WARNING: there are TWO positive sol'ns for pD1; default return plusval..." << endl;

  double returnval;

  if (plusval >= 0) {
    returnval = plusval;
  }
  else if (minusval >= 0) {
    if (!quiet_mode)
      cout << "***pD1_lab() WARNING: plusval <= 0; returning minusval for pD1..." << endl;
    returnval = minusval;
  }
  else {
    if (!quiet_mode)
      cout << "***pD1_lab() WARNING: plusval <= 0 AND minusval <= 0; returning zero for pD1..." << endl;
    returnval = 0;
  }

  return returnval;

}

bool FrameConverter::GetMasses() {

  ifstream mass_file;

  mass_file.open("/home2/kgh14d/general_sim/mass_info.txt");

  if (!mass_file) {
    cout << "***WHERE IS mass_info.txt???\n";
    return false;
  }

  string column_labels;
  int Z = 0, A = 0;
  float mass = 0, nucmass = 0;
  int success=0;
  //                 Z                A           Mass_(amu)
  mass_file >> column_labels >> column_labels >> column_labels;

  while (Z != -1) {

    mass_file >> Z >> A >> mass;
    
    nucmass = (mass - Z*RESTMASS_ELECTRON)*UTOMEV;

    if (Z == proj_lab.Z && A == proj_lab.A) {
      proj_lab.M   = nucmass;
      proj_rxncm.M = nucmass;
      proj_deccm.M = nucmass;
      success++;
    }
    if (Z == targ_lab.Z && A == targ_lab.A) {
      targ_lab.M   = nucmass;
      targ_rxncm.M = nucmass;
      targ_deccm.M = nucmass;
      success++;
    }
    if (Z == resi_lab.Z && A == resi_lab.A) {
      resi_lab.M   = nucmass;
      resi_rxncm.M = nucmass;
      resi_deccm.M = nucmass;
      success++;
    }
    if (Z == ejec_lab.Z && A == ejec_lab.A) {
      ejec_lab.M   = nucmass;
      ejec_rxncm.M = nucmass;
      ejec_deccm.M = nucmass;
      success++;
    }
    if (Z == dec1_lab.Z && A == dec1_lab.A) {
      dec1_lab.M   = nucmass;
      dec1_rxncm.M = nucmass;
      dec1_deccm.M = nucmass;
      success++;
    }
    if (Z == dec2_lab.Z && A == dec2_lab.A) {
      dec2_lab.M   = nucmass;
      dec2_rxncm.M = nucmass;
      dec2_deccm.M = nucmass;
      success++;
    }

  }

  if (success!=6)
    return false;

  return true;

}

void FrameConverter::zero_structs() {

  proj_lab.Eexc=0;  proj_rxncm.Eexc=0;  proj_deccm.Eexc=0;
  proj_lab.px=0;    proj_rxncm.px=0;    proj_deccm.px=0;
  proj_lab.py=0;    proj_rxncm.py=0;    proj_deccm.py=0;
  proj_lab.pz=0;    proj_rxncm.pz=0;    proj_deccm.pz=0;
  proj_lab.ptot=0;  proj_rxncm.ptot=0;  proj_deccm.ptot=0;
  proj_lab.KE=0;    proj_rxncm.KE=0;    proj_deccm.KE=0;
  proj_lab.theta=0; proj_rxncm.theta=0; proj_deccm.theta=0;
  proj_lab.phi=0;   proj_rxncm.phi=0;   proj_deccm.phi=0;
  proj_lab.Mtot=proj_lab.M;
  proj_rxncm.Mtot=proj_rxncm.M;
  proj_deccm.Mtot=proj_deccm.M;
  proj_lab.Etot=proj_lab.Mtot;
  proj_rxncm.Etot=proj_rxncm.Mtot;
  proj_deccm.Etot=proj_deccm.Mtot;

  targ_lab.Eexc=0;  targ_rxncm.Eexc=0;  targ_deccm.Eexc=0;
  targ_lab.px=0;    targ_rxncm.px=0;    targ_deccm.px=0;
  targ_lab.py=0;    targ_rxncm.py=0;    targ_deccm.py=0;
  targ_lab.pz=0;    targ_rxncm.pz=0;    targ_deccm.pz=0;
  targ_lab.ptot=0;  targ_rxncm.ptot=0;  targ_deccm.ptot=0;
  targ_lab.KE=0;    targ_rxncm.KE=0;    targ_deccm.KE=0;
  targ_lab.theta=0; targ_rxncm.theta=0; targ_deccm.theta=0;
  targ_lab.phi=0;   targ_rxncm.phi=0;   targ_deccm.phi=0;
  targ_lab.Mtot=targ_lab.M;
  targ_rxncm.Mtot=targ_rxncm.M;
  targ_deccm.Mtot=targ_deccm.M;
  targ_lab.Etot=targ_lab.Mtot;
  targ_rxncm.Etot=targ_rxncm.Mtot;
  targ_deccm.Etot=targ_deccm.Mtot;

  resi_lab.Eexc=0;  resi_rxncm.Eexc=0;  resi_deccm.Eexc=0;
  resi_lab.px=0;    resi_rxncm.px=0;    resi_deccm.px=0;
  resi_lab.py=0;    resi_rxncm.py=0;    resi_deccm.py=0;
  resi_lab.pz=0;    resi_rxncm.pz=0;    resi_deccm.pz=0;
  resi_lab.ptot=0;  resi_rxncm.ptot=0;  resi_deccm.ptot=0;
  resi_lab.KE=0;    resi_rxncm.KE=0;    resi_deccm.KE=0;
  resi_lab.theta=0; resi_rxncm.theta=0; resi_deccm.theta=0;
  resi_lab.phi=0;   resi_rxncm.phi=0;   resi_deccm.phi=0;
  resi_lab.Mtot=resi_lab.M;
  resi_rxncm.Mtot=resi_rxncm.M;
  resi_deccm.Mtot=resi_deccm.M;
  resi_lab.Etot=resi_lab.Mtot;
  resi_rxncm.Etot=resi_rxncm.Mtot;
  resi_deccm.Etot=resi_deccm.Mtot;

  ejec_lab.Eexc=0;  ejec_rxncm.Eexc=0;  ejec_deccm.Eexc=0;
  ejec_lab.px=0;    ejec_rxncm.px=0;    ejec_deccm.px=0;
  ejec_lab.py=0;    ejec_rxncm.py=0;    ejec_deccm.py=0;
  ejec_lab.pz=0;    ejec_rxncm.pz=0;    ejec_deccm.pz=0;
  ejec_lab.ptot=0;  ejec_rxncm.ptot=0;  ejec_deccm.ptot=0;
  ejec_lab.KE=0;    ejec_rxncm.KE=0;    ejec_deccm.KE=0;
  ejec_lab.theta=0; ejec_rxncm.theta=0; ejec_deccm.theta=0;
  ejec_lab.phi=0;   ejec_rxncm.phi=0;   ejec_deccm.phi=0;
  ejec_lab.Mtot=ejec_lab.M;
  ejec_rxncm.Mtot=ejec_rxncm.M;
  ejec_deccm.Mtot=ejec_deccm.M;
  ejec_lab.Etot=ejec_lab.Mtot;
  ejec_rxncm.Etot=ejec_rxncm.Mtot;
  ejec_deccm.Etot=ejec_deccm.Mtot;

  dec1_lab.Eexc=0;  dec1_rxncm.Eexc=0;  dec1_deccm.Eexc=0;
  dec1_lab.px=0;    dec1_rxncm.px=0;    dec1_deccm.px=0;
  dec1_lab.py=0;    dec1_rxncm.py=0;    dec1_deccm.py=0;
  dec1_lab.pz=0;    dec1_rxncm.pz=0;    dec1_deccm.pz=0;
  dec1_lab.ptot=0;  dec1_rxncm.ptot=0;  dec1_deccm.ptot=0;
  dec1_lab.KE=0;    dec1_rxncm.KE=0;    dec1_deccm.KE=0;
  dec1_lab.theta=0; dec1_rxncm.theta=0; dec1_deccm.theta=0;
  dec1_lab.phi=0;   dec1_rxncm.phi=0;   dec1_deccm.phi=0;
  dec1_lab.Mtot=dec1_lab.M;
  dec1_rxncm.Mtot=dec1_rxncm.M;
  dec1_deccm.Mtot=dec1_deccm.M;
  dec1_lab.Etot=dec1_lab.Mtot;
  dec1_rxncm.Etot=dec1_rxncm.Mtot;
  dec1_deccm.Etot=dec1_deccm.Mtot;

  dec2_lab.Eexc=0;  dec2_rxncm.Eexc=0;  dec2_deccm.Eexc=0;
  dec2_lab.px=0;    dec2_rxncm.px=0;    dec2_deccm.px=0;
  dec2_lab.py=0;    dec2_rxncm.py=0;    dec2_deccm.py=0;
  dec2_lab.pz=0;    dec2_rxncm.pz=0;    dec2_deccm.pz=0;
  dec2_lab.ptot=0;  dec2_rxncm.ptot=0;  dec2_deccm.ptot=0;
  dec2_lab.KE=0;    dec2_rxncm.KE=0;    dec2_deccm.KE=0;
  dec2_lab.theta=0; dec2_rxncm.theta=0; dec2_deccm.theta=0;
  dec2_lab.phi=0;   dec2_rxncm.phi=0;   dec2_deccm.phi=0;
  dec2_lab.Mtot=dec2_lab.M;
  dec2_rxncm.Mtot=dec2_rxncm.M;
  dec2_deccm.Mtot=dec2_deccm.M;
  dec2_lab.Etot=dec2_lab.Mtot;
  dec2_rxncm.Etot=dec2_rxncm.Mtot;
  dec2_deccm.Etot=dec2_deccm.Mtot;

  proj_vec.SetPxPyPzE(proj_lab.px,proj_lab.py,proj_lab.pz,proj_lab.Etot);
  targ_vec.SetPxPyPzE(targ_lab.px,targ_lab.py,targ_lab.pz,targ_lab.Etot);
  resi_vec.SetPxPyPzE(resi_lab.px,resi_lab.py,resi_lab.pz,resi_lab.Etot);
  ejec_vec.SetPxPyPzE(ejec_lab.px,ejec_lab.py,ejec_lab.pz,ejec_lab.Etot);
  dec1_vec.SetPxPyPzE(dec1_lab.px,dec1_lab.py,dec1_lab.pz,dec1_lab.Etot);
  dec2_vec.SetPxPyPzE(dec2_lab.px,dec2_lab.py,dec2_lab.pz,dec2_lab.Etot);

  rxn_boost_vec.SetXYZ(0,0,0);
  dec_boost_vec.SetXYZ(0,0,0);

}

void FrameConverter::QuietOn()  {quiet_mode = true;}
void FrameConverter::QuietOff() {quiet_mode = false;}
