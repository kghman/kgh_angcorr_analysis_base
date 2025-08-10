#ifndef __FRAMECONVERTER_H
#define __FRAMECONVERTER_H

#include <TVector3.h>
#include <TLorentzVector.h>

#include "singleframenucleus.h"

class FrameConverter {

 public:

  FrameConverter(int ZP,int AP,
		 int ZT,int AT,
		 int ZE,int AE,
		 int ZDL,int ADL);
  ~FrameConverter();

  //all inputs in LAB frame (targ at rest)
  bool Calculate(double beamE,
		 double ejec_theta,
		 double ejec_phi,
		 double decL_theta,
		 double decL_phi,
		 double proj_exc,
		 double targ_exc,
		 double resi_exc,
		 double ejec_exc,
		 double decL_exc,
		 double decH_exc);

  bool Calculate_Breakup(double beamE,
			 double ejec_theta,
			 double ejec_phi,
			 double proj_exc,
			 double ejec_exc,
			 double decL_exc);
  
  //all inputs in LAB frame (targ at rest)
  //EXCEPT decL theta/phi, which are in decay (resi) rest frame
  bool Calculate_DecCM(double beamE,
		       double ejec_theta,
		       double ejec_phi,
		       double decL_theta_decCM,
		       double decL_phi_decCM,
		       double proj_exc,
		       double targ_exc,
		       double resi_exc,
		       double ejec_exc,
		       double decL_exc,
		       double decH_exc);

  double LabToRxnCM_Ejec_Angle(double angle, double beamE, double proj_exc, double targ_exc, double ejec_exc, double resi_exc);
  double RxnCMToLab_Ejec_Angle(double angle, double beamE, double proj_exc, double targ_exc, double ejec_exc, double resi_exc);

  double LabToRxnCM_Ejec_Angle_Illiadis(double angle, double beamE, double proj_exc, double targ_exc, double ejec_exc, double resi_exc);
  double RxnCMToLab_Ejec_Angle_Illiadis(double angle, double beamE, double proj_exc, double targ_exc, double ejec_exc, double resi_exc);

  double LabToRxnCM_Ejec_Factor();
  double LabToRxnCM_Ejec_Factor_Incr();

  double LabToRxnCM_Resi_Factor_Incr();
  double LabToDecCM_Resi_Factor_Incr();

  double LabToRxnCM_DecL_Factor_Incr();
  double LabToDecCM_DecL_Factor_Incr();
  double LabToDecCM_DecL_Factor_ThetaOnly_Incr();

  double Q_rxn, Q_dec, Q_brk;

  singleframenucleus
    proj_lab,
    targ_lab,
    resi_lab,
    ejec_lab,
    decL_lab,
    decH_lab;

  singleframenucleus
    proj_rxncm,
    targ_rxncm,
    resi_rxncm,
    ejec_rxncm,
    decL_rxncm,
    decH_rxncm;

  singleframenucleus
    proj_deccm,
    targ_deccm,
    resi_deccm,
    ejec_deccm,
    decL_deccm,
    decH_deccm;

  singleframenucleus
    proj_brkcm,
    targ_brkcm,
    resi_brkcm,
    ejec_brkcm,
    decL_brkcm,
    decH_brkcm;

  //choose whether to display error messages
  //(like if you have some "bad" kinematic events in a loop,
  // and you don't want to see the errors each time...)
  void QuietOn();
  void QuietOff();

 private:

  TLorentzVector
    proj_vec,
    targ_vec,
    resi_vec,
    ejec_vec,
    decL_vec,
    decH_vec;

  TLorentzVector
    rxn_parent_vec;

  TVector3
    rxn_boost_vec,
    dec_boost_vec,
    brk_boost_vec;
  
  bool quiet_mode;

  double pE_lab();
  double pE_breakup_lab();
  double pDL_lab();

  bool GetMasses();

  void zero_structs();
  
};

#endif
