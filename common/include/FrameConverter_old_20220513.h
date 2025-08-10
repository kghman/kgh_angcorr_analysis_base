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
		 int ZD1,int AD1);
  ~FrameConverter();

  //all inputs in LAB frame (targ at rest)
  bool Calculate(double beamE,
		 double ejec_theta,
		 double dec1_theta,
		 double dec1_phi,
		 double proj_exc,
		 double targ_exc,
		 double resi_exc,
		 double ejec_exc,
		 double dec1_exc,
		 double dec2_exc);

  bool Calculate_Breakup(double beamE,
			 double ejec_theta,
			 double proj_exc,
			 double ejec_exc,
			 double dec1_exc);
  
  //all inputs in LAB frame (targ at rest)
  //EXCEPT dec1 theta/phi, which are in decay (resi) rest frame
  bool Calculate_DecCM(double beamE,
		       double ejec_theta,
		       double dec1_theta_decCM,
		       double dec1_phi_decCM,
		       double proj_exc,
		       double targ_exc,
		       double resi_exc,
		       double ejec_exc,
		       double dec1_exc,
		       double dec2_exc);

  double LabToRxnCM_Ejec_Angle(double angle, double beamE, double proj_exc, double targ_exc, double ejec_exc, double resi_exc);
  double RxnCMToLab_Ejec_Angle(double angle, double beamE, double proj_exc, double targ_exc, double ejec_exc, double resi_exc);

  double LabToRxnCM_Ejec_Angle_Illiadis(double angle, double beamE, double proj_exc, double targ_exc, double ejec_exc, double resi_exc);
  double RxnCMToLab_Ejec_Angle_Illiadis(double angle, double beamE, double proj_exc, double targ_exc, double ejec_exc, double resi_exc);

  double LabToRxnCM_Ejec_Factor();
  double LabToRxnCM_Ejec_Factor_Incr();

  double LabToRxnCM_Resi_Factor_Incr();
  double LabToDecCM_Resi_Factor_Incr();

  double LabToRxnCM_Dec1_Factor_Incr();
  double LabToDecCM_Dec1_Factor_Incr();
  double LabToDecCM_Dec1_Factor_ThetaOnly_Incr();

  double Q_rxn, Q_dec, Q_brk;

  singleframenucleus
    proj_lab,
    targ_lab,
    resi_lab,
    ejec_lab,
    dec1_lab,
    dec2_lab;

  singleframenucleus
    proj_rxncm,
    targ_rxncm,
    resi_rxncm,
    ejec_rxncm,
    dec1_rxncm,
    dec2_rxncm;

  singleframenucleus
    proj_deccm,
    targ_deccm,
    resi_deccm,
    ejec_deccm,
    dec1_deccm,
    dec2_deccm;

  singleframenucleus
    proj_brkcm,
    targ_brkcm,
    resi_brkcm,
    ejec_brkcm,
    dec1_brkcm,
    dec2_brkcm;

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
    dec1_vec,
    dec2_vec;

  TLorentzVector
    rxn_parent_vec;

  TVector3
    rxn_boost_vec,
    dec_boost_vec,
    brk_boost_vec;
  
  bool quiet_mode;

  double pE_lab();
  double pE_breakup_lab();
  double pD1_lab();

  bool GetMasses();

  void zero_structs();
  
};

#endif
