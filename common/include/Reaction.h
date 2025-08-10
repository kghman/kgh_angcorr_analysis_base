#ifndef __Reaction_h
#define __Reaction_h

#include <iostream>
#include <cstdlib>
#include <fstream>
#include <string>
#include <iomanip>
#include <cmath>

#include <TRandom3.h>
#include <TLorentzVector.h>

#include "constants.h"
#include "nucleusdata.h"

/*All input and output angles are in RADIANS*/

class Reaction {

 public:

  Reaction();

  Reaction(int iID,
	   int iTZ=1, int iTA=1, int iPZ=1, int iPA=1,
	   int iEZ=0, int iEA=0, int iRZ=0, int iRA=0);

  ~Reaction();

  void Initialize(int iID,
		  int iTZ=1, int iTA=1, int iPZ=1, int iPA=1,
		  int iEZ=0, int iEA=0, int iRZ=0, int iRA=0);

  /* Get Functions */
  int ID();

  int Target_Z();
  int Projectile_Z();
  int Ejectile_Z();
  int Residual_Z();

  int Target_A();
  int Projectile_A();
  int Ejectile_A();
  int Residual_A();

  double Target_Restmass();
  double Projectile_Restmass();
  double Ejectile_Restmass();
  double Residual_Restmass();

  double Target_Charge();
  double Projectile_Charge();
  double Ejectile_Charge();
  double Residual_Charge();

  double Q_Value();

  double Rutherford_XSec(double ang_cm); //input in radians, output in mb/sr

  /*Functions to return kinematics values */
  //expects an integer representing a particular reaction participant:
  // T(P,E)R = 1(2,3)4
  // defaults to ejectile
  //return units are in MeV, radians, and unit lengths
  double Get_px_lab(int    ipart_id=3);
  double Get_py_lab(int    ipart_id=3);
  double Get_pz_lab(int    ipart_id=3);
  double Get_ptot_lab(int  ipart_id=3);
  double Get_px_cm(int     ipart_id=3);
  double Get_py_cm(int     ipart_id=3);
  double Get_pz_cm(int     ipart_id=3);
  double Get_ptot_cm(int   ipart_id=3);
  double Get_E_tot_lab(int ipart_id=3);
  double Get_E_tot_cm(int  ipart_id=3);
  double Get_KE_lab(int    ipart_id=3);
  double Get_KE_cm(int     ipart_id=3);
  double Get_E_ex(int      ipart_id=3);
  double Get_invmass(int   ipart_id=3);
  double Get_theta_lab(int ipart_id=3);
  double Get_phi_lab(int   ipart_id=3);
  double Get_theta_cm(int  ipart_id=3);
  double Get_phi_cm(int    ipart_id=3);
  double Get_x(int         ipart_id=3);
  double Get_y(int         ipart_id=3);
  double Get_z(int         ipart_id=3);

  double Entrance_Channel_KEcm();
  double Exit_Channel_KEcm();

  bool Load_XSec(char *filename);

  /* Kinematics-Related Functions */
  //quantities in lab frame
  void Run_Kinematics(double Beam_E,
		      double theta_proj=0., double phi_proj=0.,
		      double targ_exc=-1,   double proj_exc=-1.,
                      double ejec_exc=-1,   double res_exc=-1);

  double Random_Ejec_Theta_CM(double thmin=0, double thmax=M_PI);
  double Random_Ejec_Theta_Lab(double thmin=0, double thmax=M_PI, double beamE=0,
			       double eexc=-1, double rexc=-1);
  double Random_Ejec_Phi(double phmin=0, double phmax=2*M_PI);

  double Random_Target_Exc_Energy();
  double Random_Projectile_Exc_Energy();
  double Random_Ejectile_Exc_Energy();
  double Random_Residual_Exc_Energy();

  double Random_Exc_Energy(int);

  double LabToCM_Angle(double angle, double beamE, double ejec_exc, double res_exc);
  double CMToLab_Angle(double angle, double beamE, double ejec_exc, double res_exc);

  double LabToCM_Factor(double angle, double beamE, double ejec_exc, double res_exc);

  void Print();

  void Reset();

 private:

  int Rxn_ID; //to label a given instantiation of this class for multiple rxns

  int TARGET_Z,     TARGET_A,
      PROJECTILE_Z, PROJECTILE_A,
      EJECTILE_Z,   EJECTILE_A,
      RESIDUAL_Z,   RESIDUAL_A;

  double RESTMASS_TARGET,
         RESTMASS_PROJECTILE,
         RESTMASS_EJECTILE,
         RESTMASS_RESIDUAL,
         RESTMASS_DECAY1,
         RESTMASS_DECAY2;

  double CHARGE_TARGET,
         CHARGE_PROJECTILE,
         CHARGE_EJECTILE,
         CHARGE_RESIDUAL,
         CHARGE_DECAY1,
         CHARGE_DECAY2;
  
  double Q_VALUE;

  //Z, A, and restmass arrays
  int Z_arr[4], A_arr[4];
  double M_arr[4];

  //to hold exc state energy, width, and relative fractions
  int   num_exc_states[4];
  float exc_state_info[4][3][100];

  //all 4 reaction participants in order: T(P,E)R
  NucleusData rxn_parts[4];

  //if given a file of diff xsec, generate probability dist.
  double **angle_weights;
  int tot_angs;

  bool xsec_loaded;

  double Get_KEcm(int, int);

  bool GetMasses();
  bool GetStates();

  void CheckRxn();

  double Get_Exc_State(float[3][100], int); //param arrays, number of states

  void NucData_Zeros(); //initialize all nucleusdata variables to zero

  TRandom3 *rand;

};

#endif
