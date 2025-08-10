#ifndef __decay_h
#define __decay_h

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

class Decay {

 public:

  Decay();
  
  Decay(int iID,
	int iPZ=1,  int iPA=1,
	int iD1Z=0, int iD1A=0, int iD2Z=0, int iD2A=0);
  
  ~Decay();

  void Initialize(int iID=1,
		  int iPZ=1,  int iPA=1,
		  int iD1Z=0, int iD1A=0, int iD2Z=0, int iD2A=0);

  int ID();

  int Parent_Z();
  int Decay1_Z();
  int Decay2_Z();

  int Parent_A();
  int Decay1_A();
  int Decay2_A();

  double Parent_Restmass();
  double Decay1_Restmass();
  double Decay2_Restmass();

  double Parent_Charge();
  double Decay1_Charge();
  double Decay2_Charge();

  double Q_Value();

  /*Functions to return kinematics values */
  //expects an integer representing a particular reaction participant:
  // P --> D1 + D2 = 1 --> 2 + 3
  // defaults to parent
  double Get_px_lab(int    ipart_id=1);
  double Get_py_lab(int    ipart_id=1);
  double Get_pz_lab(int    ipart_id=1);
  double Get_ptot_lab(int  ipart_id=1);
  double Get_px_cm(int     ipart_id=1);
  double Get_py_cm(int     ipart_id=1);
  double Get_pz_cm(int     ipart_id=1);
  double Get_ptot_cm(int   ipart_id=1);
  double Get_E_tot_lab(int ipart_id=1);
  double Get_E_tot_cm(int  ipart_id=1);
  double Get_KE_lab(int    ipart_id=1);
  double Get_KE_cm(int     ipart_id=1);
  double Get_E_ex(int      ipart_id=1);
  double Get_invmass(int   ipart_id=1);
  double Get_theta_lab(int ipart_id=1);
  double Get_phi_lab(int   ipart_id=1);
  double Get_theta_cm(int  ipart_id=1);
  double Get_phi_cm(int    ipart_id=1);
  double Get_x(int         ipart_id=1);
  double Get_y(int         ipart_id=1);
  double Get_z(int         ipart_id=1);

  double Get_decay_sep_angle_lab();
  double Get_decay_sep_angle_cm();

  double Random_Parent_Theta(double thmin=0, double thmax=M_PI);
  double Random_Parent_Phi(double phmin=0, double phmax=2*M_PI);

  double Random_Decay1_Theta_CM(double thmin=0, double thmax=M_PI);
  double Random_Decay1_Phi_CM(double phmin=0, double phmax=2*M_PI);

  double Random_Parent_Exc_Energy();
  double Random_Decay1_Exc_Energy();
  double Random_Decay2_Exc_Energy();

  double Random_Exc_Energy(int);

  //all of the following expects decay 1 angles in CM (-1 denotes random value)
  //assumes parent at rest at origin:
  void Run_Decay(double d1_theta=-1, double d1_phi=-1,
		 double P_exc=-1, double D1_exc=-1, double D2_exc=-1);
  //expects parent mom and E:
  void Run_Decay_Vector(double iP_px, double iP_py, double iP_pz, double iP_E,
			double d1_theta=-1, double d1_phi=-1,
			double D1_exc=-1, double D2_exc=-1);
  //expects parent KE and angles:
  void Run_Decay_Angles(double iP_KE, double iP_theta=-1, double iP_ph=-1,
			double d1_theta=-1, double d1_phi=-1,
			double P_exc=-1, double D1_exc=-1, double D2_exc=-1);

 private:

  int Decay_ID;

  int PARENT_Z, PARENT_A,
      DECAY1_Z, DECAY1_A,
      DECAY2_Z, DECAY2_A;

  double RESTMASS_PARENT,
         RESTMASS_DECAY1,
         RESTMASS_DECAY2;

  double CHARGE_PARENT,
         CHARGE_DECAY1,
         CHARGE_DECAY2;

  double Q_VALUE;

  int Z_arr[3], A_arr[3];
  double M_arr[4];

  int   num_exc_states[3];
  float exc_state_info[3][3][100];

  double decay_sep_angle_lab; //angle between decay products
  double decay_sep_angle_cm;

  NucleusData decay_parts[3];

  TRandom3 *rand;

  bool GetMasses();
  bool GetStates();

  void CheckDecay();

  double Get_Exc_State(float[3][100], int); //param arrays, number of states

  void NucData_Zeros(); //initialize all nucleusdata variables to zero

};

#endif
