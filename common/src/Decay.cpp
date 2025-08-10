#include "Decay.h"

using namespace std;

Decay::Decay() {

  //called when creating an array of Decay objects,
  //each of which should then be initialized via Decay::Initialize(...)

  Decay_ID = 0;
  PARENT_Z = 0;
  PARENT_A = 0;
  DECAY1_Z = 0;
  DECAY1_A = 0;
  DECAY2_Z = 0;
  DECAY2_A = 0;

  CHARGE_PARENT = 0;
  CHARGE_DECAY1 = 0;
  CHARGE_DECAY2 = 0;

  for (int k=0; k<3; k++) {
    Z_arr[k] = 0;
    A_arr[k] = 0;
    num_exc_states[k] = 0;
    for (int i=0; i<3; i++)
      for (int j=0; j<100; j++)
	exc_state_info[k][i][j] = 0;
  }

  NucData_Zeros();

  Q_VALUE = 0;

  rand = new TRandom3();
  rand->SetSeed();

}

Decay::Decay(int iID,
	     int iPZ,  int iPA,
	     int iD1Z, int iD1A, int iD2Z, int iD2A) {

  Decay_ID = iID;

  PARENT_Z = iPZ;
  PARENT_A = iPA;

  if ((iD1Z == 0 && iD1A == 0) ||
      (iD2Z == 0 && iD2A == 0)) { //default to half of parent vals (careful of odd integers...)
    DECAY1_Z = PARENT_Z/2;
    DECAY1_A = PARENT_A/2;
    DECAY2_Z = PARENT_Z/2;
    DECAY2_A = PARENT_A/2;
  }
  else {
    DECAY1_Z = iD1Z;
    DECAY1_A = iD1A;
    DECAY2_Z = iD2Z;
    DECAY2_A = iD2A;
  }

  CHARGE_PARENT = PARENT_Z*UNIT_CHARGE;
  CHARGE_DECAY1 = DECAY1_Z*UNIT_CHARGE;
  CHARGE_DECAY2 = DECAY2_Z*UNIT_CHARGE;

  Z_arr[0] = PARENT_Z;
  Z_arr[1] = DECAY1_Z;
  Z_arr[2] = DECAY2_Z;

  A_arr[0] = PARENT_A;
  A_arr[1] = DECAY1_A;
  A_arr[2] = DECAY2_A;

  CheckDecay();

  for (int k=0; k<3; k++) {
    num_exc_states[k] = 1;
    for (int i=0; i<3; i++)
      for (int j=0; j<100; j++)
	exc_state_info[k][i][j] = 0;
  }

  NucData_Zeros();

  if (!GetMasses()) cout << "***Decay() warning: problem loading rest masses\n";
  if (!GetStates()) cout << "***Decay() warning: problem loading excited states\n";

  Q_VALUE = RESTMASS_PARENT - (RESTMASS_DECAY1 + RESTMASS_DECAY2);

  rand = new TRandom3();
  rand->SetSeed();

}

void Decay::Initialize(int iID,
		       int iPZ,  int iPA,
		       int iD1Z, int iD1A, int iD2Z, int iD2A) {

  Decay_ID = iID;

  PARENT_Z = iPZ;
  PARENT_A = iPA;

  if ((iD1Z == 0 && iD1A == 0) ||
      (iD2Z == 0 && iD2A == 0)) { //default to half of parent vals (careful of odd integers...)
    DECAY1_Z = PARENT_Z/2;
    DECAY1_A = PARENT_A/2;
    DECAY2_Z = PARENT_Z/2;
    DECAY2_A = PARENT_A/2;
  }
  else {
    DECAY1_Z = iD1Z;
    DECAY1_A = iD1A;
    DECAY2_Z = iD2Z;
    DECAY2_A = iD2A;
  }

  CHARGE_PARENT = PARENT_Z*UNIT_CHARGE;
  CHARGE_DECAY1 = DECAY1_Z*UNIT_CHARGE;
  CHARGE_DECAY2 = DECAY2_Z*UNIT_CHARGE;

  Z_arr[0] = PARENT_Z;
  Z_arr[1] = DECAY1_Z;
  Z_arr[2] = DECAY2_Z;

  A_arr[0] = PARENT_A;
  A_arr[1] = DECAY1_A;
  A_arr[2] = DECAY2_A;

  CheckDecay();

  for (int k=0; k<3; k++) {
    num_exc_states[k] = 1;
    for (int i=0; i<3; i++)
      for (int j=0; j<100; j++)
	exc_state_info[k][i][j] = 0;
  }

  if (!GetMasses()) cout << "***Decay() warning: problem loading rest masses\n";
  if (!GetStates()) cout << "***Decay() warning: problem loading excited states\n";

  Q_VALUE = RESTMASS_PARENT - (RESTMASS_DECAY1 + RESTMASS_DECAY2);

}

Decay::~Decay() {

  delete rand;

}

int Decay::ID() {return Decay_ID;}

int Decay::Parent_Z() {return PARENT_Z;}
int Decay::Decay1_Z() {return DECAY1_Z;}
int Decay::Decay2_Z() {return DECAY2_Z;}

int Decay::Parent_A() {return PARENT_A;}
int Decay::Decay1_A() {return DECAY1_A;}
int Decay::Decay2_A() {return DECAY2_A;}

double Decay::Parent_Restmass() {return RESTMASS_PARENT;}
double Decay::Decay1_Restmass() {return RESTMASS_DECAY1;}
double Decay::Decay2_Restmass() {return RESTMASS_DECAY2;}

double Decay::Parent_Charge() {return CHARGE_PARENT;}
double Decay::Decay1_Charge() {return CHARGE_DECAY1;}
double Decay::Decay2_Charge() {return CHARGE_DECAY2;}

double Decay::Q_Value() {return Q_VALUE;}

double Decay::Get_decay_sep_angle_lab() {return decay_sep_angle_lab;}
double Decay::Get_decay_sep_angle_cm()  {return decay_sep_angle_cm;}

double Decay::Get_px_lab(int ipart_id)
{return (ipart_id > 0 && ipart_id <= 3) ? decay_parts[ipart_id - 1].px_lab : 0.0;}
double Decay::Get_py_lab(int ipart_id)
{return (ipart_id > 0 && ipart_id <= 3) ? decay_parts[ipart_id - 1].py_lab : 0.0;}
double Decay::Get_pz_lab(int ipart_id)
{return (ipart_id > 0 && ipart_id <= 3) ? decay_parts[ipart_id - 1].pz_lab : 0.0;}
double Decay::Get_ptot_lab(int ipart_id)
{return (ipart_id > 0 && ipart_id <= 3) ? decay_parts[ipart_id - 1].ptot_lab : 0.0;}
double Decay::Get_px_cm(int ipart_id)
{return (ipart_id > 0 && ipart_id <= 3) ? decay_parts[ipart_id - 1].px_cm : 0.0;}
double Decay::Get_py_cm(int ipart_id)
{return (ipart_id > 0 && ipart_id <= 3) ? decay_parts[ipart_id - 1].py_cm : 0.0;}
double Decay::Get_pz_cm(int ipart_id)
{return (ipart_id > 0 && ipart_id <= 3) ? decay_parts[ipart_id - 1].pz_cm : 0.0;}
double Decay::Get_ptot_cm(int ipart_id)
{return (ipart_id > 0 && ipart_id <= 3) ? decay_parts[ipart_id - 1].ptot_cm : 0.0;}
double Decay::Get_E_tot_lab(int ipart_id)
{return (ipart_id > 0 && ipart_id <= 3) ? decay_parts[ipart_id - 1].E_tot_lab : 0.0;}
double Decay::Get_E_tot_cm(int ipart_id)
{return (ipart_id > 0 && ipart_id <= 3) ? decay_parts[ipart_id - 1].E_tot_cm : 0.0;}
double Decay::Get_KE_lab(int ipart_id)
{return (ipart_id > 0 && ipart_id <= 3) ? decay_parts[ipart_id - 1].KE_lab : 0.0;}
double Decay::Get_KE_cm(int ipart_id)
{return (ipart_id > 0 && ipart_id <= 3) ? decay_parts[ipart_id - 1].KE_cm : 0.0;}
double Decay::Get_E_ex(int ipart_id)
{return (ipart_id > 0 && ipart_id <= 3) ? decay_parts[ipart_id - 1].E_ex : 0.0;}
double Decay::Get_invmass(int ipart_id)
{return (ipart_id > 0 && ipart_id <= 3) ? decay_parts[ipart_id - 1].invmass : 0.0;}
double Decay::Get_theta_lab(int ipart_id)
{return (ipart_id > 0 && ipart_id <= 3) ? decay_parts[ipart_id - 1].theta_lab : 0.0;}
double Decay::Get_phi_lab(int ipart_id)
{return (ipart_id > 0 && ipart_id <= 3) ? decay_parts[ipart_id - 1].phi_lab : 0.0;}
double Decay::Get_theta_cm(int ipart_id)
{return (ipart_id > 0 && ipart_id <= 3) ? decay_parts[ipart_id - 1].theta_cm : 0.0;}
double Decay::Get_phi_cm(int ipart_id)
{return (ipart_id > 0 && ipart_id <= 3) ? decay_parts[ipart_id - 1].phi_cm : 0.0;}
double Decay::Get_x(int ipart_id)
{return (ipart_id > 0 && ipart_id <= 3) ? decay_parts[ipart_id - 1].x : 0.0;}
double Decay::Get_y(int ipart_id)
{return (ipart_id > 0 && ipart_id <= 3) ? decay_parts[ipart_id - 1].y : 0.0;}
double Decay::Get_z(int ipart_id)
{return (ipart_id > 0 && ipart_id <= 3) ? decay_parts[ipart_id - 1].z : 0.0;}

double Decay::Random_Parent_Theta(double thmin, double thmax) {

  return rand->Uniform(thmin,thmax);

}

double Decay::Random_Parent_Phi(double phmin, double phmax) {

  return rand->Uniform(phmin,phmax);

}

double Decay::Random_Decay1_Theta_CM(double thmin, double thmax) {

  return rand->Uniform(thmin,thmax);

}

double Decay::Random_Decay1_Phi_CM(double phmin, double phmax) {

  return rand->Uniform(phmin,phmax);

}

double Decay::Random_Parent_Exc_Energy() {return Random_Exc_Energy(0);}
double Decay::Random_Decay1_Exc_Energy() {return Random_Exc_Energy(1);}
double Decay::Random_Decay2_Exc_Energy() {return Random_Exc_Energy(2);}

double Decay::Random_Exc_Energy(int particleID) {

  return Get_Exc_State(exc_state_info[particleID],num_exc_states[particleID]);

}

void Decay::Run_Decay(double d1_theta, double d1_phi,
		      double P_exc, double D1_exc, double D2_exc) {

  if (d1_theta == -1) d1_theta = Random_Decay1_Theta_CM();
  if (d1_phi   == -1) d1_phi   = Random_Decay1_Phi_CM();

  if (P_exc  == -1) P_exc  = Random_Parent_Exc_Energy();
  if (D1_exc == -1) D1_exc = Random_Decay1_Exc_Energy();
  if (D2_exc == -1) D2_exc = Random_Decay2_Exc_Energy();

  double parent_E = RESTMASS_PARENT + P_exc;
  
  Run_Decay_Vector(0, 0, 0, parent_E, d1_theta, d1_phi, D1_exc, D2_exc);

}

void Decay::Run_Decay_Vector(double iP_px, double iP_py, double iP_pz, double iP_E,
			     double d1_theta, double d1_phi,
			     double D1_exc, double D2_exc) {

  if (d1_theta == -1) d1_theta = Random_Decay1_Theta_CM();
  if (d1_phi   == -1) d1_phi   = Random_Decay1_Phi_CM();

  if (D1_exc == -1) D1_exc = Random_Decay1_Exc_Energy();
  if (D2_exc == -1) D2_exc = Random_Decay2_Exc_Energy();

  double TOTAL_RESTMASS_DECAY1 = RESTMASS_DECAY1 + D1_exc;
  double TOTAL_RESTMASS_DECAY2 = RESTMASS_DECAY2 + D2_exc;

  TLorentzVector Parent_LV(iP_px, iP_py, iP_pz, iP_E);

  decay_parts[0].px_lab    = Parent_LV.Px();
  decay_parts[0].py_lab    = Parent_LV.Py();
  decay_parts[0].pz_lab    = Parent_LV.Pz();
  decay_parts[0].ptot_lab  = Parent_LV.P();
  decay_parts[0].E_tot_lab = Parent_LV.E();
  decay_parts[0].E_ex      = Parent_LV.M() - RESTMASS_PARENT;
  double TOTAL_RESTMASS_PARENT = RESTMASS_PARENT + decay_parts[0].E_ex;
  decay_parts[0].KE_lab    = Parent_LV.E() - TOTAL_RESTMASS_PARENT;
  decay_parts[0].invmass   = Parent_LV.M();
  decay_parts[0].theta_lab = Parent_LV.Theta();
  decay_parts[0].phi_lab   = Parent_LV.Phi();
  decay_parts[0].x         = 0; //assume happens at reaction point
  decay_parts[0].y         = 0;
  decay_parts[0].z         = 0;

  TVector3 boost_V_Parent = Parent_LV.BoostVector();
  Parent_LV.Boost(-boost_V_Parent); //to CM frame
  double parent_E_cm = Parent_LV.E();

  decay_parts[0].px_cm    = Parent_LV.Px();
  decay_parts[0].py_cm    = Parent_LV.Py();
  decay_parts[0].pz_cm    = Parent_LV.Pz();
  decay_parts[0].ptot_cm  = Parent_LV.P();
  decay_parts[0].E_tot_cm = Parent_LV.E();
  decay_parts[0].KE_cm    = Parent_LV.E() - TOTAL_RESTMASS_PARENT;
  decay_parts[0].theta_cm = Parent_LV.Theta();
  decay_parts[0].phi_cm   = Parent_LV.Phi();

  //random CM angles for first decay:
  decay_parts[1].theta_cm = d1_theta;
  decay_parts[1].phi_cm   = d1_phi;

  //now calculate first decay params:
  decay_parts[1].E_tot_cm = (pow(TOTAL_RESTMASS_DECAY1,2) -
			     pow(TOTAL_RESTMASS_DECAY2,2) +
			     pow(parent_E_cm,2)) / (2*parent_E_cm);
  decay_parts[1].KE_cm   = decay_parts[1].E_tot_cm - TOTAL_RESTMASS_DECAY1;
  decay_parts[1].ptot_cm = sqrt(pow(decay_parts[1].E_tot_cm,2) - pow(TOTAL_RESTMASS_DECAY1,2));
  decay_parts[1].px_cm   = decay_parts[1].ptot_cm*sin(decay_parts[1].theta_cm)*cos(decay_parts[1].phi_cm);
  decay_parts[1].py_cm   = decay_parts[1].ptot_cm*sin(decay_parts[1].theta_cm)*sin(decay_parts[1].phi_cm);
  decay_parts[1].pz_cm   = decay_parts[1].ptot_cm*cos(decay_parts[1].theta_cm);

  TLorentzVector Decay1_LV(decay_parts[1].px_cm,
			   decay_parts[1].py_cm,
			   decay_parts[1].pz_cm,
			   decay_parts[1].E_tot_cm);

  //back to lab frame...
  Decay1_LV.Boost(boost_V_Parent);
  Parent_LV.Boost(boost_V_Parent);

  decay_parts[1].E_tot_lab = Decay1_LV.E();
  decay_parts[1].invmass   = Decay1_LV.M();
  decay_parts[1].KE_lab    = decay_parts[1].E_tot_lab - TOTAL_RESTMASS_DECAY1;
  decay_parts[1].E_ex      = decay_parts[1].invmass - RESTMASS_DECAY1;
  decay_parts[1].ptot_lab  = Decay1_LV.P();
  decay_parts[1].px_lab    = Decay1_LV.Px();
  decay_parts[1].py_lab    = Decay1_LV.Py();
  decay_parts[1].pz_lab    = Decay1_LV.Pz();
  decay_parts[1].theta_lab = Decay1_LV.Theta();
  decay_parts[1].phi_lab   = Decay1_LV.Phi();

  decay_parts[1].x = cos(decay_parts[1].phi_lab)*sin(decay_parts[1].theta_lab);
  decay_parts[1].y = sin(decay_parts[1].phi_lab)*sin(decay_parts[1].theta_lab);
  decay_parts[1].z = cos(decay_parts[1].theta_lab);

  //second decay:
  TLorentzVector Decay2_LV = Parent_LV - Decay1_LV;

  decay_parts[2].E_tot_lab = Decay2_LV.E();
  decay_parts[2].invmass   = Decay2_LV.M();
  decay_parts[2].KE_lab    = decay_parts[2].E_tot_lab - TOTAL_RESTMASS_DECAY2;
  decay_parts[2].E_ex      = decay_parts[2].invmass - RESTMASS_DECAY2;
  decay_parts[2].ptot_lab  = Decay2_LV.P();
  decay_parts[2].px_lab    = Decay2_LV.Px();
  decay_parts[2].py_lab    = Decay2_LV.Py();
  decay_parts[2].pz_lab    = Decay2_LV.Pz();
  decay_parts[2].theta_lab = Decay2_LV.Theta();
  decay_parts[2].phi_lab   = Decay2_LV.Phi();

  decay_parts[2].x = cos(decay_parts[2].phi_lab)*sin(decay_parts[2].theta_lab);
  decay_parts[2].y = sin(decay_parts[2].phi_lab)*sin(decay_parts[2].theta_lab);
  decay_parts[2].z = cos(decay_parts[2].theta_lab);

  decay_sep_angle_lab = Decay1_LV.Angle(Decay2_LV.Vect());

  //for CM values...
  Decay1_LV.Boost(-boost_V_Parent);
  Decay2_LV.Boost(-boost_V_Parent);
  decay_parts[2].E_tot_cm = Decay2_LV.E();
  decay_parts[2].KE_cm    = decay_parts[2].E_tot_cm - TOTAL_RESTMASS_DECAY2;
  decay_parts[2].ptot_cm  = Decay2_LV.P();
  decay_parts[2].px_cm    = Decay2_LV.Px();
  decay_parts[2].py_cm    = Decay2_LV.Py();
  decay_parts[2].pz_cm    = Decay2_LV.Pz();
  decay_parts[2].theta_cm = Decay2_LV.Theta();
  decay_parts[2].phi_cm   = Decay2_LV.Phi();

  decay_sep_angle_cm = Decay1_LV.Angle(Decay2_LV.Vect());

}

void Decay::Run_Decay_Angles(double iP_KE, double iP_theta, double iP_phi,
			     double d1_theta, double d1_phi,
			     double P_exc, double D1_exc, double D2_exc) {

  if (iP_theta == -1) iP_theta = Random_Parent_Theta();
  if (iP_phi   == -1) iP_phi   = Random_Parent_Phi();

  if (d1_theta == -1) d1_theta = Random_Decay1_Theta_CM();
  if (d1_phi   == -1) d1_phi   = Random_Decay1_Phi_CM();

  if (P_exc  == -1) P_exc  = Random_Parent_Exc_Energy();
  if (D1_exc == -1) D1_exc = Random_Decay1_Exc_Energy();
  if (D2_exc == -1) D2_exc = Random_Decay2_Exc_Energy();

  double TOTAL_RESTMASS_PARENT = RESTMASS_PARENT + P_exc;
  double parent_E = iP_KE + TOTAL_RESTMASS_PARENT;
  double parent_ptot = sqrt(pow(parent_E,2) - pow(TOTAL_RESTMASS_PARENT,2));
  double parent_px = parent_ptot*cos(iP_phi)*sin(iP_theta);
  double parent_py = parent_ptot*sin(iP_phi)*sin(iP_theta);
  double parent_pz = parent_ptot*cos(iP_theta);

  Run_Decay_Vector(parent_px, parent_py, parent_pz, parent_E,
		   d1_theta, d1_phi,
		   D1_exc, D2_exc);

}

bool Decay::GetMasses() {

  for (int i=0; i<3; i++) M_arr[i] = 0;

  ifstream mass_file;

  mass_file.open("/home2/kgh14d/general_sim/mass_info.txt");

  if (!mass_file) {
    cout << "***WHERE IS mass_info.txt???\n";
    return false;
  }

  TString column_labels;
  int Z = 0, A = 0;
  float mass = 0;
  //                 Z                A           Mass_(amu)
  mass_file >> column_labels >> column_labels >> column_labels;

  while (Z != -1) {

    mass_file >> Z >> A >> mass;
    
    for (int i=0; i<3; i++)
      if (Z == Z_arr[i] && A == A_arr[i])
	M_arr[i] = mass*UTOMEV;

  }

  RESTMASS_PARENT = M_arr[0];
  RESTMASS_DECAY1 = M_arr[1];
  RESTMASS_DECAY2 = M_arr[2];

  for (int i=0; i<3; i++)
    if (M_arr[i] == 0)
      return false;

  return true;

}

bool Decay::GetStates() {

  ifstream exc_state_file;
  exc_state_file.open("/home2/kgh14d/general_sim/exc_state_info.txt");
  
  bool found_states[3];
  for (int i=0; i<3; i++)
    found_states[i] = false;

  if (!exc_state_file) {
    cout << "Where's exc_state_info.txt????\n";
    return false;
  }

  int Z = 0, A = 0, num_states = 0;
  double energy = 0, width = 0, frac = 0;
  TString throwaway = "";

  //example format that it expects for header:
  /*
    Z = 4 | A = 8
    Number_of_States: 3
    E_(MeV)    Frac    Width_(MeV)
  */

  while (Z != -1) {

    exc_state_file >> throwaway >> throwaway >> Z >> throwaway
		   >> throwaway >> throwaway >> A
		   >> throwaway >> num_states
		   >> throwaway >> throwaway >> throwaway;

    for (int i=0; i<num_states; i++) {

      exc_state_file >> energy >> width >> frac;

      for (int p=0; p<3; p++) {
	if (Z == Z_arr[p] && A == A_arr[p]) {
	  exc_state_info[p][0][i] = energy;
	  exc_state_info[p][1][i] = width;
	  exc_state_info[p][2][i] = frac;
	  if (!found_states[p]) {//just so it only assigns once
	    num_exc_states[p] = num_states;
	    found_states[p] = true;
	  }
	}
      }

    }

  }

  return (found_states[0] &&
	  found_states[1] &&
	  found_states[2]);

}

void Decay::CheckDecay() {

  //check validity of decay and set decay flag appropriately
  if (DECAY1_A + DECAY2_A != PARENT_A ||
      DECAY1_Z + DECAY2_Z != PARENT_Z) {
    cout << "***WARNING: decay of "
	 << "(" << PARENT_Z << "," << PARENT_A << ") --> "
	 << "(" << DECAY1_Z << "," << DECAY1_A << ") + "
	 << "(" << DECAY2_Z << "," << DECAY2_A << ") "
	 << "not possible.\n";
  }

}

double Decay::Get_Exc_State(float particle_info[3][100], int numstates) {

  //particle_info has three columns: exc state energies, widths, and relative fractions

  double ex_current  = 0,
         current_min = 0,
         current_max = 0;

  double chosen_state = rand->Uniform(0,1); //state selection is random and uniform

  for (int i=0; i<numstates; i++) {

    current_max += particle_info[1][i];

    if (chosen_state >= current_min && chosen_state < current_max)
      return rand->Gaus(particle_info[0][i], particle_info[2][i]/2/2.355); //FWHM = 2.355*sigma

    current_min += particle_info[1][i];

  }

  return 0; //return ground state if we don't match a state,
            //i.e. the relative fractions probably don't add to one

}

void Decay::NucData_Zeros() {

  for (int i=0; i<3; i++) {
    decay_parts[i].px_lab    = 0;
    decay_parts[i].py_lab    = 0;
    decay_parts[i].pz_lab    = 0;
    decay_parts[i].ptot_lab  = 0;
    decay_parts[i].px_cm     = 0;
    decay_parts[i].py_cm     = 0;
    decay_parts[i].pz_cm     = 0;
    decay_parts[i].ptot_cm   = 0;
    decay_parts[i].E_tot_lab = 0;
    decay_parts[i].E_tot_cm  = 0;
    decay_parts[i].KE_lab    = 0;
    decay_parts[i].KE_cm     = 0;
    decay_parts[i].E_ex      = 0;
    decay_parts[i].invmass   = 0;
    decay_parts[i].theta_lab = 0;
    decay_parts[i].phi_lab   = 0;
    decay_parts[i].theta_cm  = 0;
    decay_parts[i].phi_cm    = 0;
    decay_parts[i].x         = 0;
    decay_parts[i].y         = 0;
    decay_parts[i].z         = 0;
  }

}
