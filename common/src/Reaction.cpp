#include "Reaction.h"

using namespace std;

Reaction::Reaction() {

  //called when declared as part of an array of Reaction objects,
  //each of which should then be initialized using Reaction::Initialize(...)

  Rxn_ID       = 0;
  TARGET_Z     = 0;
  TARGET_A     = 0;
  PROJECTILE_Z = 0;
  PROJECTILE_A = 0;
  EJECTILE_Z   = 0;
  EJECTILE_A   = 0;
  RESIDUAL_Z   = 0;
  RESIDUAL_A   = 0;

  CHARGE_TARGET     = 0;
  CHARGE_PROJECTILE = 0;
  CHARGE_EJECTILE   = 0;
  CHARGE_RESIDUAL   = 0;

  for (int k=0; k<4; k++) {
    Z_arr[k] = 0;
    A_arr[k] = 0;
    num_exc_states[k] = 0;
    for (int i=0; i<3; i++) {
      for (int j=0; j<100; j++) {
	exc_state_info[k][i][j] = 0;
      }
    }
  }

  tot_angs = 181;
  angle_weights = new double*[2]; //two columns, one angle, one weight
  for (int i=0; i<2; i++)
    angle_weights[i] = new double[tot_angs]; //default to 0 -- 180 degrees, 1 deg steps
  for (int i=0; i<tot_angs; i++) {
    angle_weights[0][i] = i; //angle
    angle_weights[1][i] = 1./tot_angs; //isotropic probability
  }

  xsec_loaded = false;

  NucData_Zeros(); //initializes rest of variables

  Q_VALUE = 0;

  rand = new TRandom3();
  rand->SetSeed();

}

Reaction::Reaction(int iID,
		   int iTZ, int iTA, int iPZ, int iPA,
		   int iEZ, int iEA, int iRZ, int iRA) {

  //works as follows: if nothing but ID num is given, defaults to proton-proton scattering
  //                  if only target is given, defaults to proton scattering on the target
  //                  if T and P given, defaults to P scattering on T
  //                  if T, P, and E given, runs reaction T(P,E) ; residual is calculated

  Rxn_ID = iID;

  TARGET_Z = iTZ;
  TARGET_A = iTA;

  PROJECTILE_Z = iPZ;
  PROJECTILE_A = iPA;

  if (iEZ == 0 && iEA == 0) { //assume scattering
    EJECTILE_Z = PROJECTILE_Z;
    EJECTILE_A = PROJECTILE_A;
  }
  else {
    EJECTILE_Z = iEZ;
    EJECTILE_A = iEA;
  }

  if (iRA == 0 && iRZ == 0) {
    RESIDUAL_Z = TARGET_Z + PROJECTILE_Z - EJECTILE_Z;
    RESIDUAL_A = TARGET_A + PROJECTILE_A - EJECTILE_A;
  }
  else {
    RESIDUAL_Z = iRZ;
    RESIDUAL_A = iRA;
  }

  CHARGE_TARGET     = TARGET_Z*UNIT_CHARGE;
  CHARGE_PROJECTILE = PROJECTILE_Z*UNIT_CHARGE;
  CHARGE_EJECTILE   = EJECTILE_Z*UNIT_CHARGE;
  CHARGE_RESIDUAL   = RESIDUAL_Z*UNIT_CHARGE;

  Z_arr[0] = TARGET_Z;
  Z_arr[1] = PROJECTILE_Z;
  Z_arr[2] = EJECTILE_Z;
  Z_arr[3] = RESIDUAL_Z;

  A_arr[0] = TARGET_A;
  A_arr[1] = PROJECTILE_A;
  A_arr[2] = EJECTILE_A;
  A_arr[3] = RESIDUAL_A;

  CheckRxn();

  for (int k=0; k<4; k++) {
    num_exc_states[k]=0;
    for (int i=0; i<3; i++) {
      for (int j=0; j<100; j++) {
	exc_state_info[k][i][j] = 0;
      }
    }
  }

  tot_angs = 181;
  angle_weights = new double*[2]; //two columns, one angle, one weight
  for (int i=0; i<2; i++)
    angle_weights[i] = new double[tot_angs]; //default to 0 -- 180 degrees, 1 deg steps
  for (int i=0; i<tot_angs; i++) {
    angle_weights[0][i] = i; //angle
    angle_weights[1][i] = 1./tot_angs; //isotropic probability
  }

  xsec_loaded = false;

  NucData_Zeros(); //initializes rest of variables

  if (!GetMasses()) cout << "***Reaction() warning: problem loading rest masses\n";
  if (!GetStates()) cout << "***Reaction() warning: problem loading excited states\n";

  Q_VALUE = RESTMASS_TARGET + RESTMASS_PROJECTILE - (RESTMASS_EJECTILE + RESTMASS_RESIDUAL);

  rand = new TRandom3();
  rand->SetSeed();  

}

void Reaction::Initialize(int iID,
			  int iTZ, int iTA, int iPZ, int iPA,
			  int iEZ, int iEA, int iRZ, int iRA) {

  //this member function was necessary in order to be able to create an array of Reaction
  //objects, since we technically lack a "default" constructor. In this way, multiple rxns
  //can be studied at once in the main code.

  //works as follows: if nothing but ID num is given, defaults to proton-proton scattering
  //                  if only target is given, defaults to proton scattering on the target
  //                  if T and P given, defaults to P scattering on T
  //                  if T, P, and E given, runs reaction T(P,E) ; residual is calculated

  Rxn_ID = iID;

  TARGET_Z = iTZ;
  TARGET_A = iTA;

  PROJECTILE_Z = iPZ;
  PROJECTILE_A = iPA;

  if (iEZ == 0 && iEA == 0) { //assume scattering
    EJECTILE_Z = PROJECTILE_Z;
    EJECTILE_A = PROJECTILE_A;
  }
  else {
    EJECTILE_Z = iEZ;
    EJECTILE_A = iEA;
  }

  if (iRA == 0 && iRZ == 0) {
    RESIDUAL_Z = TARGET_Z + PROJECTILE_Z - EJECTILE_Z;
    RESIDUAL_A = TARGET_A + PROJECTILE_A - EJECTILE_A;
  }
  else {
    RESIDUAL_Z = iRZ;
    RESIDUAL_A = iRA;
  }

  CHARGE_TARGET     = TARGET_Z*UNIT_CHARGE;
  CHARGE_PROJECTILE = PROJECTILE_Z*UNIT_CHARGE;
  CHARGE_EJECTILE   = EJECTILE_Z*UNIT_CHARGE;
  CHARGE_RESIDUAL   = RESIDUAL_Z*UNIT_CHARGE;

  Z_arr[0] = TARGET_Z;
  Z_arr[1] = PROJECTILE_Z;
  Z_arr[2] = EJECTILE_Z;
  Z_arr[3] = RESIDUAL_Z;

  A_arr[0] = TARGET_A;
  A_arr[1] = PROJECTILE_A;
  A_arr[2] = EJECTILE_A;
  A_arr[3] = RESIDUAL_A;

  CheckRxn();

  for (int k=0; k<4; k++) {
    num_exc_states[k]=0;
    for (int i=0; i<3; i++) {
      for (int j=0; j<100; j++) {
	exc_state_info[k][i][j] = 0;
      }
    }
  }

  if (!GetMasses()) cout << "***Reaction() warning: problem loading rest masses\n";
  if (!GetStates()) cout << "***Reaction() warning: problem loading excited states\n";

  Q_VALUE = RESTMASS_TARGET + RESTMASS_PROJECTILE - (RESTMASS_EJECTILE + RESTMASS_RESIDUAL);

}

Reaction::~Reaction() {

  for (int i=0; i<2; i++)
    delete [] angle_weights[i];
  delete [] angle_weights;

  delete rand;

}

int Reaction::ID() {return Rxn_ID;}

int Reaction::Target_Z()     {return TARGET_Z;}
int Reaction::Projectile_Z() {return PROJECTILE_Z;}
int Reaction::Ejectile_Z()   {return EJECTILE_Z;}
int Reaction::Residual_Z()   {return RESIDUAL_Z;}

int Reaction::Target_A()     {return TARGET_A;}
int Reaction::Projectile_A() {return PROJECTILE_A;}
int Reaction::Ejectile_A()   {return EJECTILE_A;}
int Reaction::Residual_A()   {return RESIDUAL_A;}

double Reaction::Target_Restmass()     {return RESTMASS_TARGET;}
double Reaction::Projectile_Restmass() {return RESTMASS_PROJECTILE;}
double Reaction::Ejectile_Restmass()   {return RESTMASS_EJECTILE;}
double Reaction::Residual_Restmass()   {return RESTMASS_RESIDUAL;}

double Reaction::Target_Charge()     {return CHARGE_TARGET;}
double Reaction::Projectile_Charge() {return CHARGE_PROJECTILE;}
double Reaction::Ejectile_Charge()   {return CHARGE_EJECTILE;}
double Reaction::Residual_Charge()   {return CHARGE_RESIDUAL;}

double Reaction::Q_Value() {return Q_VALUE;}

double Reaction::Rutherford_XSec(double ang_cm) {
  if (sin(ang_cm/2) == 0) return 0;
  double numerator = PROJECTILE_Z*TARGET_Z*FINE_STRUCTURE*HBAR_C;
  double denominator = 4*Entrance_Channel_KEcm()*pow(sin(ang_cm/2),2);
  return pow(numerator/denominator,2)*SQFMTOBARN*1000;
}

/* Kinematics Getters */
double Reaction::Get_px_lab(int ipart_id)
{return (ipart_id > 0 && ipart_id <=4 ) ? rxn_parts[ipart_id - 1].px_lab : 0.0;}
double Reaction::Get_py_lab(int ipart_id)
{return (ipart_id > 0 && ipart_id <=4 ) ? rxn_parts[ipart_id - 1].py_lab : 0.0;}
double Reaction::Get_pz_lab(int ipart_id)
{return (ipart_id > 0 && ipart_id <=4 ) ? rxn_parts[ipart_id - 1].pz_lab : 0.0;}
double Reaction::Get_ptot_lab(int ipart_id)
{return (ipart_id > 0 && ipart_id <=4 ) ? rxn_parts[ipart_id - 1].ptot_lab : 0.0;}
double Reaction::Get_px_cm(int ipart_id)
{return (ipart_id > 0 && ipart_id <=4 ) ? rxn_parts[ipart_id - 1].px_cm : 0.0;}
double Reaction::Get_py_cm(int ipart_id)
{return (ipart_id > 0 && ipart_id <=4 ) ? rxn_parts[ipart_id - 1].py_cm : 0.0;}
double Reaction::Get_pz_cm(int ipart_id)
{return (ipart_id > 0 && ipart_id <=4 ) ? rxn_parts[ipart_id - 1].pz_cm : 0.0;}
double Reaction::Get_ptot_cm(int ipart_id)
{return (ipart_id > 0 && ipart_id <=4 ) ? rxn_parts[ipart_id - 1].ptot_cm : 0.0;}
double Reaction::Get_E_tot_lab(int ipart_id)
{return (ipart_id > 0 && ipart_id <=4 ) ? rxn_parts[ipart_id - 1].E_tot_lab : 0.0;}
double Reaction::Get_E_tot_cm(int ipart_id)
{return (ipart_id > 0 && ipart_id <=4 ) ? rxn_parts[ipart_id - 1].E_tot_cm : 0.0;}
double Reaction::Get_KE_lab(int ipart_id)
{return (ipart_id > 0 && ipart_id <=4 ) ? rxn_parts[ipart_id - 1].KE_lab : 0.0;}
double Reaction::Get_KE_cm(int ipart_id)
{return (ipart_id > 0 && ipart_id <=4 ) ? rxn_parts[ipart_id - 1].KE_cm : 0.0;}
double Reaction::Get_E_ex(int ipart_id)
{return (ipart_id > 0 && ipart_id <=4 ) ? rxn_parts[ipart_id - 1].E_ex : 0.0;}
double Reaction::Get_invmass(int ipart_id)
{return (ipart_id > 0 && ipart_id <=4 ) ? rxn_parts[ipart_id - 1].invmass : 0.0;}
double Reaction::Get_theta_lab(int ipart_id)
{return (ipart_id > 0 && ipart_id <=4 ) ? rxn_parts[ipart_id - 1].theta_lab : 0.0;}
double Reaction::Get_phi_lab(int ipart_id)
{return (ipart_id > 0 && ipart_id <=4 ) ? rxn_parts[ipart_id - 1].phi_lab : 0.0;}
double Reaction::Get_theta_cm(int ipart_id)
{return (ipart_id > 0 && ipart_id <=4 ) ? rxn_parts[ipart_id - 1].theta_cm : 0.0;}
double Reaction::Get_phi_cm(int ipart_id)
{return (ipart_id > 0 && ipart_id <=4 ) ? rxn_parts[ipart_id - 1].phi_cm : 0.0;}
double Reaction::Get_x(int ipart_id)
{return (ipart_id > 0 && ipart_id <=4 ) ? rxn_parts[ipart_id - 1].x : 0.0;}
double Reaction::Get_y(int ipart_id)
{return (ipart_id > 0 && ipart_id <=4 ) ? rxn_parts[ipart_id - 1].y : 0.0;}
double Reaction::Get_z(int ipart_id)
{return (ipart_id > 0 && ipart_id <=4 ) ? rxn_parts[ipart_id - 1].z : 0.0;}

double Reaction::Entrance_Channel_KEcm() {return Get_KEcm(1,2);}
double Reaction::Exit_Channel_KEcm()     {return Get_KEcm(3,4);}

bool Reaction::Load_XSec(char *filename) {

  ifstream xsec_file;
  xsec_file.open(filename);
  if (!xsec_file) return false;

  tot_angs=0; //reinitialize to count new number of angles
  double temp1, temp2; //just during counting phase

  double tot_xsec=0; //to normalize

  while (!(xsec_file.eof())) {
    xsec_file >> temp1 >> temp2;
    tot_angs++;
  }

  //resize angle_weights:
  for (int i=0; i<2; i++) {
    delete [] angle_weights[i];
    angle_weights[i] = new double[tot_angs];
  }

  //have to close and reopen file
  xsec_file.close();
  xsec_file.open(filename);

  for (int ang=0; ang<tot_angs; ang++) {
    xsec_file >> angle_weights[0][ang] >> angle_weights[1][ang];
    tot_xsec += angle_weights[1][ang];
  }

  for (int i=0; i<tot_angs; i++)
    angle_weights[1][i] /= tot_xsec;

  xsec_file.close();

  xsec_loaded = true;

  return xsec_loaded;

}

void Reaction::Run_Kinematics(double Beam_E,
			      double theta_proj, double phi_proj,
			      double targ_exc,   double proj_exc,
			      double ejec_exc,   double res_exc) {

  NucData_Zeros(); //reset data variables

  if (targ_exc == -1) targ_exc = Random_Target_Exc_Energy();
  if (proj_exc == -1) proj_exc = Random_Projectile_Exc_Energy();
  if (ejec_exc == -1) ejec_exc = Random_Ejectile_Exc_Energy();
  if (res_exc  == -1) res_exc  = Random_Residual_Exc_Energy();

  double TOTAL_RESTMASS_TARGET     = RESTMASS_TARGET     + targ_exc;
  double TOTAL_RESTMASS_PROJECTILE = RESTMASS_PROJECTILE + proj_exc;

  TLorentzVector Target_LV(0, 0, 0, RESTMASS_TARGET+targ_exc);

  TLorentzVector Projectile_LV(0,
			       0,
			       sqrt(pow(Beam_E + TOTAL_RESTMASS_PROJECTILE,2)
				           - pow(TOTAL_RESTMASS_PROJECTILE,2)),
			       Beam_E + TOTAL_RESTMASS_PROJECTILE);

  //fill target and proj structs:
  rxn_parts[0].px_lab    = Target_LV.Px();
  rxn_parts[0].py_lab    = Target_LV.Py();
  rxn_parts[0].pz_lab    = Target_LV.Pz();
  rxn_parts[0].ptot_lab  = Target_LV.P();
  rxn_parts[0].E_tot_lab = Target_LV.E();
  rxn_parts[0].KE_lab    = Target_LV.E() - TOTAL_RESTMASS_TARGET;
  rxn_parts[0].E_ex      = Target_LV.M() - RESTMASS_TARGET;
  rxn_parts[0].invmass   = Target_LV.M();
  rxn_parts[0].theta_lab = Target_LV.Theta();
  rxn_parts[0].phi_lab   = Target_LV.Phi();
  rxn_parts[0].theta_cm  = Target_LV.Theta(); //targ at rest
  rxn_parts[0].phi_cm    = Target_LV.Phi();
  rxn_parts[0].x         = 0; //target at origin
  rxn_parts[0].y         = 0;
  rxn_parts[0].z         = 0;

  rxn_parts[1].px_lab    = Projectile_LV.Px();
  rxn_parts[1].py_lab    = Projectile_LV.Py();
  rxn_parts[1].pz_lab    = Projectile_LV.Pz();
  rxn_parts[1].ptot_lab  = Projectile_LV.P();
  rxn_parts[1].E_tot_lab = Projectile_LV.E();
  rxn_parts[1].KE_lab    = Projectile_LV.E() - TOTAL_RESTMASS_PROJECTILE;
  rxn_parts[1].E_ex      = Projectile_LV.M() - RESTMASS_PROJECTILE;
  rxn_parts[1].invmass   = Projectile_LV.M();
  rxn_parts[1].theta_lab = Projectile_LV.Theta();
  rxn_parts[1].phi_lab   = Projectile_LV.Phi();
  rxn_parts[1].theta_cm  = Projectile_LV.Theta(); //both zero for proj anyways
  rxn_parts[1].phi_cm    = Projectile_LV.Phi();
  rxn_parts[1].x         = 0; //reaction happens at origin as well
  rxn_parts[1].y         = 0;
  rxn_parts[1].z         = 0;

  //go into CM frame of "parent" nucleus (target + beam)
  TLorentzVector Parent_LV = Target_LV + Projectile_LV;
  TVector3 boost_V_parent  = Parent_LV.BoostVector();
    
  double parent_E_lab = Parent_LV.E();
  Parent_LV.Boost(-boost_V_parent);
  double parent_E_cm = Parent_LV.E();

  //grab CM values for Target and Projectile
  Target_LV.Boost(-boost_V_parent);
  Projectile_LV.Boost(-boost_V_parent);

  rxn_parts[0].px_cm    = Target_LV.Px();
  rxn_parts[0].py_cm    = Target_LV.Py();
  rxn_parts[0].pz_cm    = Target_LV.Pz();
  rxn_parts[0].ptot_cm  = Target_LV.P();
  rxn_parts[0].E_tot_cm = Target_LV.E();
  rxn_parts[0].KE_cm    = Target_LV.E() - TOTAL_RESTMASS_TARGET;

  rxn_parts[1].px_cm    = Projectile_LV.Px();
  rxn_parts[1].py_cm    = Projectile_LV.Py();
  rxn_parts[1].pz_cm    = Projectile_LV.Pz();
  rxn_parts[1].ptot_cm  = Projectile_LV.P();
  rxn_parts[1].E_tot_cm = Projectile_LV.E();
  rxn_parts[1].KE_cm    = Projectile_LV.E() - TOTAL_RESTMASS_PROJECTILE;

  const double TOTAL_RESTMASS_EJECTILE = RESTMASS_EJECTILE + ejec_exc;
  const double TOTAL_RESTMASS_RESIDUAL = RESTMASS_RESIDUAL + res_exc;

  /* KINEMATICS */

  TLorentzVector Ejectile_LV, Residual_LV;

  //convert input horizontal angles to CM frame  
  theta_proj = LabToCM_Angle(theta_proj, Beam_E, ejec_exc, res_exc);

  rxn_parts[2].theta_cm = theta_proj;
  rxn_parts[2].phi_cm   = phi_proj;
    
  //from E_parent = E_residual + E_ejectile and E^2 = p^2 + M^2 ***in parent frame
  rxn_parts[2].E_tot_cm = (pow(TOTAL_RESTMASS_EJECTILE,2) -
			   pow(TOTAL_RESTMASS_RESIDUAL,2) +
			   pow(parent_E_cm,2)) / (2*parent_E_cm);
  rxn_parts[2].KE_cm   = rxn_parts[2].E_tot_cm - TOTAL_RESTMASS_EJECTILE;
  rxn_parts[2].ptot_cm = sqrt(pow(rxn_parts[2].E_tot_cm,2) - pow(TOTAL_RESTMASS_EJECTILE,2));
  rxn_parts[2].px_cm   = rxn_parts[2].ptot_cm*sin(rxn_parts[2].theta_cm)*cos(rxn_parts[2].phi_cm);
  rxn_parts[2].py_cm   = rxn_parts[2].ptot_cm*sin(rxn_parts[2].theta_cm)*sin(rxn_parts[2].phi_cm);
  rxn_parts[2].pz_cm   = rxn_parts[2].ptot_cm*cos(rxn_parts[2].theta_cm);

  //now create vector for ejectile, and go back to lab frame
  Ejectile_LV.SetPxPyPzE(rxn_parts[2].px_cm,
			 rxn_parts[2].py_cm,
			 rxn_parts[2].pz_cm,
			 rxn_parts[2].E_tot_cm);

  //calculate residual in CM frame (and grab CM values)
  Residual_LV  = Parent_LV - Ejectile_LV;
  rxn_parts[3].px_cm    = Residual_LV.Px();
  rxn_parts[3].py_cm    = Residual_LV.Py();
  rxn_parts[3].pz_cm    = Residual_LV.Pz();
  rxn_parts[3].ptot_cm  = Residual_LV.P();
  rxn_parts[3].E_tot_cm = Residual_LV.E();
  rxn_parts[3].KE_cm    = Residual_LV.E() - TOTAL_RESTMASS_RESIDUAL;
  rxn_parts[3].theta_cm = Residual_LV.Theta();
  rxn_parts[3].phi_cm   = Residual_LV.Phi();

  //now boost everything back to lab frame for rest of parameters
  Parent_LV.Boost(boost_V_parent);
  Residual_LV.Boost(boost_V_parent);
  Ejectile_LV.Boost(boost_V_parent);

  rxn_parts[2].E_tot_lab = Ejectile_LV.E();
  rxn_parts[2].invmass   = Ejectile_LV.M();
  rxn_parts[2].KE_lab    = rxn_parts[2].E_tot_lab - TOTAL_RESTMASS_EJECTILE;
  rxn_parts[2].E_ex      = rxn_parts[2].invmass - RESTMASS_EJECTILE;
  rxn_parts[2].ptot_lab  = Ejectile_LV.P();
  rxn_parts[2].px_lab    = Ejectile_LV.Px();
  rxn_parts[2].py_lab    = Ejectile_LV.Py();
  rxn_parts[2].pz_lab    = Ejectile_LV.Pz();
  rxn_parts[2].theta_lab = Ejectile_LV.Theta();
  rxn_parts[2].phi_lab   = Ejectile_LV.Phi();

  rxn_parts[2].x = cos(rxn_parts[2].phi_lab)*sin(rxn_parts[2].theta_lab);
  rxn_parts[2].y = sin(rxn_parts[2].phi_lab)*sin(rxn_parts[2].theta_lab);
  rxn_parts[2].z = cos(rxn_parts[2].theta_lab);

  rxn_parts[3].E_tot_lab = Residual_LV.E();
  rxn_parts[3].invmass   = Residual_LV.M();
  rxn_parts[3].KE_lab    = rxn_parts[3].E_tot_lab - TOTAL_RESTMASS_RESIDUAL;
  rxn_parts[3].E_ex      = rxn_parts[3].invmass - RESTMASS_RESIDUAL;
  rxn_parts[3].ptot_lab  = Residual_LV.P();
  rxn_parts[3].px_lab    = Residual_LV.Px();
  rxn_parts[3].py_lab    = Residual_LV.Py();
  rxn_parts[3].pz_lab    = Residual_LV.Pz();
  rxn_parts[3].theta_lab = Residual_LV.Theta();
  rxn_parts[3].phi_lab   = Residual_LV.Phi();

  rxn_parts[3].x = cos(rxn_parts[3].phi_lab)*sin(rxn_parts[3].theta_lab);
  rxn_parts[3].y = sin(rxn_parts[3].phi_lab)*sin(rxn_parts[3].theta_lab);
  rxn_parts[3].z = cos(rxn_parts[3].theta_lab);

}

double Reaction::Random_Ejec_Theta_CM(double theta_min, double theta_max) {

  if (theta_max < theta_min) {
    double temptheta = theta_min;
    theta_min = theta_max;
    theta_max = temptheta;
  }

  if (xsec_loaded) { //use diff. xsec as probability distribution

    double current_min = 0,
      current_max = 0;

    double dart = rand->Uniform(theta_min/M_PI,theta_max/M_PI); //throw a "dart" within allowed theta window

    for (int i=0; i<tot_angs; i++) {
      current_max += angle_weights[1][i];
      if (dart >= current_min && dart < current_max) {
	return angle_weights[0][i];
      }
      current_min += angle_weights[1][i];
    }

    return (theta_max + theta_min)/2; //default to average

  }

  else { //assume isotropic

    return rand->Uniform(theta_min,theta_max);

  }

}

double Reaction::Random_Ejec_Theta_Lab(double theta_min, double theta_max, double Beam_E,
				       double ejec_exc, double res_exc) {

  if (ejec_exc == -1) ejec_exc = Random_Ejectile_Exc_Energy();
  if (res_exc  == -1) res_exc  = Random_Residual_Exc_Energy();

  theta_min = LabToCM_Angle(theta_min, Beam_E, ejec_exc, res_exc);
  theta_max = LabToCM_Angle(theta_max, Beam_E, ejec_exc, res_exc);

  return CMToLab_Angle(Random_Ejec_Theta_CM(theta_min, theta_max), Beam_E, ejec_exc, res_exc);

}

double Reaction::Random_Ejec_Phi(double phi_min, double phi_max) {

  if (phi_min > phi_max) {
    double temp = phi_min;
    phi_min = phi_max;
    phi_max = temp;
  }

  return rand->Uniform(phi_min, phi_max);

}

double Reaction::Random_Target_Exc_Energy() {return Random_Exc_Energy(0);}
double Reaction::Random_Projectile_Exc_Energy() {return Random_Exc_Energy(1);}
double Reaction::Random_Ejectile_Exc_Energy() {return Random_Exc_Energy(2);}
double Reaction::Random_Residual_Exc_Energy() {return Random_Exc_Energy(3);}

double Reaction::Random_Exc_Energy(int particleID) {

  return Get_Exc_State(exc_state_info[particleID],num_exc_states[particleID]);

}

double Reaction::LabToCM_Angle(double angle, double Beam_E, double Exc_Ejec, double Exc_Res) {

  //variables for easier readibility
  double
    RP  = RESTMASS_PROJECTILE,
    TRE = RESTMASS_EJECTILE + Exc_Ejec,
    TRR = RESTMASS_RESIDUAL + Exc_Res,
    Q   = Q_VALUE,
    E   = Beam_E;

  //this one I really had to twist Illiadis around...he gives the formula for CM-->lab conversion,
  //which isn't so easily invertible. Trust me, this works, even if it isn't exactly obvious.
  double velocity_ratio = sqrt(RP*TRE*E/((TRR*(TRR + TRE))*Q + TRR*(TRR + TRE - RP)*E));
  velocity_ratio = 1/(cos(angle) + sqrt(1/(velocity_ratio*velocity_ratio) - sin(angle)*sin(angle)));

  return atan2(sin(angle),cos(angle) - velocity_ratio);

}

double Reaction::CMToLab_Angle(double angle, double Beam_E, double Exc_Ejec, double Exc_Res) {

  double
    RP  = RESTMASS_PROJECTILE,
    TRE = RESTMASS_EJECTILE + Exc_Ejec,
    TRR = RESTMASS_RESIDUAL + Exc_Res,
    Q   = Q_VALUE,
    E   = Beam_E;

  double velocity_ratio = sqrt(RP*TRE*E/((TRR*(TRR + TRE))*Q + TRR*(TRR + TRE - RP)*E));

  return atan2(sin(angle),cos(angle) + velocity_ratio);
  
}

double Reaction::LabToCM_Factor(double angle_cm, double Beam_E, double Exc_Ejec, double Exc_Res) {

  double
    RP  = RESTMASS_PROJECTILE,
    TRE = RESTMASS_EJECTILE + Exc_Ejec,
    TRR = RESTMASS_RESIDUAL + Exc_Res,
    Q   = Q_VALUE,
    E   = Beam_E;

  double velocity_ratio = sqrt(RP*TRE*E/((TRR*(TRR + TRE))*Q + TRR*(TRR + TRE - RP)*E));

  double numerator = 1 + velocity_ratio*cos(angle_cm);
  double denominator = pow(sin(angle_cm),2) + pow(cos(angle_cm) + velocity_ratio,2);
  denominator = pow(denominator,1.5);
  
  return numerator/denominator;

}

void Reaction::Print() {

  cout << "********************************************\n"
       << "REACTION: "
       << "(" << Target_Z()     << "," << Target_A()     << ") + "
       << "(" << Projectile_Z() << "," << Projectile_A() << ") --> "
       << "(" << Ejectile_Z()   << "," << Ejectile_A()   << ") + "
       << "(" << Residual_Z()   << "," << Residual_A()   << endl;

  for (int i=0; i<4; i++) {

    cout << "------------------------------------------\n";

    switch(i) {
    case 0: cout << "TARGET -- ("     << Target_Z()     << "," << Target_A()     << "):\n"; break;
    case 1: cout << "PROJECTILE -- (" << Projectile_Z() << "," << Projectile_A() << "):\n"; break;
    case 2: cout << "EJECTILE -- ("   << Ejectile_Z()   << "," << Ejectile_A()   << "):\n"; break;
    case 3: cout << "RESIDUAL -- ("   << Residual_Z()   << "," << Residual_A()   << "):\n"; break;
    }

    cout << "\tPx_lab = "    << Get_px_lab(i+1)             << " MeV\n"
	 << "\tPy_lab = "    << Get_py_lab(i+1)             << " MeV\n"
	 << "\tPz_lab = "    << Get_pz_lab(i+1)             << " MeV\n"
	 << "\tPtot_lab = "  << Get_ptot_lab(i+1)           << " MeV\n"
	 << "\tPx_cm = "     << Get_px_cm(i+1)              << " MeV\n"
	 << "\tPy_cm = "     << Get_py_cm(i+1)              << " MeV\n"
	 << "\tPz_cm = "     << Get_pz_cm(i+1)              << " MeV\n"
	 << "\tPtot_cm = "   << Get_ptot_cm(i+1)            << " MeV\n"
	 << "\tEtot_lab = "  << Get_E_tot_lab(i+1)          << " MeV\n"
	 << "\tEtot_cm = "   << Get_E_tot_cm(i+1)           << " MeV\n"
	 << "\tKE_lab = "    << Get_KE_lab(i+1)             << " MeV\n"
	 << "\tKE_cm = "     << Get_KE_cm(i+1)              << " MeV\n"
	 << "\tE_exc = "     << Get_E_ex(i+1)               << " MeV\n"
	 << "\ttheta_lab = " << Get_theta_lab(i+1)*RADTODEG << " degrees\n"
	 << "\ttheta_cm = "  << Get_theta_cm(i+1)*RADTODEG  << " degrees\n"
	 << "\tphi_lab = "   << Get_phi_lab(i+1)*RADTODEG   << " degrees\n"
	 << "\tphi_cm = "    << Get_phi_cm(i+1)*RADTODEG    << " degrees\n";
  }

  cout << "------------------------------------------\n"
       << "Entrance Channel KEcm = " << Entrance_Channel_KEcm() << " MeV\n"
       << "Exit Channel KEcm = "     << Exit_Channel_KEcm()     << " MeV\n"
       << "Rutherford_XSec (CM) = "  << Rutherford_XSec(Get_theta_cm(3)) << " mb/sr\n"
       << "********************************************\n";

}

void Reaction::Reset() {NucData_Zeros();}

double Reaction::Get_KEcm(int p1, int p2) {

  return Get_KE_cm(1) + Get_KE_cm(2);

}

bool Reaction::GetMasses() {

  for (int i=0; i<4; i++) M_arr[i] = 0;

  ifstream mass_file;

  mass_file.open("/home/kghman/postgrad_work_2023/general_sim/mass_info.txt");

  if (!mass_file) {
    cout << "***WHERE IS mass_info.txt???\n";
    return false;
  }

  string column_labels;
  int Z = 0, A = 0;
  float mass = 0;
  //                 Z                A           Mass_(amu)
  mass_file >> column_labels >> column_labels >> column_labels;

  while (Z != -1) {

    mass_file >> Z >> A >> mass;
    
    for (int i=0; i<4; i++)
      if (Z == Z_arr[i] && A == A_arr[i])
	M_arr[i] = (mass - Z*RESTMASS_ELECTRON)*UTOMEV;

  }

  RESTMASS_TARGET     = M_arr[0];
  RESTMASS_PROJECTILE = M_arr[1];
  RESTMASS_EJECTILE   = M_arr[2];
  RESTMASS_RESIDUAL   = M_arr[3];

  for (int i=0; i<4; i++)
    if (M_arr[i] == 0)
      return false;

  return true;

}

bool Reaction::GetStates() {

  ifstream exc_state_file;
  exc_state_file.open("/home/kghman/postgrad_work_2023/general_sim/exc_state_info.txt");

  bool found_states[4];
  for (int i=0; i<4; i++)
    found_states[i] = false;

  if (!exc_state_file) {
    cout << "Where's exc_state_info.txt????\n";
    return false;
  }

  int Z = 0, A = 0, num_states = 0;
  double energy = 0, width = 0, frac = 0;
  string throwaway = "";

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

      for (int s=0; s<4; s++) {
	if (Z == Z_arr[s] && A == A_arr[s]) {
	  exc_state_info[s][0][i] = energy;
	  exc_state_info[s][1][i] = width;
	  exc_state_info[s][2][i] = frac;
	  if (!found_states[s]) {
	    num_exc_states[s] = num_states; //just so it only assigns once
	    found_states[s] = true;
	  }
	}
      }

    }

  }

  return (found_states[0] &&
	  found_states[1] &&
	  found_states[2] &&
	  found_states[3]);

}

void Reaction::CheckRxn() {

  //check validity of reaction and warn user
  if (TARGET_Z + PROJECTILE_Z != EJECTILE_Z + RESIDUAL_Z ||
      TARGET_A + PROJECTILE_A != EJECTILE_A + RESIDUAL_A) {
    cout << "***WARNING: reaction "
	 << "(" << TARGET_Z     << "," << TARGET_A     << ") + "
	 << "(" << PROJECTILE_Z << "," << PROJECTILE_A << ") --> "
	 << "(" << EJECTILE_Z   << "," << EJECTILE_A   << ") + "
	 << "(" << RESIDUAL_Z   << "," << RESIDUAL_A   << ") not possible. Exiting.\n";
  }

}

double Reaction::Get_Exc_State(float particle_info[3][100], int numstates) {

  //particle_info has three columns: exc state energies, widths, and relative fractions

  double ex_current  = 0,
         current_min = 0,
         current_max = 0;

  double chosen_state =  rand->Uniform(0,1); //state selection is random and uniform

  double mean, sigma;

  for (int i=0; i<numstates; i++) {

    mean = particle_info[0][i]; //centroid energy
    sigma = particle_info[2][i]/2/2.355; //FWHM = 2.355*sigma

    current_max += particle_info[1][i];

    if (chosen_state >= current_min && chosen_state < current_max)
      return rand->Gaus(mean, sigma);

    current_min += particle_info[1][i];

  }

  return 0; //return ground state if we don't match a state,
            //i.e. the relative fractions probably don't add to one

}

void Reaction::NucData_Zeros() {

  for (int i=0; i<4; i++) {
    rxn_parts[i].px_lab    = 0;
    rxn_parts[i].py_lab    = 0;
    rxn_parts[i].pz_lab    = 0;
    rxn_parts[i].ptot_lab  = 0;
    rxn_parts[i].px_cm     = 0;
    rxn_parts[i].py_cm     = 0;
    rxn_parts[i].pz_cm     = 0;
    rxn_parts[i].ptot_cm   = 0;
    rxn_parts[i].E_tot_lab = 0;
    rxn_parts[i].E_tot_cm  = 0;
    rxn_parts[i].KE_lab    = 0;
    rxn_parts[i].KE_cm     = 0;
    rxn_parts[i].E_ex      = 0;
    rxn_parts[i].invmass   = 0;
    rxn_parts[i].theta_lab = 0;
    rxn_parts[i].phi_lab   = 0;
    rxn_parts[i].theta_cm  = 0;
    rxn_parts[i].phi_cm    = 0;
    rxn_parts[i].x         = 0;
    rxn_parts[i].y         = 0;
    rxn_parts[i].z         = 0;
  }

}
