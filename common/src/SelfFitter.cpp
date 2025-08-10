#include "SelfFitter.h"
#include "mathfuncs.h"
#include <iostream>

using namespace std;

SelfFitter::SelfFitter() {

  num_states = 0;
  num_params = 3*num_states+1+2;
  toggle_interference = 1;
  toggle_relativistic = 0;

  phase = 0;
  phase_min = 0;
  phase_max = 0;

  bg_e = 0;
  bg_e_min = 0;
  bg_e_max = 0;

  bg_slope = 0;
  bg_slope_min = 0;
  bg_slope_max = 0;

  TI = 0;
  SI = 0;

  width_perc = 0;

  param_arr = 0;

  fitmin = 0;
  fitmax = 20;
  graphnum = 5000;

  bw_fit = 0;
  update_fitfunc();

}

SelfFitter::~SelfFitter() {if (param_arr!=0) delete [] param_arr;}

bool SelfFitter::Add_State(double iM, double iFWHM, double iscale, bool iis_BW, bool ifit_M, bool ifit_G) {

  num_states++;
  num_params = 3*num_states+1+2;

  fac.resize(num_states);
  M.resize(num_states);
  G.resize(num_states);
  is_BW.resize(num_states);
  fit_M.resize(num_states);
  fit_G.resize(num_states);

  min_G.resize(num_states);
  max_G.resize(num_states);

  fac[num_states-1] = iscale;
  M[num_states-1] = iM;
  G[num_states-1] = iFWHM;
  is_BW[num_states-1] = iis_BW;
  fit_M[num_states-1] = ifit_M;
  fit_G[num_states-1] = ifit_G;

  min_G[num_states-1] = iFWHM*(1-width_perc);
  max_G[num_states-1] = iFWHM*(1+width_perc);

  if (param_arr!=0) delete [] param_arr;
  param_arr = new double[num_params];
  sort_params();

  update_fitfunc();

}

bool SelfFitter::Set_Phase(double iphase, double iphase_min, double iphase_max) {

  phase = iphase;
  phase_min = iphase_min;
  phase_max = iphase_max;

}

bool SelfFitter::Set_Background(double iE_int, double islope,
				double iE_int_min, double iE_int_max,
				double islope_min, double islope_max) {

  bg_e = iE_int;
  bg_e_min = iE_int_min;
  bg_e_max = iE_int_max;
  
  bg_slope = islope;
  bg_slope_min = islope_min;
  bg_slope_max = islope_max;

}

bool SelfFitter::Set_Width_Perc(double perc) {width_perc = abs(perc);}
bool SelfFitter::Set_TI(int iTI) {TI = (iTI < num_states ? abs(iTI) : 0);}
bool SelfFitter::Set_SI(int iSI) {SI = (iSI < num_states ? abs(iSI) : 0);}

bool SelfFitter::Set_Fit_Min(double ifitmin) {fitmin = ifitmin;}
bool SelfFitter::Set_Fit_Max(double ifitmax) {fitmax = ifitmax;}
bool SelfFitter::Set_Fit_Bounds(double ifitmin, double ifitmax) {
  Set_Fit_Min(ifitmin);
  Set_Fit_Max(ifitmax);
}

int SelfFitter::Get_Num_States() {return num_states;}
int SelfFitter::Get_Num_Params() {return num_params;}

double SelfFitter::Get_Set_Fac(int istate) {return (!valid_state(istate) ? fac[istate] : 0);}
double SelfFitter::Get_Set_M(int istate)   {return (!valid_state(istate) ? M[istate] : 0);}
double SelfFitter::Get_Set_G(int istate)   {return (!valid_state(istate) ? G[istate] : 0);}
double SelfFitter::Get_Set_Phase()         {return phase;}
double SelfFitter::Get_Set_Phase_Min()     {return phase_min;}
double SelfFitter::Get_Set_Phase_Max()     {return phase_max;}
double SelfFitter::Get_Set_BG_E()          {return bg_e;}
double SelfFitter::Get_Set_BG_E_Min()      {return bg_e_min;}
double SelfFitter::Get_Set_BG_E_Max()      {return bg_e_max;}
double SelfFitter::Get_Set_BG_Slope()      {return bg_slope;}
double SelfFitter::Get_Set_BG_Slope_Min()  {return bg_slope_min;}
double SelfFitter::Get_Set_BG_Slope_Max()  {return bg_slope_max;}
double SelfFitter::Get_Set_Width_Perc()    {return width_perc;}
int    SelfFitter::Get_Set_TI()            {return TI;}
int    SelfFitter::Get_Set_SI()            {return SI;}
double SelfFitter::Get_Set_Fit_Min()       {return fitmin;}
double SelfFitter::Get_Set_Fit_Max()       {return fitmax;}

bool SelfFitter::Fit(TH1F* fithist, double ifitmin, double ifitmax) {
 
  if (num_states == 0) return false;
  if (ifitmax < ifitmin) return false;

  double fitmin_use=0, fitmax_use=0;

  if (ifitmin < 0 || ifitmax < 0) {
    fitmin_use = fitmin;
    fitmax_use = fitmax;
  }
  else {
    fitmin_use = ifitmin;
    fitmax_use = ifitmax;
  }

  fithist->Fit(bw_fit,"R","",fitmin_use,fitmax_use);

}

double SelfFitter::Get_Fit_Fac(int istate, int i) {return Get_Fit_State_Param(istate,0,i);}
double SelfFitter::Get_Fit_M(int istate, int i)   {return Get_Fit_State_Param(istate,1,i);}
double SelfFitter::Get_Fit_G(int istate, int i)   {return Get_Fit_State_Param(istate,2,i);}
double SelfFitter::Get_Fit_State_Param(int istate, int np, int i) {
  if(!valid_state(istate)) return 0;
  if (i==1) return bw_fit->GetParError(istate*3+np);
  else      return bw_fit->GetParameter(istate*3+np);
}
double SelfFitter::Get_Fit_Phase(int i) {
  if (i==1) return bw_fit->GetParError(num_params-3);
  else      return bw_fit->GetParameter(num_params-3);
}
double SelfFitter::Get_Fit_BG_E(int i) {
  if (i==1) return bw_fit->GetParError(num_params-2);
  else      return bw_fit->GetParameter(num_params-2);
}
double SelfFitter::Get_Fit_BG_Slope(int i) {
  if (i==1) return bw_fit->GetParError(num_params-1);
  else      return bw_fit->GetParameter(num_params-1);
}

double SelfFitter::Get_Fit_ChiSq() {return bw_fit->GetChisquare();}
int SelfFitter::Get_Fit_NDF()      {return bw_fit->GetNDF();}

TF1* SelfFitter::Get_Fit_TF1() {return bw_fit;}

double SelfFitter::operator()(double *var_arr, double *param_arr) {

  double return_val = 0;
  double bg_val = 0;

  double E = var_arr[0];

  double fac_op[num_states];
  double M_op[num_states];
  double G_op[num_states];
  double phase_op;
  double bg_e_op;
  double bg_slope_op;

  for (int i=0; i<num_states; i++) {
    fac_op[i] = param_arr[3*i];
    M_op[i]   = param_arr[3*i+1];
    G_op[i]   = param_arr[3*i+2];
    if      (is_BW[i] &&  toggle_relativistic) return_val += Breit_Wigner_Rel(E,M_op[i],G_op[i],fac_op[i]);
    else if (is_BW[i] && !toggle_relativistic) return_val += Breit_Wigner(E,M_op[i],G_op[i],fac_op[i]);
    else                                       return_val += Gaussian(E,M_op[i],G_op[i],fac_op[i]);
  }
  phase_op    = param_arr[num_params-3];
  bg_e_op     = param_arr[num_params-2];
  bg_slope_op = param_arr[num_params-1];
  bg_val = bg_slope_op*(E - bg_e_op);
  if (bg_val >= 0) return_val += bg_val;

  //interference term specifically for superradiance-paired states (TI and SI indices set in constructor)
  double X=0, Y=0, kT=0, kS=0, F=0;
   
  if (!toggle_relativistic) {
    //nonrelativistic form:
    X = (E-M_op[TI])*(E-M_op[SI]) + G_op[TI]*G_op[SI]/4;
    Y = (G_op[TI]/2)*(E-M_op[SI]) - (G_op[SI]/2)*(E-M_op[TI]);
    F = 2*sqrt(fac_op[TI]*fac_op[SI]*G_op[TI]*G_op[SI])/(2*M_PI)*(X*cos(phase_op) + Y*sin(phase_op))/(X*X+Y*Y);
  }
  else {
    //relativistic form:
    X = (E*E-M_op[TI]*M_op[TI])*(E*E-M_op[SI]*M_op[SI]) + M_op[TI]*M_op[SI]*G_op[TI]*G_op[SI];
    Y = M_op[TI]*G_op[TI]*(E*E-M_op[SI]*M_op[SI]) - M_op[SI]*G_op[SI]*(E*E-M_op[TI]*M_op[TI]);
    kT = k_fac(M_op[TI],G_op[TI]);
    kS = k_fac(M_op[SI],G_op[SI]);
    F = 2*sqrt(fac_op[TI]*fac_op[SI]*kT*kS)*(X*cos(phase_op) + Y*sin(phase_op))/(X*X+Y*Y);
  }

  if (toggle_interference) return_val += F;

  return return_val;

}

bool SelfFitter::Toggle_Interference() {toggle_interference = !(toggle_interference);}
bool SelfFitter::Toggle_Relativistic() {toggle_relativistic = !(toggle_relativistic);}

bool SelfFitter::valid_state(int istate) {

  if (istate >= 0 && istate < num_states) return true;
  return false;

}

bool SelfFitter::sort_params() {

  for (int i=0; i<num_params-3; i++) {
    switch (i%3) {
    case 0: param_arr[i] = fac[i/3]; break;
    case 1: param_arr[i] = M[i/3]; break;
    case 2: param_arr[i] = G[i/3]; break;
    }
  }
  param_arr[num_params-3] = phase;
  param_arr[num_params-2] = bg_e;
  param_arr[num_params-1] = bg_slope;

}

bool SelfFitter::unsort_params() {

  for (int i=0; i<num_params-3; i++) {
    switch (i%3) {
    case 0: fac_fit[i/3][0] = bw_fit->GetParameter(i); fac_fit[i/3][1] = bw_fit->GetParError(i); break;
    case 1: M_fit[i/3][0]   = bw_fit->GetParameter(i); M_fit[i/3][1]   = bw_fit->GetParError(i); break;
    case 2: G_fit[i/3][0]   = bw_fit->GetParameter(i); G_fit[i/3][1]   = bw_fit->GetParError(i); break;
    }
  }
  phase_fit[0] = bw_fit->GetParameter(num_params-3);
  phase_fit[1] = bw_fit->GetParError(num_params-3);
  bg_e_fit[0] = bw_fit->GetParameter(num_params-2);
  bg_e_fit[1] = bw_fit->GetParError(num_params-2);
  bg_slope_fit[0] = bw_fit->GetParameter(num_params-1);
  bg_slope_fit[1] = bw_fit->GetParError(num_params-1);

}

bool SelfFitter::update_fitfunc() {

  if (bw_fit != 0) delete bw_fit;

  bw_fit = new TF1("bw_fit",this,fitmin,fitmax,num_params);

  if (num_states > 0) {
    bw_fit->SetParameters(param_arr);
    bw_fit->SetNumberFitPoints(graphnum);
    bw_fit->SetNpx(graphnum);
    for (int i=0; i<num_params-3; i++) {
      switch (i%3) { //all BW params must be positive
      case 0: bw_fit->SetParLimits(i,0,1E6); break; //yield
      case 1: bw_fit->SetParLimits(i,fitmin,fitmax); //centroid mass
	if(!fit_M[i/3]) bw_fit->FixParameter(i,M[i/3]);
	break; 
      case 2: bw_fit->SetParLimits(i,min_G[i/3],max_G[i/3]);
	if(!fit_G[i/3]) bw_fit->FixParameter(i,G[i/3]);
	break; //width
      }
    }
    bw_fit->SetParLimits(num_params-3,phase_min,phase_max);
    bw_fit->SetParLimits(num_params-2,bg_e_min,bg_e_max);
    bw_fit->SetParLimits(num_params-1,bg_slope_min,bg_slope_max);
  }

}
