#ifndef __SELFFITTER_H
#define __SELFFITTER_H

#include <cmath>
#include <vector>

#include <TH1.h>
#include <TF1.h>

class SelfFitter {

 public:

  SelfFitter();
  ~SelfFitter();

  bool Add_State(double iM, double iFWHM, double iscale=1, bool iis_BW=0, bool ifit_M=1, bool ifit_G=1);
  bool Set_Phase(double iphase, double iphase_min=0, double iphase_max=2*M_PI);
  bool Set_Background(double iE_int, double islope,
		      double iE_int_min=0, double iE_int_max=20,
		      double islope_min=0, double islope_max=100);
  bool Set_Width_Perc(double perc);
  bool Set_TI(int iTI);
  bool Set_SI(int iSI);

  bool Set_Fit_Min(double ifitmin);
  bool Set_Fit_Max(double ifitmax);
  bool Set_Fit_Bounds(double ifitmin, double ifitmax);

  int Get_Num_States();
  int Get_Num_Params();

  double Get_Set_Fac(int istate);
  double Get_Set_M(int istate);
  double Get_Set_G(int istate);
  double Get_Set_Phase();
  double Get_Set_Phase_Min();
  double Get_Set_Phase_Max();
  double Get_Set_BG_E();
  double Get_Set_BG_E_Min();
  double Get_Set_BG_E_Max();
  double Get_Set_BG_Slope();
  double Get_Set_BG_Slope_Min();
  double Get_Set_BG_Slope_Max();
  double Get_Set_Width_Perc();
  int    Get_Set_TI();
  int    Get_Set_SI();
  double Get_Set_Fit_Min();
  double Get_Set_Fit_Max();

  bool Fit(TH1F*, double ifitmin=-1, double ifitmax=-1);

  //has option (i=1) to grab error instead
  double Get_Fit_Fac(int istate, int i=0);
  double Get_Fit_M(int istate, int i=0);
  double Get_Fit_G(int istate, int i=0);
  double Get_Fit_State_Param(int istate, int np, int i=0); //np=0 is fac, 1 is M, 2 is G
  double Get_Fit_Phase(int i=0);
  double Get_Fit_BG_E(int i=0);
  double Get_Fit_BG_Slope(int i=0);

  double Get_Fit_ChiSq();
  int Get_Fit_NDF();

  TF1* Get_Fit_TF1();

  double operator()(double *var_arr, double* param_arr);

  bool Toggle_Interference();
  bool Toggle_Relativistic();

 private:

  int num_states; //defaults at zero, then increments as states are added via member functions
  int num_params;
  bool toggle_interference; //default is yes
  bool toggle_relativistic; //default is no

  double phase, phase_min, phase_max; //as multiple of pi/2
  double bg_e, bg_e_min, bg_e_max;
  double bg_slope, bg_slope_min, bg_slope_max;
  std::vector<double> fac, M, G;
  std::vector<bool> is_BW, fit_M, fit_G;
  bool TI, SI; //trapped and superradiant indices

  double *param_arr;

  //two slots: one for error
  std::vector<std::vector<double>> fac_fit;
  std::vector<std::vector<double>> M_fit;
  std::vector<std::vector<double>> G_fit;
  double phase_fit[2];
  double bg_e_fit[2];
  double bg_slope_fit[2];

  double width_perc;
  std::vector<double> min_G, max_G;

  double fitmin, fitmax;
  int graphnum;

  TF1 *bw_fit;

  bool valid_state(int istate);
  bool resize(bool*, int old_size, int new_size);
  bool resize(double*, int old_size, int new_size);
  bool sort_params();
  bool unsort_params();
  bool update_fitfunc();

};

#endif
