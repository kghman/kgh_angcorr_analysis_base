#ifndef __nucleusdata_h
#define __nucleusdata_h

struct NucleusData { //to hold all relevant data about a given participant in the reaction
  double px_lab, py_lab, pz_lab, ptot_lab, 
    px_cm, py_cm, pz_cm, ptot_cm, 
    E_tot_lab, E_tot_cm, KE_lab, KE_cm, E_ex, invmass, 
    theta_cm, theta_lab, phi_cm, phi_lab,
    x, y, z;
};

#endif
