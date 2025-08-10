#ifndef __S2_DETECTOR_H
#define __S2_DETECTOR_H

//all angle inputs and outputs in radians
//distances can be whatever, PROVIDED that
//all input dists are the same

#include <TRandom3.h>
#include <cstdlib>

struct pixeldata {

  double min_theta;
  double max_theta;
  double min_costheta;
  double max_costheta;
  double min_phi;
  double max_phi;

};

class S2Detector {

 public:

  S2Detector(int nr, int nw, double inr, double outr, double z, double phirot=0);
  ~S2Detector();

  int Get_Num_Rings()  {return num_rings;}
  int Get_Num_Wedges() {return num_wedges;}

  double Get_Inner_Radius() {return inner_radius;}
  double Get_Outer_Radius() {return outer_radius;}

  double Get_Z_Dist()  {return z_dist;}
  double Get_Phi_Rot() {return phi_rot;}

  double Get_DTheta(int r) {return dtheta[r];}
  double Get_DPhi()   {return dphi;}

  std::pair<double,double> GetThetaPhi(int iringch, int iwedgech);
  std::pair<int,int> GetRingWedge(double itheta, double iphi);

  std::pair<double,double> GetCosThetaPhi(int iringch, int iwedgech);
  std::pair<int,int> GetRingWedgeCos(double icostheta, double iphi);

 private:

  int num_rings;
  int num_wedges;

  double inner_radius; //cm
  double outer_radius; //cm
  double dradius; //cm

  double inner_theta; //rad
  double outer_theta; //rad

  double inner_costheta;
  double outer_costheta;

  double z_dist; //z-dist from target, cm
  double phi_rot; //rotation of channel zero (at phi=0 by default)

  double *dtheta; //per ring (not uniform)
  double *dcostheta;
  double dphi; //per wedge (uniform)

  pixeldata **pixel;

  TRandom3 *rand;

  bool ValidRing(int r);
  bool ValidWedge(int w);

};

#endif
