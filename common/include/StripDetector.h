#ifndef __STRIP_DETECTOR_H
#define __STRIP_DETECTOR_H

// +z is along beam axis
// +y is vertically "upward" in the lab frame

//angles must be in radians, but distances can be whatever
//PROVIDED all input distances are the same

#include <TRandom3.h>
#include <cstdlib>
#include <cmath>

struct channeldata {

  double min_x, max_x; //min-->max in direction of phi
  double min_y, max_y; //same
  double min_z, max_z; //min=upstream, max=downstream

};

class StripDetector {

 public:
  
  StripDetector(int ns, int nb, double len, double wid, double cphi, double cz, double crho);
  ~StripDetector();

  std::pair<double,double> GetThetaPhi(int istripch, double istrip_ratio, double istrip_unc=-1);
  std::pair<double,double> GetCosThetaPhi(int istripch, double istrip_ratio, double istrip_unc=-1);
  
  std::pair<int,double> GetChannelRatio(double itheta, double iphi);
  std::pair<int,double> GetChannelRatioCos(double icostheta, double iphi);
  
  std::pair<double,double> GetThetaPhiBack(int istripch, int ibackch);
  std::pair<double,double> GetCosThetaPhiBack(int istripch, int ibackch);
  
  std::pair<int,double> GetStripBack(double itheta, double iphi);
  std::pair<int,double> GetStripBackCos(double icostheta, double iphi);

 private:

  const double STRIP_UNC_DEF; //% of strip length (default which can be overwritten by the functions)
  
  int num_strips;
  int num_backs;

  double length; //common to all strips, hence total
  double back_length; //assuming all equal
  double total_width;
  double strip_width; //assuming equal widths

  double center_phi; //assuming det centered above x-axis (corresponds to zero phi)
  double center_z;
  double center_rho; //perpendicular radius from axis

  channeldata *strip, **pixel; //pixel uses the back-front overlaps instead of strip ratio

  TRandom3 *rand;

  bool ValidChannel(int f);
  bool ValidBack(int b);
  bool ValidRatio(double r);

};

#endif
