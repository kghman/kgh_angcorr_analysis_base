/* Simple class to specify a geometry and check for hits within it */

#ifndef SPHERDET_H
#define SPHERDET_H

#include "constants.h"

class SpherDet {

 public:

  SpherDet();
  SpherDet(double, double, double, double); //th min,max and phi min,max
  ~SpherDet();

  bool SetThetaRange(double, double);
  bool SetPhiRange(double, double);

  bool Hit(double, double);
  
 private:

  double
    theta_min,
    theta_max,
    phi_min,
    phi_max;

};

#endif
