#include "SpherDet.h"

#include <cmath>

using namespace std;

SpherDet::SpherDet() {

  //defaults to full 4pi
  theta_min = 0;
  theta_max = M_PI;
  phi_min   = 0;
  phi_max   = 2*M_PI;

}

SpherDet::SpherDet(double itmin, double itmax, double ipmin, double ipmax) {

  if (itmin > itmax) {
    double swap = itmin;
    itmin = itmax;
    itmax = swap;
  }
  if (ipmin > ipmax) {
    double swap = ipmin;
    ipmin = ipmax;
    ipmax = swap;
  }

  theta_min = itmin;
  theta_max = itmax;
  phi_min   = ipmin;
  phi_max   = ipmax;

}

SpherDet::~SpherDet() {}

bool SpherDet::SetThetaRange(double itmin, double itmax) {

  if (itmin > itmax) {
    double swap = itmin;
    itmin = itmax;
    itmax = swap;
  }

  theta_min = itmin;
  theta_max = itmax;

}

bool SpherDet::SetPhiRange(double ipmin, double ipmax) {

  if (ipmin > ipmax) {
    double swap = ipmin;
    ipmin = ipmax;
    ipmax = swap;
  }

  theta_min = ipmin;
  theta_max = ipmax;

}

bool SpherDet::Hit(double itheta, double iphi) {

   if (itheta > theta_min && itheta < theta_max &&
         iphi > phi_min   &&   iphi < phi_max)
    return true;

  return false;

}
