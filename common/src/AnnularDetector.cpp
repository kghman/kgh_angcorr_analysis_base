#include "AnnularDetector.h"

#include <cmath>

using namespace std;

AnnularDetector::AnnularDetector(int nr, int nw, double inr, double outr, double z, double phirot) {

  num_rings = nr;
  num_wedges = nw;
  
  inner_radius = inr;
  outer_radius = outr;
  dradius = (outer_radius-inner_radius)/num_rings;

  z_dist = z;
  phi_rot = phirot;

  //can go ahead and do math just from these quantities
  inner_theta = atan2(inner_radius,z_dist);
  outer_theta = atan2(outer_radius,z_dist);

  if (inner_theta < 0) inner_theta += M_PI;
  if (outer_theta < 0) outer_theta += M_PI;

  inner_costheta = cos(inner_theta);
  outer_costheta = cos(outer_theta);

  ring_radius = new double[num_rings];
  dtheta = new double[num_rings];
  dcostheta = new double[num_rings];
  for (int r=0; r<num_rings; r++) {
    ring_radius[r] = (r==0?inner_radius+dradius/2:ring_radius[r-1]+dradius);
    dtheta[r] = abs(atan2(inner_radius+(r+1)*dradius,z_dist) - atan2(inner_radius+r*dradius,z_dist));
    dcostheta[r] = abs(cos(atan2(inner_radius+(r+1)*dradius,z_dist)) - cos(atan2(inner_radius+r*dradius,z_dist)));
  }

  dphi = 2*M_PI/num_wedges;

  pixel = new pixeldata*[num_rings];
  for (int i=0; i<num_rings; i++)
    pixel[i] = new pixeldata[num_wedges];

  for (int r=0; r<num_rings; r++) {
    for (int w=0; w<num_wedges; w++) {
      pixel[r][w].min_radius = (r==0?inner_radius:pixel[r-1][w].max_radius);
      pixel[r][w].max_radius = pixel[r][w].min_radius + dradius;
      if (z_dist >= 0) {
	pixel[r][w].min_theta = (r==0?inner_theta:pixel[r-1][w].max_theta);
	pixel[r][w].max_theta = (r==0?inner_theta+dtheta[r]:pixel[r-1][w].max_theta+dtheta[r]);
	pixel[r][w].max_costheta = (r==0?inner_costheta:pixel[r-1][w].min_costheta);
	pixel[r][w].min_costheta = (r==0?inner_costheta-dcostheta[r]:pixel[r-1][w].min_costheta-dcostheta[r]);
      }
      else {
	pixel[r][w].max_theta = (r==0?inner_theta:pixel[r-1][w].min_theta);
	pixel[r][w].min_theta = (r==0?inner_theta-dtheta[r]:pixel[r-1][w].min_theta-dtheta[r]);
	pixel[r][w].min_costheta = (r==0?inner_costheta:pixel[r-1][w].max_costheta);
	pixel[r][w].max_costheta = (r==0?inner_costheta+dcostheta[r]:pixel[r-1][w].max_costheta+dcostheta[r]);
      }
      pixel[r][w].min_phi = w*dphi + phi_rot;
      if(pixel[r][w].min_phi<0)      pixel[r][w].min_phi += 2*M_PI;
      if(pixel[r][w].min_phi>2*M_PI) pixel[r][w].min_phi -= 2*M_PI;
      pixel[r][w].max_phi = (w+1)*dphi + phi_rot;
      if(pixel[r][w].max_phi<0)      pixel[r][w].max_phi += 2*M_PI;
      if(pixel[r][w].max_phi>2*M_PI) pixel[r][w].max_phi -= 2*M_PI;
      if(pixel[r][w].max_phi<pixel[r][w].min_phi) pixel[r][w].max_phi += 2*M_PI;
    }
  }

  rand = new TRandom3();
  rand->SetSeed();

}

AnnularDetector::~AnnularDetector() {
  for (int i=0; i<num_rings; i++) delete [] pixel[i];
  delete [] pixel;
  delete [] ring_radius;
  delete [] dtheta;
  delete [] dcostheta;
  delete rand;
}

std::pair<double,double> AnnularDetector::GetThetaPhi(int iringch, int iwedgech) {

  //double theta = ((ValidRing(iringch) && ValidWedge(iwedgech)) ?
  //		  acos(rand->Uniform(cos(pixel[iringch][iwedgech].min_theta),cos(pixel[iringch][iwedgech].max_theta))) : 0);
  double phi = ((ValidRing(iringch) && ValidWedge(iwedgech)) ?
		rand->Uniform(pixel[iringch][iwedgech].min_phi,pixel[iringch][iwedgech].max_phi) : 0);
  double radius = ((ValidRing(iringch) && ValidWedge(iwedgech)) ?
		   rand->Uniform(pixel[iringch][iwedgech].min_radius,pixel[iringch][iwedgech].max_radius) : 0);

  double theta = atan2(radius,z_dist);
  
  while (phi > M_PI) phi -= 2*M_PI;

  return std::make_pair(theta,phi);

}

std::pair<double,double> AnnularDetector::GetCosThetaPhi(int iringch, int iwedgech) {

  //double costheta = ((ValidRing(iringch) && ValidWedge(iwedgech)) ?
  //		     rand->Uniform(pixel[iringch][iwedgech].min_costheta,pixel[iringch][iwedgech].max_costheta) : 0);
  double phi = ((ValidRing(iringch) && ValidWedge(iwedgech)) ?
		rand->Uniform(pixel[iringch][iwedgech].min_phi,pixel[iringch][iwedgech].max_phi) : 0);
  double radius = ((ValidRing(iringch) && ValidWedge(iwedgech)) ?
		   rand->Uniform(pixel[iringch][iwedgech].min_radius,pixel[iringch][iwedgech].max_radius) : 0);

  double theta = atan2(radius,z_dist);
  
  while (phi > M_PI) phi -= 2*M_PI;

  return std::make_pair(cos(theta),phi);

}

std::pair<int,int> AnnularDetector::GetRingWedge(double itheta, double iphi) {

  while (iphi < 0)      iphi += 2*M_PI;
  while (iphi > 2*M_PI) iphi -= 2*M_PI;

  int ring=-1, wedge=-1;
  for (int r=0; r<num_rings; r++) {
    for (int w=0; w<num_wedges; w++) {
      if (pixel[r][w].min_theta < itheta && pixel[r][w].max_theta > itheta &&
	  pixel[r][w].min_phi < iphi && pixel[r][w].max_phi > iphi)
	ring = r;
      if (pixel[r][w].min_theta < itheta && pixel[r][w].max_theta > itheta &&
	  pixel[r][w].min_phi < iphi && pixel[r][w].max_phi > iphi)
	wedge = w;
    }
  }
  return std::make_pair(ring,wedge);
}

std::pair<int,int> AnnularDetector::GetRingWedgeCos(double icostheta, double iphi) {

  while (iphi < 0)      iphi += 2*M_PI;
  while (iphi > 2*M_PI) iphi -= 2*M_PI;

  int ring=-1, wedge=-1;
  for (int r=0; r<num_rings; r++) {
    for (int w=0; w<num_wedges; w++) {
      if (pixel[r][w].min_costheta < icostheta && pixel[r][w].max_costheta > icostheta &&
	  pixel[r][w].min_phi < iphi && pixel[r][w].max_phi > iphi)
	ring = r;
      if (pixel[r][w].min_costheta < icostheta && pixel[r][w].max_costheta > icostheta &&
	  pixel[r][w].min_phi < iphi && pixel[r][w].max_phi > iphi)
	wedge = w;
    }
  }
  return std::make_pair(ring,wedge);
}

bool AnnularDetector::ValidRing(int r) {
  return ((r>=0&&r<num_rings) ? true : false);
}

bool AnnularDetector::ValidWedge(int w) {
  return ((w>=0&&w<num_wedges) ? true : false);
}
