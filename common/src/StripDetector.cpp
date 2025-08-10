#include "StripDetector.h"

#include <iostream>

using namespace std;

StripDetector::StripDetector(int ns, int nb, double len, double wid, double cphi, double cz, double crho) : STRIP_UNC_DEF(0.05) {

  num_strips = ns;
  num_backs = nb;
  
  length = fabs(len);
  back_length = length/num_backs;
  total_width = fabs(wid);
  strip_width = total_width/num_strips;

  while (cphi < 0) cphi += 2*M_PI;
  center_phi = cphi;
  center_z   = cz;
  center_rho = fabs(crho);

  strip = new channeldata[num_strips];
  for (int s=0; s<num_strips; s++) {
    strip[s].min_z = center_z - length/2;
    strip[s].max_z = center_z + length/2;
    //to get x and y, first assume we're along x s.t. center @ y=0 and x=center_rho
    strip[s].min_y = (s==0 ? -total_width/2 : strip[s-1].max_y);
    strip[s].max_y = -total_width/2 + (s+1)*strip_width;
    strip[s].min_x = center_rho;
    strip[s].max_x = center_rho;
    //when we actually call the get functions, we'll rotate by center_phi via a matrix operation
    //AFTER we've selected a uniformly random point
  }

  pixel = new channeldata*[num_backs];
  for (int b=0; b<num_backs; b++) {
    pixel[b] = new channeldata[num_strips];
    for (int s=0; s<num_strips; s++) {
      pixel[b][s].min_z = -length/2 + b*back_length + center_z;
      pixel[b][s].max_z = -length/2 + (b+1)*back_length + center_z;
      pixel[b][s].min_y = (s==0 ? -total_width/2 : strip[s-1].max_y);
      pixel[b][s].max_y = -total_width/2 + (s+1)*strip_width;
      pixel[b][s].min_x = center_rho;
      pixel[b][s].max_x = center_rho;
    }
  }
  
  rand = new TRandom3();
  rand->SetSeed();

}

StripDetector::~StripDetector() {
  for (int b=0; b<num_backs; b++)
    delete [] pixel[b];
  delete [] pixel;
  delete [] strip;
  delete rand;
}

std::pair<double,double> StripDetector::GetThetaPhi(int istripch, double istrip_ratio, double istrip_unc) {

  if (!ValidChannel(istripch) || !ValidRatio(istrip_ratio)) return std::make_pair(0,0);

  double strip_unc = (istrip_unc>=0?istrip_unc:STRIP_UNC_DEF);
  
  //recall we're still assuming phi=0 det
  //***MUST ALSO MAKE SURE we're not randomizing a "pixel" length beyond the physical edge
  //double zhit = rand->Uniform(max(istrip_ratio*(length/2)+center_z-strip_unc*length,-length/2+center_z),
  //                            min(istrip_ratio*(length/2)+center_z+strip_unc*length,+length/2+center_z));
  //alternate: randomize as gaussian (bc energy-based)
  double zhit = rand->Gaus(istrip_ratio*(length/2)+center_z,2*strip_unc*length/2.355);
  double yhit = rand->Uniform(strip[istripch].min_y,strip[istripch].max_y);
  double xhit = center_rho;

  //NOW rotate by appropriate phi
  double xhit_rot =  xhit*cos(center_phi) - yhit*sin(center_phi);
  double yhit_rot =  xhit*sin(center_phi) + yhit*cos(center_phi);

  xhit = xhit_rot;
  yhit = yhit_rot;

  double rhit = sqrt(xhit*xhit + yhit*yhit + zhit*zhit);

  double theta_hit = acos(zhit/rhit);
  double phi_hit = atan2(yhit,xhit);

  return std::make_pair(theta_hit,phi_hit);

}

std::pair<double,double> StripDetector::GetCosThetaPhi(int istripch, double istrip_ratio, double istrip_unc) {

  if (!ValidChannel(istripch) || !ValidRatio(istrip_ratio)) return std::make_pair(0,0);

  double strip_unc = (istrip_unc>=0?istrip_unc:STRIP_UNC_DEF);
  
  //recall we're still assuming phi=0 det:
  //double zhit = istrip_ratio*(length/2) + center_z + rand->Uniform(-STRIP_UNC*length,STRIP_UNC*length);
  //double zhit = rand->Uniform(max(istrip_ratio*(length/2)+center_z-strip_unc*length,-length/2+center_z),
  //			        min(istrip_ratio*(length/2)+center_z+strip_unc*length,+length/2+center_z));
  double zhit = rand->Gaus(istrip_ratio*(length/2)+center_z,2*strip_unc*length/2.355);
  double yhit = rand->Uniform(strip[istripch].min_y,strip[istripch].max_y);
  double xhit = center_rho;

  //NOW rotate by appropriate phi
  double xhit_rot =  xhit*cos(center_phi) - yhit*sin(center_phi);
  double yhit_rot =  xhit*sin(center_phi) + yhit*cos(center_phi);

  xhit = xhit_rot;
  yhit = yhit_rot;

  double rhit = sqrt(xhit*xhit + yhit*yhit + zhit*zhit);

  double costheta_hit = zhit/rhit;
  double phi_hit = atan2(yhit,xhit);

  return std::make_pair(costheta_hit,phi_hit);

}

std::pair<int,double> StripDetector::GetChannelRatio(double itheta, double iphi) {

  while (iphi < 0) iphi += 2*M_PI;

  //to make the math easier (and the same for each det), rotate the input phi
  //BACKWARD by the phi of the det, s.t. we are looking along x-axis
  iphi -= center_phi;

  if (iphi > M_PI) iphi -= 2*M_PI;

  //then we can check easily whether it even hit the detector in phi
  double det_max_phi = atan2(total_width/2,center_rho);
  double det_min_phi = -det_max_phi;
  
  if (iphi < det_min_phi || iphi > det_max_phi) return std::make_pair(-1,0);

  //for theta it's not so simple, so we have to go through the typical plane-intersect method
  
  //first thing's first: we have a fixed x for the entire detector plane:
  double xhit = center_rho;
  //thus we find the corresponding y and z for that fixed x, given the input theta and phi:
  double yhit = xhit*tan(iphi);
  double zhit = sqrt(xhit*xhit+yhit*yhit)/tan(itheta);

  int chhit=-1;
  double ratio=0;
  for (int s=0; s<num_strips; s++) {
    if (xhit >= strip[s].min_x && xhit <= strip[s].max_x &&
	yhit >= strip[s].min_y && yhit <= strip[s].max_y &&
	zhit >= strip[s].min_z && zhit <= strip[s].max_z) {
      chhit = s;
      ratio = (zhit-center_z)/(length/2);
    }
  }

  return std::make_pair(chhit,ratio);

}

std::pair<int,double> StripDetector::GetChannelRatioCos(double icostheta, double iphi) {

  while (iphi < 0) iphi += 2*M_PI;

  //to make the math easier (and the same for each det), rotate the input phi
  //BACKWARD by the phi of the det, s.t. we are looking along x-axis
  iphi -= center_phi;

  if (iphi > M_PI) iphi -= 2*M_PI;

  //then we can check easily whether it even hit the detector in phi
  double det_max_phi = atan2(total_width/2,center_rho);
  double det_min_phi = -det_max_phi;
  
  if (iphi < det_min_phi || iphi > det_max_phi) return std::make_pair(-1,0);

  //for theta it's not so simple, so we have to go through the typical plane-intersect method
  
  //first thing's first: we have a fixed x for the entire detector plane:
  double xhit = center_rho;
  //thus we find the corresponding y and z for that fixed x, given the input theta and phi:
  double yhit = xhit*tan(iphi);
  double zhit = sqrt(xhit*xhit+yhit*yhit)/tan(acos(icostheta));

  int chhit=-1;
  double ratio=0;
  for (int s=0; s<num_strips; s++) {
    if (xhit >= strip[s].min_x && xhit <= strip[s].max_x &&
	yhit >= strip[s].min_y && yhit <= strip[s].max_y &&
	zhit >= strip[s].min_z && zhit <= strip[s].max_z) {
      chhit = s;
      ratio = (zhit-center_z)/(length/2);
    }
  }

  return std::make_pair(chhit,ratio);

}

std::pair<double,double> StripDetector::GetThetaPhiBack(int istripch, int ibackch) {

  if (!ValidChannel(istripch) || !ValidBack(ibackch)) return std::make_pair(0,0);

  //recall we're still assuming phi=0 det:
  double zhit = rand->Uniform(pixel[ibackch][istripch].min_z,pixel[ibackch][istripch].max_z);
  double yhit = rand->Uniform(strip[istripch].min_y,strip[istripch].max_y);
  double xhit = center_rho;

  //NOW rotate by appropriate phi
  double xhit_rot =  xhit*cos(center_phi) - yhit*sin(center_phi);
  double yhit_rot =  xhit*sin(center_phi) + yhit*cos(center_phi);

  xhit = xhit_rot;
  yhit = yhit_rot;

  double rhit = sqrt(xhit*xhit + yhit*yhit + zhit*zhit);

  double theta_hit = acos(zhit/rhit);
  double phi_hit = atan2(yhit,xhit);

  return std::make_pair(theta_hit,phi_hit);

}

std::pair<double,double> StripDetector::GetCosThetaPhiBack(int istripch, int ibackch) {

  if (!ValidChannel(istripch) || !ValidBack(ibackch)) return std::make_pair(0,0);

  //recall we're still assuming phi=0 det:
  double zhit = rand->Uniform(pixel[ibackch][istripch].min_z,pixel[ibackch][istripch].max_z);
  double yhit = rand->Uniform(strip[istripch].min_y,strip[istripch].max_y);
  double xhit = center_rho;

  //NOW rotate by appropriate phi
  double xhit_rot =  xhit*cos(center_phi) - yhit*sin(center_phi);
  double yhit_rot =  xhit*sin(center_phi) + yhit*cos(center_phi);

  xhit = xhit_rot;
  yhit = yhit_rot;

  double rhit = sqrt(xhit*xhit + yhit*yhit + zhit*zhit);

  double costheta_hit = zhit/rhit;
  double phi_hit = atan2(yhit,xhit);

  return std::make_pair(costheta_hit,phi_hit);

}

std::pair<int,double> StripDetector::GetStripBack(double itheta, double iphi) {

  while (iphi < 0) iphi += 2*M_PI;

  //to make the math easier (and the same for each det), rotate the input phi
  //BACKWARD by the phi of the det, s.t. we are looking along x-axis
  iphi -= center_phi;

  if (iphi > M_PI) iphi -= 2*M_PI;

  //then we can check easily whether it even hit the detector in phi
  double det_max_phi = atan2(total_width/2,center_rho);
  double det_min_phi = -det_max_phi;
  
  if (iphi < det_min_phi || iphi > det_max_phi) return std::make_pair(-1,0);

  //for theta it's not so simple, so we have to go through the typical plane-intersect method
  
  //first thing's first: we have a fixed x for the entire detector plane:
  double xhit = center_rho;
  //thus we find the corresponding y and z for that fixed x, given the input theta and phi:
  double yhit = xhit*tan(iphi);
  double zhit = sqrt(xhit*xhit+yhit*yhit)/tan(itheta);

  int chhit=-1, bhit=-1;
  double ratio=0;
  for (int b=0; b<num_backs; b++) {
    for (int s=0; s<num_strips; s++) {
      if (xhit >= pixel[b][s].min_x && xhit <= pixel[b][s].max_x &&
	  yhit >= pixel[b][s].min_y && yhit <= pixel[b][s].max_y &&
	  zhit >= pixel[b][s].min_z && zhit <= pixel[b][s].max_z) {
	bhit  = b;
	chhit = s;
      }
    }
  }

  return std::make_pair(chhit,bhit);

}

std::pair<int,double> StripDetector::GetStripBackCos(double icostheta, double iphi) {

  while (iphi < 0) iphi += 2*M_PI;

  //to make the math easier (and the same for each det), rotate the input phi
  //BACKWARD by the phi of the det, s.t. we are looking along x-axis
  iphi -= center_phi;

  if (iphi > M_PI) iphi -= 2*M_PI;

  //then we can check easily whether it even hit the detector in phi
  double det_max_phi = atan2(total_width/2,center_rho);
  double det_min_phi = -det_max_phi;
  
  if (iphi < det_min_phi || iphi > det_max_phi) return std::make_pair(-1,0);

  //for theta it's not so simple, so we have to go through the typical plane-intersect method
  
  //first thing's first: we have a fixed x for the entire detector plane:
  double xhit = center_rho;
  //thus we find the corresponding y and z for that fixed x, given the input theta and phi:
  double yhit = xhit*tan(iphi);
  double zhit = sqrt(xhit*xhit+yhit*yhit)/tan(acos(icostheta));

  int chhit=-1, bhit=-1;
  double ratio=0;
  for (int b=0; b<num_backs; b++) {
    for (int s=0; s<num_strips; s++) {
      if (xhit >= pixel[b][s].min_x && xhit <= pixel[b][s].max_x &&
	  yhit >= pixel[b][s].min_y && yhit <= pixel[b][s].max_y &&
	  zhit >= pixel[b][s].min_z && zhit <= pixel[b][s].max_z) {
	bhit  = b;
	chhit = s;
      }
    }
  }

  return std::make_pair(chhit,bhit);

}

bool StripDetector::ValidChannel(int f) {
  return ((f >= 0 && f < num_strips) ? true : false);
}

bool StripDetector::ValidBack(int b) {
  return ((b >= 0 && b < num_backs) ? true : false);
}

bool StripDetector::ValidRatio(double r) {
  return ((r >= -1 && r <= 1) ? true : false);
}
