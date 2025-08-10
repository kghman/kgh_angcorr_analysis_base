#include <iostream>
#include <fstream>
#include <cstring>
#include <cmath>
#include <cstdlib>

#include <TFile.h>
#include <TH2.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TRandom3.h>

#include "AnnularDetector.h"

using namespace std;

int main(int argc, char** argv) {

  const int NUM_EVENTS = 1E7;

  TRandom3 *rand = new TRandom3();
  rand->SetSeed();

  //S1
  const int NUM_RINGS       = 16;
  const int NUM_WEDGES      = 16;
  const double INNER_RADIUS = 2.4; //cm
  const double OUTER_RADIUS = 4.8; //cm
  const double Z_POSITION   = 2.7; //cm
  const double PHI_ROTATION =-4.0*(2*M_PI)/NUM_WEDGES; //rad
  const double CUTOFF_FRAC  = 0.8125;

  /*
  //S2
  const int NUM_RINGS       = 16;
  const int NUM_WEDGES      = 16;
  const double INNER_RADIUS = 1.15; //cm
  const double OUTER_RADIUS = 3.5; //cm
  const double Z_POSITION   =-8.73; //cm
  const double PHI_ROTATION =-4.0*(2*M_PI)/NUM_WEDGES; //rad
  const double CUTOFF_FRAC  = 0.6875;
  */

  AnnularDetector det(NUM_RINGS,NUM_WEDGES,INNER_RADIUS,OUTER_RADIUS,Z_POSITION,PHI_ROTATION);

  TFile outfile("trd_outfile.root","RECREATE");

  TH2F rand_geom("rand_geom","rand_geom;#phi (rad);cos(#theta)",500,-M_PI,M_PI,500,-1,1);
  rand_geom.SetStats(0);

  TH2F det_geom("det_geom","S1 with cut top/bottom;#phi (rad);cos(#theta)",500,-M_PI,M_PI,500,-1,1);
  det_geom.SetStats(0);

  TH2F det_xy_proj("det_xy_proj","det_xy_proj;x (cm);y (cm)",500,-5,5,500,-5,5);
  det_xy_proj.SetStats(0);
  
  int hit=0, hitindet=0;
  
  for (int i=0; i<NUM_EVENTS; i++) {

    double costheta_rand = rand->Uniform(-1,1);
    double phi_rand = rand->Uniform(-M_PI,M_PI);

    rand_geom.Fill(phi_rand,costheta_rand);

    std::pair<int,int> ring_wedge_hit = det.GetRingWedgeCos(costheta_rand,phi_rand);

    if (ring_wedge_hit.first!=-1 && ring_wedge_hit.second!=-1) {

      hit++;
      
      std::pair<double,double> costheta_phi_hit = det.GetCosThetaPhi(ring_wedge_hit.first,ring_wedge_hit.second);

      double x = sqrt(pow(Z_POSITION,2) + pow(det.Get_Ring_Radius(ring_wedge_hit.first),2))*cos(costheta_phi_hit.second)*sin(acos(costheta_phi_hit.first));
      double y = sqrt(pow(Z_POSITION,2) + pow(det.Get_Ring_Radius(ring_wedge_hit.first),2))*sin(costheta_phi_hit.second)*sin(acos(costheta_phi_hit.first));
      
      if (abs(y) <= OUTER_RADIUS*CUTOFF_FRAC) {
	det_geom.Fill(costheta_phi_hit.second,costheta_phi_hit.first);
	det_xy_proj.Fill(x,y);
	hitindet++;
      }
      
    }

  }

  cout << "EVENTS HIT: " << 100.*hitindet/NUM_EVENTS << "% (" << 100.*hitindet/hit << "% of det)" << endl;
  
  outfile.cd();
  rand_geom.Write();
  det_geom.Write();
  det_xy_proj.Write();

  delete rand;

  return 0;

}
