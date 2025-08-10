#include <iostream>
#include <cmath>

#include <TROOT.h>
#include <TFile.h>
#include <TH1.h>
#include <TH2.h>
#include <TCanvas.h>
#include <TRandom3.h>

#include "constants.h"
#include "TargetLadder.h"

using namespace std;

int main(int argc, char** argv) {

  const int NUM_MC = 1E6;

  double TL_TW    = 15.72; //target width (mm)
  double TL_TH    = 3.00;  //thickness along beam direction (mm)
  double TL_TSW   = 4.08;  //target-side width (mm)
  double TL_FSW   = 2.85;  //far-side width (mm)
  double TL_VAR   = 160*DEGTORAD; //rotation around vertical (y) axis from beam direction (rad)
  double TL_MEDFT = 12.5;  //magnet edge dist from target (mm)
  double TL_MTH   = 5.00;  //magnet thickness (mm)
  double TL_MD    = 6.50;  //magnet depth (mm)
  double TL_MW    = 12.5;  //magnet width (mm)

  TargetLadder target_ladder_tilted_wmagnets(TL_TW,TL_TH,TL_TSW,TL_FSW,TL_VAR,TL_MEDFT,TL_MTH,TL_MD,TL_MW);
  TargetLadder target_ladder_untilted_wmagnets(TL_TW,TL_TH,TL_TSW,TL_FSW,M_PI,TL_MEDFT,TL_MTH,TL_MD,TL_MW);
  TargetLadder target_ladder_tilted_nomagnets(TL_TW,TL_TH,TL_TSW,TL_FSW,TL_VAR);
  TargetLadder target_ladder_untilted_nomagnets(TL_TW,TL_TH,TL_TSW,TL_FSW,M_PI);

  TH2F target_ladder_tilted_wmagnets_hist("target_ladder_tilted_wmagnets_hist","Tilted with Magnets;#phi (deg.);cos(#theta)",360,-180,180,180,-1,1);
  TH2F target_ladder_untilted_wmagnets_hist("target_ladder_untilted_wmagnets_hist","Untilted with Magnets;#phi (deg.);cos(#theta)",360,-180,180,180,-1,1);
  TH2F target_ladder_tilted_nomagnets_hist("target_ladder_tilted_nomagnets_hist","Tilted NO Magnets;#phi (deg.);cos(#theta)",360,-180,180,180,-1,1);
  TH2F target_ladder_untilted_nomagnets_hist("target_ladder_untilted_nomagnets_hist","Untilted NO Magnets;#phi (deg.);cos(#theta)",360,-180,180,180,-1,1);

  target_ladder_tilted_wmagnets_hist.SetStats(0);
  target_ladder_untilted_wmagnets_hist.SetStats(0);
  target_ladder_tilted_nomagnets_hist.SetStats(0);
  target_ladder_untilted_nomagnets_hist.SetStats(0);

  double costheta_mc=0;
  double phi_mc=0;

  TRandom3 *rand = new TRandom3();
  rand->SetSeed();

  for (int i=0; i<NUM_MC; i++) {

    cout << "\rWorking..." << (100.*i)/NUM_MC << "%";

    costheta_mc = rand->Uniform(-1,1);
    phi_mc      = rand->Uniform(-M_PI,M_PI);

    if (target_ladder_tilted_wmagnets.DoesEscape(acos(costheta_mc),phi_mc))    target_ladder_tilted_wmagnets_hist.Fill(phi_mc*RADTODEG,costheta_mc);
    if (target_ladder_untilted_wmagnets.DoesEscape(acos(costheta_mc),phi_mc))  target_ladder_untilted_wmagnets_hist.Fill(phi_mc*RADTODEG,costheta_mc);
    if (target_ladder_tilted_nomagnets.DoesEscape(acos(costheta_mc),phi_mc))   target_ladder_tilted_nomagnets_hist.Fill(phi_mc*RADTODEG,costheta_mc);
    if (target_ladder_untilted_nomagnets.DoesEscape(acos(costheta_mc),phi_mc)) target_ladder_untilted_nomagnets_hist.Fill(phi_mc*RADTODEG,costheta_mc);

  }

  cout << "done." << endl;

  TCanvas all_canv("all_canv","all_canv");
  all_canv.Divide(2,2);
  all_canv.cd(1);
  target_ladder_tilted_wmagnets_hist.Draw("colz");
  all_canv.cd(2);
  target_ladder_untilted_wmagnets_hist.Draw("colz");
  all_canv.cd(3);
  target_ladder_tilted_nomagnets_hist.Draw("colz");
  all_canv.cd(4);
  target_ladder_untilted_nomagnets_hist.Draw("colz");

  TFile outfile("tls_outfile.root","RECREATE");
  all_canv.Write();
  outfile.Close();
  
  delete rand;

  return 0;

}
