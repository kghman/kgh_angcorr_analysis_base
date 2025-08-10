#include <iostream>
#include <fstream>
#include <cmath>

#include <TROOT.h>
#include <TRandom3.h>
#include <TH2.h>
#include <TFile.h>
#include <TLatex.h>
#include <TCanvas.h>

using namespace std;

void Euler_Forward_Rot(double *outvec, double *angs, double *invec);
void Euler_Backward_Rot(double *outvec, double *angs, double *invec);

int main (int argc, char** argv) {

  const int NUM_EVENTS = 1E6;

  double euler_vec[3];
  euler_vec[0] = -M_PI/2;
  euler_vec[1] = -M_PI/2;
  euler_vec[2] =  0;

  TRandom3 *rand = new TRandom3();
  rand->SetSeed();

  TH2F randcos_cosvphi_baxisq_hist("randcos_cosvphi_baxisq_hist","BAxisQ (Random cos(#theta)));#phi (rad.);cos(#theta)",200,-M_PI,M_PI,200,-1,1);
  TH2F randcos_cosvpsi_perpq_hist("randcos_cosvpsi_perpq_hist","PerpQ (Random cos(#theta));#psi (rad.);cos(#chi)",200,-M_PI,M_PI,200,-1,1);
  TH2F randcos_thetavphi_baxisq_hist("randcos_thetavphi_baxisq_hist","BAxisQ (Random cos(#theta)));#phi (rad.);#theta (rad.)",200,-M_PI,M_PI,200,0,M_PI);
  TH2F randcos_chivpsi_perpq_hist("randcos_chivphi_perpq_hist","PerpQ (Random cos(#theta));#psi (rad.);#chi (rad.)",200,-M_PI,M_PI,200,0,M_PI);

  TH2F randang_cosvphi_baxisq_hist("randang_cosvphi_baxisq_hist","BAxisQ (Random #theta));#phi (rad.);cos(#theta)",200,-M_PI,M_PI,200,-1,1);
  TH2F randang_cosvpsi_perpq_hist("randang_cosvpsi_perpq_hist","PerpQ (Random #theta);#psi (rad.);cos(#chi)",200,-M_PI,M_PI,200,-1,1);
  TH2F randang_thetavphi_baxisq_hist("randang_thetavphi_baxisq_hist","BAxisQ (Random #theta));#phi (rad.);#theta (rad.)",200,-M_PI,M_PI,200,0,M_PI);
  TH2F randang_chivpsi_perpq_hist("randang_chivphi_perpq_hist","PerpQ (Random #theta);#psi (rad.);#chi (rad.)",200,-M_PI,M_PI,200,0,M_PI);

  for (int i=0; i<NUM_EVENTS; i++) {
    
    double costheta, theta, phi;
    double coschi, chi, psi;

    double baxisq_vec[3];
    double perpq_vec[3];

    //first randomize in cos(theta)

    costheta = rand->Uniform(-1,1);
    theta = acos(costheta);
    phi = rand->Uniform(-M_PI,M_PI);

    randcos_cosvphi_baxisq_hist.Fill(phi,costheta);
    randcos_thetavphi_baxisq_hist.Fill(phi,theta);

    /*
    baxisq_vec[0] = cos(phi)*sin(theta);
    baxisq_vec[1] = sin(phi)*sin(theta);
    baxisq_vec[2] = costheta;
    Euler_Forward_Rot(perpq_vec,euler_vec,baxisq_vec);   

    coschi = perpq_vec[2];
    chi = acos(coschi);
    psi = atan2(perpq_vec[1],perpq_vec[0]);
    */

    coschi = sin(phi)*sin(theta);
    chi = acos(coschi);
    psi = atan2(-cos(phi)*sin(theta),cos(theta));

    randcos_cosvpsi_perpq_hist.Fill(psi,coschi);
    randcos_chivpsi_perpq_hist.Fill(psi,chi);

    //then in theta

    theta = rand->Uniform(0,M_PI);
    costheta = cos(theta);
    phi = rand->Uniform(-M_PI,M_PI);

    randang_cosvphi_baxisq_hist.Fill(phi,costheta);
    randang_thetavphi_baxisq_hist.Fill(phi,theta);

    /*
    baxisq_vec[0] = cos(phi)*sin(theta);
    baxisq_vec[1] = sin(phi)*sin(theta);
    baxisq_vec[2] = costheta;
    Euler_Forward_Rot(perpq_vec,euler_vec,baxisq_vec);

    coschi = perpq_vec[2];
    chi = acos(coschi);
    psi = atan2(perpq_vec[1],perpq_vec[0]);
    */

    coschi = sin(phi)*sin(theta);
    chi = acos(coschi);
    psi = atan2(-cos(phi)*sin(theta),cos(theta));

    randang_cosvpsi_perpq_hist.Fill(psi,coschi);
    randang_chivpsi_perpq_hist.Fill(psi,chi);

  }

  randcos_cosvphi_baxisq_hist.SetStats(0);
  randcos_cosvpsi_perpq_hist.SetStats(0);
  randcos_thetavphi_baxisq_hist.SetStats(0);
  randcos_chivpsi_perpq_hist.SetStats(0);

  randang_cosvphi_baxisq_hist.SetStats(0);
  randang_cosvpsi_perpq_hist.SetStats(0);
  randang_thetavphi_baxisq_hist.SetStats(0);
  randang_chivpsi_perpq_hist.SetStats(0);

  TCanvas randcos_canv("randcos_canv","randcos_canv");
  randcos_canv.Divide(2,2);
  randcos_canv.cd(1);
  randcos_cosvphi_baxisq_hist.Draw("colz");
  randcos_canv.cd(2);
  randcos_cosvpsi_perpq_hist.Draw("colz");
  randcos_canv.cd(3);
  randcos_thetavphi_baxisq_hist.Draw("colz");
  randcos_canv.cd(4);
  randcos_chivpsi_perpq_hist.Draw("colz");

  TCanvas randang_canv("randang_canv","randang_canv");
  randang_canv.Divide(2,2);
  randang_canv.cd(1);
  randang_cosvphi_baxisq_hist.Draw("colz");
  randang_canv.cd(2);
  randang_cosvpsi_perpq_hist.Draw("colz");
  randang_canv.cd(3);
  randang_thetavphi_baxisq_hist.Draw("colz");
  randang_canv.cd(4);
  randang_chivpsi_perpq_hist.Draw("colz");

  TFile outfile("testing_solid_angle_outfile.root","RECREATE");
  randcos_canv.Write();
  randang_canv.Write();
  outfile.Close();

  delete rand;

  return 0;

}

void Euler_Forward_Rot(double *outvec, double *angs, double *invec) {

  outvec[0] = (cos(angs[0])*cos(angs[1])*cos(angs[2])-sin(angs[0])*sin(angs[2]))*invec[0]
             +(sin(angs[0])*cos(angs[1])*cos(angs[2])+cos(angs[0])*sin(angs[2]))*invec[1]
                                                    -(sin(angs[1])*cos(angs[2]))*invec[2];

  outvec[1] =-(cos(angs[0])*cos(angs[1])*sin(angs[2])+sin(angs[0])*cos(angs[2]))*invec[0]
             -(sin(angs[0])*cos(angs[1])*sin(angs[2])-cos(angs[0])*cos(angs[2]))*invec[1]
                                                    +(sin(angs[1])*sin(angs[2]))*invec[2];

  outvec[2] = cos(angs[0])*sin(angs[1])*invec[0]
             +sin(angs[0])*sin(angs[1])*invec[1]
                          +cos(angs[1])*invec[2];

}

void Euler_Backward_Rot(double *outvec, double *angs, double *invec) {

  outvec[0] = (cos(angs[0])*cos(angs[1])*cos(angs[2])-sin(angs[0])*sin(angs[2]))*invec[0]
             -(cos(angs[0])*cos(angs[1])*sin(angs[2])+sin(angs[0])*cos(angs[2]))*invec[1]
                                                    +(cos(angs[0])*sin(angs[1]))*invec[2];

  outvec[1] = (sin(angs[0])*cos(angs[1])*cos(angs[2])+cos(angs[0])*sin(angs[2]))*invec[0]
             -(sin(angs[0])*cos(angs[1])*sin(angs[2])-cos(angs[0])*cos(angs[2]))*invec[1]
                                                    +(sin(angs[0])*sin(angs[1]))*invec[2];

  outvec[2] =-sin(angs[1])*cos(angs[2])*invec[0]
             +sin(angs[1])*sin(angs[2])*invec[1]
                          +cos(angs[1])*invec[2];

}

