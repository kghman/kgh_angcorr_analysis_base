#include <TGraph.h>
#include <TMultiGraph.h>
#include <TFile.h>
#include <TTree.h>
#include <TCanvas.h>
#include <TMath.h>
#include <TH2.h>
#include <TLatex.h>

#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <cstring>

extern "C" {
#include "wiko_wiko5.h"
}

double Afunc(double Io,double If,int lambda,double delta);

using namespace std;

int main() {

  const int NUMPOINTS = 10000;
  const double I_O = 2, I_F = 2;

  double delta[NUMPOINTS];
  double A2[2][NUMPOINTS]; //one for delta > 0, one for delta < 0
  double A4[2][NUMPOINTS];

  for (int i=0; i<NUMPOINTS; i++) {
    delta[i] = (1.*i+1)/100;
    //positive delta:
    A2[0][i] = Afunc(I_O,I_F,2,delta[i]);
    A4[0][i] = Afunc(I_O,I_F,4,delta[i]);
    //negative delta:
    A2[1][i] = Afunc(I_O,I_F,2,-delta[i]);
    A4[1][i] = Afunc(I_O,I_F,4,-delta[i]);
  }

  TFile outfile("outfile_Adeltas.root","RECREATE");

  TGraph A2_gr_pos(NUMPOINTS,delta,A2[0]);
  A2_gr_pos.SetTitle("A_{2}");
  A2_gr_pos.SetName("A2_gr_pos");
  A2_gr_pos.SetLineColor(2);
  TGraph A4_gr_pos(NUMPOINTS,delta,A4[0]);
  A4_gr_pos.SetTitle("A_{4}");
  A4_gr_pos.SetName("A4_gr_pos");
  A4_gr_pos.SetLineColor(4);
  TMultiGraph pos_mg;
  pos_mg.SetTitle(Form("A_{#lambda} for Transition %f --> %f, #delta > 0;|#delta|;A_{#lambda}",I_O,I_F));
  pos_mg.SetName("pos_mg");
  pos_mg.Add(&A2_gr_pos);
  pos_mg.Add(&A4_gr_pos);

  TGraph A2_gr_neg(NUMPOINTS,delta,A2[1]);
  A2_gr_neg.SetTitle("A_{2}");
  A2_gr_neg.SetName("A2_gr_neg");
  A2_gr_neg.SetLineColor(2);
  TGraph A4_gr_neg(NUMPOINTS,delta,A4[1]);
  A4_gr_neg.SetTitle("A_{4}");
  A4_gr_neg.SetName("A4_gr_neg");
  A4_gr_neg.SetLineColor(4);
  TMultiGraph neg_mg;
  neg_mg.SetTitle(Form("A_{#lambda} for Transition %f --> %f, #delta < 0;|#delta|;A_{#lambda}",I_O,I_F));
  neg_mg.SetName("neg_mg");
  neg_mg.Add(&A2_gr_neg);
  neg_mg.Add(&A4_gr_neg);

  pos_mg.Write();
  neg_mg.Write();

  outfile.Close();

  return 0;

}

double Afunc(double Io,double If,int lambda,double delta) {

  int l = (int)TMath::Abs(Io - If)+1;
  int lp = l + 1;

  double A = wiko_f3(l,l,Io,If,lambda,lambda,0);
  A += 2*delta*wiko_f3(l,lp,Io,If,lambda,lambda,0);
  A += delta*delta*wiko_f3(lp,lp,Io,If,lambda,lambda,0);
  A /= (1+delta*delta);

  return A;

}
