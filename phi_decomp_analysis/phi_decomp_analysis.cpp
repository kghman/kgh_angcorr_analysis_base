/*

  phi_decomp_analysis.cpp:

    Analysis intended to take the output rootfiles from either bees or geebees
    and perform a transform on their 2D decay distributions from continuous phi
    to discrete q tensor indices.

    Ken H -- Mar 2022

*/

#include <iostream>
#include <fstream>
#include <cstdlib>
#include <cstring>
#include <chrono>
#include <complex>
#include <vector>

#include <TGraph.h>
#include <TFile.h>
#include <TTree.h>
#include <TCanvas.h>
#include <TMath.h>
#include <TH1.h>
#include <TH2.h>
#include <TH3.h>
#include <TLatex.h>
#include <TPaveLabel.h>
#include <TPad.h>

#include "constants.h"
#include "WikoParticleClass.h"
#include "BeesHolder.h"
#include "wiko_wiko5.hpp"
#include "FrameConverter.h"

using namespace std;

int main(int argc, char **argv) {

  if (argc!=2) {
    cout << "***WARNING: supply input file as argument." << endl;
    return 1;
  }

  ifstream input_file;
  input_file.open(argv[1]);
  if (!input_file) {
    cout << "***WARNING: input file " << argv[1] << " not found." << endl;
    return 2;
  }

  char label[100];
  char input_rootfile_name[100];
  input_file >> label >> input_rootfile_name;
  
  TFile input_rootfile(input_rootfile_name,"READ");
  if (!input_rootfile.IsOpen()) {
    cout << "***WARNING: error opening rootfile " << input_rootfile_name << endl;
    return 3;
  }

  char output_name[100];
  input_file >> label >> output_name;

  char output_rootfile_name[100];
  strcpy(output_rootfile_name,"rootfiles/");
  strcat(output_rootfile_name,output_name);
  strcat(output_rootfile_name,".root");
  
  input_file.close();
  
  double xy_angle_in, z_angle_in, W_in;
  TTree *intree = (TTree*)input_rootfile.Get("data_tree");
  intree->SetBranchAddress("xy_angle",&xy_angle_in);
  intree->SetBranchAddress("z_angle",&z_angle_in);
  intree->SetBranchAddress("W",&W_in);

  TFile output_rootfile(output_rootfile_name,"RECREATE");

  //relevant fixed quantities to pull from file
  int k_max;
  int NUM_XY_ANGS, NUM_Z_ANGS;
  TTree *intree_inputs = (TTree*)input_rootfile.Get("input_tree");
  intree_inputs->SetBranchAddress("k_max",&k_max);
  intree_inputs->SetBranchAddress("NUM_XY_ANGS",&NUM_XY_ANGS);
  intree_inputs->SetBranchAddress("NUM_Z_ANGS",&NUM_Z_ANGS);
  intree_inputs->GetEntry(0);

  double dxyang = 2.*M_PI/NUM_XY_ANGS;
  double dzang  = 2./NUM_Z_ANGS;
  
  double xy_angle_arr[NUM_XY_ANGS];
  double z_angle_arr[NUM_Z_ANGS];
  double W_arr[NUM_XY_ANGS][NUM_Z_ANGS];

  for (int i=0; i<NUM_XY_ANGS; i++)
    xy_angle_arr[i] = -M_PI + (i+0.5)*dxyang;
  for (int i=0; i<NUM_Z_ANGS; i++)
    z_angle_arr[i] = -1 + (i+0.5)*dzang;
  
  TH2F input_hist("input_hist","input_hist",NUM_XY_ANGS,-M_PI,M_PI,NUM_Z_ANGS,-1,1);

  for (int ent=0; ent<intree->GetEntries(); ent++) {

    cout << Form("\rGrabbing data from input tree (%.0f%%)...",100.*ent/intree->GetEntries()) << flush;
    
    intree->GetEntry(ent);

    input_hist.Fill(xy_angle_in,z_angle_in,W_in);

    for (int i=0; i<NUM_XY_ANGS; i++) {
      for (int j=0; j<NUM_Z_ANGS; j++) {
	if (xy_angle_in == xy_angle_arr[i] && z_angle_in == z_angle_arr[j]) {
	  W_arr[i][j] = W_in;
	  break;
	}
      }
    }
    
  }

  cout << "done." << endl
       << "Integrating decomps and writing to files...";
  
  double cos_decomp_arr[k_max+1][NUM_Z_ANGS], sin_decomp_arr[k_max+1][NUM_Z_ANGS];

  for (int i=0; i<NUM_Z_ANGS; i++) {
    for (int q=0; q<=k_max; q++) {      
      cos_decomp_arr[q][i] = 0;
      sin_decomp_arr[q][i] = 0;
      for (int j=0; j<NUM_XY_ANGS; j++) {
	cos_decomp_arr[q][i] += W_arr[j][i]*cos(q*xy_angle_arr[j])*dxyang;
	sin_decomp_arr[q][i] += W_arr[j][i]*sin(q*xy_angle_arr[j])*dxyang;
      }
    }
  }

  ofstream cos_decomp_outfile;
  char cos_decomp_outfile_name[100];
  strcpy(cos_decomp_outfile_name,"decomp_outfiles/");
  strcat(cos_decomp_outfile_name,output_name);
  strcat(cos_decomp_outfile_name,"_cosdecomp.txt");    
  cos_decomp_outfile.open(cos_decomp_outfile_name);  
  cos_decomp_outfile << "MAX_Q: " << k_max << " | TOT_ANGS: " << NUM_Z_ANGS << endl
		     << "Cos(th)   ";
  for (int q=0; q<=k_max; q++)
    cos_decomp_outfile << Form("cosq%i   ",q);
  cos_decomp_outfile << endl;
  for (int i=0; i<NUM_Z_ANGS; i++) {
    cos_decomp_outfile << z_angle_arr[i] << "   ";
    for (int q=0; q<=k_max; q++)
      cos_decomp_outfile << cos_decomp_arr[q][i] << "   ";
    cos_decomp_outfile << endl;
  }
  cos_decomp_outfile.close();

  ofstream sin_decomp_outfile;
  char sin_decomp_outfile_name[100];
  strcpy(sin_decomp_outfile_name,"decomp_outfiles/");
  strcat(sin_decomp_outfile_name,output_name);
  strcat(sin_decomp_outfile_name,"_sindecomp.txt");    
  sin_decomp_outfile.open(sin_decomp_outfile_name);  
  sin_decomp_outfile << "MAX_Q: " << k_max << " | TOT_ANGS: " << NUM_Z_ANGS << endl
		     << "Cos(th)   ";
  for (int q=0; q<=k_max; q++)
    sin_decomp_outfile << Form("sinq%i   ",q);
  sin_decomp_outfile << endl;
  for (int i=0; i<NUM_Z_ANGS; i++) {
    sin_decomp_outfile << z_angle_arr[i] << "   ";
    for (int q=0; q<=k_max; q++)
      sin_decomp_outfile << sin_decomp_arr[q][i] << "   ";
    sin_decomp_outfile << endl;
  }
  sin_decomp_outfile.close();
  
  output_rootfile.cd();
  
  input_hist.Write();
  
  input_rootfile.Close();
  output_rootfile.Close();

  cout << "done." << endl;
  
  return 0;

}
