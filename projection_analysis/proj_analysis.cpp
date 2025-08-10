/*

  proj_analysis.cpp:

    Analysis intended to take the output rootfiles from either bees or geebees
    and separate out specific 1D projections of the 2D decay distribution for
    comparison to similar data.

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

  int NUM_XY_PROJS;
  input_file >> label >> NUM_XY_PROJS
	     >> label;
  double xy_proj_arr[NUM_XY_PROJS];
  for (int i=0; i<NUM_XY_PROJS; i++)
    input_file >> xy_proj_arr[i];

  int NUM_Z_PROJS;
  input_file >> label >> NUM_Z_PROJS
	     >> label;
  double z_proj_arr[NUM_Z_PROJS];
  for (int i=0; i<NUM_Z_PROJS; i++) {
    input_file >> z_proj_arr[i];
    z_proj_arr[i]*=DEGTORAD;
  }
  
  char output_name[100];
  input_file >> label >> output_name;

  char output_rootfile_name[100];
  strcpy(output_rootfile_name,"rootfiles/");
  strcat(output_rootfile_name,output_name);
  strcat(output_rootfile_name,".root");

  char proj_outfile_name[100];
  strcpy(proj_outfile_name,"proj_outfiles/");
  strcat(proj_outfile_name,output_name);
  strcat(proj_outfile_name,".txt");
  
  input_file.close();

  TFile output_rootfile(output_rootfile_name,"RECREATE");

  cout << "Pulling data from file " << input_rootfile_name << "...";
  
  //relevant fixed quantities to pull from file
  int k_max;
  int NUM_XY_ANGS, NUM_Z_ANGS;
  TTree *intree_inputs = (TTree*)input_rootfile.Get("input_tree");
  intree_inputs->SetBranchAddress("k_max",&k_max);
  intree_inputs->SetBranchAddress("NUM_XY_ANGS",&NUM_XY_ANGS);
  intree_inputs->SetBranchAddress("NUM_Z_ANGS",&NUM_Z_ANGS);
  intree_inputs->GetEntry(0);

  int k_in, q_in;
  double akq_real_in, akq_imag_in;
  TTree *intree = (TTree*)input_rootfile.Get("coeff_tree");
  intree->SetBranchAddress("k_fill",&k_in);
  intree->SetBranchAddress("q_fill",&q_in);
  intree->SetBranchAddress("akq_real_fill",&akq_real_in);
  intree->SetBranchAddress("akq_imag_fill",&akq_imag_in);

  TH1F akq_real_in_hist("akq_real_in_hist","akq_real_in_hist",(k_max+1)*(k_max+1),0,(k_max+1)*(k_max+1)); akq_real_in_hist.SetStats(0);
  TH1F akq_imag_in_hist("akq_imag_in_hist","akq_imag_in_hist",(k_max+1)*(k_max+1),0,(k_max+1)*(k_max+1)); akq_imag_in_hist.SetStats(0);
  
  std::complex<double> **akq = new std::complex<double>*[k_max+1];
  int ent=0;
  for (int k=0; k<=k_max; k++) {
    akq[k] = new std::complex<double>[2*k+1];
    for (int q=0; q<2*k+1; q++) {
      intree->GetEntry(ent);
      akq[k][q] = std::complex<double>(akq_real_in,akq_imag_in);
      akq_real_in_hist.Fill(ent,akq_real_in);
      akq_imag_in_hist.Fill(ent,akq_imag_in);
      akq_real_in_hist.GetXaxis()->SetBinLabel(ent+1,Form("%i%i",k_in,q_in));
      akq_imag_in_hist.GetXaxis()->SetBinLabel(ent+1,Form("%i%i",k_in,q_in));
      ent++;
    }
  }
  
  double dxyang = 2.*M_PI/NUM_XY_ANGS;
  double dzang  = 2./NUM_Z_ANGS;
  
  double xy_angle_arr[NUM_XY_ANGS];
  double z_angle_arr[NUM_Z_ANGS];
  double W_arr[NUM_XY_ANGS][NUM_Z_ANGS];

  for (int i=0; i<NUM_XY_ANGS; i++)
    xy_angle_arr[i] = -M_PI + (i+0.5)*dxyang;
  for (int i=0; i<NUM_Z_ANGS; i++)
    z_angle_arr[i] = -1 + (i+0.5)*dzang;

  double W_xyproj_arr[NUM_XY_PROJS][NUM_XY_ANGS];
  double W_zproj_arr[NUM_Z_PROJS][NUM_Z_ANGS];
  
  TH2F input_hist("input_hist","input_hist",NUM_XY_ANGS,-M_PI,M_PI,NUM_Z_ANGS,-1,1);

  cout << "done." << endl
       << "Calculating projections and writing to files...";

  ofstream proj_outfile;
  proj_outfile.open(proj_outfile_name);

  double norm = 4*M_PI*real(akq[0][0]*complex_spherical_harmonic_cosine(0,0,0,0));
  
  double W_moment;
  for (int xyang=0; xyang<NUM_XY_ANGS; xyang++) {
    for (int zang=0; zang<NUM_Z_ANGS; zang++) {
      for (int k=0; k<=k_max; k++) {
	for (int q=-k; q<=k; q++) {
	  W_moment = real(akq[k][q+k]*conj(complex_spherical_harmonic_cosine(k,q,z_angle_arr[zang],xy_angle_arr[xyang])))/norm;
	  input_hist.Fill(xy_angle_arr[xyang],z_angle_arr[zang],W_moment);
	  W_arr[xyang][zang] += W_moment;
	}
      }
    }
  }

  proj_outfile << "NUM_XY_PROJS: " << NUM_XY_PROJS << endl
	       << endl;
  for (int i=0; i<NUM_XY_PROJS; i++) {
    proj_outfile << "  COS_ZANG: " << xy_proj_arr[i] << " | NUM_XY_ANGS: " << NUM_XY_ANGS << endl;
    for (int xyang=0; xyang<NUM_XY_ANGS; xyang++) {
      W_xyproj_arr[i][xyang]=0;
      for (int k=0; k<=k_max; k++) {
	for (int q=-k; q<=k; q++) {
	  W_xyproj_arr[i][xyang] += real(akq[k][q+k]*conj(complex_spherical_harmonic_cosine(k,q,xy_proj_arr[i],xy_angle_arr[xyang])))/norm;
	}
      }
      proj_outfile << "  " << xy_angle_arr[xyang]*RADTODEG << "   " << W_xyproj_arr[i][xyang] << endl;
    }
    proj_outfile << endl;
  }

  proj_outfile << "NUM_Z_PROJS: " << NUM_Z_PROJS << endl
	       << endl;
  for (int i=0; i<NUM_Z_PROJS; i++) {
    proj_outfile << "  XYANG_(deg): " << z_proj_arr[i]*RADTODEG << " | NUM_Z_ANGS: " << NUM_Z_ANGS << endl;
    for (int zang=0; zang<NUM_Z_ANGS; zang++) {
      W_zproj_arr[i][zang]=0;
      for (int k=0; k<=k_max; k++) {
	for (int q=-k; q<=k; q++) {
	  W_zproj_arr[i][zang] += real(akq[k][q+k]*conj(complex_spherical_harmonic_cosine(k,q,z_angle_arr[zang],z_proj_arr[i])))/norm;
	}
      }
      proj_outfile << "  " << z_angle_arr[zang] << "   " << W_zproj_arr[i][zang] << endl;
    }
    proj_outfile << endl;
  }

  proj_outfile.close();
  
  output_rootfile.cd();

  akq_real_in_hist.Write();
  akq_imag_in_hist.Write();
  input_hist.Write();
  
  input_rootfile.Close();
  output_rootfile.Close();

  for (int k=0; k<=k_max; k++)
    delete [] akq[k];
  delete [] akq;
  
  cout << "done." << endl;
  
  return 0;

}
