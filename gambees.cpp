/*

  bees.cpp: generates orientation tensors Bkq from fortran output fort.37 file
            then folds these with Ak coefficients and spherical harmonics to
            produce double-differential angular distributions for gamma
            decay from a state oriented by the fresco reaction.

	    Then these distributions (saved in various forms) can be used as inputs
            for simulations, projection analysis, decompositions, etc.

	    Ken H
	    Aug 2022

*/

#include <iostream>
#include <fstream>
#include <cstdlib>
#include <cstring>
#include <chrono>
#include <complex>

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
#include "WikoGammaClass.h"
#include "BeesHolder.h"
#include "wiko_wiko5.hpp"
#include "FrameConverter.h"

using namespace std;

void SetNucLabel(char *nuc_label, int Z, int A);

int main(int argc, char* argv[]) {

  if (argc != 2) {
    cout << "***WARNING: input file as argument" << endl;
    return 1;
  }
  
  ifstream input_file;
  input_file.open(argv[1]);
  if (!input_file) {
    cout << "***WARNING: cannot find input file " << argv[1] << endl;
    return 2;
  }

  chrono::high_resolution_clock::time_point t_start = chrono::high_resolution_clock::now();

  //P(T,R)E | R --> D1 + D2 (D1=gamma by default)
  int ZT,AT, ZP,AP, ZE,AE, ZR,AR, ZD1,AD1, ZD2,AD2;

  char T_label[10], P_label[10], E_label[10], R_label[10], D1_label[10], D2_label[10];
  
  //for kinematics, frame conversions, etc.
  double beam_E;
  double rxnang_lab;
  double rxnang_rxncm; //will be calculated from rxnang_lab and matched to nearest Bkq index
  double P_exc;
  double T_exc;
  double R_exc;
  double E_exc;
  double D1_exc;
  double D2_exc;

  char fort37_name[100];
  
  double Ji, Jf; //initial and final spin of decay transition
  int    Pi, Pf; //corresponding parities (+ or - 1)
  double delta21, delta31; //multipole mixing ratios
  double euler_alpha; //Euler angle (CM deg.)
  double euler_beta;  //Euler angle (CM deg.)
  double euler_gamma; //Euler angle (CM deg.)

  char canv_label[100];
  char outfile_name[100];
  char frame_label[100];
  char outfile_label[100];
  char label[100];

  input_file >> label >> ZT
	     >> label >> AT
	     >> label >> ZP
	     >> label >> AP
	     >> label >> ZE
	     >> label >> AE

	     >> label >> beam_E
	     >> label >> rxnang_lab
	     >> label >> P_exc
	     >> label >> T_exc
	     >> label >> R_exc
	     >> label >> E_exc
	     >> label >> D1_exc
	     >> label >> D2_exc

	     >> label >> fort37_name
    
             >> label >> Pi
	     >> label >> Jf
	     >> label >> Pf
	     >> label >> delta21
	     >> label >> delta31

	     >> label >> euler_alpha
	     >> label >> euler_beta
	     >> label >> euler_gamma

	     >> label;

  input_file.get(); //space
  input_file.getline(canv_label,100);

  input_file >> label >> outfile_label;

  input_file.close();

  ZR = ZP + ZT - ZE;
  AR = AP + AT - AE;
  ZD1 = 0;
  AD1 = 0;
  ZD2 = ZR;
  AD2 = AR;

  SetNucLabel(T_label,ZT,AT);
  SetNucLabel(P_label,ZP,AP);
  SetNucLabel(E_label,ZE,AE);
  SetNucLabel(R_label,ZR,AR);
  SetNucLabel(D1_label,ZD1,AD1);
  SetNucLabel(D2_label,ZD2,AD2);
  
  rxnang_lab *= DEGTORAD;

  euler_alpha *= DEGTORAD;
  euler_beta  *= DEGTORAD;
  euler_gamma *= DEGTORAD;
  
  FrameConverter fc(ZP,AP, ZT,AT, ZE,AE, ZD1,AD1);
  //calculate first only for rxn params
  fc.Calculate_DecCM(beam_E,rxnang_lab,M_PI,0,0,P_exc,T_exc,R_exc,E_exc,D1_exc,D2_exc);
  rxnang_rxncm = fc.ejec_rxncm.theta;
  
  strcpy(frame_label,"DecCM");

  BeesHolder Bkq(fort37_name,
		 euler_alpha,
		 euler_beta,
		 euler_gamma);
  if (!(Bkq.File_Loaded())) {
    cout << "***WARNING: error loading fort.37 file\n";
    return 10;
  }
  cout << "Loading and generating Bkq's...";
  if (!(Bkq.Generate_Tensors())) {
    cout << "***WARNING: error loading Bkq tensors\n";
    return 11;
  }
  cout << "done.\n";

  Ji = Bkq.J();
  
  int rxn_angle_index;
  for (int i=0; i<Bkq.Num_Angles(); i++) {
    if (rxnang_rxncm > Bkq.Angle(i) && rxnang_rxncm < Bkq.Angle(i+1)) {
      rxn_angle_index = i;
      break;
    }
  }
  
  WikoGammaClass Wiko(Pi*Ji,Pf*Jf,delta21,delta31);

  int k_max = min(Bkq.Max_Rank(),Wiko.Get_Max_Rank());

  char Pi_char = (Pi==1?'+':'-');
  char Pf_char = (Pf==1?'+':'-');

  int NUM_RXN_ANGS = Bkq.Num_Angles(); //number of reaction angles
  int NUM_Z_ANGS   = 200; //number of theta angles (also of cosines)
  int NUM_XY_ANGS  = 360; //number of phi angles
  
  const int PRINT_INPUTS = false; //usually for debugging input issues

  if (PRINT_INPUTS) {
    cout << "Read in from 37 file and input file:"        << endl
	 << "\tZT = "             << ZT                   << endl
	 << "\tAT = "             << AT                   << endl
	 << "\tZP = "             << ZP                   << endl
	 << "\tAP = "             << AP                   << endl
	 << "\tZE = "             << ZE                   << endl
	 << "\tAE = "             << AE                   << endl
	 << "\tbeam_E = "         << beam_E               << endl
	 << "\trxnang_lab = "     << rxnang_lab*RADTODEG  << endl
	 << "\tP_exc = "          << P_exc                << endl
	 << "\tT_exc = "          << T_exc                << endl
	 << "\tR_exc = "          << R_exc                << endl
	 << "\tE_exc = "          << E_exc                << endl
	 << "\tD1_exc = "         << D1_exc               << endl
	 << "\tD2_exc = "         << D2_exc               << endl
	 << "\tJi = "             << Ji                   << endl
	 << "\tPi = "             << Pi                   << endl
	 << "\tJf = "             << Jf                   << endl
	 << "\tPf = "             << Pf                   << endl
	 << "\tdelta21 = "        << delta21              << endl
      	 << "\tdelta31 = "        << delta31              << endl
	 << "\teuler_alpha = "    << euler_alpha*RADTODEG << endl
	 << "\teuler_beta  = "    << euler_beta*RADTODEG  << endl
	 << "\teuler_gamma = "    << euler_gamma*RADTODEG << endl
	 << "\tcanv_label = "     << canv_label           << endl
	 << "\tfort37_name = "    << fort37_name          << endl
	 << "\toutfile_label = "  << outfile_label        << endl;
  }

  strcpy(outfile_name,"rootfiles/");
  strcat(outfile_name,outfile_label);
  strcat(outfile_name,".root");
  
  TFile beesfile(outfile_name,"RECREATE");
  
  //input 0 --> no label
  if (strcmp(canv_label,"0") == 0)
    strcpy(canv_label,"");
  
  TPaveLabel title_pad(0.10,0.95,0.90,0.99,canv_label);

  /////// M-SUBSTATE DIST ///////////
  //--> note: this is calculated assuming rxnang=0
  //for off-axis orientation behavior examine the full density matrices
  int num_msubs = Bkq.Num_MSubs();
  double m_min = -Ji;
  double m_max =  Ji;
  double *msubs_frac = Bkq.MSubs_Integ_Frac();
  TString hist_title = Form("Residual Substate Distribution (J = %i/2^{%c});m;Fraction",(int)(2*Ji),Pi_char);
  TH1F msubs_hist("msubs_hist",hist_title,num_msubs,m_min,m_max+1);
    for (int i=0; i<num_msubs; i++) {
    double mval = i-Ji;
    msubs_hist.SetBinContent(i+1,msubs_frac[i]);
    msubs_hist.GetXaxis()->SetBinLabel(i+1,Form("%i/2",(int)(2*mval)));
  }
  msubs_hist.SetFillColor(17);
  msubs_hist.GetYaxis()->SetRangeUser(0,1);
  msubs_hist.SetStats(0);
  TCanvas msubs_canv("msubs_canv");
  TPad msubs_pad("msubs_pad","msubs_pad",0.05,0.05,0.95,0.94);
  msubs_pad.cd();
  msubs_hist.Draw("bar");
  msubs_canv.cd();
  title_pad.Draw();
  msubs_pad.Draw();
  /////////////////////////////////

  //density matrices
  TH2F densmat_real_hist(Form("densmat_real_hist"),
			 Form("Re(#rho(J = %i/2^{%c}))",(int)(2*Ji),Pi_char),
			 num_msubs,0,num_msubs,
			 num_msubs,0,num_msubs);
    densmat_real_hist.SetStats(0);
  
  TH2F densmat_imag_hist(Form("densmat_imag_hist"),
			 Form("Im(#rho(J = %i/2^{%c}))",(int)(2*Ji),Pi_char),
			 num_msubs,0,num_msubs,
			 num_msubs,0,num_msubs);
  densmat_imag_hist.SetStats(0);
  
  for (int m1=0; m1<num_msubs; m1++) {
    for (int m2=0; m2<num_msubs; m2++) {
      
      densmat_real_hist.Fill(m1,m2,Bkq.DensMat_Real(m1-Ji,m2-Ji,rxn_angle_index));
      densmat_real_hist.GetXaxis()->SetBinLabel(m1+1,Form("%.1f",m1-Ji));
      densmat_real_hist.GetYaxis()->SetBinLabel(m2+1,Form("%.1f",m2-Ji));
      
      densmat_imag_hist.Fill(m1,m2,Bkq.DensMat_Imag(m1-Ji,m2-Ji,rxn_angle_index));
      densmat_imag_hist.GetXaxis()->SetBinLabel(m1+1,Form("%.1f",m1-Ji));
      densmat_imag_hist.GetYaxis()->SetBinLabel(m2+1,Form("%.1f",m2-Ji));
      
    }
  }
  
  TCanvas densmat_real_canv("densmat_real_canv");
  TPad    densmat_real_pad("densmat_real_pad","densmat_real_pad",0.05,0.05,0.95,0.94);
  densmat_real_pad.cd();
  densmat_real_hist.Draw("text");
  densmat_real_canv.cd();
  title_pad.Draw();
  densmat_real_pad.Draw();
    
  TCanvas densmat_imag_canv("densmat_imag_canv");
  TPad    densmat_imag_pad("densmat_imag_pad","densmat_imag_pad",0.05,0.05,0.95,0.94);
  densmat_imag_pad.cd();
  densmat_imag_hist.Draw("text");
  densmat_imag_canv.cd();
  title_pad.Draw();
  densmat_imag_pad.Draw();
  
  cout << "READING FROM FILE: " << fort37_name << endl
       << "Reaction: (" << ZT << "," << AT << ") + (" << ZP << "," << AP << ") --> ("
                        << ZR << "," << AR << ") + (" << ZE << "," << AE << ") | ("
       << ZR << "," << AR << ") --> (" << ZD1 << "," << AD1 << ") + (" << ZD2 << "," << AD2 << ") @ "
       << beam_E << " MeV (lab) beam (" << ZP << "," << AP << ") energy" << endl
       << ">> Fixed lab angle = " << rxnang_lab*RADTODEG << " deg (rxncm " << rxnang_rxncm*RADTODEG << " deg, index " << rxn_angle_index << ")" << endl
       << ">> Projectile excitation = " << P_exc << " MeV" << endl
       << ">> Target excitation = " << T_exc << " MeV" << endl
       << ">> Residual excitation = " << R_exc << " MeV" << endl
       << ">> Ejectile excitation = " << E_exc << " MeV" << endl
       << ">> Light Decay excitation = " << D1_exc << " MeV" << endl
       << ">> Heavy Decay excitation = " << D2_exc << " MeV" << endl
       << "Euler angles (Rxn-CM deg.): alpha = " << euler_alpha*RADTODEG
       <<                          " | beta = " << euler_beta*RADTODEG
       <<                          " | gamma = " << euler_gamma*RADTODEG << endl
       << ">> transition from Ji = " << Ji << Pi_char << " to Jf = " << Jf << Pf_char << " (L=" << Wiko.LP() << "/L=" << Wiko.L() << " = " << delta21 << ", L=" << Wiko.LPP() << "/L=" << Wiko.L() << " = " << delta31 << ")" << endl
       << ">> MAX RANK: k = " << k_max << endl
       << "Calculating in " << frame_label << " frame" << endl
       << "----------------------------------------------------" << endl
       << ">> Creating root file " << outfile_name << endl;

  TTree input_tree("input_tree","input_tree"); //to record all inputs
  input_tree.Branch("ZT",&ZT);
  input_tree.Branch("AT",&AT);
  input_tree.Branch("ZP",&ZP);
  input_tree.Branch("AP",&AP);
  input_tree.Branch("ZE",&ZE);
  input_tree.Branch("AE",&AE);
  input_tree.Branch("ZR",&ZR);
  input_tree.Branch("AR",&AR);
  input_tree.Branch("ZD1",&ZD1);
  input_tree.Branch("AD1",&AD1);
  input_tree.Branch("ZD2",&ZD2);
  input_tree.Branch("AD2",&AD2);
  input_tree.Branch("beam_E_MeV",&beam_E);
  input_tree.Branch("rxnang_lab_rad",&rxnang_lab);
  input_tree.Branch("rxnang_rxncm_rad",&rxnang_rxncm);
  input_tree.Branch("rxn_angle_index",&rxn_angle_index);
  input_tree.Branch("P_exc_MeV",&P_exc);
  input_tree.Branch("T_exc_MeV",&T_exc);
  input_tree.Branch("R_exc_MeV",&R_exc);
  input_tree.Branch("E_exc_MeV",&E_exc);
  input_tree.Branch("D1_exc_MeV",&D1_exc);
  input_tree.Branch("D2_exc_MeV",&D2_exc);
  input_tree.Branch("Ji",&Ji);
  input_tree.Branch("Pi",&Pi);
  input_tree.Branch("Jf",&Jf);
  input_tree.Branch("Pf",&Pf);
  input_tree.Branch("delta21",&delta21);
  input_tree.Branch("delta31",&delta31);
  input_tree.Branch("k_max",&k_max);
  input_tree.Branch("euler_alpha_rad",&euler_alpha);
  input_tree.Branch("euler_beta_rad",&euler_beta);
  input_tree.Branch("euler_gamma_rad",&euler_gamma);
  input_tree.Branch("NUM_XY_ANGS",&NUM_XY_ANGS);
  input_tree.Branch("NUM_Z_ANGS",&NUM_Z_ANGS);
  input_tree.Fill();
  
  //integration steps
  double drxnang   = Bkq.DTheta();
  double dzang_cos = 2./NUM_Z_ANGS;
  double dxyang    = 360.*DEGTORAD/NUM_XY_ANGS;

  char xy_angle_label[10], z_angle_label[10];
  if (euler_alpha==0 && euler_beta==0 && euler_gamma==0) {
    strcpy(z_angle_label,"cos#theta");
    strcpy(xy_angle_label,"#phi");
  }
  else {
    strcpy(z_angle_label,"cos#chi");
    strcpy(xy_angle_label,"#psi");
  }
  
  //grab xsec and angles
  double rxnangle[NUM_RXN_ANGS], xsec[NUM_RXN_ANGS];
  for (int i=0; i<NUM_RXN_ANGS; i++) {
    rxnangle[i] = Bkq.Angle(i)*RADTODEG;
    //this procedure is to ensure angle is precisely integer or half-integer
    //for later checks
    rxnangle[i] *= 2;
    rxnangle[i] = round(rxnangle[i]);
    rxnangle[i] /= 2;
    xsec[i] = Bkq.XSec(i);
  }

  //decay angle
  double z_angle_cosine[NUM_Z_ANGS];
  for (int i=0; i<NUM_Z_ANGS; i++) {
    z_angle_cosine[i] = -1+(i+0.5)*dzang_cos;
  }

  //relative azimuthal angle
  double xy_angle[NUM_XY_ANGS];
  for (int i=0; i<NUM_XY_ANGS; i++)
    xy_angle[i] = -M_PI+(i+0.5)*dxyang;

  //the (real) ang dist itself
  double W[NUM_XY_ANGS][NUM_Z_ANGS];
  for (int i=0; i<NUM_XY_ANGS; i++)
    for (int j=0; j<NUM_Z_ANGS; j++)
      W[i][j]=0;

  TTree data_tree("data_tree","data_tree");
  double xy_angle_fill, z_angle_fill, W_fill;
  data_tree.Branch("xy_angle",&xy_angle_fill);
  data_tree.Branch("z_angle",&z_angle_fill);
  data_tree.Branch("W",&W_fill);

  TTree coeff_tree("coeff_tree","coeff_tree");
  int k_fill, q_fill;
  double akq_real_fill, akq_imag_fill;
  coeff_tree.Branch("k_fill",&k_fill);
  coeff_tree.Branch("q_fill",&q_fill);
  coeff_tree.Branch("akq_real_fill",&akq_real_fill);
  coeff_tree.Branch("akq_imag_fill",&akq_imag_fill);
  
  //graphs and canvases for Bkq only
  TGraph ***re_graphs = new TGraph**[k_max+1];
  TCanvas **re_canvs  = new TCanvas*[k_max+1];
  TPad    **re_pads   = new TPad*[k_max+1];

  TGraph ***im_graphs = new TGraph**[k_max+1];
  TCanvas **im_canvs  = new TCanvas*[k_max+1];
  TPad    **im_pads   = new TPad*[k_max+1];

  for (int i=0; i<=k_max; i++) {

    re_graphs[i] = new TGraph*[2*k_max+1];
    im_graphs[i] = new TGraph*[2*k_max+1];

    re_canvs[i] = new TCanvas(Form("B%i_canv_real",i),Form("B%i_canv_real",i));
    im_canvs[i] = new TCanvas(Form("B%i_canv_imag",i),Form("B%i_canv_imag",i));

    re_pads[i] = new TPad(Form("B%i_pad_real",i),Form("B%i_pad_real",i),0.05,0.05,0.95,0.94);
    im_pads[i] = new TPad(Form("B%i_pad_imag",i),Form("B%i_pad_imag",i),0.05,0.05,0.95,0.94);

    int rows = 0, cols = 0; //num rows and columns in which to divide canvases

    switch (i) { //can get creative w canvas divides (make it look pretty)
    case 0:  rows = 1; cols = 1; break;
    case 1:  rows = 3; cols = 1; break;
    case 2:  rows = 3; cols = 2; break;
    case 3:  rows = 4; cols = 2; break;
    case 4:  rows = 3; cols = 3; break;
    case 5:  rows = 4; cols = 3; break;
    case 6:  rows = 5; cols = 3; break;
    case 7:  rows = 5; cols = 3; break;
    case 8:  rows = 6; cols = 3; break;
    default: rows = i; cols = 3; break;
    }

    re_pads[i]->Divide(cols,rows);
    im_pads[i]->Divide(cols,rows);

  }

  for (int k=0; k<=k_max; k++) {
    for (int q=-k; q<=k; q++) {

      re_graphs[k][q+k] = new TGraph(NUM_RXN_ANGS, rxnangle, Bkq.Bkq_Real(k,q));
      re_graphs[k][q+k]->SetTitle(Form("Re(B_{%i%i});#theta_{rxn}^{CM} (deg.);Amplitude (fm^{2})",k,q));
      re_pads[k]->cd(q+k+1);
      re_graphs[k][q+k]->Draw();

      im_graphs[k][q+k] = new TGraph(NUM_RXN_ANGS, rxnangle, Bkq.Bkq_Imag(k,q));
      im_graphs[k][q+k]->SetTitle(Form("Im(B_{%i%i});#theta_{rxn}^{CM} (deg.);Amplitude (fm^{2})",k,q));
      im_pads[k]->cd(q+k+1);
      im_graphs[k][q+k]->Draw();
    
    } //end q

    re_canvs[k]->cd();
    title_pad.Draw();
    re_pads[k]->Draw();

    im_canvs[k]->cd();
    title_pad.Draw();
    im_pads[k]->Draw();

  } //end k

  //then calculate decay parts

  double norm = 4*M_PI*Wiko.Decay_Coefficient(0)*Wiko.Spherical_Harmonic(0,0,0,0);

  /* 1D HISTOGRAMS*/
  TH1F Ak_hist(Form("Ak_hist"),Form("Coupling Coefficients (A_{k})"),k_max+1,0,k_max+1);
  for (int k=0; k<=k_max; k++) {
    Ak_hist.Fill(k,Wiko.Decay_Coefficient(k));
    Ak_hist.GetXaxis()->SetBinLabel(k+1,Form("%i",k));
  }

  /* 2D HISTOGRAMS */
  TH2F DecayDist_real_hist("DecayDist_real_hist",
			   Form("W(%s,%s)_{%s} (#theta_{rxn}^{LAB} = %.1f^{o})",
				z_angle_label,xy_angle_label,frame_label,(double)(rxnang_lab*RADTODEG)),
			   NUM_XY_ANGS,-180,180,NUM_Z_ANGS,-1,1);
  DecayDist_real_hist.GetXaxis()->SetTitle(Form("%s_{%s}^{%s} (deg.)",xy_angle_label,D1_label,frame_label));
  DecayDist_real_hist.GetYaxis()->SetTitle(Form("%s_{%s}^{%s} (unitless)",z_angle_label,D1_label,frame_label));
  DecayDist_real_hist.SetStats(0);

  TH2F DecayDist_imag_hist("DecayDist_imag_hist",
			   Form("W(%s,%s)_{%s} (#theta_{rxn}^{LAB} = %.1f^{o})",
				z_angle_label,xy_angle_label,frame_label,(double)(rxnang_lab*RADTODEG)),
			   NUM_XY_ANGS,-180,180,NUM_Z_ANGS,-1,1);
  DecayDist_imag_hist.GetXaxis()->SetTitle(Form("%s_{%s}^{%s} (deg.)",xy_angle_label,D1_label,frame_label));
  DecayDist_imag_hist.GetYaxis()->SetTitle(Form("%s_{%s}^{%s} (unitless)",z_angle_label,D1_label,frame_label));
  DecayDist_imag_hist.SetStats(0);

  //variables for total integration
  double Wtot=0;
  
  /* MAIN W LOOP */
  for (int k=0; k<=k_max; k++) {
    for (int q=-k; q<=k; q++) {
      cout << "\rCalculating W(k,q) = W(" << k << "," << q << ")..." << flush;
      for (int zang=0; zang<NUM_Z_ANGS; zang++) { //decay angle	  
	for (int xyang=0; xyang<NUM_XY_ANGS; xyang++) { //relative azimuthal angle

	  /***** MAIN W CALCULATION *****/
	  double W_real = real(Wiko.Decay_Coefficient(k)*
			       Bkq.Bkq_Complex(k,q,rxn_angle_index)*
			       conj(Wiko.Complex_Spherical_Harmonic_Cosine(k,q,z_angle_cosine[zang],xy_angle[xyang])))/norm;
	  double W_imag = imag(Wiko.Decay_Coefficient(k)*
			       Bkq.Bkq_Complex(k,q,rxn_angle_index)*
			       conj(Wiko.Complex_Spherical_Harmonic_Cosine(k,q,z_angle_cosine[zang],xy_angle[xyang])))/norm;
	  /*****************************/

	  if (k%2!=0) { //***ENFORCING k=even by helicity summation in Hamilton
	    W_real=0;
	    W_imag=0;
	  }
	  
	  W[xyang][zang] += W_real;
	  
	  //total integration
	  Wtot += W_real*dzang_cos*dxyang;
	  
	  DecayDist_real_hist.Fill(xy_angle[xyang]*RADTODEG,z_angle_cosine[zang],W_real);
	  DecayDist_imag_hist.Fill(xy_angle[xyang]*RADTODEG,z_angle_cosine[zang],W_imag);
	  
	} //end xyang loop
      } //end zang loop

      k_fill = k;
      q_fill = q;
      akq_real_fill = k%2==0?real(Wiko.Decay_Coefficient(k)*Bkq.Bkq_Complex(k,q,rxn_angle_index)):0;
      akq_imag_fill = k%2==0?imag(Wiko.Decay_Coefficient(k)*Bkq.Bkq_Complex(k,q,rxn_angle_index)):0;
      coeff_tree.Fill();
      
    } //end q loop
  } //end k loop

  for (int zang=0; zang<NUM_Z_ANGS; zang++) {
    for (int xyang=0; xyang<NUM_XY_ANGS; xyang++) {
      xy_angle_fill = xy_angle[xyang];
      z_angle_fill  = z_angle_cosine[zang];
      W_fill        = W[xyang][zang];
      data_tree.Fill();
    }
  }
  
  cout << "done.\nTotal W = " << Wtot << endl;
  
  /* DATA VISUALIZATION */

  //ordinary reaction xsec (no decay)
  TGraph xsec_gr(NUM_RXN_ANGS,rxnangle,xsec);
  xsec_gr.SetTitle("Reaction Cross Section;#theta_{rxn}^{CM} (deg.);d#sigma/d#Omega (mb/sr)");
  xsec_gr.SetName("xsec_gr");
  TCanvas xsec_canv("xsec_canv");
  TPad    xsec_pad("xsec_pad","xsec_pad",0.05,0.05,0.95,0.94);
  xsec_pad.cd();
  xsec_gr.Draw();
  xsec_pad.SetLogy();
  xsec_canv.cd();
  title_pad.Draw();
  xsec_pad.Draw();

  TCanvas Ak_canv("Ak_canv");
  TPad    Ak_pad("Ak_pad","Ak_pad",0.05,0.05,0.95,0.94);
  Ak_pad.cd();
  Ak_hist.Draw("bar");
  Ak_canv.cd();
  title_pad.Draw();
  Ak_pad.Draw();
  
  TCanvas DecayDist_real_canv("DecayDist_real_canv");
  TPad    DecayDist_real_pad("DecayDist_real_pad","DecayDist_real_pad",0.05,0.05,0.95,0.94);
  DecayDist_real_pad.cd();
  DecayDist_real_hist.Draw("colz");
  DecayDist_real_canv.cd();
  title_pad.Draw();
  DecayDist_real_pad.Draw();

  TCanvas DecayDist_imag_canv("DecayDist_imag_canv");
  TPad    DecayDist_imag_pad("DecayDist_imag_pad","DecayDist_imag_pad",0.05,0.05,0.95,0.94);
  DecayDist_imag_pad.cd();
  DecayDist_imag_hist.Draw("colz");
  DecayDist_imag_canv.cd();
  title_pad.Draw();
  DecayDist_imag_pad.Draw();

  beesfile.cd();

  input_tree.Write();
  data_tree.Write();
  coeff_tree.Write();

  Ak_canv.Write();
  xsec_canv.Write();
  msubs_canv.Write();
  densmat_real_canv.Write();
  densmat_imag_canv.Write();
  for (int i=0; i<=k_max; i++) {
    re_canvs[i]->Write();
    im_canvs[i]->Write();
  }
  DecayDist_real_canv.Write();
  DecayDist_imag_canv.Write();

  beesfile.Close();
  
  for (int i=0; i<=k_max; i++) {
    delete [] re_graphs[i];
    delete [] im_graphs[i];
  }
  delete [] re_graphs;
  delete [] re_canvs;
  delete [] re_pads;
  delete [] im_graphs;
  delete [] im_canvs;
  delete [] im_pads;
  
  chrono::high_resolution_clock::time_point t_stop = chrono::high_resolution_clock::now();
  chrono::duration<double> t_span = chrono::duration_cast<chrono::duration<double>>(t_stop - t_start);

  cout << "Total time: " << t_span.count() << " seconds (" << t_span.count()/60. << " minutes)\n";

  return 0;

}

void SetNucLabel(char *nuc_label, int Z, int A) {

  if      (Z==0 && A==0) strcpy(nuc_label,"#gamma");
  else if (Z==0 && A==1) strcpy(nuc_label,"n");
  else if (Z==1 && A==1) strcpy(nuc_label,"p");
  else if (Z==1 && A==2) strcpy(nuc_label,"d");
  else if (Z==1 && A==3) strcpy(nuc_label,"t");
  else if (Z==2 && A==4) strcpy(nuc_label,"#alpha");
  else {
    switch (Z) {
    case 1:  strcpy(nuc_label,Form("^{%i}H",A));  break;
    case 2:  strcpy(nuc_label,Form("^{%i}He",A)); break;
    case 3:  strcpy(nuc_label,Form("^{%i}Li",A)); break;
    case 4:  strcpy(nuc_label,Form("^{%i}Be",A)); break;
    case 5:  strcpy(nuc_label,Form("^{%i}B",A));  break;
    case 6:  strcpy(nuc_label,Form("^{%i}C",A));  break;
    case 7:  strcpy(nuc_label,Form("^{%i}N",A));  break;
    case 8:  strcpy(nuc_label,Form("^{%i}O",A));  break;
    case 9:  strcpy(nuc_label,Form("^{%i}F",A));  break;
    case 10: strcpy(nuc_label,Form("^{%i}Ne",A)); break;
    default: strcpy(nuc_label,"xx"); break;
    }
  }
  
}
