/*

  geebees.cpp: generates and couples orientation tensors Bkq from multiple fresco
               output fort.37 files, then folds with spherical harmonics to
               produce double-differential angular distributions for particle
               decay from states oriented by the fresco reaction.

	       Then these distributions (saved in various forms) can be used as inputs
               for simulations, projection analysis, decompositions, etc.

               The intention is for this to be used with measured decays from
               overlapping unbound states.

	       ASSUMPTIONS:
	       * spin 1/2 particle decay (proton or neutron)
               * decay down to a 0+ final state of the heavy decay residual
               * states are long-lived (narrow & appropriate for weak binding approx.)

	    Ken H
	    Feb 2022

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
#include <TLegend.h>

#include "constants.h"
#include "GrandBeesHolder.h"
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

  //P(T,R)E | R --> D1 + D2
  int ZT,AT, ZP,AP, ZE,AE, ZR,AR, ZD1,AD1, ZD2,AD2;
  double MD1, MD2; //GDM inputs (masses in MeV)

  char T_label[10], P_label[10], E_label[10], R_label[10], D1_label[10], D2_label[10];
  
  //for kinematics, frame conversions, etc.
  double beam_E;
  double rxnang_lab;
  double rxnang_rxncm; //will be calculated from rxnang_lab and matched to nearest Bkq index
                       // --> since these states are NARROW and OVERLAPPING, what works for
                       //     one ought to work for all (similar kinematics)
  
  double P_exc;
  double T_exc;
  double E_exc;
  double D1_exc;
  double D2_exc;

  int NUM_STATES;

  int    *Pi; //corresponding parities (+ or - 1) (Pf=+1)
  double *Ji; //initial spin of decaying state (from fort37 file)
  double *R_exc, *R_res; //GDM class expects exc relative to threshold
  double *R_wid; //physical widths in MeV for the energy averaging
  double *BR; //branching ratio
  double *E12; //relative decay KE (calculated from the above)
  char   **fort37_name;
  char   **state_label;

  //BIG ASSUMPTIONS
  double Jf=0.;
  int    Pf=+1;
  
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
	     >> label >> ZD1
	     >> label >> AD1

	     >> label >> beam_E
	     >> label >> rxnang_lab
    
	     >> label >> P_exc
	     >> label >> T_exc
	     >> label >> E_exc
	     >> label >> D1_exc
	     >> label >> D2_exc

	     >> label >> NUM_STATES

	     >> label >> label >> label >> label >> label >> label;

  Pi    = new int[NUM_STATES];
  Ji    = new double[NUM_STATES];
  R_exc = new double[NUM_STATES];
  R_res = new double[NUM_STATES];
  R_wid = new double[NUM_STATES];
  BR    = new double[NUM_STATES];
  E12   = new double[NUM_STATES];
  fort37_name = new char*[NUM_STATES];
  state_label = new char*[NUM_STATES];
  
  for (int i=0; i<NUM_STATES; i++) {
    fort37_name[i] = new char[100];
    state_label[i] = new char[100];
    input_file >> Pi[i] >> R_exc[i] >> R_wid[i] >> BR[i] >> fort37_name[i];
    input_file.get();
    input_file.getline(state_label[i],100);
  }
  
  input_file >> label >> euler_alpha
	     >> label >> euler_beta
	     >> label >> euler_gamma
    
	     >> label;

  input_file.get(); //space
  input_file.getline(canv_label,100);

  input_file >> label >> outfile_label;

  input_file.close();

  const int PRINT_INPUTS = false; //usually for debugging input issues

  if (PRINT_INPUTS) {
    cout << "Read in from 37 file and input file:"           << endl
	 << "\tZT = "             << ZT                      << endl
	 << "\tAT = "             << AT                      << endl
	 << "\tZP = "             << ZP                      << endl
	 << "\tAP = "             << AP                      << endl
	 << "\tZE = "             << ZE                      << endl
	 << "\tAE = "             << AE                      << endl
	 << "\tZD1 = "            << ZD1                     << endl
	 << "\tAD1 = "            << AD1                     << endl
	 << "\tbeam_E = "         << beam_E                  << endl
	 << "\trxnang_lab = "     << rxnang_lab              << endl
	 << "\tP_exc = "          << P_exc                   << endl
	 << "\tT_exc = "          << T_exc                   << endl
	 << "\tE_exc = "          << E_exc                   << endl
	 << "\tD1_exc = "         << D1_exc                  << endl
	 << "\tD2_exc = "         << D2_exc                  << endl
	 << "\tNUM_STATES = "     << NUM_STATES              << endl
	 << "\tJi   Par   Exc   Wid   BR   Filename   Label" << endl;
    for (int i=0; i<NUM_STATES; i++)
      cout << "\t" << Ji[i] << "   " << Pi[i] << "   " << R_exc[i] << "   " << R_wid[i] << "   " << BR[i] << "   " << fort37_name[i] << "   " << state_label[i] << endl;
    cout << "\tJf = "             << Jf            << endl
	 << "\tPf = "             << Pf            << endl
	 << "\teuler_alpha = "    << euler_alpha   << endl
	 << "\teuler_beta  = "    << euler_beta    << endl
	 << "\teuler_gamma = "    << euler_gamma   << endl
	 << "\tcanv_label = "     << canv_label    << endl
	 << "\toutfile_label = "  << outfile_label << endl;
  }
  
  ZR = ZP + ZT - ZE;
  AR = AP + AT - AE;
  ZD2 = ZR - ZD1;
  AD2 = AR - AD1;

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
  double sepE = -fc.Q_dec;
  //calculate first only for rxn params
  //--> again, here we're assuming the states are kinematically similar enough that
  //    the same Bkq index can work for all of them
  //    BUT the decay KE E12 is subtle enough that we need to keep track of it
  for (int i=0; i<NUM_STATES; i++) {
    fc.Calculate_DecCM(beam_E,rxnang_lab,M_PI,0,0,P_exc,T_exc,R_exc[i],E_exc,D1_exc,D2_exc);
    if (i==0) { //(roughly) common to all states
      rxnang_rxncm = fc.ejec_rxncm.theta;
      MD1 = fc.decL_lab.Mtot;
      MD2 = fc.decH_lab.Mtot;
    }
    R_res[i] = R_exc[i] - sepE;
    E12[i] = fc.decL_deccm.KE + fc.decH_deccm.KE; //independent of angle in this frame
  }
  
  strcpy(frame_label,"DecCM");

  GrandBeesHolder geebees(NUM_STATES,fort37_name,
			  R_res,R_wid,BR,Pi,ZD1,ZD2,MD1,MD2,E12,
			  euler_alpha,euler_beta,euler_gamma);

  for (int i=0; i<NUM_STATES; i++)
    geebees.indiv_bees[i]->Toggle_Norm_to_B00();
  
  if (!(geebees.Files_Loaded())) {
    cout << "***WARNING: error loading fort.37 files" << endl;
    return 10;
  }
  cout << "Loading and generating Bkq's...";
  if (!(geebees.Generate_Tensors())) {
    cout << "***WARNING: error loading Bkq tensors" << endl;
    return 11;
  }

  for (int i=0; i<NUM_STATES; i++)
    Ji[i] = geebees.indiv_bees[i]->J();
  
  cout << "done." << endl;
  
  int rxn_angle_index;
  for (int i=0; i<geebees.Num_Angles(); i++) {
    if (rxnang_rxncm > geebees.Angle(i) && rxnang_rxncm < geebees.Angle(i+1)) {
      rxn_angle_index = i;
      break;
    }
  }

  int k_max = geebees.Max_Rank();

  int k_max_indiv[NUM_STATES];
  for (int i=0; i<NUM_STATES; i++)
    k_max_indiv[i] = geebees.indiv_bees[i]->Max_Rank();
  
  char *Pi_char = new char[NUM_STATES];
  for (int i=0; i<NUM_STATES; i++)
    Pi_char[i] = (Pi[i]==1?'+':'-');

  char Pf_char = (Pf==1?'+':'-');

  int NUM_RXN_ANGS = geebees.Num_Angles(); //number of reaction angles
  int NUM_Z_ANGS   = 200; //number of theta angles (also of cosines)
  int NUM_XY_ANGS  = 360; //number of phi angles
  
  strcpy(outfile_name,"rootfiles/");
  strcat(outfile_name,outfile_label);
  strcat(outfile_name,".root");
  
  TFile geebeesfile(outfile_name,"RECREATE");
  
  //input 0 --> no label
  if (strcmp(canv_label,"0") == 0)
    strcpy(canv_label,"");
  
  TPaveLabel title_pad(0.10,0.95,0.90,0.99,canv_label);

  // m-substates of individual files
  //--> note: these is calculated assuming rxnang=0
  //for off-axis orientation behavior (incl coupling) examine the full density matrices
  TH1F    **msubs_hist = new TH1F*[NUM_STATES];
  TCanvas **msubs_canv = new TCanvas*[NUM_STATES];
  TPad    **msubs_pad  = new TPad*[NUM_STATES];
  
  for (int s=0; s<NUM_STATES; s++) {
    int num_msubs = geebees.indiv_bees[s]->Num_MSubs();
    double m_min = -Ji[s];
    double m_max =  Ji[s];
    double *msubs_frac = geebees.indiv_bees[s]->MSubs_Integ_Frac();
    TString hist_title = Form("Residual Substate Distribution (J = %i/2^{%c});m;Fraction",(int)(2*Ji[s]),Pi_char[s]);
    msubs_hist[s] = new TH1F(Form("msubs_%i_hist",s),hist_title,num_msubs,m_min,m_max+1);
    for (int i=0; i<num_msubs; i++) {
      double mval = i-Ji[s];
      msubs_hist[s]->SetBinContent(i+1,msubs_frac[i]);
      msubs_hist[s]->GetXaxis()->SetBinLabel(i+1,Form("%i/2",(int)(2*mval)));
    }
    msubs_hist[s]->SetFillColor(17);
    msubs_hist[s]->GetYaxis()->SetRangeUser(0,1);
    msubs_hist[s]->SetStats(0);
    msubs_canv[s] = new TCanvas(Form("msubs_%i_canv",s));
    msubs_pad[s] = new TPad(Form("msubs_%i_pad",s),"msubs_pad",0.05,0.05,0.95,0.94);
    msubs_pad[s]->cd();
    msubs_hist[s]->Draw("bar");
    msubs_canv[s]->cd();
    title_pad.Draw();
    msubs_pad[s]->Draw();
  }

  //density matrices (blocked)
  TH2F    ***densmat_real_hist = new TH2F**[NUM_STATES];
  TCanvas ***densmat_real_canv = new TCanvas**[NUM_STATES];
  TPad    ***densmat_real_pad  = new TPad**[NUM_STATES];

  TH2F    ***densmat_imag_hist = new TH2F**[NUM_STATES];
  TCanvas ***densmat_imag_canv = new TCanvas**[NUM_STATES];
  TPad    ***densmat_imag_pad  = new TPad**[NUM_STATES];

  for (int s1=0; s1<NUM_STATES; s1++) {

    int num_m1subs = geebees.indiv_bees[s1]->Num_MSubs();
    
    densmat_real_hist[s1] = new TH2F*[NUM_STATES];
    densmat_real_canv[s1] = new TCanvas*[NUM_STATES];
    densmat_real_pad[s1]  = new TPad*[NUM_STATES];
    
    densmat_imag_hist[s1] = new TH2F*[NUM_STATES];
    densmat_imag_canv[s1] = new TCanvas*[NUM_STATES];
    densmat_imag_pad[s1]  = new TPad*[NUM_STATES];

    for (int s2=0; s2<NUM_STATES; s2++) {

      int num_m2subs = geebees.indiv_bees[s2]->Num_MSubs();
      
      densmat_real_hist[s1][s2] = new TH2F(Form("densmat_%i_%i_real_hist",s1,s2),
					   Form("Re(#rho(J = %i/2^{%c},J = %i/2^{%c}))",(int)(2*Ji[s1]),Pi_char[s1],(int)(2*Ji[s2]),Pi_char[s2]),
					   num_m1subs,0,num_m1subs,
					   num_m2subs,0,num_m2subs);
      densmat_real_hist[s1][s2]->SetStats(0);
  
      densmat_imag_hist[s1][s2] = new TH2F(Form("densmat_%i_%i_imag_hist",s1,s2),
					   Form("Im(#rho(J = %i/2^{%c},J = %i/2^{%c}))",(int)(2*Ji[s1]),Pi_char[s1],(int)(2*Ji[s2]),Pi_char[s2]),
					   num_m1subs,0,num_m1subs,
					   num_m2subs,0,num_m2subs);
      densmat_imag_hist[s1][s2]->SetStats(0);
  
      for (int m1=0; m1<num_m1subs; m1++) {
	for (int m2=0; m2<num_m2subs; m2++) {
      
	  densmat_real_hist[s1][s2]->Fill(m1,m2,geebees.DensMat_Real(s1,s2,m1-Ji[s1],m2-Ji[s2],rxn_angle_index));
	  densmat_real_hist[s1][s2]->GetXaxis()->SetBinLabel(m1+1,Form("%.1f",m1-Ji[s1]));
	  densmat_real_hist[s1][s2]->GetYaxis()->SetBinLabel(m2+1,Form("%.1f",m2-Ji[s1]));
      
	  densmat_imag_hist[s1][s2]->Fill(m1,m2,geebees.DensMat_Imag(s1,s2,m1-Ji[s1],m2-Ji[s2],rxn_angle_index));
	  densmat_imag_hist[s1][s2]->GetXaxis()->SetBinLabel(m1+1,Form("%.1f",m1-Ji[s2]));
	  densmat_imag_hist[s1][s2]->GetYaxis()->SetBinLabel(m2+1,Form("%.1f",m2-Ji[s2]));
      
	}
      }
  
      densmat_real_canv[s1][s2] = new TCanvas(Form("densmat_%i_%i_real_canv",s1,s2));
      densmat_real_pad[s1][s2] = new TPad(Form("densmat_%i_%i_real_pad",s1,s2),"densmat_real_pad",0.05,0.05,0.95,0.94);
      densmat_real_pad[s1][s2]->cd();
      densmat_real_hist[s1][s2]->Draw("text");
      densmat_real_canv[s1][s2]->cd();
      title_pad.Draw();
      densmat_real_pad[s1][s2]->Draw();
    
      densmat_imag_canv[s1][s2] = new TCanvas(Form("densmat_%i_%i_imag_canv",s1,s2));
      densmat_imag_pad[s1][s2] = new TPad(Form("densmat_%i_%i_imag_pad",s1,s2),"densmat_imag_pad",0.05,0.05,0.95,0.94);
      densmat_imag_pad[s1][s2]->cd();
      densmat_imag_hist[s1][s2]->Draw("text");
      densmat_imag_canv[s1][s2]->cd();
      title_pad.Draw();
      densmat_imag_pad[s1][s2]->Draw();

    }

  }
  
  cout << "Reaction: (" << ZT << "," << AT << ") + (" << ZP << "," << AP << ") --> ("
                        << ZR << "," << AR << ") + (" << ZE << "," << AE << ") | ("
       << ZR << "," << AR << ") --> (" << ZD1 << "," << AD1 << ") + (" << ZD2 << "," << AD2 << ") @ "
       << beam_E << " MeV (lab) beam (" << ZP << "," << AP << ") energy" << endl
       << ">> Fixed lab angle = " << rxnang_lab*RADTODEG << " deg (rxncm " << rxnang_rxncm*RADTODEG << " deg, index " << rxn_angle_index << ": " << geebees.Angle(rxn_angle_index)*RADTODEG << " deg)" << endl
       << ">> Projectile excitation = " << P_exc << " MeV" << endl
       << ">> Target excitation = " << T_exc << " MeV" << endl
       << ">> Ejectile excitation = " << E_exc << " MeV" << endl
       << ">> Light Decay excitation = " << D1_exc << " MeV" << endl
       << ">> Heavy Decay excitation = " << D2_exc << " MeV" << endl
       << ">> Final state: " << Jf << Pf_char << endl
       << ">> Number of states to couple: " << NUM_STATES << endl;
  for (int i=0; i<NUM_STATES; i++)
    cout << "    " << Ji[i] << Pi_char[i] << " @ " << R_exc[i] << " MeV (res " << R_res[i] << " MeV), width " << R_wid[i] << " MeV (BR " << BR[i] << ") | file: " << fort37_name[i] << endl;
  cout << ">> Huby-Liu weights:  ";
  for (int i=0; i<NUM_STATES; i++) {
    for (int j=0; j<NUM_STATES; j++) {
      cout << geebees.Get_Eta(i,j) << "  ";
    }
  }
  cout << endl
       << "Euler angles (Rxn-CM deg.): alpha = " << euler_alpha*RADTODEG
       <<                          " | beta = " << euler_beta*RADTODEG
       <<                          " | gamma = " << euler_gamma*RADTODEG << endl
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
  input_tree.Branch("E_exc_MeV",&E_exc);
  input_tree.Branch("D1_exc_MeV",&D1_exc);
  input_tree.Branch("D2_exc_MeV",&D2_exc);
  input_tree.Branch("Jf",&Jf);
  input_tree.Branch("Pf",&Pf);
  input_tree.Branch("k_max",&k_max);
  for (int i=0; i<NUM_STATES; i++) {
    input_tree.Branch(Form("R_exc_%i_MeV",i),&R_exc[i]);
    input_tree.Branch(Form("R_res_%i_MeV",i),&R_res[i]);
    input_tree.Branch(Form("R_wid_%i_MeV",i),&R_wid[i]);
    input_tree.Branch(Form("Ji_%i",i),&Ji[i]);
    input_tree.Branch(Form("Pi_%i",i),&Pi[i]);
    input_tree.Branch(Form("BR_%i",i),&BR[i]);
  }
  input_tree.Branch("euler_alpha_rad",&euler_alpha);
  input_tree.Branch("euler_beta_rad",&euler_beta);
  input_tree.Branch("euler_gamma_rad",&euler_gamma);
  input_tree.Branch("NUM_XY_ANGS",&NUM_XY_ANGS);
  input_tree.Branch("NUM_Z_ANGS",&NUM_Z_ANGS);
  input_tree.Fill();
  
  //integration steps (same for all states)
  double drxnang   = geebees.DTheta();
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
  
  //grab xsecs and angles
  double *rxnangle = geebees.Angle();
  double **xsec_indiv = new double*[NUM_STATES];
  double *xsec_coupled;
  //this procedure is to ensure angle is precisely integer or half-integer
  //for later checks
  for (int i=0; i<NUM_RXN_ANGS; i++) {
    rxnangle[i] *= RADTODEG;
    rxnangle[i] *= 2;
    rxnangle[i] = round(rxnangle[i]);
    rxnangle[i] /= 2;
  }
  for (int s=0; s<NUM_STATES; s++)
    xsec_indiv[s] = geebees.indiv_bees[s]->XSec();
  xsec_coupled = geebees.XSec();
  

  //decay angle
  double z_angle_cosine[NUM_Z_ANGS];
  for (int i=0; i<NUM_Z_ANGS; i++) {
    z_angle_cosine[i] = -1+(i+0.5)*dzang_cos;
  }

  //relative azimuthal angle
  double xy_angle[NUM_XY_ANGS];
  for (int i=0; i<NUM_XY_ANGS; i++)
    xy_angle[i] = -M_PI+(i+0.5)*dxyang;

  //the (real) ang dists themselves
  double W[NUM_STATES][NUM_XY_ANGS][NUM_Z_ANGS];
  double W_coupled[NUM_XY_ANGS][NUM_Z_ANGS];
  for (int i=0; i<NUM_XY_ANGS; i++) {
    for (int j=0; j<NUM_Z_ANGS; j++) {
      for (int s=0; s<NUM_STATES; s++)
	W[s][i][j]=0;
      W_coupled[i][j]=0;
    }
  }
  
  TTree **indiv_tree = new TTree*[NUM_STATES];
  double xy_angle_fill, z_angle_fill, W_fill;
  for (int i=0; i<NUM_STATES; i++) {
    indiv_tree[i] = new TTree(Form("indiv_tree_%i",i),Form("indiv_tree_%i",i));
    indiv_tree[i]->Branch("xy_angle",&xy_angle_fill);
    indiv_tree[i]->Branch("z_angle",&z_angle_fill);
    indiv_tree[i]->Branch("W",&W_fill);
  }
  TTree data_tree("data_tree","data_tree");
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
  TGraph *****re_graphs = new TGraph****[NUM_STATES];
  TCanvas ****re_canvs  = new TCanvas***[NUM_STATES];
  TPad    ****re_pads   = new TPad***[NUM_STATES];

  TGraph *****im_graphs = new TGraph****[NUM_STATES];
  TCanvas ****im_canvs  = new TCanvas***[NUM_STATES];
  TPad    ****im_pads   = new TPad***[NUM_STATES];

  for (int s1=0; s1<NUM_STATES; s1++) {
    
    re_graphs[s1] = new TGraph***[NUM_STATES];
    re_canvs[s1]  = new TCanvas**[NUM_STATES];
    re_pads[s1]   = new TPad**[NUM_STATES];

    im_graphs[s1] = new TGraph***[NUM_STATES];
    im_canvs[s1]  = new TCanvas**[NUM_STATES];
    im_pads[s1]   = new TPad**[NUM_STATES];
    
    for (int s2=0; s2<NUM_STATES; s2++) {

      re_graphs[s1][s2] = new TGraph**[k_max+1];
      re_canvs[s1][s2]  = new TCanvas*[k_max+1];
      re_pads[s1][s2]   = new TPad*[k_max+1];

      im_graphs[s1][s2] = new TGraph**[k_max+1];
      im_canvs[s1][s2]  = new TCanvas*[k_max+1];
      im_pads[s1][s2]   = new TPad*[k_max+1];
      
      for (int k=0; k<=k_max; k++) {

	re_graphs[s1][s2][k] = new TGraph*[2*k_max+1];
	im_graphs[s1][s2][k] = new TGraph*[2*k_max+1];

	re_canvs[s1][s2][k] = new TCanvas(Form("B%i_s%is%i_canv_real",k,s1,s2),Form("B%i_s%is%i_canv_real",k,s1,s2));
	im_canvs[s1][s2][k] = new TCanvas(Form("B%i_s%is%i_canv_imag",k,s1,s2),Form("B%i_s%is%i_canv_imag",k,s1,s2));

	re_pads[s1][s2][k] = new TPad(Form("B%i_s%is%i_pad_real",k,s1,s2),Form("B%i_s%is%i_pad_real",k,s1,s2),0.05,0.05,0.95,0.94);
	im_pads[s1][s2][k] = new TPad(Form("B%i_s%is%i_pad_imag",k,s1,s2),Form("B%i_s%is%i_pad_imag",k,s1,s2),0.05,0.05,0.95,0.94);

	int rows = 0, cols = 0; //num rows and columns in which to divide canvases

	switch (k) { //can get creative w canvas divides (make it look pretty)
	case 0:  rows = 1; cols = 1; break;
	case 1:  rows = 3; cols = 1; break;
	case 2:  rows = 3; cols = 2; break;
	case 3:  rows = 4; cols = 2; break;
	case 4:  rows = 3; cols = 3; break;
	case 5:  rows = 4; cols = 3; break;
	case 6:  rows = 5; cols = 3; break;
	case 7:  rows = 5; cols = 3; break;
	case 8:  rows = 6; cols = 3; break;
	default: rows = k; cols = 3; break;
	}

	re_pads[s1][s2][k]->Divide(cols,rows);
	im_pads[s1][s2][k]->Divide(cols,rows);

	for (int q=-k; q<=k; q++) {

	  re_graphs[s1][s2][k][q+k] = new TGraph(NUM_RXN_ANGS, rxnangle, geebees.Bkq_Real(s1,s2,k,q));
	  re_graphs[s1][s2][k][q+k]->SetTitle(Form("Re(B_{%i%i}(%s,%s));#theta_{rxn}^{CM} (deg.);Amplitude (fm^{2})",k,q,state_label[s1],state_label[s2]));
	  re_pads[s1][s2][k]->cd(q+k+1);
	  re_graphs[s1][s2][k][q+k]->Draw();

	  im_graphs[s1][s2][k][q+k] = new TGraph(NUM_RXN_ANGS, rxnangle, geebees.Bkq_Imag(s1,s2,k,q));
	  im_graphs[s1][s2][k][q+k]->SetTitle(Form("Im(B_{%i%i}(%s,%s));#theta_{rxn}^{CM} (deg.);Amplitude (fm^{2})",k,q,state_label[s1],state_label[s2]));
	  im_pads[s1][s2][k]->cd(q+k+1);
	  im_graphs[s1][s2][k][q+k]->Draw();
    
	} //end q

	re_canvs[s1][s2][k]->cd();
	title_pad.Draw();
	re_pads[s1][s2][k]->Draw();

	im_canvs[s1][s2][k]->cd();
	title_pad.Draw();
	im_pads[s1][s2][k]->Draw();

      } //end k
    }
  }

  //then calculate decay parts

  double indiv_norm[NUM_STATES];
  for (int i=0; i<NUM_STATES; i++)
    indiv_norm[i] = 4*M_PI*geebees.Ak(i,0)*geebees.indiv_bees[i]->Bkq_Real(0,0,rxn_angle_index)*geebees.Spherical_Harmonic(0,0,0,0);
  double coupled_norm = 4*M_PI*real(geebees.akq_Complex(0,0,rxn_angle_index))*geebees.Spherical_Harmonic(0,0,0,0);

  /* 1D HISTOGRAMS*/
  TH1F    **Ak_hist = new TH1F*[NUM_STATES];
  TCanvas **Ak_canv = new TCanvas*[NUM_STATES];
  TPad    **Ak_pad  = new TPad*[NUM_STATES];
  for (int s=0; s<NUM_STATES; s++) {
    Ak_hist[s] = new TH1F(Form("Ak_s%i_hist",s),Form("Coupling Coefficients (A_{k}) | %s",state_label[s]),geebees.indiv_bees[s]->Max_Rank()+1,0,geebees.indiv_bees[s]->Max_Rank()+1);
    for (int k=0; k<=geebees.indiv_bees[s]->Max_Rank(); k++) {
      Ak_hist[s]->Fill(k,geebees.Ak(s,k));
      Ak_hist[s]->GetXaxis()->SetBinLabel(k+1,Form("%i",k));
    }
    Ak_hist[s]->SetStats(0);
    Ak_canv[s] = new TCanvas(Form("Ak_%i_canv",s));
    Ak_pad[s]  = new TPad(Form("Ak_%i_pad",s),Form("Ak_%i_pad",s),0.05,0.05,0.95,0.94);
    Ak_pad[s]->cd();
    Ak_hist[s]->Draw();
    Ak_canv[s]->cd();
    title_pad.Draw();
    Ak_pad[s]->Draw();
  }

  TH1F akq_real_hist("akq_real_hist",Form("Re(a_{kq}) -- #theta_{rxn}^{LAB} = %.1f^{o}",(double)(rxnang_lab*RADTODEG)),(k_max+1)*(k_max+1),0,(k_max+1)*(k_max+1)); akq_real_hist.SetStats(0);
  TH1F akq_imag_hist("akq_imag_hist",Form("Im(a_{kq}) -- #theta_{rxn}^{LAB} = %.1f^{o}",(double)(rxnang_lab*RADTODEG)),(k_max+1)*(k_max+1),0,(k_max+1)*(k_max+1)); akq_imag_hist.SetStats(0);
  TCanvas akq_real_canv("akq_real_canv");
  TPad    akq_real_pad("akq_real_pad","akq_real_pad",0.05,0.05,0.95,0.94);
  TCanvas akq_imag_canv("akq_imag_canv");
  TPad    akq_imag_pad("akq_imag_pad","akq_imag_pad",0.05,0.05,0.95,0.94);
  int bin_iter=0;
  for (int k=0; k<=k_max; k++) {
    for (int q=-k; q<=k; q++) {
      k_fill = k;
      q_fill = q;
      akq_real_fill = geebees.akq_Real(k,q,rxn_angle_index);
      akq_imag_fill = geebees.akq_Imag(k,q,rxn_angle_index);
      akq_real_hist.Fill(bin_iter,akq_real_fill);
      akq_imag_hist.Fill(bin_iter,akq_imag_fill);
      akq_real_hist.GetXaxis()->SetBinLabel(bin_iter+1,Form("%i%i",k,q));
      akq_imag_hist.GetXaxis()->SetBinLabel(bin_iter+1,Form("%i%i",k,q));
      coeff_tree.Fill();
      bin_iter++;
    }
  }
  akq_real_pad.cd();
  akq_real_hist.Draw();
  akq_real_canv.cd();
  title_pad.Draw();
  akq_real_pad.Draw();
  akq_imag_pad.cd();
  akq_imag_hist.Draw();
  akq_imag_canv.cd();
  title_pad.Draw();
  akq_imag_pad.Draw();

  /* 2D HISTOGRAMS */
  TH2F **DecayDist_indiv_real_hist = new TH2F*[NUM_STATES];
  TH2F **DecayDist_indiv_imag_hist = new TH2F*[NUM_STATES];
  for (int s=0; s<NUM_STATES; s++) {
    DecayDist_indiv_real_hist[s] = new TH2F(Form("DecayDist_indiv_%i_real_hist",s),
					    Form("%s",state_label[s]),
					    NUM_XY_ANGS,-180,180,NUM_Z_ANGS,-1,1);
    DecayDist_indiv_real_hist[s]->GetXaxis()->SetTitle(Form("%s_{%s}^{%s} (deg.)",xy_angle_label,D1_label,frame_label));
    DecayDist_indiv_real_hist[s]->GetYaxis()->SetTitle(Form("%s_{%s}^{%s} (unitless)",z_angle_label,D1_label,frame_label));
    DecayDist_indiv_real_hist[s]->SetStats(0);

    DecayDist_indiv_imag_hist[s] = new TH2F(Form("DecayDist_indiv_%i_imag_hist",s),
					    Form("%s",state_label[s]),
					    NUM_XY_ANGS,-180,180,NUM_Z_ANGS,-1,1);
    DecayDist_indiv_imag_hist[s]->GetXaxis()->SetTitle(Form("%s_{%s}^{%s} (deg.)",xy_angle_label,D1_label,frame_label));
    DecayDist_indiv_imag_hist[s]->GetYaxis()->SetTitle(Form("%s_{%s}^{%s} (unitless)",z_angle_label,D1_label,frame_label));
    DecayDist_indiv_imag_hist[s]->SetStats(0);
  }
  
  TH2F DecayDist_coupled_real_hist("DecayDist_coupled_real_hist",
				   Form("Coupled"),
				   NUM_XY_ANGS,-180,180,NUM_Z_ANGS,-1,1);
  DecayDist_coupled_real_hist.GetXaxis()->SetTitle(Form("%s_{%s}^{%s} (deg.)",xy_angle_label,D1_label,frame_label));
  DecayDist_coupled_real_hist.GetYaxis()->SetTitle(Form("%s_{%s}^{%s} (unitless)",z_angle_label,D1_label,frame_label));
  DecayDist_coupled_real_hist.SetStats(0);

  TH2F DecayDist_coupled_imag_hist("DecayDist_coupled_imag_hist",
				   Form("Coupled"),
				   NUM_XY_ANGS,-180,180,NUM_Z_ANGS,-1,1);
  DecayDist_coupled_imag_hist.GetXaxis()->SetTitle(Form("%s_{%s}^{%s} (deg.)",xy_angle_label,D1_label,frame_label));
  DecayDist_coupled_imag_hist.GetYaxis()->SetTitle(Form("%s_{%s}^{%s} (unitless)",z_angle_label,D1_label,frame_label));
  DecayDist_coupled_imag_hist.SetStats(0);

  //variables for total integration
  double Wtot=0;
  
  /* MAIN W LOOPS */

  for (int s=0; s<NUM_STATES; s++) {
    Wtot=0;
    for (int k=0; k<=k_max_indiv[s]; k++) {
      for (int q=-k; q<=k; q++) {
	cout << "\rCalculating state " << s << " W(k,q) = W(" << k << "," << q << ")..." << flush;
	for (int zang=0; zang<NUM_Z_ANGS; zang++) { //decay angle	  
	  for (int xyang=0; xyang<NUM_XY_ANGS; xyang++) { //relative azimuthal angle

	    /***** MAIN W CALCULATION *****/
	    double W_real = real(geebees.Ak(s,k)*
				 geebees.indiv_bees[s]->Bkq_Complex(k,q,rxn_angle_index)*
				 conj(geebees.Complex_Spherical_Harmonic_Cosine(k,q,z_angle_cosine[zang],xy_angle[xyang])))/indiv_norm[s];
	    double W_imag = imag(geebees.Ak(s,k)*
				 geebees.indiv_bees[s]->Bkq_Complex(k,q,rxn_angle_index)*
				 conj(geebees.Complex_Spherical_Harmonic_Cosine(k,q,z_angle_cosine[zang],xy_angle[xyang])))/indiv_norm[s];
	    /*****************************/

	    W[s][xyang][zang] += W_real;
	  
	    //total integration
	    Wtot += W_real*dzang_cos*dxyang;
	  
	    DecayDist_indiv_real_hist[s]->Fill(xy_angle[xyang]*RADTODEG,z_angle_cosine[zang],W_real);
	    DecayDist_indiv_imag_hist[s]->Fill(xy_angle[xyang]*RADTODEG,z_angle_cosine[zang],W_imag);
	    
	  } //end xyang loop
	} //end zang loop
      } //end q loop
    } //end k loop

    for (int zang=0; zang<NUM_Z_ANGS; zang++) {
      for (int xyang=0; xyang<NUM_XY_ANGS; xyang++) {
	xy_angle_fill = xy_angle[xyang];
	z_angle_fill  = z_angle_cosine[zang];
	W_fill        = W[s][xyang][zang];
	indiv_tree[s]->Fill();
      }
    }

    cout << "done (total W = " << Wtot << ")" << endl;
    
  } //end s loop

  Wtot=0;
  for (int k=0; k<=k_max; k++) {
    for (int q=-k; q<=k; q++) {
      cout << "\rCalculating coupled W(k,q) = W(" << k << "," << q << ")..." << flush;
      for (int zang=0; zang<NUM_Z_ANGS; zang++) { //decay angle	  
	for (int xyang=0; xyang<NUM_XY_ANGS; xyang++) { //relative azimuthal angle

	  /***** MAIN W CALCULATION *****/
	  double W_real = real(geebees.akq_Complex(k,q,rxn_angle_index)*
			       conj(geebees.Complex_Spherical_Harmonic_Cosine(k,q,z_angle_cosine[zang],xy_angle[xyang])))/coupled_norm;
	  double W_imag = imag(geebees.akq_Complex(k,q,rxn_angle_index)*
			       conj(geebees.Complex_Spherical_Harmonic_Cosine(k,q,z_angle_cosine[zang],xy_angle[xyang])))/coupled_norm;
	  /*****************************/

	  W_coupled[xyang][zang] += W_real;
	  
	  //total integration
	  Wtot += W_real*dzang_cos*dxyang;
	  
	  DecayDist_coupled_real_hist.Fill(xy_angle[xyang]*RADTODEG,z_angle_cosine[zang],W_real);
	  DecayDist_coupled_imag_hist.Fill(xy_angle[xyang]*RADTODEG,z_angle_cosine[zang],W_imag);
	  
	} //end xyang loop
      } //end zang loop
    } //end q loop
  } //end k loop

  for (int zang=0; zang<NUM_Z_ANGS; zang++) {
    for (int xyang=0; xyang<NUM_XY_ANGS; xyang++) {
      xy_angle_fill = xy_angle[xyang];
      z_angle_fill  = z_angle_cosine[zang];
      W_fill        = W_coupled[xyang][zang];
      data_tree.Fill();
    }
  }

  cout << "done (total W = " << Wtot << ")" << endl
       << endl;
  
  /* DATA VISUALIZATION */

  //ordinary reaction xsec (no decay)
  TGraph  **xsec_indiv_gr   = new TGraph*[NUM_STATES];
  TCanvas **xsec_indiv_canv = new TCanvas*[NUM_STATES];
  TPad    **xsec_indiv_pad  = new TPad*[NUM_STATES];
  for (int s=0; s<NUM_STATES; s++) {
    xsec_indiv_gr[s] = new TGraph(NUM_RXN_ANGS,rxnangle,xsec_indiv[s]);
    xsec_indiv_gr[s]->SetTitle(";#theta_{rxn}^{CM} (deg.);d#sigma/d#Omega (mb/sr)");
    xsec_indiv_gr[s]->SetName(Form("xsec_indiv_%i_gr",s));
    xsec_indiv_gr[s]->SetLineColor(s+2);
    xsec_indiv_canv[s] = new TCanvas(Form("xsec_indiv_%i_canv",s));
    xsec_indiv_pad[s] = new TPad(Form("xsec_indiv_%i_pad",s),Form("xsec_indiv_%i_pad",s),0.05,0.05,0.95,0.94);
    xsec_indiv_pad[s]->cd();
    xsec_indiv_gr[s]->Draw();
    xsec_indiv_pad[s]->SetLogy();
    xsec_indiv_canv[s]->cd();
    title_pad.Draw();
    xsec_indiv_pad[s]->Draw();
  }

  TGraph  xsec_coupled_gr(NUM_RXN_ANGS,rxnangle,xsec_coupled);
  TCanvas xsec_coupled_canv("xsec_coupled_canv");
  TPad    xsec_coupled_pad("xsec_coupled_pad","xsec_coupled_pad",0.05,0.05,0.95,0.94);
  xsec_coupled_gr.SetTitle(";#theta_{rxn}^{CM} (deg.);d#sigma/d#Omega (mb/sr)");
  xsec_coupled_gr.SetName("xsec_coupled_gr");
  xsec_coupled_pad.cd();
  xsec_coupled_gr.Draw();
  xsec_coupled_pad.SetLogy();
  xsec_coupled_canv.cd();
  title_pad.Draw();
  xsec_coupled_pad.Draw();

  TCanvas xsec_all_canv("xsec_all_canv");
  TPad    xsec_all_pad("xsec_all_pad","xsec_all_pad",0.05,0.05,0.95,0.94);
  xsec_all_pad.cd();
  for (int s=0; s<NUM_STATES; s++)
    xsec_indiv_gr[s]->Draw(s==0?"":"same");
  xsec_coupled_gr.Draw("same");
  xsec_all_pad.SetLogy();
  TLegend *xsec_all_leg = xsec_all_pad.BuildLegend();
  xsec_all_leg->Clear();
  for (int s=0; s<NUM_STATES; s++)
    xsec_all_leg->AddEntry(Form("xsec_indiv_%i_gr",s),state_label[s],"l");
  xsec_all_leg->AddEntry("xsec_coupled_gr","Coupled","l");
  xsec_all_leg->SetFillColor(0);
  xsec_all_leg->SetBorderSize(0);
  xsec_all_leg->Draw();
  xsec_all_canv.cd();
  title_pad.Draw();
  xsec_all_pad.Draw();
  
  TCanvas **DecayDist_indiv_real_canv = new TCanvas*[NUM_STATES];
  TPad    **DecayDist_indiv_real_pad = new TPad*[NUM_STATES];
  for (int s=0; s<NUM_STATES; s++) {
    DecayDist_indiv_real_canv[s] = new TCanvas(Form("DecayDist_indiv_%i_real_canv",s));
    DecayDist_indiv_real_pad[s]  = new TPad(Form("DecayDist_indiv_%i_real_pad",s),Form("DecayDist_indiv_%i_real_pad",s),0.05,0.05,0.95,0.94);
    DecayDist_indiv_real_pad[s]->cd();
    DecayDist_indiv_real_hist[s]->Draw("colz");
    DecayDist_indiv_real_canv[s]->cd();
    title_pad.Draw();
    DecayDist_indiv_real_pad[s]->Draw();
  }

  TCanvas **DecayDist_indiv_imag_canv = new TCanvas*[NUM_STATES];
  TPad    **DecayDist_indiv_imag_pad = new TPad*[NUM_STATES];
  for (int s=0; s<NUM_STATES; s++) {
    DecayDist_indiv_imag_canv[s] = new TCanvas(Form("DecayDist_indiv_%i_imag_canv",s));
    DecayDist_indiv_imag_pad[s]  = new TPad(Form("DecayDist_indiv_%i_imag_pad",s),Form("DecayDist_indiv_%i_imag_pad",s),0.05,0.05,0.95,0.94);
    DecayDist_indiv_imag_pad[s]->cd();
    DecayDist_indiv_imag_hist[s]->Draw("colz");
    DecayDist_indiv_imag_canv[s]->cd();
    title_pad.Draw();
    DecayDist_indiv_imag_pad[s]->Draw();
  }

  TCanvas DecayDist_coupled_real_canv("DecayDist_coupled_real_canv");
  TPad    DecayDist_coupled_real_pad("DecayDist_coupled_real_pad","DecayDist_coupled_real_pad",0.05,0.05,0.95,0.94);
  DecayDist_coupled_real_pad.cd();
  DecayDist_coupled_real_hist.Draw("colz");
  DecayDist_coupled_real_canv.cd();
  title_pad.Draw();
  DecayDist_coupled_real_pad.Draw();

  TCanvas DecayDist_coupled_imag_canv("DecayDist_coupled_imag_canv");
  TPad    DecayDist_coupled_imag_pad("DecayDist_coupled_imag_pad","DecayDist_coupled_imag_pad",0.05,0.05,0.95,0.94);
  DecayDist_coupled_imag_pad.cd();
  DecayDist_coupled_imag_hist.Draw("colz");
  DecayDist_coupled_imag_canv.cd();
  title_pad.Draw();
  DecayDist_coupled_imag_pad.Draw();

  int all_cols, all_rows;
  switch (NUM_STATES+1) {
  case 2: all_cols=2; all_rows=1; break;
  case 3: all_cols=3; all_rows=1; break;
  case 4: all_cols=2; all_rows=2; break;
  case 5: all_cols=3; all_rows=2; break;
  case 6: all_cols=3; all_rows=2; break;
  default: all_cols=(NUM_STATES)/2; all_rows=(NUM_STATES)/2; break;
  }
  
  TCanvas DecayDist_all_real_canv("DecayDist_all_real_canv");
  TPad    DecayDist_all_real_pad("DecayDist_all_real_pad","DecayDist_all_real_pad",0.05,0.05,0.95,0.94);
  DecayDist_all_real_pad.Divide(all_cols,all_rows);
  for (int s=0; s<NUM_STATES; s++) {
    DecayDist_all_real_pad.cd(s+1);
    DecayDist_indiv_real_hist[s]->Draw("colz");
  }
  DecayDist_all_real_pad.cd(NUM_STATES+1);
  DecayDist_coupled_real_hist.Draw("colz");
  DecayDist_all_real_canv.cd();
  title_pad.Draw();
  DecayDist_all_real_pad.Draw();

  TCanvas DecayDist_all_imag_canv("DecayDist_all_imag_canv");
  TPad    DecayDist_all_imag_pad("DecayDist_all_imag_pad","DecayDist_all_imag_pad",0.05,0.05,0.95,0.94);
  DecayDist_all_imag_pad.Divide(all_cols,all_rows);
  for (int s=0; s<NUM_STATES; s++) {
    DecayDist_all_imag_pad.cd(s+1);
    DecayDist_indiv_imag_hist[s]->Draw("colz");
  }
  DecayDist_all_imag_pad.cd(NUM_STATES+1);
  DecayDist_coupled_imag_hist.Draw("colz");
  DecayDist_all_imag_canv.cd();
  title_pad.Draw();
  DecayDist_all_imag_pad.Draw();
  
  geebeesfile.cd();

  input_tree.Write();
  for (int s=0; s<NUM_STATES; s++)
    indiv_tree[s]->Write();
  data_tree.Write();
  coeff_tree.Write();

  for (int s=0; s<NUM_STATES; s++) {
    DecayDist_indiv_real_canv[s]->Write();
    DecayDist_indiv_imag_canv[s]->Write();
  }
  DecayDist_coupled_real_canv.Write();
  DecayDist_coupled_imag_canv.Write();
  DecayDist_all_real_canv.Write();
  DecayDist_all_imag_canv.Write();

  for (int s=0; s<NUM_STATES; s++)
    xsec_indiv_canv[s]->Write();
  xsec_coupled_canv.Write();
  xsec_all_canv.Write();

  for (int s=0; s<NUM_STATES; s++)
    Ak_canv[s]->Write();
  akq_real_canv.Write();
  akq_imag_canv.Write();

  for (int s=0; s<NUM_STATES; s++)
    msubs_canv[s]->Write();
    
  for (int s1=0; s1<NUM_STATES; s1++) {
    for (int s2=0; s2<NUM_STATES; s2++) {
      densmat_real_canv[s1][s2]->Write();
      densmat_imag_canv[s1][s2]->Write();
    }
  }
  
  for (int s1=0; s1<NUM_STATES; s1++) {
    for (int s2=0; s2<NUM_STATES; s2++) {
      for (int k=0; k<=k_max; k++) {
	re_canvs[s1][s2][k]->Write();
	im_canvs[s1][s2][k]->Write();
      }
    }
  }
  
  geebeesfile.Close();

  for (int i=0; i<NUM_STATES; i++) {

    delete [] densmat_real_hist[i];
    delete [] densmat_real_canv[i];
    delete [] densmat_real_pad[i];

    delete [] densmat_imag_hist[i];
    delete [] densmat_imag_canv[i];
    delete [] densmat_imag_pad[i];
    
  }
    
  delete [] msubs_hist;
  delete [] msubs_canv;
  delete [] msubs_pad;

  delete [] densmat_real_hist;
  delete [] densmat_real_canv;
  delete [] densmat_real_pad;

  delete [] densmat_imag_hist;
  delete [] densmat_imag_canv;
  delete [] densmat_imag_pad;

  delete [] Ak_hist;
  delete [] Ak_canv;
  delete [] Ak_pad;
  
  delete [] DecayDist_indiv_real_hist;
  delete [] DecayDist_indiv_real_canv;
  delete [] DecayDist_indiv_real_pad;
  delete [] DecayDist_indiv_imag_hist;
  delete [] DecayDist_indiv_imag_canv;
  delete [] DecayDist_indiv_imag_pad;
  
  delete [] indiv_tree;

  for (int s1=0; s1<NUM_STATES; s1++) {
    for (int s2=0; s2<NUM_STATES; s2++) {
      for (int k=0; k<=k_max; k++) {
	delete [] re_graphs[s1][s2][k];
	delete [] im_graphs[s1][s2][k];
      }
      delete [] re_graphs[s1][s2];
      delete [] re_canvs[s1][s2];
      delete [] re_pads[s1][s2];
      delete [] im_graphs[s1][s2];
      delete [] im_canvs[s1][s2];
      delete [] im_pads[s1][s2];
    }
    delete [] re_graphs[s1];
    delete [] re_canvs[s1];
    delete [] re_pads[s1];
    delete [] im_graphs[s1];
    delete [] im_canvs[s1];
    delete [] im_pads[s1];
  }
  delete [] re_graphs;
  delete [] re_canvs;
  delete [] re_pads;
  delete [] im_graphs;
  delete [] im_canvs;
  delete [] im_pads;

  delete [] xsec_indiv;
  
  delete [] Pi;
  delete [] Pi_char;
  delete [] Ji;
  delete [] R_exc;
  delete [] R_res;
  delete [] R_wid;
  delete [] BR;
  for (int i=0; i<NUM_STATES; i++) {
    delete [] fort37_name[i];
    delete [] state_label[i];
  }
  delete [] fort37_name;
  delete [] state_label;
  
  chrono::high_resolution_clock::time_point t_stop = chrono::high_resolution_clock::now();
  chrono::duration<double> t_span = chrono::duration_cast<chrono::duration<double>>(t_stop - t_start);

  cout << "Total time: " << t_span.count() << " seconds (" << t_span.count()/60. << " minutes)\n";

  return 0;

}

void SetNucLabel(char *nuc_label, int Z, int A) {

  if      (Z==0 && A==1) strcpy(nuc_label,"n");
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
