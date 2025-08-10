/* PeakFit()
 *
 * Class to make complex fits of peaks in spectra
 * Currently can provide a fit over a determined range with a preassigned number of peaks
 * in that range. Will return a reduced chi-square value as an initial test of goodness of fit
 * Current peak shapes allowed are gaussians and breit-wigner distributions
 *
 * Gordon M. -- July 2019
 */

#ifndef PEAKFIT_H
#define PEAKFIT_H

#include <TROOT.h>
#include <TH1.h>
#include <TF1.h>
#include <TFile.h>
#include <TTree.h>
#include <TCanvas.h>
#include <iostream>
#include <string>
#include <vector>
#include <fstream>

using namespace std;

class PeakFit {

  public:
    PeakFit();
    ~PeakFit();
    void createFullFunction();
    void getRanges();
    bool getHisto(string filename, string histoname);
    void saveResults(string filename);
    void fitIndividuals();
    void fitFull();
    bool drawFit();
    
  private:
    vector<TF1*> gaussians, breitwigners;
    TF1* multigaus;
    TF1* background;

    int bckgnd_order;
    int nfits;
    const double MAX_AMPLITUDE = 1e9, MAX_WIDTH = 0.5;

    Double_t BIN_WIDTH;
    Int_t nPeaks, nBW, nGaussians, totalParams;
    Float_t fullMax, fullMin;
    Double_t *params;

    Double_t chisq, r_chisq;
    Int_t ndf;

    TH1F *histo;
    TCanvas *c1;
    TFile *file;
};


#endif
