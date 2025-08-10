/* PeakFit()
 *
 * Class to make complex fits of peaks in spectra
 * Currently can provide a fit over a determined range with a preassigned number of peaks
 * in that range. Will return a reduced chi-square value as an initial test of goodness of fit
 * Current peak shapes allowed are gaussians and breit-wigner distributions
 *
 * Gordon M. -- July 2019
 */

#include "PeakFit.h"

using namespace std;

PeakFit::PeakFit() {
  params = NULL;
  nfits = 0;
  file = NULL;
  c1 = new TCanvas();
}

PeakFit::~PeakFit() {
  delete[] params;
  c1->Clear();
  c1->Close();
  if(file != NULL && file->IsOpen()) file->Close();
}

/* Pulls histogram from a file. If the histogram doesnt exist returns false 
 * and program should get terminated in main. Also gets the bin width from the original histogram
 * to be used with integration results (bin width is assumed constant)
 */
bool PeakFit::getHisto(string filename, string histoname) {
  file = new TFile(filename.c_str(), "READ");
  if(file->GetListOfKeys()->Contains(histoname.c_str())) {
    histo = (TH1F*) file->Get(histoname.c_str());
    BIN_WIDTH = histo->GetXaxis()->GetBinWidth(1);
    return true;
  } else {
    cout<<"Error in PeakFit::GetHisto!! Either the file isnt valid or the histogram name "
        <<"does not refer to a valid histogram."<<endl;
    return false;
  }
}

/* Gets the type and range of each peak (and the entire fit) from the user
 * For more complicated functions (not gaus or polN) also need initial guess of params
 */
void PeakFit::getRanges() {
  cout<<"-----------------SPS Focal Plane Spectrum Fitting Tool-----------------"<<endl;
  cout<<"Enter in the range of the histogram to fit over and the number of peaks: "<<endl;
  //histo->Draw();
  //while(c1->WaitPrimitive()) {}; //this means wait for canvas double click
  cout<<"Fit Minimum = ";
  cin>>fullMin;
  cout<<"Fit Maximum = ";
  cin>>fullMax;
  cout<<"Total number of peaks to fit: ";
  cin>>nPeaks;
  nGaussians=0; nBW=0;
  cout<<"-----------------------------------------------------------------------"<<endl;
  cout<<"Enter individual peak information..."<<endl;
  string answer, localname;
  Double_t min_i, max_i;
  Double_t amp, mean, width;
  for(int i=0; i<nPeaks; i++) {
    answer = "";
    cout<<"PEAK "<<i<<endl;
    cout<<"Is this a gaussian or a breit-wigner distribution?(g/b) ";
    cin>>answer;
    if(answer == "g") {
      cout<<"Enter in the range for gaus"<<nGaussians<<": "<<endl;
      //histo->Draw(); 
      //while(c1->WaitPrimitive()) {};
      cout<<"Min = ";
      cin>>min_i;
      cout<<"Max = ";
      cin>>max_i;
      localname = "gaus"+to_string(nGaussians);
      TF1 *f =  new TF1(localname.c_str(), "gaus",min_i, max_i);
      gaussians.push_back(f);
      nGaussians++;
    } else if(answer == "b") {
      cout<<"Enter in range for breit-wigner"<<nBW<<":"<<endl;
      //histo->Draw(); 
      //while(c1->WaitPrimitive()) {};
      cout<<"Min = ";
      cin>>min_i;
      cout<<"Max = ";
      cin>>max_i;
      cout<<"BW requires initial parameter (amplitude, mean, width) guess"<<endl;
      cout<<"Amplitude: ";
      cin>>amp;
      cout<<"Mean: ";
      cin>>mean;
      cout<<"FWHM: ";
      cin>>width;
      localname ="breit-wigner"+to_string(nBW);
      TF1 *f = new TF1(localname.c_str(), "[0]*TMath::BreitWigner(x,[1],[2])",min_i, max_i);
      f->SetParameter(0, amp);
      f->SetParameter(1, mean);
      f->SetParameter(2, width);
      f->SetParLimits(0,0,MAX_AMPLITUDE); //no negative amplitude
      f->SetParLimits(1,min_i, max_i); // mean must be within range specified
      f->SetParLimits(2,0,MAX_WIDTH); //no negative widths
      breitwigners.push_back(f);
      nBW++;
    }
  }
  cout<<"-----------------------------------------------------------------------"<<endl;
  cout<<"Enter in order for polynomial background fit: ";
  cin>>bckgnd_order;
  string bf = "pol"+to_string(bckgnd_order);
  background = new TF1("background",bf.c_str(),fullMin,fullMax);
  cout<<"Polynomial of order "<<bckgnd_order<<" chosen"<<endl;
  cout<<"-----------------------------------------------------------------------"<<endl;
  histo->GetXaxis()->SetRangeUser(fullMin, fullMax);//restrict range for fit
  totalParams = nGaussians*3+nBW*3+(bckgnd_order+1);
  params = new Double_t[totalParams];
  for(int i=bckgnd_order; i>=0; i--) {
    params[(totalParams-1)-i] = 1; //set backgnd params to 1
  }
}

/* Creates all of the TF1's to be used later. Individuals are stored in vectors
 * To make a good estimate for the full function we need both individuals and the 
 * full function
 */
void PeakFit::createFullFunction() {
  string fullname = "";
  for (int i=0; i<nGaussians; i++) {
    fullname += "gaus("+to_string(i*3)+")+";
  }
  for(int i=0; i<nBW; i++) {
    fullname += "["+to_string(nGaussians*3+i*3)+"]*TMath::BreitWigner(x,["
                +to_string(nGaussians*3+i*3+1)+"],["+to_string(nGaussians*3+i*3+2)+"])+";
  }
  fullname += "pol"+to_string(bckgnd_order)+"("+to_string(nGaussians*3+nBW*3)+")";
  multigaus = new TF1("complete_fit",fullname.c_str(),fullMin,fullMax);
  for(int i=0; i<nGaussians; i++) {
    int j = i*3;
    multigaus->SetParLimits(j,0,MAX_AMPLITUDE); //no negative amplitude
    multigaus->SetParLimits(j+2,0,MAX_WIDTH); //no negative widths
  }
  for(int i=0; i<nBW; i++) {
    int j = nGaussians*3+i*3;
    multigaus->SetParLimits(j,0,MAX_AMPLITUDE); //no negative amplitude
    multigaus->SetParLimits(j+2,0,MAX_WIDTH); //no negative widths
  }
  
}

/* Fits the individual functions and then takes the resulting parameters from the 
 * individual fits and gives them to the full fit as initial parameters. Results in a
 * much stronger fit than if initial guesses are used. (This method is particularly strong for
 * gaussians and other ROOT functions since there is no need for any user guessing)
 */
void PeakFit::fitIndividuals() {
  cout<<"---------------------Individual Fit Results----------------------------"<<endl;
  for (int i=0; i<nGaussians; i++) {
    histo->Fit(gaussians[i], "R0+");
    gaussians[i]->GetParameters(&params[i*3]);
    double max = std::max(params[i*3+1]+params[i*3+2],params[i*3+1]-params[i*3+2]);
    double min = std::min(params[i*3+1]+params[i*3+2],params[i*3+1]-params[i*3+2]);
    multigaus->SetParLimits(i*3+1, min, max); //limit mean to avoid travelling peaks
  }
  for (int i=0; i<nBW; i++) {
    histo->Fit(breitwigners[i], "R0B+");
    int bwi = nGaussians*3+i*3;
    breitwigners[i]->GetParameters(&params[bwi]);
    double max = std::max(params[bwi+1]+params[bwi+2],params[bwi+1]-params[bwi+2]);
    double min = std::min(params[bwi+1]+params[bwi+2],params[bwi+1]-params[bwi+2]);
    multigaus->SetParLimits(bwi+1, min, max);
  }
  cout<<"-----------------------------------------------------------------------"<<endl;
}

void PeakFit::fitFull() {
  cout<<"---------------------Full Results FIT"<<to_string(nfits)<<"----------------------------"<<endl;
  multigaus->SetParameters(&params[0]);
  histo->Fit(multigaus, "R0BE+");
  chisq = multigaus->GetChisquare();
  ndf  = multigaus->GetNDF();
  r_chisq = chisq/((Double_t)ndf);
  cout<<"Reduced Chi-square value: "<<r_chisq<<endl;
  cout<<"-----------------------------------------------------------------------"<<endl;
  nfits++;
}

/* Draws fit as both the single global function and the individuals with the parameters from
 * the global fit. User then has the option to either indicate desire to try the fit again or
 * to accept the current fit (to be handled by main)
 */
bool PeakFit::drawFit() {

  c1->cd();
  multigaus->GetParameters(&params[0]);
  histo->Draw();
  multigaus->SetLineColor(kBlue);
  multigaus->Draw("same");
  for(int i = 0; i<nGaussians; i++) {
    TF1 *gaus = gaussians[i];
    gaus->SetParameters(&params[i*3]);
    gaus->SetLineColor(kGreen);
    gaus->DrawF1(fullMax, fullMin, "same");
  }
  for(int i=0; i<nBW; i++) {
    TF1 *bw = breitwigners[i];
    int bwi = nGaussians*3+i*3;
    bw->SetParameters(&params[bwi]);
    bw->SetLineColor(kRed);
    bw->DrawF1(fullMax, fullMin, "same");
  }
  background->SetParameters(&params[nGaussians*3+nBW*3]);
  background->SetLineColor(kBlack);

  background->Draw("same");
  while(c1->WaitPrimitive()) {};
  string answer;
  cout<<"Did you hit STATUS=CALL LIMIT?(y/n) ";
  cin>>answer;
  //indicate needs more calls  or not. the true or false value should then be used by main to 
  //either re-run or finish
  if(answer == "y") {
    return false;
  } else if (answer == "n") {
    return true;
  } else {
    cout<<"That wasn't y or n so guess you're going again"<<endl;
    return false;
  }
}

//Stores the results of the fit in a txt file (chi square, params, and integrated area)
void PeakFit::saveResults(string filename) {
  ofstream outfile(filename);
  if(outfile.is_open()) {
    outfile<<fixed<<showpoint<<setprecision(7);
    outfile<<"Chi-square = "<<chisq<<"\t"<<"DoF = "<<ndf<<"\t"<<"Reduced chi-square = "
           <<r_chisq<<endl;
    outfile<<endl;
    outfile<<setw(10)<<"Gaussian"<<"\t"<<setw(10)<<"Amplitude"<<"\t"<<setw(10)<<"Centroid"
           <<"\t"<<setw(10)<<"Std. Dev."<<"\t"<<"Area"<<endl;
    //Note that integrals are divided by BIN_WIDTH... Numerical integration of the functions
    //is dependent on the bin width while a true integral wouldn't be impacted by the size
    //of dx since it would be assumed infinitesimal
    for(int i=0; i<nGaussians; i++) {
      TF1 *gaus = gaussians[i];
      gaus->SetParameters(&params[i*3]);
      Double_t centroid = gaus->GetParameter(1);
      Double_t amplitude = gaus->GetParameter(0);
      Double_t sigma = fabs(gaus->GetParameter(2));
      Double_t area = gaus->Integral(centroid-3*sigma, centroid+3*sigma)/BIN_WIDTH;
      outfile<<setw(10)<<gaus->GetName()<<"\t"<<amplitude<<"\t"<<centroid<<"\t"<<sigma<<"\t"
             <<area<<endl;
    }
    outfile<<endl;
    outfile<<setw(10)<<"BreitWigner"<<"\t"<<setw(10)<<"Amplitude"<<setw(10)<<"Mean"<<"\t"<<setw(10)<<"FWHM"<<"\t"
           <<setw(10)<<"Area"<<endl;
    for(int i=0; i<nBW; i++) {
      TF1 *bw = breitwigners[i];
      int bwi = nGaussians*3+i*3;
      bw->SetParameters(&params[bwi]);
      Double_t amp = bw->GetParameter(0);
      Double_t mean = bw->GetParameter(1);
      Double_t width = bw->GetParameter(2);
      Double_t area = bw->Integral(mean-width, mean+width)/BIN_WIDTH;
      outfile<<setw(10)<<bw->GetName()<<"\t"<<amp<<"\t"<<mean<<"\t"<<width<<"\t"<<area<<endl;
    }
    outfile<<endl;
    outfile<<setw(10)<<"Background"<<"\t";
    for(int i=0; i<bckgnd_order+1; i++) {
      outfile<<setw(10)<<"a"<<i<<"\t";
    }
    outfile<<endl;
    outfile<<setw(10)<<"pol"<<bckgnd_order<<"\t";
    for(int i=nGaussians*3+nBW*3; i<totalParams; i++) {
      outfile<<setw(10)<<params[i]<<"\t";
    }
    outfile<<endl;
  } else {
    cout<<"Error when writing fit results to file, could not open output file!"<<endl;
  }
  outfile.close();
}

