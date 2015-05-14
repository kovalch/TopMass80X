#include <string>
#include <vector>
#include <iostream>
#include <iomanip>

#include "TFile.h"
#include "TChain.h"
#include "TCanvas.h"
#include "TProfile.h"
#include "TH1F.h"
#include "TF1.h"
#include "TLegend.h"

#include "tdrstyle.C"

struct sample {
  std::string filename;
  bool spread;
  bool trust;
  double offset;
  TFile* file;
  TH1F* wmass;
  TH1F* tmass;
  sample(std::string f, bool s = false, bool t = false, double o = 0.)
  : filename(f), spread(s), trust(t), offset(o) {}
};

std::vector<sample> samples;
std::vector<sample>::iterator it;


void rivetMass() {
  setTDRStyle();
  gStyle->SetPadLeftMargin(0.2);
  gStyle->SetPadRightMargin(0.05);
  gStyle->SetNdivisions(505, "X");
  gStyle->SetTitleYOffset(1.75);
  gStyle->SetOptStat(0);
  gStyle->SetOptFit(0);

  // Samples
  /*
  samples.push_back(sample("TT_8TeV-amcatnlo-5f-herwigpp", true));
  samples.push_back(sample("TT_8TeV-amcatnlo-5f-pythia8", true));
  samples.push_back(sample("TT_8TeV-herwigpp", true));
  samples.push_back(sample("TT_8TeV-pythia8-noME"));
  samples.push_back(sample("TT_8TeV-pythia8", true));
  samples.push_back(sample("TT_Cluster_8TeV-sherpa2", true));
  samples.push_back(sample("TT_CUEP8M1_8TeV-powheg", true));
  samples.push_back(sample("TTJets_MSDecays_CUEP8M1_8TeV-madgraph", true));
  samples.push_back(sample("TTJets_MSDecays_CUEP8M1amc_8TeV-amcatnlo", true));
  samples.push_back(sample("TTJets_MSDecays_CUEP8M1noME_8TeV-madgraph"));
  samples.push_back(sample("TTJets_MSDecays_P11_8TeV-madgraph", true));
  samples.push_back(sample("TTJets_MSDecays_P11TeV_8TeV-madgraph"));
  samples.push_back(sample("TTJets_MSDecays_Z2_8TeV-madgraph", true));
  samples.push_back(sample("TTJets_MSDecays_Z2scaledown_8TeV-madgraph-fix"));
  samples.push_back(sample("TTJets_MSDecays_Z2scaleup_8TeV-madgraph-fix"));
  samples.push_back(sample("TT_Lund_8TeV-sherpa2"));
  samples.push_back(sample("TT_Z2_8TeV-powheg", true));
  */

  /*
  samples.push_back(sample("TTJets_MSDecays_Z2_8TeV-madgraph", true, true));
  samples.push_back(sample("TTJets_MSDecays_Z2_mass169.5_8TeV-madgraph"));
  samples.push_back(sample("TTJets_MSDecays_Z2_mass175.5_8TeV-madgraph"));
  samples.push_back(sample("TTJets_MSDecays_Z2scaledown_8TeV-madgraph-fix", true));
  samples.push_back(sample("TTJets_MSDecays_Z2scaleup_8TeV-madgraph-fix", true));

  samples.push_back(sample("TTJets_MSDecays_P11_8TeV-madgraph", true, true));
  samples.push_back(sample("TTJets_MSDecays_P11TeV_8TeV-madgraph", true));
  samples.push_back(sample("TTJets_MSDecays_P11mpiHi_8TeV-madgraph", true));
  samples.push_back(sample("TTJets_MSDecays_P11noCR_8TeV-madgraph", true));

  samples.push_back(sample("TTJets_MSDecays_CUEP8M1_8TeV-madgraph", true, true));
  samples.push_back(sample("TTJets_MSDecays_CUEP8M1gr_8TeV-madgraph", true));
  samples.push_back(sample("TTJets_MSDecays_CUEP8M1noME_8TeV-madgraph", true));

  samples.push_back(sample("TTJets_CUEP8M1amcME_8TeV-amcatnlo", true, true, 0.5));
  samples.push_back(sample("TTJets_CUEP8M1amc_8TeV-amcatnlo", true, false, 0.5));

  samples.push_back(sample("TT_8TeV-amcatnlo-5f-herwigpp-TOP13007", true));
  samples.push_back(sample("TT_8TeV-amcatnlo-5f-pythia8-TOP13007", true));
  samples.push_back(sample("TT_8TeV-amcatnlo-5f-pythia8me-TOP13007", true, true));

  samples.push_back(sample("TT_8TeV-herwigpp", true));
  samples.push_back(sample("TT_8TeV-pythia8-noME", true));
  samples.push_back(sample("TT_8TeV-pythia8", true));

  samples.push_back(sample("TT_Z2_8TeV-powheg", true, true));
  samples.push_back(sample("TT_CUEP8M1_8TeV-powheg", true, true));
  samples.push_back(sample("TT_EE3C_8TeV-powheg", true));
  samples.push_back(sample("TT_EE5C_8TeV-powheg", true));
  samples.push_back(sample("TT_EE5CnoCR_8TeV-powheg", true));
  samples.push_back(sample("TT_AUET2_8TeV-powheg", true));

  samples.push_back(sample("TT_Lund_8TeV-sherpa2", true));
  samples.push_back(sample("TT_Cluster_8TeV-sherpa2", true));
  //*/

  /*
  samples.push_back(sample("TT_mcAUET2_8TeV-mcatnlo", true));
  samples.push_back(sample("TT_mcHerwigDefault_8TeV-mcatnlo", true));
  samples.push_back(sample("TT_AUET2me_8TeV-powheg", true));
  samples.push_back(sample("TTJets_CUEP8M1amcscaleup_8TeV-amcatnlo", true));
  samples.push_back(sample("TTJets_CUEP8M1amcscaledown_8TeV-amcatnlo", true));
  samples.push_back(sample("TTJets_CUEP8M1amcMEscaleup_8TeV-amcatnlo", true));
  samples.push_back(sample("TTJets_CUEP8M1amcMEscaledown_8TeV-amcatnlo", true));
  */

  //*
  samples.push_back(sample("TT_8TeV-pythia8", true));
  samples.push_back(sample("TT_8TeV-herwigpp", true));
  samples.push_back(sample("TT_8TeV-parton-CUEP8M1noHAD"));
  samples.push_back(sample("TT_8TeV-parton-EE5CnoHAD"));
  samples.push_back(sample("TT_8TeV-parton-CUEP8M1noMPI"));
  samples.push_back(sample("TT_8TeV-parton-EE5CnoMPI"));
  samples.push_back(sample("TT_8TeV-parton-CUEP8M1"));
  samples.push_back(sample("TT_8TeV-parton-EE5C"));
  samples.push_back(sample("TT_8TeV-parton-CUEP8M1noFSR"));
  samples.push_back(sample("TT_8TeV-parton-EE5CnoFSR"));
  samples.push_back(sample("TT_8TeV-parton-CUEP8M1noISR"));
  samples.push_back(sample("TT_8TeV-parton-EE5CnoISR"));
  samples.push_back(sample("TT_8TeV-parton-CUEP8M1nothing"));
  samples.push_back(sample("TT_8TeV-parton-EE5Cnothing"));
  samples.push_back(sample("TT_8TeV-parton-CUEP8M1noME"));
  samples.push_back(sample("TT_8TeV-parton-CUEP8M1scaledown"));
  samples.push_back(sample("TT_8TeV-parton-CUEP8M1scaleup"));
  //samples.push_back(sample("TT_8TeV-parton-EE3C"));
  //*/
  
  samples.push_back(sample("TT_CUEP8M1_8TeV-powheg", true, true));
  samples.push_back(sample("TT_EE3C_8TeV-powheg", true));
  samples.push_back(sample("TT_EE5C_8TeV-powheg", true));
  samples.push_back(sample("TT_CUEP8M1parton_8TeV-powheg"));
  samples.push_back(sample("TT_EE5Cparton_8TeV-powheg"));
  samples.push_back(sample("TT_8TeV-amcatnlo-5f-herwigpp-TOP13007", true));
  samples.push_back(sample("TT_8TeV-amcatnlo-5f-pythia8-TOP13007", true));
  samples.push_back(sample("TT_8TeV-amcatnlo-5f-pythia8me-TOP13007", true, true));
  
  samples.push_back(sample("TT_ILC_0.5TeV-pythia8"));
  samples.push_back(sample("TT_ILC_0.5TeV-herwigpp"));
  samples.push_back(sample("TT_ILC_1TeV-pythia8"));
  //samples.push_back(sample("TT_ILC_2TeV-pythia8"));
  samples.push_back(sample("TT_ILC_1TeV-herwigpp"));
  samples.push_back(sample("TT_2TeV-amcatnlo-pythia8"));
  samples.push_back(sample("TT_2TeV-amcatnlo-herwigpp"));

  TH1F* wspread   = new TH1F("wspread",   "wspread",   100, -10, 10);
  TH1F* tspread   = new TH1F("tspread",   "tspread",   100, -10, 10);
  TH1F* t2dspread = new TH1F("t2dspread", "t2dspread", 100, -10, 10);

  TH1F* wtrust   = new TH1F("wtrust",   "wtrust",   100, -10, 10);
  TH1F* ttrust   = new TH1F("ttrust",   "ttrust",   100, -10, 10);
  TH1F* t2dtrust = new TH1F("t2dtrust", "t2dtrust", 100, -10, 10);

  for (it = samples.begin(); it != samples.end(); ++it) {
    it->file = new TFile((std::string("/nfs/dust/cms/user/mseidel/rivet/")+it->filename+std::string(".root")).c_str());

    //std::cout << it->filename << (it->spread ? " SPREAD " : "") << (it->trust ? " TRUST " : "") << std::endl;

    it->wmass = (TH1F*) it->file->Get("W_mass");
    it->tmass = (TH1F*) it->file->Get("t_mass_W_cut");

    TF1* wfit = new TF1("wfit", "gaus");
    TF1* tfit = new TF1("tfit", "gaus");

    it->wmass->Fit(wfit, "Q0", "", 75, 85);
    it->tmass->Fit(tfit, "Q0", "", 160, 185);

    double wref =  80.4825;
    double tref = 170.9162;

    double wshift   = wfit->GetParameter(1) - wref;
    double tshift   = tfit->GetParameter(1) - tref - it->offset;
    double t2d      = tfit->GetParameter(1) / (wfit->GetParameter(1)/wref);
    double t2dshift = t2d - tref - it->offset;

    if (it->spread) {
      wspread->Fill(wshift); tspread->Fill(tshift); t2dspread->Fill(t2dshift);
    }
    if (it->trust) {
      wtrust->Fill(wshift); ttrust->Fill(tshift); t2dtrust->Fill(t2dshift);
    }

    //printf("W Mean:  %.2lf +/- %.2lf; Sigma: %.2lf +/- %.2lf; Shift: %+.2lf \n", wfit->GetParameter(1), wfit->GetParError(1), wfit->GetParameter(2), wfit->GetParError(2), wshift);
    //printf( "t Mean: %.2lf +/- %.2lf; Sigma: %.2lf +/- %.2lf; Shift: %+.2lf  \n", tfit->GetParameter(1), tfit->GetParError(1), tfit->GetParameter(2), tfit->GetParError(2), tshift);
    //printf(" 2D mt: %.2lf;                                Shift: %+.2lf \n\n", t2d, t2dshift);
    //printf("%s & %+.2lf & %+.2lf & %+.2lf \\\\ \n\n", it->filename.c_str(), wshift, tshift, t2dshift);
    printf("%40s & %.2lf & %.2lf & %.2lf & %.2lf & %.2lf & %.2lf \\\\ \n", it->filename.c_str(), it->tmass->GetMean(), tfit->GetParameter(1), tfit->GetParameter(2), it->wmass->GetMean(), wfit->GetParameter(1), wfit->GetParameter(2));

    it->file->Close();
  }

  printf("  W Spread:  %+.2lf +/- %.2lf \n", wspread->GetMean(), wspread->GetRMS());
  printf("  t Spread:  %+.2lf +/- %.2lf \n", tspread->GetMean(), tspread->GetRMS());
  printf("t2d Spread:  %+.2lf +/- %.2lf \n", t2dspread->GetMean(), t2dspread->GetRMS());

  printf("  W Spread (trust):  %+.2lf +/- %.2lf \n", wtrust->GetMean(), wtrust->GetRMS());
  printf("  t Spread (trust):  %+.2lf +/- %.2lf \n", ttrust->GetMean(), ttrust->GetRMS());
  printf("t2d Spread (trust):  %+.2lf +/- %.2lf \n", t2dtrust->GetMean(), t2dtrust->GetRMS());
}
