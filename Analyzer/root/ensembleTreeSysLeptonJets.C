// Current state: Use largest of sys/stat

#include <vector>
#include <map>
#include <string>
#include <iostream>
#include <iomanip>

#include "TFile.h"
#include "TChain.h"
#include "TF1.h"
#include "TFitResult.h"

#include "tdrstyle.C"

struct ensemble {
  const char* file; // file selection for chain
  double size; // EFFECTIVE sample size
  double mass;
  double jes;
  double mass1d;
  double massUnc;
  double jesUnc;
  double mass1dUnc;
  TChain* chain;
  ensemble(const char* f = "temp", double s = 0, double m = -1, double j = -1, double m1 = -1)
  : file(f), size(s), mass(m), jes(j), mass1d(m1) {}
};

struct comparison {
  const char* nominal;
  const char* up;
  const char* down;
  bool correlated;
  bool active;
  comparison(const char* n = "", const char* u = "", const char* d = "", bool c = true, bool a = true)
  : nominal(n), up(u), down(d), correlated(c), active(a) {}
};

enum lepton           { kElectron, kMuon, kAll, kMuon_BReg};
std::string lepton_ [4] = { "electron", "muon", "lepton", "muon_BReg"};

//int channel = 3;
int channel = 2;

TString globalPath("/nfs/dust/cms/user/mseidel/pseudoexperiments/topmass_131012/");
//TString globalPath("/nfs/dust/test/cms/user/kirschen/BRegression_PE_NewCalibrationJan2014Applied_MCS3250MinEvtsBReg/");

double crossSection   = 230;
double peLumi         = 20000.;

void ensembleTreeSysLeptonJets(TString sPath = globalPath)
{
  std::map<std::string, ensemble> ensembles;
  
  std::map<std::string, comparison> comparisons;

  sPath += lepton_[channel]; sPath += "/";
  
  double totalMassUncertainty2 = 0;
  double totalJESUncertainty2  = 0;
  double totalMass1dUncertainty2 = 0;
  
  ///////////////////////////////////
  
  ensembles["calibration"] = ensemble("", 0, 172.5, 1., 172.5);
  ensembles["default"] = ensemble("job_*_ensemble.root", 7000000./1.75);
  //ensembles["default"] = ensemble("Summer12_TTJets1725_1.00/job_*_ensemble.root", 7000000./1.75);
  ensembles["defaultNewWeights"] = ensemble("weight.combinedWeight/job_*_ensemble.root", 7000000./1.75);
  
  ensembles["puUp"] = ensemble("weight.combinedWeight/weight.puWeight*weight.puWeightUp/job_*_ensemble.root", 7000000./1.75);
  ensembles["puDown"] = ensemble("weight.combinedWeight/weight.puWeight*weight.puWeightDown/job_*_ensemble.root", 7000000./1.75);
  
  ensembles["bFragLEP"] = ensemble("weight_frag/job_*_ensemble.root", 7000000./1.75);
  ensembles["bFragLEPHard"] = ensemble("weight_fragHard/job_*_ensemble.root", 7000000./1.75);
  ensembles["bFragLEPSoft"] = ensemble("weight_fragSoft/job_*_ensemble.root", 7000000./1.75);
  
  ensembles["bFNuUp"] = ensemble("weight_fNuUp/job_*_ensemble.root", 7000000./1.75);
  ensembles["bFNuDown"] = ensemble("weight_fNuDown/job_*_ensemble.root", 7000000./1.75);
  
  ensembles["flavorQCDUp"] = ensemble("Summer12_TTJets1725_flavor:up_FlavorQCD/job_*_ensemble.root", 7000000./1.75);
  ensembles["flavorQCDDown"] = ensemble("Summer12_TTJets1725_flavor:down_FlavorQCD/job_*_ensemble.root", 7000000./1.75);
  
  ensembles["flavorBUp"] = ensemble("Summer12_TTJets1725_flavor:up_FlavorPureBottom/job_*_ensemble.root", 7000000./1.75);
  ensembles["flavorBDown"] = ensemble("Summer12_TTJets1725_flavor:down_FlavorPureBottom/job_*_ensemble.root", 7000000./1.75);
  
  ensembles["flavorQUp"] = ensemble("Summer12_TTJets1725_flavor:up_FlavorPureQuark/job_*_ensemble.root", 7000000./1.75);
  ensembles["flavorQDown"] = ensemble("Summer12_TTJets1725_flavor:down_FlavorPureQuark/job_*_ensemble.root", 7000000./1.75);
  
  ensembles["flavorGUp"] = ensemble("Summer12_TTJets1725_flavor:up_FlavorPureGluon/job_*_ensemble.root", 7000000./1.75);
  ensembles["flavorGDown"] = ensemble("Summer12_TTJets1725_flavor:down_FlavorPureGluon/job_*_ensemble.root", 7000000./1.75);
  
  ensembles["bTagSFUp"] = ensemble("weight_bTagSFUp/job_*_ensemble.root", 7000000./1.75);
  ensembles["bTagSFDown"] = ensemble("weight_bTagSFDown/job_*_ensemble.root", 7000000./1.75);
  ensembles["misTagSFUp"] = ensemble("weight_misTagSFUp/job_*_ensemble.root", 7000000./1.75);
  ensembles["misTagSFDown"] = ensemble("weight_misTagSFDown/job_*_ensemble.root", 7000000./1.75);
  
  ensembles["topPt"] = ensemble("weight_topPt/job_*_ensemble.root", 7000000./1.75);
  
  ensembles["fSig1.0"] = ensemble("fsig_1.00/job_*_ensemble.root", 7000000./1.75);
  ensembles["fSig0.9"] = ensemble("fsig_0.90/job_*_ensemble.root", 7000000./1.75);
  
  ensembles["matchingUp"] = ensemble("Summer12_TTJets1725_matchingup/job_*_ensemble.root", 5415010./1.75);
  ensembles["matchingDown"] = ensemble("Summer12_TTJets1725_matchingdown/job_*_ensemble.root", 5476728./1.75);
  
  ensembles["scaleUp"] = ensemble("Summer12_TTJets1725_scaleup/job_*_ensemble.root", 5009488./1.75);
  ensembles["scaleDown"] = ensemble("Summer12_TTJets1725_scaledown/job_*_ensemble.root", 5387181./1.75);
  
  ensembles["jerUp"] = ensemble("Summer12_TTJets1725_jer:up/job_*_ensemble.root", 7000000./1.75);
  ensembles["jerDown"] = ensemble("Summer12_TTJets1725_jer:down/job_*_ensemble.root", 7000000./1.75);
  
  ensembles["jesUp"] = ensemble("Summer12_TTJets1725_source:up_Total/job_*_ensemble.root", 7000000./1.75);
  ensembles["jesDown"] = ensemble("Summer12_TTJets1725_source:down_Total/job_*_ensemble.root", 7000000./1.75);
  
  ensembles["jesPuUp"] = ensemble("Summer12_TTJets1725_source:up_PileUpPtBB/job_*_ensemble.root", 7000000./1.75);
  ensembles["jesPuDown"] = ensemble("Summer12_TTJets1725_source:down_PileUpPtBB/job_*_ensemble.root", 7000000./1.75);
  
  ensembles["mcatnlo"] = ensemble("Summer12_TTJets1725_mcatnlo_herwig/job_*_ensemble.root", 32852589./1.75);
  ensembles["powheg"] = ensemble("Summer12_TTJets1725_powheg/job_*_ensemble.root", 21675970./1.75);
  ensembles["powhegHerwig"] = ensemble("Summer12_TTJets1725_powheg_herwig/job_*_ensemble.root", 27684235./1.75);
  
  ensembles["defaultSC"] = ensemble("Summer12_TTJets1725_MGDecays/job_*_ensemble.root", 56000000./1.75);
  ensembles["P11"] = ensemble("Summer12_TTJets1725_MGDecays_P11/job_*_ensemble.root", 27000000./1.75);
  ensembles["P11noCR"] = ensemble("Summer12_TTJets1725_MGDecays_P11noCR/job_*_ensemble.root", 27000000./1.75);
  ensembles["P11mpiHi"] = ensemble("Summer12_TTJets1725_MGDecays_P11mpiHi/job_*_ensemble.root", 18000000./1.75);
  ensembles["P11TeV"] = ensemble("Summer12_TTJets1725_MGDecays_P11TeV/job_*_ensemble.root", 18000000./1.75);
  
  ensembles["eesUp"] = ensemble("Summer12_TTJets1725_eesShift_+1/job_*_ensemble.root", 7000000./1.75);
  ensembles["eesDown"] = ensemble("Summer12_TTJets1725_eesShift_-1/job_*_ensemble.root", 7000000./1.75);
  
  ensembles["mesUp"] = ensemble("Summer12_TTJets1725_mesShift_+1/job_*_ensemble.root", 7000000./1.75);
  ensembles["mesDown"] = ensemble("Summer12_TTJets1725_mesShift_-1/job_*_ensemble.root", 7000000./1.75);
  
  ensembles["metcl2"] = ensemble("Summer12_TTJets1725_metcl_2/job_*_ensemble.root", 7000000./1.75);
  
  ///////////////////////////////////
  
  comparisons["Pile-up (pp cross-section)       "] = comparison("defaultNewWeights", "puUp", "puDown");
  comparisons["Jet energy response (udsc)       "] = comparison("defaultNewWeights", "flavorQUp", "flavorQDown", true);//
  comparisons["Jet energy response (gluon)      "] = comparison("defaultNewWeights", "flavorGUp", "flavorGDown", true);//
  comparisons["Jet energy response (b)          "] = comparison("default", "flavorBUp", "flavorBDown", true);
  comparisons["Jet energy response (FlavorQCD)  "] = comparison("default", "flavorQCDUp", "flavorQCDDown", true, false);
  comparisons["b fragmentation                  "] = comparison("defaultNewWeights", "bFragLEP");
  comparisons["Semi-leptonic B hadron decays    "] = comparison("defaultNewWeights", "bFNuUp", "bFNuDown");
  comparisons["b-tag rate                       "] = comparison("defaultNewWeights", "bTagSFUp", "bTagSFDown");//
  comparisons["b-tag (mistag rate)              "] = comparison("defaultNewWeights", "misTagSFUp", "misTagSFDown");//
  comparisons["Top-pt reweighting               "] = comparison("defaultNewWeights", "topPt", "");//
  comparisons["ME-PS matching threshold         "] = comparison("calibration", "matchingUp", "matchingDown", false);
  comparisons["$Q^{2}$ scale                    "] = comparison("calibration", "scaleUp", "scaleDown", false);
  comparisons["Jet energy resolution            "] = comparison("default", "jerUp", "jerDown");
  comparisons["\\pt- and $\\eta$-dependent JES    "] = comparison("default", "jesUp", "jesDown");
  comparisons["Pile-up (JES)                    "] = comparison("defaultNewWeights", "jesPuUp", "jesPuDown");
  comparisons["MadGraph (no SC) vs. Powheg      "] = comparison("calibration", "powheg", "", false, false);
  comparisons["Powheg+Pythia6 vs. MC@NLO+Herwig6"] = comparison("powheg", "mcatnlo", "", false, false);
  comparisons["Powheg+Pythia6 vs. Powheg+Herwig6"] = comparison("powheg", "powhegHerwig", "", false, false);
  comparisons["MC@NLO+Herwig6 vs. Powheg+Herwig6"] = comparison("mcatnlo", "powhegHerwig", "", false, false);
  comparisons["ME generator                     "] = comparison("defaultSC", "powheg", "", false);
  comparisons["Spin correlations                "] = comparison("calibration", "defaultSC", "", false, true);
  comparisons["Pythia Z2* vs. P11               "] = comparison("defaultSC", "P11", "", false, false);
  comparisons["Color reconnection               "] = comparison("P11", "P11noCR", "", false);
  comparisons["Underlying event                 "] = comparison("P11", "P11mpiHi", "P11TeV", false);
  comparisons["Non-\\ttbar background            "] = comparison("defaultNewWeights", "fSig0.9", "fSig1.0");
  comparisons["Lepton energy scale (electron)   "] = comparison("defaultNewWeights", "eesUp", "eesDown");
  comparisons["Lepton energy scale (muon)       "] = comparison("defaultNewWeights", "mesUp", "mesDown");
  comparisons["MET                              "] = comparison("defaultNewWeights", "metcl2");
  
  
  ///////////////////////////////////
  
  std::cout << "\n### Fitting pseudo-experiments" << std::endl;
  for(std::map<std::string, ensemble>::iterator it = ensembles.begin(); it != ensembles.end(); it++) {
    std::cout << std::setiosflags(std::ios::left) << std::setw(20) << it->first;
    
    int N_PE = 0;
    if (strcmp(it->second.file, "") == 0) {
      // Calibration uncertainties
      it->second.massUnc = 0.02;
      it->second.jesUnc  = 0.;
    }
    else {
      // Get files
      it->second.chain = new TChain("tree");
      it->second.chain->Add(sPath + it->second.file);
      
      // Fit
      TF1* gaus = new TF1("gaus", "gaus");
      
      it->second.chain->Fit("gaus", "mass_mTop_JES", "mass_mTop_JES>0 & JES_mTop_JES>0 & genMass==172.5 & genJES==1", "LEMQ0");
      it->second.mass    = gaus->GetParameter(1);
      it->second.massUnc = gaus->GetParameter(2) / sqrt(it->second.size/(crossSection*peLumi));
      
      it->second.chain->Fit("gaus", "JES_mTop_JES", "mass_mTop_JES>0 & JES_mTop_JES>0 & genMass==172.5 & genJES==1", "LEMQ0");
      it->second.jes    = gaus->GetParameter(1);
      it->second.jesUnc = gaus->GetParameter(2) / sqrt(it->second.size/(crossSection*peLumi));
      
      it->second.chain->Fit("gaus", "mass_mTop", "mass_mTop>0 & genMass==172.5 & genJES==1", "LEMQ0");
      it->second.mass1d    = gaus->GetParameter(1);
      it->second.mass1dUnc = gaus->GetParameter(2) / sqrt(it->second.size/(crossSection*peLumi));
      
      N_PE = it->second.chain->GetEntries("mass_mTop_JES>0 & JES_mTop_JES>0 & genMass==172.5 & genJES==1");
    }
    
    printf("\tmass = %.2lf+/-%.2lf GeV, jes = %.3lf+/-%.3lf, mass1d = %.2lf+/-%.2lf GeV, N(PE) = %i\n", it->second.mass, it->second.massUnc, it->second.jes, it->second.jesUnc, it->second.mass1d, it->second.mass1dUnc, N_PE);
  }
  
  ///////////////////////////////////
  
  std::cout << "\n### Systematic uncertainties" << std::endl;
  std::cout << "Uncertainty name & 2D mass & JSF & 1D mass \\tabularnewline\n\\hline" << std::endl;
  for(std::map<std::string, comparison>::iterator it = comparisons.begin(); it != comparisons.end(); it++) {
    std::cout << "\\hline" <<std::endl;
    if (!it->second.active) std::cout << "(cc) ";
    std::cout << it->first;
    
    ensemble nominal  = ensembles.find(it->second.nominal)->second;
    ensemble up       = ensembles.find(it->second.up)->second;
    ensemble down;
    bool upDown = true;
    if (strcmp(it->second.down, "") == 0) {
      down   = up;
      upDown = false;
    }
    else down         = ensembles.find(it->second.down)->second;
    
    double massShift    = max(abs(nominal.mass-up.mass), abs(nominal.mass-down.mass));
    double jesShift     = max(abs(nominal.jes-up.jes), abs(nominal.jes-down.jes));
    double mass1dShift  = max(abs(nominal.mass1d-up.mass1d), abs(nominal.mass1d-down.mass1d));
    
    double massShiftUnc   = 0.;
    double jesShiftUnc    = 0.;
    double mass1dShiftUnc = 0.;
    
    if (!it->second.correlated) {
      massShiftUnc    = max(sqrt(pow(nominal.massUnc,2)+pow(up.massUnc,2)),sqrt(pow(nominal.massUnc,2)+pow(down.massUnc,2)));
      jesShiftUnc     = max(sqrt(pow(nominal.jesUnc,2)+pow(up.jesUnc,2)),sqrt(pow(nominal.jesUnc,2)+pow(down.jesUnc,2)));
      mass1dShiftUnc  = max(sqrt(pow(nominal.mass1dUnc,2)+pow(up.mass1dUnc,2)),sqrt(pow(nominal.mass1dUnc,2)+pow(down.mass1dUnc,2)));
    }
    
    printf(" & %.2lf$\\pm$%.2lf & %.3lf$\\pm$%.3lf & %.2lf$\\pm$%.2lf \\tabularnewline\n", massShift, massShiftUnc, jesShift, jesShiftUnc, mass1dShift, mass1dShiftUnc);
    
    if (upDown) {
      if (!it->second.active) std::cout << "(cc) ";
      printf("- up                              & %+.2lf & %+.3lf & %+.2lf \\tabularnewline\n", up.mass-nominal.mass, up.jes-nominal.jes, up.mass1d-nominal.mass1d);
      if (!it->second.active) std::cout << "(cc) ";
      printf("- down                            & %+.2lf & %+.3lf & %+.2lf \\tabularnewline\n", down.mass-nominal.mass, down.jes-nominal.jes, down.mass1d-nominal.mass1d);
    }
    else {
      if (!it->second.active) std::cout << "(cc) ";
      printf("- shift                           & %+.2lf & %+.3lf & %+.2lf \\tabularnewline\n", up.mass-nominal.mass, up.jes-nominal.jes, up.mass1d-nominal.mass1d);
    }
    
    if (it->second.active) {
      totalMassUncertainty2   += pow(max(massShift,   massShiftUnc  ), 2);
      totalJESUncertainty2    += pow(max(jesShift,    jesShiftUnc   ), 2);
      totalMass1dUncertainty2 += pow(max(mass1dShift, mass1dShiftUnc), 2);
    }
  }
  
  std::cout << "\n### Total uncertainties" << std::endl;
  printf("\t 2D mass = %.2lf GeV, JSF = %.3lf, 1D mass = %.2lf GeV\n", sqrt(totalMassUncertainty2), sqrt(totalJESUncertainty2), sqrt(totalMass1dUncertainty2));
}


