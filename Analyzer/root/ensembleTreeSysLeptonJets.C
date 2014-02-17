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
  
  ensembles["matchingUp"] = ensemble("Summer12_TTJets1725_matchingup/job_*_ensemble.root", 37083003./1.75);
  ensembles["matchingDown"] = ensemble("Summer12_TTJets1725_matchingdown/job_*_ensemble.root", 34053113./1.75);
  
  ensembles["scaleUp"] = ensemble("Summer12_TTJets1725_scaleup/job_*_ensemble.root", 41908271./1.75);
  ensembles["scaleDown"] = ensemble("Summer12_TTJets1725_scaledown/job_*_ensemble.root", 39286663./1.75);
  
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
  
  ///////////////////////////////////
  
  comparisons["Pile-up"] = comparison("defaultNewWeights", "puUp", "puDown");
  comparisons["JES FlavorPureQuark"] = comparison("defaultNewWeights", "flavorQUp", "flavorQDown", true, false);//
  comparisons["JES FlavorPureGluon"] = comparison("defaultNewWeights", "flavorGUp", "flavorGDown", true, false);//
  comparisons["JES FlavorPureBottom"] = comparison("default", "flavorBUp", "flavorBDown", true, false);
  comparisons["JES FlavorQCD"] = comparison("default", "flavorQCDUp", "flavorQCDDown");
  comparisons["b-fragmentation"] = comparison("defaultNewWeights", "bFragLEP");
  comparisons["Neutrino fraction"] = comparison("defaultNewWeights", "bFNuUp", "bFNuDown");
  comparisons["b-tag"] = comparison("defaultNewWeights", "bTagSFUp", "bTagSFDown");//
  comparisons["Mistag"] = comparison("defaultNewWeights", "misTagSFUp", "misTagSFDown");//
  comparisons["Top-pt"] = comparison("defaultNewWeights", "topPt", "");//
  comparisons["Matching"] = comparison("calibration", "matchingUp", "matchingDown", false);
  comparisons["Scale"] = comparison("calibration", "scaleUp", "scaleDown", false);
  comparisons["JER"] = comparison("default", "jerUp", "jerDown");
  comparisons["JES total"] = comparison("default", "jesUp", "jesDown");
  comparisons["JES PileUpPt"] = comparison("defaultNewWeights", "jesPuUp", "jesPuDown");
  comparisons["MG vs. Powheg"] = comparison("calibration", "powheg", "", false, false);
  comparisons["Powheg+Pythia vs. MCatNLO"] = comparison("powheg", "mcatnlo", "", false, false);
  comparisons["Pythia vs. Herwig"] = comparison("powheg", "powhegHerwig", "", false, false);
  comparisons["Powheg+Herwig vs. MCatNLO"] = comparison("powhegHerwig", "mcatnlo", "", false, false);
  comparisons["MG (SC) vs. Powheg"] = comparison("powheg", "defaultSC", "", false);
  comparisons["Spin correlations"] = comparison("calibration", "defaultSC", "", false, true);
  comparisons["Pythia Z2* vs. P11"] = comparison("defaultSC", "P11", "", false, false);
  comparisons["Color reconnection"] = comparison("P11", "P11noCR", "", false);
  comparisons["Underlying event"] = comparison("P11", "P11mpiHi", "P11TeV", false);
  comparisons["Background"] = comparison("defaultNewWeights", "fSig1.0", "fSig0.9");
  
  
  ///////////////////////////////////
  
  std::cout << "\n### Fitting pseudo-experiments" << std::endl;
  for(std::map<std::string, ensemble>::iterator it = ensembles.begin(); it != ensembles.end(); it++) {
    std::cout << std::setiosflags(std::ios::left) << std::setw(20) << it->first;
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
    }
    
    printf("\tmass = %.2lf+/-%.2lf GeV, jes = %.3lf+/-%.3lf, mass1d = %.2lf+/-%.2lf GeV\n", it->second.mass, it->second.massUnc, it->second.jes, it->second.jesUnc, it->second.mass1d, it->second.mass1dUnc);
  }
  
  ///////////////////////////////////
  
  std::cout << "\n### Systematic uncertainties" << std::endl;
  for(std::map<std::string, comparison>::iterator it = comparisons.begin(); it != comparisons.end(); it++) {
    if (!it->second.active) std::cout << "#";
    std::cout << std::setiosflags(std::ios::left) << std::setw(20) << it->first;
    
    ensemble nominal  = ensembles.find(it->second.nominal)->second;
    ensemble up       = ensembles.find(it->second.up)->second;
    ensemble down;
    if (strcmp(it->second.down, "") == 0) down = up;
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
    
    printf("\tmassShift = %.2lf+/-%.2lf GeV, jesShift = %.3lf+/-%.3lf, mass1dShift = %.2lf+/-%.2lf GeV\n", massShift, massShiftUnc, jesShift, jesShiftUnc, mass1dShift, mass1dShiftUnc);
    
    if (it->second.active) {
      totalMassUncertainty2   += pow(massShift, 2);
      totalJESUncertainty2    += pow(jesShift, 2);
      totalMass1dUncertainty2 += pow(mass1dShift, 2);
    }
  }
  
  std::cout << "\n### Total" << std::endl;
  printf("\tmassShift = %.2lf GeV, jesShift = %.3lf, mass1dShift = %.2lf GeV\n", sqrt(totalMassUncertainty2), sqrt(totalJESUncertainty2), sqrt(totalMass1dUncertainty2));
}


