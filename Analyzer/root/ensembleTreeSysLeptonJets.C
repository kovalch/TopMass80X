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
  double massShift;
  double massShiftUnc;
  double jesShift;
  double jesShiftUnc;
  double mass1dShift;
  double mass1dShiftUnc;
  comparison(const char* n = "", const char* u = "", const char* d = "", bool c = true, bool a = true)
  : nominal(n), up(u), down(d), correlated(c), active(a) {}
};

struct mergedcomparison {
  const char* comparison1;
  const char* comparison2;
  const char* comparison3;
  const char* comparison4;
  const char* comparison5;
  const char* comparison6;
  mergedcomparison(const char* c1 = "", const char* c2 = "", const char* c3 = "", const char* c4 = "", const char* c5 = "", const char* c6 = "")
  : comparison1(c1), comparison2(c2), comparison3(c3), comparison4(c4), comparison5(c5), comparison6(c6) {}
};

enum lepton           { kElectron, kMuon, kAll, kMuon_BReg};
std::string lepton_ [4] = { "electron", "muon", "lepton", "muon_BReg"};

//int channel = 3;
int channel = 2;

TString globalPath("/nfs/dust/cms/user/mseidel/pseudoexperiments/topmass_131012/");
TString globalPathMS("/nfs/dust/cms/user/mseidel/pseudoexperiments/topmass_140305/");
//TString globalPath("/nfs/dust/test/cms/user/kirschen/BRegression_PE_NewCalibrationJan2014Applied_MCS3250MinEvtsBReg/");

double crossSection   = 230;
double peLumi         = 20000.;

void ensembleTreeSysLeptonJets(TString sPath = globalPath)
{
  std::map<std::string, ensemble> ensembles;
  std::map<std::string, comparison> comparisons;
  std::map<std::string, mergedcomparison> mergedcomparisons;

  sPath += lepton_[channel]; sPath += "/";
  TString sPathMS(globalPathMS);
  sPathMS += lepton_[channel]; sPathMS += "/";
  
  double totalMassUncertainty2 = 0;
  double totalJESUncertainty2  = 0;
  double totalMass1dUncertainty2 = 0;
  
  ///////////////////////////////////
  
  ensembles["calibration"] = ensemble("", 0, 172.5, 1., 172.5);
  ensembles["calibrationNoUnc"] = ensemble("", 0, 172.5, 1., 172.5);
  ensembles["default"] = ensemble("job_*_ensemble.root", 7000000./1.75);
  //ensembles["default"] = ensemble("Summer12_TTJets1725_1.00/job_*_ensemble.root", 7000000./1.75);
  ensembles["defaultNewWeights"] = ensemble("weight.combinedWeight/job_*_ensemble.root", 7000000./1.75); //TODO fsig=1!
  
  ensembles["puUp"] = ensemble("weight.combinedWeight/weight.puWeight*weight.puWeightUp/job_*_ensemble.root", 7000000./1.75);
  ensembles["puDown"] = ensemble("weight.combinedWeight/weight.puWeight*weight.puWeightDown/job_*_ensemble.root", 7000000./1.75);
  
  ensembles["bFragLEP"] = ensemble("weight_frag/job_*_ensemble.root", 7000000./1.75);
  ensembles["bFragLEPHard"] = ensemble("weight_fragHard/job_*_ensemble.root", 7000000./1.75);
  ensembles["bFragLEPSoft"] = ensemble("weight_fragSoft/job_*_ensemble.root", 7000000./1.75);
  
  ensembles["bFNuUp"] = ensemble("weight_fNuUp/job_*_ensemble.root", 7000000./1.75);
  ensembles["bFNuDown"] = ensemble("weight_fNuDown/job_*_ensemble.root", 7000000./1.75);
  
  ensembles["flavorQCDUp"] = ensemble("Summer12_TTJets1725_flavor:up_FlavorQCD/job_*_ensemble.root", 62000000./1.75);
  ensembles["flavorQCDDown"] = ensemble("Summer12_TTJets1725_flavor:down_FlavorQCD/job_*_ensemble.root", 62000000./1.75);
  
  ensembles["MSflavorBUp"] = ensemble("Summer12_TTJetsMS1725_flavor:up_FlavorPureBottom/job_*_ensemble.root", 62000000./1.75);
  ensembles["MSflavorBDown"] = ensemble("Summer12_TTJetsMS1725_flavor:down_FlavorPureBottom/job_*_ensemble.root", 62000000./1.75);
  
  ensembles["MSflavorQUp"] = ensemble("Summer12_TTJetsMS1725_flavor:up_FlavorPureQuark/job_*_ensemble.root", 62000000./1.75);
  ensembles["MSflavorQDown"] = ensemble("Summer12_TTJetsMS1725_flavor:down_FlavorPureQuark/job_*_ensemble.root", 62000000./1.75);
  
  ensembles["MSflavorGUp"] = ensemble("Summer12_TTJetsMS1725_flavor:up_FlavorPureGluon/job_*_ensemble.root", 62000000./1.75);
  ensembles["MSflavorGDown"] = ensemble("Summer12_TTJetsMS1725_flavor:down_FlavorPureGluon/job_*_ensemble.root", 62000000./1.75);
  
  ensembles["MSbTagSFUp"] = ensemble("weight_bTagSFUp_bSFNewRecipe/job_*_ensemble.root", 62000000./1.75);
  ensembles["MSbTagSFDown"] = ensemble("weight_bTagSFDown_bSFNewRecipe/job_*_ensemble.root", 62000000./1.75);
  ensembles["MSmisTagSFUp"] = ensemble("weight_misTagSFUp_bSFNewRecipe/job_*_ensemble.root", 62000000./1.75);
  ensembles["MSmisTagSFDown"] = ensemble("weight_misTagSFDown_bSFNewRecipe/job_*_ensemble.root", 62000000./1.75);
  ensembles["MSbSFNewRecipe"] = ensemble("Summer12_TTJetsMS1725_bSFNewRecipe/job_*_ensemble.root", 62000000./1.75);
  
  ensembles["topPt"] = ensemble("weight_topPt/job_*_ensemble.root", 7000000./1.75);
  
  ensembles["fSig1.0"] = ensemble("fsig_1.00/job_*_ensemble.root", 7000000./1.75);
  ensembles["fSig0.9"] = ensemble("fsig_0.90/job_*_ensemble.root", 7000000./1.75);
  
  //ensembles["MSfSigUp"] = ensemble("Summer12_TTJetsMS1725_1.00_fSig_0.9748/job_*_ensemble.root", 62000000./1.75);
  ensembles["MSfSigUp"] = ensemble("Summer12_TTJetsMS1725_1.00_fSig_1.0000/job_*_ensemble.root", 62000000./1.75);
  ensembles["MSfSigDown"] = ensemble("Summer12_TTJetsMS1725_1.00_fSig_0.9280/job_*_ensemble.root", 62000000./1.75);
  
  ensembles["MSmatchingUp"] = ensemble("Summer12_TTJetsMS1725_matchingup/job_*_ensemble.root", 37000000./1.75);
  ensembles["MSmatchingDown"] = ensemble("Summer12_TTJetsMS1725_matchingdown/job_*_ensemble.root", 34000000./1.75);
  
  ensembles["MSscaleUp"] = ensemble("Summer12_TTJetsMS1725_scaleup/job_*_ensemble.root", 42000000./1.75);
  ensembles["MSscaleDown"] = ensemble("Summer12_TTJetsMS1725_scaledown/job_*_ensemble.root", 39000000./1.75);
  
  ensembles["jerUp"] = ensemble("Summer12_TTJets1725_jer:up/job_*_ensemble.root", 7000000./1.75);
  ensembles["jerDown"] = ensemble("Summer12_TTJets1725_jer:down/job_*_ensemble.root", 7000000./1.75);
  
  ensembles["MSjesUp"] = ensemble("Summer12_TTJetsMS1725_source:up_SubTotalMC/job_*_ensemble.root", 62000000./1.75);
  ensembles["MSjesDown"] = ensemble("Summer12_TTJetsMS1725_source:down_SubTotalMC/job_*_ensemble.root", 62000000./1.75);
  
  ensembles["MSjesFlavorUp"] = ensemble("Summer12_TTJetsMS1725_source:up_CorrelationGroupFlavor/job_*_ensemble.root", 62000000./1.75);
  ensembles["MSjesFlavorDown"] = ensemble("Summer12_TTJetsMS1725_source:down_CorrelationGroupFlavor/job_*_ensemble.root", 62000000./1.75);
  
  ensembles["MSjesIntercalibrationUp"] = ensemble("Summer12_TTJetsMS1725_source:up_CorrelationGroupIntercalibration/job_*_ensemble.root", 62000000./1.75);
  ensembles["MSjesIntercalibrationDown"] = ensemble("Summer12_TTJetsMS1725_source:down_CorrelationGroupIntercalibration/job_*_ensemble.root", 62000000./1.75);
  
  ensembles["MSjesMPFInSituUp"] = ensemble("Summer12_TTJetsMS1725_source:up_CorrelationGroupMPFInSitu/job_*_ensemble.root", 62000000./1.75);
  ensembles["MSjesMPFInSituDown"] = ensemble("Summer12_TTJetsMS1725_source:down_CorrelationGroupMPFInSitu/job_*_ensemble.root", 62000000./1.75);
  
  ensembles["MSjesUncorrelatedWithoutPileUpPtUp"] = ensemble("Summer12_TTJetsMS1725_source:up_CorrelationGroupUncorrelatedWithoutPileUpPt/job_*_ensemble.root", 62000000./1.75);
  ensembles["MSjesUncorrelatedWithoutPileUpPtDown"] = ensemble("Summer12_TTJetsMS1725_source:down_CorrelationGroupUncorrelatedWithoutPileUpPt/job_*_ensemble.root", 62000000./1.75);
  
  ensembles["MSjesPuBBUp"] = ensemble("Summer12_TTJetsMS1725_source:up_PileUpPtBB/job_*_ensemble.root", 62000000./1.75);
  ensembles["MSjesPuBBDown"] = ensemble("Summer12_TTJetsMS1725_source:down_PileUpPtBB/job_*_ensemble.root", 62000000./1.75);
  
  ensembles["MSjesPuECUp"] = ensemble("Summer12_TTJetsMS1725_source:up_PileUpPtEC/job_*_ensemble.root", 62000000./1.75);
  ensembles["MSjesPuECDown"] = ensemble("Summer12_TTJetsMS1725_source:down_PileUpPtEC/job_*_ensemble.root", 62000000./1.75);
  
  ensembles["MSjesPuHFUp"] = ensemble("Summer12_TTJetsMS1725_source:up_PileUpPtHF/job_*_ensemble.root", 62000000./1.75);
  ensembles["MSjesPuHFDown"] = ensemble("Summer12_TTJetsMS1725_source:down_PileUpPtHF/job_*_ensemble.root", 62000000./1.75);
  
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
  
  ensembles["metcl2"] = ensemble("Summer12_TTJets1725_metcl_2_nobg/job_*_ensemble.root", 7000000./1.75);
  
  ensembles["MSdefault"] = ensemble("Summer12_TTJetsMS1725_1.00/job_*_ensemble.root", 62000000./1.75);
  
  ensembles["MSnoLeptonSF"] = ensemble("weight_noLeptonSF/job_*_ensemble.root", 62000000./1.75);
  
  ensembles["PDFcentral"] = ensemble("", 0, 172.57, 1.000, 172.54);
  ensembles["PDFup"]      = ensemble("", 0, 172.66, 1.000, 172.59);
  ensembles["PDFdown"]    = ensemble("", 0, 172.53, 0.999, 172.50);
  
  ensembles["result"] = ensemble("", 0, 171.986, 1.01394, 173.187);
  ensembles["resultPtFix"] = ensemble("", 0, 172.040, 1.00724, 172.664);
  
  ///////////////////////////////////
  
  comparisons["Pile-up (pp cross-section)       "] = comparison("defaultNewWeights", "puUp", "puDown");
  comparisons["Jet energy response (udsc)       "] = comparison("MSdefault", "MSflavorQUp", "MSflavorQDown", true);//
  comparisons["Jet energy response (gluon)      "] = comparison("MSdefault", "MSflavorGUp", "MSflavorGDown", true);//
  comparisons["Jet energy response (b)          "] = comparison("MSdefault", "MSflavorBUp", "MSflavorBDown", true);
  comparisons["Jet energy response (FlavorQCD)  "] = comparison("default", "flavorQCDUp", "flavorQCDDown", true, false);
  comparisons["b fragmentation                  "] = comparison("defaultNewWeights", "bFragLEP");
  comparisons["Semi-leptonic B hadron decays    "] = comparison("defaultNewWeights", "bFNuUp", "bFNuDown");
  comparisons["b-tag rate                       "] = comparison("MSbSFNewRecipe", "MSbTagSFUp", "MSbTagSFDown");//
  comparisons["b-tag (mistag rate)              "] = comparison("MSbSFNewRecipe", "MSmisTagSFUp", "MSmisTagSFDown");//
  comparisons["Top-pt reweighting               "] = comparison("defaultNewWeights", "topPt", "");//
  comparisons["ME-PS matching threshold         "] = comparison("calibration", "MSmatchingUp", "MSmatchingDown", false);
  comparisons["$Q^{2}$ scale                    "] = comparison("calibration", "MSscaleUp", "MSscaleDown", false);
  comparisons["Jet energy resolution            "] = comparison("default", "jerUp", "jerDown");

  comparisons["JEC CorrelationGroupFlavor       "] = comparison("MSdefault", "MSjesFlavorUp", "MSjesFlavorDown");
  comparisons["JEC Intercalibration             "] = comparison("MSdefault", "MSjesIntercalibrationUp", "MSjesIntercalibrationDown");
  comparisons["JEC MPFInSitu                    "] = comparison("MSdefault", "MSjesMPFInSituUp", "MSjesMPFInSituDown");
  comparisons["JEC UncorrelatedWithoutPileUpPt  "] = comparison("MSdefault", "MSjesUncorrelatedWithoutPileUpPtUp", "MSjesUncorrelatedWithoutPileUpPtDown");
  
  comparisons["Pile-up (JEC PileUpPtBB)         "] = comparison("MSdefault", "MSjesPuBBUp", "MSjesPuBBDown");
  comparisons["Pile-up (JEC PileUpPtEC)         "] = comparison("MSdefault", "MSjesPuECUp", "MSjesPuECDown");
  comparisons["Pile-up (JEC PileUpPtHF)         "] = comparison("MSdefault", "MSjesPuHFUp", "MSjesPuHFDown");
  
  comparisons["MadGraph (no SC) vs. Powheg      "] = comparison("calibration", "powheg", "", false, false);
  comparisons["Powheg+Pythia6 vs. MC@NLO+Herwig6"] = comparison("powheg", "mcatnlo", "", false, false);
  comparisons["MG+Pythia6 vs. MC@NLO+Herwig6    "] = comparison("defaultSC", "mcatnlo", "", false, false);
  comparisons["Powheg+Pythia6 vs. Powheg+Herwig6"] = comparison("powheg", "powhegHerwig", "", false, false);
  comparisons["MC@NLO+Herwig6 vs. Powheg+Herwig6"] = comparison("mcatnlo", "powhegHerwig", "", false, false);
  comparisons["ME generator                     "] = comparison("defaultSC", "powheg", "", false);
  comparisons["Spin correlations                "] = comparison("calibration", "defaultSC", "", false, false);
  comparisons["Pythia Z2* vs. P11               "] = comparison("defaultSC", "P11", "", false, false);
  comparisons["Color reconnection               "] = comparison("P11", "P11noCR", "", false);
  comparisons["Underlying event                 "] = comparison("P11", "P11mpiHi", "P11TeV", false);
  comparisons["Non-\\ttbar background            "] = comparison("MSdefault", "MSfSigUp", "MSfSigDown");
  comparisons["Lepton energy scale (electron)   "] = comparison("default", "eesUp", "eesDown");
  comparisons["Lepton energy scale (muon)       "] = comparison("default", "mesUp", "mesDown");
  comparisons["MET                              "] = comparison("defaultNewWeights", "metcl2");
  comparisons["Lepton trigger/id                "] = comparison("MSdefault", "MSnoLeptonSF");
  comparisons["Calibration                      "] = comparison("calibration", "calibrationNoUnc", "", false);
  comparisons["PDF                              "] = comparison("PDFcentral", "PDFup", "PDFdown", true);
  comparisons["JEC PtFix                        "] = comparison("result", "resultPtFix", "", true);
  
  ///////////////////////////////////
  
  mergedcomparisons["JEC"] = mergedcomparison(
    "JEC CorrelationGroupFlavor       ",
    "JEC Intercalibration             ",
    "JEC MPFInSitu                    ",
    "JEC UncorrelatedWithoutPileUpPt  ",
    "JEC PtFix                        "
  );
  mergedcomparisons["JEC Flavor"] = mergedcomparison(
    "Jet energy response (udsc)       ",
    "Jet energy response (gluon)      ",
    "Jet energy response (b)          "
  );
  mergedcomparisons["JEC PileUpPt"] = mergedcomparison(
    "Pile-up (JEC PileUpPtBB)         ",
    "Pile-up (JEC PileUpPtEC)         ",
    "Pile-up (JEC PileUpPtHF)         "
  );
  mergedcomparisons["PileUp"] = mergedcomparison(
    "Pile-up (JEC PileUpPtBB)         ",
    "Pile-up (JEC PileUpPtEC)         ",
    "Pile-up (JEC PileUpPtHF)         ",
    "Pile-up (pp cross-section)       "
  );
  mergedcomparisons["b-tagging"] = mergedcomparison(
    "b-tag rate                       ",
    "b-tag (mistag rate)              "
  );
  
  ///////////////////////////////////
  
  std::cout << "\n### Fitting pseudo-experiments" << std::endl;
  for(std::map<std::string, ensemble>::iterator it = ensembles.begin(); it != ensembles.end(); it++) {
    std::cout << std::setiosflags(std::ios::left) << std::setw(20) << it->first;
    
    int N_PE = 0;
    if (it->first == "calibration") {
      // Calibration uncertainties
      it->second.massUnc = 3.62797e-02 + 6.15495e-02;
      it->second.mass1dUnc = 6.15495e-02;
      it->second.jesUnc  = 3.73474e-04 + 8.43934e-04;
    }
    else if (it->first == "calibrationNoUnc" || it->first == "PDFcentral" || it->first == "PDFup" || it->first == "PDFdown" || it->first == "result" || it->first == "resultPtFix") {
      // Calibration uncertainties
      it->second.massUnc = 0.;
      it->second.mass1dUnc = 0;
      it->second.jesUnc  = 0.;
    }
    else {
      // Get files
      it->second.chain = new TChain("tree");
      if (!it->first.compare(0,2,"MS")) {
        it->second.chain->Add(sPathMS + it->second.file);
      }
      else it->second.chain->Add(sPath + it->second.file);
      
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
    
    it->second.massShift    = max(abs(nominal.mass-up.mass), abs(nominal.mass-down.mass));
    it->second.jesShift     = max(abs(nominal.jes-up.jes), abs(nominal.jes-down.jes));
    if (it->first == "\\pt- and $\\eta$-dependent JES    ") {
      it->second.jesShift = sqrt(pow(it->second.jesShift,2)-pow(0.0138,2));
    }
    it->second.mass1dShift  = max(abs(nominal.mass1d-up.mass1d), abs(nominal.mass1d-down.mass1d));
    
    it->second.massShiftUnc   = 0.;
    it->second.jesShiftUnc    = 0.;
    it->second.mass1dShiftUnc = 0.;
    
    if (!it->second.correlated) {
      it->second.massShiftUnc    = max(sqrt(pow(nominal.massUnc,2)+pow(up.massUnc,2)),sqrt(pow(nominal.massUnc,2)+pow(down.massUnc,2)));
      it->second.jesShiftUnc     = max(sqrt(pow(nominal.jesUnc,2)+pow(up.jesUnc,2)),sqrt(pow(nominal.jesUnc,2)+pow(down.jesUnc,2)));
      it->second.mass1dShiftUnc  = max(sqrt(pow(nominal.mass1dUnc,2)+pow(up.mass1dUnc,2)),sqrt(pow(nominal.mass1dUnc,2)+pow(down.mass1dUnc,2)));
    }
    
    printf(" & %.2lf$\\pm$%.2lf & %.3lf$\\pm$%.3lf & %.2lf$\\pm$%.2lf \\tabularnewline\n", it->second.massShift, it->second.massShiftUnc, it->second.jesShift, it->second.jesShiftUnc, it->second.mass1dShift, it->second.mass1dShiftUnc);
    
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
      totalMassUncertainty2   += pow(max(it->second.massShift,   it->second.massShiftUnc  ), 2);
      totalJESUncertainty2    += pow(max(it->second.jesShift,    it->second.jesShiftUnc   ), 2);
      totalMass1dUncertainty2 += pow(max(it->second.mass1dShift, it->second.mass1dShiftUnc), 2);
    }
  }
  
  std::cout << "\n### Total uncertainties" << std::endl;
  printf("\t 2D mass = %.2lf GeV, JSF = %.3lf, 1D mass = %.2lf GeV\n", sqrt(totalMassUncertainty2), sqrt(totalJESUncertainty2-pow(0.0138, 2)), sqrt(totalMass1dUncertainty2));
  
  for(std::map<std::string, mergedcomparison>::iterator it = mergedcomparisons.begin(); it != mergedcomparisons.end(); it++) {
    std::cout << it->first;
    
    double mergedMassUncertainty2   = 0.;
    double mergedJESUncertainty2    = 0.;
    double mergedMass1dUncertainty2 = 0.;
    
                                      mergedMassUncertainty2 += pow(comparisons.find(it->second.comparison1)->second.massShift, 2);
                                      mergedMassUncertainty2 += pow(comparisons.find(it->second.comparison2)->second.massShift, 2);
    if (it->second.comparison3 != "") mergedMassUncertainty2 += pow(comparisons.find(it->second.comparison3)->second.massShift, 2);
    if (it->second.comparison4 != "") mergedMassUncertainty2 += pow(comparisons.find(it->second.comparison4)->second.massShift, 2);
    if (it->second.comparison5 != "") mergedMassUncertainty2 += pow(comparisons.find(it->second.comparison5)->second.massShift, 2);
    if (it->second.comparison6 != "") mergedMassUncertainty2 += pow(comparisons.find(it->second.comparison6)->second.massShift, 2);
    
                                      mergedJESUncertainty2 += pow(comparisons.find(it->second.comparison1)->second.jesShift, 2);
                                      mergedJESUncertainty2 += pow(comparisons.find(it->second.comparison2)->second.jesShift, 2);
    if (it->second.comparison3 != "") mergedJESUncertainty2 += pow(comparisons.find(it->second.comparison3)->second.jesShift, 2);
    if (it->second.comparison4 != "") mergedJESUncertainty2 += pow(comparisons.find(it->second.comparison4)->second.jesShift, 2);
    if (it->second.comparison5 != "") mergedJESUncertainty2 += pow(comparisons.find(it->second.comparison5)->second.jesShift, 2);
    if (it->second.comparison6 != "") mergedJESUncertainty2 += pow(comparisons.find(it->second.comparison6)->second.jesShift, 2);
    
    if (it->first == "JEC") {
      mergedJESUncertainty2 -= pow(0.0138, 2);
    }
    
                                      mergedMass1dUncertainty2 += pow(comparisons.find(it->second.comparison1)->second.mass1dShift, 2);
                                      mergedMass1dUncertainty2 += pow(comparisons.find(it->second.comparison2)->second.mass1dShift, 2);
    if (it->second.comparison3 != "") mergedMass1dUncertainty2 += pow(comparisons.find(it->second.comparison3)->second.mass1dShift, 2);
    if (it->second.comparison4 != "") mergedMass1dUncertainty2 += pow(comparisons.find(it->second.comparison4)->second.mass1dShift, 2);
    if (it->second.comparison5 != "") mergedMass1dUncertainty2 += pow(comparisons.find(it->second.comparison5)->second.mass1dShift, 2);
    if (it->second.comparison6 != "") mergedMass1dUncertainty2 += pow(comparisons.find(it->second.comparison6)->second.mass1dShift, 2);
    
    printf(" & %.2lf & %.3lf & %.2lf \\tabularnewline\n", sqrt(mergedMassUncertainty2), sqrt(mergedJESUncertainty2), sqrt(mergedMass1dUncertainty2));
  }
}


