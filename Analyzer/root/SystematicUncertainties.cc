#include "SystematicUncertainties.h"

#include <fstream>
#include <iostream>
#include <iomanip>
#include <cmath>

#include "TF1.h"

SystematicUncertainties::SystematicUncertainties()
{
  deriveSystematics();
}


void SystematicUncertainties::fillLeptonJets()
{
  enum lepton              { kElectron, kMuon, kAll, kMuon_BReg };
  std::string lepton_[4] = {"electron", "muon", "lepton", "muon_BReg"};

  //int channel = 3;
  int channel = 2;

  int selectComparisonSet = 0; //0: full set of uncertainties for Moriond; 1: only uncertainties also determined for B-Regression study; 2: minimal set used for BReg-study

  std::string path;
  if(channel==3)path = "/nfs/dust/test/cms/user/kirschen/BRegression_PE_NewCalibrationJan2014Applied_MCS3250MinEvtsBReg/";
  //else path = "/nfs/dust/cms/user/kovalch/GRID-CONTROL_JOBS/pseudoexperiments/topmass_141020/";
  else path = "/nfs/dust/cms/user/mseidel/pseudoexperiments/topmass_paper/";

  path += lepton_[channel]; path += "/";

  sample.path = path;
  sample.crossSection = 230;
  sample.peLumi = 20000.;

  sample.variables = {"mass_mTop_JES", "JES_mTop_JES", "mass_mTop"};

  sample.ensembles["ref"] = ensemble("", 0, {std::make_pair("mass_mTop_JES",std::make_pair(172.5,0.0)), std::make_pair("JES_mTop_JES",std::make_pair(1.0,0.0)), std::make_pair("mass_mTop",std::make_pair(172.5,0.0))});
  sample.ensembles["calibration"] = ensemble("", 0, {std::make_pair("mass_mTop_JES",std::make_pair(172.5,0.04)), std::make_pair("JES_mTop_JES",std::make_pair(1.0,0.0005)), std::make_pair("mass_mTop",std::make_pair(172.5,0.04))});
  sample.ensembles["stat"] = ensemble("", 0, {std::make_pair("mass_mTop_JES",std::make_pair(172.5,0.194498)), std::make_pair("JES_mTop_JES",std::make_pair(1.0,0.00180928)), std::make_pair("mass_mTop",std::make_pair(172.5,0.118054))});
  
  sample.ensembles["nobkg_pdfup"] = ensemble("", 0, {std::make_pair("mass_mTop_JES",std::make_pair(172.613,0.)), std::make_pair("JES_mTop_JES",std::make_pair(1.0007,0.)), std::make_pair("mass_mTop",std::make_pair(172.589,0.))});
  sample.ensembles["nobkg_pdfdown"] = ensemble("", 0, {std::make_pair("mass_mTop_JES",std::make_pair(172.486,0.)), std::make_pair("JES_mTop_JES",std::make_pair(0.9992,0.)), std::make_pair("mass_mTop",std::make_pair(172.501,0.))});
  
  sample.ensembles["default"] = ensemble("Summer12_TTJetsMS1725_1.00/job_*_ensemble.root", 62131965./1.2);

  sample.ensembles["puUp"] = ensemble("weight_puUp/job_*_ensemble.root", 62131965./1.2);
  sample.ensembles["puDown"] = ensemble("weight_puDown/job_*_ensemble.root", 62131965./1.2);

  sample.ensembles["bFragLEP"] = ensemble("weight_frag/job_*_ensemble.root", 62131965./1.2);
  sample.ensembles["bFragLEPHard"] = ensemble("weight_fragHard/job_*_ensemble.root", 62131965./1.2);
  sample.ensembles["bFragLEPSoft"] = ensemble("weight_fragSoft/job_*_ensemble.root", 62131965./1.2);

  sample.ensembles["bFNuUp"] = ensemble("weight_fNuUp/job_*_ensemble.root", 62131965./1.2);
  sample.ensembles["bFNuDown"] = ensemble("weight_fNuDown/job_*_ensemble.root", 62131965./1.2);

  sample.ensembles["flavorCUp"] = ensemble("Summer12_TTJetsMS1725_flavor:up_FlavorPureCharm/job_*_ensemble.root", 62131965./1.2);
  sample.ensembles["flavorCDown"] = ensemble("Summer12_TTJetsMS1725_flavor:down_FlavorPureCharm/job_*_ensemble.root", 62131965./1.2);

  sample.ensembles["flavorBUp"] = ensemble("Summer12_TTJetsMS1725_flavor:up_FlavorPureBottom/job_*_ensemble.root", 62131965./1.2);
  sample.ensembles["flavorBDown"] = ensemble("Summer12_TTJetsMS1725_flavor:down_FlavorPureBottom/job_*_ensemble.root", 62131965./1.2);

  sample.ensembles["flavorQUp"] = ensemble("Summer12_TTJetsMS1725_flavor:up_FlavorPureQuark/job_*_ensemble.root", 62131965./1.2);
  sample.ensembles["flavorQDown"] = ensemble("Summer12_TTJetsMS1725_flavor:down_FlavorPureQuark/job_*_ensemble.root", 62131965./1.2);

  sample.ensembles["flavorGUp"] = ensemble("Summer12_TTJetsMS1725_flavor:up_FlavorPureGluon/job_*_ensemble.root", 62131965./1.2);
  sample.ensembles["flavorGDown"] = ensemble("Summer12_TTJetsMS1725_flavor:down_FlavorPureGluon/job_*_ensemble.root", 62131965./1.2);

  sample.ensembles["bTagSFUp"] = ensemble("weight_bTagSFUp/job_*_ensemble.root", 62131965./1.2);
  sample.ensembles["bTagSFDown"] = ensemble("weight_bTagSFDown/job_*_ensemble.root", 62131965./1.2);
  sample.ensembles["misTagSFUp"] = ensemble("weight_misTagSFUp/job_*_ensemble.root", 62131965./1.2);
  sample.ensembles["misTagSFDown"] = ensemble("weight_misTagSFDown/job_*_ensemble.root", 62131965./1.2);
  sample.ensembles["bTagCSVUp"] = ensemble("Summer12_TTJetsMS1725_csvm_0.73/job_*_ensemble.root", 62131965./1.2);
  sample.ensembles["bTagCSVDown"] = ensemble("Summer12_TTJetsMS1725_csvm_0.63/job_*_ensemble.root", 62131965./1.2);

  sample.ensembles["nobkg"] = ensemble("nobkg/job_*_ensemble.root", 62131965./1.2);
  sample.ensembles["nobkg_topPt"] = ensemble("nobkg_weight_topPt/job_*_ensemble.root", 62131965./1.2);

  sample.ensembles["DibosonDown"] = ensemble("DibosonDown/job_*_ensemble.root", 62131965./1.2);
  sample.ensembles["DibosonUp"] = ensemble("DibosonUp/job_*_ensemble.root", 62131965./1.2);
  sample.ensembles["QCDDown"] = ensemble("QCDDown/job_*_ensemble.root", 62131965./1.2);
  sample.ensembles["QCDUp"] = ensemble("QCDUp/job_*_ensemble.root", 62131965./1.2);
  sample.ensembles["WJetsDown"] = ensemble("WJetsDown/job_*_ensemble.root", 62131965./1.2);
  sample.ensembles["WJetsUp"] = ensemble("WJetsUp/job_*_ensemble.root", 62131965./1.2);
  sample.ensembles["ZJetsDown"] = ensemble("ZJetsDown/job_*_ensemble.root", 62131965./1.2);
  sample.ensembles["ZJetsUp"] = ensemble("ZJetsUp/job_*_ensemble.root", 62131965./1.2);
  sample.ensembles["singleTopDown"] = ensemble("singleTopDown/job_*_ensemble.root", 62131965./1.2);
  sample.ensembles["singleTopUp"] = ensemble("singleTopUp/job_*_ensemble.root", 62131965./1.2);
  sample.ensembles["singleTopScaleDown"] = ensemble("singleTopScaleDown/job_*_ensemble.root", 62131965./1.2);
  sample.ensembles["singleTopScaleUp"] = ensemble("singleTopScaleUp/job_*_ensemble.root", 62131965./1.2);

  sample.ensembles["matchingUp"] = ensemble("Summer12_TTJetsMS1725_matchingup/job_*_ensemble.root", 37000000./1.2);//!
  sample.ensembles["matchingDown"] = ensemble("Summer12_TTJetsMS1725_matchingdown/job_*_ensemble.root", 34000000./1.2);//!

  sample.ensembles["scaleUp"] = ensemble("Summer12_TTJetsMS1725_scaleup/job_*_ensemble.root", 42000000./1.2);//!
  sample.ensembles["scaleDown"] = ensemble("Summer12_TTJetsMS1725_scaledown/job_*_ensemble.root", 39000000./1.2);//!

  sample.ensembles["jerUp"] = ensemble("Summer12_TTJetsMS1725_jer:1/job_*_ensemble.root", 62131965./1.2);
  sample.ensembles["jerDown"] = ensemble("Summer12_TTJetsMS1725_jer:-1/job_*_ensemble.root", 62131965./1.2);

  sample.ensembles["jesUp"] = ensemble("Summer12_TTJetsMS1725_source:up_Total/job_*_ensemble.root", 62131965./1.2);
  sample.ensembles["jesDown"] = ensemble("Summer12_TTJetsMS1725_source:down_Total/job_*_ensemble.root", 62131965./1.2);

  sample.ensembles["mcatnlo"] = ensemble("Summer12_TTJets1725_mcatnlo_herwig/job_*_ensemble.root", 32852589./1.2);
  sample.ensembles["powheg"] = ensemble("Summer12_TTJets1725_powheg/job_*_ensemble.root", 21675970./1.2);
  sample.ensembles["powhegHerwig"] = ensemble("Summer12_TTJets1725_powheg_herwig/job_*_ensemble.root", 27684235./1.2);

  sample.ensembles["defaultSC"] = ensemble("Summer12_TTJets1725_MGDecays_Z2/job_*_ensemble.root", 56000000./1.2);
  sample.ensembles["P11"] = ensemble("Summer12_TTJets1725_MGDecays_P11/job_*_ensemble.root", 27000000./1.2);//!
  sample.ensembles["P11noCR"] = ensemble("Summer12_TTJets1725_MGDecays_P11noCR/job_*_ensemble.root", 27000000./1.2);//!
  sample.ensembles["P11mpiHi"] = ensemble("Summer12_TTJets1725_MGDecays_P11mpiHi/job_*_ensemble.root", 18000000./1.2);//!
  sample.ensembles["P11TeV"] = ensemble("Summer12_TTJets1725_MGDecays_P11TeV/job_*_ensemble.root", 18000000./1.2);//!

  sample.ensembles["eesUp"] = ensemble("Summer12_TTJetsMS1725_eesShift_+1/job_*_ensemble.root", 62131965./1.2);
  sample.ensembles["eesDown"] = ensemble("Summer12_TTJetsMS1725_eesShift_-1/job_*_ensemble.root", 62131965./1.2);

  sample.ensembles["mesUp"] = ensemble("Summer12_TTJetsMS1725_mesShift_+1/job_*_ensemble.root", 62131965./1.2);
  sample.ensembles["mesDown"] = ensemble("Summer12_TTJetsMS1725_mesShift_-1/job_*_ensemble.root", 62131965./1.2);

  sample.ensembles["metUp"] = ensemble("Summer12_TTJetsMS1725_uncFactor_1.1/job_*_ensemble.root", 62131965./1.2);
  sample.ensembles["metDown"] = ensemble("Summer12_TTJetsMS1725_uncFactor_0.9/job_*_ensemble.root", 62131965./1.2);
  /////////////////////////////////////
  sample.ensembles["parj81Up"] = ensemble("Summer12_TTJets1725_FSIM_parj81_up/job_*_ensemble.root", 64000000./1.2);
  sample.ensembles["parj81Nominal"] = ensemble("Summer12_TTJets1725_FSIM_parj81_nominal/job_*_ensemble.root", 64000000./1.2);
  sample.ensembles["parj81Down"] = ensemble("Summer12_TTJets1725_FSIM_parj81_down/job_*_ensemble.root", 64000000./1.2);

  sample.ensembles["sherpa"] = ensemble("Summer12_TTJets1725_sherpa/job_*_ensemble.root", 19514018./2.*9./1.2);

  sample.ensembles["defaultSC_Z2_RD"] = ensemble("Summer12_TTJets1725_MGDecays_Z2_RD/job_*_ensemble.root", 24949110./2.*9./1.2);

  sample.ensembles["AbsoluteFlavMapUp"] = ensemble("Summer12_TTJetsMS1725_source:up_AbsoluteFlavMap/job_*_ensemble.root", 62131965./1.2);
  sample.ensembles["AbsoluteFlavMapDown"] = ensemble("Summer12_TTJetsMS1725_source:down_AbsoluteFlavMap/job_*_ensemble.root", 62131965./1.2);
  sample.ensembles["AbsoluteMPFBiasUp"] = ensemble("Summer12_TTJetsMS1725_source:up_AbsoluteMPFBias/job_*_ensemble.root", 62131965./1.2);
  sample.ensembles["AbsoluteMPFBiasDown"] = ensemble("Summer12_TTJetsMS1725_source:down_AbsoluteMPFBias/job_*_ensemble.root", 62131965./1.2);
  sample.ensembles["AbsoluteScaleUp"] = ensemble("Summer12_TTJetsMS1725_source:up_AbsoluteScale/job_*_ensemble.root", 62131965./1.2);
  sample.ensembles["AbsoluteScaleDown"] = ensemble("Summer12_TTJetsMS1725_source:down_AbsoluteScale/job_*_ensemble.root", 62131965./1.2);
  sample.ensembles["AbsoluteStatUp"] = ensemble("Summer12_TTJetsMS1725_source:up_AbsoluteStat/job_*_ensemble.root", 62131965./1.2);
  sample.ensembles["AbsoluteStatDown"] = ensemble("Summer12_TTJetsMS1725_source:down_AbsoluteStat/job_*_ensemble.root", 62131965./1.2);
  sample.ensembles["HighPtExtraUp"] = ensemble("Summer12_TTJetsMS1725_source:up_HighPtExtra/job_*_ensemble.root", 62131965./1.2);
  sample.ensembles["HighPtExtraDown"] = ensemble("Summer12_TTJetsMS1725_source:down_HighPtExtra/job_*_ensemble.root", 62131965./1.2);
  sample.ensembles["PileUpDataMCUp"] = ensemble("Summer12_TTJetsMS1725_source:up_PileUpDataMC/job_*_ensemble.root", 62131965./1.2);
  sample.ensembles["PileUpDataMCDown"] = ensemble("Summer12_TTJetsMS1725_source:down_PileUpDataMC/job_*_ensemble.root", 62131965./1.2);
  sample.ensembles["PileUpPtBBUp"] = ensemble("Summer12_TTJetsMS1725_source:up_PileUpPtBB/job_*_ensemble.root", 62131965./1.2);
  sample.ensembles["PileUpPtBBDown"] = ensemble("Summer12_TTJetsMS1725_source:down_PileUpPtBB/job_*_ensemble.root", 62131965./1.2);
  sample.ensembles["PileUpPtEC1Up"] = ensemble("Summer12_TTJetsMS1725_source:up_PileUpPtEC1/job_*_ensemble.root", 62131965./1.2);
  sample.ensembles["PileUpPtEC1Down"] = ensemble("Summer12_TTJetsMS1725_source:down_PileUpPtEC1/job_*_ensemble.root", 62131965./1.2);
  sample.ensembles["PileUpPtEC2Up"] = ensemble("Summer12_TTJetsMS1725_source:up_PileUpPtEC2/job_*_ensemble.root", 62131965./1.2);
  sample.ensembles["PileUpPtEC2Down"] = ensemble("Summer12_TTJetsMS1725_source:down_PileUpPtEC2/job_*_ensemble.root", 62131965./1.2);
  sample.ensembles["PileUpPtHFUp"] = ensemble("Summer12_TTJetsMS1725_source:up_PileUpPtHF/job_*_ensemble.root", 62131965./1.2);
  sample.ensembles["PileUpPtHFDown"] = ensemble("Summer12_TTJetsMS1725_source:down_PileUpPtHF/job_*_ensemble.root", 62131965./1.2);
  sample.ensembles["PileUpPtRefUp"] = ensemble("Summer12_TTJetsMS1725_source:up_PileUpPtRef/job_*_ensemble.root", 62131965./1.2);
  sample.ensembles["PileUpPtRefDown"] = ensemble("Summer12_TTJetsMS1725_source:down_PileUpPtRef/job_*_ensemble.root", 62131965./1.2);
  sample.ensembles["RelativeFSRUp"] = ensemble("Summer12_TTJetsMS1725_source:up_RelativeFSR/job_*_ensemble.root", 62131965./1.2);
  sample.ensembles["RelativeFSRDown"] = ensemble("Summer12_TTJetsMS1725_source:down_RelativeFSR/job_*_ensemble.root", 62131965./1.2);
  sample.ensembles["RelativeJEREC1Up"] = ensemble("Summer12_TTJetsMS1725_source:up_RelativeJEREC1/job_*_ensemble.root", 62131965./1.2);
  sample.ensembles["RelativeJEREC1Down"] = ensemble("Summer12_TTJetsMS1725_source:down_RelativeJEREC1/job_*_ensemble.root", 62131965./1.2);
  sample.ensembles["RelativeJEREC2Up"] = ensemble("Summer12_TTJetsMS1725_source:up_RelativeJEREC2/job_*_ensemble.root", 62131965./1.2);
  sample.ensembles["RelativeJEREC2Down"] = ensemble("Summer12_TTJetsMS1725_source:down_RelativeJEREC2/job_*_ensemble.root", 62131965./1.2);
  sample.ensembles["RelativeJERHFUp"] = ensemble("Summer12_TTJetsMS1725_source:up_RelativeJERHF/job_*_ensemble.root", 62131965./1.2);
  sample.ensembles["RelativeJERHFDown"] = ensemble("Summer12_TTJetsMS1725_source:down_RelativeJERHF/job_*_ensemble.root", 62131965./1.2);
  sample.ensembles["RelativePtBBUp"] = ensemble("Summer12_TTJetsMS1725_source:up_RelativePtBB/job_*_ensemble.root", 62131965./1.2);
  sample.ensembles["RelativePtBBDown"] = ensemble("Summer12_TTJetsMS1725_source:down_RelativePtBB/job_*_ensemble.root", 62131965./1.2);
  sample.ensembles["RelativePtEC1Up"] = ensemble("Summer12_TTJetsMS1725_source:up_RelativePtEC1/job_*_ensemble.root", 62131965./1.2);
  sample.ensembles["RelativePtEC1Down"] = ensemble("Summer12_TTJetsMS1725_source:down_RelativePtEC1/job_*_ensemble.root", 62131965./1.2);
  sample.ensembles["RelativePtEC2Up"] = ensemble("Summer12_TTJetsMS1725_source:up_RelativePtEC2/job_*_ensemble.root", 62131965./1.2);
  sample.ensembles["RelativePtEC2Down"] = ensemble("Summer12_TTJetsMS1725_source:down_RelativePtEC2/job_*_ensemble.root", 62131965./1.2);
  sample.ensembles["RelativePtHFUp"] = ensemble("Summer12_TTJetsMS1725_source:up_RelativePtHF/job_*_ensemble.root", 62131965./1.2);
  sample.ensembles["RelativePtHFDown"] = ensemble("Summer12_TTJetsMS1725_source:down_RelativePtHF/job_*_ensemble.root", 62131965./1.2);
  sample.ensembles["RelativeStatEC2Up"] = ensemble("Summer12_TTJetsMS1725_source:up_RelativeStatEC2/job_*_ensemble.root", 62131965./1.2);
  sample.ensembles["RelativeStatEC2Down"] = ensemble("Summer12_TTJetsMS1725_source:down_RelativeStatEC2/job_*_ensemble.root", 62131965./1.2);
  sample.ensembles["RelativeStatHFUp"] = ensemble("Summer12_TTJetsMS1725_source:up_RelativeStatHF/job_*_ensemble.root", 62131965./1.2);
  sample.ensembles["RelativeStatHFDown"] = ensemble("Summer12_TTJetsMS1725_source:down_RelativeStatHF/job_*_ensemble.root", 62131965./1.2);
  sample.ensembles["SinglePionECALUp"] = ensemble("Summer12_TTJetsMS1725_source:up_SinglePionECAL/job_*_ensemble.root", 62131965./1.2);
  sample.ensembles["SinglePionECALDown"] = ensemble("Summer12_TTJetsMS1725_source:down_SinglePionECAL/job_*_ensemble.root", 62131965./1.2);
  sample.ensembles["SinglePionHCALUp"] = ensemble("Summer12_TTJetsMS1725_source:up_SinglePionHCAL/job_*_ensemble.root", 62131965./1.2);
  sample.ensembles["SinglePionHCALDown"] = ensemble("Summer12_TTJetsMS1725_source:down_SinglePionHCAL/job_*_ensemble.root", 62131965./1.2);

  ///////////////////////////////////


  switch( selectComparisonSet )
    {
    case 0: // uncertainties not (yet) determined or considered for b-regression study
      sample.comparisons["Jet energy response (b)          "] = comparison("default", "flavorBUp", "flavorBDown", true); // not meaningful for b-regression study
      sample.comparisons["\\pt- and $\\eta$-dependent JES  "] = comparison("default", "jesUp", "jesDown", true, false); // no change expected for b-regression study
      sample.comparisons["Jet energy resolution            "] = comparison("default", "jerUp", "jerDown"); // no change expected for b-regression study
      
      sample.comparisons["Background: Diboson              "] = comparison("default", "DibosonUp", "DibosonDown");
      sample.comparisons["Background: QCD                  "] = comparison("default", "QCDUp", "QCDDown");
      sample.comparisons["Background: WJets                "] = comparison("default", "WJetsUp", "WJetsDown");
      sample.comparisons["Background: ZJets                "] = comparison("default", "ZJetsUp", "ZJetsDown");
      sample.comparisons["Background: singleTop            "] = comparison("default", "singleTopUp", "singleTopDown");
      sample.comparisons["Background: singleTop scale      "] = comparison("default", "singleTopScaleUp", "singleTopScaleDown");
      
      sample.comparisons["Lepton energy scale (electron)   "] = comparison("default", "eesUp", "eesDown");
      sample.comparisons["Lepton energy scale (muon)       "] = comparison("default", "mesUp", "mesDown");
      sample.comparisons["MET                              "] = comparison("default", "metUp", "metDown");
      
    case 1: // uncertainties determined for both variants
      sample.comparisons["Calibration                      "] = comparison("ref", "calibration", "", false);
      sample.comparisons["Stat                             "] = comparison("ref", "stat", "", false, false);
      sample.comparisons["Pile-up (pp cross-section)       "] = comparison("default", "puUp", "puDown");
      sample.comparisons["Jet energy response (uds)        "] = comparison("default", "flavorQUp", "flavorQDown", true);
      sample.comparisons["Jet energy response (c)          "] = comparison("default", "flavorCUp", "flavorCDown", true);
      sample.comparisons["Jet energy response (gluon)      "] = comparison("default", "flavorGUp", "flavorGDown", true);
      sample.comparisons["b fragmentation                  "] = comparison("default", "bFragLEP");
      sample.comparisons["Semi-leptonic B hadron decays    "] = comparison("default", "bFNuUp", "bFNuDown");
      sample.comparisons["b-tag rate                       "] = comparison("default", "bTagSFUp", "bTagSFDown");
      sample.comparisons["b-tag (mistag rate)              "] = comparison("default", "misTagSFUp", "misTagSFDown");
      sample.comparisons["b-tag disc                       "] = comparison("default", "bTagCSVUp", "bTagCSVDown", true, false);
      sample.comparisons["Top-pt reweighting               "] = comparison("nobkg", "nobkg_topPt", "");
      sample.comparisons["PDF                              "] = comparison("nobkg", "nobkg_pdfup", "nobkg_pdfdown");
      sample.comparisons["ME-PS matching threshold         "] = comparison("calibration", "matchingUp", "matchingDown", false);
      sample.comparisons["$Q^{2}$ scale                    "] = comparison("calibration", "scaleUp", "scaleDown", false);
      sample.comparisons["Powheg+Pythia6 vs. MC@NLO+Herwig6"] = comparison("powheg", "mcatnlo", "", false, false);
      sample.comparisons["Powheg+Pythia6 vs. Powheg+Herwig6"] = comparison("powheg", "powhegHerwig", "", false, false);
      sample.comparisons["MC@NLO+Herwig6 vs. Powheg+Herwig6"] = comparison("mcatnlo", "powhegHerwig", "", false, false);
      sample.comparisons["ME generator                     "] = comparison("calibration", "powheg", "", false);
      sample.comparisons["Pythia Z2* vs. P11               "] = comparison("calibration", "P11", "", false, false);
      sample.comparisons["Color reconnection               "] = comparison("P11", "P11noCR", "", false);
      sample.comparisons["Underlying event                 "] = comparison("P11", "P11mpiHi", "P11TeV", false);


      ///////
      sample.comparisons["Run dependency                   "] = comparison("defaultSC", "defaultSC_Z2_RD", "", true, false);
      sample.comparisons["FSR scale                        "] = comparison("parj81Nominal", "parj81Up", "parj81Down", false);
      sample.comparisons["sherpa                           "] = comparison("calibration", "sherpa", "", false, false);
      sample.comparisons["AbsoluteFlavMap                  "] = comparison("default", "AbsoluteFlavMapUp", "AbsoluteFlavMapDown");
      sample.comparisons["AbsoluteMPFBias                  "] = comparison("default", "AbsoluteMPFBiasUp", "AbsoluteMPFBiasDown");
      sample.comparisons["AbsoluteScale                    "] = comparison("default", "AbsoluteScaleUp", "AbsoluteScaleDown");
      sample.comparisons["AbsoluteStat                     "] = comparison("default", "AbsoluteStatUp", "AbsoluteStatDown");
      sample.comparisons["HighPtExtra                      "] = comparison("default", "HighPtExtraUp", "HighPtExtraDown");
      sample.comparisons["PileUpDataMC                     "] = comparison("default", "PileUpDataMCUp", "PileUpDataMCDown");
      sample.comparisons["PileUpPtBB                       "] = comparison("default", "PileUpPtBBUp", "PileUpPtBBDown");
      sample.comparisons["PileUpPtEC1                      "] = comparison("default", "PileUpPtEC1Up", "PileUpPtEC1Down");
      sample.comparisons["PileUpPtEC2                      "] = comparison("default", "PileUpPtEC2Up", "PileUpPtEC2Down");
      sample.comparisons["PileUpPtHF                       "] = comparison("default", "PileUpPtHFUp", "PileUpPtHFDown");
      sample.comparisons["PileUpPtRef                      "] = comparison("default", "PileUpPtRefUp", "PileUpPtRefDown");
      sample.comparisons["RelativeFSR                      "] = comparison("default", "RelativeFSRUp", "RelativeFSRDown");
      sample.comparisons["RelativeJEREC1                   "] = comparison("default", "RelativeJEREC1Up", "RelativeJEREC1Down");
      sample.comparisons["RelativeJEREC2                   "] = comparison("default", "RelativeJEREC2Up", "RelativeJEREC2Down");
      sample.comparisons["RelativeJERHF                    "] = comparison("default", "RelativeJERHFUp", "RelativeJERHFDown");
      sample.comparisons["RelativePtBB                     "] = comparison("default", "RelativePtBBUp", "RelativePtBBDown");
      sample.comparisons["RelativePtEC1                    "] = comparison("default", "RelativePtEC1Up", "RelativePtEC1Down");
      sample.comparisons["RelativePtEC2                    "] = comparison("default", "RelativePtEC2Up", "RelativePtEC2Down");
      sample.comparisons["RelativePtHF                     "] = comparison("default", "RelativePtHFUp", "RelativePtHFDown");
      sample.comparisons["RelativeStatEC2                  "] = comparison("default", "RelativeStatEC2Up", "RelativeStatEC2Down");
      sample.comparisons["RelativeStatHF                   "] = comparison("default", "RelativeStatHFUp", "RelativeStatHFDown");
      sample.comparisons["SinglePionECAL                   "] = comparison("default", "SinglePionECALUp", "SinglePionECALDown");
      sample.comparisons["SinglePionHCAL                   "] = comparison("default", "SinglePionHCALUp", "SinglePionHCALDown");


      break;
      
    case 3: // minimal set of comparison to determine differences with/ without b-regression
      sample.comparisons["Powheg+Pythia6 vs. Powheg+Herwig6"] = comparison("powheg", "powhegHerwig", "", false);
      //sample.comparisons["b fragmentation                  "] = comparison("default", "bFragLEP");
      //sample.comparisons["Semi-leptonic B hadron decays    "] = comparison("default", "bFNuUp", "bFNuDown");
      break;
    }

    sample.mergedcomparisons["JEC pt/eta                      "] = mergedcomparison({
      "AbsoluteFlavMap                  ",
      "AbsoluteMPFBias                  ",
      "AbsoluteScale                    ",
      "AbsoluteStat                     ",
      "HighPtExtra                      ",
      "PileUpDataMC                     ",
      "PileUpPtBB                       ",
      "PileUpPtEC1                      ",
      "PileUpPtEC2                      ",
      "PileUpPtHF                       ",
      "PileUpPtRef                      ",
      "RelativeFSR                      ",
      "RelativeJEREC1                   ",
      "RelativeJEREC2                   ",
      "RelativeJERHF                    ",
      "RelativePtBB                     ",
      "RelativePtEC1                    ",
      "RelativePtEC2                    ",
      "RelativePtHF                     ",
      "RelativeStatEC2                  ",
      "RelativeStatHF                   ",
      "SinglePionECAL                   ",
      "SinglePionHCAL                   "
    });
    
    sample.mergedcomparisons["CorrelationGroupMPFInSitu       "] = mergedcomparison({
      "AbsoluteMPFBias                  ",
    });
    
    sample.mergedcomparisons["CorrelationGroupInterCalibration"] = mergedcomparison({
      "RelativeFSR                      ",
    });
    
    sample.mergedcomparisons["CorrelationGroupUncorrelated    "] = mergedcomparison({
      "AbsoluteFlavMap                  ",
      "AbsoluteScale                    ",
      "AbsoluteStat                     ",
      "HighPtExtra                      ",
      "PileUpDataMC                     ",
      "PileUpPtBB                       ",
      "PileUpPtEC1                      ",
      "PileUpPtEC2                      ",
      "PileUpPtHF                       ",
      "PileUpPtRef                      ",
      "RelativeJEREC1                   ",
      "RelativeJEREC2                   ",
      "RelativeJERHF                    ",
      "RelativePtBB                     ",
      "RelativePtEC1                    ",
      "RelativePtEC2                    ",
      "RelativePtHF                     ",
      "RelativeStatEC2                  ",
      "RelativeStatHF                   ",
      "SinglePionECAL                   ",
      "SinglePionHCAL                   "
    });
    
    sample.mergedcomparisons["CorrelationGroupUncorrelatedwoPU"] = mergedcomparison({
      "AbsoluteFlavMap                  ",
      "AbsoluteScale                    ",
      "AbsoluteStat                     ",
      "HighPtExtra                      ",
      "RelativeJEREC1                   ",
      "RelativeJEREC2                   ",
      "RelativeJERHF                    ",
      "RelativePtBB                     ",
      "RelativePtEC1                    ",
      "RelativePtEC2                    ",
      "RelativePtHF                     ",
      "RelativeStatEC2                  ",
      "RelativeStatHF                   ",
      "SinglePionECAL                   ",
      "SinglePionHCAL                   "
    });
    
    sample.mergedcomparisons["CorrelationGroupUncorrelatedPU  "] = mergedcomparison({
      "PileUpDataMC                     ",
      "PileUpPtBB                       ",
      "PileUpPtEC1                      ",
      "PileUpPtEC2                      ",
      "PileUpPtHF                       ",
      "PileUpPtRef                      ",
    });
    
    sample.mergedcomparisons["CorrelationGroupFlavor          "] = mergedcomparison({
      "Jet energy response (uds)        ",
      "Jet energy response (c)          ",
      "Jet energy response (b)          ",
      "Jet energy response (gluon)      "
    }, true);
    
    sample.mergedcomparisons["CorrelationGroupbJES            "] = mergedcomparison({
      "b fragmentation                  ",
      "Semi-leptonic B hadron decays    "
    });
    
    sample.mergedcomparisons["Background                      "] = mergedcomparison({
      "Background: Diboson              ",
      "Background: QCD                  ",
      "Background: WJets                ",
      "Background: ZJets                ",
      "Background: singleTop            ",
      "Background: singleTop scale      "
    });
    
    sample.mergedcomparisons["B tagging                       "] = mergedcomparison({
      "b-tag rate                       ",
      "b-tag (mistag rate)              "
    });
    
    sample.mergedcomparisons["Lepton energy scale             "] = mergedcomparison({
      "Lepton energy scale (electron)   ",
      "Lepton energy scale (muon)       "
    });

}

void SystematicUncertainties::fillAllJets()
{
  //sample.path = "/nfs/dust/cms/user/eschliec/TopMass/topmass_140317_1201/";
  //sample.path = "/nfs/dust/cms/user/eschliec/TopMass/topmass_140401_1201/";
  //sample.path = "/nfs/dust/cms/user/eschliec/TopMass/topmass_140418_1201a/";
  //sample.path = "/nfs/dust/cms/user/eschliec/TopMass/topmass_140520_1801/";
  //sample.path = "/nfs/dust/cms/user/eschliec/TopMass/topmass_141217_1601/";
  //sample.path = "/nfs/dust/cms/user/eschliec/TopMass/topmass_150112_1601/";
  //sample.path = "/nfs/dust/cms/user/eschliec/TopMass/topmass_150114_1701/"; // old
  //sample.path = "/nfs/dust/cms/user/eschliec/TopMass/topmass_150120_1401/"; // older
  //sample.path = "/nfs/dust/cms/user/eschliec/TopMass/topmass_150122_1601/"; // oldTEST
  //sample.path = "/nfs/dust/cms/user/eschliec/TopMass/topmass_150121_1301/"; // new
  //sample.path = "/nfs/dust/cms/user/eschliec/TopMass/topmass_150121_1701/"; // newer
  //sample.path = "/nfs/dust/cms/user/eschliec/TopMass/topmass_150202_1401/"; // olderTEST
  sample.path = "/nfs/dust/cms/user/eschliec/TopMass/topmass_150211_1801/"; // older DEFAULT WITHOUT JSF-CONSTRAINT FOR PAPER
  sample.crossSection = 245.794;
  sample.peLumi = 18192.;

  double strafFaktor = 1.15;

  //sample.variables = {"mass_mTop_JES", "JES_mTop_JES", "mass_mTop"};
  //sample.variables = {"mass_mTop_JES_fSig", "JES_mTop_JES_fSig", "mass_mTop_fSig"};
  //sample.variables = {"mass_mTop_JES_fSig_fCP", "JES_mTop_JES_fSig_fCP", "mass_mTop"};
  sample.variables = {"mass_mTop_JES_fSig_fCP", "JES_mTop_JES_fSig_fCP", "mass_mTop_fSig_fCP"};  // DEFAULT, fSig & fCP are free
  //sample.variables = {"mass_mTop_JES_fSig", "JES_mTop_JES_fSig", "mass_mTop_fSig"};  // only fSig is free
  //sample.variables = {"mass_mTop_JES_fSig_fCP", "JES_mTop_JES_fSig_fCP", "mass_mTop_fSig"};
  //sample.variables = {"mass_mTop_JES_fSig_fCP", "mass_mTop_JES_fSig", "mass_mTop_JES"};
  //sample.variables = {"fCP_mTop_JES_fSig_fCP", "fCP_mTop_JES_fCP", "fCP_mTop_fCP"};
  //sample.variables = {"fSig_mTop_JES_fSig_fCP", "fSig_mTop_JES_fSig", "fSig_mTop_fSig"};
  //sample.variables = {"fCP_mTop_JES_fSig_fCP", "JES_mTop_JES_fSig_fCP", "fCP_mTop_fSig_fCP"};
  //sample.variables = {"fSig_mTop_JES_fSig_fCP", "JES_mTop_JES_fSig_fCP", "fSig_mTop_fSig_fCP"};

  sample.ensembles["calibration"] = ensemble("", 0, {std::make_pair(sample.variables[0],std::make_pair(172.5,0.0544808)), std::make_pair(sample.variables[1],std::make_pair(1.0,0.00046571)), std::make_pair(sample.variables[2],std::make_pair(172.5,0.0542347))});
  sample.ensembles["calibrationDummy"] = ensemble("", 0, {std::make_pair(sample.variables[0],std::make_pair(172.5,0.)), std::make_pair(sample.variables[1],std::make_pair(1.0,0.)), std::make_pair(sample.variables[2],std::make_pair(172.5,0.))});
  sample.ensembles["default"] = ensemble("Summer12_TTJetsMS1725_1.00_alljets/job_*_ensemble.root", 62131965./strafFaktor);

  sample.ensembles["puUp"  ] = ensemble("Z2_S12_PU_Up/job_*_ensemble.root"  , 62131965./strafFaktor);
  sample.ensembles["puDown"] = ensemble("Z2_S12_PU_Down/job_*_ensemble.root", 62131965./strafFaktor);

  sample.ensembles["bFragLEP"    ] = ensemble("Z2_S12_FRAG/job_*_ensemble.root"     , 62131965./strafFaktor);
  sample.ensembles["bFragLEPHard"] = ensemble("Z2_S12_FRAG_Hard/job_*_ensemble.root", 62131965./strafFaktor);
  sample.ensembles["bFragLEPSoft"] = ensemble("Z2_S12_FRAG_Soft/job_*_ensemble.root", 62131965./strafFaktor);

  sample.ensembles["bFNuUp"  ] = ensemble("Z2_S12_NU_Up/job_*_ensemble.root"  , 62131965./strafFaktor);
  sample.ensembles["bFNuDown"] = ensemble("Z2_S12_NU_Down/job_*_ensemble.root", 62131965./strafFaktor);

  sample.ensembles["flavorBUp"  ] = ensemble("Z2_S12_BJES_Up/job_*_ensemble.root"  , 62131965./strafFaktor);
  sample.ensembles["flavorBDown"] = ensemble("Z2_S12_BJES_Down/job_*_ensemble.root", 62131965./strafFaktor);

  sample.ensembles["flavorCUp"  ] = ensemble("Z2_S12_CJES_Up/job_*_ensemble.root"  , 62131965./strafFaktor);
  sample.ensembles["flavorCDown"] = ensemble("Z2_S12_CJES_Down/job_*_ensemble.root", 62131965./strafFaktor);

  sample.ensembles["flavorGUp"  ] = ensemble("Z2_S12_GJES_Up/job_*_ensemble.root"  , 62131965./strafFaktor);
  sample.ensembles["flavorGDown"] = ensemble("Z2_S12_GJES_Down/job_*_ensemble.root", 62131965./strafFaktor);

  sample.ensembles["flavorQUp"  ] = ensemble("Z2_S12_LJES_Up/job_*_ensemble.root"  , 62131965./strafFaktor);
  sample.ensembles["flavorQDown"] = ensemble("Z2_S12_LJES_Down/job_*_ensemble.root", 62131965./strafFaktor);

  sample.ensembles["bTagDiscUp"  ] = ensemble("Z2_S12_BTAG_Disc_Up/job_*_ensemble.root"  , 62131965./strafFaktor);
  sample.ensembles["bTagDiscDown"] = ensemble("Z2_S12_BTAG_Disc_Down/job_*_ensemble.root", 62131965./strafFaktor);

  sample.ensembles["bTagSFUp"    ] = ensemble("Z2_S12_BTAG_Up/job_*_ensemble.root"  , 62131965./strafFaktor);
  sample.ensembles["bTagSFDown"  ] = ensemble("Z2_S12_BTAG_Down/job_*_ensemble.root", 62131965./strafFaktor);
  sample.ensembles["misTagSFUp"  ] = ensemble("Z2_S12_MTAG_Up/job_*_ensemble.root"  , 62131965./strafFaktor);
  sample.ensembles["misTagSFDown"] = ensemble("Z2_S12_MTAG_Down/job_*_ensemble.root", 62131965./strafFaktor);

  sample.ensembles["topPt"] = ensemble("Z2_S12_TopPt/job_*_ensemble.root", 62131965./strafFaktor);

  sample.ensembles["fSigUp"  ] = ensemble("Z2_S12_FSIG_Up/job_*_ensemble.root"  , 62131965./strafFaktor);
  sample.ensembles["fSigDown"] = ensemble("Z2_S12_FSIG_Down/job_*_ensemble.root", 62131965./strafFaktor);

  sample.ensembles["matchingUp"  ] = ensemble("Z2_S12_Matching_Up/job_*_ensemble.root"  , 37083003./strafFaktor);
  sample.ensembles["matchingDown"] = ensemble("Z2_S12_Matching_Down/job_*_ensemble.root", 34053113./strafFaktor);

  sample.ensembles["scaleUp"  ] = ensemble("Z2_S12_Scale_Up/job_*_ensemble.root"  , 41908271./strafFaktor);
  sample.ensembles["scaleDown"] = ensemble("Z2_S12_Scale_Down/job_*_ensemble.root", 39286663./strafFaktor);

  sample.ensembles["jerUp"  ] = ensemble("Z2_S12_JER_Up/job_*_ensemble.root"  , 62131965./strafFaktor);
  sample.ensembles["jerDown"] = ensemble("Z2_S12_JER_Down/job_*_ensemble.root", 62131965./strafFaktor);

  sample.ensembles["jer3Up"     ] = ensemble("Z2_S12_OLDJER3_Up/job_*_ensemble.root"     , 62131965./strafFaktor);
  sample.ensembles["jer3Central"] = ensemble("Z2_S12_OLDJER3_Central/job_*_ensemble.root", 62131965./strafFaktor);
  sample.ensembles["jer3Down"   ] = ensemble("Z2_S12_OLDJER3_Down/job_*_ensemble.root"   , 62131965./strafFaktor);

  //sample.ensembles["jesFlaUp"  ] = ensemble("Z2_S12_CORFlaJES_Up/job_*_ensemble.root"  , 62131965./strafFaktor);
  //sample.ensembles["jesFlaDown"] = ensemble("Z2_S12_CORFlaJES_Down/job_*_ensemble.root", 62131965./strafFaktor);
  //sample.ensembles["jesIntUp"  ] = ensemble("Z2_S12_CORIntJES_Up/job_*_ensemble.root"  , 62131965./strafFaktor);
  //sample.ensembles["jesIntDown"] = ensemble("Z2_S12_CORIntJES_Down/job_*_ensemble.root", 62131965./strafFaktor);
  //sample.ensembles["jesMPFUp"  ] = ensemble("Z2_S12_CORMPFJES_Up/job_*_ensemble.root"  , 62131965./strafFaktor);
  //sample.ensembles["jesMPFDown"] = ensemble("Z2_S12_CORMPFJES_Down/job_*_ensemble.root", 62131965./strafFaktor);
  //sample.ensembles["jesUncUp"  ] = ensemble("Z2_S12_CORUncJES_Up/job_*_ensemble.root"  , 62131965./strafFaktor);
  //sample.ensembles["jesUncDown"] = ensemble("Z2_S12_CORUncJES_Down/job_*_ensemble.root", 62131965./strafFaktor);

  sample.ensembles["jesAbsoluteMPFBiasUp"  ] = ensemble("Z2_S12_JES_AbsoluteMPFBias_Up/job_*_ensemble.root"  , 62131965./strafFaktor);
  sample.ensembles["jesAbsoluteMPFBiasDown"] = ensemble("Z2_S12_JES_AbsoluteMPFBias_Down/job_*_ensemble.root", 62131965./strafFaktor);
  sample.ensembles["jesAbsoluteScaleUp"    ] = ensemble("Z2_S12_JES_AbsoluteScale_Up/job_*_ensemble.root"  , 62131965./strafFaktor);
  sample.ensembles["jesAbsoluteScaleDown"  ] = ensemble("Z2_S12_JES_AbsoluteScale_Down/job_*_ensemble.root", 62131965./strafFaktor);
  sample.ensembles["jesAbsoluteStatUp"     ] = ensemble("Z2_S12_JES_AbsoluteStat_Up/job_*_ensemble.root"  , 62131965./strafFaktor);
  sample.ensembles["jesAbsoluteStatDown"   ] = ensemble("Z2_S12_JES_AbsoluteStat_Down/job_*_ensemble.root", 62131965./strafFaktor);
  sample.ensembles["jesFragmentationUp"    ] = ensemble("Z2_S12_JES_Fragmentation_Up/job_*_ensemble.root"  , 62131965./strafFaktor);
  sample.ensembles["jesFragmentationDown"  ] = ensemble("Z2_S12_JES_Fragmentation_Down/job_*_ensemble.root", 62131965./strafFaktor);
  sample.ensembles["jesPileUpDataMCUp"     ] = ensemble("Z2_S12_JES_PileUpDataMC_Up/job_*_ensemble.root"  , 62131965./strafFaktor);
  sample.ensembles["jesPileUpDataMCDown"   ] = ensemble("Z2_S12_JES_PileUpDataMC_Down/job_*_ensemble.root", 62131965./strafFaktor);
  sample.ensembles["jesPileUpPtBBUp"       ] = ensemble("Z2_S12_JES_PileUpPtBB_Up/job_*_ensemble.root"  , 62131965./strafFaktor);
  sample.ensembles["jesPileUpPtBBDown"     ] = ensemble("Z2_S12_JES_PileUpPtBB_Down/job_*_ensemble.root", 62131965./strafFaktor);
  sample.ensembles["jesPileUpPtEC1Up"      ] = ensemble("Z2_S12_JES_PileUpPtEC1_Up/job_*_ensemble.root"  , 62131965./strafFaktor);
  sample.ensembles["jesPileUpPtEC1Down"    ] = ensemble("Z2_S12_JES_PileUpPtEC1_Down/job_*_ensemble.root", 62131965./strafFaktor);
  sample.ensembles["jesPileUpPtRefUp"      ] = ensemble("Z2_S12_JES_PileUpPtRef_Up/job_*_ensemble.root"  , 62131965./strafFaktor);
  sample.ensembles["jesPileUpPtRefDown"    ] = ensemble("Z2_S12_JES_PileUpPtRef_Down/job_*_ensemble.root", 62131965./strafFaktor);
  sample.ensembles["jesRelativeFSRUp"      ] = ensemble("Z2_S12_JES_RelativeFSR_Up/job_*_ensemble.root"  , 62131965./strafFaktor);
  sample.ensembles["jesRelativeFSRDown"    ] = ensemble("Z2_S12_JES_RelativeFSR_Down/job_*_ensemble.root", 62131965./strafFaktor);
  sample.ensembles["jesRelativeStatFSRUp"  ] = ensemble("Z2_S12_JES_RelativeStatFSR_Up/job_*_ensemble.root"  , 62131965./strafFaktor);
  sample.ensembles["jesRelativeStatFSRDown"] = ensemble("Z2_S12_JES_RelativeStatFSR_Down/job_*_ensemble.root", 62131965./strafFaktor);
  sample.ensembles["jesRelativeJEREC1Up"   ] = ensemble("Z2_S12_JES_RelativeJEREC1_Up/job_*_ensemble.root"  , 62131965./strafFaktor);
  sample.ensembles["jesRelativeJEREC1Down" ] = ensemble("Z2_S12_JES_RelativeJEREC1_Down/job_*_ensemble.root", 62131965./strafFaktor);
  sample.ensembles["jesRelativePtBBUp"     ] = ensemble("Z2_S12_JES_RelativePtBB_Up/job_*_ensemble.root"  , 62131965./strafFaktor);
  sample.ensembles["jesRelativePtBBDown"   ] = ensemble("Z2_S12_JES_RelativePtBB_Down/job_*_ensemble.root", 62131965./strafFaktor);
  sample.ensembles["jesRelativePtEC1Up"    ] = ensemble("Z2_S12_JES_RelativePtEC1_Up/job_*_ensemble.root"  , 62131965./strafFaktor);
  sample.ensembles["jesRelativePtEC1Down"  ] = ensemble("Z2_S12_JES_RelativePtEC1_Down/job_*_ensemble.root", 62131965./strafFaktor);
  sample.ensembles["jesSinglePionECALUp"   ] = ensemble("Z2_S12_JES_SinglePionECAL_Up/job_*_ensemble.root"  , 62131965./strafFaktor);
  sample.ensembles["jesSinglePionECALDown" ] = ensemble("Z2_S12_JES_SinglePionECAL_Down/job_*_ensemble.root", 62131965./strafFaktor);
  sample.ensembles["jesSinglePionHCALUp"   ] = ensemble("Z2_S12_JES_SinglePionHCAL_Up/job_*_ensemble.root"  , 62131965./strafFaktor);
  sample.ensembles["jesSinglePionHCALDown" ] = ensemble("Z2_S12_JES_SinglePionHCAL_Down/job_*_ensemble.root", 62131965./strafFaktor);

  sample.ensembles["jesTotalUp"  ] = ensemble("Z2_S12_JES_Total_Up/job_*_ensemble.root"  , 62131965./strafFaktor);
  sample.ensembles["jesTotalDown"] = ensemble("Z2_S12_JES_Total_Down/job_*_ensemble.root", 62131965./strafFaktor);

  sample.ensembles["triggerUp"   ] = ensemble("Z2_S12_TRIGGER_Up/job_*_ensemble.root"  , 62131965./strafFaktor);
  sample.ensembles["triggerDown" ] = ensemble("Z2_S12_TRIGGER_Down/job_*_ensemble.root", 62131965./strafFaktor);

  sample.ensembles["mcatnlo"     ] = ensemble("Z2_S12_MCNLO/job_*_ensemble.root" , 32852589./strafFaktor);
  sample.ensembles["powheg"      ] = ensemble("Z2_S12_POWHEG/job_*_ensemble.root", 21675970./strafFaktor);
  sample.ensembles["powhegHerwig"] = ensemble("Z2_S12_POWHER/job_*_ensemble.root", 27684235./strafFaktor);

  sample.ensembles["amcatnloPy8FAST"] = ensemble("Z2_S12_aMCNLO_Py8/job_*_ensemble.root" , 18280992./strafFaktor);
  sample.ensembles["amcatnloHerFAST"] = ensemble("Z2_S12_aMCNLO_Her/job_*_ensemble.root" , 22734440./strafFaktor);
  //sample.ensembles["Pythia8FAST"    ] = ensemble("Z2_S12_FSIM_Py8/job_*_ensemble.root"   ,  8660000./strafFaktor);

  //sample.ensembles["parj81CenFAST" ] = ensemble("Z2_S12_FSIM_parj81_Cen/job_*_ensemble.root"  , 60507383./strafFaktor);
  //sample.ensembles["parj81UpFAST"  ] = ensemble("Z2_S12_FSIM_parj81_Down/job_*_ensemble.root" , 63957354./strafFaktor);
  //sample.ensembles["parj81DownFAST"] = ensemble("Z2_S12_FSIM_parj81_Up/job_*_ensemble.root"   , 64220697./strafFaktor);

  sample.ensembles["defaultSC"] = ensemble("Z2_S12_Z2/job_*_ensemble.root"      , 41761265./6./6.*9.*9./strafFaktor);
  sample.ensembles["P11"      ] = ensemble("Z2_S12_P11/job_*_ensemble.root"     , 11651739./6./6.*9.*9./strafFaktor);
  sample.ensembles["P11noCR"  ] = ensemble("Z2_S12_P11NoCR/job_*_ensemble.root" , 11919063./6./6.*9.*9./strafFaktor);
  sample.ensembles["P11mpiHi" ] = ensemble("Z2_S12_P11mpiHi/job_*_ensemble.root",  7953758./6./6.*9.*9./strafFaktor);
  sample.ensembles["P11TeV"   ] = ensemble("Z2_S12_P11TeV/job_*_ensemble.root"  ,  7946264./6./6.*9.*9./strafFaktor);

  sample.ensembles["RD"] = ensemble("Z2_S12_Z2_RD/job_*_ensemble.root", 31228390./6./6.*9.*9.);

  sample.ensembles["bkgMwDef" ] = ensemble("Z2_S12_Background_Shape_mW_def/job_*_ensemble.root", 62131965./strafFaktor);
  sample.ensembles["bkgMwUp"  ] = ensemble("Z2_S12_Background_Shape_mW_Up/job_*_ensemble.root", 62131965./strafFaktor);
  sample.ensembles["bkgMwDown"] = ensemble("Z2_S12_Background_Shape_mW_Down/job_*_ensemble.root", 62131965./strafFaktor);
  sample.ensembles["bkgMtDef" ] = ensemble("Z2_S12_Background_Shape_mT_def/job_*_ensemble.root", 62131965./strafFaktor);
  sample.ensembles["bkgMtUp"  ] = ensemble("Z2_S12_Background_Shape_mT_Up/job_*_ensemble.root", 62131965./strafFaktor);
  sample.ensembles["bkgMtDown"] = ensemble("Z2_S12_Background_Shape_mT_Down/job_*_ensemble.root", 62131965./strafFaktor);

  //sample.ensembles["PDFDown"] = ensemble("", 0, {std::make_pair(sample.variables[0],std::make_pair(172.322,0.0)), std::make_pair(sample.variables[1],std::make_pair(1.0007,0.0)), std::make_pair(sample.variables[2],std::make_pair(172.337,0.0))}); // DUMMY OLD VALUES
  //sample.ensembles["PDFUp"  ] = ensemble("", 0, {std::make_pair(sample.variables[0],std::make_pair(172.347,0.0)), std::make_pair(sample.variables[1],std::make_pair(1.0010,0.0)), std::make_pair(sample.variables[2],std::make_pair(172.362,0.0))}); // DUMMY OLD VALUES

  ///////////////////////////////////

  sample.comparisons["Calibration                      "] = comparison("calibration", "calibrationDummy", "", false);
  sample.comparisons["Pile-up (pp cross-section)       "] = comparison("default", "puUp", "puDown", true);
  sample.comparisons["Jet energy response (uds)        "] = comparison("default", "flavorQUp", "flavorQDown", true);
  sample.comparisons["Jet energy response (gluon)      "] = comparison("default", "flavorGUp", "flavorGDown", true);
  sample.comparisons["Jet energy response (c)          "] = comparison("default", "flavorCUp", "flavorCDown", true);
  sample.comparisons["Jet energy response (b)          "] = comparison("default", "flavorBUp", "flavorBDown", true);
  sample.comparisons["b fragmentation                  "] = comparison("default", "bFragLEP", "", true);
  sample.comparisons["Semi-leptonic B hadron decays    "] = comparison("default", "bFNuUp", "bFNuDown", true);
  sample.comparisons["b-tag disc                       "] = comparison("default", "bTagDiscUp", "bTagDiscDown", true, false);
  sample.comparisons["b-tag rate                       "] = comparison("default", "bTagSFUp", "bTagSFDown", true);
  sample.comparisons["b-tag (mistag rate)              "] = comparison("default", "misTagSFUp", "misTagSFDown", true);
  sample.comparisons["Trigger                          "] = comparison("default", "triggerUp", "triggerDown", true);
  sample.comparisons["Top-\\pt reweighting             "] = comparison("default", "topPt", "", true);
  sample.comparisons["ME-PS matching threshold         "] = comparison("default", "matchingUp", "matchingDown", false);
  sample.comparisons["$Q^{2}$ scale                    "] = comparison("default", "scaleUp", "scaleDown", false);
  sample.comparisons["Jet energy resolution            "] = comparison("default", "jerUp" , "jerDown" , true);
  sample.comparisons["Jet energy resolution OLD        "] = comparison("jer3Central", "jer3Up" , "jer3Down" , true, false);
  //sample.comparisons["JES MPFInSitu                    "] = comparison("default", "jesMPFUp", "jesMPFDown", true, false);
  //sample.comparisons["JES Flavor                       "] = comparison("default", "jesFlaUp", "jesFlaDown", true, false);
  //sample.comparisons["JES Intercalibration             "] = comparison("default", "jesIntUp", "jesIntDown", true, false);
  //sample.comparisons["JES Uncorrelated                 "] = comparison("default", "jesUncUp", "jesUncDown", true, false);

  sample.comparisons["JES AbsoluteMPFBias              "] = comparison("default", "jesAbsoluteMPFBiasUp", "jesAbsoluteMPFBiasDown", true);
  sample.comparisons["JES AbsoluteScale                "] = comparison("default", "jesAbsoluteScaleUp", "jesAbsoluteScaleDown", true);
  sample.comparisons["JES AbsoluteStat                 "] = comparison("default", "jesAbsoluteStatUp", "jesAbsoluteStatDown", true);
  sample.comparisons["JES HighPtExtra                  "] = comparison("default", "jesHighPtExtraUp", "jesHighPtExtraDown", true);
  sample.comparisons["JES PileUpDataMC                 "] = comparison("default", "jesPileUpDataMCUp", "jesPileUpDataMCDown", true);
  sample.comparisons["JES PileUpPtBB                   "] = comparison("default", "jesPileUpPtBBUp", "jesPileUpPtBBDown", true);
  sample.comparisons["JES PileUpPtEC1                  "] = comparison("default", "jesPileUpPtEC1Up", "jesPileUpPtEC1Down", true);
  sample.comparisons["JES PileUpPtRef                  "] = comparison("default", "jesPileUpPtRefUp", "jesPileUpPtRefDown", true);
  sample.comparisons["JES RelativeFSR                  "] = comparison("default", "jesRelativeFSRUp", "jesRelativeFSRDown", true);
  sample.comparisons["JES RelativeJEREC1               "] = comparison("default", "jesRelativeJEREC1Up", "jesRelativeJEREC1Down", true);
  sample.comparisons["JES RelativePtBB                 "] = comparison("default", "jesRelativePtBBUp", "jesRelativePtBBDown", true);
  sample.comparisons["JES RelativePtEC1                "] = comparison("default", "jesRelativePtEC1Up", "jesRelativePtEC1Down", true);
  sample.comparisons["JES SinglePionECAL               "] = comparison("default", "jesSinglePionECALUp", "jesSinglePionECALDown", true);
  sample.comparisons["JES SinglePionHCAL               "] = comparison("default", "jesSinglePionHCALUp", "jesSinglePionHCALDown", true);

  sample.comparisons["JES Total                        "] = comparison("default", "jesTotalUp", "jesTotalDown", true, false);

  sample.comparisons["MadGraph vs. MC@NLO+Herwig6      "] = comparison("default", "mcatnlo", "", false, false);
  sample.comparisons["Powheg+Pythia6 vs. MC@NLO+Herwig6"] = comparison("powheg", "mcatnlo", "", false, false);
  sample.comparisons["Powheg+Pythia6 vs. Powheg+Herwig6"] = comparison("powheg", "powhegHerwig", "", false, false);
  sample.comparisons["MC@NLO+Herwig6 vs. Powheg+Herwig6"] = comparison("mcatnlo", "powhegHerwig", "", false, false);
  sample.comparisons["ME generator                     "] = comparison("default", "powheg", "", false);
  sample.comparisons["Spin correlations                "] = comparison("default", "defaultSC", "", false, false);
  sample.comparisons["Pythia Z2* vs. P11               "] = comparison("defaultSC", "P11", "", false, false);
  sample.comparisons["Color reconnection               "] = comparison("P11", "P11noCR", "", false);
  sample.comparisons["Underlying event                 "] = comparison("P11", "P11mpiHi", "P11TeV", false);
  sample.comparisons["Non-\\ttbar background \\fsig    "] = comparison("default", "fSigUp", "fSigDown", true);
  sample.comparisons["Non-\\ttbar background shape mW  "] = comparison("bkgMwDef", "bkgMwUp", "bkgMwDown", true);
  sample.comparisons["Non-\\ttbar background shape mT  "] = comparison("bkgMtDef", "bkgMtUp", "bkgMtDown", true);
  sample.comparisons["Run dependent                    "] = comparison("defaultSC", "RD", "", false, false);
  //sample.comparisons["PDF                              "] = comparison("default", "PDFDown", "PDFUp", true, true);
  //sample.comparisons["FSR                              "] = comparison("parj81CenFAST", "parj81UpFAST", "parj81DownFAST", true, true);
  sample.comparisons["aMCNLO+Pythia8 vs. aMCNLO+Herwig+"] = comparison("amcatnloPy8FAST", "amcatnloHerFAST", "", false, false);

  sample.mergedcomparisons["JEC pt/eta                      "] = mergedcomparison({
      "JES AbsoluteMPFBias              ",
      "JES AbsoluteScale                ",
      "JES AbsoluteStat                 ",
      "JES HighPtExtra                  ",
      "JES PileUpDataMC                 ",
      "JES PileUpPtBB                   ",
      "JES PileUpPtEC1                  ",
      "JES PileUpPtRef                  ",
      "JES RelativeFSR                  ",
      "JES RelativeStatFSR              ",
      "JES RelativeJEREC1               ",
      "JES RelativePtBB                 ",
      "JES RelativePtEC1                ",
      "JES SinglePionECAL               ",
      "JES SinglePionHCAL               "
    });
    
  //sample.mergedcomparisons["JEC pt/eta from groups          "] = mergedcomparison(
  //    {"JES Intercalibration             ",
  //     "JES MPFInSitu                    ",
  //     "JES Uncorrelated                 "}
  //);
  sample.mergedcomparisons["CorrelationGroupMPFInSitu       "] = mergedcomparison({
      "JES AbsoluteMPFBias              ",
    });
    
  sample.mergedcomparisons["CorrelationGroupInterCalibration"] = mergedcomparison({
      "JES RelativeFSR                  ",
    });
    
  sample.mergedcomparisons["CorrelationGroupUncorrelated    "] = mergedcomparison({
      "JES AbsoluteScale                ",
      "JES AbsoluteStat                 ",
      "JES HighPtExtra                  ",
      "JES PileUpDataMC                 ",
      "JES PileUpPtBB                   ",
      "JES PileUpPtEC1                  ",
      "JES PileUpPtRef                  ",
      "JES RelativeJEREC1               ",
      "JES RelativePtBB                 ",
      "JES RelativePtEC1                ",
      "JES SinglePionECAL               ",
      "JES SinglePionHCAL               "
      "JES RelativeStatFSR              ",
    });
    
  sample.mergedcomparisons["CorrelationGroupUncorrelatedwoPU"] = mergedcomparison({
      "JES AbsoluteScale                ",
      "JES AbsoluteStat                 ",
      "JES HighPtExtra                  ",
      "JES RelativeJEREC1               ",
      "JES RelativePtBB                 ",
      "JES RelativePtEC1                ",
      "JES SinglePionECAL               ",
      "JES SinglePionHCAL               "
      "JES RelativeStatFSR              ",
    });
    
  sample.mergedcomparisons["CorrelationGroupUncorrelatedPU  "] = mergedcomparison({
      "JES PileUpDataMC                 ",
      "JES PileUpPtBB                   ",
      "JES PileUpPtEC1                  ",
      "JES PileUpPtRef                  ",
    });
    
  sample.mergedcomparisons["JEC Flavor                      "] = mergedcomparison(
      {"Jet energy response (uds)        ",
       "Jet energy response (gluon)      ",
       "Jet energy response (c)          ",
       "Jet energy response (b)          "}, true
  );
  sample.mergedcomparisons["CorrelationGroupbJES            "] = mergedcomparison({
      "b fragmentation                  ",
      "Semi-leptonic B hadron decays    "
    });
    
  sample.mergedcomparisons["JES PileUpPt                    "] = mergedcomparison(
      {"JES PileUpPtBB                   ",
       "JES PileUpPtEC1                  ",
       "JES PileUpPtRef                  "}
  );
  sample.mergedcomparisons["PileUp                          "] = mergedcomparison(
      {"JES PileUpDataMC                 ",
       "JES PileUpPtBB                   ",
       "JES PileUpPtEC1                  ",
       "JES PileUpPtRef                  ",
       "Pile-up (pp cross-section)       "}
  );
  sample.mergedcomparisons["b-tagging                       "] = mergedcomparison(
      {"b-tag rate                       ",
       "b-tag (mistag rate)              "}
  );
  sample.mergedcomparisons["Non-\\ttbar background shape    "] = mergedcomparison(
      {"Non-\\ttbar background shape mW  ",
       "Non-\\ttbar background shape mT  "}
  );
  sample.mergedcomparisons["Non-\\ttbar background          "] = mergedcomparison(
      {"Non-\\ttbar background \\fsig    ",
       "Non-\\ttbar background shape mW  ",
       "Non-\\ttbar background shape mT  "}
  );

}

void SystematicUncertainties::deriveSystematics()
{
  //fillLeptonJets();
  fillAllJets();

  std::map<std::string, double> totalUncertainties2;
  for(auto& var : sample.variables)
    totalUncertainties2[var] = 0;

  std::cout << "\n### Fitting pseudo-experiments" << std::endl;
  for(std::map<std::string, ensemble>::iterator it = sample.ensembles.begin(); it != sample.ensembles.end(); ++it) {
    std::cout << std::setiosflags(std::ios::left) << std::setw(20) << it->first;

    int N_PE = 0;
    if (strcmp(it->second.file.c_str(), "") == 0){}
    else{
      // Get files
      it->second.chain = new TChain("tree");
      it->second.chain->Add((sample.path + it->second.file).c_str());

      // Fit
      TF1* gaus = new TF1("gaus", "gaus");

      for(auto& var : sample.variables){
        std::string selection = getSelectionFromVariable(var);
        it->second.chain->Fit("gaus", var.c_str(), selection.c_str(), "LEMQ0");
        it->second.values[var].first  = gaus->GetParameter(1);
        it->second.values[var].second = gaus->GetParameter(2) / sqrt(it->second.size/(sample.crossSection*sample.peLumi));
      }

      std::string selection = sample.variables[0] + std::string(">0 & ") + sample.variables[1] + std::string(">0 & genMass==172.5 & genJES==1");
      N_PE = it->second.chain->GetEntries(selection.c_str());
    }
    printf("\tmass = %.2lf+/-%.2lf GeV, jes = %.3lf+/-%.3lf, mass1d = %.2lf+/-%.2lf GeV, N(PE) = %i\n", it->second.values[sample.variables[0]].first, it->second.values[sample.variables[0]].second, it->second.values[sample.variables[1]].first, it->second.values[sample.variables[1]].second, it->second.values[sample.variables[2]].first, it->second.values[sample.variables[2]].second, N_PE);
  }

  ///////////////////////////////////

  std::ofstream myfile;
  char buffer[999];
  myfile.open("systematicUncertainties.txt", std::ios::out | std::ios::app);

  std::cout << "\n### Systematic uncertainties" << std::endl;
  std::cout << "Uncertainty name \t \t \t  2D mass \t  JSF \t \t    1D mass" << std::endl;
  myfile    << "Uncertainty name & 2D mass & JSF & 1D mass \\tabularnewline\n\\hline" << "\n";
  for(std::map<std::string, comparison>::iterator it = sample.comparisons.begin(); it != sample.comparisons.end(); it++) {
    myfile    << "\\hline" << "\n";
    if (!it->second.active) std::cout << "(cc) ";
    if (!it->second.active) myfile    << "(cc) ";
    std::cout << it->first;
    myfile    << it->first;

    ensemble nominal  = sample.ensembles.find(it->second.nominal)->second;
    ensemble up       = sample.ensembles.find(it->second.up     )->second;
    ensemble down;
    bool upDown = true;
    if (strcmp(it->second.down.c_str(), "") == 0) {
      down   = up;
      upDown = false;
    }
    else down = sample.ensembles.find(it->second.down)->second;

    //std::map<std::string, double> shifts;
    for(auto& var : sample.variables) {
      it->second.shifts[var] = std::max(std::abs(nominal.values[var].first-up.values[var].first), std::abs(nominal.values[var].first-down.values[var].first));
      it->second.shiftsUp  [var] =   up.values[var].first-nominal.values[var].first;
      it->second.shiftsDown[var] = down.values[var].first-nominal.values[var].first;
    }

    //std::map<std::string, double> shiftUncs;
    for(auto& var : sample.variables)
      it->second.shiftUncs[var] = 0.;

    if (!it->second.correlated) {
      for(auto& var : sample.variables)
        it->second.shiftUncs[var] = std::max(sqrt(pow(nominal.values[var].second,2)+pow(up.values[var].second,2)),sqrt(pow(nominal.values[var].second,2)+pow(down.values[var].second,2)));
    }

    sprintf(buffer," & %.2lf$\\pm$%.2lf & %.3lf$\\pm$%.3lf & %.2lf$\\pm$%.2lf \\tabularnewline\n", it->second.shifts[sample.variables[0]], it->second.shiftUncs[sample.variables[0]], it->second.shifts[sample.variables[1]], it->second.shiftUncs[sample.variables[1]], it->second.shifts[sample.variables[2]], it->second.shiftUncs[sample.variables[2]]);
    printf(" \t  %.2lf +/- %.2lf   %.3lf +/- %.3lf   %.2lf +/- %.2lf \n", it->second.shifts[sample.variables[0]], it->second.shiftUncs[sample.variables[0]], it->second.shifts[sample.variables[1]], it->second.shiftUncs[sample.variables[1]], it->second.shifts[sample.variables[2]], it->second.shiftUncs[sample.variables[2]]);
    myfile << buffer;

    if (upDown) {
      if (!it->second.active) std::cout << "(cc) ";
      sprintf(buffer, "- up                              & %+.2lf & %+.3lf & %+.2lf \\tabularnewline\n", up.values[sample.variables[0]].first-nominal.values[sample.variables[0]].first, up.values[sample.variables[1]].first-nominal.values[sample.variables[1]].first, up.values[sample.variables[2]].first-nominal.values[sample.variables[2]].first);
      printf("- up                              \t %+.2lf \t \t %+.3lf \t   %+.2lf \n", up.values[sample.variables[0]].first-nominal.values[sample.variables[0]].first, up.values[sample.variables[1]].first-nominal.values[sample.variables[1]].first, up.values[sample.variables[2]].first-nominal.values[sample.variables[2]].first);
      myfile << buffer;
      if (!it->second.active) std::cout << "(cc) ";
      sprintf(buffer, "- down                            & %+.2lf & %+.3lf & %+.2lf \\tabularnewline\n", down.values[sample.variables[0]].first-nominal.values[sample.variables[0]].first, down.values[sample.variables[1]].first-nominal.values[sample.variables[1]].first, down.values[sample.variables[2]].first-nominal.values[sample.variables[2]].first);
      printf("- down                            \t %+.2lf \t \t %+.3lf \t   %+.2lf \n", down.values[sample.variables[0]].first-nominal.values[sample.variables[0]].first, down.values[sample.variables[1]].first-nominal.values[sample.variables[1]].first, down.values[sample.variables[2]].first-nominal.values[sample.variables[2]].first);
      myfile << buffer;
    }
    else {
      if (!it->second.active) std::cout << "(cc) ";
      sprintf(buffer, "- shift                           & %+.2lf & %+.3lf & %+.2lf \\tabularnewline\n", up.values[sample.variables[0]].first-nominal.values[sample.variables[0]].first, up.values[sample.variables[1]].first-nominal.values[sample.variables[1]].first, up.values[sample.variables[2]].first-nominal.values[sample.variables[2]].first);
      printf("- shift                           \t %+.2lf \t \t %+.3lf \t   %+.2lf \n", up.values[sample.variables[0]].first-nominal.values[sample.variables[0]].first, up.values[sample.variables[1]].first-nominal.values[sample.variables[1]].first, up.values[sample.variables[2]].first-nominal.values[sample.variables[2]].first);
      myfile << buffer;
    }

    if (it->second.active) {
      for(auto& var : sample.variables){
        totalUncertainties2[var] += pow(std::max(it->second.shifts[var], it->second.shiftUncs[var]), 2);
      }
    }
  }

  std::cout << "\n### Total systematic uncertainty" << std::endl;
  myfile    << "\\hline\n\\hline\nTotal";
  sprintf(buffer, "& %.2lf & %.3lf & %.2lf \n", sqrt(totalUncertainties2[sample.variables[0]]), sqrt(totalUncertainties2[sample.variables[1]]), sqrt(totalUncertainties2[sample.variables[2]]));
  printf("\t 2D mass = %.2lf GeV, JSF = %.3lf, 1D mass = %.2lf GeV\n\n", sqrt(totalUncertainties2[sample.variables[0]]), sqrt(totalUncertainties2[sample.variables[1]]), sqrt(totalUncertainties2[sample.variables[2]]));
  myfile << buffer;

  for(std::map<std::string, mergedcomparison>::const_iterator it = sample.mergedcomparisons.begin(); it != sample.mergedcomparisons.end(); ++it) {
    std::cout << it->first;
    myfile    << it->first;

    std::map<std::string, double> mergedUncertainties2;
    std::map<std::string, double> mergedUncertaintiesUp;
    std::map<std::string, double> mergedUncertaintiesDown;
    for(auto& var : sample.variables) {
      mergedUncertainties2   [var] = 0;
      mergedUncertaintiesUp  [var] = 0;
      mergedUncertaintiesDown[var] = 0;
    }

    for(std::vector<std::string>::const_iterator it2 = it->second.comparisons.begin(); it2 != it->second.comparisons.end(); ++it2){
      for(auto& var : sample.variables){
        mergedUncertainties2[var] += pow(sample.comparisons.find(*it2)->second.shifts[var],2);
        // Add linearly, statistical uncertainties neglected
        mergedUncertaintiesUp  [var] += sample.comparisons.find(*it2)->second.shiftsUp  [var];
        mergedUncertaintiesDown[var] += sample.comparisons.find(*it2)->second.shiftsDown[var];
        
        // Subtract expected JEC shift
        //if (it->first == "JEC") {
        //  if(it2->substr(0,3) == "JES"){
        //    mergedUncertainties2[var] -= pow(0.0138, 2);
        //  }
        //}
      }
    }
    sprintf(buffer, " & %.2lf & %.3lf & %.2lf \\tabularnewline\n", sqrt(mergedUncertainties2[sample.variables[0]]), sqrt(mergedUncertainties2[sample.variables[1]]), sqrt(mergedUncertainties2[sample.variables[2]]));
    printf(" \t  %.2lf             %.3lf             %.2lf           \n", sqrt(mergedUncertainties2[sample.variables[0]]), sqrt(mergedUncertainties2[sample.variables[1]]), sqrt(mergedUncertainties2[sample.variables[2]]));
    myfile << buffer;
    
    if (it->second.linear) {
      for(auto& var : sample.variables) {
        mergedUncertaintiesUp[var] = std::abs(mergedUncertaintiesUp[var]);
	mergedUncertaintiesDown[var] = std::abs(mergedUncertaintiesDown[var]);
      }
      sprintf(buffer, "- linear & %.2lf & %.3lf & %.2lf \\tabularnewline\n", std::max(mergedUncertaintiesUp[sample.variables[0]], mergedUncertaintiesDown[sample.variables[0]]), std::max(mergedUncertaintiesUp[sample.variables[1]], mergedUncertaintiesDown[sample.variables[1]]), std::max(mergedUncertaintiesUp[sample.variables[2]], mergedUncertaintiesDown[sample.variables[2]]));
      printf("- linear                        \t  %.2lf             %.3lf             %.2lf           \n", std::max(mergedUncertaintiesUp[sample.variables[0]], mergedUncertaintiesDown[sample.variables[0]]), std::max(mergedUncertaintiesUp[sample.variables[1]], mergedUncertaintiesDown[sample.variables[1]]), std::max(mergedUncertaintiesUp[sample.variables[2]], mergedUncertaintiesDown[sample.variables[2]]));
      myfile << buffer;
      
      for(auto& var : sample.variables) {
        totalUncertainties2[var] -= mergedUncertainties2[var];
        totalUncertainties2[var] += pow(std::max(mergedUncertaintiesUp[var], mergedUncertaintiesDown[var]),2);
      }
    }
  }
  
  std::cout << "\n### Total systematic uncertainty (after linear merges)" << std::endl;
  myfile    << "\\hline\n\\hline\nTotal";
  sprintf(buffer, "& %.2lf & %.3lf & %.2lf \n", sqrt(totalUncertainties2[sample.variables[0]]), sqrt(totalUncertainties2[sample.variables[1]]), sqrt(totalUncertainties2[sample.variables[2]]));
  printf("\t 2D mass = %.2lf GeV, JSF = %.3lf, 1D mass = %.2lf GeV\n\n", sqrt(totalUncertainties2[sample.variables[0]]), sqrt(totalUncertainties2[sample.variables[1]]), sqrt(totalUncertainties2[sample.variables[2]]));
  myfile << buffer;

  myfile.close();
}

std::string SystematicUncertainties::getSelectionFromVariable(std::string& var)
{
  std::string selection = "genMass==172.5 & genJES==1";

  std::string ext = var.substr(var.find_first_of('_'));
  std::vector<std::string> freeVars;
  std::string tmp = ext;
  while(tmp.find_last_of('_') != std::string::npos){
    std::string tmp2 = tmp.substr(tmp.find_last_of('_')+1);
    if(tmp2 == "mTop") tmp2 = "mass";
    freeVars.push_back(tmp2);
    tmp = tmp.substr(0,tmp.find_last_of('_'));
  }
  for(auto& free : freeVars){
    selection += " & " + free+ext + ">0";
  }
  return selection;
}


