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
  int channel = 1;

  int selectComparisonSet = 4; //0: full set of uncertainties for Moriond; 1: only uncertainties also determined for B-Regression study; 2: minimal set used for BReg-study: 76x Test of LHE weight scale uncertainties

  std::string path;
  if(channel==3)path = "/nfs/dust/test/cms/user/kirschen/BRegression_PE_NewCalibrationJan2014Applied_MCS3250MinEvtsBReg/";
  //else path = "/nfs/dust/cms/user/kovalch/GRID-CONTROL_JOBS/pseudoexperiments/topmass_141020/";
  //else path = "/nfs/dust/cms/user/mseidel/pseudoexperiments/topmass_paper/";
  else path = "/nfs/dust/cms/user/garbersc/TopMass/2015D_TemplateCalibration/76x/";

  path += lepton_[channel]; path += "_calibrated/2015_JESVariations/76x/";

  sample.path = path;
  sample.crossSection = 832;
  sample.peLumi = 2000.;

  sample.variables = {"mass_mTop_JES", "JES_mTop_JES", "mass_mTop"};

  //TODO whats the second double in the pair?
  sample.ensembles["ref"] = ensemble("", 0, {std::make_pair("mass_mTop_JES",std::make_pair(172.5,0.0)), std::make_pair("JES_mTop_JES",std::make_pair(1.0,0.0)), std::make_pair("mass_mTop",std::make_pair(172.5,0.0))});
  sample.ensembles["calibration"] = ensemble("", 0, {std::make_pair("mass_mTop_JES",std::make_pair(172.5,0.04)), std::make_pair("JES_mTop_JES",std::make_pair(1.0,0.0005)), std::make_pair("mass_mTop",std::make_pair(172.5,0.04))});
  sample.ensembles["stat"] = ensemble("", 0, {std::make_pair("mass_mTop_JES",std::make_pair(172.5,0.194498)), std::make_pair("JES_mTop_JES",std::make_pair(1.0,0.00180928)), std::make_pair("mass_mTop",std::make_pair(172.5,0.118054))});
  
  // PDF default
  /*
  sample.ensembles["nobkg_pdfup"] = ensemble("", 0, {std::make_pair("mass_mTop_JES",std::make_pair(172.613,0.)), std::make_pair("JES_mTop_JES",std::make_pair(1.0007,0.)), std::make_pair("mass_mTop",std::make_pair(172.589,0.))});
  sample.ensembles["nobkg_pdfdown"] = ensemble("", 0, {std::make_pair("mass_mTop_JES",std::make_pair(172.486,0.)), std::make_pair("JES_mTop_JES",std::make_pair(0.9992,0.)), std::make_pair("mass_mTop",std::make_pair(172.501,0.))});
  //*/
  // PDF JSF constraint
  //*  //FIXME check this nmbrs
  sample.ensembles["nobkg_pdfup"] = ensemble("", 0, {std::make_pair("mass_mTop_JES",std::make_pair(172.574,0.)), std::make_pair("JES_mTop_JES",std::make_pair(1.0004,0.)), std::make_pair("mass_mTop",std::make_pair(172.574,0.))});
  sample.ensembles["nobkg_pdfdown"] = ensemble("", 0, {std::make_pair("mass_mTop_JES",std::make_pair(172.512,0.)), std::make_pair("JES_mTop_JES",std::make_pair(0.9996,0.)), std::make_pair("mass_mTop",std::make_pair(172.561,0.))});
  //*/
  
  sample.ensembles["default"] = ensemble("RunIISpring15MiniAODv2_76x_asymptotic_v12_ext3-v1_TT_TuneCUETP8M1_powheg-pythia8_MS1725_patJetsUpdated_1.00/job_*_ensemble.root", 97994442); //FIXME event nrs have to be effective, atm no weights are considered

  //sample.ensembles["LHEscaleWeightID_1001"] = ensemble("RunIISpring15MiniAODv2_76x_asymptotic_v12_ext3-v1_TT_TuneCUETP8M1_powheg-pythia8_MS1725_1.00_newBG_scaledWidthLHEEventWeight_0/job_*_ensemble.root", 97994442);
  //sample.ensembles["LHEscaleWeightID_1002"] = ensemble("RunIISpring15MiniAODv2_76x_asymptotic_v12_ext3-v1_TT_TuneCUETP8M1_powheg-pythia8_MS1725_1.00_newBG_scaledWidthLHEEventWeight_1/job_*_ensemble.root", 97994442);
  //sample.ensembles["LHEscaleWeightID_1003"] = ensemble("RunIISpring15MiniAODv2_76x_asymptotic_v12_ext3-v1_TT_TuneCUETP8M1_powheg-pythia8_MS1725_1.00_newBG_scaledWidthLHEEventWeight_2/job_*_ensemble.root", 97994442);
  //sample.ensembles["LHEscaleWeightID_1004"] = ensemble("RunIISpring15MiniAODv2_76x_asymptotic_v12_ext3-v1_TT_TuneCUETP8M1_powheg-pythia8_MS1725_1.00_newBG_scaledWidthLHEEventWeight_3/job_*_ensemble.root", 97994442);
  //sample.ensembles["LHEscaleWeightID_1005"] = ensemble("RunIISpring15MiniAODv2_76x_asymptotic_v12_ext3-v1_TT_TuneCUETP8M1_powheg-pythia8_MS1725_1.00_newBG_scaledWidthLHEEventWeight_4/job_*_ensemble.root", 97994442);
  //sample.ensembles["LHEscaleWeightID_1006"] = ensemble("RunIISpring15MiniAODv2_76x_asymptotic_v12_ext3-v1_TT_TuneCUETP8M1_powheg-pythia8_MS1725_1.00_newBG_scaledWidthLHEEventWeight_5/job_*_ensemble.root", 97994442);
  //sample.ensembles["LHEscaleWeightID_1007"] = ensemble("RunIISpring15MiniAODv2_76x_asymptotic_v12_ext3-v1_TT_TuneCUETP8M1_powheg-pythia8_MS1725_1.00_newBG_scaledWidthLHEEventWeight_6/job_*_ensemble.root", 97994442);
  //sample.ensembles["LHEscaleWeightID_1008"] = ensemble("RunIISpring15MiniAODv2_76x_asymptotic_v12_ext3-v1_TT_TuneCUETP8M1_powheg-pythia8_MS1725_1.00_newBG_scaledWidthLHEEventWeight_7/job_*_ensemble.root", 97994442);
  //sample.ensembles["LHEscaleWeightID_1009"] = ensemble("RunIISpring15MiniAODv2_76x_asymptotic_v12_ext3-v1_TT_TuneCUETP8M1_powheg-pythia8_MS1725_1.00_newBG_scaledWidthLHEEventWeight_8/job_*_ensemble.root", 97994442);

  sample.ensembles["scaleUp"] = ensemble("RunII_TT_TuneCUETP8M1_powheg-pythia8_76X_mcRun2_asymptotic_v12_ext3_scaleup_patJetsUpdated__1.00_allCut_scaled/job_*_ensemble.root", 38507969);
  sample.ensembles["scaleDown"] = ensemble("RunII_TT_TuneCUETP8M1_powheg-pythia8_76X_mcRun2_asymptotic_v12_ext3_scaledown_patJetsUpdated__1.00_allCut_scaled/job_*_ensemble.root", 39461147);

  sample.ensembles["puUp"] = ensemble("RunII_TT_TuneCUETP8M1_powheg-pythia8_76X_mcRun2_asymptotic_v12_ext3_wPUVarWeights_1.00_allCut_puUp/job_*_ensemble.root", 97994442);
  sample.ensembles["puDown"] = ensemble("RunII_TT_TuneCUETP8M1_powheg-pythia8_76X_mcRun2_asymptotic_v12_ext3_wPUVarWeights_1.00_allCut_puDown/job_*_ensemble.root", 97994442);
  sample.ensembles["puDef"] = ensemble("RunII_TT_TuneCUETP8M1_powheg-pythia8_76X_mcRun2_asymptotic_v12_ext3_wPUVarWeights_1.00_allCut_puDef/job_*_ensemble.root", 97994442);

  sample.ensembles["halfBG"] = ensemble("RunII_TT_TuneCUETP8M1_powheg-pythia8_76X_mcRun2_asymptotic_v12_ext3_wPUVarWeights_1.00_allCut_halfBG/job_*_ensemble.root", 97994442);
  sample.ensembles["doubleBG"] = ensemble("RunII_TT_TuneCUETP8M1_powheg-pythia8_76X_mcRun2_asymptotic_v12_ext3_wPUVarWeights_1.00_allCut_doubleBG/job_*_ensemble.root", 97994442);

  sample.ensembles["amcatnlo"] = ensemble("RunII_TT_TuneCUETP8M1_powheg-amcatnlo_76X_mcRun2_asymptotic_v12_v1_patJetsUpdated_1.00_allCut/job_*_ensemble.root", 9373386);
  sample.ensembles["amcatnloFXFX"] = ensemble("RunII_TT_TuneCUETP8M1_amcatnloFXFX-pythia8_76X_mcRun2_asymptotic_v12-v1_patJetsUpdated__1.00_allCut/job_*_ensemble.root", 9373386);
  sample.ensembles["amcatnloFXFX Q2 up"] = ensemble("RunII_TT_TuneCUETP8M1_amcatnloFXFX-pythia8_76X_mcRun2_asymptotic_v12-v1_scaleup_patJetsUpdated__1.00_allCut/job_*_ensemble.root", 9373386);
  sample.ensembles["amcatnloFXFX Q2 down"] = ensemble("RunII_TT_TuneCUETP8M1_amcatnloFXFX-pythia8_76X_mcRun2_asymptotic_v12-v1_scaledown_patJetsUpdated__1.00_allCut/job_*_ensemble.root", 9373386);


  sample.ensembles["flavor:up light"] = ensemble("RunII_TT_TuneCUETP8M1_powheg-pythia8_76X_mcRun2_asymptotic_v12_ext3_patJetsUpdated_flavor:up_FlavorPureQuark_allCut/job_*_ensemble.root", 97994442);
  sample.ensembles["flavor:down light"] = ensemble("RunII_TT_TuneCUETP8M1_powheg-pythia8_76X_mcRun2_asymptotic_v12_ext3_patJetsUpdated_flavor:down_FlavorPureQuark_allCut/job_*_ensemble.root", 97994442);
  sample.ensembles["flavor:up gluon"] = ensemble("RunII_TT_TuneCUETP8M1_powheg-pythia8_76X_mcRun2_asymptotic_v12_ext3_patJetsUpdated_flavor:up_FlavorPureGluon_allCut/job_*_ensemble.root", 97994442);
  sample.ensembles["flavor:down gluon"] = ensemble("RunII_TT_TuneCUETP8M1_powheg-pythia8_76X_mcRun2_asymptotic_v12_ext3_patJetsUpdated_flavor:down_FlavorPureGluon_allCut/job_*_ensemble.root", 97994442);
  sample.ensembles["flavor:up charm"] = ensemble("RunII_TT_TuneCUETP8M1_powheg-pythia8_76X_mcRun2_asymptotic_v12_ext3_patJetsUpdated_flavor:up_FlavorPureCharm_allCut/job_*_ensemble.root", 97994442);
  sample.ensembles["flavor:down charm"] = ensemble("RunII_TT_TuneCUETP8M1_powheg-pythia8_76X_mcRun2_asymptotic_v12_ext3_patJetsUpdated_flavor:down_FlavorPureCharm_allCut/job_*_ensemble.root", 97994442);
  sample.ensembles["flavor:up bottom"] = ensemble("RunII_TT_TuneCUETP8M1_powheg-pythia8_76X_mcRun2_asymptotic_v12_ext3_patJetsUpdated_flavor:up_FlavorPureBottom_allCut/job_*_ensemble.root", 97994442);
  sample.ensembles["flavor:down bottom"] = ensemble("RunII_TT_TuneCUETP8M1_powheg-pythia8_76X_mcRun2_asymptotic_v12_ext3_patJetsUpdated_flavor:down_FlavorPureBottom_allCut/job_*_ensemble.root", 97994442);

  sample.ensembles["source:up_TotalNoFlavorNoTime"] = ensemble("RunII_TT_TuneCUETP8M1_powheg-pythia8_76X_mcRun2_asymptotic_v12_ext3_patJetsUpdated_source:up_TotalNoFlavorNoTime_allCut/job_*_ensemble.root", 97994442);
  sample.ensembles["source:down_TotalNoFlavorNoTime"] = ensemble("RunII_TT_TuneCUETP8M1_powheg-pythia8_76X_mcRun2_asymptotic_v12_ext3_patJetsUpdated_source:down_TotalNoFlavorNoTime_allCut/job_*_ensemble.root", 97994442);

  sample.ensembles["source:up_CorrelationGroupIntercalibration"] = ensemble("RunII_TT_TuneCUETP8M1_powheg-pythia8_76X_mcRun2_asymptotic_v12_ext3_patJetsUpdated_source:up_CorrelationGroupIntercalibration_allCut/job_*_ensemble.root", 97994442);
  sample.ensembles["source:down_CorrelationGroupIntercalibration"] = ensemble("RunII_TT_TuneCUETP8M1_powheg-pythia8_76X_mcRun2_asymptotic_v12_ext3_patJetsUpdated_source:down_CorrelationGroupIntercalibration_allCut/job_*_ensemble.root", 97994442);
  sample.ensembles["source:up_CorrelationGroupUncorrelated"] = ensemble("RunII_TT_TuneCUETP8M1_powheg-pythia8_76X_mcRun2_asymptotic_v12_ext3_patJetsUpdated_source:up_CorrelationGroupUncorrelated_allCut/job_*_ensemble.root", 97994442);
  sample.ensembles["source:down_CorrelationGroupUncorrelated"] = ensemble("RunII_TT_TuneCUETP8M1_powheg-pythia8_76X_mcRun2_asymptotic_v12_ext3_patJetsUpdated_source:down_CorrelationGroupUncorrelated_allCut/job_*_ensemble.root", 97994442);
  sample.ensembles["source:up_CorrelationGroupMPFInSitu"] = ensemble("RunII_TT_TuneCUETP8M1_powheg-pythia8_76X_mcRun2_asymptotic_v12_ext3_patJetsUpdated_source:up_CorrelationGroupMPFInSitu_allCut/job_*_ensemble.root", 97994442);
  sample.ensembles["source:down_CorrelationGroupMPFInSitu"] = ensemble("RunII_TT_TuneCUETP8M1_powheg-pythia8_76X_mcRun2_asymptotic_v12_ext3_patJetsUpdated_source:down_CorrelationGroupMPFInSitu_allCut/job_*_ensemble.root", 97994442);

  sample.ensembles["JER correction"] = ensemble("RunII_TT_TuneCUETP8M1_powheg-pythia8_76X_mcRun2_asymptotic_v12_ext3_sysJER_1.00_sigRes0.0_allCut/job_*_ensemble.root", 97994442);
  sample.ensembles["JER up"] = ensemble("RunII_TT_TuneCUETP8M1_powheg-pythia8_76X_mcRun2_asymptotic_v12_ext3_sysJER_1.00_sigRes1.0_allCut/job_*_ensemble.root", 97994442);
  sample.ensembles["JER down"] = ensemble("RunII_TT_TuneCUETP8M1_powheg-pythia8_76X_mcRun2_asymptotic_v12_ext3_sysJER_1.00_sigRes-1.0_allCut/job_*_ensemble.root", 97994442);

  sample.ensembles["mpiOff"] = ensemble("RunII_TT_TuneCUETP8M1mpiOFF_powheg-pythia8_76X_mcRun2_asymptotic_v12_v1_1.00_allCut/job_*_ensemble.root", 24523000);
  sample.ensembles["noCR"] = ensemble("RunII_TT_TuneCUETP8M1noCR_powheg-pythia8_76X_mcRun2_asymptotic_v12-v1_1.00_allCut/job_*_ensemble.root", 24478000);


/*
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
  sample.ensembles["bTagCSVUp"] = ensemble("Summer12_TTJetsMS1725_csvm_0.71/job_*_ensemble.root", 62131965./1.2);
  sample.ensembles["bTagCSVDown"] = ensemble("Summer12_TTJetsMS1725_csvm_0.655/job_*_ensemble.root", 62131965./1.2);

  sample.ensembles["nobkg"] = ensemble("nobkg/job_*_ensemble.root", 62131965./1.2);
  sample.ensembles["nobkg_topPt"] = ensemble("nobkg_weight_topPt/job_*_ensemble.root", 62131965./1.2);
  sample.ensembles["topPt"] = ensemble("weight_topPt_fix/job_*_ensemble.root", 62131965./1.2);

  sample.ensembles["DibosonDown"] = ensemble("DibosonDown/job_*_ensemble.root", 62131965./1.2);
  sample.ensembles["DibosonUp"] = ensemble("DibosonUp/job_*_ensemble.root", 62131965./1.2);
  sample.ensembles["QCDDown"] = ensemble("QCDDown/job_*_ensemble.root", 62131965./1.2);
  sample.ensembles["QCDUp"] = ensemble("QCDUp/job_*_ensemble.root", 62131965./1.2);
  sample.ensembles["WJetsDown"] = ensemble("WJetsDown/job_*_ensemble.root", 62131965./1.2);
  sample.ensembles["WJetsUp"] = ensemble("WJetsUp/job_*_ensemble.root", 62131965./1.2);
  sample.ensembles["WJetsScaleDown"] = ensemble("WJetsScaleDown/job_*_ensemble.root", 62131965./1.2);
  sample.ensembles["WJetsScaleUp"] = ensemble("WJetsScaleUp/job_*_ensemble.root", 62131965./1.2);
  sample.ensembles["WJetsMatchingDown"] = ensemble("WJetsMatchingDown/job_*_ensemble.root", 62131965./1.2);
  sample.ensembles["WJetsMatchingUp"] = ensemble("WJetsMatchingUp/job_*_ensemble.root", 62131965./1.2);
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
  sample.ensembles["powheg1"] = ensemble("Summer12_TTJets1725_powheg1/job_*_ensemble.root", 21675970./1.2);
  sample.ensembles["powheg2"] = ensemble("Summer12_TTJets1725_powheg2_pythia6/job_*_ensemble.root", 21675970./1.2);
  sample.ensembles["powheg2P11C"] = ensemble("Summer12_TTJets1725_powheg2_P11C/job_*_ensemble.root", 9500000./1.2);
  sample.ensembles["powheg2Pythia8"] = ensemble("Summer12_TTJets1725_powheg2_pythia8/job_*_ensemble.root", 21675970./1.2);
  sample.ensembles["powhegHerwig"] = ensemble("Summer12_TTJets1725_powheg_herwig/job_*_ensemble.root", 27684235./1.2);
  sample.ensembles["powheg1Herwig"] = ensemble("Summer12_TTJets1725_powheg2_herwig/job_*_ensemble.root", 27684235./1.2);
  sample.ensembles["powheg2Herwig"] = ensemble("Summer12_TTJets1725_powheg2_herwig/job_*_ensemble.root", 27684235./1.2);

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
  sample.ensembles["fastparj81Up"] = ensemble("Summer12_TTJets1725_FSIM_parj81_up/job_*_ensemble.root", 64000000./1.2);
  sample.ensembles["fastNominal"] = ensemble("Summer12_TTJets1725_FSIM_parj81_nominal/job_*_ensemble.root", 64000000./1.2);
  sample.ensembles["fastparj81Down"] = ensemble("Summer12_TTJets1725_FSIM_parj81_down/job_*_ensemble.root", 64000000./1.2);
  
  sample.ensembles["fastPythia8"] = ensemble("Summer12_TT1725_FSIM_pythia8/job_*_ensemble.root", 64000000./1.2);
  sample.ensembles["fastamcHerwigppCorr"] = ensemble("Summer12_TT1725_amcatnlo_herwigpp:pythia8/job_*_ensemble.root", 64000000./1.2);
  sample.ensembles["fastamcHerwigpp"] = ensemble("Summer12_TT1725_amcatnlo_herwigpp/job_*_ensemble.root", 64000000./1.2);
  sample.ensembles["fastamcPythia8"] = ensemble("Summer12_TT1725_amcatnlo_pythia8/job_*_ensemble.root", 64000000./1.2);
  sample.ensembles["fastSherpaCluster"] = ensemble("Summer12_TT1725_sherpa2_cluster/job_*_ensemble.root", 64000000./1.2);
  sample.ensembles["fastSherpaLund"] = ensemble("Summer12_TT1725_sherpa2_lund/job_*_ensemble.root", 64000000./1.2);

  sample.ensembles["sherpa"] = ensemble("Summer12_TTJets1725_sherpa/job_*_ensemble.root", 19514018./2.*9./1.2);

  sample.ensembles["defaultSC_Z2_RD"] = ensemble("Summer12_TTJets1725_MGDecays_Z2_RD/job_*_ensemble.root", 24949110./2.*9./1.2);

  sample.ensembles["AbsoluteFlavMapUp"] = ensemble("Summer12_TTJetsMS1725_source:up_AbsoluteFlavMap_Winter14_V8/job_*_ensemble.root", 62131965./1.2);
  sample.ensembles["AbsoluteFlavMapDown"] = ensemble("Summer12_TTJetsMS1725_source:down_AbsoluteFlavMap_Winter14_V8/job_*_ensemble.root", 62131965./1.2);
  sample.ensembles["AbsoluteMPFBiasUp"] = ensemble("Summer12_TTJetsMS1725_source:up_AbsoluteMPFBias_Winter14_V8/job_*_ensemble.root", 62131965./1.2);
  sample.ensembles["AbsoluteMPFBiasDown"] = ensemble("Summer12_TTJetsMS1725_source:down_AbsoluteMPFBias_Winter14_V8/job_*_ensemble.root", 62131965./1.2);
  sample.ensembles["AbsoluteScaleUp"] = ensemble("Summer12_TTJetsMS1725_source:up_AbsoluteScale_Winter14_V8/job_*_ensemble.root", 62131965./1.2);
  sample.ensembles["AbsoluteScaleDown"] = ensemble("Summer12_TTJetsMS1725_source:down_AbsoluteScale_Winter14_V8/job_*_ensemble.root", 62131965./1.2);
  sample.ensembles["AbsoluteStatUp"] = ensemble("Summer12_TTJetsMS1725_source:up_AbsoluteStat_Winter14_V8/job_*_ensemble.root", 62131965./1.2);
  sample.ensembles["AbsoluteStatDown"] = ensemble("Summer12_TTJetsMS1725_source:down_AbsoluteStat_Winter14_V8/job_*_ensemble.root", 62131965./1.2);
  sample.ensembles["FragmentationUp"] = ensemble("Summer12_TTJetsMS1725_source:up_Fragmentation_Winter14_V8/job_*_ensemble.root", 62131965./1.2);
  sample.ensembles["FragmentationDown"] = ensemble("Summer12_TTJetsMS1725_source:down_Fragmentation_Winter14_V8/job_*_ensemble.root", 62131965./1.2);
  sample.ensembles["PileUpDataMCUp"] = ensemble("Summer12_TTJetsMS1725_source:up_PileUpDataMC_Winter14_V8/job_*_ensemble.root", 62131965./1.2);
  sample.ensembles["PileUpDataMCDown"] = ensemble("Summer12_TTJetsMS1725_source:down_PileUpDataMC_Winter14_V8/job_*_ensemble.root", 62131965./1.2);
  sample.ensembles["PileUpPtBBUp"] = ensemble("Summer12_TTJetsMS1725_source:up_PileUpPtBB_Winter14_V8/job_*_ensemble.root", 62131965./1.2);
  sample.ensembles["PileUpPtBBDown"] = ensemble("Summer12_TTJetsMS1725_source:down_PileUpPtBB_Winter14_V8/job_*_ensemble.root", 62131965./1.2);
  sample.ensembles["PileUpPtEC1Up"] = ensemble("Summer12_TTJetsMS1725_source:up_PileUpPtEC1_Winter14_V8/job_*_ensemble.root", 62131965./1.2);
  sample.ensembles["PileUpPtEC1Down"] = ensemble("Summer12_TTJetsMS1725_source:down_PileUpPtEC1_Winter14_V8/job_*_ensemble.root", 62131965./1.2);
  sample.ensembles["PileUpPtEC2Up"] = ensemble("Summer12_TTJetsMS1725_source:up_PileUpPtEC2_Winter14_V8/job_*_ensemble.root", 62131965./1.2);
  sample.ensembles["PileUpPtEC2Down"] = ensemble("Summer12_TTJetsMS1725_source:down_PileUpPtEC2_Winter14_V8/job_*_ensemble.root", 62131965./1.2);
  sample.ensembles["PileUpPtHFUp"] = ensemble("Summer12_TTJetsMS1725_source:up_PileUpPtHF_Winter14_V8/job_*_ensemble.root", 62131965./1.2);
  sample.ensembles["PileUpPtHFDown"] = ensemble("Summer12_TTJetsMS1725_source:down_PileUpPtHF_Winter14_V8/job_*_ensemble.root", 62131965./1.2);
  sample.ensembles["PileUpPtRefUp"] = ensemble("Summer12_TTJetsMS1725_source:up_PileUpPtRef_Winter14_V8/job_*_ensemble.root", 62131965./1.2);
  sample.ensembles["PileUpPtRefDown"] = ensemble("Summer12_TTJetsMS1725_source:down_PileUpPtRef_Winter14_V8/job_*_ensemble.root", 62131965./1.2);
  sample.ensembles["RelativeFSRUp"] = ensemble("Summer12_TTJetsMS1725_source:up_RelativeFSR_Winter14_V8/job_*_ensemble.root", 62131965./1.2);
  sample.ensembles["RelativeFSRDown"] = ensemble("Summer12_TTJetsMS1725_source:down_RelativeFSR_Winter14_V8/job_*_ensemble.root", 62131965./1.2);
  sample.ensembles["RelativeJEREC1Up"] = ensemble("Summer12_TTJetsMS1725_source:up_RelativeJEREC1_Winter14_V8/job_*_ensemble.root", 62131965./1.2);
  sample.ensembles["RelativeJEREC1Down"] = ensemble("Summer12_TTJetsMS1725_source:down_RelativeJEREC1_Winter14_V8/job_*_ensemble.root", 62131965./1.2);
  sample.ensembles["RelativeJEREC2Up"] = ensemble("Summer12_TTJetsMS1725_source:up_RelativeJEREC2_Winter14_V8/job_*_ensemble.root", 62131965./1.2);
  sample.ensembles["RelativeJEREC2Down"] = ensemble("Summer12_TTJetsMS1725_source:down_RelativeJEREC2_Winter14_V8/job_*_ensemble.root", 62131965./1.2);
  sample.ensembles["RelativeJERHFUp"] = ensemble("Summer12_TTJetsMS1725_source:up_RelativeJERHF_Winter14_V8/job_*_ensemble.root", 62131965./1.2);
  sample.ensembles["RelativeJERHFDown"] = ensemble("Summer12_TTJetsMS1725_source:down_RelativeJERHF_Winter14_V8/job_*_ensemble.root", 62131965./1.2);
  sample.ensembles["RelativePtBBUp"] = ensemble("Summer12_TTJetsMS1725_source:up_RelativePtBB_Winter14_V8/job_*_ensemble.root", 62131965./1.2);
  sample.ensembles["RelativePtBBDown"] = ensemble("Summer12_TTJetsMS1725_source:down_RelativePtBB_Winter14_V8/job_*_ensemble.root", 62131965./1.2);
  sample.ensembles["RelativePtEC1Up"] = ensemble("Summer12_TTJetsMS1725_source:up_RelativePtEC1_Winter14_V8/job_*_ensemble.root", 62131965./1.2);
  sample.ensembles["RelativePtEC1Down"] = ensemble("Summer12_TTJetsMS1725_source:down_RelativePtEC1_Winter14_V8/job_*_ensemble.root", 62131965./1.2);
  sample.ensembles["RelativePtEC2Up"] = ensemble("Summer12_TTJetsMS1725_source:up_RelativePtEC2_Winter14_V8/job_*_ensemble.root", 62131965./1.2);
  sample.ensembles["RelativePtEC2Down"] = ensemble("Summer12_TTJetsMS1725_source:down_RelativePtEC2_Winter14_V8/job_*_ensemble.root", 62131965./1.2);
  sample.ensembles["RelativePtHFUp"] = ensemble("Summer12_TTJetsMS1725_source:up_RelativePtHF_Winter14_V8/job_*_ensemble.root", 62131965./1.2);
  sample.ensembles["RelativePtHFDown"] = ensemble("Summer12_TTJetsMS1725_source:down_RelativePtHF_Winter14_V8/job_*_ensemble.root", 62131965./1.2);
  sample.ensembles["RelativeStatFSRUp"] = ensemble("Summer12_TTJetsMS1725_source:up_RelativeStatFSR_Winter14_V8/job_*_ensemble.root", 62131965./1.2);
  sample.ensembles["RelativeStatFSRDown"] = ensemble("Summer12_TTJetsMS1725_source:down_RelativeStatFSR_Winter14_V8/job_*_ensemble.root", 62131965./1.2);
  sample.ensembles["RelativeStatEC2Up"] = ensemble("Summer12_TTJetsMS1725_source:up_RelativeStatEC2_Winter14_V8/job_*_ensemble.root", 62131965./1.2);
  sample.ensembles["RelativeStatEC2Down"] = ensemble("Summer12_TTJetsMS1725_source:down_RelativeStatEC2_Winter14_V8/job_*_ensemble.root", 62131965./1.2);
  sample.ensembles["RelativeStatHFUp"] = ensemble("Summer12_TTJetsMS1725_source:up_RelativeStatHF_Winter14_V8/job_*_ensemble.root", 62131965./1.2);
  sample.ensembles["RelativeStatHFDown"] = ensemble("Summer12_TTJetsMS1725_source:down_RelativeStatHF_Winter14_V8/job_*_ensemble.root", 62131965./1.2);
  sample.ensembles["SinglePionECALUp"] = ensemble("Summer12_TTJetsMS1725_source:up_SinglePionECAL_Winter14_V8/job_*_ensemble.root", 62131965./1.2);
  sample.ensembles["SinglePionECALDown"] = ensemble("Summer12_TTJetsMS1725_source:down_SinglePionECAL_Winter14_V8/job_*_ensemble.root", 62131965./1.2);
  sample.ensembles["SinglePionHCALUp"] = ensemble("Summer12_TTJetsMS1725_source:up_SinglePionHCAL_Winter14_V8/job_*_ensemble.root", 62131965./1.2);
  sample.ensembles["SinglePionHCALDown"] = ensemble("Summer12_TTJetsMS1725_source:down_SinglePionHCAL_Winter14_V8/job_*_ensemble.root", 62131965./1.2);
*/
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
      sample.comparisons["Background: WJets scale          "] = comparison("default", "WJetsScaleUp", "WJetsScaleDown");
      sample.comparisons["Background: WJets matching       "] = comparison("default", "WJetsMatchingUp", "WJetsMatchingDown");
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
      sample.comparisons["b-tag rate                       "] = comparison("default", "bTagSFUp", "bTagSFDown", true, false);
      sample.comparisons["b-tag (mistag rate)              "] = comparison("default", "misTagSFUp", "misTagSFDown", true, false);
      sample.comparisons["b-tag disc                       "] = comparison("default", "bTagCSVUp", "bTagCSVDown");
      //sample.comparisons["Top-pt reweighting (nobkg)       "] = comparison("nobkg", "nobkg_topPt", "", true, false);
      sample.comparisons["Top-pt reweighting               "] = comparison("default", "topPt", "");
      sample.comparisons["PDF                              "] = comparison("nobkg", "nobkg_pdfup", "nobkg_pdfdown");
      sample.comparisons["ME-PS matching threshold         "] = comparison("calibration", "matchingUp", "matchingDown", false);
      sample.comparisons["$Q^{2}$ scale                    "] = comparison("calibration", "scaleUp", "scaleDown", false);
      sample.comparisons["Powheg+Pythia6 vs. MC@NLO+Herwig6"] = comparison("powheg", "mcatnlo", "", false, false);
      sample.comparisons["Powheg+Pythia6 vs. Powheg+Herwig6"] = comparison("powheg", "powhegHerwig", "", false, false);
      sample.comparisons["MC@NLO+Herwig6 vs. Powheg+Herwig6"] = comparison("mcatnlo", "powhegHerwig", "", false, false);
      sample.comparisons["ME generator                     "] = comparison("calibration", "powheg", "", false);
      sample.comparisons["Pythia Z2* vs. P11               "] = comparison("calibration", "P11", "", false, false);
      //sample.comparisons["D0 crosscheck                    "] = comparison("calibration", "powhegP11C", "", false, false);
      sample.comparisons["Color reconnection               "] = comparison("P11", "P11noCR", "", false);
      sample.comparisons["Underlying event                 "] = comparison("P11", "P11mpiHi", "P11TeV", false);


      ///////
      sample.comparisons["Run dependency                   "] = comparison("defaultSC", "defaultSC_Z2_RD", "", true, false);
      //sample.comparisons["FSR scale                        "] = comparison("fastNominal", "fastparj81Up", "fastparj81Down", false);
      sample.comparisons["sherpa                           "] = comparison("calibration", "sherpa", "", false, false);
      sample.comparisons["AbsoluteFlavMap                  "] = comparison("default", "AbsoluteFlavMapUp", "AbsoluteFlavMapDown");
      sample.comparisons["AbsoluteMPFBias                  "] = comparison("default", "AbsoluteMPFBiasUp", "AbsoluteMPFBiasDown");
      sample.comparisons["AbsoluteScale                    "] = comparison("default", "AbsoluteScaleUp", "AbsoluteScaleDown");
      sample.comparisons["AbsoluteStat                     "] = comparison("default", "AbsoluteStatUp", "AbsoluteStatDown");
      sample.comparisons["Fragmentation                    "] = comparison("default", "FragmentationUp", "FragmentationDown");
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
      sample.comparisons["RelativeStatFSR                  "] = comparison("default", "RelativeStatFSRUp", "RelativeStatFSRDown");
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

    case 4:
    //	sample.comparisons["Factorisation                   "] = comparison("default", "LHEscaleWeightID_1002", "LHEscaleWeightID_1003", true,false);
    //	sample.comparisons["Renormalisation                 "] = comparison("default", "LHEscaleWeightID_1004", "LHEscaleWeightID_1007", true,false);
    	//sample.comparisons["$Q^{2}$ scale via weight        "] = comparison("default", "LHEscaleWeightID_1005", "LHEscaleWeightID_1009", true,false);
    	//sample.comparisons["$Q^{2}$ scale   active?         "] = comparison("default", "LHEscaleWeightID_1005", "LHEscaleWeightID_1009", true);
    	//sample.comparisons["$Q^{2}$ scale   uncor?          "] = comparison("default", "LHEscaleWeightID_1005", "LHEscaleWeightID_1009", false);
    	//sample.comparisons["$Q^{2}$ scale  active and uncot?"] = comparison("default", "LHEscaleWeightID_1005", "LHEscaleWeightID_1009", false);
    //	sample.comparisons["$Q^{2}$ scale to 1001 weight    "] = comparison("LHEscaleWeightID_1001", "LHEscaleWeightID_1005", "LHEscaleWeightID_1009", true,false);
    	sample.comparisons["$Q^{2}$ scale         "] = comparison("default", "scaleUp", "scaleDown", true,false);
    	//sample.comparisons["Pileup, wrong BG      "] = comparison("default", "puUp", "puDown", true,false);
    	sample.comparisons["rough BG overestimate "] = comparison("default", "doubleBG", "halfBG", true,false);
    	sample.comparisons["Pileup, no BG         "] = comparison("puDef", "puUp", "puDown", true,false);
    	//sample.comparisons["oldschool $Q^{2}$ scale uncor?  "] = comparison("default", "scaleUp", "scaleDown", false,false);
    	//sample.comparisons["aMCatNLO, same configuration    "] = comparison("default", "amcatnlo", "", false,false);
    	//sample.comparisons["aMCatNLOFXFX, same configuration"] = comparison("default", "amcatnloFXFX", "", false,false);
    	//sample.comparisons["aMCatNLOFXFX $Q^{2}$ scale      "] = comparison("amcatnloFXFX", "amcatnloFXFX Q2 up", "amcatnloFXFX Q2 down", false,false);
    	sample.comparisons["flavor light                    "] = comparison("default", "flavor:up light", "flavor:down light", true,false);
    	//sample.comparisons["calibration                     "] = comparison("default", "calibration", "", true,false);
    	sample.comparisons["flavor gluon                    "] = comparison("default", "flavor:up gluon", "flavor:down gluon", true,false);
    	sample.comparisons["flavor charm                    "] = comparison("default", "flavor:up charm", "flavor:down charm", true,false);
    	sample.comparisons["flavor bottom                   "] = comparison("default", "flavor:up bottom", "flavor:down bottom", true,false);
       	//sample.comparisons["source:up_TotalNoFlavorNoTime   "] = comparison("default", "source:up_TotalNoFlavorNoTime", "source:down_TotalNoFlavorNoTime", true,false);
       	sample.comparisons["JEC : Intercalibration          "] = comparison("default", "source:up_CorrelationGroupIntercalibration", "source:down_CorrelationGroupIntercalibration", true,false);
       	sample.comparisons["JEC : MPFInSitu                 "] = comparison("default", "source:up_CorrelationGroupMPFInSitu", "source:down_CorrelationGroupMPFInSitu", true,false);
       	sample.comparisons["JEC : Uncorrelated              "] = comparison("default", "source:up_CorrelationGroupUncorrelated", "source:down_CorrelationGroupUncorrelated", true,false);
       	sample.comparisons["JER correction                  "] = comparison("default", "JER correction", "", true,false);
       	sample.comparisons["JER                             "] = comparison("JER correction", "JER up", "JER down", true,false);
       	//sample.comparisons["Underlying Event                "] = comparison("default", "mpiOff", "", false,false);
       	//sample.comparisons["Color reconnection modelling    "] = comparison("default", "noCR", "", false,false);
       	sample.comparisons["Underlying Event, corrected JER "] = comparison("JER correction", "mpiOff", "", false,false);
       	sample.comparisons["ColorReconnectionModellingCorJER"] = comparison("JER correction", "noCR", "", false,false);
    	break;
    }

  sample.mergedcomparisons["JEC: Flavor           "] = mergedcomparison({
          "flavor light                    ",
          "flavor gluon                    ",
  		  "flavor charm                    ",
		  "flavor bottom                   "
        },true);

  sample.mergedcomparisons["JEC quad. sum         "] = mergedcomparison({
	  "JEC : Intercalibration          ",
	  "JEC : MPFInSitu                 ",
	  "JEC : Uncorrelated              "
  });

  sample.mergedcomparisons["nonlin-merged         "] = mergedcomparison({
	  "$Q^{2}$ scale         ",
	  "Pileup, no BG         ",
	  "JER                             ",
	  "Underlying Event, corrected JER ",
	  "ColorReconnectionModellingCorJER",
      "flavor light                    ",
      "flavor gluon                    ",
	"flavor charm                    ",
	  "flavor bottom                   "
  });


 /* sample.mergedcomparisons["all weight tests         "] = mergedcomparison({
        "Factorisation                   ",
        "Renormalisation                 ",
		"$Q^{2}$ scale                   ",
		//"$Q^{2}$ scale anti-correlated   "
      });

  sample.mergedcomparisons["weight to gen scale       "] = mergedcomparison({
		"$Q^{2}$ scale                   ",
		"$Q^{2}$ scale on gen lvl        "
      });
*/
/*    sample.mergedcomparisons["JEC pt/eta                      "] = mergedcomparison({
      "AbsoluteFlavMap                  ",
      "AbsoluteMPFBias                  ",
      "AbsoluteScale                    ",
      "AbsoluteStat                     ",
      "Fragmentation                    ",
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
      "Fragmentation                    ",
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
      "RelativeStatFSR                  ",
      "RelativeStatEC2                  ",
      "RelativeStatHF                   ",
      "SinglePionECAL                   ",
      "SinglePionHCAL                   "
    });
    
    sample.mergedcomparisons["CorrelationGroupUncorrelatedwoPU"] = mergedcomparison({
      "AbsoluteFlavMap                  ",
      "AbsoluteScale                    ",
      "AbsoluteStat                     ",
      "Fragmentation                    ",
      "RelativeJEREC1                   ",
      "RelativeJEREC2                   ",
      "RelativeJERHF                    ",
      "RelativePtBB                     ",
      "RelativePtEC1                    ",
      "RelativePtEC2                    ",
      "RelativePtHF                     ",
      "RelativeStatFSR                  ",
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
      "Background: WJets scale          ",
      "Background: WJets matching       ",
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
*/
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
  //sample.path =
  //"/nfs/dust/cms/user/eschliec/TopMass/topmass_150211_1802/"; //
  //older DEFAULT WITH JSF-CONSTRAINT FOR PAPER
  //sample.path="/nfs/dust/cms/user/stadie/TopMass/syscHSt0/";
  //sample.path="/nfs/dust/cms/user/stadie/TopMass/sysct0/";
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
  //sample.variables = {"fSig_mTop_JES_fSig", "JES_mTop_JES_fSig", "fSig_mTop_fSig"};

  sample.ensembles["calibration"] = ensemble("", 0, {std::make_pair(sample.variables[0],std::make_pair(172.5,0.0616207)), std::make_pair(sample.variables[1],std::make_pair(1.0,0.000531145)), std::make_pair(sample.variables[2],std::make_pair(172.5,0.0607915))});
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

  //sample.ensembles["bTagDiscUp"  ] = ensemble("Z2_S12_BTAG_Disc_Up/job_*_ensemble.root"  , 62131965./strafFaktor);
  //sample.ensembles["bTagDiscDown"] = ensemble("Z2_S12_BTAG_Disc_Down/job_*_ensemble.root", 62131965./strafFaktor);

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

  //sample.ensembles["amcatnloPy8FAST"] = ensemble("Z2_S12_aMCNLO_Py8/job_*_ensemble.root" , 18280992./strafFaktor);
  //sample.ensembles["amcatnloHerFAST"] = ensemble("Z2_S12_aMCNLO_Her/job_*_ensemble.root" , 22734440./strafFaktor);
  //sample.ensembles["Pythia8FAST"    ] = ensemble("Z2_S12_FSIM_Py8/job_*_ensemble.root"   ,  8660000./strafFaktor);

  //sample.ensembles["parj81CenFAST" ] = ensemble("Z2_S12_FSIM_parj81_Cen/job_*_ensemble.root"  , 60507383./strafFaktor);
  //sample.ensembles["parj81UpFAST"  ] = ensemble("Z2_S12_FSIM_parj81_Down/job_*_ensemble.root" , 63957354./strafFaktor);
  //sample.ensembles["parj81DownFAST"] = ensemble("Z2_S12_FSIM_parj81_Up/job_*_ensemble.root"   , 64220697./strafFaktor);

  //sample.ensembles["defaultSC"] = ensemble("Z2_S12_Z2/job_*_ensemble.root"      , 41761265./6./6.*9.*9./strafFaktor);
  sample.ensembles["P11"      ] = ensemble("Z2_S12_P11/job_*_ensemble.root"     , 11651739./6./6.*9.*9./strafFaktor);
  sample.ensembles["P11noCR"  ] = ensemble("Z2_S12_P11NoCR/job_*_ensemble.root" , 11919063./6./6.*9.*9./strafFaktor);
  sample.ensembles["P11mpiHi" ] = ensemble("Z2_S12_P11mpiHi/job_*_ensemble.root",  7953758./6./6.*9.*9./strafFaktor);
  sample.ensembles["P11TeV"   ] = ensemble("Z2_S12_P11TeV/job_*_ensemble.root"  ,  7946264./6./6.*9.*9./strafFaktor);

  sample.ensembles["RD"] = ensemble("Z2_S12_Z2_RD/job_*_ensemble.root", 31228390./6./6.*9.*9.);

  sample.ensembles["bkgMwDef" ] = ensemble("Z2_S12_Background_Shape_mW_def/job_*_ensemble.root", 62131965./strafFaktor);
  sample.ensembles["bkgMwUp"  ] = ensemble("Z2_S12_Background_Shape_mW_Up/job_*_ensemble.root", 62131965./strafFaktor);
  sample.ensembles["bkgMwDown"] = ensemble("Z2_S12_Background_Shape_mW_Down/job_*_ensemble.root", 62131965./strafFaktor);
  sample.ensembles["bkgMwAlt" ] = ensemble("Z2_S12_Background_Shape_mW_alt/job_*_ensemble.root", 62131965./strafFaktor);
  sample.ensembles["bkgMtDef" ] = ensemble("Z2_S12_Background_Shape_mT_def/job_*_ensemble.root", 62131965./strafFaktor);
  sample.ensembles["bkgMtUp"  ] = ensemble("Z2_S12_Background_Shape_mT_Up/job_*_ensemble.root", 62131965./strafFaktor);
  sample.ensembles["bkgMtDown"] = ensemble("Z2_S12_Background_Shape_mT_Down/job_*_ensemble.root", 62131965./strafFaktor);
  sample.ensembles["bkgMtAlt"]  = ensemble("Z2_S12_Background_Shape_mT_alt/job_*_ensemble.root", 62131965./strafFaktor);


  //external PDF uncertainties
  //2D
  //sample.ensembles["PDFref"] =  ensemble("Z2_S12_ABS_JES_100_172_5/pdf_weight/148/job_*_ensemble.root", 62131965./strafFaktor);
  //sample.ensembles["PDFDown"] = ensemble("", 0, {std::make_pair(sample.variables[0],std::make_pair(172.528,0.0)), std::make_pair(sample.variables[1],std::make_pair(0.9970,0.0)), std::make_pair(sample.variables[2],std::make_pair(172.332,0.0))});
  //sample.ensembles["PDFUp"  ] = ensemble("", 0, {std::make_pair(sample.variables[0],std::make_pair(172.619,0.0)), std::make_pair(sample.variables[1],std::make_pair(0.9978,0.0)), std::make_pair(sample.variables[2],std::make_pair(172.406,0.0))});
  //constraint 2D
  sample.ensembles["PDFref"] =  ensemble("Z2_S12_ABS_JES_100_172_5/pdf_weight/148/job_*_ensemble.root", 62131965./strafFaktor);
  sample.ensembles["PDFDown"] = ensemble("", 0, {std::make_pair(sample.variables[0],std::make_pair(172.372,0.0)), std::make_pair(sample.variables[1],std::make_pair(0.9994,0.0)), std::make_pair(sample.variables[2],std::make_pair(172.333,0.0))});
  sample.ensembles["PDFUp"  ] = ensemble("", 0, {std::make_pair(sample.variables[0],std::make_pair(172.435,0.0)), std::make_pair(sample.variables[1],std::make_pair(0.9996,0.0)), std::make_pair(sample.variables[2],std::make_pair(172.391,0.0))});

  ///////////////////////////////////

  sample.comparisons["Calibration                      "] = comparison("calibration", "calibrationDummy", "", false);
  sample.comparisons["Pile-up (pp cross-section)       "] = comparison("default", "puUp", "puDown", true);
  sample.comparisons["Jet energy response (uds)        "] = comparison("default", "flavorQUp", "flavorQDown", true);
  sample.comparisons["Jet energy response (gluon)      "] = comparison("default", "flavorGUp", "flavorGDown", true);
  sample.comparisons["Jet energy response (c)          "] = comparison("default", "flavorCUp", "flavorCDown", true);
  sample.comparisons["Jet energy response (b)          "] = comparison("default", "flavorBUp", "flavorBDown", true);
  sample.comparisons["b fragmentation                  "] = comparison("default", "bFragLEP", "", true);
  sample.comparisons["Semi-leptonic B hadron decays    "] = comparison("default", "bFNuUp", "bFNuDown", true);
  //sample.comparisons["b-tag disc                       "] = comparison("default", "bTagDiscUp", "bTagDiscDown", true, false);
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
  sample.comparisons["JES Fragmentation                "] = comparison("default", "jesFragmentationUp", "jesFragmentationDown", true);
  sample.comparisons["JES PileUpDataMC                 "] = comparison("default", "jesPileUpDataMCUp", "jesPileUpDataMCDown", true);
  sample.comparisons["JES PileUpPtBB                   "] = comparison("default", "jesPileUpPtBBUp", "jesPileUpPtBBDown", true);
  sample.comparisons["JES PileUpPtEC1                  "] = comparison("default", "jesPileUpPtEC1Up", "jesPileUpPtEC1Down", true);
  sample.comparisons["JES PileUpPtRef                  "] = comparison("default", "jesPileUpPtRefUp", "jesPileUpPtRefDown", true);
  sample.comparisons["JES RelativeFSR                  "] = comparison("default", "jesRelativeFSRUp", "jesRelativeFSRDown", true);
  sample.comparisons["JES RelativeStatFSR              "] = comparison("default", "jesRelativeStatFSRUp", "jesRelativeStatFSRDown", true);
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
  //sample.comparisons["Spin correlations                "] = comparison("default", "defaultSC", "", false, false);
  //sample.comparisons["Pythia Z2* vs. P11               "] = comparison("defaultSC", "P11", "", false, false);
  sample.comparisons["Color reconnection               "] = comparison("P11", "P11noCR", "", false);
  sample.comparisons["Underlying event                 "] = comparison("P11", "P11mpiHi", "P11TeV", false);
  sample.comparisons["Non-\\ttbar background \\fsig    "] = comparison("default", "fSigUp", "fSigDown", true);
  sample.comparisons["Non-\\ttbar background shape mW  "] = comparison("bkgMwDef", "bkgMwUp", "bkgMwDown", true);
  sample.comparisons["Non-\\ttbar background shape mT  "] = comparison("bkgMtDef", "bkgMtUp", "bkgMtDown", true);
  sample.comparisons["Run dependent                    "] = comparison("default", "RD", "", false, false);
  sample.comparisons["PDF                              "] = comparison("PDFref", "PDFDown", "PDFUp", true, true);
  //sample.comparisons["FSR                              "] = comparison("parj81CenFAST", "parj81UpFAST", "parj81DownFAST", true, true);
  //sample.comparisons["aMCNLO+Pythia8 vs. aMCNLO+Herwig+"] = comparison("amcatnloPy8FAST", "amcatnloHerFAST", "", false, false);

  sample.mergedcomparisons["JEC pt/eta                      "] = mergedcomparison({
      "JES AbsoluteMPFBias              ",
      "JES AbsoluteScale                ",
      "JES AbsoluteStat                 ",
      "JES Fragmentation                ",
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
      "JES Fragmentation                ",
      "JES PileUpDataMC                 ",
      "JES PileUpPtBB                   ",
      "JES PileUpPtEC1                  ",
      "JES PileUpPtRef                  ",
      "JES RelativeJEREC1               ",
      "JES RelativePtBB                 ",
      "JES RelativePtEC1                ",
      "JES SinglePionECAL               ",
      "JES SinglePionHCAL               ",
      "JES RelativeStatFSR              ",
    });
    
  sample.mergedcomparisons["CorrelationGroupUncorrelatedwoPU"] = mergedcomparison({
      "JES AbsoluteScale                ",
      "JES AbsoluteStat                 ",
      "JES Fragmentation                ",
      "JES RelativeJEREC1               ",
      "JES RelativePtBB                 ",
      "JES RelativePtEC1                ",
      "JES SinglePionECAL               ",
      "JES SinglePionHCAL               ",
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
  fillLeptonJets();
  //fillAllJets();

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

    	// std::cout<<"debug "<<(sample.path + it->second.file)<<std::endl;

      it->second.chain = new TChain("tree");
      it->second.chain->Add((sample.path + it->second.file).c_str());
      // Fit
      TF1* gaus = new TF1("gaus", "gaus");

      for(auto& var : sample.variables){
        std::string selection = getSelectionFromVariable(var);
     //   std::cout<<"debug var:"<<var<<"    selection:"<<selection<<std::endl;
        it->second.chain->Fit("gaus", var.c_str(), selection.c_str(), "LEMQ0");
        it->second.values[var].first  = gaus->GetParameter(1);
        it->second.values[var].second = gaus->GetParameter(2) / sqrt(it->second.size/(sample.crossSection*sample.peLumi));
      }

      std::string selection = sample.variables[0] + std::string(">0 & ") + sample.variables[1] + std::string(">0 & genMass==172.5 & genJES==1");
      N_PE = it->second.chain->GetEntries(selection.c_str());
    }
    printf("\tmass = %.3lf+/-%.3lf GeV, jes = %.3lf+/-%.3lf, mass1d = %.3lf+/-%.3lf GeV, N(PE) = %i\n", it->second.values[sample.variables[0]].first, it->second.values[sample.variables[0]].second, it->second.values[sample.variables[1]].first, it->second.values[sample.variables[1]].second, it->second.values[sample.variables[2]].first, it->second.values[sample.variables[2]].second, N_PE);
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
    else{
    	down = sample.ensembles.find(it->second.down)->second;
    }

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
    std::map<std::string, double> mergedUncertaintiesUpSign;
    std::map<std::string, double> mergedUncertaintiesDown;
    for(auto& var : sample.variables) {
      mergedUncertainties2   [var] = 0;
      mergedUncertaintiesUp  [var] = 0;
      mergedUncertaintiesUpSign[var] = 0;
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
    
    for(auto& var : sample.variables) {
      mergedUncertaintiesUpSign[var] = mergedUncertaintiesUp[var]/std::abs(mergedUncertaintiesUp[var]);
    }
    
    sprintf(buffer, " & %.2lf & %.3lf & %.2lf \\tabularnewline\n", sqrt(mergedUncertainties2[sample.variables[0]]), sqrt(mergedUncertainties2[sample.variables[1]]), sqrt(mergedUncertainties2[sample.variables[2]]));
    printf(" \t  %+.2lf             %+.3lf             %+.2lf           \n",
      mergedUncertaintiesUpSign[sample.variables[0]]*sqrt(mergedUncertainties2[sample.variables[0]]), 
      mergedUncertaintiesUpSign[sample.variables[1]]*sqrt(mergedUncertainties2[sample.variables[1]]), 
      mergedUncertaintiesUpSign[sample.variables[2]]*sqrt(mergedUncertainties2[sample.variables[2]]));
    myfile << buffer;
    
    if (it->second.linear) {
      for(auto& var : sample.variables) {
        mergedUncertaintiesUp  [var] = std::abs(mergedUncertaintiesUp  [var]);
        mergedUncertaintiesDown[var] = std::abs(mergedUncertaintiesDown[var]);
      }
      sprintf(buffer, "- linear & %.2lf & %.3lf & %.2lf \\tabularnewline\n",
        mergedUncertaintiesUpSign[sample.variables[0]]*std::max(mergedUncertaintiesUp[sample.variables[0]], mergedUncertaintiesDown[sample.variables[0]]),
        mergedUncertaintiesUpSign[sample.variables[1]]*std::max(mergedUncertaintiesUp[sample.variables[1]], mergedUncertaintiesDown[sample.variables[1]]),
        mergedUncertaintiesUpSign[sample.variables[2]]*std::max(mergedUncertaintiesUp[sample.variables[2]], mergedUncertaintiesDown[sample.variables[2]]));
      printf("- linear                        \t  %+.2lf             %+.3lf             %+.2lf           \n",
        mergedUncertaintiesUpSign[sample.variables[0]]*std::max(mergedUncertaintiesUp[sample.variables[0]], mergedUncertaintiesDown[sample.variables[0]]),
        mergedUncertaintiesUpSign[sample.variables[1]]*std::max(mergedUncertaintiesUp[sample.variables[1]], mergedUncertaintiesDown[sample.variables[1]]),
        mergedUncertaintiesUpSign[sample.variables[2]]*std::max(mergedUncertaintiesUp[sample.variables[2]], mergedUncertaintiesDown[sample.variables[2]]));
      myfile << buffer;
      
      for(auto& var : sample.variables) {
        totalUncertainties2[var] -= mergedUncertainties2[var];
        totalUncertainties2[var] += pow(std::max(mergedUncertaintiesUp[var], mergedUncertaintiesDown[var]),2);
      }
    }
  }
  
  std::cout << "\n### Total systematic uncertainty (after linear merges)" << std::endl;
  myfile    << "\\hline\n\\hline\nTotal";
  sprintf(buffer, "& %.2lf & %.3lf & %.2lf \\tabularnewline \n \\hline \n \\hline \n", sqrt(totalUncertainties2[sample.variables[0]]), sqrt(totalUncertainties2[sample.variables[1]]), sqrt(totalUncertainties2[sample.variables[2]]));
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


