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
  int channel = 3;

  int selectComparisonSet = 0; //0: full set of uncertainties for Moriond; 1: only uncertainties also determined for B-Regression study; 2: minimal set used for BReg-study

  std::string path;
  if(channel==3||channel==1)path = "/nfs/dust/test/cms/user/kirschen/BRegression_PE_NewCalibrationJan2014Applied_MCS3250MinEvtsBReg/";
  else path = "/nfs/dust/cms/user/mseidel/pseudoexperiments/topmass_131012/";

  path += lepton_[channel]; path += "/";

  sample.path = path;
  sample.crossSection = 230;
  sample.peLumi = 20000.;

  sample.variables = {"mass_mTop_JES", "JES_mTop_JES", "mass_mTop"};

  sample.ensembles["calibration"] = ensemble("", 0, {std::make_pair("mass_mTop_JES",std::make_pair(172.5,0.02)), std::make_pair("JES_mTop_JES",std::make_pair(1.0,0.0)), std::make_pair("mass_mTop",std::make_pair(172.5,0.01))});
  if(channel==3)sample.ensembles["default"] = ensemble("Summer12_TTJets1725_1.00/job_*_ensemble.root", 7000000./1.75);
  else sample.ensembles["default"] = ensemble("job_*_ensemble.root", 7000000./1.75);
  sample.ensembles["defaultNewWeights"] = ensemble("weight.combinedWeight/job_*_ensemble.root", 7000000./1.75);

  sample.ensembles["puUp"] = ensemble("weight.combinedWeight/weight.puWeight*weight.puWeightUp/job_*_ensemble.root", 7000000./1.75);
  sample.ensembles["puDown"] = ensemble("weight.combinedWeight/weight.puWeight*weight.puWeightDown/job_*_ensemble.root", 7000000./1.75);

  sample.ensembles["bFragLEP"] = ensemble("weight_frag/job_*_ensemble.root", 7000000./1.75);
  sample.ensembles["bFragLEPHard"] = ensemble("weight_fragHard/job_*_ensemble.root", 7000000./1.75);
  sample.ensembles["bFragLEPSoft"] = ensemble("weight_fragSoft/job_*_ensemble.root", 7000000./1.75);

  sample.ensembles["bFNuUp"] = ensemble("weight_fNuUp/job_*_ensemble.root", 7000000./1.75);
  sample.ensembles["bFNuDown"] = ensemble("weight_fNuDown/job_*_ensemble.root", 7000000./1.75);

  sample.ensembles["flavorQCDUp"] = ensemble("Summer12_TTJets1725_flavor:up_FlavorQCD/job_*_ensemble.root", 7000000./1.75);
  sample.ensembles["flavorQCDDown"] = ensemble("Summer12_TTJets1725_flavor:down_FlavorQCD/job_*_ensemble.root", 7000000./1.75);

  sample.ensembles["flavorBUp"] = ensemble("Summer12_TTJets1725_flavor:up_FlavorPureBottom/job_*_ensemble.root", 7000000./1.75);
  sample.ensembles["flavorBDown"] = ensemble("Summer12_TTJets1725_flavor:down_FlavorPureBottom/job_*_ensemble.root", 7000000./1.75);

  sample.ensembles["flavorQUp"] = ensemble("Summer12_TTJets1725_flavor:up_FlavorPureQuark/job_*_ensemble.root", 7000000./1.75);
  sample.ensembles["flavorQDown"] = ensemble("Summer12_TTJets1725_flavor:down_FlavorPureQuark/job_*_ensemble.root", 7000000./1.75);

  sample.ensembles["flavorGUp"] = ensemble("Summer12_TTJets1725_flavor:up_FlavorPureGluon/job_*_ensemble.root", 7000000./1.75);
  sample.ensembles["flavorGDown"] = ensemble("Summer12_TTJets1725_flavor:down_FlavorPureGluon/job_*_ensemble.root", 7000000./1.75);

  sample.ensembles["bTagSFUp"] = ensemble("weight_bTagSFUp/job_*_ensemble.root", 7000000./1.75);
  sample.ensembles["bTagSFDown"] = ensemble("weight_bTagSFDown/job_*_ensemble.root", 7000000./1.75);
  sample.ensembles["misTagSFUp"] = ensemble("weight_misTagSFUp/job_*_ensemble.root", 7000000./1.75);
  sample.ensembles["misTagSFDown"] = ensemble("weight_misTagSFDown/job_*_ensemble.root", 7000000./1.75);

  sample.ensembles["topPt"] = ensemble("weight_topPt/job_*_ensemble.root", 7000000./1.75);

  sample.ensembles["fSig1.0"] = ensemble("fsig_1.00/job_*_ensemble.root", 7000000./1.75);
  sample.ensembles["fSig0.9"] = ensemble("fsig_0.90/job_*_ensemble.root", 7000000./1.75);

  sample.ensembles["matchingUp"] = ensemble("Summer12_TTJets1725_matchingup/job_*_ensemble.root", 5415010./1.75);
  sample.ensembles["matchingDown"] = ensemble("Summer12_TTJets1725_matchingdown/job_*_ensemble.root", 5476728./1.75);

  sample.ensembles["scaleUp"] = ensemble("Summer12_TTJets1725_scaleup/job_*_ensemble.root", 5009488./1.75);
  sample.ensembles["scaleDown"] = ensemble("Summer12_TTJets1725_scaledown/job_*_ensemble.root", 5387181./1.75);

  sample.ensembles["jerUp"] = ensemble("Summer12_TTJets1725_jer:up/job_*_ensemble.root", 7000000./1.75);
  sample.ensembles["jerDown"] = ensemble("Summer12_TTJets1725_jer:down/job_*_ensemble.root", 7000000./1.75);

  sample.ensembles["jesUp"] = ensemble("Summer12_TTJets1725_source:up_Total/job_*_ensemble.root", 7000000./1.75);
  sample.ensembles["jesDown"] = ensemble("Summer12_TTJets1725_source:down_Total/job_*_ensemble.root", 7000000./1.75);

  sample.ensembles["jesPuUp"] = ensemble("Summer12_TTJets1725_source:up_PileUpPtBB/job_*_ensemble.root", 7000000./1.75);
  sample.ensembles["jesPuDown"] = ensemble("Summer12_TTJets1725_source:down_PileUpPtBB/job_*_ensemble.root", 7000000./1.75);

  sample.ensembles["mcatnlo"] = ensemble("Summer12_TTJets1725_mcatnlo_herwig/job_*_ensemble.root", 32852589./1.75);
  sample.ensembles["powheg"] = ensemble("Summer12_TTJets1725_powheg/job_*_ensemble.root", 21675970./1.75);
  sample.ensembles["powhegHerwig"] = ensemble("Summer12_TTJets1725_powheg_herwig/job_*_ensemble.root", 27684235./1.75);

  sample.ensembles["defaultSC"] = ensemble("Summer12_TTJets1725_MGDecays/job_*_ensemble.root", 56000000./1.75);
  sample.ensembles["P11"] = ensemble("Summer12_TTJets1725_MGDecays_P11/job_*_ensemble.root", 27000000./1.75);
  sample.ensembles["P11noCR"] = ensemble("Summer12_TTJets1725_MGDecays_P11noCR/job_*_ensemble.root", 27000000./1.75);
  sample.ensembles["P11mpiHi"] = ensemble("Summer12_TTJets1725_MGDecays_P11mpiHi/job_*_ensemble.root", 18000000./1.75);
  sample.ensembles["P11TeV"] = ensemble("Summer12_TTJets1725_MGDecays_P11TeV/job_*_ensemble.root", 18000000./1.75);

  sample.ensembles["eesUp"] = ensemble("Summer12_TTJets1725_eesShift_+1/job_*_ensemble.root", 7000000./1.75);
  sample.ensembles["eesDown"] = ensemble("Summer12_TTJets1725_eesShift_-1/job_*_ensemble.root", 7000000./1.75);

  sample.ensembles["mesUp"] = ensemble("Summer12_TTJets1725_mesShift_+1/job_*_ensemble.root", 7000000./1.75);
  sample.ensembles["mesDown"] = ensemble("Summer12_TTJets1725_mesShift_-1/job_*_ensemble.root", 7000000./1.75);

  sample.ensembles["metcl2"] = ensemble("Summer12_TTJets1725_metcl_2/job_*_ensemble.root", 7000000./1.75);

  ///////////////////////////////////


  switch( selectComparisonSet )
    {
    case 0: // uncertainties not (yet) determined or considered for b-regression study
      sample.comparisons["Jet energy response (b)          "] = comparison("default", "flavorBUp", "flavorBDown", true); // not meaningful for b-regression study
      sample.comparisons["Jet energy response (FlavorQCD)  "] = comparison("default", "flavorQCDUp", "flavorQCDDown", true, false);// not meaningful for b-regression study
      sample.comparisons["\\pt- and $\\eta$-dependent JES  "] = comparison("default", "jesUp", "jesDown"); // no change expected for b-regression study
      sample.comparisons["Jet energy resolution            "] = comparison("default", "jerUp", "jerDown"); // no change expected for b-regression study
      sample.comparisons["Pile-up (JES)                    "] = comparison("defaultNewWeights", "jesPuUp", "jesPuDown");
      sample.comparisons["Non-\\ttbar background           "] = comparison("defaultNewWeights", "fSig0.9", "fSig1.0");
      sample.comparisons["Lepton energy scale (electron)   "] = comparison("defaultNewWeights", "eesUp", "eesDown");
      sample.comparisons["Lepton energy scale (muon)       "] = comparison("defaultNewWeights", "mesUp", "mesDown");
      sample.comparisons["MET                              "] = comparison("defaultNewWeights", "metcl2");
      sample.comparisons["Spin correlations                "] = comparison("calibration", "defaultSC", "", false, true); // not relevant anymore when using MadSpin samples consistently, not considered as common for both anymore
      
    case 1: // uncertainties determined for both variants
      sample.comparisons["Pile-up (pp cross-section)       "] = comparison("defaultNewWeights", "puUp", "puDown");
      sample.comparisons["Jet energy response (udsc)       "] = comparison("defaultNewWeights", "flavorQUp", "flavorQDown", true);//
      sample.comparisons["Jet energy response (gluon)      "] = comparison("defaultNewWeights", "flavorGUp", "flavorGDown", true);//
      sample.comparisons["b fragmentation                  "] = comparison("defaultNewWeights", "bFragLEP");
      sample.comparisons["Semi-leptonic B hadron decays    "] = comparison("defaultNewWeights", "bFNuUp", "bFNuDown");
      sample.comparisons["b-tag rate                       "] = comparison("defaultNewWeights", "bTagSFUp", "bTagSFDown");//
      sample.comparisons["b-tag (mistag rate)              "] = comparison("defaultNewWeights", "misTagSFUp", "misTagSFDown");//
      sample.comparisons["Top-pt reweighting               "] = comparison("defaultNewWeights", "topPt", "");//
      sample.comparisons["ME-PS matching threshold         "] = comparison("calibration", "matchingUp", "matchingDown", false);
      sample.comparisons["$Q^{2}$ scale                    "] = comparison("calibration", "scaleUp", "scaleDown", false);
      sample.comparisons["MadGraph (no SC) vs. Powheg      "] = comparison("calibration", "powheg", "", false, false);
      sample.comparisons["Powheg+Pythia6 vs. MC@NLO+Herwig6"] = comparison("powheg", "mcatnlo", "", false, false);
      sample.comparisons["Powheg+Pythia6 vs. Powheg+Herwig6"] = comparison("powheg", "powhegHerwig", "", false, false);
      sample.comparisons["MC@NLO+Herwig6 vs. Powheg+Herwig6"] = comparison("mcatnlo", "powhegHerwig", "", false, false);
      sample.comparisons["ME generator                     "] = comparison("defaultSC", "powheg", "", false);
      sample.comparisons["Pythia Z2* vs. P11               "] = comparison("defaultSC", "P11", "", false, false);
      sample.comparisons["Color reconnection               "] = comparison("P11", "P11noCR", "", false);
      sample.comparisons["Underlying event                 "] = comparison("P11", "P11mpiHi", "P11TeV", false);
      break;
      
    case 3: // minimal set of comparison to determine differences with/ without b-regression
      sample.comparisons["Powheg+Pythia6 vs. Powheg+Herwig6"] = comparison("powheg", "powhegHerwig", "", false);
      sample.comparisons["b fragmentation                  "] = comparison("defaultNewWeights", "bFragLEP");
      sample.comparisons["Semi-leptonic B hadron decays    "] = comparison("defaultNewWeights", "bFNuUp", "bFNuDown");
      break;
    }

    


}

void SystematicUncertainties::fillAllJets()
{
  //sample.path = "/nfs/dust/cms/user/eschliec/TopMass/topmass_140317_1201/";
  //sample.path = "/nfs/dust/cms/user/eschliec/TopMass/topmass_140401_1201/";
  //sample.path = "/nfs/dust/cms/user/eschliec/TopMass/topmass_140418_1201a/";
  //sample.path = "/nfs/dust/cms/user/eschliec/TopMass/topmass_140520_1801/";
  //sample.path = "/nfs/dust/cms/user/eschliec/TopMass/topmass_140911_1201/";
  //sample.path = "/nfs/dust/cms/user/eschliec/TopMass/topmass_140916_1201/";
  //sample.path = "/nfs/dust/cms/user/eschliec/TopMass/topmass_140917_1701/";
  //sample.path = "/nfs/dust/cms/user/eschliec/TopMass/topmass_141006_1401/";
  //sample.path = "/nfs/dust/cms/user/eschliec/TopMass/topmass_140918_1601/";
  //sample.path = "/nfs/dust/cms/user/eschliec/TopMass/topmass_141006_1701/";
  //sample.path = "/nfs/dust/cms/user/eschliec/TopMass/topmass_141006_2001/";
  //sample.path = "/nfs/dust/cms/user/eschliec/TopMass/topmass_140918_1401/";
  //sample.path = "/nfs/dust/cms/user/eschliec/TopMass/topmass_141007_1301/";
  //sample.path = "/nfs/dust/cms/user/eschliec/TopMass/topmass_141007_1501/";
  sample.path = "/nfs/dust/cms/user/eschliec/TopMass/topmass_141017_1601/";
  sample.crossSection = 245.794;
  sample.peLumi = 18192.;

  //sample.variables = {"mass_mTop_JES", "JES_mTop_JES", "mass_mTop"};
  //sample.variables = {"mass_mTop_JES_fSig", "JES_mTop_JES_fSig", "mass_mTop_fSig"};
  //sample.variables = {"mass_mTop_JES_fSig_fCP", "JES_mTop_JES_fSig_fCP", "mass_mTop"};
  sample.variables = {"mass_mTop_JES_fSig_fCP", "JES_mTop_JES_fSig_fCP", "mass_mTop_fSig_fCP"};  // DEFAULT
  //sample.variables = {"mass_mTop_JES_fSig_fCP", "JES_mTop_JES_fSig_fCP", "mass_mTop_fSig"};
  //sample.variables = {"mass_mTop_JES_fSig_fCP", "mass_mTop_JES_fSig", "mass_mTop_JES"};
  //sample.variables = {"fCP_mTop_JES_fSig_fCP", "fCP_mTop_JES_fCP", "fCP_mTop_fCP"};
  //sample.variables = {"fSig_mTop_JES_fSig_fCP", "fSig_mTop_JES_fSig", "fSig_mTop_fSig"};
  //sample.variables = {"fCP_mTop_JES_fSig_fCP", "JES_mTop_JES_fSig_fCP", "fCP_mTop_fSig_fCP"};
  //sample.variables = {"fSig_mTop_JES_fSig_fCP", "JES_mTop_JES_fSig_fCP", "fSig_mTop_fSig_fCP"};

  sample.ensembles["calibration"] = ensemble("", 0, {std::make_pair(sample.variables[0],std::make_pair(172.5,0.0544808)), std::make_pair(sample.variables[1],std::make_pair(1.0,0.00046571)), std::make_pair(sample.variables[2],std::make_pair(172.5,0.0542347))});
  sample.ensembles["calibrationDummy"] = ensemble("", 0, {std::make_pair(sample.variables[0],std::make_pair(172.5,0.)), std::make_pair(sample.variables[1],std::make_pair(1.0,0.)), std::make_pair(sample.variables[2],std::make_pair(172.5,0.))});
  sample.ensembles["default"] = ensemble("Z2_S12_ABS_JES_100_172_5/job_*_ensemble.root", 62131965./1.8);

  sample.ensembles["puUp"  ] = ensemble("Z2_S12_PU_Up/job_*_ensemble.root"  , 62131965./1.8);
  sample.ensembles["puDown"] = ensemble("Z2_S12_PU_Down/job_*_ensemble.root", 62131965./1.8);

  sample.ensembles["bFragLEP"    ] = ensemble("Z2_S12_FRAG/job_*_ensemble.root"     , 62131965./1.8);
  sample.ensembles["bFragLEPHard"] = ensemble("Z2_S12_FRAG_Hard/job_*_ensemble.root", 62131965./1.8);
  sample.ensembles["bFragLEPSoft"] = ensemble("Z2_S12_FRAG_Soft/job_*_ensemble.root", 62131965./1.8);

  sample.ensembles["bFNuUp"  ] = ensemble("Z2_S12_NU_Up/job_*_ensemble.root"  , 62131965./1.8);
  sample.ensembles["bFNuDown"] = ensemble("Z2_S12_NU_Down/job_*_ensemble.root", 62131965./1.8);

  sample.ensembles["flavorBUp"  ] = ensemble("Z2_S12_BJES_Up/job_*_ensemble.root"  , 62131965./1.8);
  sample.ensembles["flavorBDown"] = ensemble("Z2_S12_BJES_Down/job_*_ensemble.root", 62131965./1.8);

  sample.ensembles["flavorQUp"  ] = ensemble("Z2_S12_LJES_Up/job_*_ensemble.root"  , 62131965./1.8);
  sample.ensembles["flavorQDown"] = ensemble("Z2_S12_LJES_Down/job_*_ensemble.root", 62131965./1.8);

  sample.ensembles["flavorGUp"  ] = ensemble("Z2_S12_GJES_Up/job_*_ensemble.root"  , 62131965./1.8);
  sample.ensembles["flavorGDown"] = ensemble("Z2_S12_GJES_Down/job_*_ensemble.root", 62131965./1.8);

  sample.ensembles["bTagSFUp"    ] = ensemble("Z2_S12_BTAG_Up/job_*_ensemble.root"  , 62131965./1.8);
  sample.ensembles["bTagSFDown"  ] = ensemble("Z2_S12_BTAG_Down/job_*_ensemble.root", 62131965./1.8);
  sample.ensembles["misTagSFUp"  ] = ensemble("Z2_S12_MTAG_Up/job_*_ensemble.root"  , 62131965./1.8);
  sample.ensembles["misTagSFDown"] = ensemble("Z2_S12_MTAG_Down/job_*_ensemble.root", 62131965./1.8);

  sample.ensembles["topPt"] = ensemble("Z2_S12_TopPt/job_*_ensemble.root", 62131965./1.8);

  sample.ensembles["fSigUp"  ] = ensemble("Z2_S12_FSIG_Up/job_*_ensemble.root"  , 62131965./1.8);
  sample.ensembles["fSigDown"] = ensemble("Z2_S12_FSIG_Down/job_*_ensemble.root", 62131965./1.8);

  //sample.ensembles["shape"   ] = ensemble("Z2_S12_BackgroundSys/job_*_ensemble.root", 62131965./1.8);
  //sample.ensembles["shapeAlt"] = ensemble("../topmass_140401_1201a/Z2_S12_BackgroundSys2/job_*_ensemble.root", 62131965./1.8);
  //sample.ensembles["bkg0x"   ] = ensemble("../topmass_140401_1201d/Z2_S12_ABS_JES_100_172_5/job_*_ensemble.root", 62131965./1.8);
  //sample.ensembles["bkg2x"   ] = ensemble("../topmass_140401_1201e/Z2_S12_BKG_2X/job_*_ensemble.root", 62131965./1.8);
  //sample.ensembles["shapeSid"] = ensemble("../topmass_140520_1805/Z2_S12_ABS_JES_100_172_5/job_*_ensemble.root", 62131965./1.8);

  //sample.ensembles["shapeDef"  ] = ensemble("../topmass_140520_1808/Z2_S12_ABS_JES_100_172_5/job_*_ensemble.root", 62131965./1.8); // identical to default, as hoped
  //sample.ensembles["shapeTUp"  ] = ensemble("../topmass_140520_1808/Z2_S12_BackgroundShapeT_Up/job_*_ensemble.root", 62131965./1.8);
  //sample.ensembles["shapeTDown"] = ensemble("../topmass_140520_1808/Z2_S12_BackgroundShapeT_Down/job_*_ensemble.root", 62131965./1.8);
  //sample.ensembles["shapeWUp"  ] = ensemble("../topmass_140520_1808/Z2_S12_BackgroundShapeW_Up/job_*_ensemble.root", 62131965./1.8);
  //sample.ensembles["shapeWDown"] = ensemble("../topmass_140520_1808/Z2_S12_BackgroundShapeW_Down/job_*_ensemble.root", 62131965./1.8);

  //sample.ensembles["shapeDef2" ] = ensemble("../topmass_140520_1809/Z2_S12_ABS_JES_100_172_5/job_*_ensemble.root", 62131965./1.8); // identical to default, as hoped
  //sample.ensembles["shapeVar"  ] = ensemble("../topmass_140520_1809/Z2_S12_BackgroundShapeSys/job_*_ensemble.root", 62131965./1.8);
  //sample.ensembles["shapeTEST1"] = ensemble("../topmass_140520_1809a/Z2_S12_BackgroundShapeSys1/job_*_ensemble.root", 62131965./1.8);
  //sample.ensembles["shapeTEST2"] = ensemble("../topmass_140520_1809a/Z2_S12_BackgroundShapeSys2/job_*_ensemble.root", 62131965./1.8);

  //sample.ensembles["shape0"] = ensemble("../topmass_140520_1814/Z2_S12_BackgroundShapeSys0/job_*_ensemble.root", 62131965./1.8);
  //sample.ensembles["shape1"] = ensemble("../topmass_140520_1814/Z2_S12_BackgroundShapeSys1/job_*_ensemble.root", 62131965./1.8);
  //sample.ensembles["shape2"] = ensemble("../topmass_140520_1814/Z2_S12_BackgroundShapeSys2/job_*_ensemble.root", 62131965./1.8);
  //sample.ensembles["shape3"] = ensemble("../topmass_140520_1814/Z2_S12_BackgroundShapeSys3/job_*_ensemble.root", 62131965./1.8);
  //sample.ensembles["shape4"] = ensemble("../topmass_140520_1814/Z2_S12_BackgroundShapeSys4/job_*_ensemble.root", 62131965./1.8);
  //sample.ensembles["shape5"] = ensemble("../topmass_140520_1814/Z2_S12_BackgroundShapeSys5/job_*_ensemble.root", 62131965./1.8);
  //sample.ensembles["shape6"] = ensemble("../topmass_140520_1814/Z2_S12_BackgroundShapeSys6/job_*_ensemble.root", 62131965./1.8);
  //sample.ensembles["shape7"] = ensemble("../topmass_140520_1814/Z2_S12_BackgroundShapeSys7/job_*_ensemble.root", 62131965./1.8);
  //sample.ensembles["shape8"] = ensemble("../topmass_140520_1814/Z2_S12_BackgroundShapeSys8/job_*_ensemble.root", 62131965./1.8);

  sample.ensembles["matchingUp"  ] = ensemble("Z2_S12_Matching_Up/job_*_ensemble.root"  , 37083003./1.8);
  sample.ensembles["matchingDown"] = ensemble("Z2_S12_Matching_Down/job_*_ensemble.root", 34053113./1.8);

  sample.ensembles["scaleUp"  ] = ensemble("Z2_S12_Scale_Up/job_*_ensemble.root"  , 41908271./1.8);
  sample.ensembles["scaleDown"] = ensemble("Z2_S12_Scale_Down/job_*_ensemble.root", 39286663./1.8);

  sample.ensembles["jerUp"  ] = ensemble("Z2_S12_JER_Up/job_*_ensemble.root"  , 62131965./1.8);
  sample.ensembles["jerDown"] = ensemble("Z2_S12_JER_Down/job_*_ensemble.root", 62131965./1.8);

  sample.ensembles["jer3Up"     ] = ensemble("Z2_S12_NEWJER3_Up/job_*_ensemble.root"     , 62131965./1.8);
  sample.ensembles["jer3Central"] = ensemble("Z2_S12_NEWJER3_Central/job_*_ensemble.root", 62131965./1.8);
  sample.ensembles["jer3Down"   ] = ensemble("Z2_S12_NEWJER3_Down/job_*_ensemble.root"   , 62131965./1.8);

  sample.ensembles["jes1Up"  ] = ensemble("Z2_S12_COR1JES_Up/job_*_ensemble.root"  , 62131965./1.8);
  sample.ensembles["jes1Down"] = ensemble("Z2_S12_COR1JES_Down/job_*_ensemble.root", 62131965./1.8);
  sample.ensembles["jes2Up"  ] = ensemble("Z2_S12_COR2JES_Up/job_*_ensemble.root"  , 62131965./1.8);
  sample.ensembles["jes2Down"] = ensemble("Z2_S12_COR2JES_Down/job_*_ensemble.root", 62131965./1.8);
  sample.ensembles["jes3Up"  ] = ensemble("Z2_S12_COR3JES_Up/job_*_ensemble.root"  , 62131965./1.8);
  sample.ensembles["jes3Down"] = ensemble("Z2_S12_COR3JES_Down/job_*_ensemble.root", 62131965./1.8);
  sample.ensembles["jes4Up"  ] = ensemble("Z2_S12_COR4JES_Up/job_*_ensemble.root"  , 62131965./1.8);
  sample.ensembles["jes4Down"] = ensemble("Z2_S12_COR4JES_Down/job_*_ensemble.root", 62131965./1.8);

  //sample.ensembles["triggerUp"   ] = ensemble("Z2_S12_TRIGGERJES_Up/job_*_ensemble.root"  , 62131965./1.8);
  //sample.ensembles["triggerDown" ] = ensemble("Z2_S12_TRIGGERJES_Down/job_*_ensemble.root", 62131965./1.8);
  //sample.ensembles["trigger2Up"  ] = ensemble("Z2_S12_TRIGGERJES2_Up/job_*_ensemble.root"  , 62131965./1.8);
  //sample.ensembles["trigger2Down"] = ensemble("Z2_S12_TRIGGERJES2_Down/job_*_ensemble.root", 62131965./1.8);
  //sample.ensembles["trigger3Up"  ] = ensemble("Z2_S12_TRIGGERJES3_Up/job_*_ensemble.root"  , 62131965./1.8);
  //sample.ensembles["trigger3Down"] = ensemble("Z2_S12_TRIGGERJES3_Down/job_*_ensemble.root", 62131965./1.8);
  //sample.ensembles["trigger4Up"  ] = ensemble("Z2_S12_TRIGGERJES4_Up/job_*_ensemble.root"  , 62131965./1.8);
  //sample.ensembles["trigger4Down"] = ensemble("Z2_S12_TRIGGERJES4_Down/job_*_ensemble.root", 62131965./1.8);
  //sample.ensembles["triggerKosta"] = ensemble("../topmass_140520_1802/Z2_S12_ABS_JES_100_172_5/job_*_ensemble.root", 62131965./1.8);
  //sample.ensembles["triggerTEST" ] = ensemble("../topmass_140520_1803/Z2_S12_ABS_JES_100_172_5/job_*_ensemble.root", 62131965./1.8);
  //sample.ensembles["triggerTEST2"] = ensemble("../topmass_140520_1804/Z2_S12_ABS_JES_100_172_5/job_*_ensemble.root", 62131965./1.8);
  //sample.ensembles["triggerTEST3"] = ensemble("../topmass_140520_1806/Z2_S12_ABS_JES_100_172_5/job_*_ensemble.root", 62131965./1.8); //should be identical to triggerTEST2 and/or default !!!!
  //sample.ensembles["triggerTEST4"] = ensemble("../topmass_140520_1807/Z2_S12_ABS_JES_100_172_5/job_*_ensemble.root", 62131965./1.8); //identical to triggerTEST3 !!!!

  sample.ensembles["jesPuBBUp"  ] = ensemble("Z2_S12_PUPTBBJES_Up/job_*_ensemble.root"  , 62131965./1.8);
  sample.ensembles["jesPuBBDown"] = ensemble("Z2_S12_PUPTBBJES_Down/job_*_ensemble.root", 62131965./1.8);
  sample.ensembles["jesPuECUp"  ] = ensemble("Z2_S12_PUPTECJES_Up/job_*_ensemble.root"  , 62131965./1.8);
  sample.ensembles["jesPuECDown"] = ensemble("Z2_S12_PUPTECJES_Down/job_*_ensemble.root", 62131965./1.8);

  sample.ensembles["mcatnlo"     ] = ensemble("Z2_S12_MCNLO/job_*_ensemble.root" , 32852589./1.8);
  sample.ensembles["powheg"      ] = ensemble("Z2_S12_POWHEG/job_*_ensemble.root", 21675970./1.8);
  sample.ensembles["powhegHerwig"] = ensemble("Z2_S12_POWHER/job_*_ensemble.root", 27684235./1.8);

  sample.ensembles["defaultSC"] = ensemble("Z2_S12_Z2/job_*_ensemble.root"      , 41761265./6./6.*9.*9./1.8);
  sample.ensembles["P11"      ] = ensemble("Z2_S12_P11/job_*_ensemble.root"     , 11651739./6./6.*9.*9./1.8);
  sample.ensembles["P11noCR"  ] = ensemble("Z2_S12_P11NoCR/job_*_ensemble.root" , 11919063./6./6.*9.*9./1.8);
  sample.ensembles["P11mpiHi" ] = ensemble("Z2_S12_P11mpiHi/job_*_ensemble.root",  7953758./6./6.*9.*9./1.8);
  sample.ensembles["P11TeV"   ] = ensemble("Z2_S12_P11TeV/job_*_ensemble.root"  ,  7946264./6./6.*9.*9./1.8);

  //sample.ensembles["RD"] = ensemble("../topmass_140401_1201b/Z2_S12_RD1/job_*_ensemble.root", 31228390./6./6.*9.*9.);

  //sample.ensembles["PDFDown"] = ensemble("", 0, {std::make_pair(sample.variables[0],std::make_pair(172.322,0.0)), std::make_pair(sample.variables[1],std::make_pair(1.0007,0.0)), std::make_pair(sample.variables[2],std::make_pair(172.337,0.0))});
  //sample.ensembles["PDFUp"  ] = ensemble("", 0, {std::make_pair(sample.variables[0],std::make_pair(172.347,0.0)), std::make_pair(sample.variables[1],std::make_pair(1.0010,0.0)), std::make_pair(sample.variables[2],std::make_pair(172.362,0.0))});

  ///////////////////////////////////

  sample.comparisons["Calibration                      "] = comparison("calibration", "calibrationDummy", "", false);
  sample.comparisons["Pile-up (pp cross-section)       "] = comparison("default", "puUp", "puDown", true);
  sample.comparisons["Jet energy response (udsc)       "] = comparison("default", "flavorQUp", "flavorQDown", true);
  sample.comparisons["Jet energy response (gluon)      "] = comparison("default", "flavorGUp", "flavorGDown", true);
  sample.comparisons["Jet energy response (b)          "] = comparison("default", "flavorBUp", "flavorBDown", true);
  //sample.comparisons["Jet energy response (FlavorQCD)  "] = comparison("default", "flavorQCDUp", "flavorQCDDown", true, false);
  sample.comparisons["b fragmentation                  "] = comparison("default", "bFragLEP", "", true);
  sample.comparisons["Semi-leptonic B hadron decays    "] = comparison("default", "bFNuUp", "bFNuDown", true);
  sample.comparisons["b-tag rate                       "] = comparison("default", "bTagSFUp", "bTagSFDown", true);
  sample.comparisons["b-tag (mistag rate)              "] = comparison("default", "misTagSFUp", "misTagSFDown", true);
  //sample.comparisons["Trigger                          "] = comparison("default", "triggerUp", "triggerDown", true, false);
  //sample.comparisons["Trigger2                         "] = comparison("default", "trigger2Up", "trigger2Down", true);
  //sample.comparisons["Trigger3                         "] = comparison("default", "trigger3Up", "trigger3Down", true, false);
  //sample.comparisons["Trigger4                         "] = comparison("default", "trigger4Up", "trigger4Down", true, false);
  //sample.comparisons["TriggerKostas                    "] = comparison("default", "triggerKosta", "", true, false);
  //sample.comparisons["TriggerTEST                      "] = comparison("default", "triggerTEST", "", true, false);
  //sample.comparisons["TriggerTEST2                     "] = comparison("default", "triggerTEST2", "", true, false);
  //sample.comparisons["TriggerTEST Bug?                 "] = comparison("default", "triggerTEST3", "triggerTEST4", true, false);
  sample.comparisons["Top-\\pt reweighting             "] = comparison("default", "topPt", "", true);
  sample.comparisons["ME-PS matching threshold         "] = comparison("default", "matchingUp", "matchingDown", false);
  sample.comparisons["$Q^{2}$ scale                    "] = comparison("default", "scaleUp", "scaleDown", false);
  sample.comparisons["Jet energy resolution            "] = comparison("default", "jerUp" , "jerDown" , true);
  sample.comparisons["Jet energy resolution NEW        "] = comparison("jer3Central", "jer3Up" , "jer3Down" , true, false);
  sample.comparisons["JES MPFInSitu                    "] = comparison("default", "jes1Up", "jes1Down", true);
  sample.comparisons["JES Flavor                       "] = comparison("default", "jes2Up", "jes2Down", true);
  sample.comparisons["JES Intercalibration             "] = comparison("default", "jes3Up", "jes3Down", true);
  sample.comparisons["JES Uncorrelated                 "] = comparison("default", "jes4Up", "jes4Down", true);
  sample.comparisons["Pile-up (JES) BB                 "] = comparison("default", "jesPuBBUp", "jesPuBBDown", true);
  sample.comparisons["Pile-up (JES) EC                 "] = comparison("default", "jesPuECUp", "jesPuECDown", true);
  sample.comparisons["MadGraph (no SC) vs. Powheg      "] = comparison("default", "powheg", "", false, false);
  sample.comparisons["MadGraph vs. MC@NLO+Herwig6      "] = comparison("default", "mcatnlo", "", false, false);
  sample.comparisons["Powheg+Pythia6 vs. MC@NLO+Herwig6"] = comparison("powheg", "mcatnlo", "", false, false);
  sample.comparisons["Powheg+Pythia6 vs. Powheg+Herwig6"] = comparison("powheg", "powhegHerwig", "", false, false);
  sample.comparisons["MC@NLO+Herwig6 vs. Powheg+Herwig6"] = comparison("mcatnlo", "powhegHerwig", "", false, false);
  sample.comparisons["ME generator                     "] = comparison("defaultSC", "powheg", "", false);
  sample.comparisons["Spin correlations                "] = comparison("default", "defaultSC", "", false, false);
  sample.comparisons["Pythia Z2* vs. P11               "] = comparison("defaultSC", "P11", "", false, false);
  sample.comparisons["Color reconnection               "] = comparison("P11", "P11noCR", "", false);
  sample.comparisons["Underlying event                 "] = comparison("P11", "P11mpiHi", "P11TeV", false);
  sample.comparisons["Non-\\ttbar background \\fsig    "] = comparison("default", "fSigUp", "fSigDown", true);
  //sample.comparisons["Non-\\ttbar background \\fsig 2  "] = comparison("default", "bkg0x", "bkg2x", true, false);
  //sample.comparisons["Non-\\ttbar background shape     "] = comparison("default", "shape", "", true, false);
  //sample.comparisons["Non-\\ttbar background shapeAlt  "] = comparison("default", "shapeAlt", "", true, false);
  //sample.comparisons["Non-\\ttbar background shapeSide "] = comparison("default", "shapeSid", "", true, false);
  //sample.comparisons["Non-\\ttbar background shapeT    "] = comparison("default", "shapeTUp", "shapeTDown", true, false);
  //sample.comparisons["Non-\\ttbar background shapeW    "] = comparison("default", "shapeWUp", "shapeWDown", true, false);
  //sample.comparisons["Non-\\ttbar background shapeVar  "] = comparison("default", "shapeVar", "", true, false);
  //sample.comparisons["Non-\\ttbar background shapeTEST "] = comparison("shapeTEST1", "shapeTEST2", "", true, false);
  //sample.comparisons["Non-\\ttbar background shape1    "] = comparison("shape0", "shape1", "shape2", true, false);
  //sample.comparisons["Non-\\ttbar background shape2    "] = comparison("shape0", "shape3", "shape4", true, false);
  //sample.comparisons["Non-\\ttbar background shape3    "] = comparison("shape0", "shape5", "shape6", true, true);
  //sample.comparisons["Non-\\ttbar background shape4    "] = comparison("shape0", "shape7", "shape8", true, true);
  //sample.comparisons["Run dependent                    "] = comparison("defaultSC", "RD", "", false, false);
  //sample.comparisons["PDF                              "] = comparison("default", "PDFDown", "PDFUp", true, true);

  sample.mergedcomparisons["JEC"] = mergedcomparison(
      {"JES Flavor                       ",
       "JES Intercalibration             ",
       "JES MPFInSitu                    ",
       "JES Uncorrelated                 "}
  );
  sample.mergedcomparisons["JEC Flavor"] = mergedcomparison(
      {"Jet energy response (udsc)       ",
       "Jet energy response (gluon)      ",
       "Jet energy response (b)          "}
  );
  sample.mergedcomparisons["JES PileUpPt"] = mergedcomparison(
      {"Pile-up (JES) BB                 ",
       "Pile-up (JES) EC                 "}
  );
  sample.mergedcomparisons["PileUp"] = mergedcomparison(
      {"Pile-up (JES) BB                 ",
       "Pile-up (JES) EC                 ",
       "Pile-up (pp cross-section)       "}
  );
  sample.mergedcomparisons["b-tagging"] = mergedcomparison(
      {"b-tag rate                       ",
       "b-tag (mistag rate)              "}
  );
  //sample.mergedcomparisons["Non-\\ttbar background shape"] = mergedcomparison(
  //    {"Non-\\ttbar background shape1    ",
  //     "Non-\\ttbar background shape2    "}
  //);
  //sample.mergedcomparisons["Non-\\ttbar background"] = mergedcomparison(
  //    {"Non-\\ttbar background \\fsig    ",
  //     "Non-\\ttbar background shape1    ",
  //     "Non-\\ttbar background shape2    "}
  //);

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
    for(auto& var : sample.variables)
      it->second.shifts[var] = std::max(std::abs(nominal.values[var].first-up.values[var].first), std::abs(nominal.values[var].first-down.values[var].first));

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

  std::cout << "\n### Total uncertainties" << std::endl;
  myfile    << "\\hline\n\\hline\nTotal uncertainties";
  sprintf(buffer, "& %.2lf & %.3lf & %.2lf \n", sqrt(totalUncertainties2[sample.variables[0]]), sqrt(totalUncertainties2[sample.variables[1]]), sqrt(totalUncertainties2[sample.variables[2]]));
  printf("\t 2D mass = %.2lf GeV, JSF = %.3lf, 1D mass = %.2lf GeV\n", sqrt(totalUncertainties2[sample.variables[0]]), sqrt(totalUncertainties2[sample.variables[1]]), sqrt(totalUncertainties2[sample.variables[2]]));
  myfile << buffer;

  for(std::map<std::string, mergedcomparison>::const_iterator it = sample.mergedcomparisons.begin(); it != sample.mergedcomparisons.end(); ++it) {
    std::cout << it->first;
    myfile    << it->first;

    std::map<std::string, double> mergedUncertainties2;
    for(auto& var : sample.variables)
      mergedUncertainties2[var] = 0;

    for(std::vector<std::string>::const_iterator it2 = it->second.comparisons.begin(); it2 != it->second.comparisons.end(); ++it2){
      for(auto& var : sample.variables){
        mergedUncertainties2[var] += pow(sample.comparisons.find(*it2)->second.shifts[var],2);
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
  }


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


