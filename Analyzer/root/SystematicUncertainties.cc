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

  std::string path("/nfs/dust/cms/user/mseidel/pseudoexperiments/topmass_131012/");
  //std::string path("/nfs/dust/test/cms/user/kirschen/BRegression_PE_NewCalibrationJan2014Applied_MCS3250MinEvtsBReg/");

  path += lepton_[channel]; path += "/";

  sample.path = path;
  sample.crossSection = 230;
  sample.peLumi = 20000.;

  sample.variables = {"mass_mTop_JES", "JES_mTop_JES", "mass_mTop"};

  sample.ensembles["calibration"] = ensemble("", 0, {std::make_pair("mass_mTop_JES",std::make_pair(172.5,0.02)), std::make_pair("JES_mTop_JES",std::make_pair(1.0,0.0)), std::make_pair("mass_mTop",std::make_pair(172.5,0.01))});
  sample.ensembles["default"] = ensemble("job_*_ensemble.root", 7000000./1.75);
  //sample.ensembles["default"] = ensemble("Summer12_TTJets1725_1.00/job_*_ensemble.root", 7000000./1.75);
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

  sample.comparisons["Pile-up (pp cross-section)       "] = comparison("defaultNewWeights", "puUp", "puDown");
  sample.comparisons["Jet energy response (udsc)       "] = comparison("defaultNewWeights", "flavorQUp", "flavorQDown", true);//
  sample.comparisons["Jet energy response (gluon)      "] = comparison("defaultNewWeights", "flavorGUp", "flavorGDown", true);//
  sample.comparisons["Jet energy response (b)          "] = comparison("default", "flavorBUp", "flavorBDown", true);
  sample.comparisons["Jet energy response (FlavorQCD)  "] = comparison("default", "flavorQCDUp", "flavorQCDDown", true, false);
  sample.comparisons["b fragmentation                  "] = comparison("defaultNewWeights", "bFragLEP");
  sample.comparisons["Semi-leptonic B hadron decays    "] = comparison("defaultNewWeights", "bFNuUp", "bFNuDown");
  sample.comparisons["b-tag rate                       "] = comparison("defaultNewWeights", "bTagSFUp", "bTagSFDown");//
  sample.comparisons["b-tag (mistag rate)              "] = comparison("defaultNewWeights", "misTagSFUp", "misTagSFDown");//
  sample.comparisons["Top-pt reweighting               "] = comparison("defaultNewWeights", "topPt", "");//
  sample.comparisons["ME-PS matching threshold         "] = comparison("calibration", "matchingUp", "matchingDown", false);
  sample.comparisons["$Q^{2}$ scale                    "] = comparison("calibration", "scaleUp", "scaleDown", false);
  sample.comparisons["Jet energy resolution            "] = comparison("default", "jerUp", "jerDown");
  sample.comparisons["\\pt- and $\\eta$-dependent JES  "] = comparison("default", "jesUp", "jesDown");
  sample.comparisons["Pile-up (JES)                    "] = comparison("defaultNewWeights", "jesPuUp", "jesPuDown");
  sample.comparisons["MadGraph (no SC) vs. Powheg      "] = comparison("calibration", "powheg", "", false, false);
  sample.comparisons["Powheg+Pythia6 vs. MC@NLO+Herwig6"] = comparison("powheg", "mcatnlo", "", false, false);
  sample.comparisons["Powheg+Pythia6 vs. Powheg+Herwig6"] = comparison("powheg", "powhegHerwig", "", false, false);
  sample.comparisons["MC@NLO+Herwig6 vs. Powheg+Herwig6"] = comparison("mcatnlo", "powhegHerwig", "", false, false);
  sample.comparisons["ME generator                     "] = comparison("defaultSC", "powheg", "", false);
  sample.comparisons["Spin correlations                "] = comparison("calibration", "defaultSC", "", false, true);
  sample.comparisons["Pythia Z2* vs. P11               "] = comparison("defaultSC", "P11", "", false, false);
  sample.comparisons["Color reconnection               "] = comparison("P11", "P11noCR", "", false);
  sample.comparisons["Underlying event                 "] = comparison("P11", "P11mpiHi", "P11TeV", false);
  sample.comparisons["Non-\\ttbar background           "] = comparison("defaultNewWeights", "fSig0.9", "fSig1.0");
  sample.comparisons["Lepton energy scale (electron)   "] = comparison("defaultNewWeights", "eesUp", "eesDown");
  sample.comparisons["Lepton energy scale (muon)       "] = comparison("defaultNewWeights", "mesUp", "mesDown");
  sample.comparisons["MET                              "] = comparison("defaultNewWeights", "metcl2");
}

void SystematicUncertainties::fillAllJets()
{
  sample.path = "/nfs/dust/cms/user/eschliec/TopMass/";
  sample.crossSection = 245.794;
  sample.peLumi = 18192.;

  //sample.variables = {"mass_mTop_JES", "JES_mTop_JES", "mass_mTop"};
  sample.variables = {"mass_mTop_JES_fSig_fCP", "JES_mTop_JES_fSig_fCP", "mass_mTop"};

  sample.ensembles["calibration"] = ensemble("", 0, {std::make_pair(sample.variables[0],std::make_pair(172.5,0.0585859)), std::make_pair(sample.variables[1],std::make_pair(1.0,0.000489326)), std::make_pair(sample.variables[2],std::make_pair(172.5,0.058894))});
  sample.ensembles["default"] = ensemble("topmass_140305_1201/Z2_S12_ABS_JES_100_172_5_MadSpin_sig/MASS_172.5/JES_1.00/SHAPE_1.00/PERMU_0.0/FSIG_0.767995487/WEIGHT_weight.combinedWeight/job_*_ensemble.root", 62131965./1.8);

  sample.ensembles["puUp"  ] = ensemble("topmass_140305_1201a/Z2_S12_PU_Up/job_*_ensemble.root", 62131965./1.8);
  sample.ensembles["puDown"] = ensemble("topmass_140305_1201a/Z2_S12_PU_Down/job_*_ensemble.root", 62131965./1.8);

  sample.ensembles["bFragLEP"    ] = ensemble("topmass_140305_1201e/Z2_S12_FRAG/job_*_ensemble.root", 62131965./1.8);
  sample.ensembles["bFragLEPHard"] = ensemble("topmass_140305_1201e/Z2_S12_FRAG_Hard/job_*_ensemble.root", 62131965./1.8);
  sample.ensembles["bFragLEPSoft"] = ensemble("topmass_140305_1201e/Z2_S12_FRAG_Soft/job_*_ensemble.root", 62131965./1.8);

  sample.ensembles["bFNuUp"  ] = ensemble("topmass_140305_1201e/Z2_S12_NU_Up/job_*_ensemble.root", 62131965./1.8);
  sample.ensembles["bFNuDown"] = ensemble("topmass_140305_1201e/Z2_S12_NU_Down/job_*_ensemble.root", 62131965./1.8);

  //sample.ensembles["flavorQCDUp"] = ensemble("Summer12_TTJets1725_flavor:up_FlavorQCD/job_*_ensemble.root", 62131965./1.8);
  //sample.ensembles["flavorQCDDown"] = ensemble("Summer12_TTJets1725_flavor:down_FlavorQCD/job_*_ensemble.root", 62131965./1.8);

  sample.ensembles["flavorBUp"  ] = ensemble("topmass_140305_1201/Z2_S12_BJES_Up_sig/MASS_172.5/JES_1.00/SHAPE_1.00/PERMU_0.0/FSIG_0.7685747868/WEIGHT_weight.combinedWeight/job_*_ensemble.root", 62131965./1.8);
  sample.ensembles["flavorBDown"] = ensemble("topmass_140305_1201/Z2_S12_BJES_Down_sig/MASS_172.5/JES_1.00/SHAPE_1.00/PERMU_0.0/FSIG_0.7673455689/WEIGHT_weight.combinedWeight/job_*_ensemble.root", 62131965./1.8);

  //sample.ensembles["flavorQUp"] = ensemble("Summer12_TTJets1725_flavor:up_FlavorPureQuark/job_*_ensemble.root", 62131965./1.8);
  //sample.ensembles["flavorQDown"] = ensemble("Summer12_TTJets1725_flavor:down_FlavorPureQuark/job_*_ensemble.root", 62131965./1.8);

  //sample.ensembles["flavorGUp"] = ensemble("Summer12_TTJets1725_flavor:up_FlavorPureGluon/job_*_ensemble.root", 62131965./1.8);
  //sample.ensembles["flavorGDown"] = ensemble("Summer12_TTJets1725_flavor:down_FlavorPureGluon/job_*_ensemble.root", 62131965./1.8);

  sample.ensembles["bTagSFUp"    ] = ensemble("topmass_140305_1201a/Z2_S12_BTAG_Up/job_*_ensemble.root", 62131965./1.8);
  sample.ensembles["bTagSFDown"  ] = ensemble("topmass_140305_1201a/Z2_S12_BTAG_Down/job_*_ensemble.root", 62131965./1.8);
  sample.ensembles["misTagSFUp"  ] = ensemble("topmass_140305_1201a/Z2_S12_MTAG_Up/job_*_ensemble.root", 62131965./1.8);
  sample.ensembles["misTagSFDown"] = ensemble("topmass_140305_1201a/Z2_S12_MTAG_Down/job_*_ensemble.root", 62131965./1.8);

  sample.ensembles["topPt"] = ensemble("topmass_140305_1201e/Z2_S12_TopPt/job_*_ensemble.root", 62131965./1.8);

  sample.ensembles["fSigUp"  ] = ensemble("topmass_140305_1201a/Z2_S12_FSIG_Up/job_*_ensemble.root", 62131965./1.8);
  sample.ensembles["fSigDown"] = ensemble("topmass_140305_1201a/Z2_S12_FSIG_Down/job_*_ensemble.root", 62131965./1.8);

  sample.ensembles["matchingUp"  ] = ensemble("topmass_140305_1201/Z2_S12_Matching_Up_MadSpin_sig/MASS_172.5/JES_1.00/SHAPE_1.00/PERMU_0.0/FSIG_0.6322077386/WEIGHT_weight.combinedWeight/job_*_ensemble.root", 37083003./1.8);
  sample.ensembles["matchingDown"] = ensemble("topmass_140305_1201/Z2_S12_Matching_Down_MadSpin_sig/MASS_172.5/JES_1.00/SHAPE_1.00/PERMU_0.0/FSIG_0.6615011853/WEIGHT_weight.combinedWeight/job_*_ensemble.root", 34053113./1.8);

  sample.ensembles["scaleUp"  ] = ensemble("topmass_140305_1201/Z2_S12_Scale_Up_MadSpin_sig/MASS_172.5/JES_1.00/SHAPE_1.00/PERMU_0.0/FSIG_0.7123135044/WEIGHT_weight.combinedWeight/job_*_ensemble.root", 41908271./1.8);
  sample.ensembles["scaleDown"] = ensemble("topmass_140305_1201/Z2_S12_Scale_Down_MadSpin_sig/MASS_172.5/JES_1.00/SHAPE_1.00/PERMU_0.0/FSIG_0.635286828/WEIGHT_weight.combinedWeight/job_*_ensemble.root", 39286663./1.8);

  sample.ensembles["jerUp"  ] = ensemble("topmass_140305_1201/Z2_S12_JER_Up_sig/MASS_172.5/JES_1.00/SHAPE_1.00/PERMU_0.0/FSIG_0.7532557569/WEIGHT_weight.combinedWeight/job_*_ensemble.root", 62131965./1.8);
  sample.ensembles["jerDown"] = ensemble("topmass_140305_1201/Z2_S12_JER_Down_sig/MASS_172.5/JES_1.00/SHAPE_1.00/PERMU_0.0/FSIG_0.7835384716/WEIGHT_weight.combinedWeight/job_*_ensemble.root", 62131965./1.8);

  sample.ensembles["jes1Up"  ] = ensemble("topmass_140305_1201/Z2_S12_COR1JES_Up_sig/MASS_172.5/JES_1.00/SHAPE_1.00/PERMU_0.0/FSIG_0.7672252189/WEIGHT_weight.combinedWeight/job_*_ensemble.root", 62131965./1.8);
  sample.ensembles["jes1Down"] = ensemble("topmass_140305_1201/Z2_S12_COR1JES_Down_sig/MASS_172.5/JES_1.00/SHAPE_1.00/PERMU_0.0/FSIG_0.7695591532/WEIGHT_weight.combinedWeight/job_*_ensemble.root", 62131965./1.8);
  sample.ensembles["jes2Up"  ] = ensemble("topmass_140305_1201/Z2_S12_COR2JES_Up_sig/MASS_172.5/JES_1.00/SHAPE_1.00/PERMU_0.0/FSIG_0.7748250179/WEIGHT_weight.combinedWeight/job_*_ensemble.root", 62131965./1.8);
  sample.ensembles["jes2Down"] = ensemble("topmass_140305_1201/Z2_S12_COR2JES_Down_sig/MASS_172.5/JES_1.00/SHAPE_1.00/PERMU_0.0/FSIG_0.7601530052/WEIGHT_weight.combinedWeight/job_*_ensemble.root", 62131965./1.8);
  sample.ensembles["jes3Up"  ] = ensemble("topmass_140305_1201/Z2_S12_COR3JES_Up_sig/MASS_172.5/JES_1.00/SHAPE_1.00/PERMU_0.0/FSIG_0.76795735/WEIGHT_weight.combinedWeight/job_*_ensemble.root", 62131965./1.8);
  sample.ensembles["jes3Down"] = ensemble("topmass_140305_1201/Z2_S12_COR3JES_Down_sig/MASS_172.5/JES_1.00/SHAPE_1.00/PERMU_0.0/FSIG_0.7681779725/WEIGHT_weight.combinedWeight/job_*_ensemble.root", 62131965./1.8);
  sample.ensembles["jes4Up"  ] = ensemble("topmass_140305_1201/Z2_S12_COR4JES_Up_sig/MASS_172.5/JES_1.00/SHAPE_1.00/PERMU_0.0/FSIG_0.7835827954/WEIGHT_weight.combinedWeight/job_*_ensemble.root", 62131965./1.8);
  sample.ensembles["jes4Down"] = ensemble("topmass_140305_1201/Z2_S12_COR4JES_Down_sig/MASS_172.5/JES_1.00/SHAPE_1.00/PERMU_0.0/FSIG_0.7514697972/WEIGHT_weight.combinedWeight/job_*_ensemble.root", 62131965./1.8);

  sample.ensembles["triggerUp"  ] = ensemble("topmass_140305_1201/Z2_S12_TRIGGERJES_Up_sig/MASS_172.5/JES_1.00/SHAPE_1.00/PERMU_0.0/FSIG_1.0000000000/WEIGHT_weight.combinedWeight/job_*_ensemble.root", 62131965./1.8);
  sample.ensembles["triggerDown"] = ensemble("topmass_140305_1201/Z2_S12_TRIGGERJES_Down_sig/MASS_172.5/JES_1.00/SHAPE_1.00/PERMU_0.0/FSIG_0.4383580517/WEIGHT_weight.combinedWeight/job_*_ensemble.root", 62131965./1.8);

  sample.ensembles["jesPuBBUp"  ] = ensemble("topmass_140305_1201d/Z2_S12_PUPTBBJES_Up/job_*_ensemble.root", 62131965./1.8);
  sample.ensembles["jesPuBBDown"] = ensemble("topmass_140305_1201d/Z2_S12_PUPTBBJES_Down/job_*_ensemble.root", 62131965./1.8);
  sample.ensembles["jesPuECUp"  ] = ensemble("topmass_140305_1201d/Z2_S12_PUPTECJES_Up/job_*_ensemble.root", 62131965./1.8);
  sample.ensembles["jesPuECDown"] = ensemble("topmass_140305_1201d/Z2_S12_PUPTECJES_Down/job_*_ensemble.root", 62131965./1.8);

  sample.ensembles["mcatnlo"     ] = ensemble("topmass_140305_1201/Z2_S12_MCNLO_sig/MASS_172.5/JES_1.00/SHAPE_1.00/PERMU_0.0/FSIG_0.5557773815/WEIGHT_weight.combinedWeight/job_*_ensemble.root", 32852589./1.8);
  sample.ensembles["powheg"      ] = ensemble("topmass_140305_1201/Z2_S12_POWHEG_sig/MASS_172.5/JES_1.00/SHAPE_1.00/PERMU_0.0/FSIG_0.7029775376/WEIGHT_weight.combinedWeight/job_*_ensemble.root", 21675970./1.8);
  sample.ensembles["powhegHerwig"] = ensemble("topmass_140305_1201/Z2_S12_POWHER_sig/MASS_172.5/JES_1.00/SHAPE_1.00/PERMU_0.0/FSIG_0.8471413429/WEIGHT_weight.combinedWeight/job_*_ensemble.root", 27684235./1.8);

  sample.ensembles["defaultSC"] = ensemble("topmass_140305_1201/Z2_S12_Z2_sig/MASS_172.5/JES_1.00/SHAPE_1.00/PERMU_0.0/FSIG_0.7916110679/WEIGHT_weight.combinedWeight/job_*_ensemble.root", 41761265./6./6.*9.*9./1.8);
  sample.ensembles["P11"      ] = ensemble("topmass_140305_1201/Z2_S12_P11_sig/MASS_172.5/JES_1.00/SHAPE_1.00/PERMU_0.0/FSIG_0.7238466061/WEIGHT_weight.combinedWeight/job_*_ensemble.root", 11651739./6./6.*9.*9./1.8);
  sample.ensembles["P11noCR"  ] = ensemble("topmass_140305_1201/Z2_S12_P11NoCR_sig/MASS_172.5/JES_1.00/SHAPE_1.00/PERMU_0.0/FSIG_0.7413235903/WEIGHT_weight.combinedWeight/job_*_ensemble.root", 11919063./6./6.*9.*9./1.8);
  sample.ensembles["P11mpiHi" ] = ensemble("topmass_140305_1201/Z2_S12_P11mpiHi_sig/MASS_172.5/JES_1.00/SHAPE_1.00/PERMU_0.0/FSIG_0.7353096816/WEIGHT_weight.combinedWeight/job_*_ensemble.root", 7953758./6./6.*9.*9./1.8);
  sample.ensembles["P11TeV"   ] = ensemble("topmass_140305_1201/Z2_S12_P11TeV_sig/MASS_172.5/JES_1.00/SHAPE_1.00/PERMU_0.0/FSIG_0.7431802926/WEIGHT_weight.combinedWeight/job_*_ensemble.root", 7946264./6./6.*9.*9./1.8);

  ///////////////////////////////////

  sample.comparisons["Pile-up (pp cross-section)       "] = comparison("default", "puUp", "puDown", true);
  //sample.comparisons["Jet energy response (udsc)       "] = comparison("default", "flavorQUp", "flavorQDown", true);
  //sample.comparisons["Jet energy response (gluon)      "] = comparison("default", "flavorGUp", "flavorGDown", true);
  sample.comparisons["Jet energy response (b)          "] = comparison("default", "flavorBUp", "flavorBDown", true);
  //sample.comparisons["Jet energy response (FlavorQCD)  "] = comparison("default", "flavorQCDUp", "flavorQCDDown", true, false);
  sample.comparisons["b fragmentation                  "] = comparison("default", "bFragLEP", "", true);
  sample.comparisons["Semi-leptonic B hadron decays    "] = comparison("default", "bFNuUp", "bFNuDown", true);
  sample.comparisons["b-tag rate                       "] = comparison("default", "bTagSFUp", "bTagSFDown", true);
  sample.comparisons["b-tag (mistag rate)              "] = comparison("default", "misTagSFUp", "misTagSFDown", true);
  sample.comparisons["Trigger                          "] = comparison("default", "triggerUp", "triggerDown", true);
  sample.comparisons["Top-pt reweighting               "] = comparison("default", "topPt", "", true);
  sample.comparisons["ME-PS matching threshold         "] = comparison("default", "matchingUp", "matchingDown", false);
  sample.comparisons["$Q^{2}$ scale                    "] = comparison("default", "scaleUp", "scaleDown", false);
  sample.comparisons["Jet energy resolution            "] = comparison("default", "jerUp" , "jerDown" , true);
  sample.comparisons["\\pt- and $\\eta$-dependent JES1 "] = comparison("default", "jes1Up", "jes1Down", true);
  sample.comparisons["\\pt- and $\\eta$-dependent JES2 "] = comparison("default", "jes2Up", "jes2Down", true);
  sample.comparisons["\\pt- and $\\eta$-dependent JES3 "] = comparison("default", "jes3Up", "jes3Down", true);
  sample.comparisons["\\pt- and $\\eta$-dependent JES4 "] = comparison("default", "jes4Up", "jes4Down", true);
  sample.comparisons["Pile-up (JES) BB                 "] = comparison("default", "jesPuBBUp", "jesPuBBDown", true);
  sample.comparisons["Pile-up (JES) EC                 "] = comparison("default", "jesPuECUp", "jesPuECDown", true);
  sample.comparisons["MadGraph (no SC) vs. Powheg      "] = comparison("default", "powheg", "", false, false);
  sample.comparisons["Powheg+Pythia6 vs. MC@NLO+Herwig6"] = comparison("powheg", "mcatnlo", "", false, false);
  sample.comparisons["Powheg+Pythia6 vs. Powheg+Herwig6"] = comparison("powheg", "powhegHerwig", "", false, true);
  sample.comparisons["MC@NLO+Herwig6 vs. Powheg+Herwig6"] = comparison("mcatnlo", "powhegHerwig", "", false, false);
  sample.comparisons["ME generator                     "] = comparison("defaultSC", "powheg", "", false);
  sample.comparisons["Spin correlations                "] = comparison("default", "defaultSC", "", false, true);
  sample.comparisons["Pythia Z2* vs. P11               "] = comparison("defaultSC", "P11", "", false, false);
  sample.comparisons["Color reconnection               "] = comparison("P11", "P11noCR", "", false);
  sample.comparisons["Underlying event                 "] = comparison("P11", "P11mpiHi", "P11TeV", false);
  sample.comparisons["Non-\\ttbar background fSig      "] = comparison("default", "fSigUp", "fSigDown", true);
  //sample.comparisons["Non-\\ttbar background shape     "] = comparison("default", "", "", true);
}

void SystematicUncertainties::deriveSystematics()
{
  //fillLeptonJets();
  fillAllJets();

  std::map<std::string, double> totalUncertainties2;
  for(auto& var : sample.variables)
    totalUncertainties2[var] = 0;

  std::cout << "\n### Fitting pseudo-experiments" << std::endl;
  for(std::map<std::string, ensemble>::iterator it = sample.ensembles.begin(); it != sample.ensembles.end(); it++) {
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
  //myfile    << "\n### Systematic uncertainties" << "\n";
  std::cout << "Uncertainty name \t \t \t  2D mass \t  JSF \t \t    1D mass" << std::endl;
  myfile    << "Uncertainty name & 2D mass & JSF & 1D mass \\tabularnewline\n\\hline" << "\n";
  for(std::map<std::string, comparison>::iterator it = sample.comparisons.begin(); it != sample.comparisons.end(); it++) {
    myfile    << "\\hline" << "\n";
    if (!it->second.active) std::cout << "(cc) ";
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
    else down         = sample.ensembles.find(it->second.down)->second;

    std::map<std::string, double> shifts;
    for(auto& var : sample.variables)
      shifts[var] = std::max(std::abs(nominal.values[var].first-up.values[var].first), std::abs(nominal.values[var].first-down.values[var].first));

    std::map<std::string, double> shiftUncs;
    for(auto& var : sample.variables)
      shiftUncs[var] = 0.;

    if (!it->second.correlated) {
      for(auto& var : sample.variables)
        shiftUncs[var] = std::max(sqrt(pow(nominal.values[var].second,2)+pow(up.values[var].second,2)),sqrt(pow(nominal.values[var].second,2)+pow(down.values[var].second,2)));
    }

    sprintf(buffer," & %.2lf$\\pm$%.2lf & %.3lf$\\pm$%.3lf & %.2lf$\\pm$%.2lf \\tabularnewline\n", shifts[sample.variables[0]], shiftUncs[sample.variables[0]], shifts[sample.variables[1]], shiftUncs[sample.variables[1]], shifts[sample.variables[2]], shiftUncs[sample.variables[2]]);
    printf(" \t  %.2lf +/- %.2lf   %.3lf +/- %.3lf   %.2lf +/- %.2lf \n", shifts[sample.variables[0]], shiftUncs[sample.variables[0]], shifts[sample.variables[1]], shiftUncs[sample.variables[1]], shifts[sample.variables[2]], shiftUncs[sample.variables[2]]);
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
        totalUncertainties2[var] += pow(std::max(shifts[var], shiftUncs[var]), 2);
      }
    }
  }

  std::cout << "\n### Total uncertainties" << std::endl;
  myfile    << "\\hline\n\\hline\nTotal uncertainties";
  sprintf(buffer, "& %.2lf & %.3lf & %.2lf \n", sqrt(totalUncertainties2[sample.variables[0]]), sqrt(totalUncertainties2[sample.variables[1]]), sqrt(totalUncertainties2[sample.variables[2]]));
  printf("\t 2D mass = %.2lf GeV, JSF = %.3lf, 1D mass = %.2lf GeV\n", sqrt(totalUncertainties2[sample.variables[0]]), sqrt(totalUncertainties2[sample.variables[1]]), sqrt(totalUncertainties2[sample.variables[2]]));
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


