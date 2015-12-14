#include "TemplateDerivation.h"

#include "TAxis.h"
#include "TCanvas.h"
#include "TF1.h"
#include "TFile.h"
#include "TLatex.h"
#include "TLegend.h"
#include "TROOT.h"
#include "TString.h"
#include "TStyle.h"
#include "TLorentzVector.h"
#include "TClonesArray.h"
#include "TEventList.h"
#include "TTreeFormula.h"

#include "RooAddition.h"
#include "RooAddPdf.h"
#include "RooBernstein.h"
#include "RooBifurGauss.h"
#include "RooCategory.h"
#include "RooChebychev.h"
#include "RooDataHist.h"
#include "RooDataSet.h"
#include "RooFitResult.h"
#include "RooGamma.h"
#include "RooCBShape.h"
#include "RooGaussian.h"
#include "RooLandau.h"
#include "RooMinimizer.h"
#include "RooMinuit.h"
#include "RooNLLVar.h"
#include "RooPlot.h"
#include "RooProdPdf.h"
#include "RooRealVar.h"
#include "RooSimultaneous.h"
#include "RooSimWSTool.h"
#include "RooVoigtian.h"
#include "RooWorkspace.h"
#include "RooFactoryWSTool.h"

#include "TopMass/TopEventTree/interface/JetEvent.h"
#include "TopMass/TopEventTree/interface/TopEvent.h"
#include "TopMass/TopEventTree/interface/WeightEvent.h"

#include "ProgramOptionsReader.h"

#include "iostream"
#include "fstream"

#include <boost/algorithm/string/split.hpp>
#include <boost/algorithm/string/classification.hpp>
#include <boost/format.hpp>

typedef ProgramOptionsReader po;

TemplateDerivation::TemplateDerivation()
    : selection_(po::GetOption<std::string>("analysisConfig.selection")),
      samplePath_(po::GetOption<std::string>("analysisConfig.samplePath")),
      fVar1_(po::GetOption<std::string>("analysisConfig.var1")),
      fVar2_(po::GetOption<std::string>("analysisConfig.var2")),
      fVar3_(po::GetOption<std::string>("analysisConfig.var3")),
      fWeight_(po::GetOption<std::string>("weight")),
      fChannel_(po::GetOption<std::string>("channel")),
      activeBranches_(
          po::GetOption<std::string>("analysisConfig.activeBranches")),
      maxPermutations_(po::GetOption<int>("analysisConfig.maxPermutations")),
      maxMtop_(po::GetOption<double>("templates.maxTopMass")),
      workspace_(0) {
  if (!strncmp(fChannel_, "alljets", 7))
    channelID_ = kAllJets;
  else if (!strncmp(fChannel_, "muon", 4))
    channelID_ = kMuonJets;
  else if (!strncmp(fChannel_, "electron", 8))
    channelID_ = kElectronJets;
  else if (!strncmp(fChannel_, "lepton", 6))
    channelID_ = kLeptonJets;
  else
    UnknownChannelAbort();

  jesValues_ = {0.96, 0.98, 1.00, 1.02, 1.04};
  massValues_ = {166.5, 169.5, 171.5, 172.5, 173.5, 175.5, 178.5};

  rooFitTopMass();
}

TemplateDerivation::~TemplateDerivation() {}

std::string TemplateDerivation::constructFileName(double mass, double jes) {
  std::string fileName = str(boost::format("Summer12_TTJetsMS%1%%2%_%3$3.2f") %
                             (int)mass % (int)((mass - int(mass)) * 10) % jes);
  if (channelID_ == kAllJets) {
    fileName += "_alljets.root";
  }
  return fileName;
}

std::string TemplateDerivation::constructTemplateName(double mass, double jes) {
  std::string name = str(boost::format("jes%3%%4$02imass%1%%2%") % (int)mass %
                         (int)((mass - int(mass)) * 10) % (int)jes %
                         round((jes - int(jes)) * 100));
  return name;
}

std::vector<RooDataSet *> TemplateDerivation::createDataSets(
    const RooArgSet &varSet) {
  TFile *tmpFile = TFile::Open("tmpFileRooFitTopMass.root", "RECREATE");
  std::vector<RooDataSet *> dataset;
  for (auto mass : massValues_) {
    for (auto jes : jesValues_) {
      std::string fileName = constructFileName(mass, jes);
      std::cout << "Creating RooDataSet for: " << fileName;
      TChain *chain = new TChain("analyzeKinFit/eventTree");
      std::cout << " (nFiles: " << chain->Add(samplePath_ + TString(fileName))
                << ") ";  // << std::endl;
      tmpFile->cd();
      TTree *tree = modifiedTree(chain);  //, minComboType, maxComboType);
      // needed to avoid a crash
      double co, /*pr,*/ mt, mw, cw;
      tree->SetBranchAddress("comboType", &co);
      // tree->SetBranchAddress("prob", &pr);
      tree->SetBranchAddress("topMass", &mt);
      tree->SetBranchAddress("meanWMass", &mw);
      tree->SetBranchAddress("combinedWeight", &cw);
      TString name = constructName("dataset", dataset.size());
      dataset.push_back(new RooDataSet(name, name, varSet,
                                       RooFit::Import(*tree),
                                       RooFit::WeightVar("combinedWeight")));
      // if(binnedTemplates){
      //  for(unsigned h=0; h<1; ++h) {
      //    name = "hist_"; name += iTempl;
      //    hist[h][iTempl] = new RooDataHist(name, name, mTop,
      //    *dataset[iTempl]);
      //  }
      //}
      workspace_->import(*dataset.back());
    }
  }
  tmpFile->Close();
  return dataset;
}

void TemplateDerivation::addTemplateFunction(int templType, int comboType,
                                             int nComboTypes, RooRealVar *mTop,
                                             RooRealVar *JES) {
  int index = nComboTypes * templType + comboType;

  RooRealVar MassOffset = RooRealVar("MassOffset", "MassOffset", 0.);
  MassOffset.setConstant(kTRUE);
  RooRealVar MassSlopeMass = RooRealVar("MassSlopeMass", "MassSlopeMass", 0.);
  MassSlopeMass.setConstant(kTRUE);
  RooRealVar MassSlopeJES = RooRealVar("MassSlopeJES", "MassSlopeJES", 0.);
  MassSlopeJES.setConstant(kTRUE);
  RooRealVar MassSlopeMassJES =
      RooRealVar("MassSlopeMassJES", "MassSlopeMassJES", 0.);
  MassSlopeMassJES.setConstant(kTRUE);
  RooRealVar JESOffset = RooRealVar("JESOffset", "JESOffset", 0.);
  JESOffset.setConstant(kTRUE);
  RooRealVar JESSlopeMass = RooRealVar("JESSlopeMass", "JESSlopeMass", 0.);
  JESSlopeMass.setConstant(kTRUE);
  RooRealVar JESSlopeJES = RooRealVar("JESSlopeJES", "JESSlopeJES", 0.);
  JESSlopeJES.setConstant(kTRUE);
  RooRealVar JESSlopeMassJES =
      RooRealVar("JESSlopeMassJES", "JESSlopeMassJES", 0.);
  JESSlopeMassJES.setConstant(kTRUE);
  RooRealVar MassOffset2 = RooRealVar("MassOffset2", "MassOffset2", 0.);
  MassOffset2.setConstant(kTRUE);
  RooRealVar MassSlopeMass2 =
      RooRealVar("MassSlopeMass2", "MassSlopeMass2", 0.);
  MassSlopeMass2.setConstant(kTRUE);
  RooRealVar MassSlopeJES2 = RooRealVar("MassSlopeJES2", "MassSlopeJES2", 0.);
  MassSlopeJES2.setConstant(kTRUE);
  RooRealVar MassSlopeMassJES2 =
      RooRealVar("MassSlopeMassJES2", "MassSlopeMassJES2", 0.);
  MassSlopeMassJES2.setConstant(kTRUE);
  RooRealVar JESOffset2 = RooRealVar("JESOffset2", "JESOffset2", 0.);
  JESOffset2.setConstant(kTRUE);
  RooRealVar JESSlopeMass2 = RooRealVar("JESSlopeMass2", "JESSlopeMass2", 0.);
  JESSlopeMass2.setConstant(kTRUE);
  RooRealVar JESSlopeJES2 = RooRealVar("JESSlopeJES2", "JESSlopeJES2", 0.);
  JESSlopeJES2.setConstant(kTRUE);
  RooRealVar JESSlopeMassJES2 =
      RooRealVar("JESSlopeMassJES2", "JESSlopeMassJES2", 0.);
  JESSlopeMassJES2.setConstant(kTRUE);
  // 0 = JES, 1 = mTop, 2 = Offset, 3 = SlopeMass, 4 = SlopeJES, 5 =
  // SlopeMassJES
  TString formula_mTop =
      "@1+@2+@3*(@1-172.5)+@4*(@0-1.0)+@5*(@1-172.5)*(@0-1.0)";
  TString formula_JES =
      "@0+@2+@3*(@1-172.5)+@4*(@0-1.0)+@5*(@1-172.5)*(@0-1.0)";
  TString name = "mTop_intermediate_";
  name += index;
  RooFormulaVar mTop_intermediate =
      RooFormulaVar(name, name, formula_mTop,
                    RooArgSet(*JES, *mTop, MassOffset, MassSlopeMass,
                              MassSlopeJES, MassSlopeMassJES));
  name = "mJES_intermediate_";
  name += index;
  RooFormulaVar JES_intermediate = RooFormulaVar(
      name, name, formula_JES, RooArgSet(*JES, *mTop, JESOffset, JESSlopeMass,
                                         JESSlopeJES, JESSlopeMassJES));
  name = "mTop_corrected_";
  name += index;
  RooFormulaVar mTop_corrected = RooFormulaVar(
      name, name, formula_mTop,
      RooArgSet(JES_intermediate, mTop_intermediate, MassOffset2,
                MassSlopeMass2, MassSlopeJES2, MassSlopeMassJES2));
  name = "mJES_corrected_";
  name += index;
  RooFormulaVar JES_corrected =
      RooFormulaVar(name, name, formula_JES,
                    RooArgSet(JES_intermediate, mTop_intermediate, JESOffset2,
                              JESSlopeMass2, JESSlopeJES2, JESSlopeMassJES2));
  // double a[] = {84.6333,  0.0205799 , 16.0275, -0.280052,
  //              4.76909,  0.018934  , 2.18012, -0.12839 ,
  //              7.08441, -0.00812571, -6.2732,  0.22637 };
  // iniPar = std::vector<double>(a,a+sizeof(a)/sizeof(double));

  RooRealVar *topMass =
      new RooRealVar("topMass", "m_{t}^{fit}", 100., maxMtop_, "GeV");
  RooRealVar *meanWMass =
      new RooRealVar("meanWMass", "m_{W}^{rec}", 50., 300., "GeV");
  RooRealVar *var = 0;
  TString comboRangeName = "";
  std::vector<RooRealVar *> par;
  std::vector<RooFormulaVar *> alpha;
  std::vector<double> iniPar;

  if (channelID_ == kAllJets) {
    if (comboType == 0) {
      comboRangeName = "CP";
    }
    // else if(comboType == 1){ comboRangeName = "R4"; }
    // else if(comboType == 2){ comboRangeName = "R5"; }
    else if (comboType == 1) {
      comboRangeName = "WP";
    }

    if (templType == 0) {
      if (comboType == 0) {  // mTop, correct
        iniPar = {172.399, 0.989847,  81.225,  0.758306,
                  7.92651, 0.0730216, 5.03205, 0.0241954};
      } else if (comboType == 1) {  // mTop, wrong
        iniPar = {190.009, 0.23107,   123.714,  -2.41287,   173.325, 1.11063,
                  80.0112, 1.32969,   26.2827,  -0.0172727, 35.0065, -0.442854,
                  10.0152, 0.0656903, -4.72985, -0.108491};
      } else if (comboType == 2) {  // mTop, unmatched
        iniPar = {178.883, 0.306507,  40.7573, 0.445759,  174.681, 1.06829,
                  83.8601, 0.280563,  22.065,  0.0542991, 4.93823, 0.3491,
                  9.74,    0.0548014, 2.31277, -0.181041};
      }
    } else if (templType == 1) {
      if (comboType == 0) {  // mW, correct
        iniPar = {84.4535, -0.0154884, 91.1495,  -0.00943573,
                  5.20726, -0.0115877, 23.4646,  -0.0249863,
                  6.75636, 0.0019251,  -19.8341, 0.0724597};
      } else if (comboType == 1) {  // mW, wrong
        iniPar = {87.2093,    -0.00584533, 20.9648,  -0.0912001, 3.91928,
                  0.00521933, 3.8585,      -0.28899, 5.9744,     0.00351353,
                  -5.70177,   -0.231188,   3.80803,  0.00289055, -11.4316,
                  -0.264811,  4.89211,     4.89211,  0.25,       0.5};
      } else if (comboType == 2) {  // mW, unmatched
        iniPar = {86.6185, -0.0184471, 32.1242,  0.0815632,
                  6.48169, -0.0108763, 7.37265,  0.0139849,
                  6.48169, -0.0108763, 7.37265,  0.0139849,
                  7.91363, 0.00815331, -17.3423, -0.00100213};
      }
    }

    for (unsigned p = 0; p < iniPar.size(); ++p) {
      // std::cout << iniPar[p] << ", ";
      name = constructName("par", index, p);
      RooRealVar *myVar = new RooRealVar(name, name, iniPar[p]);
      myVar->setConstant(kFALSE);
      if (p == 16) myVar->setVal(5.0);
      if (p == 16) myVar->setConstant();
      // if(p == 17) myVar->setVal(7.5);
      if (p == 17) myVar->setVal(5.0);
      if (p == 17) myVar->setConstant();

      if (p == 18) myVar->setVal(0.25);
      if (p == 18) myVar->setConstant();
      if (p == 19) myVar->setVal(0.5);
      if (p == 19) myVar->setConstant();
      //}

      par.push_back(myVar);
    }
    std::cout << std::endl;
    name = "sig_";
    name += index;
    TString varName;
    if (templType == 0) {
      var = topMass;
      if (comboType == 0) {
        fillAlpha(alpha, index, RooArgSet(*par[0], *par[1], *par[2], *par[3],
                                          mTop_corrected, JES_corrected));
        fillAlpha(alpha, index, RooArgSet(*par[4], *par[5], *par[6], *par[7],
                                          mTop_corrected, JES_corrected));
        varName = constructName("width", index);
        RooRealVar *width = new RooRealVar(varName, varName, 2.0);
        width->setConstant(kTRUE);
        RooVoigtian voigt =
            RooVoigtian(name, name, *var, *alpha[0], *width, *alpha[1]);
        workspace_->import(voigt);
      } else if (comboType == 1 || comboType == 2) {
        varName = constructName("ratio", index);
        RooRealVar *ratio = new RooRealVar(varName, varName, 0.0, 1.0);
        if (comboType == 1) ratio->setVal(0.700556);  // 0.648584); //0.662784);
                                                      // //0.820546);
        else if (comboType == 2)
          ratio->setVal(0.685704);  // 0.774588); //0.758690);
        // ratio->setConstant(kTRUE);
        ratio->setConstant(kFALSE);
        fillAlpha(alpha, index, RooArgSet(*par[0], *par[1], *par[2], *par[3],
                                          mTop_corrected, JES_corrected));
        fillAlpha(alpha, index, RooArgSet(*par[8], *par[9], *par[10], *par[11],
                                          mTop_corrected, JES_corrected));
        varName = constructName("landau", index);
        RooLandau *landau =
            new RooLandau(varName, varName, *var, *alpha[0], *alpha[1]);
        fillAlpha(alpha, index, RooArgSet(*par[4], *par[5], *par[6], *par[7],
                                          mTop_corrected, JES_corrected));
        fillAlpha(alpha, index,
                  RooArgSet(*par[12], *par[13], *par[14], *par[15],
                            mTop_corrected, JES_corrected));
        varName = constructName("width", index);
        RooRealVar *width = new RooRealVar(varName, varName, 2.0);
        width->setConstant(kTRUE);
        varName = "voigt_";
        varName += index;
        // RooVoigtian *voigt = new RooVoigtian(varName, varName, *var,
        // *alpha[2], *width, *alpha[3]);
        RooGaussian *voigt =
            new RooGaussian(varName, varName, *var, *alpha[2], *alpha[3]);
        RooAddPdf *add = new RooAddPdf(name, name, *landau, *voigt, *ratio);
        workspace_->import(*add);
      }
    } else if (templType == 1) {
      var = meanWMass;
      if (comboType == 0) {  // mW, correct
        fillAlpha(alpha, index, RooArgSet(*par[0], *par[1], *par[2], *par[3],
                                          mTop_corrected, JES_corrected));
        fillAlpha(alpha, index, RooArgSet(*par[4], *par[5], *par[6], *par[7],
                                          mTop_corrected, JES_corrected));
        fillAlpha(alpha, index, RooArgSet(*par[8], *par[9], *par[10], *par[11],
                                          mTop_corrected, JES_corrected));
        RooBifurGauss *asymGaus = new RooBifurGauss(name, name, *var, *alpha[0],
                                                    *alpha[1], *alpha[2]);
        workspace_->import(*asymGaus);
      } else {
        // fillAlpha(alpha, h, RooArgSet(*par[ 0], *par[ 1], *par[ 2],
        // *par[ 3], mTop_corrected, JES_corrected), "-7.5");
        fillAlpha(alpha, index,
                  RooArgSet(*par[0], *par[1], *par[2], *par[3], mTop_corrected,
                            JES_corrected, *par[16]),
                  "-@6");
        fillAlpha(alpha, index, RooArgSet(*par[0], *par[1], *par[2], *par[3],
                                          mTop_corrected, JES_corrected));
        // fillAlpha(alpha, h, RooArgSet(*par[ 0], *par[ 1], *par[ 2],
        // *par[ 3], mTop_corrected, JES_corrected), "+7.5");
        fillAlpha(alpha, index,
                  RooArgSet(*par[0], *par[1], *par[2], *par[3], mTop_corrected,
                            JES_corrected, *par[17]),
                  "+@6");
        fillAlpha(alpha, index, RooArgSet(*par[4], *par[5], *par[6], *par[7],
                                          mTop_corrected, JES_corrected));
        fillAlpha(alpha, index, RooArgSet(*par[8], *par[9], *par[10], *par[11],
                                          mTop_corrected, JES_corrected));
        fillAlpha(alpha, index,
                  RooArgSet(*par[12], *par[13], *par[14], *par[15],
                            mTop_corrected, JES_corrected));
        varName = constructName("gaus1", index);
        // RooRealVar *ratio1 = new RooRealVar("ratioGaus1","",0.25);
        // RooRealVar *ratio2 = new RooRealVar("ratioGaus2","",0.5 );
        // RooRealVar *ratio3 = new RooRealVar("ratioGaus3","",0.25);
        RooGaussian *gaus1 =
            new RooGaussian(varName, varName, *var, *alpha[0], *alpha[3]);
        varName = constructName("gaus2", index);
        RooGaussian *gaus2 =
            new RooGaussian(varName, varName, *var, *alpha[1], *alpha[4]);
        varName = constructName("gaus3", index);
        RooGaussian *gaus3 =
            new RooGaussian(varName, varName, *var, *alpha[2], *alpha[5]);
        RooAddPdf *add = new RooAddPdf(
            name, name, RooArgList(*gaus1, *gaus2, *gaus3),
            RooArgList(*par[18],
                       *par[19]) /*RooArgList(*ratio1, *ratio2, *ratio3)*/);
        workspace_->import(*add);
      }
    }
  }
  unsigned int nPars = numVariables(constructName("par", index));
  assert(par.size() == nPars);
  std::cout << "old: ";
  for (unsigned int i = 0; i < nPars; ++i) {
    std::cout << iniPar[i];
    if (i != nPars - 1) std::cout << "|";
  }
  std::cout << std::endl;
}

RooFitResult *TemplateDerivation::fitTemplate(int templType, int comboType,
                                              int nComboTypes) {
  RooSimWSTool *simWST = new RooSimWSTool(*workspace_);
  const int index = nComboTypes * templType + comboType;
  TString name = constructName("sim", index);
  TString name_pdf = constructName("sig", index);
  simWST->build(name, name_pdf, RooFit::SplitParam("JES", "cat_templ"),
                RooFit::SplitParam("mTop", "cat_templ"));

  RooArgSet nllSet;
  int templateIndex = 0;
  for (auto mass : massValues_) {
    for (auto jes : jesValues_) {
      RooAbsData *dataset =
          workspace_->data(constructName("dataset", templateIndex));
      std::string templateName = constructTemplateName(mass, jes);
      workspace_->var((std::string("JES_") + templateName).c_str())
          ->setVal(jes);
      workspace_->var((std::string("mTop_") + templateName).c_str())
          ->setVal(mass);
      TString name = constructName("nll", index, templateIndex);
      TString name_pdf = "sig_";
      name_pdf += index;
      name_pdf += "_" + templateName;
      // if(binnedTemplates)
      //  nll[h][t] = new RooNLLVar(name, name,
      //  *workspace_->pdf(name_pdf), *hist[h][t], RooFit::NumCPU(1),
      //  RooFit::Range(100., 550.));
      // else
      TString comboRangeName;
      if (channelID_ == kAllJets) {
        if (comboType == 0) {
          comboRangeName = "CP";
        }
        // else if(comboType == 1){ comboRangeName = "R4"; }
        // else if(comboType == 2){ comboRangeName = "R5"; }
        else if (comboType == 1) {
          comboRangeName = "WP";
        }
      }
      RooDataSet *reducedDataset =
          (RooDataSet *)dataset->reduce(RooFit::CutRange(comboRangeName));
      reducedDataset->SetName(constructName("dataset", index, templateIndex));
      std::cout << reducedDataset->GetName() << std::endl;
      workspace_->import(*reducedDataset);
      std::cout << "Entries in " << name << ": " << reducedDataset->numEntries()
                << std::endl;
      RooNLLVar *nll =
          new RooNLLVar(name, name, *workspace_->pdf(name_pdf), *reducedDataset,
                        RooFit::NumCPU(1), RooFit::Range("mTopFitRange"));
      nllSet.add(*nll);
      ++templateIndex;
    }
  }
  TString nllName = constructName("nllSim", index);
  RooAddition *nllSim = new RooAddition(nllName, nllName, nllSet);
  // workspace_->import(*nllSim);
  gStyle->SetOptTitle(1);
  // perform the fit

  RooMinuit minuit(*nllSim);
  minuit.setStrategy(2);
  minuit.migrad();
  minuit.improve();
  minuit.hesse();
  // minuit.minos();
  /*
   RooMinimizer minuit(*nllSim[h]);
   minuit.setStrategy(2);
   minuit.minimize("Genetic");
   //minuit.minimize("GSLMultiMin","BFGS2");
   //minuit.minimize("GSLMultiMin","SteepestDescent");
   */
  RooFitResult *result = minuit.save();
  workspace_->import(*result);
  return result;
}

void TemplateDerivation::plotResult(int templType, int comboType,
                                    int nComboTypes) {
  const int index = nComboTypes * templType + comboType;
  RooFitResult *result =
      (RooFitResult *)workspace_->genobj(constructName("nllSim", index));

  RooRealVar *var = (templType == 0) ? workspace_->var("topMass")
                                     : workspace_->var("meanWMass");
  RooPlot *frame = 0;
  TCanvas *canvas = new TCanvas("canvas", "canvas", 10, 10, 600, 600);
  TString figName = "";
  int templateIndex = 0;
  for (auto mass : massValues_) {
    for (auto jes : jesValues_) {
      // frame = var_mass[index]->frame();
      if (templType == 0 && comboType == 0)
        frame = var->frame(RooFit::Range(100., 215.));
      else if (templType == 0)
        frame = var->frame(RooFit::Range(100., maxMtop_));
      else if (templType == 1)
        frame = var->frame(RooFit::Range(50., 120.));

      RooAbsData *reducedDataset =
          workspace_->data(constructName("dataset", index, templateIndex));
      // if(! binnedTemplates) {
      reducedDataset->statOn(frame, RooFit::Layout(.6, .9, .9));
      reducedDataset->plotOn(frame);
      //}
      // else {
      //  hist[h][t]->statOn(frame, RooFit::Layout(.6, .9, .9));
      //  hist[h][t]->plotOn(frame);
      //}
      TString name_pdf = "sig_";
      name_pdf += index;
      name_pdf += "_" + constructTemplateName(mass, jes);
      std::cout << "PDF Name: " << name_pdf << std::endl;
      workspace_->pdf(name_pdf)->plotOn(frame, RooFit::FillColor(kGray),
                                        RooFit::VisualizeError(*result));
      workspace_->pdf(name_pdf)->plotOn(frame, RooFit::LineColor(kRed + 1));
      frame->SetMinimum(.0);
      std::string title =
          str(boost::format("MC for mass=%1$4.1f GeV and JSF= %2$3.2f") % mass %
              jes);
      frame->SetTitle(title.c_str());
      frame->Draw();
      figName = "template_";
      figName += templType;
      figName += "_";
      figName += comboType;
      figName += "_";
      figName += templateIndex;
      TString outDir = "plot/calibration";
      canvas->Print(outDir + "/" + figName + ".eps");
      canvas->Print(outDir + "/catalog.ps(");
      templateIndex++;
    }
  }
  // const int redPalette[5] = {kRed-2, kRed-6, kRed+1, kRed-8,
  // kRed+3};
  // const int redPalette[9] = {kRed-9, kRed-8, kRed-7, kRed-6,
  // kRed-5,
  // kRed-4, kRed-3, kRed-2, kRed-1};
  // show functions for different JES values
  // frame = var_mass[index]->frame();
  if (templType == 0 && comboType == 0)
    frame = var->frame(RooFit::Range(100., 215.));
  else if (templType == 0)
    frame = var->frame(RooFit::Range(100., maxMtop_));
  else if (templType == 1)
    frame = var->frame(RooFit::Range(50., 120.));
  // for(unsigned t=20; t<25; ++templateIndex) {
  for (auto jes : jesValues_) {
    TString name_pdf = "sig_";
    name_pdf += index;
    name_pdf += "_" + constructTemplateName(172.5, jes);
    // std::cout << name_pdf << std::endl;
    // workspace_->pdf(name_pdf)->plotOn(frame,
    // RooFit::FillColor(kGray), RooFit::VisualizeError(*result));
    workspace_->pdf(name_pdf)->plotOn(frame, RooFit::LineColor(kRed));
  }
  frame->SetMinimum(.0);
  frame->GetYaxis()->SetTitle("Probability");
  frame->SetTitle("JES = 0.96, 0.98, 1.00, 1.02, 1.04 (mTop = 172.5 GeV)");
  frame->Draw();
  figName = constructName("funcFamily_jes", templType, comboType);
  TString outDir = "plot/calibration";
  canvas->Print(outDir + "/" + figName + ".eps");
  canvas->Print(outDir + "/catalog.ps");
  // show functions for different mTop values
  // frame = var_mass[index]->frame();
  if (templType == 0 && comboType == 0)
    frame = var->frame(RooFit::Range(100., 215.));
  else if (templType == 0)
    frame = var->frame(RooFit::Range(100., maxMtop_));
  else if (templType == 1)
    frame = var->frame(RooFit::Range(50., 120.));
  // name_pdf = "sig_"; name_pdf += h; name_pdf += "_"+templ[2];
  // workspace[h]->pdf(name_pdf)->plotOn(frame,
  // RooFit::LineColor(kRed+1));
  for (auto mass : massValues_) {
    TString name_pdf = "sig_";
    name_pdf += index;
    name_pdf += "_" + constructTemplateName(mass, 1.00);
    // std::cout << name_pdf << std::endl;
    // workspace_->pdf(name_pdf)->plotOn(frame,
    // RooFit::FillColor(kGray), RooFit::VisualizeError(*result));
    workspace_->pdf(name_pdf)->plotOn(frame, RooFit::LineColor(kRed));
  }
  frame->SetMinimum(.0);
  frame->GetYaxis()->SetTitle("Probability");
  // frame->SetTitle("mTop = 161.5, 163.5, 166.5, 169.5, 172.5, 175.5,
  // 178.5, 181.5, 184.5 GeV (JES = 1.00)");
  frame->SetTitle(
      "mTop = 166.5, 169.5, 171.5, 172.5, 173.5, 175.5, 178.5 GeV "
      "(JES = "
      "1.00)");
  frame->Draw();
  figName = "funcFamily_mass_";
  figName += templType;
  figName += "_";
  figName += comboType;
  canvas->Print(outDir + "/" + figName + ".eps");
  canvas->Print(outDir + "/catalog.ps");
  gStyle->SetOptTitle(0);
  // plot different alpha_i as a function of JES
  char yTitle[99];
  for (unsigned i = 0; i < numVariables(constructName("alpha", index)); ++i) {
    TString jesName("JES_");
    jesName += constructTemplateName(massValues_[0], jesValues_[0]);
    frame = workspace_->var(jesName)->frame(RooFit::Range(0.9, 1.1));
    TString name = constructName("alpha", index, i);
    name += "_" + constructTemplateName(massValues_[0], jesValues_[0]);
    if (workspace_->function(name)) {
      workspace_->function(name)->plotOn(frame, RooFit::FillColor(kGray),
                                         RooFit::VisualizeError(*result));
      workspace_->function(name)->plotOn(frame, RooFit::Name("sig"),
                                         RooFit::LineColor(kRed + 1),
                                         RooFit::FillColor(kRed + 1));
    }
    frame->SetMinimum(.0);
    if (index == 0)
      sprintf(yTitle, "#alpha_{%i} (m_{top} = 172.5 GeV)", i);
    else
      sprintf(yTitle, "#alpha_{%i}", i);
    frame->GetYaxis()->SetTitle(yTitle);
    frame->GetXaxis()->SetTitle("JES");
    frame->Draw();
    double yLegHigh = 0.31;
    double yLegLow = yLegHigh - 2 * 0.06;
    TLegend legend(0.72, yLegLow, 0.92, yLegHigh);
    legend.SetBorderSize(0);
    legend.SetFillStyle(0);
    legend.SetTextFont(42);
    // legend.AddEntry("sig", "P ("+var_mass[h]->getTitle()+")", "F");
    legend.AddEntry("sig", "P (" + var->getTitle() + ")", "F");
    legend.Draw();
    figName = "alpha_";
    figName += i;
    figName += "_jes_";
    figName += templType;
    figName += "_";
    figName += comboType;
    canvas->Print(outDir + "/" + figName + ".eps");
    canvas->Print(outDir + "/catalog.ps");
    // plot different alpha_i as a function of mTop
    TString mtopName("mTop_");
    mtopName += constructTemplateName(massValues_[2], jesValues_[0]);
    frame = workspace_->var(mtopName)->frame(RooFit::Range(150., 200.));
    name = "alpha_";
    name += index;
    name += "_";
    name += i;
    name += "_" + constructTemplateName(massValues_[2], jesValues_[0]);
    if (workspace_->function(name)) {
      workspace_->function(name)->plotOn(frame, RooFit::FillColor(kGray),
                                         RooFit::VisualizeError(*result));
      workspace_->function(name)->plotOn(frame, RooFit::LineColor(kRed + 1));
    }
    frame->SetMinimum(.0);
    sprintf(yTitle, "#alpha_{%i} (JES = 1.00)", i);
    frame->GetYaxis()->SetTitle(yTitle);
    frame->GetXaxis()->SetTitle("m_{top} (GeV)");
    frame->Draw();
    legend.Draw();
    figName = "alpha_";
    figName += i;
    figName += "_mass_";
    figName += templType;
    figName += "_";
    figName += comboType;
    canvas->Print(outDir + "/" + figName + ".eps");
    canvas->Print(outDir + "/catalog.ps)");
  }
}

void TemplateDerivation::printResult(int templType, int comboType,
                                     int nComboTypes) {
  const int index = nComboTypes * templType + comboType;
  std::cout << "q[0] + q[1] * (x[0]-172.5) + (q[2] + q[3] * "
               "(x[0]-172.5)) * (x[1]-1.)" << std::endl;
  std::cout << "new: ";
  std::string vars_helper = "";
  std::stringstream vars(vars_helper,
                         std::stringstream::in | std::stringstream::out);
  for (unsigned int i = 0; i < numVariables(constructName("par", index)); ++i) {
    RooRealVar *curPar = workspace_->var(constructName("par", index, i));
    curPar->setConstant(kTRUE);
    vars << curPar->getVal();
    if (i != numVariables(constructName("par", index) - 1)) vars << "|";
  }
  if (templType == 0 && comboType == 1)
    vars << "|" << workspace_->var("ratio_1")->getVal();
  if (templType == 0 && comboType == 2)
    vars << "|" << workspace_->var("ratio_2")->getVal();
  std::cout << vars.str();
  std::cout << std::endl;

  std::ofstream myfile;
  TString outDir = "plot/calibration";
  myfile.open(outDir + "/variables.txt", std::ios::out | std::ios::app);
  myfile << "Template Type: " << templType << ", Combo Type: " << comboType
         << "\n";
  myfile << vars.str() << "\n";
  myfile.close();
}

void TemplateDerivation::rooFitTopMass() {
  const int nTemplTypes =
      2;  // number of different distributions, e.g., mTop & mW
  // const int nComboTypes = 3; // number of different permutation types,
  // e.g.,
  // correct, wrong, unmatched
  const int nComboTypes =
      (channelID_ == kAllJets) ? 2 : 3;  // number of different
                                         // permutation types, e.g.,
                                         // correct, rest
  RooRealVar *comboTypeVar =
      new RooRealVar("comboType", "comboType", -10., 20., "");
  // RooRealVar prob           = RooRealVar("prob" ,"P(#chi^{2})",  0.,
  // 1.,"");
  RooRealVar *MTOP = new RooRealVar("topMass", "m_{templateIndex}^{fit}", 100.,
                                    maxMtop_, "GeV");
  RooRealVar *meanMW =
      new RooRealVar("meanWMass", "m_{W}^{rec}", 50., 300., "GeV");
  RooRealVar *combinedWeight =
      new RooRealVar("combinedWeight", "weight", 0., 100., "");

  if (channelID_ == kAllJets) {
    comboTypeVar->setRange("CP", 0.9, 1.1);
    // comboTypeVar.setRange("UN",-9.9,-0.1);
    // comboTypeVar.setRange("WP",1.9,10.1);
    comboTypeVar->setRange("WP", 1.9, 20.1);
  }
  MTOP->setRange("mTopFitRange", 100., maxMtop_);
  meanMW->setRange("mWFitRange", 50., 300.);

  RooArgSet varSet = RooArgSet(*comboTypeVar, /*prob,*/ *MTOP, *meanMW,
                               *combinedWeight, "varSet");

  workspace_ = new RooWorkspace("workspaceMtop", "workspaceMtop");
  std::string categories = "";
  for (auto mass : massValues_) {
    for (auto jes : jesValues_) {
      categories += constructTemplateName(mass, jes);
      categories += ',';
    }
  }
  categories.pop_back();  // remove last ','
  workspace_->factory().createCategory("cat_templ", categories.c_str());
  /// create datasets and import for later use
  std::vector<RooDataSet *> dataset = createDataSets(varSet);

  RooRealVar *JES = new RooRealVar("JES", "JES", 1.0, 0.8, 1.2);
  RooRealVar *mTop =
      new RooRealVar("mTop", "mTop", 172.5, 100., maxMtop_, "GeV");
  JES->setConstant(kTRUE);
  mTop->setConstant(kTRUE);

  for (int templType = 0; templType < nTemplTypes; ++templType) {
    for (int comboType = 0; comboType < nComboTypes; ++comboType) {
      addTemplateFunction(templType, comboType, nComboTypes, mTop, JES);
      fitTemplate(templType, comboType, nComboTypes);
      plotResult(templType, comboType, nComboTypes);
      printResult(templType, comboType, nComboTypes);
    }
  }
  workspace_->writeToFile("RooWorkspace_TEST.root");
}

// because RooFit only likes plain trees with standard data types (int,
// float,
// double, ...)
// the original tree has to be adapted for the new content
TTree *TemplateDerivation::modifiedTree(TChain *tree, int minComboType,
                                        int maxComboType, bool isData) {
  tree->SetBranchStatus("*", 0);
  std::vector<std::string> vActiveBanches;
  boost::split(vActiveBanches, activeBranches_, boost::is_any_of("|"));
  for (const auto &branch : vActiveBanches) {
    tree->SetBranchStatus(branch.c_str(), 1);
  }

  TTreeFormula *f1 = new TTreeFormula("f1", fVar1_, tree);
  TTreeFormula *f2 = new TTreeFormula("f2", fVar2_, tree);
  TTreeFormula *f3 = new TTreeFormula("f3", fVar3_, tree);
  TTreeFormula *weight = new TTreeFormula("weight", fWeight_, tree);
  TTreeFormula *sel = new TTreeFormula("sel", selection_, tree);
  TTreeFormula *combo = new TTreeFormula("combo", "top.combinationType", tree);

  double topMass, meanWMass, combinedWeight, comboType;
  TTree *newTree = new TTree("tree", "tree");
  newTree->Branch("topMass", &topMass, "topMass/D");
  newTree->Branch("meanWMass", &meanWMass, "meanWMass/D");
  if (!isData) {
    newTree->Branch("combinedWeight", &combinedWeight, "combinedWeight/D");
    newTree->Branch("comboType", &comboType, "comboType/D");
  }

  int selected = 0;
  for (int i = 0;; ++i) {
    long entry = tree->LoadTree(i);
    if (entry < 0) break;
    if (entry == 0) {
      f1->UpdateFormulaLeaves();
      f2->UpdateFormulaLeaves();
      f3->UpdateFormulaLeaves();
      weight->UpdateFormulaLeaves();
      sel->UpdateFormulaLeaves();
      combo->UpdateFormulaLeaves();
    }
    if (!f1->GetNdata()) continue;
    if (!f2->GetNdata()) continue;
    if (!f3->GetNdata()) continue;
    if (!weight->GetNdata()) continue;
    if (!sel->GetNdata()) continue;
    if (!combo->GetNdata()) continue;
    int filledPermutations = 0;
    for (int j = 0, l = std::min(maxPermutations_, sel->GetNdata()); j < l;
         ++j) {
      if (!sel->EvalInstance(j)) continue;
      topMass = f1->EvalInstance(j);
      meanWMass = f2->EvalInstance(j);
      combinedWeight = f3->EvalInstance(j) * weight->EvalInstance(j);
      if (!isData) {
        comboType = combo->EvalInstance(j);
        if (comboType < 0) comboType = std::abs(comboType) + 10.;
      }
      newTree->Fill();
      filledPermutations++;
    }
    if (filledPermutations) ++selected;
  }

  std::cout << ": " << selected << " events" << std::endl;

  delete tree;
  delete f1;
  delete f2;
  delete f3;
  delete weight;
  delete sel;
  delete combo;

  return newTree;
}

void TemplateDerivation::UnknownChannelAbort() {
  std::cout << "Channel name *" << fChannel_
            << "* not know! Aborting program execution!" << std::endl;
  exit(1);
}

void TemplateDerivation::fillAlpha(std::vector<RooFormulaVar *> &alpha, int &h,
                                   RooArgSet argSet, std::string add) {
  TString formula_alpha = "@0+@1*(@4-172.5)+(@2+@3*(@4-172.5))*(@5-1.)";
  if (add.size() > 0) formula_alpha += add.c_str();
  std::cout << formula_alpha << std::endl;
  TString varName = constructName("alpha", h, alpha.size());
  alpha.push_back(new RooFormulaVar(varName, varName, formula_alpha, argSet));
}

TString TemplateDerivation::constructName(const std::string &name, int i) {
  TString temp(name);
  temp += "_";
  temp += i;
  return temp;
}

TString TemplateDerivation::constructName(const std::string &name, int i,
                                          int j) {
  TString temp(name);
  temp += "_";
  temp += i;
  temp += "_";
  temp += j;
  return temp;
}

unsigned int TemplateDerivation::numVariables(TString startName) const {
  return workspace_->allVars().selectByName(startName + "_*")->getSize();
}
