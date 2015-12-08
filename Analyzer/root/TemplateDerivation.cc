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
      doCalibration_(true),
      fitBackground_(false) {
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

  rooFitTopMass_();
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

std::vector<RooDataSet *> TemplateDerivation::createDataSets(
    std::vector<double> masses, std::vector<double> JESes,
    const RooArgSet &varSet) {
  TFile *tmpFile = TFile::Open("tmpFileRooFitTopMass.root", "RECREATE");
  std::vector<RooDataSet *> dataset;
  for (auto mass : masses) {
    for (auto jes : JESes) {
      std::string fileName = constructFileName(mass, jes);
      std::cout << "Creating RooDataSet for: " << fileName;
      TChain *chain = new TChain("analyzeKinFit/eventTree");
      std::cout << " (nFiles: " << chain->Add(samplePath_ + TString(fileName))
                << ") ";  // << std::endl;
      tmpFile->cd();
      TTree *tree = modifiedTree_(chain);  //, minComboType, maxComboType);
      int iTempl = dataset.size();
      TString name = "dataset_";
      name += iTempl;
      // needed to avoid a crash
      double co, /*pr,*/ mt, mw, cw;
      tree->SetBranchAddress("comboType", &co);
      // tree->SetBranchAddress("prob", &pr);
      tree->SetBranchAddress("topMass", &mt);
      tree->SetBranchAddress("meanWMass", &mw);
      tree->SetBranchAddress("combinedWeight", &cw);
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
    }
  }
  tmpFile->Close();
  return dataset;
}

void TemplateDerivation::rooFitTopMass_() {
  RooWorkspace *workspace[1];

  const int nTemplTypes =
      2;  // number of different distributions, e.g., mTop & mW
  // const int nComboTypes = 3; // number of different permutation types, e.g.,
  // correct, wrong, unmatched
  const int nComboTypes =
      (channelID_ == kAllJets)
          ? 2
          : 3;  // number of different permutation types, e.g., correct, rest
  const int nPDFs = nTemplTypes * nComboTypes;

  // bool binnedTemplates = false;

  std::vector<double> jesValues = {0.96, 0.98, 1.00, 1.02, 1.04};
  // const double mTop_templ[nMasses] = {161.5, 163.5, 166.5, 169.5, 172.5,
  // 175.5, 178.5, 181.5, 184.5};
  std::vector<double> massValues = {166.5, 169.5, 171.5, 172.5,
                                    173.5, 175.5, 178.5};

  unsigned int nTemplates = jesValues.size() * massValues.size();
  const unsigned iTemplateJES[] = {0, 1, 2, 3, 4, 0, 1, 2, 3, 4, 0, 1,
                                   2, 3, 4, 0, 1, 2, 3, 4, 0, 1, 2, 3,
                                   4, 0, 1, 2, 3, 4, 0, 1, 2, 3, 4};

  const unsigned iTemplateMass[] = {0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 2, 2,
                                    2, 2, 2, 3, 3, 3, 3, 3, 4, 4, 4, 4,
                                    4, 5, 5, 5, 5, 5, 6, 6, 6, 6, 6};

  TString templ[nTemplates] = {
      //"jes096mass1615","jes098mass1615","jes100mass1615","jes102mass1615","jes104mass1615",
      //"jes096mass1635","jes098mass1635","jes100mass1635","jes102mass1635","jes104mass1635",
      "jes096mass1665", "jes098mass1665", "jes100mass1665", "jes102mass1665",
      "jes104mass1665", "jes096mass1695", "jes098mass1695", "jes100mass1695",
      "jes102mass1695", "jes104mass1695", "jes096mass1715", "jes098mass1715",
      "jes100mass1715", "jes102mass1715", "jes104mass1715", "jes096mass1725",
      "jes098mass1725", "jes100mass1725", "jes102mass1725", "jes104mass1725",
      "jes096mass1735", "jes098mass1735", "jes100mass1735", "jes102mass1735",
      "jes104mass1735", "jes096mass1755", "jes098mass1755", "jes100mass1755",
      "jes102mass1755", "jes104mass1755", "jes096mass1785", "jes098mass1785",
      "jes100mass1785", "jes102mass1785", "jes104mass1785"  //,
                                                            //"jes096mass1815","jes098mass1815","jes100mass1815","jes102mass1815","jes104mass1815",
                                                            //"jes096mass1845","jes098mass1845","jes100mass1845","jes102mass1845","jes104mass1845"
  };
  TString templateSettings[nTemplates] = {
      //"MC for JES = 0.96, mTop = 161.5 GeV","MC for JES = 0.98, mTop = 161.5
      // GeV","MC for JES = 1.00, mTop = 161.5 GeV","MC for JES = 1.02, mTop =
      // 161.5 GeV","MC for JES = 1.04, mTop = 161.5 GeV",
      //"MC for JES = 0.96, mTop = 163.5 GeV","MC for JES = 0.98, mTop = 163.5
      // GeV","MC for JES = 1.00, mTop = 163.5 GeV","MC for JES = 1.02, mTop =
      // 163.5 GeV","MC for JES = 1.04, mTop = 163.5 GeV",
      "MC for JES = 0.96, mTop = 166.5 GeV",
      "MC for JES = 0.98, mTop = 166.5 GeV",
      "MC for JES = 1.00, mTop = 166.5 GeV",
      "MC for JES = 1.02, mTop = 166.5 GeV",
      "MC for JES = 1.04, mTop = 166.5 GeV",
      "MC for JES = 0.96, mTop = 169.5 GeV",
      "MC for JES = 0.98, mTop = 169.5 GeV",
      "MC for JES = 1.00, mTop = 169.5 GeV",
      "MC for JES = 1.02, mTop = 169.5 GeV",
      "MC for JES = 1.04, mTop = 169.5 GeV",
      "MC for JES = 0.96, mTop = 171.5 GeV",
      "MC for JES = 0.98, mTop = 171.5 GeV",
      "MC for JES = 1.00, mTop = 171.5 GeV",
      "MC for JES = 1.02, mTop = 171.5 GeV",
      "MC for JES = 1.04, mTop = 171.5 GeV",
      "MC for JES = 0.96, mTop = 172.5 GeV",
      "MC for JES = 0.98, mTop = 172.5 GeV",
      "MC for JES = 1.00, mTop = 172.5 GeV",
      "MC for JES = 1.02, mTop = 172.5 GeV",
      "MC for JES = 1.04, mTop = 172.5 GeV",
      "MC for JES = 0.96, mTop = 173.5 GeV",
      "MC for JES = 0.98, mTop = 173.5 GeV",
      "MC for JES = 1.00, mTop = 173.5 GeV",
      "MC for JES = 1.02, mTop = 173.5 GeV",
      "MC for JES = 1.04, mTop = 173.5 GeV",
      "MC for JES = 0.96, mTop = 175.5 GeV",
      "MC for JES = 0.98, mTop = 175.5 GeV",
      "MC for JES = 1.00, mTop = 175.5 GeV",
      "MC for JES = 1.02, mTop = 175.5 GeV",
      "MC for JES = 1.04, mTop = 175.5 GeV",
      "MC for JES = 0.96, mTop = 178.5 GeV",
      "MC for JES = 0.98, mTop = 178.5 GeV",
      "MC for JES = 1.00, mTop = 178.5 GeV",
      "MC for JES = 1.02, mTop = 178.5 GeV",
      "MC for JES = 1.04, mTop = 178.5 GeV"  //,
      //"MC for JES = 0.96, mTop = 181.5 GeV","MC for JES = 0.98, mTop = 181.5
      // GeV","MC for JES = 1.00, mTop = 181.5 GeV","MC for JES = 1.02, mTop =
      // 181.5 GeV","MC for JES = 1.04, mTop = 181.5 GeV",
      //"MC for JES = 0.96, mTop = 184.5 GeV","MC for JES = 0.98, mTop = 184.5
      // GeV","MC for JES = 1.00, mTop = 184.5 GeV","MC for JES = 1.02, mTop =
      // 184.5 GeV","MC for JES = 1.04, mTop = 184.5 GeV"
  };

  // RooDataHist* hist[1][nTemplates];

  RooRealVar comboTypeVar = RooRealVar("comboType", "comboType", -10., 20., "");
  // RooRealVar prob           = RooRealVar("prob"          ,"P(#chi^{2})",  0.,
  // 1.,"");
  RooRealVar MTOP = RooRealVar("topMass", "m_{t}^{fit}", 100., maxMtop_, "GeV");
  RooRealVar meanMW = RooRealVar("meanWMass", "m_{W}^{rec}", 50., 300., "GeV");
  RooRealVar combinedWeight =
      RooRealVar("combinedWeight", "weight", 0., 100., "");

  if (channelID_ == kAllJets) {
    comboTypeVar.setRange("R1", 0.9, 1.1);
    // comboTypeVar.setRange("R4",-9.9,-0.1);
    // comboTypeVar.setRange("R5",1.9,10.1);
    comboTypeVar.setRange("R2", 1.9, 20.1);
  }
  MTOP.setRange("mTopFitRange", 100., maxMtop_);
  meanMW.setRange("mWFitRange", 50., 300.);

  RooArgSet varSet =
      RooArgSet(comboTypeVar, /*prob,*/ MTOP, meanMW, combinedWeight, "varSet");

  /// create datasets for later use
  std::vector<RooDataSet *> dataset =
      createDataSets(massValues, jesValues, varSet);
  RooDataSet *reducedDataset[nTemplates];

  // RooAddPdf *topAdd = 0;
  // RooAddPdf   *wAdd = 0;

  TString name = "";
  TString outDir = "plot/calibration";

  if (doCalibration_) {
    RooCategory *cat_templ = new RooCategory("cat_templ", "cat_templ");
    for (unsigned t = 0; t < nTemplates; ++t) cat_templ->defineType(templ[t]);

    workspace[0] = new RooWorkspace("workspaceMtop", "workspaceMtop");
    workspace[0]->import(*cat_templ);

    RooSimWSTool *simWST[nPDFs];
    RooSimultaneous *sim[nPDFs];
    RooNLLVar *nll[nPDFs][nTemplates];
    RooAddition *nllSim[nPDFs];

    RooRealVar *JES = new RooRealVar("JES", "JES", 1.0, 0.8, 1.2);
    RooRealVar *mTop =
        new RooRealVar("mTop", "mTop", 172.5, 100., maxMtop_, "GeV");
    JES->setConstant(kTRUE);
    mTop->setConstant(kTRUE);

    for (int templType = 0; templType < nTemplTypes; ++templType) {
      for (int comboType = 0; comboType < nComboTypes; ++comboType) {
        int h = nComboTypes * templType + comboType;

        RooRealVar MassOffset = RooRealVar("MassOffset", "MassOffset", 0.);
        MassOffset.setConstant(kTRUE);
        RooRealVar MassSlopeMass =
            RooRealVar("MassSlopeMass", "MassSlopeMass", 0.);
        MassSlopeMass.setConstant(kTRUE);
        RooRealVar MassSlopeJES =
            RooRealVar("MassSlopeJES", "MassSlopeJES", 0.);
        MassSlopeJES.setConstant(kTRUE);
        RooRealVar MassSlopeMassJES =
            RooRealVar("MassSlopeMassJES", "MassSlopeMassJES", 0.);
        MassSlopeMassJES.setConstant(kTRUE);
        RooRealVar JESOffset = RooRealVar("JESOffset", "JESOffset", 0.);
        JESOffset.setConstant(kTRUE);
        RooRealVar JESSlopeMass =
            RooRealVar("JESSlopeMass", "JESSlopeMass", 0.);
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
        RooRealVar MassSlopeJES2 =
            RooRealVar("MassSlopeJES2", "MassSlopeJES2", 0.);
        MassSlopeJES2.setConstant(kTRUE);
        RooRealVar MassSlopeMassJES2 =
            RooRealVar("MassSlopeMassJES2", "MassSlopeMassJES2", 0.);
        MassSlopeMassJES2.setConstant(kTRUE);
        RooRealVar JESOffset2 = RooRealVar("JESOffset2", "JESOffset2", 0.);
        JESOffset2.setConstant(kTRUE);
        RooRealVar JESSlopeMass2 =
            RooRealVar("JESSlopeMass2", "JESSlopeMass2", 0.);
        JESSlopeMass2.setConstant(kTRUE);
        RooRealVar JESSlopeJES2 =
            RooRealVar("JESSlopeJES2", "JESSlopeJES2", 0.);
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
        name = "mTop_intermediate_";
        name += h;
        RooFormulaVar mTop_intermediate =
            RooFormulaVar(name, name, formula_mTop,
                          RooArgSet(*JES, *mTop, MassOffset, MassSlopeMass,
                                    MassSlopeJES, MassSlopeMassJES));
        name = "mJES_intermediate_";
        name += h;
        RooFormulaVar JES_intermediate =
            RooFormulaVar(name, name, formula_JES,
                          RooArgSet(*JES, *mTop, JESOffset, JESSlopeMass,
                                    JESSlopeJES, JESSlopeMassJES));
        name = "mTop_corrected_";
        name += h;
        RooFormulaVar mTop_corrected = RooFormulaVar(
            name, name, formula_mTop,
            RooArgSet(JES_intermediate, mTop_intermediate, MassOffset2,
                      MassSlopeMass2, MassSlopeJES2, MassSlopeMassJES2));
        name = "mJES_corrected_";
        name += h;
        RooFormulaVar JES_corrected = RooFormulaVar(
            name, name, formula_JES,
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
            comboRangeName = "R1";
          }
          // else if(comboType == 1){ comboRangeName = "R4"; }
          // else if(comboType == 2){ comboRangeName = "R5"; }
          else if (comboType == 1) {
            comboRangeName = "R2";
          }

          if (templType == 0) {
            if (comboType == 0) {  // mTop, correct
              iniPar = {172.399, 0.989847,  81.225,  0.758306,
                        7.92651, 0.0730216, 5.03205, 0.0241954};
            } else if (comboType == 1) {  // mTop, wrong
              iniPar = {190.009, 0.23107,    123.714,  -2.41287,
                        173.325, 1.11063,    80.0112,  1.32969,
                        26.2827, -0.0172727, 35.0065,  -0.442854,
                        10.0152, 0.0656903,  -4.72985, -0.108491};
            } else if (comboType == 2) {  // mTop, unmatched
              iniPar = {178.883, 0.306507,  40.7573, 0.445759,
                        174.681, 1.06829,   83.8601, 0.280563,
                        22.065,  0.0542991, 4.93823, 0.3491,
                        9.74,    0.0548014, 2.31277, -0.181041};
            }
          } else if (templType == 1) {
            if (comboType == 0) {  // mW, correct
              iniPar = {84.4535, -0.0154884, 91.1495,  -0.00943573,
                        5.20726, -0.0115877, 23.4646,  -0.0249863,
                        6.75636, 0.0019251,  -19.8341, 0.0724597};
            } else if (comboType == 1) {  // mW, wrong
              iniPar = {87.2093, -0.00584533, 20.9648,  -0.0912001,
                        3.91928, 0.00521933,  3.8585,   -0.28899,
                        5.9744,  0.00351353,  -5.70177, -0.231188,
                        3.80803, 0.00289055,  -11.4316, -0.264811,
                        4.89211, 4.89211,     0.25,     0.5};
            } else if (comboType == 2) {  // mW, unmatched
              iniPar = {86.6185, -0.0184471, 32.1242,  0.0815632,
                        6.48169, -0.0108763, 7.37265,  0.0139849,
                        6.48169, -0.0108763, 7.37265,  0.0139849,
                        7.91363, 0.00815331, -17.3423, -0.00100213};
            }
          }

          for (unsigned p = 0; p < iniPar.size(); ++p) {
            // std::cout << iniPar[p] << ", ";
            name = "par_";
            name += h;
            name += "_";
            name += p;
            RooRealVar *myVar = new RooRealVar(name, name, iniPar[p]);
            myVar->setConstant(kFALSE);
            /*
             // for WP JES ONLY
             if(templType == 1 && comboType == 1){
             if(p ==  0) myVar->setRange(86.5,88);
             //if(p ==  1) myVar->setRange(-0.01,0.01);
             if(p ==  1) myVar->setVal(0.0);
             if(p ==  1) myVar->setConstant();
             if(p ==  2) myVar->setRange(18,24);
             //if(p ==  3) myVar->setRange(-0.5,0.5);
             if(p ==  3) myVar->setVal(0.0);
             if(p ==  3) myVar->setConstant();

             if(p ==  4) myVar->setRange(2,6);
             //if(p ==  5) myVar->setRange(-0.01,0.01);
             if(p ==  5) myVar->setVal(0.0);
             if(p ==  5) myVar->setConstant();
             if(p ==  6) myVar->setRange(0,10);
             //if(p ==  7) myVar->setRange(-0.5,0.5);
             if(p ==  7) myVar->setVal(0.0);
             if(p ==  7) myVar->setConstant();

             if(p ==  8) myVar->setRange(4,9);
             //if(p ==  9) myVar->setRange(-0.01,0.01);
             if(p ==  9) myVar->setVal(0.0);
             if(p ==  9) myVar->setConstant();
             if(p == 10) myVar->setRange(-10,0);
             //if(p == 11) myVar->setRange(-0.5,0.5);
             if(p == 11) myVar->setVal(0.0);
             if(p == 11) myVar->setConstant();

             if(p == 12) myVar->setRange(2,6);
             //if(p == 13) myVar->setRange(-0.01,0.01);
             if(p == 13) myVar->setVal(0.0);
             if(p == 13) myVar->setConstant();
             if(p == 14) myVar->setRange(-20,-5);
             //if(p == 15) myVar->setRange(-0.5,0.5);
             if(p == 15) myVar->setVal(0.0);
             if(p == 15) myVar->setConstant();

             if(p == 16) myVar->setRange(2,9);
             if(p == 17) myVar->setRange(2,9);

             if(p == 18) myVar->setRange(0,0.25);
             if(p == 19) myVar->setRange(0.25,0.75);
             }
             // for CP JES ONLY
             if(templType == 1 && comboType == 0){
             */
            // if(p == 16) myVar->setVal(7.5);
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
          name += h;
          TString varName;
          if (templType == 0) {
            var = topMass;
            if (comboType == 0) {
              fillAlpha(alpha, h, RooArgSet(*par[0], *par[1], *par[2], *par[3],
                                            mTop_corrected, JES_corrected));
              fillAlpha(alpha, h, RooArgSet(*par[4], *par[5], *par[6], *par[7],
                                            mTop_corrected, JES_corrected));
              varName = "width_";
              varName += h;
              RooRealVar *width = new RooRealVar(varName, varName, 2.0);
              width->setConstant(kTRUE);
              RooVoigtian voigt =
                  RooVoigtian(name, name, *var, *alpha[0], *width, *alpha[1]);
              workspace[0]->import(voigt);
            } else if (comboType == 1 || comboType == 2) {
              varName = "ratio_";
              varName += h;
              RooRealVar *ratio = new RooRealVar(varName, varName, 0.0, 1.0);
              if (comboType == 1)
                ratio->setVal(0.700556);  // 0.648584); //0.662784);
                                          // //0.820546);
              else if (comboType == 2)
                ratio->setVal(0.685704);  // 0.774588); //0.758690);
              // ratio->setConstant(kTRUE);
              ratio->setConstant(kFALSE);
              fillAlpha(alpha, h, RooArgSet(*par[0], *par[1], *par[2], *par[3],
                                            mTop_corrected, JES_corrected));
              fillAlpha(alpha, h,
                        RooArgSet(*par[8], *par[9], *par[10], *par[11],
                                  mTop_corrected, JES_corrected));
              varName = "landau_";
              varName += h;
              RooLandau *landau =
                  new RooLandau(varName, varName, *var, *alpha[0], *alpha[1]);
              fillAlpha(alpha, h, RooArgSet(*par[4], *par[5], *par[6], *par[7],
                                            mTop_corrected, JES_corrected));
              fillAlpha(alpha, h,
                        RooArgSet(*par[12], *par[13], *par[14], *par[15],
                                  mTop_corrected, JES_corrected));
              varName = "width_";
              varName += h;
              RooRealVar *width = new RooRealVar(varName, varName, 2.0);
              width->setConstant(kTRUE);
              varName = "voigt_";
              varName += h;
              // RooVoigtian *voigt = new RooVoigtian(varName, varName, *var,
              // *alpha[2], *width, *alpha[3]);
              RooGaussian *voigt =
                  new RooGaussian(varName, varName, *var, *alpha[2], *alpha[3]);
              RooAddPdf *add =
                  new RooAddPdf(name, name, *landau, *voigt, *ratio);
              workspace[0]->import(*add);
            }
          } else if (templType == 1) {
            var = meanWMass;
            if (comboType == 0) {  // mW, correct
              fillAlpha(alpha, h, RooArgSet(*par[0], *par[1], *par[2], *par[3],
                                            mTop_corrected, JES_corrected));
              fillAlpha(alpha, h, RooArgSet(*par[4], *par[5], *par[6], *par[7],
                                            mTop_corrected, JES_corrected));
              fillAlpha(alpha, h,
                        RooArgSet(*par[8], *par[9], *par[10], *par[11],
                                  mTop_corrected, JES_corrected));
              RooBifurGauss *asymGaus = new RooBifurGauss(
                  name, name, *var, *alpha[0], *alpha[1], *alpha[2]);
              workspace[0]->import(*asymGaus);
            } else {
              // fillAlpha(alpha, h, RooArgSet(*par[ 0], *par[ 1], *par[ 2],
              // *par[ 3], mTop_corrected, JES_corrected), "-7.5");
              fillAlpha(alpha, h,
                        RooArgSet(*par[0], *par[1], *par[2], *par[3],
                                  mTop_corrected, JES_corrected, *par[16]),
                        "-@6");
              fillAlpha(alpha, h, RooArgSet(*par[0], *par[1], *par[2], *par[3],
                                            mTop_corrected, JES_corrected));
              // fillAlpha(alpha, h, RooArgSet(*par[ 0], *par[ 1], *par[ 2],
              // *par[ 3], mTop_corrected, JES_corrected), "+7.5");
              fillAlpha(alpha, h,
                        RooArgSet(*par[0], *par[1], *par[2], *par[3],
                                  mTop_corrected, JES_corrected, *par[17]),
                        "+@6");
              fillAlpha(alpha, h, RooArgSet(*par[4], *par[5], *par[6], *par[7],
                                            mTop_corrected, JES_corrected));
              fillAlpha(alpha, h,
                        RooArgSet(*par[8], *par[9], *par[10], *par[11],
                                  mTop_corrected, JES_corrected));
              fillAlpha(alpha, h,
                        RooArgSet(*par[12], *par[13], *par[14], *par[15],
                                  mTop_corrected, JES_corrected));
              varName = "gaus1_";
              varName += h;
              // RooRealVar *ratio1 = new RooRealVar("ratioGaus1","",0.25);
              // RooRealVar *ratio2 = new RooRealVar("ratioGaus2","",0.5 );
              // RooRealVar *ratio3 = new RooRealVar("ratioGaus3","",0.25);
              RooGaussian *gaus1 =
                  new RooGaussian(varName, varName, *var, *alpha[0], *alpha[3]);
              varName = "gaus2_";
              varName += h;
              RooGaussian *gaus2 =
                  new RooGaussian(varName, varName, *var, *alpha[1], *alpha[4]);
              varName = "gaus3_";
              varName += h;
              RooGaussian *gaus3 =
                  new RooGaussian(varName, varName, *var, *alpha[2], *alpha[5]);
              RooAddPdf *add = new RooAddPdf(
                  name, name, RooArgList(*gaus1, *gaus2, *gaus3),
                  RooArgList(
                      *par[18],
                      *par[19]) /*RooArgList(*ratio1, *ratio2, *ratio3)*/);
              workspace[0]->import(*add);
            }
          }
        }
        simWST[h] = new RooSimWSTool(*workspace[0]);
        name = "sim_";
        name += h;
        TString name_pdf = "sig_";
        name_pdf += h;
        sim[h] = simWST[h]->build(name, name_pdf,
                                  RooFit::SplitParam("JES", "cat_templ"),
                                  RooFit::SplitParam("mTop", "cat_templ"));
        for (unsigned t = 0; t < nTemplates; ++t) {
          workspace[0]->var("JES_" + templ[t])->setVal(
              jesValues[iTemplateJES[t]]);
          workspace[0]->var("mTop_" + templ[t])->setVal(
              massValues[iTemplateMass[t]]);
          name = "nll_";
          name += h;
          name += "_";
          name += t;
          name_pdf = "sig_";
          name_pdf += h;
          name_pdf += "_" + templ[t];
          // if(binnedTemplates)
          //  nll[h][t] = new RooNLLVar(name, name,
          //  *workspace[0]->pdf(name_pdf), *hist[h][t], RooFit::NumCPU(1),
          //  RooFit::Range(100., 550.));
          // else
          reducedDataset[t] = (RooDataSet *)dataset[t]->reduce(
              RooFit::CutRange(comboRangeName));
          TString name_dataset = "dataset_";
          name_dataset += h;
          name_dataset += "_";
          name_dataset += t;
          reducedDataset[t]->SetName(name_dataset);
          std::cout << reducedDataset[t]->GetName() << std::endl;
          workspace[0]->import(*reducedDataset[t]);
          std::cout << "Entries in " << name << ": "
                    << reducedDataset[t]->numEntries() << std::endl;
          nll[h][t] = new RooNLLVar(name, name, *workspace[0]->pdf(name_pdf),
                                    *reducedDataset[t], RooFit::NumCPU(1),
                                    RooFit::Range("mTopFitRange"));
        }
        name = "nllSim_";
        name += h;
        RooArgSet argSet;
        for (unsigned t = 0; t < nTemplates; ++t) argSet.add(*nll[h][t]);
        nllSim[h] = new RooAddition(name, name, argSet);
        RooFitResult *result;
        gStyle->SetOptTitle(1);
        // perform the fit

        RooMinuit minuit(*nllSim[h]);
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
        result = minuit.save();
        workspace[0]->import(*result);

        /// control plots

        RooPlot *frame = 0;
        TCanvas *canvas = new TCanvas("canvas", "canvas", 10, 10, 600, 600);
        TString figName = "";
        for (unsigned t = 0; t < nTemplates; ++t) {
          // frame = var_mass[h]->frame();
          if (templType == 0 && comboType == 0)
            frame = var->frame(RooFit::Range(100., 215.));
          else if (templType == 0)
            frame = var->frame(RooFit::Range(100., maxMtop_));
          else if (templType == 1)
            frame = var->frame(RooFit::Range(50., 120.));
          // if(! binnedTemplates) {
          reducedDataset[t]->statOn(frame, RooFit::Layout(.6, .9, .9));
          reducedDataset[t]->plotOn(frame);
          //}
          // else {
          //  hist[h][t]->statOn(frame, RooFit::Layout(.6, .9, .9));
          //  hist[h][t]->plotOn(frame);
          //}
          name_pdf = "sig_";
          name_pdf += h;
          name_pdf += "_" + templ[t];
          std::cout << "PDF Name: " << name_pdf << std::endl;
          workspace[0]->pdf(name_pdf)->plotOn(frame, RooFit::FillColor(kGray),
                                              RooFit::VisualizeError(*result));
          workspace[0]->pdf(name_pdf)->plotOn(frame,
                                              RooFit::LineColor(kRed + 1));
          frame->SetMinimum(.0);
          frame->SetTitle(templateSettings[t]);
          frame->Draw();
          figName = "template_";
          figName += templType;
          figName += "_";
          figName += comboType;
          figName += "_";
          figName += t;
          canvas->Print(outDir + "/" + figName + ".eps");
          canvas->Print(outDir + "/catalog.ps(");
        }
        // const int redPalette[5] = {kRed-2, kRed-6, kRed+1, kRed-8, kRed+3};
        // const int redPalette[9] = {kRed-9, kRed-8, kRed-7, kRed-6, kRed-5,
        // kRed-4, kRed-3, kRed-2, kRed-1};
        // show functions for different JES values
        // frame = var_mass[h]->frame();
        if (templType == 0 && comboType == 0)
          frame = var->frame(RooFit::Range(100., 215.));
        else if (templType == 0)
          frame = var->frame(RooFit::Range(100., maxMtop_));
        else if (templType == 1)
          frame = var->frame(RooFit::Range(50., 120.));
        // for(unsigned t=20; t<25; ++t) {
        for (unsigned t = 15; t < 20; ++t) {
          name_pdf = "sig_";
          name_pdf += h;
          name_pdf += "_" + templ[t];
          // std::cout << name_pdf << std::endl;
          // workspace[0]->pdf(name_pdf)->plotOn(frame,
          // RooFit::FillColor(kGray), RooFit::VisualizeError(*result));
          workspace[0]->pdf(name_pdf)->plotOn(frame, RooFit::LineColor(kRed));
        }
        frame->SetMinimum(.0);
        frame->GetYaxis()->SetTitle("Probability");
        frame->SetTitle(
            "JES = 0.96, 0.98, 1.00, 1.02, 1.04 (mTop = 172.5 GeV)");
        frame->Draw();
        figName = "funcFamily_jes_";
        figName += templType;
        figName += "_";
        figName += comboType;
        canvas->Print(outDir + "/" + figName + ".eps");
        canvas->Print(outDir + "/catalog.ps");
        // show functions for different mTop values
        // frame = var_mass[h]->frame();
        if (templType == 0 && comboType == 0)
          frame = var->frame(RooFit::Range(100., 215.));
        else if (templType == 0)
          frame = var->frame(RooFit::Range(100., maxMtop_));
        else if (templType == 1)
          frame = var->frame(RooFit::Range(50., 120.));
        // name_pdf = "sig_"; name_pdf += h; name_pdf += "_"+templ[2];
        // workspace[h]->pdf(name_pdf)->plotOn(frame,
        // RooFit::LineColor(kRed+1));
        for (unsigned t = 1, i = 0; t < nTemplates; t += 5, ++i) {
          name_pdf = "sig_";
          name_pdf += h;
          name_pdf += "_" + templ[t];
          // std::cout << name_pdf << std::endl;
          // workspace[0]->pdf(name_pdf)->plotOn(frame,
          // RooFit::FillColor(kGray), RooFit::VisualizeError(*result));
          workspace[0]->pdf(name_pdf)->plotOn(frame, RooFit::LineColor(kRed));
        }
        frame->SetMinimum(.0);
        frame->GetYaxis()->SetTitle("Probability");
        // frame->SetTitle("mTop = 161.5, 163.5, 166.5, 169.5, 172.5, 175.5,
        // 178.5, 181.5, 184.5 GeV (JES = 1.00)");
        frame->SetTitle(
            "mTop = 166.5, 169.5, 171.5, 172.5, 173.5, 175.5, 178.5 GeV (JES = "
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
        for (unsigned i = 0; i < alpha.size(); ++i) {
          frame = workspace[0]->var("JES_" + templ[0])->frame(
              RooFit::Range(0.9, 1.1));
          name = "alpha_";
          name += h;
          name += "_";
          name += i;
          if (workspace[0]->function(name + "_" + templ[0])) {
            workspace[0]->function(name + "_" + templ[0])->plotOn(
                frame, RooFit::FillColor(kGray),
                RooFit::VisualizeError(*result));
            workspace[0]->function(name + "_" + templ[0])->plotOn(
                frame, RooFit::Name("sig"), RooFit::LineColor(kRed + 1),
                RooFit::FillColor(kRed + 1));
          }
          frame->SetMinimum(.0);
          if (h == 0)
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
          frame = workspace[0]->var("mTop_" + templ[2])->frame(
              RooFit::Range(150., 200.));
          name = "alpha_";
          name += h;
          name += "_";
          name += i;
          if (workspace[0]->function(name + "_" + templ[2])) {
            workspace[0]->function(name + "_" + templ[2])->plotOn(
                frame, RooFit::FillColor(kGray),
                RooFit::VisualizeError(*result));
            workspace[0]->function(name + "_" + templ[2])->plotOn(
                frame, RooFit::LineColor(kRed + 1));
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

        //////////////////////////////////////////////////////////////////////////////////////////////////////

        std::cout << "q[0] + q[1] * (x[0]-172.5) + (q[2] + q[3] * "
                     "(x[0]-172.5)) * (x[1]-1.)" << std::endl;
        std::cout << "new: ";
        std::string vars_helper = "";
        std::stringstream vars(vars_helper,
                               std::stringstream::in | std::stringstream::out);
        for (unsigned int i = 0; i < par.size(); ++i) {
          TString varName = "par_";
          varName += h;
          varName += "_";
          varName += i;
          RooRealVar *curPar = workspace[0]->var(varName);
          curPar->setConstant(kTRUE);
          vars << curPar->getVal();
          if (i != par.size() - 1) vars << "|";
        }
        if (templType == 0 && comboType == 1)
          vars << "|" << workspace[0]->var("ratio_1")->getVal();
        if (templType == 0 && comboType == 2)
          vars << "|" << workspace[0]->var("ratio_2")->getVal();
        std::cout << vars.str();
        std::cout << std::endl;

        std::ofstream myfile;
        myfile.open(outDir + "/variables.txt", std::ios::out | std::ios::app);
        myfile << "Template Type: " << templType
               << ", Combo Type: " << comboType << "\n";
        myfile << vars.str() << "\n";
        myfile.close();
        std::cout << "old: ";
        for (unsigned int i = 0; i < iniPar.size(); ++i) {
          std::cout << iniPar[i];
          if (i != iniPar.size() - 1) std::cout << "|";
        }
        std::cout << std::endl;
        //}
      }
    }
    workspace[0]->writeToFile("RooWorkspace_TEST.root");
  } else {
    TFile *calibrationFile = TFile::Open("RooWorkspace_TEST.root");
    workspace[0] =
        (RooWorkspace *)calibrationFile->Get("workspaceMtop")->Clone();
  }

  /*

   RooRealVar *fCP = new RooRealVar("fCP", "fCP", 0.254136);
   fCP->setConstant(kTRUE);
   RooRealVar *fWP = new RooRealVar("fWP", "fWP",
   6.90155e-06+0.000212875+0.00280065+0.0244401+0.0921562);
   fWP->setConstant(kTRUE);
   RooRealVar *fUN = new RooRealVar("fUN", "fUN",
   0.00541135+3.66988e-05+0.0570029+0.000743396+0.563053);
   fUN->setConstant(kTRUE);
   RooArgSet permutationFractions =
   RooArgSet(*fCP,*fWP,*fUN,"permutationFractions");

   RooAbsPdf *topCP = workspace[0]->pdf("sig_0");
   RooAbsPdf *topWP = workspace[0]->pdf("sig_1");
   RooAbsPdf *topUN = workspace[0]->pdf("sig_2");
   topCP->setNormRange("mTopFitRange");
   topWP->setNormRange("mTopFitRange");
   topUN->setNormRange("mTopFitRange");
   RooArgSet topPDFs =  RooArgSet(*topCP,*topWP,*topUN,"topPDFs");
   topAdd = new RooAddPdf("mTopPDF", "mTopPDF", topPDFs, permutationFractions);
   topAdd->setNormRange("mTopFitRange");
   //workspace[0]->import(*topAdd);

   RooAbsPdf *wCP = workspace[0]->pdf("sig_3");
   RooAbsPdf *wWP = workspace[0]->pdf("sig_4");
   RooAbsPdf *wUN = workspace[0]->pdf("sig_5");
   wCP->setNormRange("mTopFitRange");
   wWP->setNormRange("mTopFitRange");
   wUN->setNormRange("mTopFitRange");
   RooArgSet wPDFs =  RooArgSet(*wCP,*wWP,*wUN,"wPDFs");
   wAdd = new RooAddPdf("mWPDF", "mWPDF", wPDFs, permutationFractions);
   wAdd->setNormRange("mTopFitRange");
   //workspace[0]->import(*wAdd);

   */

  //  if(doMeasurement_){
  //    TFile * workspaceFile = TFile::Open("RooWorkspace_TEST.root");
  //    workspace[0] = (RooWorkspace*)workspaceFile->Get("workspaceMtop");
  //
  //    //MassOffset       = *workspace[0]->var("MassOffset"      ); MassOffset
  //    .setVal( 1.05709e-02);
  //    //MassSlopeMass    = *workspace[0]->var("MassSlopeMass"   );
  //    MassSlopeMass   .setVal( 6.41444e-02);
  //    //MassSlopeJES     = *workspace[0]->var("MassSlopeJES"    );
  //    MassSlopeJES    .setVal(-4.61594e+00);
  //    //MassSlopeMassJES = *workspace[0]->var("MassSlopeMassJES");
  //    MassSlopeMassJES.setVal( 0.         );
  //    //JESOffset        = *workspace[0]->var("JESOffset"       ); JESOffset
  //    .setVal(-3.76801e-03);
  //    //JESSlopeMass     = *workspace[0]->var("JESSlopeMass"    );
  //    JESSlopeMass    .setVal(-1.50418e-04);
  //    //JESSlopeJES      = *workspace[0]->var("JESSlopeJES"     ); JESSlopeJES
  //    .setVal( 4.60147e-02);
  //    //JESSlopeMassJES  = *workspace[0]->var("JESSlopeMassJES" );
  //    JESSlopeMassJES .setVal( 0.         );
  //
  //    //workspace[0]->var("MassOffset"      )->setVal( 1.05709e-02);
  //    //workspace[0]->var("MassSlopeMass"   )->setVal( 6.41444e-02);
  //    //workspace[0]->var("MassSlopeJES"    )->setVal(-4.61594e+00);
  //    //workspace[0]->var("MassSlopeMassJES")->setVal( 0.         );
  //    //workspace[0]->var("JESOffset"       )->setVal(-3.76801e-03);
  //    //workspace[0]->var("JESSlopeMass"    )->setVal(-1.50418e-04);
  //    //workspace[0]->var("JESSlopeJES"     )->setVal( 4.60147e-02);
  //    //workspace[0]->var("JESSlopeMassJES" )->setVal( 0.         );
  //
  //    //workspace[0]->var("MassOffset"      )->setVal( -0.71);
  //    //workspace[0]->var("MassSlopeMass"   )->setVal(  0.  );
  //    //workspace[0]->var("MassSlopeJES"    )->setVal( 12.59);
  //    //workspace[0]->var("MassSlopeMassJES")->setVal(  0.  );
  //    //workspace[0]->var("JESOffset"       )->setVal( 0.0199);
  //    //workspace[0]->var("JESSlopeMass"    )->setVal( 0.    );
  //    //workspace[0]->var("JESSlopeJES"     )->setVal( 0.8281);
  //    //workspace[0]->var("JESSlopeMassJES" )->setVal( 0.    );
  //
  //    //mTop_corrected   =
  //    *(RooFormulaVar*)workspace[0]->var("mTop_corrected");
  //    //JES_corrected    = *(RooFormulaVar*)workspace[0]->var(
  //    "JES_corrected");
  //
  //    topAdd = (RooAddPdf*)workspace[0]->pdf("mTopPDF");
  //    wAdd   = (RooAddPdf*)workspace[0]->pdf("mWPDF"  );
  //
  //    varSet  = RooArgSet(/*prob,*///MTOP,meanMW,"varSet");
  //
  //    TString fileName = "writeFullHadTree_data_2011.root";
  //    //TString fileName = "Z2_F11_172_5_sig.root";
  //    //TFile* file = TFile::Open(samplePath_+fileName);
  //    //TTree* oldTree = (TTree*)file->Get("FullHadTreeWriter/tree");
  //    TChain* oldTree = new TChain("FullHadTreeWriter/tree");
  //    oldTree->Add(samplePath_+fileName);
  //    tmpFile->cd();
  //    bool isData = fileName.Contains("data") ? true : false;
  //    TTree* tree = modifiedTree_(oldTree, -10, 10, isData);
  //
  //    RooDataSet* data = new
  //    RooDataSet("data","data",varSet,RooFit::Import(*tree));//,RooFit::WeightVar("prob"));
  //
  //    name = "";
  //
  //    RooAbsPdf* model  = workspace[0]->pdf("model");
  //    RooAbsPdf* topBKG = workspace[0]->pdf("topBKG");
  //    RooAbsPdf* wBKG   = workspace[0]->pdf("wBKG");
  //
  //    RooAbsReal* fSig = workspace[0]->var("fSig");
  //
  //    RooAbsReal* nll = model->createNLL(*data);
  //    RooMinuit mini(*nll);
  //    mini.setStrategy(2);
  //    mini.migrad();
  //    mini.improve();
  //    mini.hesse();
  //
  //    TCanvas* canvas1 = new TCanvas("canvas1", "canvas1", 1, 1, 600, 600);
  //    canvas1->cd();
  //    RooPlot* frame1 = MTOP.frame(RooFit::Range(100,350));
  //    data  ->statOn (frame1, RooFit::Layout(.5, .9, .9));
  //    model ->paramOn(frame1, RooFit::Layout(.5, .9, .7));
  //    data  ->plotOn(frame1);
  //    model ->plotOn(frame1, RooFit::LineColor(8));
  //    topAdd->plotOn(frame1, RooFit::LineColor(kRed) , RooFit::Normalization(
  //    fSig->getVal()));
  //    topBKG->plotOn(frame1, RooFit::LineColor(kBlue),
  //    RooFit::Normalization(1.-fSig->getVal()));
  //    //workspace[0]->pdf("sig_0")   ->plotOn(frame1,
  //    RooFit::LineColor(kRed+1),
  //    RooFit::Normalization(1./*fSig.getVal()*workspace[0]->var("fCP")->getVal()*/));
  //    //workspace[0]->pdf("sig_1")   ->plotOn(frame1,
  //    RooFit::LineColor(kRed+2),
  //    RooFit::Normalization(1./*fSig.getVal()*workspace[0]->var("fWP")->getVal()*/));
  //    //workspace[0]->pdf("sig_2")   ->plotOn(frame1,
  //    RooFit::LineColor(kRed+3),
  //    RooFit::Normalization(1./*fSig.getVal()*workspace[0]->var("fUN")->getVal()*/));
  //    frame1->Draw();
  //
  //    //TCanvas* canvas2 = new TCanvas("canvas2", "canvas2", 641, 1, 600,
  //    600);
  //    //canvas2->cd();
  //    //
  //    //RooPlot* frame2 = fSig.frame();
  //    ////nll->plotOn(frame2,RooFit::ShiftToZero());
  //    //
  //    //RooAbsReal* pll_fSig = nll->createProfile(fSig);
  //    //pll_fSig->plotOn(frame2,RooFit::LineColor(kRed));
  //    //frame2->Draw();
  //    //
  //    //TCanvas* canvas3 = new TCanvas("canvas3", "canvas3", 641, 1, 600,
  //    600);
  //    //canvas3->cd();
  //    //
  //    //RooPlot* frame3 =
  //    workspace[0]->var("mTop")->frame();//RooFit::Bins(100),RooFit::Range(172,175));
  //    ////nll->plotOn(frame3,RooFit::ShiftToZero());
  //    //
  //    //RooAbsReal* pll_mTop = nll->createProfile(*workspace[0]->var("mTop"));
  //    //pll_mTop->plotOn(frame3,RooFit::LineColor(kRed));
  //    //frame3->Draw();
  //    //
  //    //TCanvas* canvas4 = new TCanvas("canvas4", "canvas4", 641, 1, 600,
  //    600);
  //    //canvas4->cd();
  //    //
  //    //RooPlot* frame4 =
  //    workspace[0]->var("JES")->frame();//RooFit::Bins(100),RooFit::Range(172,175));
  //    ////nll->plotOn(frame4,RooFit::ShiftToZero());
  //    //
  //    //RooAbsReal* pll_JES = nll->createProfile(*workspace[0]->var("JES"));
  //    //pll_JES->plotOn(frame4,RooFit::LineColor(kRed));
  //    //frame4->Draw();
  //
  //    TCanvas* canvas5 = new TCanvas("canvas5", "canvas5", 641, 1, 600, 600);
  //    canvas5->cd();
  //    RooPlot* frame5 = meanMW.frame(RooFit::Range(60,140));
  //    data ->statOn (frame5, RooFit::Layout(.5, .9, .9));
  //    model->paramOn(frame5, RooFit::Layout(.5, .9, .7));
  //    data ->plotOn(frame5);
  //    model->plotOn(frame5, RooFit::LineColor(8));
  //    wAdd ->plotOn(frame5, RooFit::LineColor(kRed) , RooFit::Normalization(
  //    fSig->getVal()));
  //    wBKG ->plotOn(frame5, RooFit::LineColor(kBlue),
  //    RooFit::Normalization(1.-fSig->getVal()));
  //    //workspace[0]->pdf("sig_3")   ->plotOn(frame5,
  //    RooFit::LineColor(kRed+1),
  //    RooFit::Normalization(1./*fSig.getVal()*workspace[0]->var("fCP")->getVal()*/));
  //    //workspace[0]->pdf("sig_4")   ->plotOn(frame5,
  //    RooFit::LineColor(kRed+2),
  //    RooFit::Normalization(1./*fSig.getVal()*workspace[0]->var("fWP")->getVal()*/));
  //    //workspace[0]->pdf("sig_5")   ->plotOn(frame5,
  //    RooFit::LineColor(kRed+3),
  //    RooFit::Normalization(1./*fSig.getVal()*workspace[0]->var("fUN")->getVal()*/));
  //    frame5->Draw();
  //
  //    //tmpFile->Close();
  //  }
}

// because RooFit only likes plain trees with standard data types (int, float,
// double, ...)
// the original tree has to be adapted for the new content
TTree *TemplateDerivation::modifiedTree_(TChain *tree, int minComboType,
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

TTree *TemplateDerivation::modifiedTree_(TChain *tree, int comboType) {
  return modifiedTree_(tree, comboType, comboType);
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
  TString varName = "alpha_";
  varName += h;
  varName += "_";
  varName += alpha.size();
  alpha.push_back(new RooFormulaVar(varName, varName, formula_alpha, argSet));
}
