#include "calibration.h"

//#include "TArrow.h"
#include "TAxis.h"
#include "TCanvas.h"
#include "TF1.h"
//#include "TF2.h"
#include "TFile.h"
//#include "TGraphErrors.h"
//#include "TGraph2DErrors.h"
//#include "TH1F.h"
#include "TLatex.h"
#include "TLegend.h"
//#include "TPaveText.h"
//#include "TMath.h"
//#include "TRandom3.h"
#include "TROOT.h"
#include "TString.h"
#include "TStyle.h"
//#include "TSystem.h"
#include "TLorentzVector.h"
#include "TClonesArray.h"
#include "TEventList.h"

#include "RooAddition.h"
#include "RooAddPdf.h"
#include "RooBifurGauss.h"
#include "RooCategory.h"
#include "RooDataHist.h"
#include "RooDataSet.h"
#include "RooFitResult.h"
#include "RooGamma.h"
#include "RooCBShape.h"
#include "RooGaussian.h"
#include "RooLandau.h"
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
#include "sstream"

typedef ProgramOptionsReader po;

TopMassCalibration::TopMassCalibration() :
    bTagEff_(0),
    cTagEff_(0),
    lTagEff_(0),
    //samplePath_("/scratch/hh/dust/naf/cms/user/eschliec/TopMass/19/"),
    //samplePath_("/scratch/hh/lustre/cms/user/eschliec/TopMass/19/"),
    selection_  (po::GetOption<std::string>("analysisConfig.selection")),
    samplePath_ (po::GetOption<std::string>("analysisConfig.samplePath")),
    fChannel_   (po::GetOption<std::string>("channel")),
    doCalibration_(true),
    fitBackground_(false),
    doMeasurement_(false)
{
  if      (!strcmp(fChannel_, "alljets" )) channelID_ = kAllJets;
  else if (!strcmp(fChannel_, "muon"    )) channelID_ = kMuonJets;
  else if (!strcmp(fChannel_, "electron")) channelID_ = kElectronJets;
  else if (!strcmp(fChannel_, "lepton"  )) channelID_ = kLeptonJets;
  else UnknownChannelAbort();

  rooFitTopMass_();
}

TopMassCalibration::~TopMassCalibration()
{
  delete bTagEff_;
  delete cTagEff_;
  delete lTagEff_;
}

void
TopMassCalibration::rooFitTopMass_()
{
  TFile *tmpFile = TFile::Open("tmpFileRooFitTopMass.root","RECREATE");
  RooWorkspace *workspace[1];
  workspace[0] = new RooWorkspace("workspaceMtop", "workspaceMtop");

  const int nTemplTypes = 2; // number of different distributions, e.g., mTop & mW
  const int nComboTypes = 3; // number of different permutation types, e.g., correct, wrong, unmatched
  const int nPDFs = nTemplTypes*nComboTypes;

  //bool binnedTemplates = false;

  TFile *bTagFile = TFile::Open(samplePath_+TString("bTagFile.root"));
  gROOT->cd();
  bTagEff_ = (TH2F*)bTagFile->Get("histb")->Clone();
  cTagEff_ = (TH2F*)bTagFile->Get("histc")->Clone();
  lTagEff_ = (TH2F*)bTagFile->Get("histl")->Clone();
  bTagFile->Close();
  tmpFile->cd();

  /// definitions for fit
  const unsigned short nJES = 5; // number of JES to be used for calibration
  const unsigned short nMasses = 9; // number of top masses to be used for calibration
  const unsigned short nTemplates = nJES*nMasses;

  const double  jes_templ[nJES] = {0.96, 0.98, 1.00, 1.02, 1.04};
  const double mTop_templ[nMasses] = {161.5, 163.5, 166.5, 169.5, 172.5, 175.5, 178.5, 181.5, 184.5};

  const unsigned iTemplateJES [] = {0,1,2,3,4,
                                    0,1,2,3,4,
                                    0,1,2,3,4,
                                    0,1,2,3,4,
                                    0,1,2,3,4,
                                    0,1,2,3,4,
                                    0,1,2,3,4,
                                    0,1,2,3,4,
                                    0,1,2,3,4};

  const unsigned iTemplateMass[] = {0,0,0,0,0,
                                    1,1,1,1,1,
                                    2,2,2,2,2,
                                    3,3,3,3,3,
                                    4,4,4,4,4,
                                    5,5,5,5,5,
                                    6,6,6,6,6,
                                    7,7,7,7,7,
                                    8,8,8,8,8};

  TString templ[nTemplates] = {"jes096mass1615","jes098mass1615","jes100mass1615","jes102mass1615","jes104mass1615",
                               "jes096mass1635","jes098mass1635","jes100mass1635","jes102mass1635","jes104mass1635",
                               "jes096mass1665","jes098mass1665","jes100mass1665","jes102mass1665","jes104mass1665",
                               "jes096mass1695","jes098mass1695","jes100mass1695","jes102mass1695","jes104mass1695",
                               "jes096mass1725","jes098mass1725","jes100mass1725","jes102mass1725","jes104mass1725",
                               "jes096mass1755","jes098mass1755","jes100mass1755","jes102mass1755","jes104mass1755",
                               "jes096mass1785","jes098mass1785","jes100mass1785","jes102mass1785","jes104mass1785",
                               "jes096mass1815","jes098mass1815","jes100mass1815","jes102mass1815","jes104mass1815",
                               "jes096mass1845","jes098mass1845","jes100mass1845","jes102mass1845","jes104mass1845"};

  TString templateSettings[nTemplates] = {"MC for JES = 0.96, mTop = 161.5 GeV","MC for JES = 0.98, mTop = 161.5 GeV","MC for JES = 1.00, mTop = 161.5 GeV","MC for JES = 1.02, mTop = 161.5 GeV","MC for JES = 1.04, mTop = 161.5 GeV",
                                          "MC for JES = 0.96, mTop = 163.5 GeV","MC for JES = 0.98, mTop = 163.5 GeV","MC for JES = 1.00, mTop = 163.5 GeV","MC for JES = 1.02, mTop = 163.5 GeV","MC for JES = 1.04, mTop = 163.5 GeV",
                                          "MC for JES = 0.96, mTop = 166.5 GeV","MC for JES = 0.98, mTop = 166.5 GeV","MC for JES = 1.00, mTop = 166.5 GeV","MC for JES = 1.02, mTop = 166.5 GeV","MC for JES = 1.04, mTop = 166.5 GeV",
                                          "MC for JES = 0.96, mTop = 169.5 GeV","MC for JES = 0.98, mTop = 169.5 GeV","MC for JES = 1.00, mTop = 169.5 GeV","MC for JES = 1.02, mTop = 169.5 GeV","MC for JES = 1.04, mTop = 169.5 GeV",
                                          "MC for JES = 0.96, mTop = 172.5 GeV","MC for JES = 0.98, mTop = 172.5 GeV","MC for JES = 1.00, mTop = 172.5 GeV","MC for JES = 1.02, mTop = 172.5 GeV","MC for JES = 1.04, mTop = 172.5 GeV",
                                          "MC for JES = 0.96, mTop = 175.5 GeV","MC for JES = 0.98, mTop = 175.5 GeV","MC for JES = 1.00, mTop = 175.5 GeV","MC for JES = 1.02, mTop = 175.5 GeV","MC for JES = 1.04, mTop = 175.5 GeV",
                                          "MC for JES = 0.96, mTop = 178.5 GeV","MC for JES = 0.98, mTop = 178.5 GeV","MC for JES = 1.00, mTop = 178.5 GeV","MC for JES = 1.02, mTop = 178.5 GeV","MC for JES = 1.04, mTop = 178.5 GeV",
                                          "MC for JES = 0.96, mTop = 181.5 GeV","MC for JES = 0.98, mTop = 181.5 GeV","MC for JES = 1.00, mTop = 181.5 GeV","MC for JES = 1.02, mTop = 181.5 GeV","MC for JES = 1.04, mTop = 181.5 GeV",
                                          "MC for JES = 0.96, mTop = 184.5 GeV","MC for JES = 0.98, mTop = 184.5 GeV","MC for JES = 1.00, mTop = 184.5 GeV","MC for JES = 1.02, mTop = 184.5 GeV","MC for JES = 1.04, mTop = 184.5 GeV"};

  RooCategory* cat_templ = new RooCategory("cat_templ", "cat_templ");
  for(unsigned t=0; t<nTemplates; ++t)
    cat_templ->defineType(templ[t]);

  workspace[0]->import(*cat_templ);

  /// create datasets for later use
  RooDataSet *dataset[nTemplates];
  RooDataSet *reducedDataset[nTemplates];
  //RooDataHist* hist[1][nTemplates];

  RooRealVar comboTypeVar   = RooRealVar("comboType"     ,"comboType"  ,  -10.,   6.9,"");
  //RooRealVar prob           = RooRealVar("prob"          ,"P(#chi^{2})",  0.,   1.,"");
  RooRealVar MTOP           = RooRealVar("topMass"       ,"m_{t}^{fit}",100., 550.,"GeV");
  RooRealVar meanMW         = RooRealVar("meanWMass"     ,"m_{W}^{rec}", 50., 300.,"GeV");
  RooRealVar combinedWeight = RooRealVar("combinedWeight","weight"     ,  0., 100.,"");

  if(channelID_ == kAllJets){
    comboTypeVar.setRange("R1",0.9,1.1);
    comboTypeVar.setRange("R4",-9.9,-0.1);
    comboTypeVar.setRange("R5",1.9,6.1);
  }
  MTOP  .setRange("mTopFitRange",100.,550.);
  meanMW.setRange("mWFitRange", 50.,300.);

  RooArgSet varSet = RooArgSet(comboTypeVar,/*prob,*/MTOP,meanMW,combinedWeight,"varSet");

  RooAddPdf *topAdd = 0;
  RooAddPdf   *wAdd = 0;

  TString name = "";

  if(doCalibration_){

    for(unsigned int iMass = 0; iMass < nMasses; ++iMass){
      for(unsigned int iJES = 0; iJES < nJES; ++iJES){
        TString fileName;
        if(channelID_ == kAllJets){
          //fileName += "Z2_F11_ABS_JES";
          //fileName += "Z2_S12_ABS_JES";
          fileName += "TopMassTreeWriter_02_MC01/Z2_S12_ABS_JES";
          if     (iJES  == 0) fileName += "_096_";
          else if(iJES  == 1) fileName += "_098_";
          else if(iJES  == 2) fileName += "_100_";
          else if(iJES  == 3) fileName += "_102_";
          else if(iJES  == 4) fileName += "_104_";
          if     (iMass == 0) fileName += "161_5";
          else if(iMass == 1) fileName += "163_5";
          else if(iMass == 2) fileName += "166_5";
          else if(iMass == 3) fileName += "169_5";
          else if(iMass == 4) fileName += "172_5";
          else if(iMass == 5) fileName += "175_5";
          else if(iMass == 6) fileName += "178_5";
          else if(iMass == 7) fileName += "181_5";
          else if(iMass == 8) fileName += "184_5";
          fileName += "_sig/*.root";
          //fileName += "_sig.root";
        }
        std::cout << "Creating RooDataSet for: " << fileName;

        //TFile *file = TFile::Open(samplePath_+TString(fileName));
        //TTree *tree = (TTree*)file->Get("analyzeKinFit/eventTree");
        TChain *chain = new TChain("analyzeKinFit/eventTree");
        if(iMass == 4) {
          TString whichJES;
          if     (iJES  == 0) whichJES += "_096_";
          else if(iJES  == 1) whichJES += "_098_";
          else if(iJES  == 2) whichJES += "_100_";
          else if(iJES  == 3) whichJES += "_102_";
          else if(iJES  == 4) whichJES += "_104_";
          std::cout << " (nFiles: " << chain->Add(samplePath_+TString("TopMassTreeWriter_02_MC01/Z2_S12_Had1_ABS_JES")+whichJES+TString("172_5_sig/*.root"))
                                     + chain->Add(samplePath_+TString("TopMassTreeWriter_02_MC01/Z2_S12_Had2_ABS_JES")+whichJES+TString("172_5_sig/*.root"))
                                     + chain->Add(samplePath_+TString("TopMassTreeWriter_02_MC01/Z2_S12_Semi_ABS_JES")+whichJES+TString("172_5_sig/*.root"))
                                     + chain->Add(samplePath_+TString("TopMassTreeWriter_02_MC01/Z2_S12_Lept_ABS_JES")+whichJES+TString("172_5_sig/*.root"))
                                    << ")" << std::endl;
        }
        else {
          std::cout << " (nFiles: " << chain->Add(samplePath_+TString(fileName)) << ")" << std::endl;
        }
        //TTree *tree = (TTree*)file->Get("FullHadTreeWriter/tree");
        tmpFile->cd();
        TTree* tree = modifiedTree_(chain); //, minComboType, maxComboType);
        //file->Close();
        int iTempl = iMass*nJES+iJES;
        name = "dataset_"; name += iTempl;

        // needed to avoid a crash
        double co, /*pr,*/ mt, mw, cw;
        tree->SetBranchAddress("comboType", &co);
        //tree->SetBranchAddress("prob", &pr);
        tree->SetBranchAddress("topMass", &mt);
        tree->SetBranchAddress("meanWMass", &mw);
        tree->SetBranchAddress("combinedWeight", &cw);

        dataset[iTempl] = new RooDataSet(name,name,varSet,RooFit::Import(*tree),RooFit::WeightVar("combinedWeight"));
        //if(binnedTemplates){
        //  for(unsigned h=0; h<1; ++h) {
        //    name = "hist_"; name += iTempl;
        //    hist[h][iTempl] = new RooDataHist(name, name, mTop, *dataset[iTempl]);
        //  }
        //}
      }
    }

    RooSimWSTool* simWST[nPDFs];
    RooSimultaneous* sim[nPDFs];
    RooNLLVar* nll[nPDFs][nTemplates];
    RooAddition* nllSim[nPDFs];

    RooRealVar* JES  = new RooRealVar("JES" , "JES" , 1.0, 0.8, 1.2);
    RooRealVar* mTop = new RooRealVar("mTop", "mTop", 172.5, 100., 550., "GeV");
    JES ->setConstant(kTRUE);
    mTop->setConstant(kTRUE);

    for(int templType = 0; templType < nTemplTypes; ++templType){
      for(int comboType = 0; comboType < nComboTypes; ++comboType){

        int h = nComboTypes*templType + comboType;

        RooRealVar MassOffset        = RooRealVar("MassOffset"       ,"MassOffset"       , 0.); MassOffset       .setConstant(kTRUE);
        RooRealVar MassSlopeMass     = RooRealVar("MassSlopeMass"    ,"MassSlopeMass"    , 0.); MassSlopeMass    .setConstant(kTRUE);
        RooRealVar MassSlopeJES      = RooRealVar("MassSlopeJES"     ,"MassSlopeJES"     , 0.); MassSlopeJES     .setConstant(kTRUE);
        RooRealVar MassSlopeMassJES  = RooRealVar("MassSlopeMassJES" ,"MassSlopeMassJES" , 0.); MassSlopeMassJES .setConstant(kTRUE);
        RooRealVar JESOffset         = RooRealVar("JESOffset"        ,"JESOffset"        , 0.); JESOffset        .setConstant(kTRUE);
        RooRealVar JESSlopeMass      = RooRealVar("JESSlopeMass"     ,"JESSlopeMass"     , 0.); JESSlopeMass     .setConstant(kTRUE);
        RooRealVar JESSlopeJES       = RooRealVar("JESSlopeJES"      ,"JESSlopeJES"      , 0.); JESSlopeJES      .setConstant(kTRUE);
        RooRealVar JESSlopeMassJES   = RooRealVar("JESSlopeMassJES"  ,"JESSlopeMassJES"  , 0.); JESSlopeMassJES  .setConstant(kTRUE);
        RooRealVar MassOffset2       = RooRealVar("MassOffset2"      ,"MassOffset2"      , 0.); MassOffset2      .setConstant(kTRUE);
        RooRealVar MassSlopeMass2    = RooRealVar("MassSlopeMass2"   ,"MassSlopeMass2"   , 0.); MassSlopeMass2   .setConstant(kTRUE);
        RooRealVar MassSlopeJES2     = RooRealVar("MassSlopeJES2"    ,"MassSlopeJES2"    , 0.); MassSlopeJES2    .setConstant(kTRUE);
        RooRealVar MassSlopeMassJES2 = RooRealVar("MassSlopeMassJES2","MassSlopeMassJES2", 0.); MassSlopeMassJES2.setConstant(kTRUE);
        RooRealVar JESOffset2        = RooRealVar("JESOffset2"       ,"JESOffset2"       , 0.); JESOffset2       .setConstant(kTRUE);
        RooRealVar JESSlopeMass2     = RooRealVar("JESSlopeMass2"    ,"JESSlopeMass2"    , 0.); JESSlopeMass2    .setConstant(kTRUE);
        RooRealVar JESSlopeJES2      = RooRealVar("JESSlopeJES2"     ,"JESSlopeJES2"     , 0.); JESSlopeJES2     .setConstant(kTRUE);
        RooRealVar JESSlopeMassJES2  = RooRealVar("JESSlopeMassJES2" ,"JESSlopeMassJES2" , 0.); JESSlopeMassJES2 .setConstant(kTRUE);
        // 0 = JES, 1 = mTop, 2 = JESOffset, 3 = JESSlopeMass, 4 = JESSlopeJES, 5 = JESSlopeMassJES
        TString formula_mTop = "@1+@2+@3*(@1-172.5)+@4*(@0-1.0)+@5*(@1-172.5)*(@0-1.0)";
        TString formula_JES  = "@0+@2+@3*(@1-172.5)+@4*(@0-1.0)+@5*(@1-172.5)*(@0-1.0)";
        name = "mTop_intermediate_"; name += h;
        RooFormulaVar mTop_intermediate = RooFormulaVar(name, name, formula_mTop, RooArgSet(*JES, *mTop, MassOffset, MassSlopeMass, MassSlopeJES, MassSlopeMassJES));
        name = "mJES_intermediate_"; name += h;
        RooFormulaVar  JES_intermediate = RooFormulaVar(name, name, formula_JES , RooArgSet(*JES, *mTop,  JESOffset,  JESSlopeMass,  JESSlopeJES,  JESSlopeMassJES));
        name = "mTop_corrected_"; name += h;
        RooFormulaVar mTop_corrected = RooFormulaVar(name, name, formula_mTop, RooArgSet(JES_intermediate, mTop_intermediate, MassOffset2, MassSlopeMass2, MassSlopeJES2, MassSlopeMassJES2));
        name = "mJES_corrected_"; name += h;
        RooFormulaVar  JES_corrected = RooFormulaVar(name, name, formula_JES , RooArgSet(JES_intermediate, mTop_intermediate,  JESOffset2,  JESSlopeMass2,  JESSlopeJES2,  JESSlopeMassJES2));

        RooRealVar *topMass   = new RooRealVar("topMass"  , "m_{t}^{fit}", 100.,  550., "GeV");
        RooRealVar *meanWMass = new RooRealVar("meanWMass", "m_{W}^{rec}",  50.,  300., "GeV");
        RooRealVar *var = 0;
        TString comboRangeName = "";
        std::vector<RooRealVar*> par;
        std::vector<RooFormulaVar*> alpha;
        std::vector<double> iniPar;

        if(channelID_ == kAllJets){
          if     (comboType == 0){ comboRangeName = "R1"; }
          else if(comboType == 1){ comboRangeName = "R4"; }
          else if(comboType == 2){ comboRangeName = "R5"; }

          if(templType == 0){
            if (comboType == 0){ // mTop, correct
              //double a[] = {172.328, 0.982351 , 85.7756 , 0.871109,
              //              7.17985, 0.0687401,  6.68148, 0.281139};
              //iniPar = std::vector<double>(a,a+sizeof(a)/sizeof(double));
              iniPar = {172.478, 0.980786 , 79.2462 ,  0.724549,
                        7.88763, 0.0719042,  6.25116, -0.0799414};
            }
            else if (comboType == 1){ // mTop, wrong
              double a[] = {182.221,  0.180998 , 35.0042 ,  0.701448,
                            172.565,  1.0841   , 99.7299 ,  0.226163,
                            24.3286, -0.0251085, -6.24156,  0.392888,
                            8.87901, -0.072422 ,  7.45875, -1.6491};
              iniPar = std::vector<double>(a,a+sizeof(a)/sizeof(double));
              //iniPar = {171.534 , 0.36392  , 50.5469 ,  0.885869,
              //          174.476 , 0.967762 , 84.2048 , -0.469726,
              //           20.2423, 0.0825402, 15.4528 ,  1.20711,
              //           12.4386, 0.0652505,  7.62482,  0.000278128};
            }
            else if(comboType == 2){ // mTop, unmatched
              //double a[] = {178.206, 0.26746    , 35.9815  , 1.50419 ,
              //              174.467, 0.991727   , 88.1759  , 0.991246,
              //              21.9042, 0.0375821  ,  3.86448 , 1.34186 ,
              //               8.9953, 0.000352339,  0.885277, 0.126578};
              //iniPar = std::vector<double>(a,a+sizeof(a)/sizeof(double));
              iniPar = {178.883, 0.306507 , 40.7573 ,  0.445759,
                        174.681, 1.06829  , 83.8601 ,  0.280563,
                         22.065, 0.0542991,  4.93823,  0.3491,
                          9.74 , 0.0548014,  2.31277, -0.181041};
            }
          }
          else if(templType == 1){
            if (comboType == 0){ // mW, correct
              //double a[] = {84.4714 , -0.00288904, 57.2879 ,  0.221738 ,
              //               4.25403,  0.00163248, 11.8082 ,  0.0886389,
              //               5.43406, -0.00640861, -6.93172, -0.0742648};
              //iniPar = std::vector<double>(a,a+sizeof(a)/sizeof(double));
              iniPar = {84.5861, -0.00713988, 82.7429, -0.196859,
                        5.20718, -0.000562235, 21.8792, -0.13701,
                        6.76429, 0.000500302, -17.8822, 0.0366887};
            }
            else if (comboType == 1){ // mW, wrong
              //double a[] = {84.059 ,  0.00363846, 27.6651, -0.0122909,
              //              4.28601, -0.0106779 ,  5.2756, -0.18915  ,
              //              6.2658 , -0.0133067 , -6.2384,  0.0742743};
              //iniPar = std::vector<double>(a,a+sizeof(a)/sizeof(double));
              iniPar = {85.2743, -0.0165617, 26.6677, -0.332143,
                        5.58215, -0.0131054, 6.19793, -0.200832,
                        7.56877, 0.00863095, -9.94996, 0.100932};
            }
            else if(comboType == 2){ // mW, unmatched
              //double a[] = {84.6333,  0.0205799 , 16.0275, -0.280052,
              //              4.76909,  0.018934  , 2.18012, -0.12839 ,
              //              7.08441, -0.00812571, -6.2732,  0.22637 };
              //iniPar = std::vector<double>(a,a+sizeof(a)/sizeof(double));
              iniPar = {86.6185, -0.0184471, 32.1242, 0.0815632,
                        6.48169, -0.0108763, 7.37265, 0.0139849,
                        7.91363, 0.00815331, -17.3423, -0.00100213};
            }
          }

          for(unsigned p=0; p<iniPar.size(); ++p){
            //std::cout << iniPar[p] << ", ";
            name = "par_"; name += h; name += "_"; name += p;
            RooRealVar *myVar = new RooRealVar(name, name, iniPar[p]); myVar->setConstant(kFALSE);
            par.push_back(myVar);
          }
          std::cout << std::endl;
          name = "sig_"; name += h;
          TString varName;
          if(templType == 0){
            var = topMass;
            if(comboType == 0) {
              fillAlpha(alpha, h, RooArgSet(*par[0], *par[1], *par[2], *par[3], mTop_corrected, JES_corrected));
              fillAlpha(alpha, h, RooArgSet(*par[4], *par[5], *par[6], *par[7], mTop_corrected, JES_corrected));
              varName = "width_"; varName += h;
              RooRealVar *width = new RooRealVar(varName, varName, 2.0);
              width->setConstant(kTRUE);
              RooVoigtian voigt = RooVoigtian(name, name, *var, *alpha[0], *width, *alpha[1]);
              workspace[0]->import(voigt);
            }
            else if(comboType == 1 || comboType == 2) {
              varName = "ratio_"; varName += h;
              RooRealVar *ratio = new RooRealVar(varName, varName, 0.0, 1.0);
              if     (comboType == 1) ratio->setVal(0.662784); //0.820546);
              else if(comboType == 2) ratio->setVal(0.774588); //0.758690);
              ratio->setConstant(kTRUE);
              //ratio->setConstant(kFALSE);
              fillAlpha(alpha, h, RooArgSet(*par[0], *par[1], *par[ 2], *par[ 3], mTop_corrected, JES_corrected));
              fillAlpha(alpha, h, RooArgSet(*par[8], *par[9], *par[10], *par[11], mTop_corrected, JES_corrected));
              varName = "landau_"; varName += h;
              RooLandau *landau = new RooLandau(varName, varName, *var, *alpha[0], *alpha[1]);
              fillAlpha(alpha, h, RooArgSet(*par[ 4], *par[ 5], *par[ 6], *par[ 7], mTop_corrected, JES_corrected));
              fillAlpha(alpha, h, RooArgSet(*par[12], *par[13], *par[14], *par[15], mTop_corrected, JES_corrected));
              varName = "width_"; varName += h;
              RooRealVar *width = new RooRealVar(varName, varName, 2.0);
              width->setConstant(kTRUE);
              varName = "voigt_"; varName += h;
              RooVoigtian *voigt = new RooVoigtian(varName, varName, *var, *alpha[2], *width, *alpha[3]);
              RooAddPdf *add = new RooAddPdf(name, name, *landau, *voigt, *ratio);
              workspace[0]->import(*add);
            }
          }
          else if(templType == 1){
            var = meanWMass;
            fillAlpha(alpha, h, RooArgSet(*par[0], *par[1], *par[ 2], *par[ 3], mTop_corrected, JES_corrected));
            fillAlpha(alpha, h, RooArgSet(*par[4], *par[5], *par[ 6], *par[ 7], mTop_corrected, JES_corrected));
            fillAlpha(alpha, h, RooArgSet(*par[8], *par[9], *par[10], *par[11], mTop_corrected, JES_corrected));
            RooBifurGauss *asymGaus = new RooBifurGauss(name, name, *var, *alpha[0], *alpha[1], *alpha[2]);
            workspace[0]->import(*asymGaus);
          }
        }
        simWST[h] = new RooSimWSTool(*workspace[0]);
        name = "sim_"; name += h;
        TString name_pdf = "sig_"; name_pdf += h;
        sim[h] = simWST[h]->build(name, name_pdf,
            RooFit::SplitParam("JES" ,"cat_templ"),
            RooFit::SplitParam("mTop","cat_templ"));
        for(unsigned t=0; t<nTemplates; ++t) {
          workspace[0]->var("JES_" +templ[t])->setVal(jes_templ [iTemplateJES [t]]);
          workspace[0]->var("mTop_"+templ[t])->setVal(mTop_templ[iTemplateMass[t]]);
          name = "nll_"; name += h; name += "_"; name += t;
          name_pdf = "sig_"; name_pdf += h; name_pdf += "_" + templ[t];
          //if(binnedTemplates)
          //  nll[h][t] = new RooNLLVar(name, name, *workspace[0]->pdf(name_pdf), *hist[h][t], RooFit::NumCPU(1), RooFit::Range(100., 550.));
          //else
          reducedDataset[t] = (RooDataSet*)dataset[t]->reduce(RooFit::CutRange(comboRangeName));
          std::cout << "Entries in " << name << ": " << reducedDataset[t]->numEntries() << std::endl;
          nll[h][t] = new RooNLLVar(name, name, *workspace[0]->pdf(name_pdf), *reducedDataset[t], RooFit::NumCPU(1), RooFit::Range("mTopFitRange"));
        }
        name = "nllSim_"; name += h;
        RooArgSet argSet;
        for(unsigned t=0; t<nTemplates; ++t)
          argSet.add(*nll[h][t]);
        nllSim[h] = new RooAddition(name, name, argSet);
        RooFitResult* result;
        gStyle->SetOptTitle(1);
        // perform the fit
        RooMinuit minuit(*nllSim[h]);
        minuit.setStrategy(2);
        minuit.migrad();
        minuit.improve();
        minuit.hesse();
        //minuit.minos();
        result = minuit.save();
        workspace[0]->import(*result);

        /// control plots
        TString outDir = "plot/calibration";

        RooPlot* frame = 0;
        TCanvas* canvas = new TCanvas("canvas", "canvas", 10, 10, 600, 600);
        TString figName = "";
        for(unsigned t=0; t<nTemplates; ++t) {
          //frame = var_mass[h]->frame();
          if     (templType == 0 && comboType == 0) frame = var->frame(RooFit::Range(100., 250.));
          else if(templType == 0)                   frame = var->frame(RooFit::Range(100., 550.));
          else if(templType == 1)                   frame = var->frame(RooFit::Range( 50., 120.));
          //if(! binnedTemplates) {
          reducedDataset[t]->statOn(frame, RooFit::Layout(.6, .9, .9));
          reducedDataset[t]->plotOn(frame);
          //}
          //else {
          //  hist[h][t]->statOn(frame, RooFit::Layout(.6, .9, .9));
          //  hist[h][t]->plotOn(frame);
          //}
          name_pdf = "sig_"; name_pdf += h; name_pdf += "_"+templ[t];
          std::cout << "PDF Name: " << name_pdf << std::endl;
          workspace[0]->pdf(name_pdf)->plotOn(frame, RooFit::FillColor(kGray), RooFit::VisualizeError(*result));
          workspace[0]->pdf(name_pdf)->plotOn(frame, RooFit::LineColor(kRed+1));
          frame->SetMinimum(.0);
          frame->SetTitle(templateSettings[t]);
          frame->Draw();
          figName = "template_"; figName += templType; figName += "_"; figName += comboType; figName += "_"; figName += t;
          canvas->Print(outDir + "/" +  figName + ".eps");
          canvas->Print(outDir + "/catalog.ps(");
        }
        //const int redPalette[5] = {kRed-2, kRed-6, kRed+1, kRed-8, kRed+3};
        //const int redPalette[9] = {kRed-9, kRed-8, kRed-7, kRed-6, kRed-5, kRed-4, kRed-3, kRed-2, kRed-1};
        // show functions for different JES values
        //frame = var_mass[h]->frame();
        if     (templType == 0 && comboType == 0) frame = var->frame(RooFit::Range(100., 250.));
        else if(templType == 0)                   frame = var->frame(RooFit::Range(100., 550.));
        else if(templType == 1)                   frame = var->frame(RooFit::Range( 50., 120.));
        for(unsigned t=20; t<25; ++t) {
          name_pdf = "sig_"; name_pdf += h; name_pdf += "_"+templ[t];
          //std::cout << name_pdf << std::endl;
          //workspace[0]->pdf(name_pdf)->plotOn(frame, RooFit::FillColor(kGray), RooFit::VisualizeError(*result));
          workspace[0]->pdf(name_pdf)->plotOn(frame, RooFit::LineColor(kRed));
        }
        frame->SetMinimum(.0);
        frame->GetYaxis()->SetTitle("Probability");
        frame->SetTitle("JES = 0.96, 0.98, 1.00, 1.02, 1.04 (mTop = 172.5 GeV)");
        frame->Draw();
        figName = "funcFamily_jes_"; figName += templType; figName += "_"; figName += comboType;
        canvas->Print(outDir + "/" + figName + ".eps");
        canvas->Print(outDir + "/catalog.ps");
        // show functions for different mTop values
        //frame = var_mass[h]->frame();
        if     (templType == 0 && comboType == 0) frame = var->frame(RooFit::Range(100., 250.));
        else if(templType == 0)                   frame = var->frame(RooFit::Range(100., 550.));
        else if(templType == 1)                   frame = var->frame(RooFit::Range( 50., 120.));
        //name_pdf = "sig_"; name_pdf += h; name_pdf += "_"+templ[2];
        //workspace[h]->pdf(name_pdf)->plotOn(frame, RooFit::LineColor(kRed+1));
        for(unsigned t=1, i=0; t<nTemplates; t+=5, ++i) {
          name_pdf = "sig_"; name_pdf += h; name_pdf += "_"+templ[t];
          //std::cout << name_pdf << std::endl;
          //workspace[0]->pdf(name_pdf)->plotOn(frame, RooFit::FillColor(kGray), RooFit::VisualizeError(*result));
          workspace[0]->pdf(name_pdf)->plotOn(frame, RooFit::LineColor(kRed));
        }
        frame->SetMinimum(.0);
        frame->GetYaxis()->SetTitle("Probability");
        frame->SetTitle("mTop = 161.5, 163.5, 166.5, 169.5, 172.5, 175.5, 178.5, 181.5, 184.5 GeV (JES = 1.00)");
        frame->Draw();
        figName = "funcFamily_mass_"; figName += templType; figName += "_"; figName += comboType;
        canvas->Print(outDir + "/" + figName + ".eps");
        canvas->Print(outDir + "/catalog.ps");
        gStyle->SetOptTitle(0);
        // plot different alpha_i as a function of JES
        char yTitle[99];
        for(unsigned i=0; i<alpha.size(); ++i) {
          frame = workspace[0]->var("JES_"+templ[0])->frame(RooFit::Range(0.9, 1.1));
          name = "alpha_"; name += h; name += "_"; name += i;
          if(workspace[0]->function(name+"_"+templ[0])) {
            workspace[0]->function(name+"_"+templ[0])->plotOn(frame, RooFit::FillColor(kGray),
                RooFit::VisualizeError(*result));
            workspace[0]->function(name+"_"+templ[0])->plotOn(frame, RooFit::Name("sig"),
                RooFit::LineColor(kRed+1),
                RooFit::FillColor(kRed+1));
          }
          frame->SetMinimum(.0);
          if(h==0)
            sprintf(yTitle, "#alpha_{%i} (m_{top} = 172.5 GeV)", i);
          else
            sprintf(yTitle, "#alpha_{%i}", i);
          frame->GetYaxis()->SetTitle(yTitle);
          frame->GetXaxis()->SetTitle("JES");
          frame->Draw();
          double yLegHigh = 0.31;
          double yLegLow = yLegHigh-2*0.06;
          TLegend legend(0.72, yLegLow, 0.92, yLegHigh);
          legend.SetBorderSize(0);
          legend.SetFillStyle(0);
          legend.SetTextFont(42);
          //legend.AddEntry("sig", "P ("+var_mass[h]->getTitle()+")", "F");
          legend.AddEntry("sig", "P ("+var->getTitle()+")", "F");
          legend.Draw();
          figName = "alpha_"; figName += i; figName += "_jes_"; figName += templType; figName += "_"; figName += comboType;
          canvas->Print(outDir + "/" +  figName + ".eps");
          canvas->Print(outDir + "/catalog.ps");
          // plot different alpha_i as a function of mTop
          frame = workspace[0]->var("mTop_"+templ[2])->frame(RooFit::Range(150., 200.));
          name = "alpha_"; name += h; name += "_"; name += i;
          if(workspace[0]->function(name+"_"+templ[2])) {
            workspace[0]->function(name+"_"+templ[2])->plotOn(frame, RooFit::FillColor(kGray),
                RooFit::VisualizeError(*result));
            workspace[0]->function(name+"_"+templ[2])->plotOn(frame, RooFit::LineColor(kRed+1));
          }
          frame->SetMinimum(.0);
          sprintf(yTitle, "#alpha_{%i} (JES = 1.00)", i);
          frame->GetYaxis()->SetTitle(yTitle);
          frame->GetXaxis()->SetTitle("m_{top} (GeV)");
          frame->Draw();
          legend.Draw();
          figName = "alpha_"; figName += i; figName += "_mass_"; figName += templType; figName += "_"; figName += comboType;
          canvas->Print(outDir + "/" +  figName + ".eps");
          canvas->Print(outDir + "/catalog.ps)");
        }


        //////////////////////////////////////////////////////////////////////////////////////////////////////

        std::cout << "q[0] + q[1] * (x[0]-172.5) + (q[2] + q[3] * (x[0]-172.5)) * (x[1]-1.)" << std::endl;
        std::cout <<  "new: ";
        std::string vars_helper = "";
        std::stringstream vars(vars_helper, std::stringstream::in|std::stringstream::out);
        for(unsigned int i = 0; i < par.size(); ++i){
          TString varName = "par_"; varName += h; varName += "_"; varName += i;
          RooRealVar* curPar = workspace[0]->var(varName);
          curPar->setConstant(kTRUE);
          vars << curPar->getVal();
          if(i != par.size()-1) vars << ", ";
        }
        std::cout << vars.str();
        std::cout << std::endl;

        std::ofstream myfile;
        myfile.open(outDir + "/variables.txt", std::ios::out | std::ios::app);
        myfile << "Template Type: " << templType << ", Combo Type: " << comboType << "\n";
        myfile << vars.str() << "\n";
        myfile.close();
        std::cout << "old: ";
        for(unsigned int i = 0; i < iniPar.size(); ++i){
          std::cout << iniPar[i];
          if(i != iniPar.size()-1) std::cout << ", ";
        }
        std::cout << std::endl;
        //}
      }
    }

    RooRealVar *fCP = new RooRealVar("fCP", "fCP", 2.79599875211715698e-01);
    fCP->setConstant(kTRUE);
    RooRealVar *fWP = new RooRealVar("fWP", "fWP", 4.28881868720054626e-03+1.13365799188613892e-03+2.13942706584930420e-01);
    fWP->setConstant(kTRUE);
    RooRealVar *fUN = new RooRealVar("fUN", "fUN", 5.01034975051879883e-01);
    fUN->setConstant(kTRUE);
    RooArgSet permutationFractions = RooArgSet(*fCP,*fWP,*fUN,"permutationFractions");

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

    RooDataSet *BKG = 0;
    if(fitBackground_){
      std::cout << "Creating BKG dataset" << std::endl;

      TString fileName = "QCDEstimationMix_2011_NEW_skimmed2.root";
      //TFile* file = TFile::Open(samplePath_+fileName);
      //TTree* oldTree = (TTree*)file->Get("tree");
      TChain* oldTree = new TChain("tree");
      oldTree->Add(samplePath_+fileName);
      tmpFile->cd();
      TTree* tree = modifiedTree_(oldTree);
      //TTree* tree = modifiedTree_(oldTree, -1, 10, true);

      BKG = new RooDataSet("BKG","BKG",varSet,RooFit::Import(*tree));//,RooFit::WeightVar("prob"));

      std::cout << "Creating BKG PDFs" << std::endl;
    }

    RooRealVar *topBKGGammaNorm  = new RooRealVar("topBKGGammaNorm" , "topBKGGammaNorm" ,  0.848667);
    RooRealVar *topBKGGammaGamma = new RooRealVar("topBKGGammaGamma", "topBKGGammaGamma",  4.14752 );
    RooRealVar *topBKGGammaBeta  = new RooRealVar("topBKGGammaBeta" , "topBKGGammaBeta" , 33.7793  );
    RooRealVar *topBKGGammaMu    = new RooRealVar("topBKGGammaMu"   , "topBKGGammaMu"   , 89.6645  );
    RooGamma *topBKGGamma = new RooGamma("topBKGGamma","topBKGGamma",MTOP,*topBKGGammaGamma,*topBKGGammaBeta,*topBKGGammaMu);
    topBKGGamma->setNormRange("mTopFitRange");

    RooRealVar *topBKGLandauMean  = new RooRealVar("topBKGLandauMean" , "topBKGLandauMean" , 217.114 );
    RooRealVar *topBKGLandauSigma = new RooRealVar("topBKGLandauSigma", "topBKGLandauSigma",  27.6076);
    RooLandau *topBKGLandau = new RooLandau("topBKGLandau","topBKGLandau",MTOP,*topBKGLandauMean,*topBKGLandauSigma);
    topBKGLandau->setNormRange("mTopFitRange");

    RooAddPdf *topBKG = new RooAddPdf("topBKG", "topBKG", *topBKGGamma, *topBKGLandau, *topBKGGammaNorm);
    topBKG->setNormRange("mTopFitRange");

    if(fitBackground_){
      topBKGGammaNorm  ->setConstant(kFALSE);
      topBKGGammaGamma ->setConstant(kFALSE);
      topBKGGammaBeta  ->setConstant(kFALSE);
      topBKGGammaMu    ->setConstant(kFALSE);
      topBKGLandauMean ->setConstant(kFALSE);
      topBKGLandauSigma->setConstant(kFALSE);

      RooAbsReal* nllTopBKG = topBKG->createNLL(*BKG);
      RooMinuit miniTopBKG(*nllTopBKG);
      miniTopBKG.setStrategy(2);
      miniTopBKG.migrad();
      miniTopBKG.improve();
      miniTopBKG.hesse();
      RooFitResult* result = miniTopBKG.save();
      workspace[0]->import(*result);

      topBKGGammaNorm  ->setConstant(kTRUE);
      topBKGGammaGamma ->setConstant(kTRUE);
      topBKGGammaBeta  ->setConstant(kTRUE);
      topBKGGammaMu    ->setConstant(kTRUE);
      topBKGLandauMean ->setConstant(kTRUE);
      topBKGLandauSigma->setConstant(kTRUE);
    }
    //workspace[0]->import(topBKG);

    RooRealVar *wBKGMean       = new RooRealVar("wBKGMean"      , "wBKGMean"      , 86.9853 );
    RooRealVar *wBKGSigmaLeft  = new RooRealVar("wBKGSigmaLeft" , "wBKGSigmaLeft" ,  5.78569);
    RooRealVar *wBKGSigmaRight = new RooRealVar("wBKGSigmaRight", "wBKGSigmaRight",  7.12755);
    RooBifurGauss *wBKG = new RooBifurGauss("wBKG", "wBKG", meanMW, *wBKGMean, *wBKGSigmaLeft, *wBKGSigmaRight);
    wBKG->setNormRange("mTopFitRange");

    if(fitBackground_){
      wBKGMean      ->setConstant(kFALSE);
      wBKGSigmaLeft ->setConstant(kFALSE);
      wBKGSigmaRight->setConstant(kFALSE);

      RooAbsReal* nllWBKG = wBKG->createNLL(*BKG);
      RooMinuit miniWBKG(*nllWBKG);
      miniWBKG.setStrategy(2);
      miniWBKG.migrad();
      miniWBKG.improve();
      miniWBKG.hesse();
      RooFitResult* result = miniWBKG.save();
      workspace[0]->import(*result);

      wBKGMean      ->setConstant(kTRUE);
      wBKGSigmaLeft ->setConstant(kTRUE);
      wBKGSigmaRight->setConstant(kTRUE);
    }
    //workspace[0]->import(wBKG);

    RooRealVar *fSig = new RooRealVar("fSig", "fSig", 0.539, 0, 1);

    RooAddPdf *topModel = new RooAddPdf("topModel", "topModel", *topAdd, *topBKG, *fSig);
    RooAddPdf *wModel   = new RooAddPdf("wModel"  , "wModel"  , *wAdd  , *wBKG  , *fSig);
    topModel->setNormRange("mTopFitRange");
    wModel  ->setNormRange("mTopFitRange");
    //workspace[0]->import(topModel);
    //workspace[0]->import(wModel);

    RooProdPdf *model = new RooProdPdf("model", "model", *topModel, *wModel);
    model->setNormRange("mTopFitRange");
    workspace[0]->import(*model);

    workspace[0]->var("JES" )->setConstant(kFALSE);
    workspace[0]->var("mTop")->setConstant(kFALSE);
    workspace[0]->var("fSig")->setConstant(kTRUE);

    workspace[0]->writeToFile("RooWorkspace_TEST.root");
  }

  if(doMeasurement_){
    TFile * workspaceFile = TFile::Open("RooWorkspace_TEST.root");
    workspace[0] = (RooWorkspace*)workspaceFile->Get("workspaceMtop");

    //MassOffset       = *workspace[0]->var("MassOffset"      ); MassOffset      .setVal( 1.05709e-02);
    //MassSlopeMass    = *workspace[0]->var("MassSlopeMass"   ); MassSlopeMass   .setVal( 6.41444e-02);
    //MassSlopeJES     = *workspace[0]->var("MassSlopeJES"    ); MassSlopeJES    .setVal(-4.61594e+00);
    //MassSlopeMassJES = *workspace[0]->var("MassSlopeMassJES"); MassSlopeMassJES.setVal( 0.         );
    //JESOffset        = *workspace[0]->var("JESOffset"       ); JESOffset       .setVal(-3.76801e-03);
    //JESSlopeMass     = *workspace[0]->var("JESSlopeMass"    ); JESSlopeMass    .setVal(-1.50418e-04);
    //JESSlopeJES      = *workspace[0]->var("JESSlopeJES"     ); JESSlopeJES     .setVal( 4.60147e-02);
    //JESSlopeMassJES  = *workspace[0]->var("JESSlopeMassJES" ); JESSlopeMassJES .setVal( 0.         );

    //workspace[0]->var("MassOffset"      )->setVal( 1.05709e-02);
    //workspace[0]->var("MassSlopeMass"   )->setVal( 6.41444e-02);
    //workspace[0]->var("MassSlopeJES"    )->setVal(-4.61594e+00);
    //workspace[0]->var("MassSlopeMassJES")->setVal( 0.         );
    //workspace[0]->var("JESOffset"       )->setVal(-3.76801e-03);
    //workspace[0]->var("JESSlopeMass"    )->setVal(-1.50418e-04);
    //workspace[0]->var("JESSlopeJES"     )->setVal( 4.60147e-02);
    //workspace[0]->var("JESSlopeMassJES" )->setVal( 0.         );

    //workspace[0]->var("MassOffset"      )->setVal( -0.71);
    //workspace[0]->var("MassSlopeMass"   )->setVal(  0.  );
    //workspace[0]->var("MassSlopeJES"    )->setVal( 12.59);
    //workspace[0]->var("MassSlopeMassJES")->setVal(  0.  );
    //workspace[0]->var("JESOffset"       )->setVal( 0.0199);
    //workspace[0]->var("JESSlopeMass"    )->setVal( 0.    );
    //workspace[0]->var("JESSlopeJES"     )->setVal( 0.8281);
    //workspace[0]->var("JESSlopeMassJES" )->setVal( 0.    );

    //mTop_corrected   = *(RooFormulaVar*)workspace[0]->var("mTop_corrected");
    //JES_corrected    = *(RooFormulaVar*)workspace[0]->var( "JES_corrected");

    topAdd = (RooAddPdf*)workspace[0]->pdf("mTopPDF");
    wAdd   = (RooAddPdf*)workspace[0]->pdf("mWPDF"  );

    varSet  = RooArgSet(/*prob,*/MTOP,meanMW,"varSet");

    TString fileName = "writeFullHadTree_data_2011.root";
    //TString fileName = "Z2_F11_172_5_sig.root";
    //TFile* file = TFile::Open(samplePath_+fileName);
    //TTree* oldTree = (TTree*)file->Get("FullHadTreeWriter/tree");
    TChain* oldTree = new TChain("FullHadTreeWriter/tree");
    oldTree->Add(samplePath_+fileName);
    tmpFile->cd();
    bool isData = fileName.Contains("data") ? true : false;
    TTree* tree = modifiedTree_(oldTree, -10, 10, isData);

    RooDataSet* data = new RooDataSet("data","data",varSet,RooFit::Import(*tree));//,RooFit::WeightVar("prob"));

    name = "";

    RooAbsPdf* model  = workspace[0]->pdf("model");
    RooAbsPdf* topBKG = workspace[0]->pdf("topBKG");
    RooAbsPdf* wBKG   = workspace[0]->pdf("wBKG");

    RooAbsReal* fSig = workspace[0]->var("fSig");

    RooAbsReal* nll = model->createNLL(*data);
    RooMinuit mini(*nll);
    mini.setStrategy(2);
    mini.migrad();
    mini.improve();
    mini.hesse();

    TCanvas* canvas1 = new TCanvas("canvas1", "canvas1", 1, 1, 600, 600);
    canvas1->cd();
    RooPlot* frame1 = MTOP.frame(RooFit::Range(100,350));
    data  ->statOn (frame1, RooFit::Layout(.5, .9, .9));
    model ->paramOn(frame1, RooFit::Layout(.5, .9, .7));
    data  ->plotOn(frame1);
    model ->plotOn(frame1, RooFit::LineColor(8));
    topAdd->plotOn(frame1, RooFit::LineColor(kRed) , RooFit::Normalization(   fSig->getVal()));
    topBKG->plotOn(frame1, RooFit::LineColor(kBlue), RooFit::Normalization(1.-fSig->getVal()));
    //workspace[0]->pdf("sig_0")   ->plotOn(frame1, RooFit::LineColor(kRed+1), RooFit::Normalization(1./*fSig.getVal()*workspace[0]->var("fCP")->getVal()*/));
    //workspace[0]->pdf("sig_1")   ->plotOn(frame1, RooFit::LineColor(kRed+2), RooFit::Normalization(1./*fSig.getVal()*workspace[0]->var("fWP")->getVal()*/));
    //workspace[0]->pdf("sig_2")   ->plotOn(frame1, RooFit::LineColor(kRed+3), RooFit::Normalization(1./*fSig.getVal()*workspace[0]->var("fUN")->getVal()*/));
    frame1->Draw();

    //TCanvas* canvas2 = new TCanvas("canvas2", "canvas2", 641, 1, 600, 600);
    //canvas2->cd();
    //
    //RooPlot* frame2 = fSig.frame();
    ////nll->plotOn(frame2,RooFit::ShiftToZero());
    //
    //RooAbsReal* pll_fSig = nll->createProfile(fSig);
    //pll_fSig->plotOn(frame2,RooFit::LineColor(kRed));
    //frame2->Draw();
    //
    //TCanvas* canvas3 = new TCanvas("canvas3", "canvas3", 641, 1, 600, 600);
    //canvas3->cd();
    //
    //RooPlot* frame3 = workspace[0]->var("mTop")->frame();//RooFit::Bins(100),RooFit::Range(172,175));
    ////nll->plotOn(frame3,RooFit::ShiftToZero());
    //
    //RooAbsReal* pll_mTop = nll->createProfile(*workspace[0]->var("mTop"));
    //pll_mTop->plotOn(frame3,RooFit::LineColor(kRed));
    //frame3->Draw();
    //
    //TCanvas* canvas4 = new TCanvas("canvas4", "canvas4", 641, 1, 600, 600);
    //canvas4->cd();
    //
    //RooPlot* frame4 = workspace[0]->var("JES")->frame();//RooFit::Bins(100),RooFit::Range(172,175));
    ////nll->plotOn(frame4,RooFit::ShiftToZero());
    //
    //RooAbsReal* pll_JES = nll->createProfile(*workspace[0]->var("JES"));
    //pll_JES->plotOn(frame4,RooFit::LineColor(kRed));
    //frame4->Draw();

    TCanvas* canvas5 = new TCanvas("canvas5", "canvas5", 641, 1, 600, 600);
    canvas5->cd();
    RooPlot* frame5 = meanMW.frame(RooFit::Range(60,140));
    data ->statOn (frame5, RooFit::Layout(.5, .9, .9));
    model->paramOn(frame5, RooFit::Layout(.5, .9, .7));
    data ->plotOn(frame5);
    model->plotOn(frame5, RooFit::LineColor(8));
    wAdd ->plotOn(frame5, RooFit::LineColor(kRed) , RooFit::Normalization(   fSig->getVal()));
    wBKG ->plotOn(frame5, RooFit::LineColor(kBlue), RooFit::Normalization(1.-fSig->getVal()));
    //workspace[0]->pdf("sig_3")   ->plotOn(frame5, RooFit::LineColor(kRed+1), RooFit::Normalization(1./*fSig.getVal()*workspace[0]->var("fCP")->getVal()*/));
    //workspace[0]->pdf("sig_4")   ->plotOn(frame5, RooFit::LineColor(kRed+2), RooFit::Normalization(1./*fSig.getVal()*workspace[0]->var("fWP")->getVal()*/));
    //workspace[0]->pdf("sig_5")   ->plotOn(frame5, RooFit::LineColor(kRed+3), RooFit::Normalization(1./*fSig.getVal()*workspace[0]->var("fUN")->getVal()*/));
    frame5->Draw();

    //tmpFile->Close();
  }
}

//// return the PU weights for the different samples
//double
//TopMassCalibration::calcPUWeight_(enum enumForPUWeights sample, short nPU)
//{
//  // ----
//  // default weight to be used for fall11 samples with S6 PU scenario
//  // ----
//  // 3.54fb-1 !!! precale weighted !!! 73.5mb
//  //double weightsPUFall11[] = {0.0953411, 0.0209421, 0.0389577, 0.355278, 1.3007, 1.8981, 1.95124, 1.81828, 1.63994, 1.55631, 1.54116, 1.44885, 1.24065, 0.915466, 0.580092, 0.350426, 0.209707, 0.113744, 0.0509568, 0.0183767, 0.0054363, 0.00141399, 0.000360394, 0.000106747, 3.96548e-05, 1.65259e-05, 6.57705e-06, 2.32168e-06, 7.15271e-07, 1.99901e-07, 5.67217e-08, 1.86344e-08, 7.44984e-09, 3.25689e-09, 1.39639e-09, 5.58487e-10, 2.00145e-10, 6.60545e-11, 1.89381e-11, 4.82982e-12, 1.12634e-12, 2.21706e-13, 4.06965e-14, 6.1115e-15, 9.20171e-16, 1.05937e-16, 1.38254e-17, 1.14592e-18, 1.13769e-19, 43889.5};
//  // 3.54fb-1 !!! precale weighted !!! 68mb
//  double weightsPUFall11[] = {0.0956436, 0.0219907, 0.0713245, 0.735605, 1.85831, 2.19452, 2.04027, 1.76898, 1.58702, 1.53704, 1.45119, 1.24648, 0.904423, 0.552217, 0.319224, 0.179453, 0.0864062, 0.0326249, 0.0095542, 0.00229068, 0.000501167, 0.000121599, 3.80191e-05, 1.3964e-05, 4.89459e-06, 1.48493e-06, 3.82906e-07, 8.94646e-08, 2.18024e-08, 6.46154e-09, 2.31321e-09, 8.61013e-10, 3.0218e-10, 9.44128e-11, 2.57874e-11, 6.17332e-12, 1.27288e-12, 2.34488e-13, 3.65422e-14, 4.94041e-15, 5.96064e-16, 5.92555e-17, 5.36326e-18, 3.87753e-19, 2.74435e-20, 1.45014e-21, 8.48143e-23, 3.07617e-24, 1.3049e-25, 58771.4};
//  // ----
//  // fall11 weight with S6 PU scenario shifted up by 5%
//  // ----
//  // 3.54fb-1 !!! precale weighted !!! 68mb
//  double weightsPUFall11Plus05[] = {0.0954358, 0.020682, 0.0484447, 0.47369, 1.50053, 2.01825, 1.99328, 1.80561, 1.61769, 1.54968, 1.51784, 1.3867, 1.12852, 0.773361, 0.467097, 0.277645, 0.157531, 0.0761859, 0.02937, 0.00905285, 0.00233558, 0.000560472, 0.00014654, 4.85625e-05, 1.90022e-05, 7.38414e-06, 2.54485e-06, 7.59523e-07, 2.01995e-07, 5.28236e-08, 1.59273e-08, 5.8689e-09, 2.42327e-09, 9.83602e-10, 3.67109e-10, 1.23672e-10, 3.6646e-11, 9.87495e-12, 2.28817e-12, 4.67308e-13, 8.65061e-14, 1.34004e-14, 1.91938e-15, 2.23009e-16, 2.57593e-17, 2.25592e-18, 2.2207e-19, 1.37666e-20, 1.01363e-21, 50283.5};
//
//  // ----
//  // fall11 weight with S6 PU scenario shifted down by 5%
//  // ----
//  // 3.54fb-1 !!! precale weighted !!! 68mb
//  double weightsPUFall11Minus05[] = {0.0962068, 0.0254305, 0.10942, 1.09375, 2.23819, 2.34231, 2.05248, 1.72372, 1.56859, 1.50191, 1.34308, 1.04368, 0.658236, 0.373059, 0.206629, 0.0997198, 0.0369724, 0.0103136, 0.00227993, 0.000454437, 0.000101297, 3.00008e-05, 1.01979e-05, 3.2025e-06, 8.38403e-07, 1.8544e-07, 3.81727e-08, 8.90286e-09, 2.62892e-09, 8.75387e-10, 2.83582e-10, 8.20241e-11, 2.07262e-11, 4.47805e-12, 8.2373e-13, 1.30017e-13, 1.73391e-14, 2.0282e-15, 1.9709e-16, 1.63193e-17, 1.18443e-18, 6.95731e-20, 3.6548e-21, 1.50639e-22, 5.9703e-24, 1.73528e-25, 5.48351e-27, 1.0555e-28, 2.33406e-30, 73177.8};
//
//  // ----
//  // default weight to be used for summer11 samples with S4 PU scenario
//  // ----
//  // 3.54 fb-1 !!! precale weighted !!! 73.5mb
//  //double weightsPUSummer11[] = {0.015788, 0.168999, 0.398126, 0.720134, 1.04936, 1.31397, 1.48381, 1.5613, 1.57478, 1.54969, 1.50736, 1.46098, 1.41641, 1.3734, 1.33783, 1.30218, 1.26742, 1.23224, 1.19459, 1.15701, 1.11803, 1.07132, 1.02626, 0.982213, 0.936123, 0.886997, 0.840387, 0.783362, 0.7366, 0.697965, 0.645112, 0.586587, 0.536603, 0.514893, 0.46368};
//  // 3.54fb-1 !!! precale weighted !!! 68mb
//  double weightsPUSummer11[] = {0.022778, 0.231718, 0.517403, 0.888432, 1.23288, 1.47585, 1.59957, 1.62107, 1.57897, 1.50264, 1.41361, 1.32374, 1.23758, 1.15455, 1.07948, 1.00629, 0.936187, 0.868571, 0.802406, 0.739714, 0.679652, 0.618656, 0.562462, 0.510472, 0.460945, 0.413436, 0.370475, 0.326335, 0.289734, 0.259019, 0.225714, 0.193379, 0.166592, 0.150473, 0.127517};
//  // ----
//  // summer11 weight with S6 PU scenario shifted up by 5%
//  // ----
//  // 3.54fb-1 !!! precale weighted !!! 68mb
//  double weightsPUSummer11Plus05[] = {0.0181301, 0.190491, 0.439904, 0.780434, 1.11676, 1.37519, 1.52946, 1.58712, 1.58037, 1.53624, 1.47625, 1.4131, 1.35212, 1.29288, 1.24081, 1.18892, 1.13828, 1.08789, 1.03619, 0.985577, 0.93491, 0.879107, 0.82611, 0.775364, 0.724452, 0.67272, 0.624431, 0.570059, 0.524814, 0.486735, 0.440211, 0.391576, 0.350349, 0.32874, 0.289457};
//
//  // ----
//  // summer11 weight with S6 PU scenario shifted down by 5%
//  // ----
//  // 3.54fb-1 !!! precale weighted !!! 68mb
//  double weightsPUSummer11Minus05[] = {0.0285702, 0.28113, 0.606747, 1.0079, 1.35549, 1.57602, 1.66289, 1.64388, 1.56397, 1.45448, 1.33669, 1.22153, 1.11292, 1.01022, 0.917654, 0.829964, 0.748293, 0.672149, 0.600689, 0.535311, 0.475158, 0.417591, 0.366348, 0.320643, 0.279063, 0.241112, 0.208012, 0.176312, 0.150554, 0.12939, 0.108351, 0.0891757, 0.0737799, 0.0639892, 0.0520639};
//
//  switch(sample){
//  case kFall11:
//    return (nPU < int(sizeof(weightsPUFall11)/sizeof(double)) && nPU > -1) ? weightsPUFall11[nPU] : 0. ;
//  case kFall11Plus05:
//    return (nPU < int(sizeof(weightsPUFall11Plus05)/sizeof(double)) && nPU > -1) ? weightsPUFall11Plus05[nPU] : 0. ;
//  case kFall11Minus05:
//    return (nPU < int(sizeof(weightsPUFall11Minus05)/sizeof(double)) && nPU > -1) ? weightsPUFall11Minus05[nPU] : 0. ;
//  case kSummer11:
//    return (nPU < int(sizeof(weightsPUSummer11)/sizeof(double)) && nPU > -1) ? weightsPUSummer11[nPU] : 0. ;
//  case kSummer11Plus05:
//    return (nPU < int(sizeof(weightsPUSummer11Plus05)/sizeof(double)) && nPU > -1) ? weightsPUSummer11Plus05[nPU] : 0. ;
//  case kSummer11Minus05:
//    return (nPU < int(sizeof(weightsPUSummer11Minus05)/sizeof(double)) && nPU > -1) ? weightsPUSummer11Minus05[nPU] : 0. ;
//  case kFall10:
//    return 1.;
//  case kSummer12:
//    return -1.;
//  default:
//    return -1.;
//  }
//}
//
//// calculate the probability of b-tagging one event with 2 b-tags
//double
//TopMassCalibration::eventBTagProbability_(std::vector<double> &oneMinusBEffies, std::vector<double> &oneMinusBMistags){
//  double bTaggingEfficiency = 1.;
//  double tmp = 1.;
//
//  //std::cout << bTaggingEfficiency << std::endl;
//
//  // subtract probability that no jet is tagged
//  for(std::vector<double>::const_iterator eff = oneMinusBEffies.begin(); eff != oneMinusBEffies.end(); ++eff)
//    tmp *= (*eff);
//  for(std::vector<double>::const_iterator mis = oneMinusBMistags.begin(); mis != oneMinusBMistags.end(); ++mis)
//    tmp *= (*mis);
//  bTaggingEfficiency -= tmp;
//
//  //std::cout << bTaggingEfficiency << std::endl;
//
//  // subtract probability that 1 bJet is tagged
//  for(std::vector<double>::const_iterator eff = oneMinusBEffies.begin(); eff != oneMinusBEffies.end(); ++eff){
//    tmp = 1.-(*eff);
//    for(std::vector<double>::const_iterator eff2 = oneMinusBEffies.begin(); eff2 != oneMinusBEffies.end(); ++eff2)
//      if(eff != eff2) tmp *= (*eff2);
//    for(std::vector<double>::const_iterator mis = oneMinusBMistags.begin(); mis != oneMinusBMistags.end(); ++mis)
//      tmp *= (*mis);
//    bTaggingEfficiency -= tmp;
//  }
//
//  //std::cout << bTaggingEfficiency << std::endl;
//
//  // subtract probability that 1 non-bJet is tagged
//  for(std::vector<double>::const_iterator mis = oneMinusBMistags.begin(); mis != oneMinusBMistags.end(); ++mis){
//    tmp = 1.-(*mis);
//    for(std::vector<double>::const_iterator eff = oneMinusBEffies.begin(); eff != oneMinusBEffies.end(); ++eff)
//      tmp *= (*eff);
//    for(std::vector<double>::const_iterator mis2 = oneMinusBMistags.begin(); mis2 != oneMinusBMistags.end(); ++mis2)
//      if(mis != mis2) tmp *= (*mis2);
//    bTaggingEfficiency -= tmp;
//  }
//
//  //std::cout << bTaggingEfficiency << std::endl;
//
//  return bTaggingEfficiency;
//}
//
//double
//TopMassCalibration::calcBTagWeight_(int Njet, float *bTag, short *pdgId, TClonesArray *jets)
//{
//  double bTaggingEfficiency = 0., bTaggingEfficiency_scaled = 0.;
//  //double bTaggingEfficiency_scaled_EffUp = 0., bTaggingEfficiency_scaled_EffDown = 0., bTaggingEfficiency_scaled_MisUp = 0., bTaggingEfficiency_scaled_MisDown = 0.;
//  double pt, eta, eff, effyScale_pt; //, effVariation_pt, misTagScale_pt, misVariation_pt;
//
//  std::vector<double> oneMinusBEffies(0) , oneMinusBEffies_scaled(0) ; //, oneMinusBEffies_scaled_EffUp(0) , oneMinusBEffies_scaled_EffDown(0) ;
//  std::vector<double> oneMinusBMistags(0), oneMinusBMistags_scaled(0); //, oneMinusBMistags_scaled_MisUp(0), oneMinusBMistags_scaled_MisDown(0);
//  for(int i = 0; i < Njet; ++i){
//    pt  = ((TLorentzVector*)jets->At(i))->Pt();
//    eta = ((TLorentzVector*)jets->At(i))->Eta();
//
//    if(pt > 670.)
//      effyScale_pt    = 0.901615*((1.+(0.552628*670.))/(1.+(0.547195*670.)));
//    if(pt < 30.)
//      effyScale_pt    = 0.901615*((1.+(0.552628*30.))/(1.+(0.547195*30.)));
//    else
//      effyScale_pt    = 0.901615*((1.+(0.552628*pt))/(1.+(0.547195*pt)));
//    //effVariation_pt = bTagEffScaleFactor->GetBinError(bTagEffScaleFactor->FindBin(pt));
//
//    if(pdgId[i] == 5 || pdgId[i] == -5){
//      eff = bTagEff_->GetBinContent(bTagEff_->FindBin(pt,std::abs(eta)));
//      oneMinusBEffies               .push_back(1.- eff);
//      oneMinusBEffies_scaled        .push_back(1.-(eff* effyScale_pt));
//      //if(isDefaultSample){
//      //  oneMinusBEffies_scaled_EffUp  .push_back(1.-(eff*(effyScale_pt+effVariation_pt)));
//      //  oneMinusBEffies_scaled_EffDown.push_back(1.-(eff*(effyScale_pt-effVariation_pt)));
//      //}
//    }
//    else if(pdgId[i] == 4 || pdgId[i] == -4){
//      eff = cTagEff_->GetBinContent(cTagEff_->FindBin(pt,std::abs(eta)));
//      oneMinusBMistags               .push_back(1.- eff);
//      oneMinusBMistags_scaled        .push_back(1.-(eff* effyScale_pt));
//      //if(pt<240){
//      //	oneMinusBMistags_scaled_MisUp  .push_back(1.-(eff*(effyScale_pt+(2*effVariation_pt))));
//      //	oneMinusBMistags_scaled_MisDown.push_back(1.-(eff*(effyScale_pt-(2*effVariation_pt))));
//      //}
//      //else{
//      //	oneMinusBMistags_scaled_MisUp  .push_back(1.-(eff*(effyScale_pt+effVariation_pt)));
//      //	oneMinusBMistags_scaled_MisDown.push_back(1.-(eff*(effyScale_pt-effVariation_pt)));
//      //}
//    }
//    else{
//      eff = lTagEff_->GetBinContent(lTagEff_->FindBin(pt,std::abs(eta)));
//      oneMinusBMistags               .push_back(1.- eff);
//      oneMinusBMistags_scaled        .push_back(1.-(eff* (((0.948463+(0.00288102*pt))+(-7.98091e-06*(pt*pt)))+(5.50157e-09*(pt*(pt*pt)))) ));
//      //if(isDefaultSample){
//      //  oneMinusBMistags_scaled_MisUp  .push_back(1.-(eff* (((0.997077+(0.00473953*pt))+(-1.34985e-05*(pt*pt)))+(1.0032e-08*(pt*(pt*pt)))) ));
//      //  oneMinusBMistags_scaled_MisDown.push_back(1.-(eff* (((0.899715+(0.00102278*pt))+(-2.46335e-06*(pt*pt)))+(9.71143e-10*(pt*(pt*pt)))) ));
//      //}
//    }
//  }
//  bTaggingEfficiency        = eventBTagProbability_(oneMinusBEffies       , oneMinusBMistags       );
//  bTaggingEfficiency_scaled = eventBTagProbability_(oneMinusBEffies_scaled, oneMinusBMistags_scaled);
//  //if(isDefaultSample){
//  //  bTaggingEfficiency_scaled_EffUp   += eventBTagProbability_(oneMinusBEffies_scaled_EffUp  , oneMinusBMistags_scaled        );
//  //  bTaggingEfficiency_scaled_EffDown += eventBTagProbability_(oneMinusBEffies_scaled_EffDown, oneMinusBMistags_scaled        );
//  //  bTaggingEfficiency_scaled_MisUp   += eventBTagProbability_(oneMinusBEffies_scaled        , oneMinusBMistags_scaled_MisUp  );
//  //  bTaggingEfficiency_scaled_MisDown += eventBTagProbability_(oneMinusBEffies_scaled        , oneMinusBMistags_scaled_MisDown);
//  //}
//  //std::cout << bTaggingEfficiency << " " << bTaggingEfficiency_scaled << " " << bTaggingEfficiency_scaled/bTaggingEfficiency << std::endl;
//  return bTaggingEfficiency_scaled/bTaggingEfficiency;
//}
//TTree*
//TopMassCalibration::modifiedTree_(TTree * tree, int minComboType, int maxComboType, bool isData)
//{
//  tree->SetBranchStatus("*", 0);
//  tree->SetBranchStatus("probs", 1);
//  tree->SetBranchStatus("dRbb", 1);
//  tree->SetBranchStatus("topMasses", 1);
//  tree->SetBranchStatus("w1Mass", 1);
//  tree->SetBranchStatus("w2Mass", 1);
//  tree->SetBranchStatus("Njet", 1);
//  tree->SetBranchStatus("jets", 1);
//  tree->SetBranchStatus("bTag_CSV", 1);
//  if(!isData){
//    tree->SetBranchStatus("comboTypes", 1);
//    tree->SetBranchStatus("partonFlavour", 1);
//  }
//  //tree->SetBranchStatus("", 1);
//
//  const int kMAX = 10000;
//  double * probs     = new double[kMAX];
//  double * topMasses = new double[kMAX];
//  double * w1Masses  = new double[kMAX];
//  double * w2Masses  = new double[kMAX];
//  float dRbbs;
//  unsigned short * comboTypes = new unsigned short[kMAX];
//  short nPU = -1;
//
//  const int kMAX_(50);
//  TString bTagAlgo_ = "CSV";
//  int Njet;
//  tree->SetBranchStatus ("Njet", 1);
//  tree->SetBranchAddress("Njet", &Njet );
//  float bTag[kMAX_];
//  tree->SetBranchStatus (TString("bTag_")+bTagAlgo_, 1);
//  tree->SetBranchAddress(TString("bTag_")+bTagAlgo_, &bTag );
//  short pdgId[kMAX_];
//  if(!isData){
//    tree->SetBranchStatus ("partonFlavour", 1);       // pdgId  *or*  partonFlavour
//    tree->SetBranchAddress("partonFlavour", &pdgId); //  pdgId  *or*  partonFlavour
//  }
//  TClonesArray * jets = new TClonesArray("TLorentzVector");
//  tree->SetBranchStatus("jets", 1);
//  tree->GetBranch("jets")->SetAutoDelete(kFALSE);
//  tree->SetBranchAddress("jets", &jets);
//
//  double PUWeight = 1.;
//  double MCWeight = 1.;
//  double BTagWeight = 1.;
//
//  tree->SetBranchAddress("probs",probs);
//  tree->SetBranchAddress("dRbb",&dRbbs);
//  tree->SetBranchAddress("topMasses",topMasses);
//  tree->SetBranchAddress("w1Mass",w1Masses);
//  tree->SetBranchAddress("w2Mass",w2Masses);
//  if(!isData){
//    tree->SetBranchAddress("comboTypes",comboTypes);
//  }
//  //tree->SetBranchAddress("",);
//
//  // event weights
//  TString filename = tree->GetCurrentFile()->GetName();
//  enumForPUWeights whichSample = kFall10;
//  if(filename.Contains("S11")) whichSample = kSummer11;
//  else if(filename.Contains("F11")) whichSample = kFall11;
//
//  if(whichSample == kSummer11){
//    tree->SetBranchStatus("nPU", 1);
//    tree->SetBranchAddress("nPU", &nPU);
//  }
//  else if(whichSample == kFall11){
//    tree->SetBranchStatus("nPUTru", 1);
//    tree->SetBranchAddress("nPUTru", &nPU);
//  }
//  tree->SetBranchStatus("MCweight", 1);
//  tree->SetBranchAddress("MCweight", &MCWeight);
//
//
//  double prob,dRbb,topMass,meanWMass,combinedWeight,comboType; //w1Mass,w2Mass;
//  TTree * newTree = new TTree("tree","tree");
//  newTree->Branch("prob", &prob, "prob/D");
//  newTree->Branch("topMass", &topMass, "topMass/D");
//  newTree->Branch("meanWMass", &meanWMass, "meanWMass/D");
//  if(!isData){
//    newTree->Branch("combinedWeight", &combinedWeight, "combinedWeight/D");
//    newTree->Branch("comboType"     , &comboType     , "comboType/D");
//  }
//  //newTree->Branch("dRbb", &dRbb, "dRbb/D");
//  //newTree->Branch("w1Mass", &w1Mass, "w1Mass/D");
//  //newTree->Branch("w2Mass", &w2Mass, "w2Mass/D");
//  //newTree->Branch("", &);
//
//  for(int i = 0, l = tree->GetEntries(); i < l; ++i){
//    tree->GetEntry(i);
//    if(probs[0] < 0.09) continue;
//    if(dRbbs < 1.5) continue;
//    if(topMasses[0] < 100.0) continue;
//    if(topMasses[0] > 550.0) continue;
//    if(!isData && comboTypes[0] < minComboType) continue;
//    if(!isData && comboTypes[0] > maxComboType) continue;
//    prob = probs[0];
//    dRbb = dRbbs;
//    topMass = topMasses[0];
//    meanWMass = (w1Masses[0]+w2Masses[0])/2.;
//    if(!isData){
//      PUWeight = calcPUWeight_(whichSample, nPU);
//      BTagWeight = calcBTagWeight_(Njet, bTag, pdgId, jets);
//      combinedWeight = probs[0] * PUWeight * MCWeight * BTagWeight;
//      comboType = comboTypes[0];
//    }
//    //combinedWeight = 1.;
//    //w1Mass = w1Masses[0];
//    //w2Mass = w2Masses[0];
//    newTree->Fill();
//  }
//  return newTree;
//}

// because RooFit only likes plain trees with standard data types (int, float, double, ...)
// the original tree has to be adapted for the new content
TTree*
TopMassCalibration::modifiedTree_(TChain *tree, int minComboType, int maxComboType, bool isData)
{
  //TTree* treeCopy = tree->CopyTree(selection_,"",1000);
  tree->Draw(">>selectedEntries", selection_); //,"",100000);
  TEventList *selectedEntries = (TEventList*)gDirectory->Get("selectedEntries");

  tree->SetBranchStatus("*", 0);
  //tree->SetBranchStatus("jet.*"   , 1);
  tree->SetBranchStatus("top.*"   , 1);
  tree->SetBranchStatus("weight.*", 1);

  //JetEvent    *jetEvent    = new JetEvent();
  TopEvent    *topEvent    = new TopEvent();
  WeightEvent *weightEvent = new WeightEvent();

  //tree->SetBranchAddress("jet."   , &jetEvent);
  tree->SetBranchAddress("top."   , &topEvent);
  tree->SetBranchAddress("weight.", &weightEvent);

  double /*prob,dRbb,*/topMass,meanWMass,combinedWeight,comboType; //w1Mass,w2Mass;
  TTree *newTree = new TTree("tree","tree");
  //newTree->Branch("prob", &prob, "prob/D");
  newTree->Branch("topMass", &topMass, "topMass/D");
  newTree->Branch("meanWMass", &meanWMass, "meanWMass/D");
  if(!isData){
    newTree->Branch("combinedWeight", &combinedWeight, "combinedWeight/D");
    newTree->Branch("comboType"     , &comboType     , "comboType/D");
  }
  //newTree->Branch("dRbb", &dRbb, "dRbb/D");
  //newTree->Branch("w1Mass", &w1Mass, "w1Mass/D");
  //newTree->Branch("w2Mass", &w2Mass, "w2Mass/D");
  //newTree->Branch("", &);

  //for(int i = 0, l = tree->GetEntries(); i < l; ++i){
  for(int i = 0, l = selectedEntries->GetN(); i < l; ++i){
    //jetEvent->init();
    topEvent->init();
    weightEvent->init();
    tree->GetEntry(selectedEntries->GetEntry(i));
    //selectedEntries->Next();
    //if(!isData && topEvent->combinationType[0] < minComboType) continue;
    //if(!isData && topEvent->combinationType[0] > maxComboType) continue;
    //prob = topEvent->fitProb[0];
    //dRbb = topEvent->fitB1[0].DeltaR(topEvent->fitB2[0]);
    topMass = topEvent->fitTop1[0].M();
    meanWMass = (topEvent->recoW1[0].M()+topEvent->recoW2[0].M())/2.;
    if(!isData){
      //double PUWeight   = weightEvent->puWeight;
      //double BTagWeight = weightEvent->bTagEffWeight; // calcBTagWeight_(Njet, bTag, pdgId, jets);
      //double MCWeight   = weightEvent->mcWeight;
      //combinedWeight = prob * PUWeight * MCWeight * BTagWeight;
      combinedWeight = /*prob * */weightEvent->combinedWeight;
      comboType = topEvent->combinationType[0];
    }
    //combinedWeight = 1.;
    //w1Mass = w1Masses[0];
    //w2Mass = w2Masses[0];
    //std::cout << comboType << " " << topMass << " " << meanWMass << std::endl;
    newTree->Fill();
  }
  return newTree;
}

TTree*
TopMassCalibration::modifiedTree_(TChain *tree, int minComboType, int maxComboType)
{
  return modifiedTree_(tree, minComboType, maxComboType, false);
}
TTree*
TopMassCalibration::modifiedTree_(TChain *tree, int comboType)
{
  return modifiedTree_(tree, comboType, comboType);
}
TTree*
TopMassCalibration::modifiedTree_(TChain *tree)
{
  return modifiedTree_(tree, -10, 10);
}

void TopMassCalibration::UnknownChannelAbort(){
  std::cout << "Channel name *" << fChannel_ << "* not know! Aborting program execution!" << std::endl;
  exit(1);
}

void TopMassCalibration::fillAlpha(std::vector<RooFormulaVar*>& alpha, int& h, RooArgSet argSet){
  static TString formula_alpha = "@0+@1*(@4-172.5)+(@2+@3*(@4-172.5))*(@5-1.)";
  TString varName = "alpha_"; varName += h; varName += "_"; varName += alpha.size();
  alpha.push_back(new RooFormulaVar(varName, varName, formula_alpha, argSet));
}
