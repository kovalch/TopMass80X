#include "RooFitTemplateAnalyzer.h"
#include "Helper.h"

#include <iomanip>

#include "TH1D.h"
#include "TH2D.h"
#include "TF2.h"
#include "TFile.h"
#include "TCanvas.h"
#include "TMath.h"
#include "TROOT.h"
#include "TStyle.h"
#include "TLegend.h"
#include "TDirectory.h"
#include "TSystem.h"
#include "TColor.h"
#include "TRandom3.h"

//#include "RooAddition.h"
#include "RooAddPdf.h"
#include "RooBifurGauss.h"
#include "RooCategory.h"
//#include "RooDataHist.h"
#include "RooDataSet.h"
#include "RooFitResult.h"
#include "RooFormulaVar.h"
#include "RooGamma.h"
//#include "RooCBShape.h"
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

void RooFitTemplateAnalyzer::Analyze(const TString& cuts, int i, int j, po::variables_map vm) {
  Scan(cuts, i, j, vm, "mTop_JES_fSig");
  Scan(cuts, i, j, vm, "mTop_JES");
  Scan(cuts, i, j, vm, "mTop");
  //Scan(cuts, i, j, vm, "mTop fSig");
  //Scan(cuts, i, j, vm, "mTop JES fSig");
}

void RooFitTemplateAnalyzer::Scan(const TString& cuts, int i, int j, po::variables_map vm, TString variables) 
{
  /*
  //TString samplePath("/scratch/hh/dust/naf/cms/user/eschliec/TopMass/19/");
  //TString samplePath("/scratch/hh/lustre/cms/user/eschliec/TopMass/19/");

  //bool doCalibration = false;
  //bool fitBackground = false;

  //TFile * tmpFile = TFile::Open(TString(gSystem->Getenv("TMPDIR"))+TString("/tmpFileRooFitTemplateAnalyzer.root"),"RECREATE");
  //workspace = loadedWorkspace;

  //RooWorkspace* workspace;
  //workspace = new RooWorkspace("workspaceMtop", "workspaceMtop");
  */

  const int nTemplTypes = 2;
  const int nComboTypes = 3;

  /*
  //const int nPDFs = nTemplTypes*nComboTypes;
  //
  ///// definitions for fit
  //tmpFile->cd();
  //const unsigned short nTemplates = 45;
  //
  //TString templ[nTemplates] = {"jes096mass1615","jes098mass1615","jes100mass1615","jes102mass1615","jes104mass1615",
  //			       "jes096mass1635","jes098mass1635","jes100mass1635","jes102mass1635","jes104mass1635",
  //			       "jes096mass1665","jes098mass1665","jes100mass1665","jes102mass1665","jes104mass1665",
  //			       "jes096mass1695","jes098mass1695","jes100mass1695","jes102mass1695","jes104mass1695",
  //			       "jes096mass1725","jes098mass1725","jes100mass1725","jes102mass1725","jes104mass1725",
  //			       "jes096mass1755","jes098mass1755","jes100mass1755","jes102mass1755","jes104mass1755",
  //			       "jes096mass1785","jes098mass1785","jes100mass1785","jes102mass1785","jes104mass1785",
  //			       "jes096mass1815","jes098mass1815","jes100mass1815","jes102mass1815","jes104mass1815",
  //			       "jes096mass1845","jes098mass1845","jes100mass1845","jes102mass1845","jes104mass1845"};
  //
  //TString templateSettings[nTemplates] = {"MC for JES = 0.96, mTop = 161.5 GeV","MC for JES = 0.98, mTop = 161.5 GeV","MC for JES = 1.00, mTop = 161.5 GeV","MC for JES = 1.02, mTop = 161.5 GeV","MC for JES = 1.04, mTop = 161.5 GeV",
  //					  "MC for JES = 0.96, mTop = 163.5 GeV","MC for JES = 0.98, mTop = 163.5 GeV","MC for JES = 1.00, mTop = 163.5 GeV","MC for JES = 1.02, mTop = 163.5 GeV","MC for JES = 1.04, mTop = 163.5 GeV",
  //					  "MC for JES = 0.96, mTop = 166.5 GeV","MC for JES = 0.98, mTop = 166.5 GeV","MC for JES = 1.00, mTop = 166.5 GeV","MC for JES = 1.02, mTop = 166.5 GeV","MC for JES = 1.04, mTop = 166.5 GeV",
  //					  "MC for JES = 0.96, mTop = 169.5 GeV","MC for JES = 0.98, mTop = 169.5 GeV","MC for JES = 1.00, mTop = 169.5 GeV","MC for JES = 1.02, mTop = 169.5 GeV","MC for JES = 1.04, mTop = 169.5 GeV",
  //					  "MC for JES = 0.96, mTop = 172.5 GeV","MC for JES = 0.98, mTop = 172.5 GeV","MC for JES = 1.00, mTop = 172.5 GeV","MC for JES = 1.02, mTop = 172.5 GeV","MC for JES = 1.04, mTop = 172.5 GeV",
  //					  "MC for JES = 0.96, mTop = 175.5 GeV","MC for JES = 0.98, mTop = 175.5 GeV","MC for JES = 1.00, mTop = 175.5 GeV","MC for JES = 1.02, mTop = 175.5 GeV","MC for JES = 1.04, mTop = 175.5 GeV",
  //					  "MC for JES = 0.96, mTop = 178.5 GeV","MC for JES = 0.98, mTop = 178.5 GeV","MC for JES = 1.00, mTop = 178.5 GeV","MC for JES = 1.02, mTop = 178.5 GeV","MC for JES = 1.04, mTop = 178.5 GeV",
  //					  "MC for JES = 0.96, mTop = 181.5 GeV","MC for JES = 0.98, mTop = 181.5 GeV","MC for JES = 1.00, mTop = 181.5 GeV","MC for JES = 1.02, mTop = 181.5 GeV","MC for JES = 1.04, mTop = 181.5 GeV",
  //					  "MC for JES = 0.96, mTop = 184.5 GeV","MC for JES = 0.98, mTop = 184.5 GeV","MC for JES = 1.00, mTop = 184.5 GeV","MC for JES = 1.02, mTop = 184.5 GeV","MC for JES = 1.04, mTop = 184.5 GeV"};
  //
  //RooCategory* cat_templ = new RooCategory("cat_templ", "cat_templ");
  //for(unsigned t=0; t<nTemplates; ++t)
  //  cat_templ->defineType(templ[t]);
  //
  //workspace->import(*cat_templ);
  //
  //const unsigned short nJES = 5;
  //const unsigned short nMasses = 9;
  //
  //const double jes_templ [nJES] = {0.96, 0.98, 1.00, 1.02, 1.04};
  //const double mTop_templ[nMasses] = {161.5, 163.5, 166.5, 169.5, 172.5, 175.5, 178.5, 181.5, 184.5};
  //
  //const unsigned iTemplateJES [] = {0,1,2,3,4,
  //				    0,1,2,3,4,
  //				    0,1,2,3,4,
  //				    0,1,2,3,4,
  //				    0,1,2,3,4,
  //				    0,1,2,3,4,
  //				    0,1,2,3,4,
  //				    0,1,2,3,4,
  //				    0,1,2,3,4};
  //
  //const unsigned iTemplateMass[] = {0,0,0,0,0,
  //				    1,1,1,1,1,
  //				    2,2,2,2,2,
  //				    3,3,3,3,3,
  //				    4,4,4,4,4,
  //				    5,5,5,5,5,
  //				    6,6,6,6,6,
  //				    7,7,7,7,7,
  //				    8,8,8,8,8};
  //
  ///// create datasets for later use
  //RooDataSet* dataset[nTemplates];
  //RooDataSet* reducedDataset[nTemplates];
  */

  RooRealVar comboTypeVar   = RooRealVar("comboType"     ,"comboType"  ,  0.,   6.,"");
  RooRealVar prob           = RooRealVar("prob"          ,"P(#chi^{2})",  0.,   1.,"");
  RooRealVar MTOP           = RooRealVar("topMass"       ,"m_{t}^{fit}",100., 550.,"GeV");
  RooRealVar meanMW         = RooRealVar("meanWMass"     ,"m_{W}^{rec}",  0.,7000.,"GeV");
  RooRealVar combinedWeight = RooRealVar("combinedWeight","weight"     ,  0., 100.,"");
  
  //comboTypeVar.setRange("R1",0.9,1.1);
  //comboTypeVar.setRange("R4",1.9,4.1);
  //comboTypeVar.setRange("R5",4.9,5.1);
  MTOP.setRange("mTopFitRange",100.,550.);
  
  RooAddPdf* topAdd = 0;
  RooAddPdf*   wAdd = 0;
  
  TString name = "";

  /*
  //RooArgSet varSet = RooArgSet(comboTypeVar,prob,MTOP,meanMW,combinedWeight,"varSet");
  //
  //if(doCalibration){
  //
  //  for(unsigned int iMass = 0; iMass < nMasses; ++iMass){
  //    for(unsigned int iJES = 0; iJES < nJES; ++iJES){
  //	TString fileName = "Z2_";
  //	//fileName += ((iMass == 6 || iMass == 7) ? "S11_" : "F11_" );
  //	fileName += "F11_";
  //	if(iJES == 0)
  //	  fileName += "ABS_JES_096_";
  //	else if(iJES == 1)
  //	  fileName += "ABS_JES_098_";
  //	else if(iJES == 2)
  //	  fileName += "";
  //	else if(iJES == 3)
  //	  fileName += "ABS_JES_102_";
  //	else if(iJES == 4)
  //	  fileName += "ABS_JES_104_";
  //	if     (iMass==0) fileName += "161_5";
  //	else if(iMass==1) fileName += "163_5";
  //	else if(iMass==2) fileName += "166_5";
  //	else if(iMass==3) fileName += "169_5";
  //	else if(iMass==4) fileName += "172_5";
  //	else if(iMass==5) fileName += "175_5";
  //	else if(iMass==6) fileName += "178_5";
  //	else if(iMass==7) fileName += "181_5";
  //	else if(iMass==8) fileName += "184_5";
  //	fileName += "_sig.root";
  //
  //	std::cout << "Creating RooDataSet for: " << fileName << std::endl;
  //
  //	TString path(samplePath);
  //	TFile * file = TFile::Open(path+TString(fileName));
  //	TTree * tree = (TTree*)file->Get("FullHadTreeWriter/tree");
  //	tmpFile->cd();
  //	tree = modifiedTree(tree); //, minComboType, maxComboType);
  //	file->Close();
  //	int iTempl = iMass*nJES+iJES;
  //	name = "dataset_"; name += iTempl;
  //
  //	// needed to avoid a crash
  //	double co, pr, mt, mw, cw;
  //	tree->SetBranchAddress("comboType", &co);
  //	tree->SetBranchAddress("prob", &pr);
  //	tree->SetBranchAddress("topMass", &mt);
  //	tree->SetBranchAddress("meanWMass", &mw);
  //	tree->SetBranchAddress("combinedWeight", &cw);
  //    
  //	dataset[iTempl] = new RooDataSet(name,name,varSet,RooFit::Import(*tree),RooFit::WeightVar("combinedWeight"));
  //    }
  //  }
  //
  //  RooSimWSTool* simWST[nPDFs];
  //  RooSimultaneous* sim[nPDFs];
  //  RooNLLVar* nll[nPDFs][nTemplates];
  //  RooAddition* nllSim[nPDFs];
  //
  //  RooRealVar* JES  = new RooRealVar("JES" , "JES" , 1.0, 0.8, 1.2);
  //  RooRealVar* mTop = new RooRealVar("mTop", "mTop", 172.5, 100., 550., "GeV");
  //  JES ->setConstant(kTRUE);
  //  mTop->setConstant(kTRUE);
  //
  //  RooRealVar MassOffset       = RooRealVar("MassOffset"      ,"MassOffset"      , 0.); MassOffset      .setConstant(kTRUE);
  //  RooRealVar MassSlopeMass    = RooRealVar("MassSlopeMass"   ,"MassSlopeMass"   , 0.); MassSlopeMass   .setConstant(kTRUE);
  //  RooRealVar MassSlopeJES     = RooRealVar("MassSlopeJES"    ,"MassSlopeJES"    , 0.); MassSlopeJES    .setConstant(kTRUE);
  //  RooRealVar MassSlopeMassJES = RooRealVar("MassSlopeMassJES","MassSlopeMassJES", 0.); MassSlopeMassJES.setConstant(kTRUE);
  //  RooRealVar JESOffset        = RooRealVar("JESOffset"       ,"JESOffset"       , 0.); JESOffset       .setConstant(kTRUE);
  //  RooRealVar JESSlopeMass     = RooRealVar("JESSlopeMass"    ,"JESSlopeMass"    , 0.); JESSlopeMass    .setConstant(kTRUE);
  //  RooRealVar JESSlopeJES      = RooRealVar("JESSlopeJES"     ,"JESSlopeJES"     , 0.); JESSlopeJES     .setConstant(kTRUE);
  //  RooRealVar JESSlopeMassJES  = RooRealVar("JESSlopeMassJES" ,"JESSlopeMassJES" , 0.); JESSlopeMassJES .setConstant(kTRUE);
  //  // 0 = JES, 1 = mTop, 2 = JESOffset, 3 = JESSlopeMass, 4 = JESSlopeJES, 5 = JESSlopeMassJES
  //  RooFormulaVar mTop_corrected_1 = RooFormulaVar("mTop_corrected_1", "mTop_corrected_1", "@1+@2+@3*(@1-172.5)+@4*(@0-1.0)+@5*(@1-172.5)*(@0-1.0)", RooArgSet(*JES, *mTop, MassOffset, MassSlopeMass, MassSlopeJES, MassSlopeMassJES));
  //  RooFormulaVar  JES_corrected_1 = RooFormulaVar( "JES_corrected_1",  "JES_corrected_1", "@0+@2+@3*(@1-172.5)+@4*(@0-1.0)+@5*(@1-172.5)*(@0-1.0)", RooArgSet(*JES, *mTop,  JESOffset,  JESSlopeMass,  JESSlopeJES,  JESSlopeMassJES));
  //  RooFormulaVar mTop_corrected_2 = RooFormulaVar("mTop_corrected_2", "mTop_corrected_2", "@1+@2+@3*(@1-172.5)+@4*(@0-1.0)+@5*(@1-172.5)*(@0-1.0)", RooArgSet(*JES, *mTop, MassOffset, MassSlopeMass, MassSlopeJES, MassSlopeMassJES));
  //  RooFormulaVar  JES_corrected_2 = RooFormulaVar( "JES_corrected_2",  "JES_corrected_2", "@0+@2+@3*(@1-172.5)+@4*(@0-1.0)+@5*(@1-172.5)*(@0-1.0)", RooArgSet(*JES, *mTop,  JESOffset,  JESSlopeMass,  JESSlopeJES,  JESSlopeMassJES));
  //  RooFormulaVar mTop_corrected_3 = RooFormulaVar("mTop_corrected_3", "mTop_corrected_3", "@1+@2+@3*(@1-172.5)+@4*(@0-1.0)+@5*(@1-172.5)*(@0-1.0)", RooArgSet(*JES, *mTop, MassOffset, MassSlopeMass, MassSlopeJES, MassSlopeMassJES));
  //  RooFormulaVar  JES_corrected_3 = RooFormulaVar( "JES_corrected_3",  "JES_corrected_3", "@0+@2+@3*(@1-172.5)+@4*(@0-1.0)+@5*(@1-172.5)*(@0-1.0)", RooArgSet(*JES, *mTop,  JESOffset,  JESSlopeMass,  JESSlopeJES,  JESSlopeMassJES));
  //  RooFormulaVar mTop_corrected_4 = RooFormulaVar("mTop_corrected_4", "mTop_corrected_4", "@1+@2+@3*(@1-172.5)+@4*(@0-1.0)+@5*(@1-172.5)*(@0-1.0)", RooArgSet(*JES, *mTop, MassOffset, MassSlopeMass, MassSlopeJES, MassSlopeMassJES));
  //  RooFormulaVar  JES_corrected_4 = RooFormulaVar( "JES_corrected_4",  "JES_corrected_4", "@0+@2+@3*(@1-172.5)+@4*(@0-1.0)+@5*(@1-172.5)*(@0-1.0)", RooArgSet(*JES, *mTop,  JESOffset,  JESSlopeMass,  JESSlopeJES,  JESSlopeMassJES));
  //  RooFormulaVar mTop_corrected_5 = RooFormulaVar("mTop_corrected_5", "mTop_corrected_5", "@1+@2+@3*(@1-172.5)+@4*(@0-1.0)+@5*(@1-172.5)*(@0-1.0)", RooArgSet(*JES, *mTop, MassOffset, MassSlopeMass, MassSlopeJES, MassSlopeMassJES));
  //  RooFormulaVar  JES_corrected_5 = RooFormulaVar( "JES_corrected_5",  "JES_corrected_5", "@0+@2+@3*(@1-172.5)+@4*(@0-1.0)+@5*(@1-172.5)*(@0-1.0)", RooArgSet(*JES, *mTop,  JESOffset,  JESSlopeMass,  JESSlopeJES,  JESSlopeMassJES));
  //  RooFormulaVar mTop_corrected_6 = RooFormulaVar("mTop_corrected_6", "mTop_corrected_6", "@1+@2+@3*(@1-172.5)+@4*(@0-1.0)+@5*(@1-172.5)*(@0-1.0)", RooArgSet(*JES, *mTop, MassOffset, MassSlopeMass, MassSlopeJES, MassSlopeMassJES));
  //  RooFormulaVar  JES_corrected_6 = RooFormulaVar( "JES_corrected_6",  "JES_corrected_6", "@0+@2+@3*(@1-172.5)+@4*(@0-1.0)+@5*(@1-172.5)*(@0-1.0)", RooArgSet(*JES, *mTop,  JESOffset,  JESSlopeMass,  JESSlopeJES,  JESSlopeMassJES));
  //
  //  for(int templType = 0; templType < nTemplTypes; ++templType){
  //    for(int comboType = 0; comboType < nComboTypes; ++comboType){
  //
  //	//int minComboType = minComboType_;
  //	//int maxComboType = maxComboType_;
  //	//int templType = templType_;
  //
  //	int h = nComboTypes*templType + comboType;
  //
  //	TString comboRangeName = "";
  //	int minComboType = -1;
  //	int maxComboType = -1;
  //	if(comboType == 0){
  //	  minComboType = 1;
  //	  maxComboType = 1;
  //	  comboRangeName = "R1";
  //	}
  //	else if(comboType == 1){
  //	  minComboType = 2;
  //	  maxComboType = 4;
  //	  comboRangeName = "R4";
  //	}
  //	else if(comboType == 2){
  //	  minComboType = 5;
  //	  maxComboType = 5;
  //	  comboRangeName = "R5";
  //	}
  //
  //	const unsigned nAlpha = (templType == 0 && minComboType == 1 && maxComboType == 1) ? 2 : 3;
  //
  //	std::vector<double> iniPar;
  //	if(templType == 0){
  //	  if (minComboType == 1 && maxComboType == 1){ // mTop, correct
  //	    double a[] = {172.351, 0.994749 , 85.8624 , 0.819629 ,
  //			  7.31249, 0.0549339,  6.80791, 0.0648612};
  //	    iniPar = std::vector<double>(a,a+sizeof(a)/sizeof(double));
  //	  }
  //	  else if (minComboType == 2 && maxComboType == 4){ // mTop, missing+wrong
  //	    double a[] = {177.506,  0.72069   , 72.7239 , 0.,
  //			  23.0854,  0.141175  ,  4.41442, 0.,
  //			  8.21835, -0.00774364, -2.46418, 0.};
  //	    //double a[] = {160.,  0.7 ,  80. , 0.,
  //	    //		     50., -0.5 , -20. , 0.,
  //	    //		    100., -30. , -2.46418, 0.};
  //	    iniPar = std::vector<double>(a,a+sizeof(a)/sizeof(double));
  //	  }
  //	  else if(minComboType == 5 && maxComboType == 5){ // mTop, unmatched
  //	    double a[] = {176.206,  0.697207 , 69.1265, 0.,
  //			  22.2574,  0.18319  , 16.9791, 0.,
  //			  9.61471, -0.0348847, 1.59068, 0.};
  //	    //double a[] = {160.,  0.7 ,  80. , 0.,
  //	    //		     50., -0.5 , -20. , 0.,
  //	    //		    100., -30. , -2.46418, 0.};
  //	    iniPar = std::vector<double>(a,a+sizeof(a)/sizeof(double));
  //	  }
  //	}
  //	else if(templType == 1){
  //	  if (minComboType == 1 && maxComboType == 1){ // mW, correct
  //	    double a[] = {84.4836 , -0.00235611, 57.7472 ,  0.260235,
  //			  4.3087 , -0.00226071, 12.085  ,  0.108488,
  //			  5.43718, -0.0112025 , -7.38103, -0.125601};
  //	    iniPar = std::vector<double>(a,a+sizeof(a)/sizeof(double));
  //	  }
  //	  else if (minComboType == 2 && maxComboType == 4){ // mW, missing+wrong
  //	    double a[] = {84.0586 ,  0.0314863 , 28.6706 , -0.524549,
  //			  4.31361,  0.00693565,  5.86459, -0.361966,
  //			  6.27985, -0.0230866 , -6.57438,  0.40124};
  //	    iniPar = std::vector<double>(a,a+sizeof(a)/sizeof(double));
  //	  }
  //	  else if(minComboType == 5 && maxComboType == 5){ // mW, unmatched
  //	    double a[] = {84.7383 , 0.00308977, 16.0727 ,  0.38017,
  //			  4.84492, 0.00678078,  2.43329,  0.14005,
  //			  7.04282, 0.00211119, -6.0521 , -0.17391};
  //	    iniPar = std::vector<double>(a,a+sizeof(a)/sizeof(double));
  //	  }                                                                                                                                                     
  //	}
  //	for(unsigned p=0; p<4*nAlpha; p++)
  //	  std::cout << iniPar[p] << ", ";
  //	std::cout << std::endl;
  //
  //	RooRealVar* par[nPDFs][4*nAlpha];
  //	RooFormulaVar* alpha[nPDFs][nAlpha];
  //
  //	//for(unsigned h=0; h<1; ++h) {
  //	for(unsigned p=0; p<4*nAlpha; ++p) {
  //	  name = "par_"; name += h; name += "_"; name += p;
  //	  par[h][p] = new RooRealVar(name, name, iniPar[p]);
  //	  //par[h][p] = new RooRealVar(name, name, 0.);
  //	  //par[h][p]->setVal(iniPar[p]);
  //	  par[h][p]->setConstant(kFALSE);
  //	}
  //	for(unsigned i=0; i<nAlpha; ++i) {
  //	  name = "alpha_"; name += h; name += "_"; name += i;
  //	  if(h==0){
  //	    alpha[h][i] = new RooFormulaVar(name, name, "@0+@1*(@4-172.5)+(@2+@3*(@4-172.5))*(@5-1.)",
  //					    RooArgSet(*par[h][4*i], *par[h][4*i+1], *par[h][4*i+2], *par[h][4*i+3], mTop_corrected_1, JES_corrected_1));
  //	  }
  //	  else if(h==1){
  //	    alpha[h][i] = new RooFormulaVar(name, name, "@0+@1*(@4-172.5)+(@2+@3*(@4-172.5))*(@5-1.)",
  //					    RooArgSet(*par[h][4*i], *par[h][4*i+1], *par[h][4*i+2], *par[h][4*i+3], mTop_corrected_2, JES_corrected_2));
  //	  }
  //	  else if(h==2){
  //	    alpha[h][i] = new RooFormulaVar(name, name, "@0+@1*(@4-172.5)+(@2+@3*(@4-172.5))*(@5-1.)",
  //					    RooArgSet(*par[h][4*i], *par[h][4*i+1], *par[h][4*i+2], *par[h][4*i+3], mTop_corrected_3, JES_corrected_3));
  //	  }
  //	  else if(h==3){
  //	    alpha[h][i] = new RooFormulaVar(name, name, "@0+@1*(@4-172.5)+(@2+@3*(@4-172.5))*(@5-1.)",
  //					    RooArgSet(*par[h][4*i], *par[h][4*i+1], *par[h][4*i+2], *par[h][4*i+3], mTop_corrected_4, JES_corrected_4));
  //	  }
  //	  else if(h==4){
  //	    alpha[h][i] = new RooFormulaVar(name, name, "@0+@1*(@4-172.5)+(@2+@3*(@4-172.5))*(@5-1.)",
  //					    RooArgSet(*par[h][4*i], *par[h][4*i+1], *par[h][4*i+2], *par[h][4*i+3], mTop_corrected_5, JES_corrected_5));
  //	  }
  //	  else{
  //	    alpha[h][i] = new RooFormulaVar(name, name, "@0+@1*(@4-172.5)+(@2+@3*(@4-172.5))*(@5-1.)",
  //					    RooArgSet(*par[h][4*i], *par[h][4*i+1], *par[h][4*i+2], *par[h][4*i+3], mTop_corrected_6, JES_corrected_6));
  //	  }
  //	}
  //
  //	name = "sig_"; name += h;
  //	RooRealVar * topMass   = new RooRealVar("topMass"  , "m_{t}^{fit}", 100.,  550., "GeV");
  //	RooRealVar * meanWMass = new RooRealVar("meanWMass", "m_{W}^{rec}",   0., 7000., "GeV");
  //	RooRealVar * var = 0;
  //	if(templType == 0){
  //	  var = topMass;
  //	  if(minComboType == 1 && maxComboType == 1) {
  //	    TString varName = "width_"; varName += h;
  //	    RooRealVar* VoigtWidth = new RooRealVar(varName, varName, 2.0);
  //	    VoigtWidth->setConstant(kTRUE);
  //	    RooVoigtian* voigt = new RooVoigtian(name, name, *var, *alpha[h][0], *VoigtWidth, *alpha[h][1]);
  //	    workspace->import(*voigt);
  //	  }
  //	  else if(minComboType == 2 && maxComboType == 4) {
  //	    //RooRealVar* a = new RooRealVar("a", "a", 15.);
  //	    //a->setConstant(kTRUE);
  //	    //RooRealVar* n = new RooRealVar("n", "n", 5.);
  //	    //n->setConstant(kTRUE);
  //	    //RooCBShape * cb = new RooCBShape(name, name, *var, *alpha[h][0], *alpha[h][1], *a, *n);
  //	    //workspace[h]->import(*cb);
  //	    TString varName = "ratio_"; varName += h;
  //	    RooRealVar* ratio = new RooRealVar(varName, varName, 0.9);
  //	    ratio->setConstant(kTRUE);
  //	    varName = "landau_"; varName += h;
  //	    RooLandau*   landau   = new RooLandau  (varName, varName, *var, *alpha[h][0], *alpha[h][1]);
  //	    varName = "gaussian_"; varName += h;
  //	    RooGaussian* gaussian = new RooGaussian(varName, varName, *var, *alpha[h][0], *alpha[h][2]);
  //	    RooAddPdf* add = new RooAddPdf(name, name, *landau, *gaussian, *ratio);
  //	    workspace->import(*add);
  //	  }
  //	  else if(minComboType == 5 && maxComboType == 5) {
  //	    //RooRealVar* a = new RooRealVar("a", "a", 15.);
  //	    //a->setConstant(kTRUE);
  //	    //RooRealVar* n = new RooRealVar("n", "n", 5.);
  //	    //n->setConstant(kTRUE);
  //	    //RooCBShape * cb = new RooCBShape(name, name, *var, *alpha[h][0], *alpha[h][1], *a, *n);
  //	    //workspace[h]->import(*cb);
  //	    TString varName = "ratio_"; varName += h;
  //	    RooRealVar* ratio = new RooRealVar(varName, varName, 0.8);
  //	    ratio->setConstant(kTRUE);
  //	    varName = "landau_"; varName += h;
  //	    RooLandau*   landau   = new RooLandau  (varName, varName, *var, *alpha[h][0], *alpha[h][1]);
  //	    varName = "gaussian_"; varName += h;
  //	    RooGaussian* gaussian = new RooGaussian(varName, varName, *var, *alpha[h][0], *alpha[h][2]);
  //	    RooAddPdf* add = new RooAddPdf(name, name, *landau, *gaussian, *ratio);
  //	    workspace->import(*add);
  //	  }
  //	}
  //	else if(templType == 1){
  //	  var = meanWMass;
  //	  RooBifurGauss* asymGaus = new RooBifurGauss(name, name, *var, *alpha[h][0], *alpha[h][1], *alpha[h][2]);
  //	  workspace->import(*asymGaus);
  //	}
  //
  //	simWST[h] = new RooSimWSTool(*workspace);
  //	name = "sim_"; name += h;
  //	TString name_pdf = "sig_"; name_pdf += h;
  //	sim[h] = simWST[h]->build(name, name_pdf,
  //				  RooFit::SplitParam("JES" ,"cat_templ"),
  //				  RooFit::SplitParam("mTop","cat_templ"));
  //	for(unsigned t=0; t<nTemplates; ++t) {
  //	  workspace->var("JES_" +templ[t])->setVal(jes_templ [iTemplateJES [t]]);
  //	  workspace->var("mTop_"+templ[t])->setVal(mTop_templ[iTemplateMass[t]]);
  //	  name = "nll_"; name += h; name += "_"; name += t;
  //	  name_pdf = "sig_"; name_pdf += h; name_pdf += "_" + templ[t];
  //	  //if(binnedTemplates)
  //	  //  nll[h][t] = new RooNLLVar(name, name, *workspace->pdf(name_pdf), *hist[h][t], RooFit::NumCPU(1), RooFit::Range(100., 550.));
  //	  //else
  //	  reducedDataset[t] = (RooDataSet*)dataset[t]->reduce(RooFit::CutRange(comboRangeName));
  //	  nll[h][t] = new RooNLLVar(name, name, *workspace->pdf(name_pdf), *reducedDataset[t], RooFit::NumCPU(1), RooFit::Range("mTopFitRange"));
  //	}
  //	name = "nllSim_"; name += h;
  //	RooArgSet argSet;
  //	for(unsigned t=0; t<nTemplates; ++t)
  //	  argSet.add(*nll[h][t]);
  //	nllSim[h] = new RooAddition(name, name, argSet);
  //	RooFitResult* result;
  //	gStyle->SetOptTitle(1);
  //	// perform the fit
  //	RooMinuit minuit(*nllSim[h]);
  //	//minuit.setStrategy(0);
  //	//minuit.simplex();
  //	minuit.setStrategy(2);
  //	minuit.migrad();
  //	minuit.improve();
  //	minuit.hesse();
  //	//minuit.minos();
  //	result = minuit.save();
  //	workspace->import(*result);
  //    
  //	/// control plots
  //	TString outDir = "rooFitTopMassPlots";
  //    
  //	RooPlot* frame = 0;
  //	TCanvas* canvas = new TCanvas("canvas", "canvas", 10, 10, 600, 600);
  //	TString figName = "";
  //	for(unsigned t=0; t<nTemplates; ++t) {
  //	  //frame = var_mass[h]->frame();
  //	  if     (templType == 0) frame = var->frame(RooFit::Range(100., 250.));
  //	  else if(templType == 1) frame = var->frame(RooFit::Range( 50., 120.));
  //	  //if(! binnedTemplates) {
  //	  reducedDataset[t]->statOn(frame, RooFit::Layout(.6, .9, .9));
  //	  reducedDataset[t]->plotOn(frame);
  //	  //}
  //	  //else {
  //	  //  hist[h][t]->statOn(frame, RooFit::Layout(.6, .9, .9));
  //	  //  hist[h][t]->plotOn(frame);
  //	  //}
  //	  name_pdf = "sig_"; name_pdf += h; name_pdf += "_"+templ[t];
  //	  workspace->pdf(name_pdf)->plotOn(frame, RooFit::FillColor(kGray), RooFit::VisualizeError(*result));
  //	  workspace->pdf(name_pdf)->plotOn(frame, RooFit::LineColor(kRed+1));
  //	  frame->SetMinimum(.0);
  //	  frame->SetTitle(templateSettings[t]);
  //	  frame->Draw();
  //	  figName = "template_"; figName += templType; figName += "_"; figName += maxComboType; figName += "_"; figName += t;
  //	  canvas->Print(outDir + "/" +  figName + ".eps");
  //	  canvas->Print(outDir + "/catalog.ps(");
  //	}
  //	//const int redPalette[5] = {kRed-2, kRed-6, kRed+1, kRed-8, kRed+3};
  //	const int redPalette[9] = {kRed-9, kRed-8, kRed-7, kRed-6, kRed-5, kRed-4, kRed-3, kRed-2, kRed-1};
  //	// show functions for different JES values
  //	//frame = var_mass[h]->frame();
  //	if     (templType == 0) frame = var->frame(RooFit::Range(100., 250.));
  //	else if(templType == 1) frame = var->frame(RooFit::Range( 50., 120.));
  //	for(unsigned t=20; t<25; ++t) {
  //	  name_pdf = "sig_"; name_pdf += h; name_pdf += "_"+templ[t];
  //	  //std::cout << name_pdf << std::endl;
  //	  //workspace->pdf(name_pdf)->plotOn(frame, RooFit::FillColor(kGray), RooFit::VisualizeError(*result));
  //	  workspace->pdf(name_pdf)->plotOn(frame, RooFit::LineColor(kRed));
  //	}
  //	frame->SetMinimum(.0);
  //	frame->GetYaxis()->SetTitle("Probability");
  //	frame->SetTitle("JES = 0.96, 0.98, 1.00, 1.02, 1.04 (mTop = 172.5 GeV)");
  //	frame->Draw();
  //	figName = "funcFamily_jes_"; figName += templType; figName += "_"; figName += maxComboType;
  //	canvas->Print(outDir + "/" + figName + ".eps");
  //	canvas->Print(outDir + "/catalog.ps");
  //	// show functions for different mTop values
  //	//frame = var_mass[h]->frame();
  //	if     (templType == 0) frame = var->frame(RooFit::Range(100., 250.));
  //	else if(templType == 1) frame = var->frame(RooFit::Range( 50., 120.));
  //	//name_pdf = "sig_"; name_pdf += h; name_pdf += "_"+templ[2];
  //	//workspace[h]->pdf(name_pdf)->plotOn(frame, RooFit::LineColor(kRed+1));
  //	for(unsigned t=1, i=0; t<nTemplates; t+=5, ++i) {
  //	  name_pdf = "sig_"; name_pdf += h; name_pdf += "_"+templ[t];
  //	  //std::cout << name_pdf << std::endl;
  //	  //workspace->pdf(name_pdf)->plotOn(frame, RooFit::FillColor(kGray), RooFit::VisualizeError(*result));
  //	  workspace->pdf(name_pdf)->plotOn(frame, RooFit::LineColor(kRed));
  //	}
  //	frame->SetMinimum(.0);
  //	frame->GetYaxis()->SetTitle("Probability");
  //	frame->SetTitle("mTop = 161.5, 163.5, 166.5, 169.5, 172.5, 175.5, 178.5, 181.5, 184.5 GeV (JES = 1.00)");
  //	frame->Draw();
  //	figName = "funcFamily_mass_"; figName += templType; figName += "_"; figName += maxComboType;
  //	canvas->Print(outDir + "/" + figName + ".eps");
  //	canvas->Print(outDir + "/catalog.ps");
  //	gStyle->SetOptTitle(0);
  //	// plot different alpha_i as a function of JES
  //	char yTitle[99];
  //	for(unsigned i=0; i<nAlpha; ++i) {
  //	  frame = workspace->var("JES_"+templ[0])->frame(RooFit::Range(0.9, 1.1));
  //	  name = "alpha_"; name += h; name += "_"; name += i;
  //	  if(workspace->function(name+"_"+templ[0])) {
  //	    workspace->function(name+"_"+templ[0])->plotOn(frame, RooFit::FillColor(kGray),
  //							      RooFit::VisualizeError(*result));
  //	    workspace->function(name+"_"+templ[0])->plotOn(frame, RooFit::Name("sig"),
  //							      RooFit::LineColor(kRed+1),
  //							      RooFit::FillColor(kRed+1));
  //	  }
  //	  frame->SetMinimum(.0);
  //	  if(h==0)
  //	    sprintf(yTitle, "#alpha_{%i} (m_{top} = 172.5 GeV)", i);
  //	  else
  //	    sprintf(yTitle, "#alpha_{%i}", i);
  //	  frame->GetYaxis()->SetTitle(yTitle);
  //	  frame->GetXaxis()->SetTitle("JES");
  //	  frame->Draw();
  //	  double yLegHigh = 0.31;
  //	  double yLegLow = yLegHigh-2*0.06;
  //	  TLegend legend(0.72, yLegLow, 0.92, yLegHigh);
  //	  legend.SetBorderSize(0);
  //	  legend.SetFillStyle(0);
  //	  legend.SetTextFont(42);
  //	  //legend.AddEntry("sig", "P ("+var_mass[h]->getTitle()+")", "F");
  //	  legend.AddEntry("sig", "P ("+var->getTitle()+")", "F");
  //	  legend.Draw();
  //	  figName = "alpha_"; figName += i; figName += "_jes_"; figName += templType; figName += "_"; figName += maxComboType;
  //	  canvas->Print(outDir + "/" +  figName + ".eps");
  //	  canvas->Print(outDir + "/catalog.ps");
  //	  // plot different alpha_i as a function of mTop
  //	  frame = workspace->var("mTop_"+templ[2])->frame(RooFit::Range(150., 200.));
  //	  name = "alpha_"; name += h; name += "_"; name += i;
  //	  if(workspace->function(name+"_"+templ[2])) {
  //	    workspace->function(name+"_"+templ[2])->plotOn(frame, RooFit::FillColor(kGray),
  //							      RooFit::VisualizeError(*result));
  //	    workspace->function(name+"_"+templ[2])->plotOn(frame, RooFit::LineColor(kRed+1));
  //	  }
  //	  frame->SetMinimum(.0);
  //	  sprintf(yTitle, "#alpha_{%i} (JES = 1.00)", i);
  //	  frame->GetYaxis()->SetTitle(yTitle);
  //	  frame->GetXaxis()->SetTitle("m_{top} (GeV)");
  //	  frame->Draw();
  //	  legend.Draw();
  //	  figName = "alpha_"; figName += i; figName += "_mass_"; figName += templType; figName += "_"; figName += maxComboType;
  //	  canvas->Print(outDir + "/" +  figName + ".eps");
  //	  canvas->Print(outDir + "/catalog.ps)");
  //	}
  //
  //    
  //	//////////////////////////////////////////////////////////////////////////////////////////////////////
  //
  //	std::cout << "q[0] + q[1] * (x[0]-172.5) + (q[2] + q[3] * (x[0]-172.5)) * (x[1]-1.)" << std::endl;
  //	std::cout <<  "new: ";
  //	std::string vars_helper = "";
  //	std::stringstream vars(vars_helper, std::stringstream::in|std::stringstream::out);
  //	for(unsigned int i = 0; i < nAlpha*4; ++i){
  //	  TString varName = "par_"; varName += h; varName += "_"; varName += i;
  //	  vars << workspace->var(varName)->getVal();
  //	  if(i != nAlpha*4-1) vars << ", ";
  //	}
  //	std::cout << vars.str();
  //	std::cout << std::endl;
  //
  //	std::ofstream myfile;
  //	myfile.open(outDir + "/variables.txt", ios::out | ios::app);
  //	myfile << "Template Type: " << templType << ", Min Combo Type: " << minComboType << ", Max Combo Type: " << maxComboType << "\n";
  //	myfile << vars.str() << "\n";
  //	myfile.close();
  //	std::cout << "old: ";
  //	if(templType == 0){
  //	  if     (minComboType == 1 && maxComboType == 1)
  //	    std::cout << "172.351, 0.994749, 85.8624, 0.819629, 7.31249, 0.0549339, 6.80791, 0.0648612";
  //	  else if(minComboType == 2 && maxComboType == 4)
  //	    std::cout << "177.506, 0.72069, 72.7239, 0., 23.0854, 0.141175, 4.41442, 0., 8.21835, -0.00774364, -2.46418, 0.";
  //	  else if(minComboType == 5 && maxComboType == 5)
  //	    std::cout << "176.206, 0.697207, 69.1265, 0., 22.2574, 0.18319, 16.9791, 0., 9.61471, -0.0348847, 1.59068, 0.";
  //	}
  //	else if(templType == 1){
  //	  if     (minComboType == 1 && maxComboType == 1)
  //	    std::cout << "84.4836, -0.00235611, 57.7472, 0.260235, 4.3087, -0.00226071, 12.085, 0.108488, 5.43718, -0.0112025, -7.38103, -0.125601";
  //	  else if(minComboType == 2 && maxComboType == 4)
  //	    std::cout << "84.0586, 0.0314863, 28.6706, -0.524549, 4.31361, 0.00693565, 5.86459, -0.361966, 6.27985, -0.0230866, -6.57438, 0.40124";
  //	  else if(minComboType == 5 && maxComboType == 5)
  //	    std::cout << "84.7383, 0.00308977, 16.0727, 0.38017, 4.84492, 0.00678078, 2.43329, 0.14005, 7.04282, 0.00211119, -6.0521, -0.17391";
  //	}
  //	std::cout << std::endl;
  //	//}
  //    }
  //  }
  //
  //  RooRealVar fCP = RooRealVar("fCP", "fCP", 2.79599875211715698e-01);
  //  fCP.setConstant(kTRUE);
  //  RooRealVar fWP = RooRealVar("fWP", "fWP", 4.28881868720054626e-03+1.13365799188613892e-03+2.13942706584930420e-01);
  //  fWP.setConstant(kTRUE);
  //  RooRealVar fUN = RooRealVar("fUN", "fUN", 5.01034975051879883e-01);
  //  fUN.setConstant(kTRUE);
  //  RooArgSet permutationFractions = RooArgSet(fCP,fWP,fUN,"permutationFractions");
  // 
  //  RooAbsPdf* topCP = workspace->pdf("sig_0");
  //  RooAbsPdf* topWP = workspace->pdf("sig_1");
  //  RooAbsPdf* topUN = workspace->pdf("sig_2");
  //  topCP->setNormRange("mTopFitRange");
  //  topWP->setNormRange("mTopFitRange");
  //  topUN->setNormRange("mTopFitRange");
  //  RooArgSet topPDFs =  RooArgSet(*topCP,*topWP,*topUN,"topPDFs");
  //  topAdd = new RooAddPdf("mTopPDF", "mTopPDF", topPDFs, permutationFractions);
  //  topAdd->setNormRange("mTopFitRange");
  //  workspace->import(*topAdd);
  //
  //  RooAbsPdf* wCP = workspace->pdf("sig_3");
  //  RooAbsPdf* wWP = workspace->pdf("sig_4");
  //  RooAbsPdf* wUN = workspace->pdf("sig_5");
  //  wCP->setNormRange("mTopFitRange");
  //  wWP->setNormRange("mTopFitRange");
  //  wUN->setNormRange("mTopFitRange");
  //  RooArgSet wPDFs =  RooArgSet(*wCP,*wWP,*wUN,"wPDFs");
  //  wAdd = new RooAddPdf("mWPDF", "mWPDF", wPDFs, permutationFractions);
  //  wAdd->setNormRange("mTopFitRange");
  //  workspace->import(*wAdd);
  //
  //  workspace->writeToFile("PUT_A_REASONABLE_NAME.HERE");
  //}

  //TFile * workspaceFile = TFile::Open("RooFitCalibration.root");
  //workspace = (RooWorkspace*)workspaceFile->Get("workspaceMtop");
  */

  /*
  //MassOffset       = *workspace->var("MassOffset"      ); MassOffset      .setVal( 1.05709e-02);
  //MassSlopeMass    = *workspace->var("MassSlopeMass"   ); MassSlopeMass   .setVal( 6.41444e-02);
  //MassSlopeJES     = *workspace->var("MassSlopeJES"    ); MassSlopeJES    .setVal(-4.61594e+00);
  //MassSlopeMassJES = *workspace->var("MassSlopeMassJES"); MassSlopeMassJES.setVal( 0.         );
  //JESOffset        = *workspace->var("JESOffset"       ); JESOffset       .setVal(-3.76801e-03);
  //JESSlopeMass     = *workspace->var("JESSlopeMass"    ); JESSlopeMass    .setVal(-1.50418e-04);
  //JESSlopeJES      = *workspace->var("JESSlopeJES"     ); JESSlopeJES     .setVal( 4.60147e-02);
  //JESSlopeMassJES  = *workspace->var("JESSlopeMassJES" ); JESSlopeMassJES .setVal( 0.         );
  
  //mTop_corrected   = *(RooFormulaVar*)workspace->var("mTop_corrected");
  //JES_corrected    = *(RooFormulaVar*)workspace->var( "JES_corrected");
  */

  /*
  // no calibration
  workspace->var("MassOffset"      )->setVal(0.0);
  workspace->var("MassSlopeMass"   )->setVal(0.0);
  workspace->var("MassSlopeJES"    )->setVal(0.0);
  workspace->var("MassSlopeMassJES")->setVal(0.0);

  workspace->var("JESOffset"       )->setVal(0.0);
  workspace->var("JESSlopeMass"    )->setVal(0.0);
  workspace->var("JESSlopeJES"     )->setVal(0.0);
  workspace->var("JESSlopeMassJES" )->setVal(0.0);
  */

  // parametrization of mean and constant width
  workspace->var("MassOffset"      )->setVal( 8.43067e-01);
  workspace->var("MassSlopeMass"   )->setVal( 1.11208e-02);
  workspace->var("MassSlopeJES"    )->setVal(-8.57688e+00);
  workspace->var("MassSlopeMassJES")->setVal(-3.33429e-01);

  workspace->var("JESOffset"       )->setVal(-1.06016e-02);
  workspace->var("JESSlopeMass"    )->setVal( 3.05550e-05);
  workspace->var("JESSlopeJES"     )->setVal( 8.41141e-02);
  workspace->var("JESSlopeMassJES" )->setVal( 3.29863e-03);

  /*
  // parametrization of mean and width independent
  workspace->var("MassOffset"      )->setVal( 4.29624e-01+1.01432e-01);
  workspace->var("MassSlopeMass"   )->setVal( 8.00953e-02);
  workspace->var("MassSlopeJES"    )->setVal(-8.77033e+00);
  workspace->var("MassSlopeMassJES")->setVal(-4.08061e-01);
  
  workspace->var("JESOffset"       )->setVal(-9.38332e-03-1.21815e-03);
  workspace->var("JESSlopeMass"    )->setVal(-2.67342e-04);
  workspace->var("JESSlopeJES"     )->setVal( 1.07114e-01);
  workspace->var("JESSlopeMassJES" )->setVal( 2.72959e-03);

  // parametrization of width as const*mean
  workspace->var("MassOffset"      )->setVal( 8.74877e-01+1.48996e-01);
  workspace->var("MassSlopeMass"   )->setVal( 1.14649e-02);
  workspace->var("MassSlopeJES"    )->setVal(-1.18615e+01);
  workspace->var("MassSlopeMassJES")->setVal(-3.44788e-01);

  workspace->var("JESOffset"       )->setVal(-1.14916e-02-1.45748e-03);
  workspace->var("JESSlopeMass"    )->setVal( 2.42972e-05);
  workspace->var("JESSlopeJES"     )->setVal( 1.19464e-01);
  workspace->var("JESSlopeMassJES" )->setVal( 3.07011e-03);
  */
  
  topAdd = (RooAddPdf*)workspace->pdf("mTopPDF");
  wAdd   = (RooAddPdf*)workspace->pdf("mWPDF"  );
  
  RooArgSet varSet = RooArgSet(prob,MTOP,meanMW,"varSet");
  
  RooDataSet data = RooDataSet("data","data",varSet,RooFit::Import(*fTree),RooFit::WeightVar("prob"));
  
  name = "";
  for(int templType = 0; templType < nTemplTypes; ++templType){
    for(int comboType = 0; comboType < nComboTypes; ++comboType){
      int h = nComboTypes*templType + comboType;
      //const unsigned nAlpha = (templType == 0 && comboType == 0) ? 2 : 3;
      const unsigned nAlpha = (templType == 0 && comboType == 0) ? 2 : ((templType == 1) ? 3 : 4);
      const unsigned n4ParVars = (templType == 0 && comboType != 0) ? 2 : 1;
      const unsigned nPar   = 4*n4ParVars + nAlpha-n4ParVars;
      for(unsigned p=0; p<nPar; ++p) {
  	name = "par_"; name += h; name += "_"; name += p;
  	workspace->var(name)->setConstant(true);
  	//workspace->var(name)->removeError();
  	//workspace->var(name)->removeAsymError();
      }
    }
  }
  
  RooRealVar topBKGGammaNorm  = RooRealVar("topBKGGammaNorm" , "topBKGGammaNorm" ,  0.848667);
  RooRealVar topBKGGammaGamma = RooRealVar("topBKGGammaGamma", "topBKGGammaGamma",  4.14752 );
  RooRealVar topBKGGammaBeta  = RooRealVar("topBKGGammaBeta" , "topBKGGammaBeta" , 33.7793  );
  RooRealVar topBKGGammaMu    = RooRealVar("topBKGGammaMu"   , "topBKGGammaMu"   , 89.6645  );
  RooGamma   topBKGGamma      = RooGamma  ("topBKGGamma"     , "topBKGGamma"     , MTOP, topBKGGammaGamma, topBKGGammaBeta, topBKGGammaMu);
  topBKGGamma.setNormRange("mTopFitRange");
  
  RooRealVar topBKGLandauMean  = RooRealVar("topBKGLandauMean" , "topBKGLandauMean" , 217.114 );
  RooRealVar topBKGLandauSigma = RooRealVar("topBKGLandauSigma", "topBKGLandauSigma",  27.6076);
  RooLandau  topBKGLandau      = RooLandau ("topBKGLandau"     , "topBKGLandau"     , MTOP, topBKGLandauMean, topBKGLandauSigma);
  topBKGLandau.setNormRange("mTopFitRange");
  
  RooAddPdf topBKG = RooAddPdf("topBKG", "topBKG", topBKGGamma, topBKGLandau, topBKGGammaNorm);
  topBKG.setNormRange("mTopFitRange");
  
  RooRealVar wBKGMean       = RooRealVar("wBKGMean"      , "wBKGMean"      , 86.9853 );
  RooRealVar wBKGSigmaLeft  = RooRealVar("wBKGSigmaLeft" , "wBKGSigmaLeft" ,  5.78569);
  RooRealVar wBKGSigmaRight = RooRealVar("wBKGSigmaRight", "wBKGSigmaRight",  7.12755);
  RooBifurGauss wBKG        = RooBifurGauss("wBKG", "wBKG", meanMW, wBKGMean, wBKGSigmaLeft, wBKGSigmaRight);
  wBKG.setNormRange("mTopFitRange");

  /*
  //RooDataSet* BKG = 0;
  //if(fitBackground){
  //  std::cout << "Creating BKG dataset" << std::endl;
  //
  //  fileName = "QCDEstimationMix_2011_NEW_skimmed2.root";
  //  file = TFile::Open(path+fileName);
  //  oldTree = (TTree*)file->Get("tree");
  //  tmpFile->cd();
  //  tree = modifiedTree(oldTree, -10, 10, isData);
  //
  //  BKG = new RooDataSet("BKG","BKG",varSet,RooFit::Import(*tree),RooFit::WeightVar("prob"));
  //
  //  std::cout << "Creating BKG PDFs" << std::endl;
  //
  //  topBKGGammaNorm  .setConstant(kFALSE);
  //  topBKGGammaGamma .setConstant(kFALSE);
  //  topBKGGammaBeta  .setConstant(kFALSE);
  //  topBKGGammaMu    .setConstant(kFALSE);
  //  topBKGLandauMean .setConstant(kFALSE);
  //  topBKGLandauSigma.setConstant(kFALSE); 
  //
  //  RooAbsReal* nllTopBKG = topBKG->createNLL(*BKG);
  //  RooMinuit miniTopBKG = RooMinuit(*nllTopBKG);
  //  miniTopBKG.setStrategy(2);
  //  miniTopBKG.migrad();
  //  miniTopBKG.improve();
  //  miniTopBKG.hesse();
  //
  //  topBKGGammaNorm  .setConstant(kTRUE);
  //  topBKGGammaGamma .setConstant(kTRUE);
  //  topBKGGammaBeta  .setConstant(kTRUE);
  //  topBKGGammaMu    .setConstant(kTRUE);
  //  topBKGLandauMean .setConstant(kTRUE);
  //  topBKGLandauSigma.setConstant(kTRUE);
  //
  //  wBKGMean      .setConstant(kFALSE);
  //  wBKGSigmaLeft .setConstant(kFALSE);
  //  wBKGSigmaRight.setConstant(kFALSE);
  //
  //  RooAbsReal* nllWBKG = wBKG->createNLL(*BKG);
  //  RooMinuit miniWBKG = RooMinuit(*nllWBKG);
  //  miniWBKG.setStrategy(2);
  //  miniWBKG.migrad();
  //  miniWBKG.improve();
  //  miniWBKG.hesse();
  //
  //  wBKGMean      .setConstant(kTRUE);
  //  wBKGSigmaLeft .setConstant(kTRUE);
  //  wBKGSigmaRight.setConstant(kTRUE);
  //}
  */

  RooRealVar fSig = RooRealVar("fSig", "fSig", 0.539, 0, 1);
  //fSig.setVal(1);
  
  RooAddPdf topModel = RooAddPdf("topModel", "topModel", *topAdd, topBKG, fSig);
  RooAddPdf wModel   = RooAddPdf("wModel"  , "wModel"  , *wAdd  , wBKG  , fSig);
  topModel.setNormRange("mTopFitRange");
  wModel  .setNormRange("mTopFitRange");
  
  RooProdPdf model = RooProdPdf("model", "model", topModel, wModel);
  model.setNormRange("mTopFitRange");
  
  workspace->var("JES" )->setVal(1.0);
  workspace->var("mTop")->setVal(172.5);

  //workspace->var("JES" )->removeError();
  //workspace->var("JES" )->removeAsymError();
  //workspace->var("mTop")->removeError();
  //workspace->var("mTop")->removeAsymError();

  workspace->var("JES" )->setError(0.0);
  workspace->var("mTop")->setError(0.0);
  workspace->var("JES" )->setAsymError(0.0,0.0);
  workspace->var("mTop")->setAsymError(0.0,0.0);

  if( variables.Contains("JES") ) workspace->var("JES" )->setConstant(false);
  else                            workspace->var("JES" )->setConstant(true);
  if( variables.Contains("mTop")) workspace->var("mTop")->setConstant(false);
  else                            workspace->var("mTop")->setConstant(true);
  if(!variables.Contains("fSig")) fSig.setConstant(true);
  else                            fSig.setConstant(false);
  
  RooAbsReal* nll = model.createNLL(data,RooFit::CloneData(false));
  RooMinuit mini(*nll);
  mini.setStrategy(2);
  mini.migrad();
  mini.migrad();
  //mini.improve();
  mini.hesse();
  
  RooFitResult* fitResult = mini.save();
  int fitStatus = fitResult->status();
  delete fitResult;

  std::cout << "FIT STATUS: " << fitStatus << std::endl;

  if(variables.Contains("fSig") && variables.Contains("JES") && variables.Contains("mTop")){
    if(fitStatus == 0){
      ffSig          = fSig.getVal();
      ffSigError     = fSig.getError();
      fMassfSig      = workspace->var("mTop")->getVal();
      fMassfSigError = workspace->var("mTop")->getError();
      fJESfSig       = workspace->var("JES" )->getVal();
      fJESfSigError  = workspace->var("JES" )->getError();
    }
    else{
      ffSig          = -1.;
      ffSigError     = -1.;
      fMassfSig      = -1.;
      fMassfSigError = -1.;
      fJESfSig       = -1.;
      fJESfSigError  = -1.;
    }
  }
  if(!variables.Contains("fSig") && variables.Contains("JES") && variables.Contains("mTop")){
    if(fitStatus == 0){
      fMass      = workspace->var("mTop")->getVal();
      fMassError = workspace->var("mTop")->getError();
      fJES       = workspace->var("JES" )->getVal();
      fJESError  = workspace->var("JES" )->getError();
    }
    else{
      fMass      = -1.;
      fMassError = -1.;
      fJES       = -1.;
      fJESError  = -1.;
    }
  }
  if(!variables.Contains("fSig") && !variables.Contains("JES") && variables.Contains("mTop")){
    if(fitStatus == 0){
      fMassConstJES      = workspace->var("mTop")->getVal();
      fMassConstJESError = workspace->var("mTop")->getError();
    }
    else{
      fMassConstJES      = -1.;
      fMassConstJESError = -1.;
    }
  }
  
  if(!vm["task"].as<std::string>().compare("sm")){
    TCanvas* canvas1 = new TCanvas("canvas1", "canvas1", 1, 1, 600, 600);
    canvas1->cd();
    RooPlot* frame1 = MTOP.frame(RooFit::Range(100,350));
    data  .statOn (frame1, RooFit::Layout(.5, .9, .9));
    model .paramOn(frame1, RooFit::Layout(.5, .9, .7));
    data   .plotOn(frame1);
    model  .plotOn(frame1, RooFit::LineColor(8));
    topAdd->plotOn(frame1, RooFit::LineColor(kRed) , RooFit::Normalization(   fSig.getVal()));
    topBKG .plotOn(frame1, RooFit::LineColor(kBlue), RooFit::Normalization(1.-fSig.getVal()));
    //workspace->pdf("sig_0")   ->plotOn(frame1, RooFit::LineColor(kRed+1), RooFit::Normalization(1./*fSig.getVal()*workspace->var("fCP")->getVal()*/));
    //workspace->pdf("sig_1")   ->plotOn(frame1, RooFit::LineColor(kRed+2), RooFit::Normalization(1./*fSig.getVal()*workspace->var("fWP")->getVal()*/));
    //workspace->pdf("sig_2")   ->plotOn(frame1, RooFit::LineColor(kRed+3), RooFit::Normalization(1./*fSig.getVal()*workspace->var("fUN")->getVal()*/));
    frame1->Draw();
    TString path("plot/RooFitTemplate/mTop_"); path+= fIdentifier; path += "_"; path += i; path += "_"; path += j; path += "_"; path += variables; path += ".eps";
    canvas1->Print(path);
  
    /*
    //TCanvas* canvas2 = new TCanvas("canvas2", "canvas2", 641, 1, 600, 600);
    //canvas2->cd();
    //
    //RooPlot* frame2 = fSig.frame();
    ////nll.plotOn(frame2,RooFit::ShiftToZero());
    //
    //RooAbsReal* pll_fSig = nll->createProfile(fSig);
    //pll_fSig->plotOn(frame2,RooFit::LineColor(kRed));
    //frame2->Draw();
    //
    //TCanvas* canvas3 = new TCanvas("canvas3", "canvas3", 641, 1, 600, 600);
    //canvas3->cd();
    //
    //RooPlot* frame3 = workspace->var("mTop")->frame();//RooFit::Bins(100),RooFit::Range(172,175));
    ////nll.plotOn(frame3,RooFit::ShiftToZero());
    //
    //RooAbsReal* pll_mTop = nll->createProfile(*workspace->var("mTop"));
    //pll_mTop->plotOn(frame3,RooFit::LineColor(kRed));
    //frame3->Draw();
    //
    //TCanvas* canvas4 = new TCanvas("canvas4", "canvas4", 641, 1, 600, 600);
    //canvas4->cd();
    //
    //RooPlot* frame4 = workspace->var("JES")->frame();//RooFit::Bins(100),RooFit::Range(172,175));
    ////nll.plotOn(frame4,RooFit::ShiftToZero());
    //
    //RooAbsReal* pll_JES = nll->createProfile(*workspace->var("JES"));
    //pll_JES->plotOn(frame4,RooFit::LineColor(kRed));
    //frame4->Draw();
    */

    TCanvas* canvas5 = new TCanvas("canvas5", "canvas5", 641, 1, 600, 600);
    canvas5->cd();
    RooPlot* frame5 = meanMW.frame(RooFit::Range(60,140));
    data .statOn (frame5, RooFit::Layout(.5, .9, .9));
    model.paramOn(frame5, RooFit::Layout(.5, .9, .7));
    data  .plotOn(frame5);
    model .plotOn(frame5, RooFit::LineColor(8));
    wAdd ->plotOn(frame5, RooFit::LineColor(kRed) , RooFit::Normalization(   fSig.getVal()));
    wBKG  .plotOn(frame5, RooFit::LineColor(kBlue), RooFit::Normalization(1.-fSig.getVal()));
    //workspace->pdf("sig_3")   ->plotOn(frame5, RooFit::LineColor(kRed+1), RooFit::Normalization(1./*fSig.getVal()*workspace->var("fCP")->getVal()*/));
    //workspace->pdf("sig_4")   ->plotOn(frame5, RooFit::LineColor(kRed+2), RooFit::Normalization(1./*fSig.getVal()*workspace->var("fWP")->getVal()*/));
    //workspace->pdf("sig_5")   ->plotOn(frame5, RooFit::LineColor(kRed+3), RooFit::Normalization(1./*fSig.getVal()*workspace->var("fUN")->getVal()*/));
    frame5->Draw();
    path = "plot/RooFitTemplate/mW_"; path+= fIdentifier; path += "_"; path += i; path += "_"; path += j; path += "_"; path += variables; path += ".eps";
    canvas5->Print(path);
  
    delete canvas1;
    delete canvas5;
    delete frame1;
    delete frame5;
  }

  delete nll;
  //delete workspace;

  //tmpFile->Close();
  //delete tmpFile;

  std::cout << "RooFitTemplateAnalyzer done" << std::endl;
}

RooWorkspace* RooFitTemplateAnalyzer::workspace(0);

RooFitTemplateAnalyzer::RooFitTemplateAnalyzer(const TString& identifier, TTree* tree) : MassAnalyzer(identifier, tree)
{
  if(!workspace){
    //std::cout << "!!! Reading in ROOWORKSPACE from file !!!" << std::endl;
    TFile * workspaceFile = TFile::Open("/afs/naf.desy.de/user/e/eschliec/wd/TopMass/RooFitCalibration.root");
    //workspaceFile->ls();
    gROOT->cd();
    workspace = (RooWorkspace*)((RooWorkspace*)workspaceFile->Get("workspaceMtop"))->Clone();
    workspaceFile->Close();
    delete workspaceFile;
  }
}

RooFitTemplateAnalyzer::~RooFitTemplateAnalyzer()
{
  //delete workspace;
}
