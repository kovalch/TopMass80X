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
      jsfValues_({0.96, 0.98, 1.00, 1.02, 1.04}),
      massValues_({166.5, 169.5, 171.5, 172.5, 173.5, 175.5, 178.5}),
      // jsfValues_({0.96, 1.00, 1.04}),
      // massValues_({166.5, 172.5, 178.5}),
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

std::string TemplateDerivation::templateName(double mass, double jes) {
  std::string name = str(boost::format("jes%3%%4$02imass%1%%2%") % (int)mass %
                         (int)((mass - int(mass)) * 10) % (int)jes %
                         round((jes - int(jes)) * 100));
  return name;
}

std::vector<RooDataSet *> TemplateDerivation::createDataSets(
    const RooArgSet *varSet) {
  TFile *tmpFile = TFile::Open("tmpFileRooFitTopMass.root", "RECREATE");
  std::vector<RooDataSet *> dataset;
  for (auto mass : massValues_) {
    for (auto jes : jsfValues_) {
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
      tree->SetBranchAddress("fitTopMass", &mt);
      tree->SetBranchAddress("meanWMass", &mw);
      tree->SetBranchAddress("combinedWeight", &cw);
      TString name = templateName(mass, jes);
      dataset.push_back(new RooDataSet(name, name, *varSet,
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

void TemplateDerivation::addTemplateFunction(const std::string &varType,
                                             const std::string &comboType,
                                             RooRealVar *mTop,
                                             RooRealVar *JSF) {
  std::cout << "build template for " << varType << " and " << comboType << '\n';
  std::vector<RooRealVar *> par;
  std::vector<double> iniPar;

  if (varType == "mTop") {
    if (comboType == "CP") {
      iniPar = {172.399, 0.989847,  81.225,  0.758306,
                7.92651, 0.0730216, 5.03205, 0.0241954};
    } else if (comboType == "WP") {
      iniPar = {190.009, 0.23107,   123.714,  -2.41287,   173.325, 1.11063,
                80.0112, 1.32969,   26.2827,  -0.0172727, 35.0065, -0.442854,
                10.0152, 0.0656903, -4.72985, -0.108491,  0.700556};
    } else if (comboType == "UN") {
      iniPar = {178.883, 0.306507,  40.7573, 0.445759,  174.681, 1.06829,
                83.8601, 0.280563,  22.065,  0.0542991, 4.93823, 0.3491,
                9.74,    0.0548014, 2.31277, -0.181041};
    }
  }
  if (varType == "mW") {
    if (comboType == "CP") {  // mW, correct
      iniPar = {84.4535, -0.0154884, 91.1495, -0.00943573, 5.20726,  -0.0115877,
                23.4646, -0.0249863, 6.75636, 0.0019251,   -19.8341, 0.0724597};
    } else if (comboType == "WP") {  // mW, wrong
      iniPar = {87.2093,    -0.00584533, 20.9648,  -0.0912001, 3.91928,
                0.00521933, 3.8585,      -0.28899, 5.9744,     0.00351353,
                -5.70177,   -0.231188,   3.80803,  0.00289055, -11.4316,
                -0.264811,  5.0,         5.0,      0.25,       0.5};
    } else if (comboType == "UN") {  // mW, unmatched
      iniPar = {86.6185, -0.0184471, 32.1242,  0.0815632,  6.48169, -0.0108763,
                7.37265, 0.0139849,  6.48169,  -0.0108763, 7.37265, 0.0139849,
                7.91363, 0.00815331, -17.3423, -0.00100213};
    }
  }
  // init parameters
  for (unsigned int i = 0; i < iniPar.size(); ++i) {
    TString name = addInt(catName("par", varType, comboType), i);
    RooRealVar myVar(name, name, iniPar[i]);
    myVar.setVal(iniPar[i]);
    myVar.setConstant(kFALSE);
    workspace_->import(myVar);
    par.push_back(workspace_->var(name));
  }

  if (varType == "mTop") {
    if (comboType == "CP") {
      createAlpha("alpha_mTop_CP_0",
                  RooArgSet(*par[0], *par[1], *par[2], *par[3], *mTop, *JSF));
      createAlpha("alpha_mTop_CP_1",
                  RooArgSet(*par[4], *par[5], *par[6], *par[7], *mTop, *JSF));
      workspace_->factory(
          "Voigtian::pdf_mTop_CP(fitTopMass,alpha_mTop_CP_0,width_mTop_CP[2],"
          "alpha_mTop_CP_1)");
    } else if (comboType == "WP") {
      createAlpha("alpha_mTop_WP_0",
                  RooArgSet(*par[0], *par[1], *par[2], *par[3], *mTop, *JSF));
      createAlpha("alpha_mTop_WP_1",
                  RooArgSet(*par[8], *par[9], *par[10], *par[11], *mTop, *JSF));
      workspace_->factory(
          "Landau::pdf_mTop_WP1(fitTopMass,alpha_mTop_WP_0,alpha_mTop_WP_1)");

      createAlpha("alpha_mTop_WP_2",
                  RooArgSet(*par[4], *par[5], *par[6], *par[7], *mTop, *JSF));
      createAlpha("alpha_mTop_WP_3", RooArgSet(*par[12], *par[13], *par[14],
                                               *par[15], *mTop, *JSF));
      workspace_->factory(
          "Gaussian::pdf_mTop_WP2(fitTopMass,alpha_mTop_WP_2,alpha_"
          "mTop_WP_3)");
      workspace_->factory(
          "SUM::pdf_mTop_WP(par_mTop_WP_16*pdf_mTop_WP1,pdf_mTop_"
          "WP2)");
    } else if (comboType == "UN") {
    }
  }
  if (varType == "mW") {
    if (comboType == "CP") {  // mW, correct
      createAlpha("alpha_mW_CP_0",
                  RooArgSet(*par[0], *par[1], *par[2], *par[3], *mTop, *JSF));
      createAlpha("alpha_mw_CP_1",
                  RooArgSet(*par[4], *par[5], *par[6], *par[7], *mTop, *JSF));
      createAlpha("alpha_mW_CP_2",
                  RooArgSet(*par[8], *par[9], *par[10], *par[11], *mTop, *JSF));
      workspace_->factory(
          "BifurGauss::pdf_mW_CP(meanWMass,alpha_mW_CP_0,alpha_mw_CP_1,alpha_"
          "mW_CP_2)");
    } else if (comboType == "WP") {  // mW, wrong
      createAlpha("alpha_mW_WP_0", RooArgSet(*par[0], *par[1], *par[2], *par[3],
                                             *mTop, *JSF, *par[16]),
                  "-@6");
      createAlpha("alpha_mW_WP_1",
                  RooArgSet(*par[0], *par[1], *par[2], *par[3], *mTop, *JSF));
      createAlpha("alpha_mW_WP_2", RooArgSet(*par[0], *par[1], *par[2], *par[3],
                                             *mTop, *JSF, *par[17]),
                  "+@6");
      createAlpha("alpha_mW_WP_3",
                  RooArgSet(*par[4], *par[5], *par[6], *par[7], *mTop, *JSF));
      createAlpha("alpha_mW_WP_4",
                  RooArgSet(*par[8], *par[9], *par[10], *par[11], *mTop, *JSF));
      createAlpha("alpha_mW_WP_5", RooArgSet(*par[12], *par[13], *par[14],
                                             *par[15], *mTop, *JSF));
      workspace_->factory(
          "Gaussian::pdf_mW_WP1(meanWMass,alpha_mW_WP_0,alpha_mW_WP_3)");
      workspace_->factory(
          "Gaussian::pdf_mW_WP2(meanWMass,alpha_mW_WP_1,alpha_mW_WP_4)");

      workspace_->factory(
          "Gaussian::pdf_mW_WP3(meanWMass,alpha_mW_WP_2,alpha_mW_WP_5)");
      workspace_->factory(
          "SUM::pdf_mW_WP(par_mW_WP_18*pdf_mW_WP1,par_mW_WP_19*pdf_mW_WP2, "
          "pdf_mW_WP3)");
      par[16]->setConstant(true);
      par[17]->setConstant(true);
      par[18]->setConstant(true);
      par[19]->setConstant(true);
    } else if (comboType == "UN") {  // mW, unmatched
    }
  }
  std::cout << "old: ";
  for (auto val : iniPar) {
    std::cout << val;
    std::cout << "|";
  }
  std::cout << std::endl;
}

RooFitResult *TemplateDerivation::fitTemplate(const std::string &varType,
                                              const std::string &comboType) {
  RooSimWSTool *simWST = new RooSimWSTool(*workspace_);
  std::string pdfName = catName("pdf", varType, comboType);
  std::string simName = pdfName + "_sim";
  RooSimultaneous *sim =
      simWST->build(simName.c_str(), pdfName.c_str(),
                    RooFit::SplitParam("JSF", "calibPoints"),
                    RooFit::SplitParam("mTop", "calibPoints"));

  // RooArgSet nllSet;
  // int templateIndex = 0;
  std::map<std::string, RooDataSet *> datamap;
  for (auto mass : massValues_) {
    for (auto jsf : jsfValues_) {
      std::string templname = templateName(mass, jsf);
      std::string name = catName(templname, varType, comboType);
      RooAbsData *dataset = workspace_->data(templname.c_str());
      workspace_->var(catName("JSF", templname).c_str())->setVal(jsf);
      workspace_->var(catName("mTop", templname).c_str())->setVal(mass);
      RooDataSet *reducedDataset =
          (RooDataSet *)dataset->reduce(RooFit::CutRange(comboType.c_str()));
      reducedDataset->SetName(name.c_str());
      std::cout << reducedDataset->GetName() << std::endl;
      workspace_->import(*reducedDataset);
      std::cout << "Entries in " << name << ": " << reducedDataset->numEntries()
                << std::endl;
      datamap[templname] = reducedDataset;
    }
  }
  std::string dataName = catName("data", varType, comboType);

  RooDataSet combData(dataName.c_str(), dataName.c_str(),
                      *(workspace_->set("varSet")),
                      RooFit::Index(*(workspace_->cat("calibPoints"))),
                      RooFit::Import(datamap));
  RooFitResult *result = sim->fitTo(
      combData, RooFit::Strategy(2), RooFit::NumCPU(4, RooFit::SimComponents),
      RooFit::Range("mTopFitRange"), RooFit::Save(true),
      RooFit::SumW2Error(true));

  std::cout << "result:" << result << '\n';
  result->Print();
  workspace_->import(*result);
  return result;
}

void TemplateDerivation::plotResult(const std::string &varType,
                                    const std::string &comboType) {
  std::string pdfName = catName("pdf", varType, comboType);
  std::string simName = pdfName + "_sim";
  RooFitResult *result = (RooFitResult *)workspace_->genobj(simName.c_str());
  RooRealVar *var = (varType == "mTop") ? workspace_->var("fitTopMass")
                                        : workspace_->var("meanWMass");
  RooPlot *frame = 0;

  TCanvas *canvas = new TCanvas("canvas", "canvas", 10, 10, 600, 600);
  for (auto mass : massValues_) {
    for (auto jsf : jsfValues_) {
      if (varType == "mTop") {
        if (comboType == "CP")
          frame = var->frame(RooFit::Range(100., 215.));
        else
          frame = var->frame(RooFit::Range(100., maxMtop_));
      } else if (varType == "mW")
        frame = var->frame(RooFit::Range(50., 120.));

      std::string templName = templateName(mass, jsf);
      std::string name = catName(templName, varType, comboType);

      RooAbsData *reducedDataset = workspace_->data(name.c_str());
      reducedDataset->statOn(frame, RooFit::Layout(.6, .9, .9));
      reducedDataset->plotOn(frame);

      std::string tempname = catName(pdfName, templName);
      std::cout << "PDF Name: " << tempname << std::endl;
      workspace_->pdf(tempname.c_str())->plotOn(
          frame, RooFit::FillColor(kGray), RooFit::VisualizeError(*result));
      workspace_->pdf(tempname.c_str())
          ->plotOn(frame, RooFit::LineColor(kRed + 1));
      frame->SetMinimum(.0);
      std::string title =
          str(boost::format("MC for mass=%1$4.1f GeV and JSF= %2$3.2f") % mass %
              jsf);
      frame->SetTitle(title.c_str());
      frame->Draw();
      std::string figName = catName("template", varType, comboType, templName);
      TString outDir = "plot/calibration";
      canvas->Print(outDir + "/" + figName + ".eps");
      canvas->Print(outDir + "/catalog.ps(");
    }
  }
  // const int redPalette[5] = {kRed-2, kRed-6, kRed+1, kRed-8,
  // kRed+3};
  // const int redPalette[9] = {kRed-9, kRed-8, kRed-7, kRed-6,
  // kRed-5,
  // kRed-4, kRed-3, kRed-2, kRed-1};
  // show functions for different JES values
  // frame = var_mass[index]->frame();

  // for(unsigned t=20; t<25; ++templateIndex) {
  if (varType == "mTop") {
    if (comboType == "CP")
      frame = var->frame(RooFit::Range(100., 215.));
    else
      frame = var->frame(RooFit::Range(100., maxMtop_));
  } else if (varType == "mW")
    frame = var->frame(RooFit::Range(50., 120.));
  for (auto jsf : jsfValues_) {
    std::string templName = templateName(172.5, jsf);
    std::string pdfname = catName("pdf", varType, comboType, templName);
    workspace_->pdf(pdfname.c_str())->plotOn(frame, RooFit::LineColor(kRed));
  }
  frame->SetMinimum(.0);
  frame->GetYaxis()->SetTitle("Probability");
  frame->SetTitle("JES = 0.96, 0.98, 1.00, 1.02, 1.04 (mTop = 172.5 GeV)");
  frame->Draw();
  std::string figName = catName("funcFamily", "jsf", varType, comboType);
  TString outDir = "plot/calibration";
  canvas->Print(outDir + "/" + figName + ".eps");
  canvas->Print(outDir + "/catalog.ps");
  // show functions for different mTop values
  // frame = var_mass[index]->frame();
  // name_pdf = "sig_"; name_pdf += h; name_pdf += "_"+templ[2];
  // workspace[h]->pdf(name_pdf)->plotOn(frame,
  // RooFit::LineColor(kRed+1));
  if (varType == "mTop") {
    if (comboType == "CP")
      frame = var->frame(RooFit::Range(100., 215.));
    else
      frame = var->frame(RooFit::Range(100., maxMtop_));
  } else if (varType == "mW")
    frame = var->frame(RooFit::Range(50., 120.));
  for (auto mass : massValues_) {
    std::string templName = templateName(mass, 1.0);
    std::string pdfname = catName("pdf", varType, comboType, templName);
    workspace_->pdf(pdfname.c_str())->plotOn(frame, RooFit::LineColor(kRed));
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
  figName = "funcFamily_mass_" + varType + "_" + comboType;
  canvas->Print(outDir + "/" + figName + ".eps");
  canvas->Print(outDir + "/catalog.ps");
  gStyle->SetOptTitle(0);
  // plot different alpha_i as a function of JES
  char yTitle[99];
  std::string alphaName = catName("alpha", varType, comboType);
  for (unsigned i = 0; i < numVariables(alphaName); ++i) {
    TString jesName("JSF_");
    jesName += templateName(massValues_[0], jsfValues_[0]);
    frame = workspace_->var(jesName)->frame(RooFit::Range(0.9, 1.1));
    TString name = addInt(alphaName, i);
    name += "_" + templateName(massValues_[0], jsfValues_[0]);
    if (workspace_->function(name)) {
      workspace_->function(name)->plotOn(frame, RooFit::FillColor(kGray),
                                         RooFit::VisualizeError(*result));
      workspace_->function(name)->plotOn(frame, RooFit::Name("sig"),
                                         RooFit::LineColor(kRed + 1),
                                         RooFit::FillColor(kRed + 1));
    }
    frame->SetMinimum(.0);
    if (varType == "mTop")
      sprintf(yTitle, "#alpha_{%i} (m_{top} = 172.5 GeV)", i);
    else
      sprintf(yTitle, "#alpha_{%i}", i);
    frame->GetYaxis()->SetTitle(yTitle);
    frame->GetXaxis()->SetTitle("JSF");
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
    figName += varType;
    figName += "_";
    figName += comboType;
    canvas->Print(outDir + "/" + figName + ".eps");
    canvas->Print(outDir + "/catalog.ps");
    // plot different alpha_i as a function of mTop
    TString mtopName("mTop_");
    mtopName += templateName(massValues_[2], jsfValues_[0]);
    frame = workspace_->var(mtopName)->frame(RooFit::Range(150., 200.));
    name = catName(addInt(catName("alpha", varType, comboType), i),
                   templateName(massValues_[2], jsfValues_[0]));
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
    figName += varType;
    figName += "_";
    figName += comboType;
    canvas->Print(outDir + "/" + figName + ".eps");
    canvas->Print(outDir + "/catalog.ps)");
  }
}

void TemplateDerivation::printResult(const std::string &varType,
                                     const std::string &comboType) {
  std::cout << "q[0] + q[1] * (x[0]-172.5) + (q[2] + q[3] * "
               "(x[0]-172.5)) * (x[1]-1.)" << std::endl;
  std::cout << "new: ";
  std::string vars_helper = "";
  std::stringstream vars(vars_helper,
                         std::stringstream::in | std::stringstream::out);
  std::string parName = catName("par", varType, comboType);
  for (unsigned int i = 0; i < numVariables(parName); ++i) {
    RooRealVar *curPar = workspace_->var(addInt(parName, i).c_str());
    curPar->setConstant(kTRUE);
    vars << curPar->getVal();
    if (i != numVariables(parName) - 1) vars << "|";
  }
  std::cout << vars.str();
  std::cout << std::endl;

  std::ofstream myfile;
  TString outDir = "plot/calibration";
  myfile.open(outDir + "/variables.txt", std::ios::out | std::ios::app);
  myfile << "Template Type: " << varType << ", Combo Type: " << comboType
         << "\n";
  myfile << vars.str() << "\n";
  myfile.close();
}

void TemplateDerivation::run() {
  workspace_ = new RooWorkspace("workspaceMtop", "workspaceMtop");
  workspace_->factory("comboType[-10,20]");
  workspace_->factory("fitTopMass[100,400]");
  workspace_->var("fitTopMass")->setMax(maxMtop_);
  workspace_->factory("meanWMass[50,300]");
  workspace_->factory("combinedWeight[-100,100]");

  workspace_->var("fitTopMass")->setRange("mTopFitRange", 100, maxMtop_);

  std::vector<std::string> permutationTypes({"CP", "WP", "UN"});
  if (channelID_ == kAllJets) {
    workspace_->var("comboType")->setRange("CP", 0.9, 1.1);
    workspace_->var("comboType")->setRange("WP", 1.9, 20.1);
    permutationTypes = {"CP", "WP"};
  } else {
    workspace_->var("comboType")->setRange("CP", 0.9, 1.1);
    workspace_->var("comboType")->setRange("UN", 5.1, 20.1);
    workspace_->var("comboType")->setRange("WP", 1.9, 5.1);
  }
  std::vector<std::string> variables({"mTop", "mW"});
  workspace_->factory().createCategory("varType", "mTop,mW");
  if (channelID_ == kAllJets) {
    workspace_->factory().createCategory("permType", "CP,WP");
  } else {
    workspace_->factory().createCategory("permType", "CP,WP,UN");
  }
  std::string categories = "";
  for (auto mass : massValues_) {
    for (auto jes : jsfValues_) {
      categories += templateName(mass, jes);
      categories += ',';
    }
  }
  categories.pop_back();  // remove last ','
  workspace_->factory().createCategory("calibPoints", categories.c_str());

  workspace_->factory("JSF[1.0, 0.8, 1.2]");
  workspace_->factory("mTop[172.5,100.,400]");
  workspace_->var("JSF")->setConstant(kTRUE);
  workspace_->var("mTop")->setConstant(kTRUE);

  // import data

  // 0 = JES, 1 = mTop, 2 = Offset, 3 = SlopeMass, 4 = SlopeJES, 5 =
  // SlopeMassJES
  //  workspace_->factory(
  // "expr::mTop_intermediate('@1+@2+@3*(@1-172.5)+@4*(@0-1.0)+@5*(@1-172.5)*("
  //      "@0-1.0)',JES,mTop,"
  // "MassOffset[0],MassSlopeMass[0],MassSlopeJES[0],MassSlopeMassJES[0])");
  //  workspace_->factory(
  // "expr::JES_intermediate('@0+@2+@3*(@1-172.5)+@4*(@0-1.0)+@5*(@1-172.5)*(@"
  //      "0-1.0)', JES, mTop,"
  //      "JESOffset[0], JESSlopeMass[0],JESSlopeJES[0],
  //      JESSlopeMassJES[0])");
  //  RooRealVar *mTop = (RooRealVar *)workspace_->factory(
  // "expr::mTop_corrected('@1+@2+@3*(@1-172.5)+@4*(@0-1.0)+@5*(@1-172.5)*(@0-"
  //      "1.0)',JES_intermediate,mTop_intermediate,MassOffset2[0],"
  //      "MassSlopeMass2[0],MassSlopeJES2[0],"
  //      "MassSlopeMassJES2[0])");
  //  RooRealVar *JES = (RooRealVar *)workspace_->factory(
  // "expr::JES_corrected('@0+@2+@3*(@1-172.5)+@4*(@0-1.0)+@5*(@1-172.5)*(@0-"
  //      "1.0)',JES_intermediate,mTop_intermediate,JESOffset2[0], "
  //      "JESSlopeMass2[0],JESSlopeJES2[0], JESSlopeMassJES2[0])");

  for (auto variable : variables) {
    for (auto permutationType : permutationTypes) {
      addTemplateFunction(variable, permutationType, workspace_->var("mTop"),
                          workspace_->var("JSF"));
    }
  }
  /// create datasets and import for later use
  workspace_->defineSet("varSet",
                        "comboType,fitTopMass,meanWMass,combinedWeight");
  std::vector<RooDataSet *> datasets =
      createDataSets(workspace_->set("varSet"));
  // fitTemplate("mTop", "WP");
  // plotResult("mTop", "WP");
  // printResult("mTop", "WP");

  // workspace_->Print();
  for (auto variable : variables) {
    for (auto permutationType : permutationTypes) {
      fitTemplate(variable, permutationType);
      plotResult(variable, permutationType);
    }
  }
  // fitTemplate("mTop", "CP");
  //  for (int templType = 0; templType < nTemplTypes; ++templType) {
  //    for (int comboType = 0; comboType < nComboTypes; ++comboType) {
  //      fitTemplate(templType, comboType, nComboTypes);
  //      // plotResult(templType, comboType, nComboTypes);
  //      printResult(templType, comboType, nComboTypes);
  //    }
  //  }

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
  newTree->Branch("fitTopMass", &topMass, "fitTopMass/D");
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

RooFormulaVar *TemplateDerivation::createAlpha(const std::string &name,
                                               const RooArgSet &argSet,
                                               const std::string &addTerm) {
  std::string cmd = "expr::";
  cmd += name + "('@0+@1*(@4-172.5)+(@2+@3*(@4-172.5))*(@5-1.)" + addTerm +
         "'," + argSet.contentsString() + ")";
  RooFormulaVar *form = (RooFormulaVar *)workspace_->factory(cmd.c_str());
  return form;
}

std::string TemplateDerivation::addInt(const std::string &name, int i) {
  return str(boost::format(name + "_%1%") % i);
}

std::string TemplateDerivation::catName(const std::string &name,
                                        const std::string &cat1,
                                        const std::string &cat2,
                                        const std::string &cat3) {
  std::string temp(name);
  temp += "_" + cat1;
  if (cat2 != "") {
    temp += "_" + cat2;
  }
  if (cat3 != "") {
    temp += "_" + cat3;
  }
  return temp;
}

unsigned int TemplateDerivation::numVariables(
    const std::string &startName) const {
  std::string temp = startName + "_";
  return workspace_->allVars().selectByName(temp.c_str())->getSize();
}
