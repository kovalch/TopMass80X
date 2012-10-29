#include "TopMass.h"

#include <cmath>
#include <fstream>
#include <map>
#include <string>
#include <time.h>

#include "TCanvas.h"
#include "TFile.h"
#include "TGraphErrors.h"
#include "TLegend.h"
#include "TMultiGraph.h"
#include "TPaveStats.h"
#include "TStyle.h"
#include "TSystem.h"

#include "Helper.h"
#include "ProgramOptionsReader.h"

typedef ProgramOptionsReader po;

TopMass::TopMass() :
  fMethod_(po::GetOption<std::string>("method")),
  fBins_  (po::GetOption<int        >("bins")),
  fLumi_  (po::GetOption<double     >("lumi"))
{
  // check existence of a temp directory and create one if not available
  TString tempDir(gSystem->Getenv("TMPDIR"));
  if(tempDir.IsNull() || !tempDir.Length()){
    tempDir = gSystem->GetFromPipe("mktemp -d");
    gSystem->Setenv("TMPDIR", tempDir);
  }
  std::cout << "Directory to be used for temporary files: " << tempDir << std::endl;

  // set environment variables needed for LHAPDF
  gSystem->Setenv("LHAPATH", "/afs/naf.desy.de/user/e/eschliec/wd/LHAPDF/share/lhapdf/PDFsets");

  //QuickCalibration();
  //LoadXML();
  //QuickSystematics();
  
  if (!po::GetOption<std::string>("task").compare("cal")) {
    QuickCalibration();
  }
  
  else if (!po::GetOption<std::string>("task").compare("pe")) {
    WriteEnsembleTest();
  }
  
  else if (!po::GetOption<std::string>("task").compare("sm")) {
    Analysis* analysis = new Analysis();
    analysis->Analyze();
  }
  
  else if (!po::GetOption<std::string>("task").compare("hc")) {
  
    /*  Systematic samples
    analyses.push_back(new Analysis("1725_flavordown", "/scratch/hh/current/cms/user/mseidel/Summer11_TTJets1725_flavor:down/analyzeTop.root", fMethod, fBins, 0));
    analyses.push_back(new Analysis("1725_flavorup", "/scratch/hh/current/cms/user/mseidel/Summer11_TTJets1725_flavor:up/analyzeTop.root", fMethod, fBins, 0));
    analyses.push_back(new Analysis("1725_jesdown", "/scratch/hh/current/cms/user/mseidel/Summer11_TTJets1725_jes:down/analyzeTop.root", fMethod, fBins, 0));
    analyses.push_back(new Analysis("1725_jesup", "/scratch/hh/current/cms/user/mseidel/Summer11_TTJets1725_jes:up/analyzeTop.root", fMethod, fBins, 0));
    analyses.push_back(new Analysis("1725_resdown", "/scratch/hh/current/cms/user/mseidel/Summer11_TTJets1725_res:down/analyzeTop.root", fMethod, fBins, 0));
    analyses.push_back(new Analysis("1725_resup", "/scratch/hh/current/cms/user/mseidel/Summer11_TTJets1725_res:up/analyzeTop.root", fMethod, fBins, 0));
    //*/
    
    /*  Uncorrelated systematic samples
    analyses.push_back(new Analysis("1725_scaledown", "/scratch/hh/current/cms/user/mseidel/Summer11_TTJets1725_scaledown/analyzeTop.root", fMethod, fBins, 0));
    analyses.push_back(new Analysis("1725_scaleup", "/scratch/hh/current/cms/user/mseidel/Summer11_TTJets1725_scaleup/analyzeTop.root", fMethod, fBins, 0));
    analyses.push_back(new Analysis("1725_matchingdown", "/scratch/hh/current/cms/user/mseidel/Summer11_TTJets1725_matchingdown/analyzeTop.root", fMethod, fBins, 0));
    analyses.push_back(new Analysis("1725_matchingup", "/scratch/hh/current/cms/user/mseidel/Summer11_TTJets1725_matchingup/analyzeTop.root", fMethod, fBins, 0));
    //*/
    
    /*  Spring11
    analyses.push_back(new Analysis("Spring11_D6T", "/scratch/hh/current/cms/user/mseidel/Spring11_TTJets1725_D6T/analyzeTop.root", fMethod, fBins, 0));
    analyses.push_back(new Analysis("Spring11_Z2", "/scratch/hh/current/cms/user/mseidel/Spring11_TTJets1725_Z2/analyzeTop.root", fMethod, fBins, 0));
    //*/
    
    for (unsigned int i = 0; i < analyses.size(); i++) {
      analyses[i]->Analyze();
    }

    /* DATA, full 2011
    Analysis* aRun2011 = new Analysis("Run2011", "/scratch/hh/current/cms/user/mseidel/Run2011/analyzeTop.root", fMethod, fBins, 0);
    Measure(aRun2011);
    //*/
  
  }
}


void TopMass::WriteEnsembleTest() {
  time_t start, end;
  time(&start);
  time(&end);
  
  //if (readCalibration) LoadXML();
  
  int nEnsembles = po::GetOption<int>("number");
  int n = 0;
  
  double genMass, genJES, genfSig;

  std::map<TString, double> values;

  bool treeCreated = false;
  TTree* tree = 0;
  TFile* ensembleFile = new TFile("ensemble.root", "UPDATE");
  ensembleFile->GetObject("tree", tree);
  
  Analysis* analysis = new Analysis();
  
  while (difftime(end, start) < po::GetOption<double>("walltime") * 60 && n < nEnsembles) {
    std::cout << "\n- - - - - - - - - - " << n << " - - - - - - - - - -\n" << std::endl;
    
    analysis->Analyze();

    const std::map<TString, TH2F*> histograms = analysis->GetH2s();

    if(!treeCreated){

      for(std::map<TString, TH2F*>::const_iterator hist = histograms.begin(); hist != histograms.end(); ++hist){
        values[hist->first] = -1;
      }

      if (!tree) {
        ensembleFile->cd();
        tree = new TTree("tree", "tree");
        tree->Branch("genMass", &genMass, "genMass/D");
        tree->Branch("genJES" , &genJES , "genJES/D" );
        tree->Branch("genfSig", &genfSig, "genfSig/D");

        for(std::map<TString, TH2F*>::const_iterator hist = histograms.begin(); hist != histograms.end(); ++hist){
          TString leafName = hist->first +TString("/D");
          tree->Branch(hist->first, &values[hist->first], leafName);
        }
      }
      else {
        tree->SetBranchAddress("genMass", &genMass);
        tree->SetBranchAddress("genJES" , &genJES );
        tree->SetBranchAddress("genfSig", &genfSig);

        for(std::map<TString, TH2F*>::const_iterator hist = histograms.begin(); hist != histograms.end(); ++hist){
          tree->SetBranchAddress(hist->first, &values[hist->first]);
        }
      }
      treeCreated = true;
      std::cout << "branching finished" << std::endl;
    }

    for (int i = 0; i < fBins_; i++) {
      for (int j = 0; j < fBins_; j++) {
        genMass   = po::GetOption<double>("mass");
        genJES    = po::GetOption<double>("jes" );
        genfSig   = po::GetOption<double>("fsig");

        for(std::map<TString, TH2F*>::const_iterator hist = histograms.begin(); hist != histograms.end(); ++hist){
          if(hist->first.Contains("Pull")){
            TString histName = hist->first; histName.ReplaceAll("_Pull","");
            TString histNameError = histName + TString("_Error");
            double val = histograms.find(histName     )->second->GetCellContent(i+1, j+1);
            double err = histograms.find(histNameError)->second->GetCellContent(i+1, j+1);
            double gen = 0;
            if     (histName.BeginsWith("mass")) gen = genMass;
            else if(histName.BeginsWith("JES" )) gen = genJES;
            else if(histName.BeginsWith("fSig")) gen = genfSig;
            values[hist->first] = (val - gen)/err;
          }
          else{
            values[hist->first] = hist->second->GetCellContent(i+1, j+1);
          }
        }
        tree->Fill();
      }
    }
    
    std::cout << "fill finished" << std::endl;
    
    ++n;
    time(&end);
  }
  delete analysis;

  tree->GetCurrentFile()->Write();
  ensembleFile->Close();
}

void TopMass::QuickCalibration() {

  std::cerr << "Please don't use TopMass::QuickCalibration, as it is not usable currently!!!" << std::endl;
  return;

  enum styles          { kDown, kNominal, kUp};
  int color_   [ 3 ] = { kRed+1, kBlue+1, kGreen+1};
  int marker_  [ 3 ] = { 23, 20, 22};

  Helper* helper = new Helper(fBins_);
  helper->SetTDRStyle();  

  std::vector< std::vector<TH2F*> > hMass;
  std::vector< std::vector<TH2F*> > hMassError;

  std::vector< std::vector<TH2F*> > hJES;
  std::vector< std::vector<TH2F*> > hJESError;

  for (int i = 0; i < 3; i++) {
    calibrationAnalyses.push_back(std::vector<Analysis*>());
    hMass              .push_back(std::vector<TH2F*>());
    hMassError         .push_back(std::vector<TH2F*>());
    hJES               .push_back(std::vector<TH2F*>());
    hJESError          .push_back(std::vector<TH2F*>());
  }

  //TString path = "/scratch/hh/current/cms/user/eschliec/TopMass/19/";
  TString path = "/scratch/hh/dust/naf/cms/user/eschliec/TopMass/19/";

  Analysis* a1615_096 = new Analysis("Z2_S11_ABS_JES_Down_161_5_sig", path+TString("Z2_S11_ABS_JES_Down_161_5_sig.root"), fMethod_, fBins_, 0.0);
  Analysis* a1635_096 = new Analysis("Z2_S11_ABS_JES_Down_163_5_sig", path+TString("Z2_S11_ABS_JES_Down_163_5_sig.root"), fMethod_, fBins_, 0.0);
  Analysis* a1665_096 = new Analysis("Z2_S11_ABS_JES_Down_166_5_sig", path+TString("Z2_S11_ABS_JES_Down_166_5_sig.root"), fMethod_, fBins_, 0.0);
  Analysis* a1695_096 = new Analysis("Z2_S11_ABS_JES_Down_169_5_sig", path+TString("Z2_S11_ABS_JES_Down_169_5_sig.root"), fMethod_, fBins_, 0.0);
  Analysis* a1725_096 = new Analysis("Z2_S11_ABS_JES_Down_172_5_sig", path+TString("Z2_S11_ABS_JES_Down_172_5_sig.root"), fMethod_, fBins_, 0.0);
  Analysis* a1755_096 = new Analysis("Z2_S11_ABS_JES_Down_175_5_sig", path+TString("Z2_S11_ABS_JES_Down_175_5_sig.root"), fMethod_, fBins_, 0.0);
  Analysis* a1785_096 = new Analysis("Z2_S11_ABS_JES_Down_178_5_sig", path+TString("Z2_S11_ABS_JES_Down_178_5_sig.root"), fMethod_, fBins_, 0.0);
  Analysis* a1815_096 = new Analysis("Z2_S11_ABS_JES_Down_181_5_sig", path+TString("Z2_S11_ABS_JES_Down_181_5_sig.root"), fMethod_, fBins_, 0.0);
  Analysis* a1845_096 = new Analysis("Z2_S11_ABS_JES_Down_184_5_sig", path+TString("Z2_S11_ABS_JES_Down_184_5_sig.root"), fMethod_, fBins_, 0.0);

  Analysis* a1615_100 = new Analysis("Z2_S11_161_5_sig", path+TString("Z2_S11_161_5_sig.root"), fMethod_, fBins_, 0.0);
  Analysis* a1635_100 = new Analysis("Z2_S11_163_5_sig", path+TString("Z2_S11_163_5_sig.root"), fMethod_, fBins_, 0.0);
  Analysis* a1665_100 = new Analysis("Z2_S11_166_5_sig", path+TString("Z2_S11_166_5_sig.root"), fMethod_, fBins_, 0.0);
  Analysis* a1695_100 = new Analysis("Z2_S11_169_5_sig", path+TString("Z2_S11_169_5_sig.root"), fMethod_, fBins_, 0.0);
  Analysis* a1725_100 = new Analysis("Z2_S11_172_5_sig", path+TString("Z2_S11_172_5_sig.root"), fMethod_, fBins_, 0.0);
  Analysis* a1755_100 = new Analysis("Z2_S11_175_5_sig", path+TString("Z2_S11_175_5_sig.root"), fMethod_, fBins_, 0.0);
  Analysis* a1785_100 = new Analysis("Z2_S11_178_5_sig", path+TString("Z2_S11_178_5_sig.root"), fMethod_, fBins_, 0.0);
  Analysis* a1815_100 = new Analysis("Z2_S11_181_5_sig", path+TString("Z2_S11_181_5_sig.root"), fMethod_, fBins_, 0.0);
  Analysis* a1845_100 = new Analysis("Z2_S11_184_5_sig", path+TString("Z2_S11_184_5_sig.root"), fMethod_, fBins_, 0.0);

  Analysis* a1615_104 = new Analysis("Z2_S11_ABS_JES_Up_161_5_sig", path+TString("Z2_S11_ABS_JES_Up_161_5_sig.root"), fMethod_, fBins_, 0.0);
  Analysis* a1635_104 = new Analysis("Z2_S11_ABS_JES_Up_163_5_sig", path+TString("Z2_S11_ABS_JES_Up_163_5_sig.root"), fMethod_, fBins_, 0.0);
  Analysis* a1665_104 = new Analysis("Z2_S11_ABS_JES_Up_166_5_sig", path+TString("Z2_S11_ABS_JES_Up_166_5_sig.root"), fMethod_, fBins_, 0.0);
  Analysis* a1695_104 = new Analysis("Z2_S11_ABS_JES_Up_169_5_sig", path+TString("Z2_S11_ABS_JES_Up_169_5_sig.root"), fMethod_, fBins_, 0.0);
  Analysis* a1725_104 = new Analysis("Z2_S11_ABS_JES_Up_172_5_sig", path+TString("Z2_S11_ABS_JES_Up_172_5_sig.root"), fMethod_, fBins_, 0.0);
  Analysis* a1755_104 = new Analysis("Z2_S11_ABS_JES_Up_175_5_sig", path+TString("Z2_S11_ABS_JES_Up_175_5_sig.root"), fMethod_, fBins_, 0.0);
  Analysis* a1785_104 = new Analysis("Z2_S11_ABS_JES_Up_178_5_sig", path+TString("Z2_S11_ABS_JES_Up_178_5_sig.root"), fMethod_, fBins_, 0.0);
  Analysis* a1815_104 = new Analysis("Z2_S11_ABS_JES_Up_181_5_sig", path+TString("Z2_S11_ABS_JES_Up_181_5_sig.root"), fMethod_, fBins_, 0.0);
  Analysis* a1845_104 = new Analysis("Z2_S11_ABS_JES_Up_184_5_sig", path+TString("Z2_S11_ABS_JES_Up_184_5_sig.root"), fMethod_, fBins_, 0.0);

  calibrationAnalyses[0].push_back(a1615_096);
  calibrationAnalyses[0].push_back(a1635_096);
  calibrationAnalyses[0].push_back(a1665_096);
  calibrationAnalyses[0].push_back(a1695_096);
  calibrationAnalyses[0].push_back(a1725_096);
  calibrationAnalyses[0].push_back(a1755_096);
  calibrationAnalyses[0].push_back(a1785_096);
  calibrationAnalyses[0].push_back(a1815_096);
  calibrationAnalyses[0].push_back(a1845_096);

  calibrationAnalyses[1].push_back(a1615_100);
  calibrationAnalyses[1].push_back(a1635_100);
  calibrationAnalyses[1].push_back(a1665_100);
  calibrationAnalyses[1].push_back(a1695_100);
  calibrationAnalyses[1].push_back(a1725_100);
  calibrationAnalyses[1].push_back(a1755_100);
  calibrationAnalyses[1].push_back(a1785_100);
  calibrationAnalyses[1].push_back(a1815_100);
  calibrationAnalyses[1].push_back(a1845_100);

  calibrationAnalyses[2].push_back(a1615_104);
  calibrationAnalyses[2].push_back(a1635_104);
  calibrationAnalyses[2].push_back(a1665_104);
  calibrationAnalyses[2].push_back(a1695_104);
  calibrationAnalyses[2].push_back(a1725_104);
  calibrationAnalyses[2].push_back(a1755_104);
  calibrationAnalyses[2].push_back(a1785_104);
  calibrationAnalyses[2].push_back(a1815_104);
  calibrationAnalyses[2].push_back(a1845_104);

  std::cout << "vectors filled" << std::endl;

  TCanvas* canvasFit = new TCanvas("canvasFit", "hadronic top h2Mass", 500, 500);

  const unsigned short numberOfMasses = 9;
  double genMass[]      = {161.5, 163.5, 166.5, 169.5, 172.5, 175.5, 178.5, 181.5, 184.5};
  double genMassError[] = {1e-6, 1e-6, 1e-6, 1e-6, 1e-6, 1e-6, 1e-6, 1e-6, 1e-6};
  //double genMass[]      = {172.5, 175.5, 178.5, 181.5};
  //double genMassError[] = {1e-6, 1e-6, 1e-6, 1e-6};

  double genJES[]       = {0.96, 1.00, 1.04};
  //double genJESError[]  = {1e-6, 1e-6, 1e-6};

  double hadTopMass[3][numberOfMasses];
  double hadTopMassBias[3][numberOfMasses];
  double hadTopMassError[3][numberOfMasses];

  //double JES[3][numberOfMasses];
  double JESBias[3][numberOfMasses];
  double JESError[3][numberOfMasses];

  //*
  for(int iJES = 0; iJES < 3; iJES++) {
    for(int iMass = 0; iMass < numberOfMasses; iMass++) {
      calibrationAnalyses[iJES][iMass]->Analyze();

      hMass[iJES].push_back(calibrationAnalyses[iJES][iMass]->GetH2("mass"));
      hMassError[iJES].push_back(calibrationAnalyses[iJES][iMass]->GetH2("massError"));

      hJES[iJES].push_back(calibrationAnalyses[iJES][iMass]->GetH2("JES"));
      hJESError[iJES].push_back(calibrationAnalyses[iJES][iMass]->GetH2("JESError"));
    }
  }
  //*/

  /* TEST
  calibrationAnalyses[1][2]->Analyze();
  for(int iJES = 0; iJES < 3; iJES++) {
    for(int iMass = 0; iMass < numberOfMasses; iMass++) {
      hMass[iJES].push_back(calibrationAnalyses[1][2]->GetH2Mass());
      hMassError[iJES].push_back(calibrationAnalyses[1][2]->GetH2MassError());

      hJES[iJES].push_back(calibrationAnalyses[1][2]->GetH2JES());
      hJESError[iJES].push_back(calibrationAnalyses[1][2]->GetH2JESError());
    }
  }
  //*/

  canvasFit->cd();

  std::vector<TGraphErrors*> gMass;
  std::vector<TGraphErrors*> gJES;

  gStyle->SetOptFit(1);

  for (int i = 0; i < fBins_; i++) {
    for (int j = 0; j < fBins_; j++) {
      //*
      if (   hMass     [0][0]->GetCellContent(i+1, j+1) > 0 && hMass     [2][2]->GetCellContent(i+1, j+1) > 0
          && hMassError[0][0]->GetCellContent(i+1, j+1) > 0 && hMassError[2][2]->GetCellContent(i+1, j+1) > 0) {
        for (int iJES = 0; iJES < 3; iJES++) {
          for (int iMass = 0; iMass < numberOfMasses; iMass++) {
            hadTopMass[iJES][iMass]      = genMass[iMass];
            hadTopMassBias[iJES][iMass]  = hMass[iJES][iMass]->GetCellContent(i+1, j+1) - genMass[iMass];
            hadTopMassError[iJES][iMass] = hMassError[iJES][iMass]->GetCellContent(i+1, j+1);

            //JES[iJES][iMass]      = genJES[iJES];
            JESBias[iJES][iMass]  = hJES[iJES][iMass]->GetCellContent(i+1, j+1) - genJES[iJES];
            JESError[iJES][iMass] = hJESError[iJES][iMass]->GetCellContent(i+1, j+1);
          }
        }
        //*/

        TMultiGraph *mgJES = new TMultiGraph();
        mgJES->SetTitle(";m_{t,gen};JES_{meas}-JES_{gen}");

        for (int iJES = 0; iJES < 3; iJES++) {
          double genMassMod[numberOfMasses];
          for (int i = 0; i < numberOfMasses; i++) {
            genMassMod[i] = genMass[i] + 0.2*(iJES-1);
          }

          gJES.push_back(new TGraphErrors(numberOfMasses, genMassMod, JESBias[iJES], genMassError, JESError[iJES]));

          gJES[iJES]->SetMarkerStyle(marker_[iJES]);
          gJES[iJES]->SetMarkerColor(color_ [iJES]);
          gJES[iJES]->SetLineColor  (color_ [iJES]);

          TString sFit("[0]+(x-172.5-"); sFit += 0.2*(iJES-1); sFit += ")*[1]";

          TF1* linearFit = new TF1("linearFit", sFit);
          linearFit->SetParLimits(1, -1, 1);
          linearFit->SetParNames("offset", "slope");
          linearFit->SetLineColor(color_[iJES]);
          gJES[iJES]->Fit("linearFit", "EM");

          mgJES->Add(gJES[iJES]);

          /*
	    for (int l = 0; l < 2; l++) {
	    fCalibFitParameter[i][j][l] = linearFit->GetParameter(l);
	    fCalibFitParError[i][j][l]  = linearFit->GetParError(l);
	    }
           */
        }

        mgJES->SetMinimum(-0.05);
        mgJES->SetMaximum( 0.05);

        mgJES->Draw("AP");

        canvasFit->Update();

        TPaveStats* stats0 = (TPaveStats*) gJES[0]->GetListOfFunctions()->FindObject("stats");
        stats0->SetTextColor(color_[0]);
        stats0->SetX1NDC(0.16);
        stats0->SetY1NDC(0.7);
        stats0->SetX2NDC(0.4233);
        stats0->SetY2NDC(0.825);

        TPaveStats* stats1 = (TPaveStats*) gJES[1]->GetListOfFunctions()->FindObject("stats");
        stats1->SetTextColor(color_[1]);
        stats1->SetX1NDC(0.4233);
        stats1->SetY1NDC(0.7);
        stats1->SetX2NDC(0.6867);
        stats1->SetY2NDC(0.825);

        TPaveStats* statsGlobal = (TPaveStats*) gJES[2]->GetListOfFunctions()->FindObject("stats");

        TPaveStats* stats2 = (TPaveStats*) statsGlobal->Clone();
        stats2->SetTextColor(color_[2]);
        stats2->SetX1NDC(0.6867);
        stats2->SetY1NDC(0.7);
        stats2->SetX2NDC(0.95);
        stats2->SetY2NDC(0.825);
        stats2->Draw();

        canvasFit->Modified();

        TLegend *leg = new TLegend(0.16, 0.825, 0.555, 0.95);
        leg->SetFillColor(kWhite);
        leg->SetBorderSize(1);
        leg->AddEntry( gJES[0], "JES=0.96", "LEP");
        leg->AddEntry( gJES[1], "JES=1.00", "LEP");
        leg->AddEntry( gJES[2], "JES=1.04", "LEP");
        leg->Draw();

        TF1* constFit = new TF1("constFit", "[0]");
        constFit->SetParNames("offset");
        constFit->SetLineColor(kBlack);
        constFit->SetLineWidth(3);
        constFit->SetLineStyle(7);
        mgJES->Fit("constFit", "EM");

        statsGlobal->SetX1NDC(0.555);
        statsGlobal->SetY1NDC(0.825);
        statsGlobal->SetX2NDC(0.95);
        statsGlobal->SetY2NDC(0.95);

        TString path("plot/"); path += fMethod_; path += "/ensembletest/"; path += "fit_JES_"; path += i; path += "_"; path += j; path += ".eps";
        canvasFit->Print(path);

        canvasFit->Clear();

        TMultiGraph *mgMass = new TMultiGraph();
        mgMass->SetTitle("m_{t} bias;m_{t,gen};m_{t,meas}-m_{t,gen}");

        for (int iJES = 0; iJES < 3; iJES++) {
          double genMassMod[numberOfMasses];
          for (int i = 0; i < numberOfMasses; i++) {
            genMassMod[i] = genMass[i] + 0.2*(iJES-1);
          }

          gMass.push_back(new TGraphErrors(numberOfMasses, genMassMod, hadTopMassBias[iJES], genMassError, hadTopMassError[iJES]));

          gMass[iJES]->SetMarkerStyle(marker_[iJES]);
          gMass[iJES]->SetMarkerColor(color_ [iJES]);
          gMass[iJES]->SetLineColor  (color_ [iJES]);

          TString sFit("[0]+(x-172.5-"); sFit += 0.2*(iJES-1); sFit += ")*[1]";

          TF1* linearFit = new TF1("linearFit", sFit);
          linearFit->SetParLimits(1, -1, 1);
          linearFit->SetParNames("offset", "slope");
          linearFit->SetLineColor(color_[iJES]);
          gMass[iJES]->Fit("linearFit", "EM");

          mgMass->Add(gMass[iJES]);

          /*
	    for (int l = 0; l < 2; l++) {
	    fCalibFitParameter[i][j][l] = linearFit->GetParameter(l);
	    fCalibFitParError[i][j][l]  = linearFit->GetParError(l);
	    }
           */
        }

        //mgJES->GetXaxis()->SetLimits( 0.95, 1.05);
        //mgJES->GetYaxis()->SetLimits(-0.05, 0.05);

        mgMass->SetMinimum(-5);
        mgMass->SetMaximum( 5);

        mgMass->Draw("AP");

        canvasFit->Update();
        TPaveStats* statsMass0 = (TPaveStats*) gMass[0]->GetListOfFunctions()->FindObject("stats");
        statsMass0->SetTextColor(color_[0]);
        statsMass0->SetX1NDC(0.16);
        statsMass0->SetY1NDC(0.7);
        statsMass0->SetX2NDC(0.4233);
        statsMass0->SetY2NDC(0.825);
        TPaveStats* statsMass1 = (TPaveStats*) gMass[1]->GetListOfFunctions()->FindObject("stats");
        statsMass1->SetTextColor(color_[1]);
        statsMass1->SetX1NDC(0.4233);
        statsMass1->SetY1NDC(0.7);
        statsMass1->SetX2NDC(0.6867);
        statsMass1->SetY2NDC(0.825);
        TPaveStats* statsMassGlobal = (TPaveStats*) gMass[2]->GetListOfFunctions()->FindObject("stats");
        TPaveStats* statsMass2 = (TPaveStats*) statsMassGlobal->Clone();
        statsMass2->SetTextColor(color_[2]);
        statsMass2->SetX1NDC(0.6867);
        statsMass2->SetY1NDC(0.7);
        statsMass2->SetX2NDC(0.95);
        statsMass2->SetY2NDC(0.825);
        statsMass2->Draw();
        canvasFit->Modified();

        TLegend *legMass = new TLegend(0.16, 0.825, 0.555, 0.95);
        legMass->SetFillColor(kWhite);
        legMass->SetBorderSize(1);
        legMass->AddEntry( gMass[0], "JES=0.96", "LEP");
        legMass->AddEntry( gMass[1], "JES=1.00", "LEP");
        legMass->AddEntry( gMass[2], "JES=1.04", "LEP");
        legMass->Draw();

        mgMass->Fit("constFit", "EM");

        statsMassGlobal->SetX1NDC(0.555);
        statsMassGlobal->SetY1NDC(0.825);
        statsMassGlobal->SetX2NDC(0.95);
        statsMassGlobal->SetY2NDC(0.95);

        path = "plot/"; path += fMethod_; path += "/ensembletest/"; path += "fit_Mass_"; path += i; path += "_"; path += j; path += ".eps";
        canvasFit->Print(path);
      }
    }
  }

  /*
  Measure(a1665);
  Measure(a1725);
  Measure(a1785);
   */
}


void TopMass::QuickSystematics() {

  std::cerr << "Please don't use TopMass::QuickSystematics, as it is not usable currently!!!" << std::endl;
  return;

  Analysis* a1725 = new Analysis("1725", "/scratch/hh/current/cms/user/mseidel/Summer11_TTJets1725_1.00/analyzeTop.root", fMethod_, fBins_, 100000);
  Analysis* a1725_jes_up = new Analysis("1725_jes_up", "/scratch/hh/current/cms/user/mseidel/Summer11_TTJets1725-S4_flavor:up/analyzeTop.root", fMethod_, fBins_, 100000);
  Analysis* a1725_jes_down = new Analysis("1725_jes_down", "/scratch/hh/current/cms/user/mseidel/Summer11_TTJets1725-S4_flavor:down/analyzeTop.root", fMethod_, fBins_, 100000);
  
  a1725_jes_down->Analyze();
  a1725         ->Analyze();
  a1725_jes_up  ->Analyze();
  
  TH2F* hMassJESdown = a1725_jes_down->GetH2("mass");
  TH2F* hMassJESnorm = a1725         ->GetH2("mass");
  TH2F* hMassJESup   = a1725_jes_up  ->GetH2("mass");
  
  hMassJESup->Add(hMassJESnorm, -1.000000);
  hMassJESnorm->Add(hMassJESdown, -1.000000);
  
  hMassJESup->SetTitle("JES up MassError");
  hMassJESnorm->SetTitle("JES down MassError");
  
  TCanvas* canvas = new TCanvas("canvas", "Hadronic top mass JES error", 1000, 500);
  
  canvas->Divide(2,1);
  
  canvas->cd(1);
  hMassJESnorm->Draw("COLZ,TEXT");
  hMassJESnorm->SetAxisRange(0.05, 2, "Z");
  
  //*
  canvas->cd(2);
  hMassJESup->Draw("COLZ,TEXT");
  hMassJESup->SetAxisRange(0.05, 5, "Z");
  //*/
  
  TString path("plot/"); path += fMethod_; path += "_jeserror.eps";
  canvas->Print(path);
}

/*
void TopMass::LoadXML() {
  TString xmlFileName;
  if (fexists("/scratch/hh/lustre/cms/user/mseidel/calibration.xml")) {
    xmlFileName = "/scratch/hh/lustre/cms/user/mseidel/calibration.xml";
  }
  else xmlFileName = "calibration.xml";

  xml::XMLDocument doc(xmlFileName);
  //bool loadOkay = doc.LoadFile();
  
  xml::XMLElement *pRoot, *pParm;
  
  pRoot = doc.FirstChildElement( "calibration" );
  
  pParm = pRoot->FirstChildElement("bin");
  while ( pParm ) {
    int i = -1, j = -1;
    double p0, p0error, p1, p1error;
    pParm->QueryIntAttribute("binx", &i);
    pParm->QueryIntAttribute("biny", &j);
    pParm->QueryDoubleAttribute("p0", &p0);
    pParm->QueryDoubleAttribute("p0error", &p0error);
    pParm->QueryDoubleAttribute("p1", &p1);
    pParm->QueryDoubleAttribute("p1error", &p1error);
    
    fCalibFitParameter[i][j][0] = p0;
    fCalibFitParameter[i][j][1] = p1;
    fCalibFitParError[i][j][0]  = p0error;
    fCalibFitParError[i][j][1]  = p1error;
    
    pParm = pParm->NextSiblingElement( "bin" );
  }
}
*/

bool TopMass::fexists(const char *filename) {
  ifstream ifile(filename);
  return ifile;
}
