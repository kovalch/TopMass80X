#include "TopMass.h"

TopMass::TopMass(po::variables_map vm) {
  fMethod = vm["method"].as<std::string>();
  fBins   = vm["bins"].as<int>();
  fLumi   = vm["lumi"].as<double>();
  
  //QuickCalibration();
  //LoadXML();
  //QuickSystematics();
  
  if (!vm["task"].as<std::string>().compare("cal")) {
    QuickCalibration(vm);
  }
  
  else if (!vm["task"].as<std::string>().compare("pe")) {
    WriteEnsembleTest(vm);
  }
  
  else if (!vm["task"].as<std::string>().compare("sm")) {
    Analysis* analysis = new Analysis(vm);
    analysis->Analyze(vm);
  }
  
  else if (!vm["task"].as<std::string>().compare("hc")) {
  
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
    
    for (int i = 0; i < analyses.size(); i++) {
      analyses[i]->Analyze(vm);
    }

    /* DATA, full 2011
    Analysis* aRun2011 = new Analysis("Run2011", "/scratch/hh/current/cms/user/mseidel/Run2011/analyzeTop.root", fMethod, fBins, 0);
    Measure(aRun2011);
    //*/
  
  }
}


void TopMass::WriteEnsembleTest(po::variables_map vm) {
  time_t start, end;
  time(&start);
  time(&end);
  
  //if (readCalibration) LoadXML();
  
  int nEnsembles = vm["number"].as<int>();
  int n = 0;
  
  double genMass, mass, massError, massPull, genJES, JES, JESError, JESPull;
  TFile* ensembleFile;
  TTree* tree;
  
  ensembleFile = new TFile("ensemble.root", "UPDATE");
  ensembleFile->GetObject("tree", tree);
  
  if (!tree) {
    tree = new TTree("tree", "tree");
    tree->Branch("genMass", &genMass, "genMass/D");
    tree->Branch("mass", &mass, "mass/D");
    tree->Branch("massError", &massError, "massError/D");
    tree->Branch("massPull", &massPull, "massPull/D");
    tree->Branch("genJES", &genJES, "genJES/D");
    tree->Branch("JES", &JES, "JES/D");
    tree->Branch("JESError", &JESError, "JESError/D");
    tree->Branch("JESPull", &JESPull, "JESPull/D");
  }
  else {
    tree->SetBranchAddress("genMass", &genMass);
    tree->SetBranchAddress("mass", &mass);
    tree->SetBranchAddress("massError", &massError);
    tree->SetBranchAddress("massPull", &massPull);
    tree->SetBranchAddress("genJES", &genJES);
    tree->SetBranchAddress("JES", &JES);
    tree->SetBranchAddress("JESError", &JESError);
    tree->SetBranchAddress("JESPull", &JESPull);
  }
  
  Analysis* analysis = new Analysis(vm);
  
  while (difftime(end, start) < vm["walltime"].as<int>() * 60 && n < nEnsembles) {
    std::cout << "\n- - - - - - - - - - " << n << " - - - - - - - - - -\n" << std::endl;
    
    analysis->Analyze(vm);
    
    std::cout << "branching finished" << std::endl;
    
    for (int i = 0; i < fBins; i++) {
      for (int j = 0; j < fBins; j++) {
        genMass   = vm["mass"].as<double>();
        mass      = analysis->GetH2Mass()->GetCellContent(i+1, j+1);
        massError = analysis->GetH2MassError()->GetCellContent(i+1, j+1);
        massPull  = (mass - genMass)/massError;
        
        genJES    = vm["jes" ].as<double>();
        JES       = analysis->GetH2JES()->GetCellContent(i+1, j+1);
        JESError  = analysis->GetH2JESError()->GetCellContent(i+1, j+1);
        JESPull   = (JES - genJES)/JESError;
        
        tree->Fill();
      }
    }
    
    std::cout << "fill finished" << std::endl;
    
    n++;
    time(&end);
  }
  tree->GetCurrentFile()->Write();
  ensembleFile->Close();
}

void TopMass::QuickCalibration(po::variables_map vm) {
	enum styles          { kDown, kNominal, kUp};
	int color_   [ 3 ] = { kRed+1, kBlue+1, kGreen+1};
	int marker_  [ 3 ] = { 23, 20, 22};
	
  Helper* helper = new Helper(fBins);
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
  
	
  Analysis* a1615_096 = new Analysis("1615_096", "/scratch/hh/current/cms/user/mseidel/Summer11_TTJets1615_0.96/analyzeTop.root", fMethod, fBins, 0);
	Analysis* a1635_096 = new Analysis("1635_096", "/scratch/hh/current/cms/user/mseidel/Summer11_TTJets1635_0.96/analyzeTop.root", fMethod, fBins, 0);
	Analysis* a1665_096 = new Analysis("1665_096", "/scratch/hh/current/cms/user/mseidel/Summer11_TTJets1665_0.96/analyzeTop.root", fMethod, fBins, 0);
	Analysis* a1695_096 = new Analysis("1695_096", "/scratch/hh/current/cms/user/mseidel/Summer11_TTJets1695_0.96/analyzeTop.root", fMethod, fBins, 0);
  Analysis* a1725_096 = new Analysis("1725_096", "/scratch/hh/current/cms/user/mseidel/Summer11_TTJets1725_0.96/analyzeTop.root", fMethod, fBins, 0);
	Analysis* a1755_096 = new Analysis("1755_096", "/scratch/hh/current/cms/user/mseidel/Summer11_TTJets1755_0.96/analyzeTop.root", fMethod, fBins, 0);
  Analysis* a1785_096 = new Analysis("1785_096", "/scratch/hh/current/cms/user/mseidel/Summer11_TTJets1785_0.96/analyzeTop.root", fMethod, fBins, 0);
	Analysis* a1815_096 = new Analysis("1815_096", "/scratch/hh/current/cms/user/mseidel/Summer11_TTJets1815_0.96/analyzeTop.root", fMethod, fBins, 0);
	Analysis* a1845_096 = new Analysis("1845_096", "/scratch/hh/current/cms/user/mseidel/Summer11_TTJets1845_0.96/analyzeTop.root", fMethod, fBins, 0);
  
  Analysis* a1615_100 = new Analysis("1615_100", "/scratch/hh/current/cms/user/mseidel/Summer11_TTJets1615_1.00/analyzeTop.root", fMethod, fBins, 0);
	Analysis* a1635_100 = new Analysis("1635_100", "/scratch/hh/current/cms/user/mseidel/Summer11_TTJets1635_1.00/analyzeTop.root", fMethod, fBins, 0);
	Analysis* a1665_100 = new Analysis("1665_100", "/scratch/hh/current/cms/user/mseidel/Summer11_TTJets1665_1.00/analyzeTop.root", fMethod, fBins, 0);
	Analysis* a1695_100 = new Analysis("1695_100", "/scratch/hh/current/cms/user/mseidel/Summer11_TTJets1695_1.00/analyzeTop.root", fMethod, fBins, 0);
  Analysis* a1725_100 = new Analysis("1725_100", "/scratch/hh/current/cms/user/mseidel/Summer11_TTJets1725_1.00/analyzeTop.root", fMethod, fBins, 0);
	Analysis* a1755_100 = new Analysis("1755_100", "/scratch/hh/current/cms/user/mseidel/Summer11_TTJets1755_1.00/analyzeTop.root", fMethod, fBins, 0);
  Analysis* a1785_100 = new Analysis("1785_100", "/scratch/hh/current/cms/user/mseidel/Summer11_TTJets1785_1.00/analyzeTop.root", fMethod, fBins, 0);
	Analysis* a1815_100 = new Analysis("1815_100", "/scratch/hh/current/cms/user/mseidel/Summer11_TTJets1815_1.00/analyzeTop.root", fMethod, fBins, 0);
	Analysis* a1845_100 = new Analysis("1845_100", "/scratch/hh/current/cms/user/mseidel/Summer11_TTJets1845_1.00/analyzeTop.root", fMethod, fBins, 0);
  
  Analysis* a1615_104 = new Analysis("1615_104", "/scratch/hh/current/cms/user/mseidel/Summer11_TTJets1615_1.04/analyzeTop.root", fMethod, fBins, 0);
	Analysis* a1635_104 = new Analysis("1635_104", "/scratch/hh/current/cms/user/mseidel/Summer11_TTJets1635_1.04/analyzeTop.root", fMethod, fBins, 0);
	Analysis* a1665_104 = new Analysis("1665_104", "/scratch/hh/current/cms/user/mseidel/Summer11_TTJets1665_1.04/analyzeTop.root", fMethod, fBins, 0);
	Analysis* a1695_104 = new Analysis("1695_104", "/scratch/hh/current/cms/user/mseidel/Summer11_TTJets1695_1.04/analyzeTop.root", fMethod, fBins, 0);
  Analysis* a1725_104 = new Analysis("1725_104", "/scratch/hh/current/cms/user/mseidel/Summer11_TTJets1725_1.04/analyzeTop.root", fMethod, fBins, 0);
	Analysis* a1755_104 = new Analysis("1755_104", "/scratch/hh/current/cms/user/mseidel/Summer11_TTJets1755_1.04/analyzeTop.root", fMethod, fBins, 0);
  Analysis* a1785_104 = new Analysis("1785_104", "/scratch/hh/current/cms/user/mseidel/Summer11_TTJets1785_1.04/analyzeTop.root", fMethod, fBins, 0);
	Analysis* a1815_104 = new Analysis("1815_104", "/scratch/hh/current/cms/user/mseidel/Summer11_TTJets1815_1.04/analyzeTop.root", fMethod, fBins, 0);
	Analysis* a1845_104 = new Analysis("1845_104", "/scratch/hh/current/cms/user/mseidel/Summer11_TTJets1845_1.04/analyzeTop.root", fMethod, fBins, 0);
	  
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
  
  double genMass[]      = {161.5, 163.5, 166.5, 169.5, 172.5, 175.5, 178.5, 181.5, 184.5};
  double genMassError[] = {1e-6, 1e-6, 1e-6, 1e-6, 1e-6, 1e-6, 1e-6, 1e-6, 1e-6};
  
  double genJES[]       = {0.96, 1.00, 1.04};
  double genJESError[]  = {1e-6, 1e-6, 1e-6};
  
  double hadTopMass[3][9];
  double hadTopMassBias[3][9];
  double hadTopMassError[3][9];
  
  double JES[3][9];
  double JESBias[3][9];
  double JESError[3][9];
  
  //*
  for(int iJES = 0; iJES < 3; iJES++) {
    for(int iMass = 0; iMass < 9; iMass++) {
      calibrationAnalyses[iJES][iMass]->Analyze(vm);
			
      hMass[iJES].push_back(calibrationAnalyses[iJES][iMass]->GetH2Mass());
      hMassError[iJES].push_back(calibrationAnalyses[iJES][iMass]->GetH2MassError());
			
			hJES[iJES].push_back(calibrationAnalyses[iJES][iMass]->GetH2JES());
      hJESError[iJES].push_back(calibrationAnalyses[iJES][iMass]->GetH2JESError());
    }
  }
  //*/
  
  /* TEST
  calibrationAnalyses[1][2]->Analyze();
  for(int iJES = 0; iJES < 3; iJES++) {
    for(int iMass = 0; iMass < 9; iMass++) {
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
  
  for (int i = 0; i < fBins; i++) {
    for (int j = 0; j < fBins; j++) {
      //*
      if (   hMass[0][0]->GetCellContent(i+1, j+1) > 0      && hMass[2][2]->GetCellContent(i+1, j+1) > 0
          && hMassError[0][0]->GetCellContent(i+1, j+1) > 0 && hMassError[2][2]->GetCellContent(i+1, j+1) > 0) {
        for (int iJES = 0; iJES < 3; iJES++) {
          for (int iMass = 0; iMass < 9; iMass++) {
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
					double genMassMod[9];
					for (int i = 0; i < 9; i++) {
						genMassMod[i] = genMass[i] + 0.2*(iJES-1);
					}
					
        	gJES.push_back(new TGraphErrors(9, genMassMod, JESBias[iJES], genMassError, JESError[iJES]));
					
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
				
				TString path("plot/"); path += fMethod; path += "/ensembletest/"; path += "fit_JES_"; path += i; path += "_"; path += j; path += ".eps";
        canvasFit->Print(path);
        
        canvasFit->Clear();
        
        TMultiGraph *mgMass = new TMultiGraph();
				mgMass->SetTitle("m_{t} bias;m_{t,gen};m_{t,meas}-m_{t,gen}");
								
				for (int iJES = 0; iJES < 3; iJES++) {
					double genMassMod[9];
					for (int i = 0; i < 9; i++) {
						genMassMod[i] = genMass[i] + 0.2*(iJES-1);
					}
					
        	gMass.push_back(new TGraphErrors(9, genMassMod, hadTopMassBias[iJES], genMassError, hadTopMassError[iJES]));
					
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
				
				path = "plot/"; path += fMethod; path += "/ensembletest/"; path += "fit_Mass_"; path += i; path += "_"; path += j; path += ".eps";
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


void TopMass::QuickSystematics(po::variables_map vm) {

  Analysis* a1725 = new Analysis("1725", "/scratch/hh/current/cms/user/mseidel/Summer11_TTJets1725_1.00/analyzeTop.root", fMethod, fBins, 100000);
  Analysis* a1725_jes_up = new Analysis("1725_jes_up", "/scratch/hh/current/cms/user/mseidel/Summer11_TTJets1725-S4_flavor:up/analyzeTop.root", fMethod, fBins, 100000);
  Analysis* a1725_jes_down = new Analysis("1725_jes_down", "/scratch/hh/current/cms/user/mseidel/Summer11_TTJets1725-S4_flavor:down/analyzeTop.root", fMethod, fBins, 100000);
  
  a1725_jes_down->Analyze(vm);
  a1725         ->Analyze(vm);
  a1725_jes_up  ->Analyze(vm);
  
  TH2F* hMassJESdown = a1725_jes_down->GetH2Mass();
  TH2F* hMassJESnorm = a1725         ->GetH2Mass();
  TH2F* hMassJESup   = a1725_jes_up  ->GetH2Mass();
  
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
  
  TString path("plot/"); path += fMethod; path += "_jeserror.eps";
  canvas->Print(path);
}

void TopMass::LoadXML() {
  TString xmlFileName;
  if (fexists("/scratch/hh/lustre/cms/user/mseidel/calibration.xml")) {
    xmlFileName = "/scratch/hh/lustre/cms/user/mseidel/calibration.xml";
  }
  else xmlFileName = "calibration.xml";

  TiXmlDocument doc(xmlFileName);
  bool loadOkay = doc.LoadFile();
  
  TiXmlElement *pRoot, *pParm;
  
  pRoot = doc.FirstChildElement( "calibration" );
  
  pParm = pRoot->FirstChildElement("bin");
  while ( pParm ) {
    int i, j;
    double p0, p0error, p1, p1error;
    pParm->Attribute("binx", &i);
    pParm->Attribute("biny", &j);
    pParm->Attribute("p0", &p0);
    pParm->Attribute("p0error", &p0error);
    pParm->Attribute("p1", &p1);
    pParm->Attribute("p1error", &p1error);
    
    fCalibFitParameter[i][j][0] = p0;
    fCalibFitParameter[i][j][1] = p1;
    fCalibFitParError[i][j][0]  = p0error;
    fCalibFitParError[i][j][1]  = p1error;
    
    pParm = pParm->NextSiblingElement( "bin" );
  }
}

int main(int ac, char** av)
{
  // Declare the supported options.
  po::options_description desc("Allowed options");
  desc.add_options()
    ("help,h", "produce help message")
    ("method,m", po::value<std::string>()->default_value("Ideogram"),
      "Top mass measurement method\n"
      "  GenMatch: \tGaussian fit of correct permutations (MC only)\n"
      "  Ideogram: \tJoint likelihood fit of mt and JES given data sample"
    )
    ("task,t", po::value<std::string>()->default_value("sm"),
      "Task to be done\n"
      "  sm: \tSingle measurement based on input parameters\n"
      "  pe: \tPerform pseudo-experiments, additional settings necessary\n"
      "  hc: \tPerform any hardcoded tasks in TopMass constructor"
    )
    ("input,i", po::value<std::string>()->default_value("Summer11_TTJets1725_1.00", "1725_1.00"),
      "Identifier of input file to be analyzed")
    ("lepton,l", po::value<std::string>()->default_value("all"),
      "Lepton+jets channel\n"
      "  electron \t\n"
      "  muon \t\n"
      "  all: \telectron+muon"
    )
    ("bins,b", po::value<int>()->default_value(1), "Number of phasespace bins")
    //("weight,w", po::value<int>()->default_value(0), "Weight type -> think, document")
    ("weight,w", po::value<std::string>()->default_value("muWeight*bWeight*PUWeight"),
      "Event weight used in pseudo-experiments")
    ("fsig,f", po::value<double>()->default_value(1.0), "Signal fraction used in pseudo-experiments")
    ("bdisc,B", po::value<double>()->default_value(0.679), "Threshold for b-jets")
    ("mass,M", po::value<double>()->default_value(172.5), "Input top mass for pseudo-experiments")
    ("jes,J", po::value<double>()->default_value(1.0), "Input JES for pseudo-experiments")
    ("lumi,L", po::value<double>()->default_value(4700.0), "Luminosity for each pseudo-experiment")
    ("number,N", po::value<int>()->default_value(10000), "Number of pseudo-experiments per job")
    ("walltime,W", po::value<int>()->default_value(10), "set walltime limit for pseudo-experiments in minutes")
  ;

  po::variables_map vm;
  po::store(po::parse_command_line(ac, av, desc), vm);
  po::notify(vm);    

  if (vm.count("help")) {
      std::cout << desc << "\n";
      return 1;
  }
  
  TopMass* top = new TopMass(vm);
}
