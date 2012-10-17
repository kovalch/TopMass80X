#include "Analysis.h"

#include "TSystem.h"
#include "TEventList.h"

#include "LHAPDF/LHAPDF.h"
#include <string>

Analysis::Analysis(po::variables_map vm) :
  samplePath("/scratch/hh/lustre/cms/user/eschliec/TopMass/19/"),
  //samplePath"/scratch/hh/dust/naf/cms/user/eschliec/TopMass/19/"),
  fIdentifier(vm["input"].as<std::string>()),
  fMethod    (vm["method"].as<std::string>()),
  fBins      (vm["bins"  ].as<int>()),
  fLumi      (vm["lumi"  ].as<double>()),
  fTree(0),
  tTree(0),
  tTreeBkg(0),
  tempFilePath(gSystem->Getenv("TMPDIR")),
  tempFile(0)
{
  fFile += samplePath;
  fFile += fIdentifier;
  fFile += ".root";

  fChain = new TChain("FullHadTreeWriter/tree");
  fChain->Add(fFile);

  CreateHistos();
}

Analysis::Analysis(TString identifier, TString file, TString method, int bins, double lumi) :
  samplePath("/scratch/hh/lustre/cms/user/eschliec/TopMass/19/"),
  //samplePath"/scratch/hh/dust/naf/cms/user/eschliec/TopMass/19/"),
  fIdentifier(identifier), fMethod(method), fBins(bins), fLumi(lumi),
  fTree(0), tTree(0), tTreeBkg(0),
  tempFilePath(gSystem->Getenv("TMPDIR")),
  tempFile(0)
{
  fFile += samplePath;
  fFile += fIdentifier;
  fFile += ".root";

  fChain = new TChain("FullHadTreeWriter/tree");
  fChain->Add(fFile);

  CreateHistos();
}

Analysis::~Analysis()
{
  delete fChain;

  delete fTree;
  delete tTree;
  delete tTreeBkg;
  delete hEntries;
  delete hMass;
  delete hMassError;
  //delete hMassSigma;
  delete hJES;
  delete hJESError;
  delete hMassConstJES;
  delete hMassConstJESError;
  //delete hMassConstJESSigma;
  delete hFSig;
  delete hFSigError;
  delete hMassfSig;
  delete hMassfSigError;
  //delete hMassfSigSigma;
  delete hJESfSig;
  delete hJESfSigError;
  delete bTagEff;
  delete cTagEff;
  delete lTagEff;
  //delete hMassCalibrated;
  //delete hMassErrorCalibrated;
  tempFile->Close();
  delete tempFile;
}

void Analysis::Analyze(po::variables_map vm) {

  std::cout << "Analyze " << fIdentifier << " with method " << fMethod << std::endl;

  CreateRandomSubset();

  MassAnalyzer* fAnalyzer = 0;

  if (!strcmp(fMethod, "GenMatch")) {
    fAnalyzer = new GenMatchAnalyzer(fIdentifier, fTree);
  }
  else if (!strcmp(fMethod, "MVA")) {
    fAnalyzer = new MVAAnalyzer(fIdentifier, fTree);
  }
  else if (!strcmp(fMethod, "Ideogram")) {
    fAnalyzer = new IdeogramAnalyzer(fIdentifier, fTree);
  }
  else if (!strcmp(fMethod, "RooFit")) {
    fAnalyzer = new RooFitTemplateAnalyzer(fIdentifier, fTree);
  }
  else {
    return;
  }

  Helper* helper = new Helper(fBins);
  helper->SetTDRStyle();
  delete helper;

  TCanvas* canvas = new TCanvas("canvas", "Top mass", 900, 600);
  canvas->cd();

  double smearBins = 1.;
  double rangeX = 10.;
  double rangeY = 10.;

  double smearX = smearBins/fBins*rangeX;
  double smearY = smearBins/fBins*rangeY;

  double minEntries = 25;

  /*if (!strcmp(fMethod, "Ideogram")) {
    minEntries = 1500;
  }*/

  TString observableX = "dRbb";
  TString observableY = "dRbb";

  for(int i = 0; i < fBins; i++) {
    for(int j = 0; j < fBins; j++) {
      // calculate cuts
      std::stringstream stream;
      stream << (rangeX/fBins)*(i)-smearX << "<" << observableX << "&"
             << observableX << "<" << (rangeX/fBins)*(i+1)+smearX << " & "
             << rangeY/fBins*(j)-smearY << "<" << observableY << "&"
             << observableY << "<" << rangeY/fBins*(j+1)+smearY;
      TString cuts(stream.str());

      if (!strcmp(fMethod, "GenMatch")) {
        cuts += " & comboType == 1";
      }
      else if (!strcmp(fMethod, "MVA")) {
        cuts += " & mvaDisc > 0";
      }
      else if (!strcmp(fMethod, "Ideogram")) {
//        //cuts += " & target == 1";
//        //cuts += " & (target == 0 | target == 1)";
//        //cuts += " & Njet == 6";
//        cuts += " & probs[0] > 0.09";
//        cuts += " & dRbb > 1.5";
//        //cuts += " & nVertex >= 5";
//        //cuts += " & nVertex <= 5";
//        cuts += " & topMasses[0] > 100. & topMasses[0] < 550.";
      }
      else if (!strcmp(fMethod, "RooFit")) {
      }

      int entries = fTree->GetEntries(cuts);
      std::cout << cuts << std::endl;
      std::cout << entries << std::endl;

      hEntries->SetCellContent(i+1, j+1, entries);

      if (entries > minEntries) {
        fAnalyzer->Analyze(cuts, i, j, vm);

        hMass     ->SetCellContent(i+1, j+1, fAnalyzer->GetMass());
        hMassError->SetCellContent(i+1, j+1, fAnalyzer->GetMassError());
        //hMassSigma->SetCellContent(i+1, j+1, fAnalyzer->GetMassSigma());
        hJES      ->SetCellContent(i+1, j+1, fAnalyzer->GetJES());
        hJESError ->SetCellContent(i+1, j+1, fAnalyzer->GetJESError());
        hMassConstJES     ->SetCellContent(i+1, j+1, fAnalyzer->GetMassConstJES());
        hMassConstJESError->SetCellContent(i+1, j+1, fAnalyzer->GetMassConstJESError());
        //hMassConstJESSigma->SetCellContent(i+1, j+1, fAnalyzer->GetMassConstJESSigma());
        hFSig     ->SetCellContent(i+1, j+1, fAnalyzer->GetFSig());
        hFSigError->SetCellContent(i+1, j+1, fAnalyzer->GetFSigError());
        hMassfSig     ->SetCellContent(i+1, j+1, fAnalyzer->GetMassfSig());
        hMassfSigError->SetCellContent(i+1, j+1, fAnalyzer->GetMassfSigError());
        //hMassfSigSigma->SetCellContent(i+1, j+1, fAnalyzer->GetMassfSigSigma());
        hJESfSig      ->SetCellContent(i+1, j+1, fAnalyzer->GetJESfSig());
        hJESfSigError ->SetCellContent(i+1, j+1, fAnalyzer->GetJESfSigError());

        std::cout << std::endl;
        std::cout << "Measured 1D mass: " << fAnalyzer->GetMassConstJES() << " +/- "
                  << fAnalyzer->GetMassConstJESError() << " GeV" << std::endl;
        std::cout << std::endl;
        std::cout << "Measured 2D mass: " << fAnalyzer->GetMass() << " +/- "
                  << fAnalyzer->GetMassError() << " GeV" << std::endl;
        std::cout << "Measured 2D JES: " << fAnalyzer->GetJES() << " +/- "
                  << fAnalyzer->GetJESError() << " " << std::endl;
        std::cout << std::endl;
        std::cout << "Measured 3D mass: " << fAnalyzer->GetMassfSig() << " +/- "
                  << fAnalyzer->GetMassfSigError() << " GeV" << std::endl;
        std::cout << "Measured 3D JES: " << fAnalyzer->GetJESfSig() << " +/- "
                  << fAnalyzer->GetJESfSigError() << " " << std::endl;
        std::cout << "Measured 3D fSig: " << fAnalyzer->GetFSig() << " +/- "
                  << fAnalyzer->GetFSigError() << " GeV" << std::endl;
        std::cout << std::endl;
      }
    }
  }

  canvas->Clear();
  canvas->Divide(2,2);

  canvas->cd(1);
  hEntries->Draw("COLZ");

  canvas->cd(2);
  hMass->Draw("COLZ,TEXT");
  hMass->SetAxisRange(hMass->GetMinimum(0.05), hMass->GetMaximum(), "Z");

  canvas->cd(3);
  hMassError->Draw("COLZ,TEXT");
  hMassError->SetAxisRange(0.05, 5, "Z");

  //canvas->cd(4);
  //hMassSigma->Draw("COLZ,TEXT");
  //hMassSigma->SetAxisRange(hMassSigma->GetMinimum(0.05), hMassSigma->GetMaximum(), "Z");

  TString path("plot/"); path += fMethod; path += "_"; path += fIdentifier; path += ".eps";
  canvas->Print(path);

  delete canvas;
  delete fAnalyzer;
}

void Analysis::CreateHistos() {
  Helper* helper = new Helper(fBins);

  hEntries = helper->GetH2("Entries");

  hMass = helper->GetH2("Mass");
  hMassError = helper->GetH2("MassError");
  //hMassSigma = helper->GetH2("MassSigma");
  hJES = helper->GetH2("JES");
  hJESError = helper->GetH2("JESError");
  
  hMassConstJES = helper->GetH2("MassConstJES");
  hMassConstJESError = helper->GetH2("MassConstJESError");
  //hMassConstJESSigma = helper->GetH2("MassConstJESSigma");
  
  hFSig = helper->GetH2("fSig");
  hFSigError = helper->GetH2("fSigError");
  hMassfSig = helper->GetH2("MassfSig");
  hMassfSigError = helper->GetH2("MassfSigError");
  //hMassfSigSigma = helper->GetH2("MassfSigSigma");
  hJESfSig = helper->GetH2("JESfSig");
  hJESfSigError = helper->GetH2("JESfSigError");
  
  //hMassCalibrated = helper->GetH2("Mass (Calibrated)");
  //hMassErrorCalibrated = helper->GetH2("MassError (Calibrated)");

  TFile * bTagFile = TFile::Open(samplePath+TString("bTagFile.root"));
  gROOT->cd();
  bTagEff = (TH2F*)bTagFile->Get("histb")->Clone();
  cTagEff = (TH2F*)bTagFile->Get("histc")->Clone();
  lTagEff = (TH2F*)bTagFile->Get("histl")->Clone();
  bTagFile->Close();

  delete bTagFile;
  delete helper;
}

void Analysis::CreateRandomSubset() {
  std::cout << "Create random subset..." << std::endl;
  fChain->SetBranchStatus("*",0);
  fChain->SetBranchStatus("nCombos", 1);
  fChain->SetBranchStatus("comboTypes", 1);
  fChain->SetBranchStatus("topMasses", 1);
  fChain->SetBranchStatus("topMass", 1);
  fChain->SetBranchStatus("w1Mass", 1);
  fChain->SetBranchStatus("w2Mass", 1);
  fChain->SetBranchStatus("probs", 1);
  fChain->SetBranchStatus("prob", 1);
  fChain->SetBranchStatus("dRbb", 1);
  
  fChain->SetBranchStatus("Njet", 1);
  fChain->SetBranchStatus("jets", 1);
  //fChain->SetBranchStatus("bTag_CSV", 1);
  fChain->SetBranchStatus("partonFlavour", 1);

  fChain->SetBranchStatus("runNumber", 1);
  fChain->SetBranchStatus("luminosityBlockNumber", 1);
  fChain->SetBranchStatus("eventNumber", 1);

  fChain->SetBranchStatus("dRbb", 1);

  fChain->SetBranchStatus("nPU", 1);
  fChain->SetBranchStatus("nPUTru", 1);

  fChain->SetBranchStatus("MCweight", 1);
  
  if(fFile.Contains("PDF")){
    fChain->SetBranchStatus("id1", 1);
    fChain->SetBranchStatus("id2", 1);
    fChain->SetBranchStatus("x1", 1);
    fChain->SetBranchStatus("x2", 1);
    fChain->SetBranchStatus("Q", 1);
  }
  const short kMAXCombos = 12000;
  const short kMAXNJets  = 50;

  double CombinedWeight = -1.;
  double meanWMass      = -1.;
  double * topMasses    = new double[kMAXCombos];
  double topMass        = -1.;
  double * w1Masses     = new double[kMAXCombos];
  double * w2Masses     = new double[kMAXCombos];
  double * probs        = new double[kMAXCombos];
  double prob           = -1.;
  int Njet              = -1;
  TClonesArray * jets   = new TClonesArray("TLorentzVector");
  //float * bTag = new float[kMAXNJets];
  float dRbb            = -1.;
  unsigned int nCombos  =  0;
  unsigned short * comboTypes = new unsigned short[kMAXCombos];
  short * pdgId = new short[kMAXNJets];
  unsigned int runNumber             = 0;
  unsigned int luminosityBlockNumber = 0;
  unsigned int eventNumber           = 0;
  double MCweight = 0.;
  short nPU = -1;
  short nPUTru = -1;
  double x1 = -1.,x2 = -1.;
  float Q = -1.;
  int id1 = 0, id2 = 0;  

  TString sel = "probs[0]>0.09 && dRbb>1.5 && topMasses[0] > 100. && topMasses[0] < 550.";
  TString triggerUncertainty = "";
  if(fIdentifier.Contains("jet4_02")) triggerUncertainty += " & jets[3].Pt() > 62.";
  if(fIdentifier.Contains("jet4_05")) triggerUncertainty += " & jets[3].Pt() > 65.";
  if(fIdentifier.Contains("jet4_10")) triggerUncertainty += " & jets[3].Pt() > 70.";
  if(fIdentifier.Contains("jet5_02")) triggerUncertainty += " & jets[4].Pt() > 52.";
  if(fIdentifier.Contains("jet5_05")) triggerUncertainty += " & jets[4].Pt() > 55.";
  if(fIdentifier.Contains("jet5_10")) triggerUncertainty += " & jets[4].Pt() > 60.";
  if(fIdentifier.Contains("jet6_02")) triggerUncertainty += " & jets[5].Pt() > 42.";
  if(fIdentifier.Contains("jet6_05")) triggerUncertainty += " & jets[5].Pt() > 45.";
  if(fIdentifier.Contains("jet6_10")) triggerUncertainty += " & jets[5].Pt() > 50.";

  if (fLumi!=0) {
    tempFile = TFile::Open(tempFilePath+TString("/tempTree.root"), "READ");
    if(tempFile && !tempFile->IsZombie()){
      tTree    = (TTree*)tempFile->Get(TString("fullTree_")+fIdentifier);
      tTreeBkg = (TTree*)tempFile->Get("fullTreeBkg");
    }
    if(!tempFile || tempFile->IsZombie() || !tTree || !tTreeBkg){
      tempFile = TFile::Open(tempFilePath+TString("/tempTree.root"), "RECREATE", "", 0);
      std::cout << "Creating: " << tempFile->GetName() << std::endl;

      if(!tTree){
	fChain->Draw(">>selectedEvents",sel+triggerUncertainty,"goff");
	TEventList* selectedEvents = (TEventList*)gDirectory->Get("selectedEvents");
      
	fChain->SetBranchAddress("nCombos", &nCombos);
	fChain->SetBranchAddress("comboTypes", comboTypes);
	fChain->SetBranchAddress("topMasses", topMasses);
	fChain->SetBranchAddress("topMass", &topMass);
	fChain->SetBranchAddress("w1Mass", w1Masses);
	fChain->SetBranchAddress("w2Mass", w2Masses);
	fChain->SetBranchAddress("probs", probs);
	fChain->SetBranchAddress("prob", &prob);
	fChain->SetBranchAddress("dRbb", &dRbb);
  
	fChain->SetBranchAddress("Njet", &Njet);
	fChain->SetBranchAddress("jets", &jets);
	//fChain->SetBranchAddress("bTag_CSV", bTag);
	fChain->SetBranchAddress("partonFlavour", pdgId);
  
	fChain->SetBranchAddress("runNumber", &runNumber);
	fChain->SetBranchAddress("luminosityBlockNumber", &luminosityBlockNumber);
	fChain->SetBranchAddress("eventNumber", &eventNumber);
  
	fChain->SetBranchAddress("dRbb", &dRbb);
  
	fChain->SetBranchAddress("nPU", &nPU);
	fChain->SetBranchAddress("nPUTru", &nPUTru);
  
	fChain->SetBranchAddress("MCweight", &MCweight);
  
	if(fFile.Contains("PDF")){
	  fChain->SetBranchAddress("id1", &id1);
	  fChain->SetBranchAddress("id2", &id2);
	  fChain->SetBranchAddress("x1", &x1);
	  fChain->SetBranchAddress("x2", &x2);
	  fChain->SetBranchAddress("Q", &Q);
	}

	tTree = fChain->CloneTree(0);
	tTree->SetName(TString("fullTree_")+fIdentifier);
      
	for(int idx = 0 , l = selectedEvents->GetN(); idx < l; ++idx){
	  fChain->GetEntry(selectedEvents->GetEntry(idx));
	  nCombos = 1;
	  tTree->Fill();
	}
	AddWeights(tTree);
	tTree->Write(0,TObject::kOverwrite);
	delete tTree;
      }

      if(!tTreeBkg){
	TFile * fileBkg = TFile::Open(samplePath+TString("QCDEstimationMix_2011_skimmed2.root"));
	if(TString(fileBkg->GetName()).Contains("_skimmed")) tTreeBkg = (TTree*)fileBkg->Get("tree");
	else  tTreeBkg = (TTree*)fileBkg->Get("analyzeFullHadEventMixer/tree");
    
	tTreeBkg->SetBranchStatus("*",0);
	tTreeBkg->SetBranchStatus("topMasses", 1);
	tTreeBkg->SetBranchStatus("topMass", 1);
	tTreeBkg->SetBranchStatus("w1Mass", 1);
	tTreeBkg->SetBranchStatus("w2Mass", 1);
	tTreeBkg->SetBranchStatus("probs", 1);
	tTreeBkg->SetBranchStatus("prob", 1);
	tTreeBkg->SetBranchStatus("dRbb", 1);
    
	tTreeBkg->SetBranchStatus("jets", 1);
    
	tTreeBkg->SetBranchStatus("runNumber", 1);
	tTreeBkg->SetBranchStatus("luminosityBlockNumber", 1);
	tTreeBkg->SetBranchStatus("eventNumber", 1);
    
	tempFile->cd();
    
	TString selBkg = sel; selBkg.ReplaceAll("dRbb","dRbb[0]");
	TTree* tempTreeBkg = tTreeBkg->CopyTree(selBkg);
	tempTreeBkg->SetName("fullTreeBkg");
	AddWeights(tempTreeBkg,true);
	tempTreeBkg->Write(0,TObject::kOverwrite);
	delete tempTreeBkg;

	delete tTreeBkg;
	fileBkg->Close();
	delete fileBkg;
      }

      tempFile->Close();
      tempFile = TFile::Open(tempFilePath+TString("/tempTree.root"),"READ");
      tTree    = (TTree*)tempFile->Get(TString("fullTree_")+fIdentifier);
      tTreeBkg = (TTree*)tempFile->Get("fullTreeBkg");
    }

    TRandom3* myRandom = new TRandom3(0);
    std::cout << "Random seed: " << myRandom->GetSeed() << std::endl;

    tTree   ->SetCacheSize(10000000000);
    tTreeBkg->SetCacheSize(10000000000);

    tTree->SetBranchStatus("*",0);
    tTree->SetBranchStatus("CombinedWeight",1);
    tTree->SetBranchStatus("meanWMass",1);
    tTree->SetBranchStatus("topMasses",1);
    tTree->SetBranchStatus("topMass",1);
    tTree->SetBranchStatus("w1Mass",1);
    tTree->SetBranchStatus("w2Mass",1);
    tTree->SetBranchStatus("probs",1);
    tTree->SetBranchStatus("prob",1);
    tTree->SetBranchStatus("dRbb",1);
    tTree->SetBranchStatus("jets",1);
    tTree->SetBranchStatus("nCombos",1);
    tTree->SetBranchStatus("comboTypes",1);
    tTree->SetBranchStatus("runNumber", 1);
    tTree->SetBranchStatus("luminosityBlockNumber", 1);
    tTree->SetBranchStatus("eventNumber", 1);

    tTree->SetBranchAddress("CombinedWeight",&CombinedWeight);
    tTree->SetBranchAddress("meanWMass",&meanWMass);
    tTree->SetBranchAddress("topMasses",topMasses);
    tTree->SetBranchAddress("topMass",&topMass);
    tTree->SetBranchAddress("w1Mass",w1Masses);
    tTree->SetBranchAddress("w2Mass",w2Masses);
    tTree->SetBranchAddress("probs",probs);
    tTree->SetBranchAddress("prob",&prob);
    tTree->SetBranchAddress("dRbb",&dRbb);
    tTree->GetBranch("jets")->SetAutoDelete(kFALSE);
    tTree->SetBranchAddress("jets", &jets);
    tTree->SetBranchAddress("nCombos",&nCombos);
    tTree->SetBranchAddress("comboTypes",comboTypes);
    tTree->SetBranchAddress("runNumber", &runNumber);
    tTree->SetBranchAddress("luminosityBlockNumber", &luminosityBlockNumber);
    tTree->SetBranchAddress("eventNumber", &eventNumber);

    tTreeBkg->SetBranchAddress("meanWMass",&meanWMass);
    tTreeBkg->SetBranchAddress("topMasses", topMasses);
    tTreeBkg->SetBranchAddress("topMass", &topMass);
    tTreeBkg->SetBranchAddress("w1Mass", w1Masses);
    tTreeBkg->SetBranchAddress("w2Mass", w2Masses);
    tTreeBkg->SetBranchAddress("probs", probs);
    tTreeBkg->SetBranchAddress("prob", &prob);
    tTreeBkg->SetBranchAddress("dRbb", &dRbb);
    tTreeBkg->GetBranch("jets")->SetAutoDelete(kFALSE);
    tTreeBkg->SetBranchAddress("jets", &jets);
    tTreeBkg->SetBranchAddress("nCombos",&nCombos);
    
    tTreeBkg->SetBranchAddress("runNumber", &runNumber);
    tTreeBkg->SetBranchAddress("luminosityBlockNumber", &luminosityBlockNumber);
    tTreeBkg->SetBranchAddress("eventNumber", &eventNumber);
    
    double maxMCWeight = tTree->GetMaximum("CombinedWeight");

    if (maxMCWeight == -1) { std::cout << "Running over data?" << std::endl; }
    
    double fSig = 0.539; // 0.504; // 
    if(fFile.Contains("fSig_Up"))
      fSig += 0.10;
    else if(fFile.Contains("fSig_Down"))
      fSig -= 0.10;
    else if(fFile.Contains("BJES_Up"))
      fSig *= 1.0124;
    else if(fFile.Contains("BJES_Down"))
      fSig *= 0.9895;
    else if(fFile.Contains("JES_Up"))
      fSig *= 1.1042;
    else if(fFile.Contains("JES_Down"))
      fSig *= 0.8967;
    else if(fFile.Contains("JER_Up"))
      fSig *= 0.9816;
    else if(fFile.Contains("JER_Down"))
      fSig *= 1.0192;
    else if(fFile.Contains("BTAG_Up"))
      fSig *= 1.0509;
    else if(fFile.Contains("BTAG_Down"))
      fSig *= 0.9565;
    else if(fFile.Contains("Scale_Up"))
      fSig *= 0.8726;
    else if(fFile.Contains("Scale_Down"))
      fSig *= 1.1184;
    else if(fFile.Contains("Matching_Up"))
      fSig *= 0.9718;
    else if(fFile.Contains("P11_NoCR"))
      fSig *= 1.0294;
    else if(fFile.Contains("P11mpiHi"))
      fSig *= 1.0202;
    else if(fFile.Contains("P11TeV"))
      fSig *= 1.0234;
    else if(fFile.Contains("jet4_02"))
      fSig *= 0.9256;
    else if(fFile.Contains("jet5_02"))
      fSig *= 0.9431;
    else if(fFile.Contains("jet6_02"))
      fSig *= 0.9305;

    int permsMC  = tTree   ->GetEntries();
    int permsBkg = tTreeBkg->GetEntries();
    int eventsPE = myRandom->Poisson(2410./3544.844*fLumi); // add poisson 
    int eventsDrawn = 0;
    
    tempFile->ReOpen("UPDATE");
    fTree = tTree->CloneTree(0);
    fTree->SetName("tree");
    // remove unnecessary branches from fTree
    fTree->SetBranchStatus("CombinedWeight",0);
    fTree->SetBranchStatus("jets",0);
    //fTree->SetBranchStatus("runNumber", 0);
    //fTree->SetBranchStatus("luminosityBlockNumber", 0);
    //fTree->SetBranchStatus("eventNumber", 0);

    int signalDrawn = 0, backgroundDrawn = 0;
    if(fLumi>0) {
      //int drawCounter = 0;
      while (eventsDrawn < eventsPE) {
	double fSigRndm = myRandom->Uniform(0.,1.);
	if(fSigRndm < fSig){
	  int oldEventsDrawn = eventsDrawn;
	  do {
	    int drawn = myRandom->Integer(permsMC);
	    //++drawCounter;
	    tTree->GetEntry(drawn);
	    //std::cout << CombinedWeight << " " << maxMCWeight << std::endl;
	    if (CombinedWeight > myRandom->Uniform(0., maxMCWeight)) {
	      // reduce size of tree before filling
	      nCombos = 1;
	      //if(fTree->Fill() == -1) std::cout << eventsDrawn << " " << meanWMass << " " << topMasses[0] << " " << topMass << " " << w1Masses[0] << " " << w2Masses[0] << " " << probs[0] << " " << prob << " " << dRbb << " " << nCombos << " " << comboTypes[0] << " " << runNumber << " " << luminosityBlockNumber << " " << eventNumber << " " << std::endl;
	      fTree->Fill();

	      ++eventsDrawn;
	      ++signalDrawn;
	    }
	  }
	  while(oldEventsDrawn == eventsDrawn);
	}
	else{
	  int drawn = myRandom->Integer(permsBkg);
	  //++drawCounter;
	  tTreeBkg->GetEntry(drawn);
	  // reduce size of tree before filling
	  nCombos = 1;
	  // set missing variable for background tree
	  comboTypes[0] = 0;
	  //if(fTree->Fill() == -1) std::cout << eventsDrawn << " " << meanWMass << " " << topMasses[0] << " " << topMass << " " << w1Masses[0] << " " << w2Masses[0] << " " << probs[0] << " " << prob << " " << dRbb << " " << nCombos << " " << comboTypes[0] << " " << runNumber << " " << luminosityBlockNumber << " " << eventNumber << " " << std::endl;
	  fTree->Fill();
	  ++eventsDrawn;
	  ++backgroundDrawn;
	}
      }
      //std::cout << "DRAWCOUNTER: " << drawCounter << std::endl;
    }
    else {
      for(int ev = 0; ev < permsMC; ++ev){
	tTree->GetEntry(ev);
	fTree->Fill();
	++signalDrawn;
      }
      std::cout << "wanted / available events: " << permsMC*((1.-fSig)/fSig) << " / " << permsBkg << std::endl;
      for(int ev = 0; ev < permsMC*((1.-fSig)/fSig); ++ev){
	int drawn = myRandom->Integer(permsBkg);
      	tTreeBkg->GetEntry(drawn);
      	fTree->Fill();
	++backgroundDrawn;
      }
    }
    fTree->Write(0,TObject::kOverwrite);
    std::cout << "Events drawn: " << eventsDrawn << " (sig: " << signalDrawn << ", bkg: " << backgroundDrawn << ") -> fSig: " << double(signalDrawn)/double(signalDrawn+backgroundDrawn) << " (default: " << fSig << ")" << std::endl;
    delete fTree;
    delete tTree;
    delete tTreeBkg;
    tempFile->Close();
    tempFile = TFile::Open(tempFilePath+TString("/tempTree.root"));
    fTree = (TTree*)tempFile->Get("tree");
    delete myRandom;
  }
  else {
    tempFile = new TFile(tempFilePath+TString("/tempTree.root"), "RECREATE");
    tempFile->cd();
    tTree = fChain->CopyTree(sel);
    AddWeights(tTree,true);
    tempFile->Write(0,TObject::kOverwrite);
    delete tTree;
    tempFile->Close();
    tempFile = TFile::Open(tempFilePath+TString("/tempTree.root"));
    fTree = (TTree*)tempFile->Get("tree");
  }
  delete[] topMasses;
  delete[] w1Masses;
  delete[] w2Masses;
  delete[] probs;
  delete[] comboTypes;
  delete[] pdgId;
  jets->Delete();
  std::cout << "Created random subset." << std::endl;
}

void Analysis::AddWeights(TTree* tempTree, bool isData) { //, enumForPUWeights whichSample, int whichPDF) {
  std::cout << "Adapting " << tempTree->GetName();
  if(!isData) std::cout << " and adding combined weight";
  std::cout << " ..." << std::endl;

  enumForPUWeights whichSample = kFall10;
  if(!isData){
    if     (fFile.Contains("S11"        )) whichSample = kSummer11;
    else if(fFile.Contains("FAST"       )) whichSample = kFall10;
    else if(fFile.Contains("F11_PU_Up"  )) whichSample = kFall11Plus05;
    else if(fFile.Contains("F11_PU_Down")) whichSample = kFall11Minus05;
    else if(fFile.Contains("F11"        )) whichSample = kFall11;
  }
  int whichPDFUncertainty = -1;
  bool upVariation = true; // value does not matter as long as _whichPDFUncertainty_==-1
  if(!isData){
    if(fFile.Contains("PDF")){
      TString pdfUncert = fFile;
      upVariation = pdfUncert.Contains("Up") ? true : false;
      pdfUncert = TString(pdfUncert(pdfUncert.Index("PDF",TString::kIgnoreCase)+3,2).Data());
      whichPDFUncertainty = pdfUncert.Atoi();
    }
    std::cout << "whichPDFUncertainty: " << whichPDFUncertainty << " (" << fFile << ")" << std::endl;
  }

  const int kMAX_(50);
  TString bTagAlgo_ = "CSV";
  int Njet;
  if(!isData) tempTree->SetBranchStatus ("Njet", 1);
  if(!isData) tempTree->SetBranchAddress("Njet", &Njet);
  //float bTag[kMAX_];
  //if(!isData) tempTree->SetBranchStatus (TString("bTag_")+bTagAlgo_, 1);
  //if(!isData) tempTree->SetBranchAddress(TString("bTag_")+bTagAlgo_, &bTag );
  short pdgId[kMAX_];
  if(!isData) tempTree->SetBranchStatus ("partonFlavour", 1);       // pdgId  *or*  partonFlavour
  if(!isData) tempTree->SetBranchAddress("partonFlavour", &pdgId); //  pdgId  *or*  partonFlavour
  TClonesArray * jets = new TClonesArray("TLorentzVector");
  if(!isData) tempTree->SetBranchStatus("jets", 1);
  if(!isData) tempTree->GetBranch("jets")->SetAutoDelete(kFALSE);
  if(!isData) tempTree->SetBranchAddress("jets", &jets);
  //unsigned int nCombos  = 0;
  //tempTree->SetBranchStatus ("nCombos", 1);
  //tempTree->SetBranchAddress("nCombos", &nCombos);
  double * w1Masses  = new double[12000];
  double * w2Masses  = new double[12000];
  tempTree->SetBranchStatus("w1Mass", 1);
  tempTree->SetBranchStatus("w2Mass", 1);
  tempTree->SetBranchAddress("w1Mass",w1Masses);
  tempTree->SetBranchAddress("w2Mass",w2Masses);

  //tempTree->SetBranchStatus("*",0);
  double meanWMass      = -1.;
  double CombinedWeight = -1.;
  double PUWeight       = -1.;
  double MCWeight       =  1.;
  double PDFWeight      =  1.;
  double BTagWeight     =  1.;
  short nPU = -1;
  double x1,x2;
  float Q;
  int id1,id2;
  TBranch * br = 0;
  if(!isData) br = tempTree->Branch("CombinedWeight", &CombinedWeight, "CombinedWeight/D");
  TBranch * br5 = tempTree->Branch("meanWMass", &meanWMass, "meanWMass/D");
  //TBranch * br1 = tempTree->Branch("PUWeight"      , &PUWeight      , "PUWeight/D");
  ////TBranch * br2 = tempTree->Branch("MCWeight"      , &MCWeight      , "MCWeight/D");
  //TBranch * br3 = tempTree->Branch("PDFWeight"     , &PDFWeight     , "PDFWeight/D");
  //TBranch * br4 = tempTree->Branch("BTagWeight"    , &BTagWeight    , "BTagWeight/D");
  if(!isData){
    if(whichSample == kSummer11 || whichSample == kSummer11Plus05 || whichSample == kSummer11Minus05){
      tempTree->SetBranchStatus("nPUTru", 0);
      tempTree->SetBranchAddress("nPU", &nPU);
    }
    else if(whichSample == kFall11 || whichSample == kFall11Plus05 || whichSample == kFall11Minus05){
      tempTree->SetBranchStatus("nPU", 0);
      tempTree->SetBranchAddress("nPUTru", &nPU);
    }
    if(whichPDFUncertainty >= 0){
      tempTree->SetBranchAddress("id1", &id1);
      tempTree->SetBranchAddress("id2", &id2);
      tempTree->SetBranchAddress("x1", &x1);
      tempTree->SetBranchAddress("x2", &x2);
      tempTree->SetBranchAddress("Q", &Q);
    }
    tempTree->SetBranchAddress("MCweight", &MCWeight);
    //tempTree->SetBranchAddress("PUweight", &PUWeight);
  }
  for(int i = 0, l = tempTree->GetEntries(); i < l; ++i){
    tempTree->GetEntry(i);
    if(!isData) {
      PUWeight = calcPUWeight_(whichSample, nPU);
      //PUWeight = 1.;
      BTagWeight = calcBTagWeight_(Njet, pdgId, jets);
      if(whichPDFUncertainty >= 0) PDFWeight = calcPDFWeight_(whichPDFUncertainty, upVariation, x1, id1, Q, x2, id2);
      CombinedWeight = PUWeight * BTagWeight * MCWeight * PDFWeight;
      //std::cout << CombinedWeight << " ";
      br->Fill();
    }
    meanWMass = (w1Masses[0]+w2Masses[0])/2.0;
    br5->Fill();
    //br1->Fill();
    ////br2->Fill();
    //br3->Fill();
    //br4->Fill();
  }
  jets->Delete();
  delete[] w1Masses;
  delete[] w2Masses;
  std::cout << "Adapted  " << tempTree->GetName();
  if(!isData) std::cout << " and added  combined weight";
  std::cout << " ..." << std::endl;
}

void Analysis::AdaptBkgTree(TTree* tempTreeBkg) {
  std::cout << "ADAPTING BACKGROUND TREE" << std::endl;
  
  double CombinedWeight = 1.;
  //const unsigned int kMAXCombo = 12000;
  //float* dRbbOLD = new float[kMAXCombo];
  //float dRbb;
  TBranch * br1 = tempTreeBkg->Branch("CombinedWeight", &CombinedWeight, "CombinedWeight/D");
  //TBranch * br2 = tempTreeBkg->Branch("dRbbNEW", &dRbb, "dRbb/F");
  //tempTreeBkg->SetBranchAddress("dRbb", dRbbOLD);
  for(int i = 0; i < tempTreeBkg->GetEntries(); ++i){
    tempTreeBkg->GetEntry(i);
    //dRbb = dRbbOLD[0];
    br1->Fill();
    //br2->Fill();
  }
  //tempTreeBkg->SetBranchStatus("dRbb", 0);
  //br2->SetName("dRbb");

  std::cout << "ADAPTED  BACKGROUND TREE" << std::endl;
}

TH2F* Analysis::GetH2Mass() {
  return hMass;
}

TH2F* Analysis::GetH2MassError() {
  return hMassError;
}

//TH2F* Analysis::GetH2MassSigma() {
//  return hMassSigma;
//}

TH2F* Analysis::GetH2JES() {
  return hJES;
}

TH2F* Analysis::GetH2JESError() {
  return hJESError;
}

TH2F* Analysis::GetH2MassConstJES() {
  return hMassConstJES;
}

TH2F* Analysis::GetH2MassConstJESError() {
  return hMassConstJESError;
}

//TH2F* Analysis::GetH2MassConstJESSigma() {
//  return hMassConstJESSigma;
//}

TH2F* Analysis::GetH2FSig() {
  return hFSig;
}

TH2F* Analysis::GetH2FSigError() {
  return hFSigError;
}

TH2F* Analysis::GetH2MassfSig() {
  return hMassfSig;
}

TH2F* Analysis::GetH2MassfSigError() {
  return hMassfSigError;
}

//TH2F* Analysis::GetH2MassfSigSigma() {
//  return hMassfSigSigma;
//}

TH2F* Analysis::GetH2JESfSig() {
  return hJESfSig;
}

TH2F* Analysis::GetH2JESfSigError() {
  return hJESfSigError;
}

//TH2F* Analysis::GetH2MassCalibrated() {
//  return hMassCalibrated;
//}
//
//TH2F* Analysis::GetH2MassErrorCalibrated() {
//  return hMassErrorCalibrated;
//}

TString Analysis::GetIdentifier() {
  return fIdentifier;
}

// return the PU weights for the different samples
double Analysis::calcPUWeight_(enum enumForPUWeights sample, short nPU)
{
  // ----
  // default weight to be used for fall11 samples with S6 PU scenario
  // ----
  // 3.54fb-1 !!! precale weighted !!! 73.5mb
  //double weightsPUFall11[] = {0.0953411, 0.0209421, 0.0389577, 0.355278, 1.3007, 1.8981, 1.95124, 1.81828, 1.63994, 1.55631, 1.54116, 1.44885, 1.24065, 0.915466, 0.580092, 0.350426, 0.209707, 0.113744, 0.0509568, 0.0183767, 0.0054363, 0.00141399, 0.000360394, 0.000106747, 3.96548e-05, 1.65259e-05, 6.57705e-06, 2.32168e-06, 7.15271e-07, 1.99901e-07, 5.67217e-08, 1.86344e-08, 7.44984e-09, 3.25689e-09, 1.39639e-09, 5.58487e-10, 2.00145e-10, 6.60545e-11, 1.89381e-11, 4.82982e-12, 1.12634e-12, 2.21706e-13, 4.06965e-14, 6.1115e-15, 9.20171e-16, 1.05937e-16, 1.38254e-17, 1.14592e-18, 1.13769e-19, 43889.5};
  // 3.54fb-1 !!! precale weighted !!! 68mb
  double weightsPUFall11[] = {0.0956436, 0.0219907, 0.0713245, 0.735605, 1.85831, 2.19452, 2.04027, 1.76898, 1.58702, 1.53704, 1.45119, 1.24648, 0.904423, 0.552217, 0.319224, 0.179453, 0.0864062, 0.0326249, 0.0095542, 0.00229068, 0.000501167, 0.000121599, 3.80191e-05, 1.3964e-05, 4.89459e-06, 1.48493e-06, 3.82906e-07, 8.94646e-08, 2.18024e-08, 6.46154e-09, 2.31321e-09, 8.61013e-10, 3.0218e-10, 9.44128e-11, 2.57874e-11, 6.17332e-12, 1.27288e-12, 2.34488e-13, 3.65422e-14, 4.94041e-15, 5.96064e-16, 5.92555e-17, 5.36326e-18, 3.87753e-19, 2.74435e-20, 1.45014e-21, 8.48143e-23, 3.07617e-24, 1.3049e-25}; //, 58771.4};
  // ----
  // fall11 weight with S6 PU scenario shifted up by 5%
  // ----
  // 3.54fb-1 !!! precale weighted !!! 68mb
  double weightsPUFall11Plus05[] = {0.0954358, 0.020682, 0.0484447, 0.47369, 1.50053, 2.01825, 1.99328, 1.80561, 1.61769, 1.54968, 1.51784, 1.3867, 1.12852, 0.773361, 0.467097, 0.277645, 0.157531, 0.0761859, 0.02937, 0.00905285, 0.00233558, 0.000560472, 0.00014654, 4.85625e-05, 1.90022e-05, 7.38414e-06, 2.54485e-06, 7.59523e-07, 2.01995e-07, 5.28236e-08, 1.59273e-08, 5.8689e-09, 2.42327e-09, 9.83602e-10, 3.67109e-10, 1.23672e-10, 3.6646e-11, 9.87495e-12, 2.28817e-12, 4.67308e-13, 8.65061e-14, 1.34004e-14, 1.91938e-15, 2.23009e-16, 2.57593e-17, 2.25592e-18, 2.2207e-19, 1.37666e-20, 1.01363e-21}; //, 50283.5};

  // ----
  // fall11 weight with S6 PU scenario shifted down by 5%
  // ----
  // 3.54fb-1 !!! precale weighted !!! 68mb
  double weightsPUFall11Minus05[] = {0.0962068, 0.0254305, 0.10942, 1.09375, 2.23819, 2.34231, 2.05248, 1.72372, 1.56859, 1.50191, 1.34308, 1.04368, 0.658236, 0.373059, 0.206629, 0.0997198, 0.0369724, 0.0103136, 0.00227993, 0.000454437, 0.000101297, 3.00008e-05, 1.01979e-05, 3.2025e-06, 8.38403e-07, 1.8544e-07, 3.81727e-08, 8.90286e-09, 2.62892e-09, 8.75387e-10, 2.83582e-10, 8.20241e-11, 2.07262e-11, 4.47805e-12, 8.2373e-13, 1.30017e-13, 1.73391e-14, 2.0282e-15, 1.9709e-16, 1.63193e-17, 1.18443e-18, 6.95731e-20, 3.6548e-21, 1.50639e-22, 5.9703e-24, 1.73528e-25, 5.48351e-27, 1.0555e-28, 2.33406e-30}; //, 73177.8};

  // ----
  // default weight to be used for summer11 samples with S4 PU scenario
  // ----
  // 3.54 fb-1 !!! precale weighted !!! 73.5mb
  //double weightsPUSummer11[] = {0.015788, 0.168999, 0.398126, 0.720134, 1.04936, 1.31397, 1.48381, 1.5613, 1.57478, 1.54969, 1.50736, 1.46098, 1.41641, 1.3734, 1.33783, 1.30218, 1.26742, 1.23224, 1.19459, 1.15701, 1.11803, 1.07132, 1.02626, 0.982213, 0.936123, 0.886997, 0.840387, 0.783362, 0.7366, 0.697965, 0.645112, 0.586587, 0.536603, 0.514893, 0.46368};
  // 3.54fb-1 !!! precale weighted !!! 68mb
  double weightsPUSummer11[] = {0.022778, 0.231718, 0.517403, 0.888432, 1.23288, 1.47585, 1.59957, 1.62107, 1.57897, 1.50264, 1.41361, 1.32374, 1.23758, 1.15455, 1.07948, 1.00629, 0.936187, 0.868571, 0.802406, 0.739714, 0.679652, 0.618656, 0.562462, 0.510472, 0.460945, 0.413436, 0.370475, 0.326335, 0.289734, 0.259019, 0.225714, 0.193379, 0.166592, 0.150473, 0.127517};
  // ----
  // summer11 weight with S6 PU scenario shifted up by 5%
  // ----
  // 3.54fb-1 !!! precale weighted !!! 68mb
  double weightsPUSummer11Plus05[] = {0.0181301, 0.190491, 0.439904, 0.780434, 1.11676, 1.37519, 1.52946, 1.58712, 1.58037, 1.53624, 1.47625, 1.4131, 1.35212, 1.29288, 1.24081, 1.18892, 1.13828, 1.08789, 1.03619, 0.985577, 0.93491, 0.879107, 0.82611, 0.775364, 0.724452, 0.67272, 0.624431, 0.570059, 0.524814, 0.486735, 0.440211, 0.391576, 0.350349, 0.32874, 0.289457};

  // ----
  // summer11 weight with S6 PU scenario shifted down by 5%
  // ----
  // 3.54fb-1 !!! precale weighted !!! 68mb
  double weightsPUSummer11Minus05[] = {0.0285702, 0.28113, 0.606747, 1.0079, 1.35549, 1.57602, 1.66289, 1.64388, 1.56397, 1.45448, 1.33669, 1.22153, 1.11292, 1.01022, 0.917654, 0.829964, 0.748293, 0.672149, 0.600689, 0.535311, 0.475158, 0.417591, 0.366348, 0.320643, 0.279063, 0.241112, 0.208012, 0.176312, 0.150554, 0.12939, 0.108351, 0.0891757, 0.0737799, 0.0639892, 0.0520639};

  switch(sample){
  case kFall11:
    return (nPU < int(sizeof(weightsPUFall11)/sizeof(double)) && nPU > -1) ? weightsPUFall11[nPU] : 0. ;
  case kFall11Plus05:
    return (nPU < int(sizeof(weightsPUFall11Plus05)/sizeof(double)) && nPU > -1) ? weightsPUFall11Plus05[nPU] : 0. ;
  case kFall11Minus05:
    return (nPU < int(sizeof(weightsPUFall11Minus05)/sizeof(double)) && nPU > -1) ? weightsPUFall11Minus05[nPU] : 0. ;
  case kSummer11:
    return (nPU < int(sizeof(weightsPUSummer11)/sizeof(double)) && nPU > -1) ? weightsPUSummer11[nPU] : 0. ;
  case kSummer11Plus05:
    return (nPU < int(sizeof(weightsPUSummer11Plus05)/sizeof(double)) && nPU > -1) ? weightsPUSummer11Plus05[nPU] : 0. ;
  case kSummer11Minus05:
    return (nPU < int(sizeof(weightsPUSummer11Minus05)/sizeof(double)) && nPU > -1) ? weightsPUSummer11Minus05[nPU] : 0. ;
  case kFall10:
    return 1.;
  default:
    return -1.;
  }
}

double Analysis::calcPDFWeight_(int whichPDFUncertainty, bool upVariation, double x1, int id1, float Q, double x2, int id2)
{
  if(whichPDFUncertainty < 0)
    return 1.;

  const std::string NAME = "cteq6mE";
  LHAPDF::initPDFSetM(1, NAME, LHAPDF::LHGRID);
  LHAPDF::initPDFSetM(2, NAME, LHAPDF::LHGRID);

  LHAPDF::initPDFM(1, 0);
  if(upVariation) LHAPDF::initPDFM(2, 2*whichPDFUncertainty-1);
  else            LHAPDF::initPDFM(2, 2*whichPDFUncertainty);
  //std::cout << whichPDFUncertainty << ": " << LHAPDF::xfxM(2, x1, Q, id1) << " " << LHAPDF::xfxM(2, x2, Q, id2) << " " << LHAPDF::xfxM(1, x1, Q, id1) << " " << LHAPDF::xfxM(1, x2, Q, id2) << std::endl;
  return (LHAPDF::xfxM(2, x1, Q, id1)*LHAPDF::xfxM(2, x2, Q, id2)/(LHAPDF::xfxM(1, x1, Q, id1)*LHAPDF::xfxM(1, x2, Q, id2)));
}

double Analysis::calcBTagWeight_(int Njet, short * pdgId, TClonesArray * jets)
{
  double bTaggingEfficiency = 0., bTaggingEfficiency_scaled = 0.;
  //double bTaggingEfficiency_scaled_EffUp = 0., bTaggingEfficiency_scaled_EffDown = 0., bTaggingEfficiency_scaled_MisUp = 0., bTaggingEfficiency_scaled_MisDown = 0.;
  double pt, eta, eff, effyScale_pt; //, effVariation_pt, misTagScale_pt, misVariation_pt;

  std::vector<double> oneMinusBEffies(0) , oneMinusBEffies_scaled(0) ; //, oneMinusBEffies_scaled_EffUp(0) , oneMinusBEffies_scaled_EffDown(0) ;
  std::vector<double> oneMinusBMistags(0), oneMinusBMistags_scaled(0); //, oneMinusBMistags_scaled_MisUp(0), oneMinusBMistags_scaled_MisDown(0);
  for(int i = 0; i < Njet; ++i){
    pt  = ((TLorentzVector*)jets->At(i))->Pt();
    eta = ((TLorentzVector*)jets->At(i))->Eta();

    if(pt > 670.)
      effyScale_pt    = 0.901615*((1.+(0.552628*670.))/(1.+(0.547195*670.)));
    if(pt < 30.)
      effyScale_pt    = 0.901615*((1.+(0.552628*30.))/(1.+(0.547195*30.)));
    else
      effyScale_pt    = 0.901615*((1.+(0.552628*pt))/(1.+(0.547195*pt)));
    //effVariation_pt = bTagEffScaleFactor->GetBinError(bTagEffScaleFactor->FindBin(pt));

    if(pdgId[i] == 5 || pdgId[i] == -5){
      eff = bTagEff->GetBinContent(bTagEff->FindBin(pt,std::abs(eta)));
      oneMinusBEffies               .push_back(1.- eff);
      oneMinusBEffies_scaled        .push_back(1.-(eff* effyScale_pt));
      //if(isDefaultSample){
      //  oneMinusBEffies_scaled_EffUp  .push_back(1.-(eff*(effyScale_pt+effVariation_pt)));
      //  oneMinusBEffies_scaled_EffDown.push_back(1.-(eff*(effyScale_pt-effVariation_pt)));
      //}
    }
    else if(pdgId[i] == 4 || pdgId[i] == -4){
      eff = cTagEff->GetBinContent(cTagEff->FindBin(pt,std::abs(eta)));
      oneMinusBMistags               .push_back(1.- eff);
      oneMinusBMistags_scaled        .push_back(1.-(eff* effyScale_pt));
      //if(pt<240){
      //	oneMinusBMistags_scaled_MisUp  .push_back(1.-(eff*(effyScale_pt+(2*effVariation_pt))));
      //	oneMinusBMistags_scaled_MisDown.push_back(1.-(eff*(effyScale_pt-(2*effVariation_pt))));
      //}
      //else{
      //	oneMinusBMistags_scaled_MisUp  .push_back(1.-(eff*(effyScale_pt+effVariation_pt)));
      //	oneMinusBMistags_scaled_MisDown.push_back(1.-(eff*(effyScale_pt-effVariation_pt)));
      //}
    }
    else{
      eff = lTagEff->GetBinContent(lTagEff->FindBin(pt,std::abs(eta)));
      oneMinusBMistags               .push_back(1.- eff);
      oneMinusBMistags_scaled        .push_back(1.-(eff* (((0.948463+(0.00288102*pt))+(-7.98091e-06*(pt*pt)))+(5.50157e-09*(pt*(pt*pt)))) ));
      //if(isDefaultSample){
      //  oneMinusBMistags_scaled_MisUp  .push_back(1.-(eff* (((0.997077+(0.00473953*pt))+(-1.34985e-05*(pt*pt)))+(1.0032e-08*(pt*(pt*pt)))) ));
      //  oneMinusBMistags_scaled_MisDown.push_back(1.-(eff* (((0.899715+(0.00102278*pt))+(-2.46335e-06*(pt*pt)))+(9.71143e-10*(pt*(pt*pt)))) ));
      //}
    }
  }
  bTaggingEfficiency        = eventBTagProbability_(oneMinusBEffies       , oneMinusBMistags       );
  bTaggingEfficiency_scaled = eventBTagProbability_(oneMinusBEffies_scaled, oneMinusBMistags_scaled);

  //if(bTaggingEfficiency_scaled/bTaggingEfficiency > 5.){
  //  std::cout << "Ratio: " << bTaggingEfficiency_scaled/bTaggingEfficiency << " !!!" << std::endl;
  //  std::cout << "Eff for: " << Njet << std::endl;
  //  for(int i = 0; i < Njet; ++i)
  //    std::cout << pdgId[i] << ", ";
  //  std::cout << std::endl;
  //  std::cout << "(1-eff " << effyScale_pt << "): ";
  //  for(std::vector<double>::const_iterator i = oneMinusBEffies.begin(); i != oneMinusBEffies.end(); ++i)
  //    std::cout << (*i) << ", ";
  //  std::cout << std::endl;
  //  std::cout << "(1-eff scaled): ";
  //  for(std::vector<double>::const_iterator i = oneMinusBEffies_scaled.begin(); i != oneMinusBEffies_scaled.end(); ++i)
  //    std::cout << (*i) << ", ";
  //  std::cout << std::endl;
  //  std::cout << "(1-mis): ";
  //  for(std::vector<double>::const_iterator i = oneMinusBMistags.begin(); i != oneMinusBMistags.end(); ++i)
  //    std::cout << (*i) << ", ";
  //  std::cout << std::endl;
  //  std::cout << "(1-mis scaled): ";
  //  for(std::vector<double>::const_iterator i = oneMinusBMistags_scaled.begin(); i != oneMinusBMistags_scaled.end(); ++i)
  //    std::cout << (*i) << ", ";
  //  std::cout << std::endl;
  //  bTaggingEfficiency        = eventBTagProbability_(oneMinusBEffies       , oneMinusBMistags       , true);
  //  std::cout << "Scaled Eff: " << std::endl;
  //  bTaggingEfficiency_scaled = eventBTagProbability_(oneMinusBEffies_scaled, oneMinusBMistags_scaled, true);
  //}

  //if(isDefaultSample){
  //  bTaggingEfficiency_scaled_EffUp   += eventBTagProbability_(oneMinusBEffies_scaled_EffUp  , oneMinusBMistags_scaled        );
  //  bTaggingEfficiency_scaled_EffDown += eventBTagProbability_(oneMinusBEffies_scaled_EffDown, oneMinusBMistags_scaled        );
  //  bTaggingEfficiency_scaled_MisUp   += eventBTagProbability_(oneMinusBEffies_scaled        , oneMinusBMistags_scaled_MisUp  );
  //  bTaggingEfficiency_scaled_MisDown += eventBTagProbability_(oneMinusBEffies_scaled        , oneMinusBMistags_scaled_MisDown);
  //}
  //std::cout << bTaggingEfficiency << " " << bTaggingEfficiency_scaled << " " << bTaggingEfficiency_scaled/bTaggingEfficiency << std::endl;
  return bTaggingEfficiency_scaled/bTaggingEfficiency;
}

// calculate the probability of b-tagging one event with 2 b-tags
double
Analysis::eventBTagProbability_(std::vector<double> &oneMinusBEffies, std::vector<double> &oneMinusBMistags, bool verbose){
  double bTaggingEfficiency = 1.;
  double tmp = 1.;

  if(verbose) std::cout << bTaggingEfficiency << std::endl;

  // subtract probability that no jet is tagged
  for(std::vector<double>::const_iterator eff = oneMinusBEffies.begin(); eff != oneMinusBEffies.end(); ++eff)
    tmp *= (*eff);
  for(std::vector<double>::const_iterator mis = oneMinusBMistags.begin(); mis != oneMinusBMistags.end(); ++mis)
    tmp *= (*mis);
  bTaggingEfficiency -= tmp;

  if(verbose) std::cout << bTaggingEfficiency << std::endl;

  // subtract probability that 1 bJet is tagged
  for(std::vector<double>::const_iterator eff = oneMinusBEffies.begin(); eff != oneMinusBEffies.end(); ++eff){
    tmp = 1.-(*eff);
    for(std::vector<double>::const_iterator eff2 = oneMinusBEffies.begin(); eff2 != oneMinusBEffies.end(); ++eff2)
      if(eff != eff2) tmp *= (*eff2);
    for(std::vector<double>::const_iterator mis = oneMinusBMistags.begin(); mis != oneMinusBMistags.end(); ++mis)
      tmp *= (*mis);
    bTaggingEfficiency -= tmp;
  }

  if(verbose) std::cout << bTaggingEfficiency << std::endl;

  // subtract probability that 1 non-bJet is tagged
  for(std::vector<double>::const_iterator mis = oneMinusBMistags.begin(); mis != oneMinusBMistags.end(); ++mis){
    tmp = 1.-(*mis);
    for(std::vector<double>::const_iterator eff = oneMinusBEffies.begin(); eff != oneMinusBEffies.end(); ++eff)
      tmp *= (*eff);
    for(std::vector<double>::const_iterator mis2 = oneMinusBMistags.begin(); mis2 != oneMinusBMistags.end(); ++mis2)
      if(mis != mis2) tmp *= (*mis2);
    bTaggingEfficiency -= tmp;
  }

  if(verbose) std::cout << bTaggingEfficiency << std::endl;

  return bTaggingEfficiency;
}
