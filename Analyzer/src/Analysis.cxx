#include "Analysis.h"

#include "TSystem.h"
#include "TEventList.h"

#include "LHAPDF/LHAPDF.h"
#include <string>
#include <stdlib.h>

Analysis::Analysis(po::variables_map vm) :
  _samplePath("/scratch/hh/lustre/cms/user/eschliec/TopMass/19/"),
  //samplePath"/scratch/hh/dust/naf/cms/user/eschliec/TopMass/19/"),
  _fIdentifier(vm["input"].as<std::string>()),
  _fMethod    (vm["method"].as<std::string>()),
  _fBins      (vm["bins"  ].as<int>()),
  _fLumi      (vm["lumi"  ].as<double>()),
  _fTree(0),
  _tTree(0),
  _tTreeBkg(0),
  _tempFilePath(gSystem->Getenv("TMPDIR")),
  _tempFile(0)
{
  _fFile += _samplePath;
  _fFile += _fIdentifier;
  _fFile += ".root";

  _fChain = new TChain("FullHadTreeWriter/tree");
  _fChain->Add(_fFile);

  CreateHistos();
  ReadConfigFromXMLFile();
}

Analysis::Analysis(TString identifier, TString file, TString method, int bins, double lumi) :
  _samplePath("/scratch/hh/lustre/cms/user/eschliec/TopMass/19/"),
  //samplePath"/scratch/hh/dust/naf/cms/user/eschliec/TopMass/19/"),
  _fIdentifier(identifier), _fMethod(method), _fBins(bins), _fLumi(lumi),
  _fTree(0), _tTree(0), _tTreeBkg(0),
  _tempFilePath(gSystem->Getenv("TMPDIR")),
  _tempFile(0)
{
  _fFile += _samplePath;
  _fFile += _fIdentifier;
  _fFile += ".root";

  _fChain = new TChain("FullHadTreeWriter/tree");
  _fChain->Add(_fFile);

  CreateHistos();
  ReadConfigFromXMLFile();
}

Analysis::~Analysis()
{
  delete _fChain;

  delete _fTree;
  delete _tTree;
  delete _tTreeBkg;
  delete _hEntries;
  delete _hMass;
  delete _hMassError;
  delete _hJES;
  delete _hJESError;
  delete _hMassConstJES;
  delete _hMassConstJESError;
  delete _hFSig;
  delete _hFSigError;
  delete _hMassfSig;
  delete _hMassfSigError;
  delete _hJESfSig;
  delete _hJESfSigError;
  delete _bTagEff;
  delete _cTagEff;
  delete _lTagEff;
  _tempFile->Close();
  delete _tempFile;
}

void Analysis::Analyze(po::variables_map vm) {

  std::cout << "Analyze " << _fIdentifier << " with method " << _fMethod << std::endl;

  CreateRandomSubset();

  MassAnalyzer* fAnalyzer = 0;

  if (!strcmp(_fMethod, "GenMatch")) {
    fAnalyzer = new GenMatchAnalyzer(_fIdentifier, _fTree);
  }
  else if (!strcmp(_fMethod, "MVA")) {
    fAnalyzer = new MVAAnalyzer(_fIdentifier, _fTree);
  }
  else if (!strcmp(_fMethod, "Ideogram")) {
    fAnalyzer = new IdeogramAnalyzer(_fIdentifier, _fTree);
  }
  else if (!strcmp(_fMethod, "RooFit")) {
    fAnalyzer = new RooFitTemplateAnalyzer(_fIdentifier, _fTree);
  }
  else {
    return;
  }

  Helper* helper = new Helper(_fBins);
  helper->SetTDRStyle();
  delete helper;

  TCanvas* canvas = new TCanvas("canvas", "Top mass", 900, 600);
  canvas->cd();

  double smearBins = 1.;
  double rangeX = 10.;
  double rangeY = 10.;

  double smearX = smearBins/_fBins*rangeX;
  double smearY = smearBins/_fBins*rangeY;

  double minEntries = 25;

  TString observableX = "dRbb";
  TString observableY = "dRbb";

  for(int i = 0; i < _fBins; i++) {
    for(int j = 0; j < _fBins; j++) {
      // calculate cuts
      std::stringstream stream;
      stream << (rangeX/_fBins)*(i)-smearX << "<" << observableX << "&"
             << observableX << "<" << (rangeX/_fBins)*(i+1)+smearX << " & "
             << rangeY/_fBins*(j)-smearY << "<" << observableY << "&"
             << observableY << "<" << rangeY/_fBins*(j+1)+smearY;
      TString cuts(stream.str());

      if (!strcmp(_fMethod, "GenMatch")) {
        cuts += " & comboType == 1";
      }
      else if (!strcmp(_fMethod, "MVA")) {
        cuts += " & mvaDisc > 0";
      }
      else if (!strcmp(_fMethod, "Ideogram")) {
//        //cuts += " & target == 1";
//        //cuts += " & (target == 0 | target == 1)";
//        //cuts += " & Njet == 6";
//        cuts += " & probs[0] > 0.09";
//        cuts += " & dRbb > 1.5";
//        //cuts += " & nVertex >= 5";
//        //cuts += " & nVertex <= 5";
//        cuts += " & topMasses[0] > 100. & topMasses[0] < 550.";
      }
      else if (!strcmp(_fMethod, "RooFit")) {
      }

      int entries = _fTree->GetEntries(cuts);
      std::cout << cuts << std::endl;
      std::cout << entries << std::endl;

      _hEntries->SetCellContent(i+1, j+1, entries);

      if (entries > minEntries) {
        fAnalyzer->Analyze(cuts, i, j, vm);

        _hMass     ->SetCellContent(i+1, j+1, fAnalyzer->GetMass());
        _hMassError->SetCellContent(i+1, j+1, fAnalyzer->GetMassError());
        _hJES      ->SetCellContent(i+1, j+1, fAnalyzer->GetJES());
        _hJESError ->SetCellContent(i+1, j+1, fAnalyzer->GetJESError());
        _hMassConstJES     ->SetCellContent(i+1, j+1, fAnalyzer->GetMassConstJES());
        _hMassConstJESError->SetCellContent(i+1, j+1, fAnalyzer->GetMassConstJESError());
        _hFSig     ->SetCellContent(i+1, j+1, fAnalyzer->GetFSig());
        _hFSigError->SetCellContent(i+1, j+1, fAnalyzer->GetFSigError());
        _hMassfSig     ->SetCellContent(i+1, j+1, fAnalyzer->GetMassfSig());
        _hMassfSigError->SetCellContent(i+1, j+1, fAnalyzer->GetMassfSigError());
        _hJESfSig      ->SetCellContent(i+1, j+1, fAnalyzer->GetJESfSig());
        _hJESfSigError ->SetCellContent(i+1, j+1, fAnalyzer->GetJESfSigError());

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
  _hEntries->Draw("COLZ");

  canvas->cd(2);
  _hMass->Draw("COLZ,TEXT");
  _hMass->SetAxisRange(_hMass->GetMinimum(0.05), _hMass->GetMaximum(), "Z");

  canvas->cd(3);
  _hMassError->Draw("COLZ,TEXT");
  _hMassError->SetAxisRange(0.05, 5, "Z");

  TString path("plot/"); path += _fMethod; path += "_"; path += _fIdentifier; path += ".eps";
  canvas->Print(path);

  delete canvas;
  delete fAnalyzer;
}

void Analysis::CreateHistos() {
  Helper* helper = new Helper(_fBins);

  _hEntries = helper->GetH2("Entries");

  _hMass = helper->GetH2("Mass");
  _hMassError = helper->GetH2("MassError");
  _hJES = helper->GetH2("JES");
  _hJESError = helper->GetH2("JESError");
  
  _hMassConstJES = helper->GetH2("MassConstJES");
  _hMassConstJESError = helper->GetH2("MassConstJESError");
  
  _hFSig = helper->GetH2("fSig");
  _hFSigError = helper->GetH2("fSigError");
  _hMassfSig = helper->GetH2("MassfSig");
  _hMassfSigError = helper->GetH2("MassfSigError");
  _hJESfSig = helper->GetH2("JESfSig");
  _hJESfSigError = helper->GetH2("JESfSigError");
  
  TFile * bTagFile = TFile::Open(_samplePath+TString("bTagFile.root"));
  gROOT->cd();
  _bTagEff = (TH2F*)bTagFile->Get("histb")->Clone();
  _cTagEff = (TH2F*)bTagFile->Get("histc")->Clone();
  _lTagEff = (TH2F*)bTagFile->Get("histl")->Clone();
  bTagFile->Close();

  delete bTagFile;
  delete helper;
}

void Analysis::CreateRandomSubset() {
  std::cout << "Create random subset..." << std::endl;
  SetBranchStatuses(_fChain);
  /*
  _fChain->SetBranchStatus("*",0);
  _fChain->SetBranchStatus("nCombos", 1);
  _fChain->SetBranchStatus("comboTypes", 1);
  _fChain->SetBranchStatus("topMasses", 1);
  _fChain->SetBranchStatus("topMass", 1);
  _fChain->SetBranchStatus("w1Mass", 1);
  _fChain->SetBranchStatus("w2Mass", 1);
  _fChain->SetBranchStatus("probs", 1);
  _fChain->SetBranchStatus("prob", 1);
  _fChain->SetBranchStatus("dRbb", 1);
  
  _fChain->SetBranchStatus("Njet", 1);
  _fChain->SetBranchStatus("jets", 1);
  //_fChain->SetBranchStatus("bTag_CSV", 1);
  _fChain->SetBranchStatus("partonFlavour", 1);

  _fChain->SetBranchStatus("runNumber", 1);
  _fChain->SetBranchStatus("luminosityBlockNumber", 1);
  _fChain->SetBranchStatus("eventNumber", 1);

  _fChain->SetBranchStatus("nPU", 1);
  _fChain->SetBranchStatus("nPUTru", 1);

  _fChain->SetBranchStatus("MCweight", 1);
  
  if(_fFile.Contains("PDF")){
    _fChain->SetBranchStatus("id1", 1);
    _fChain->SetBranchStatus("id2", 1);
    _fChain->SetBranchStatus("x1", 1);
    _fChain->SetBranchStatus("x2", 1);
    _fChain->SetBranchStatus("Q", 1);
  }
  */
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

  TString triggerUncertainty = "";
  if(_fIdentifier.Contains("jet4_02")) triggerUncertainty += " & jets[3].Pt() > 62.";
  if(_fIdentifier.Contains("jet4_05")) triggerUncertainty += " & jets[3].Pt() > 65.";
  if(_fIdentifier.Contains("jet4_10")) triggerUncertainty += " & jets[3].Pt() > 70.";
  if(_fIdentifier.Contains("jet5_02")) triggerUncertainty += " & jets[4].Pt() > 52.";
  if(_fIdentifier.Contains("jet5_05")) triggerUncertainty += " & jets[4].Pt() > 55.";
  if(_fIdentifier.Contains("jet5_10")) triggerUncertainty += " & jets[4].Pt() > 60.";
  if(_fIdentifier.Contains("jet6_02")) triggerUncertainty += " & jets[5].Pt() > 42.";
  if(_fIdentifier.Contains("jet6_05")) triggerUncertainty += " & jets[5].Pt() > 45.";
  if(_fIdentifier.Contains("jet6_10")) triggerUncertainty += " & jets[5].Pt() > 50.";

  if (_fLumi!=0) {
    _tempFile = TFile::Open(_tempFilePath+TString("/tempTree.root"), "READ");
    if(_tempFile && !_tempFile->IsZombie()){
      _tTree    = (TTree*)_tempFile->Get(TString("fullTree_")+_fIdentifier);
      _tTreeBkg = (TTree*)_tempFile->Get("fullTreeBkg");
      std::cout << "Reading: " << _tempFile->GetName() << std::endl;
    }
    else{
      _tempFile = TFile::Open(_tempFilePath+TString("/tempTree.root"), "RECREATE", "", 0);
      std::cout << "Creating: " << _tempFile->GetName() << std::endl;
    }
    if(!_tTree || !_tTreeBkg){
      _tempFile->ReOpen("UPDATE");
      std::cout << "Updating: " << _tempFile->GetName() << std::endl;

      if(!_tTree){
        _fChain->Draw(">>selectedEvents",_selection+triggerUncertainty,"goff");
        TEventList* selectedEvents = (TEventList*)gDirectory->Get("selectedEvents");

        _fChain->SetBranchAddress("nCombos", &nCombos);
        _fChain->SetBranchAddress("comboTypes", comboTypes);
        _fChain->SetBranchAddress("topMasses", topMasses);
        _fChain->SetBranchAddress("topMass", &topMass);
        _fChain->SetBranchAddress("w1Mass", w1Masses);
        _fChain->SetBranchAddress("w2Mass", w2Masses);
        _fChain->SetBranchAddress("probs", probs);
        _fChain->SetBranchAddress("prob", &prob);
        _fChain->SetBranchAddress("dRbb", &dRbb);

        _fChain->SetBranchAddress("Njet", &Njet);
        _fChain->SetBranchAddress("jets", &jets);
        //fChain->SetBranchAddress("bTag_CSV", bTag);
        _fChain->SetBranchAddress("partonFlavour", pdgId);

        _fChain->SetBranchAddress("runNumber", &runNumber);
        _fChain->SetBranchAddress("luminosityBlockNumber", &luminosityBlockNumber);
        _fChain->SetBranchAddress("eventNumber", &eventNumber);

        _fChain->SetBranchAddress("dRbb", &dRbb);

        _fChain->SetBranchAddress("nPU", &nPU);
        _fChain->SetBranchAddress("nPUTru", &nPUTru);

        _fChain->SetBranchAddress("MCweight", &MCweight);

        if(_fFile.Contains("PDF")){
          _fChain->SetBranchAddress("id1", &id1);
          _fChain->SetBranchAddress("id2", &id2);
          _fChain->SetBranchAddress("x1", &x1);
          _fChain->SetBranchAddress("x2", &x2);
          _fChain->SetBranchAddress("Q", &Q);
        }

        _tTree = _fChain->CloneTree(0);
        _tTree->SetName(TString("fullTree_")+_fIdentifier);

        for(int idx = 0 , l = selectedEvents->GetN(); idx < l; ++idx){
          _fChain->GetEntry(selectedEvents->GetEntry(idx));
          nCombos = 1;
          _tTree->Fill();
        }
        AddWeights(_tTree);
        _tTree->Write(0,TObject::kOverwrite);
        delete _tTree;
      }

      if(!_tTreeBkg){
        TFile * fileBkg = TFile::Open(_samplePath+TString("QCDEstimationMix_2011_skimmed2.root"));
        if(TString(fileBkg->GetName()).Contains("_skimmed")) _tTreeBkg = (TTree*)fileBkg->Get("tree");
        else  _tTreeBkg = (TTree*)fileBkg->Get("analyzeFullHadEventMixer/tree");

        _tTreeBkg->SetBranchStatus("*",0);
        _tTreeBkg->SetBranchStatus("topMasses", 1);
        _tTreeBkg->SetBranchStatus("topMass", 1);
        _tTreeBkg->SetBranchStatus("w1Mass", 1);
        _tTreeBkg->SetBranchStatus("w2Mass", 1);
        _tTreeBkg->SetBranchStatus("probs", 1);
        _tTreeBkg->SetBranchStatus("prob", 1);
        _tTreeBkg->SetBranchStatus("dRbb", 1);

        _tTreeBkg->SetBranchStatus("jets", 1);

        _tTreeBkg->SetBranchStatus("runNumber", 1);
        _tTreeBkg->SetBranchStatus("luminosityBlockNumber", 1);
        _tTreeBkg->SetBranchStatus("eventNumber", 1);

        _tempFile->cd();

        TString selBkg = _selection; selBkg.ReplaceAll("dRbb","dRbb[0]");
        TTree* tempTreeBkg = _tTreeBkg->CopyTree(selBkg);
        tempTreeBkg->SetName("fullTreeBkg");
        AddWeights(tempTreeBkg,true);
        tempTreeBkg->Write(0,TObject::kOverwrite);
        delete tempTreeBkg;

        delete _tTreeBkg;
        fileBkg->Close();
        delete fileBkg;
      }

      _tempFile->Close();
      _tempFile = TFile::Open(_tempFilePath+TString("/tempTree.root"),"READ");
      _tTree    = (TTree*)_tempFile->Get(TString("fullTree_")+_fIdentifier);
      _tTreeBkg = (TTree*)_tempFile->Get("fullTreeBkg");
    }

    TRandom3* myRandom = new TRandom3(0);
    std::cout << "Random seed: " << myRandom->GetSeed() << std::endl;

    _tTree   ->SetCacheSize(10000000000);
    _tTreeBkg->SetCacheSize(10000000000);

    _tTree->SetBranchStatus("*",0);
    _tTree->SetBranchStatus("CombinedWeight",1);
    _tTree->SetBranchStatus("meanWMass",1);
    _tTree->SetBranchStatus("topMasses",1);
    _tTree->SetBranchStatus("topMass",1);
    _tTree->SetBranchStatus("w1Mass",1);
    _tTree->SetBranchStatus("w2Mass",1);
    _tTree->SetBranchStatus("probs",1);
    _tTree->SetBranchStatus("prob",1);
    _tTree->SetBranchStatus("dRbb",1);
    _tTree->SetBranchStatus("jets",1);
    _tTree->SetBranchStatus("nCombos",1);
    _tTree->SetBranchStatus("comboTypes",1);
    _tTree->SetBranchStatus("runNumber", 1);
    _tTree->SetBranchStatus("luminosityBlockNumber", 1);
    _tTree->SetBranchStatus("eventNumber", 1);

    _tTree->SetBranchAddress("CombinedWeight",&CombinedWeight);
    _tTree->SetBranchAddress("meanWMass",&meanWMass);
    _tTree->SetBranchAddress("topMasses",topMasses);
    _tTree->SetBranchAddress("topMass",&topMass);
    _tTree->SetBranchAddress("w1Mass",w1Masses);
    _tTree->SetBranchAddress("w2Mass",w2Masses);
    _tTree->SetBranchAddress("probs",probs);
    _tTree->SetBranchAddress("prob",&prob);
    _tTree->SetBranchAddress("dRbb",&dRbb);
    _tTree->GetBranch("jets")->SetAutoDelete(kFALSE);
    _tTree->SetBranchAddress("jets", &jets);
    _tTree->SetBranchAddress("nCombos",&nCombos);
    _tTree->SetBranchAddress("comboTypes",comboTypes);
    _tTree->SetBranchAddress("runNumber", &runNumber);
    _tTree->SetBranchAddress("luminosityBlockNumber", &luminosityBlockNumber);
    _tTree->SetBranchAddress("eventNumber", &eventNumber);

    _tTreeBkg->SetBranchAddress("meanWMass",&meanWMass);
    _tTreeBkg->SetBranchAddress("topMasses", topMasses);
    _tTreeBkg->SetBranchAddress("topMass", &topMass);
    _tTreeBkg->SetBranchAddress("w1Mass", w1Masses);
    _tTreeBkg->SetBranchAddress("w2Mass", w2Masses);
    _tTreeBkg->SetBranchAddress("probs", probs);
    _tTreeBkg->SetBranchAddress("prob", &prob);
    _tTreeBkg->SetBranchAddress("dRbb", &dRbb);
    _tTreeBkg->GetBranch("jets")->SetAutoDelete(kFALSE);
    _tTreeBkg->SetBranchAddress("jets", &jets);
    _tTreeBkg->SetBranchAddress("nCombos",&nCombos);
    
    _tTreeBkg->SetBranchAddress("runNumber", &runNumber);
    _tTreeBkg->SetBranchAddress("luminosityBlockNumber", &luminosityBlockNumber);
    _tTreeBkg->SetBranchAddress("eventNumber", &eventNumber);
    
    double maxMCWeight = _tTree->GetMaximum("CombinedWeight");

    if (maxMCWeight == -1) { std::cout << "Running over data?" << std::endl; }
    
    double fSig = 0.504; // 0.539; //
    if(_fFile.Contains("fSig_Up"))
      fSig += 0.10;
    else if(_fFile.Contains("fSig_Down"))
      fSig -= 0.10;
    else if(_fFile.Contains("BJES_Up"))
      fSig *= 1.0124;
    else if(_fFile.Contains("BJES_Down"))
      fSig *= 0.9895;
    else if(_fFile.Contains("JES_Up"))
      fSig *= 1.1042;
    else if(_fFile.Contains("JES_Down"))
      fSig *= 0.8967;
    else if(_fFile.Contains("JER_Up"))
      fSig *= 0.9816;
    else if(_fFile.Contains("JER_Down"))
      fSig *= 1.0192;
    else if(_fFile.Contains("BTAG_Up"))
      fSig *= 1.0509;
    else if(_fFile.Contains("BTAG_Down"))
      fSig *= 0.9565;
    else if(_fFile.Contains("Scale_Up"))
      fSig *= 0.8726;
    else if(_fFile.Contains("Scale_Down"))
      fSig *= 1.1184;
    else if(_fFile.Contains("Matching_Up"))
      fSig *= 0.9718;
    else if(_fFile.Contains("P11_NoCR"))
      fSig *= 1.0294;
    else if(_fFile.Contains("P11mpiHi"))
      fSig *= 1.0202;
    else if(_fFile.Contains("P11TeV"))
      fSig *= 1.0234;
    else if(_fFile.Contains("jet4_02"))
      fSig *= 0.9256;
    else if(_fFile.Contains("jet5_02"))
      fSig *= 0.9431;
    else if(_fFile.Contains("jet6_02"))
      fSig *= 0.9305;

    int permsMC  = _tTree   ->GetEntries();
    int permsBkg = _tTreeBkg->GetEntries();
    int eventsPE = myRandom->Poisson(2410./3544.844*_fLumi); // add poisson 
    int eventsDrawn = 0;
    
    _tempFile->ReOpen("UPDATE");
    _fTree = _tTree->CloneTree(0);
    _fTree->SetName("tree");
    // remove unnecessary branches from fTree
    _fTree->SetBranchStatus("CombinedWeight",0);
    _fTree->SetBranchStatus("jets",0);
    //fTree->SetBranchStatus("runNumber", 0);
    //fTree->SetBranchStatus("luminosityBlockNumber", 0);
    //fTree->SetBranchStatus("eventNumber", 0);

    int signalDrawn = 0, backgroundDrawn = 0;
    if(_fLumi>0) {
      //int drawCounter = 0;
      while (eventsDrawn < eventsPE) {
        double fSigRndm = myRandom->Uniform(0.,1.);
        if(fSigRndm < fSig){
          int oldEventsDrawn = eventsDrawn;
          do {
            int drawn = myRandom->Integer(permsMC);
            //++drawCounter;
            _tTree->GetEntry(drawn);
            //std::cout << CombinedWeight << " " << maxMCWeight << std::endl;
            if (CombinedWeight > myRandom->Uniform(0., maxMCWeight)) {
              // reduce size of tree before filling
              nCombos = 1;
              //if(fTree->Fill() == -1) std::cout << eventsDrawn << " " << meanWMass << " " << topMasses[0] << " " << topMass << " " << w1Masses[0] << " " << w2Masses[0] << " " << probs[0] << " " << prob << " " << dRbb << " " << nCombos << " " << comboTypes[0] << " " << runNumber << " " << luminosityBlockNumber << " " << eventNumber << " " << std::endl;
              _fTree->Fill();

              ++eventsDrawn;
              ++signalDrawn;
            }
          }
          while(oldEventsDrawn == eventsDrawn);
        }
        else{
          int drawn = myRandom->Integer(permsBkg);
          //++drawCounter;
          _tTreeBkg->GetEntry(drawn);
          // reduce size of tree before filling
          nCombos = 1;
          // set missing variable for background tree
          comboTypes[0] = 0;
          //if(fTree->Fill() == -1) std::cout << eventsDrawn << " " << meanWMass << " " << topMasses[0] << " " << topMass << " " << w1Masses[0] << " " << w2Masses[0] << " " << probs[0] << " " << prob << " " << dRbb << " " << nCombos << " " << comboTypes[0] << " " << runNumber << " " << luminosityBlockNumber << " " << eventNumber << " " << std::endl;
          _fTree->Fill();
          ++eventsDrawn;
          ++backgroundDrawn;
        }
      }
      //std::cout << "DRAWCOUNTER: " << drawCounter << std::endl;
    }
    else {
      for(int ev = 0; ev < permsMC; ++ev){
        _tTree->GetEntry(ev);
        _fTree->Fill();
        ++signalDrawn;
      }
      std::cout << "wanted / available events: " << permsMC*((1.-fSig)/fSig) << " / " << permsBkg << std::endl;
      for(int ev = 0; ev < permsMC*((1.-fSig)/fSig); ++ev){
        int drawn = myRandom->Integer(permsBkg);
        _tTreeBkg->GetEntry(drawn);
        _fTree->Fill();
        ++backgroundDrawn;
      }
    }
    _fTree->Write(0,TObject::kOverwrite);
    std::cout << "Events drawn: " << eventsDrawn << " (sig: " << signalDrawn << ", bkg: " << backgroundDrawn << ") -> fSig: " << double(signalDrawn)/double(signalDrawn+backgroundDrawn) << " (default: " << fSig << ")" << std::endl;
    delete _fTree;
    delete _tTree;
    delete _tTreeBkg;
    _tempFile->Close();
    _tempFile = TFile::Open(_tempFilePath+TString("/tempTree.root"));
    _fTree = (TTree*)_tempFile->Get("tree");
    delete myRandom;
  }
  else {
    _tempFile = new TFile(_tempFilePath+TString("/tempTree.root"), "RECREATE");
    _tempFile->cd();
    _tTree = _fChain->CopyTree(_selection);
    AddWeights(_tTree,true);
    _tempFile->Write(0,TObject::kOverwrite);
    delete _tTree;
    _tempFile->Close();
    _tempFile = TFile::Open(_tempFilePath+TString("/tempTree.root"));
    _fTree = (TTree*)_tempFile->Get("tree");
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

void Analysis::AddWeights(TTree* tempTree, bool isData) {
  std::cout << "Adapting " << tempTree->GetName();
  if(!isData) std::cout << " and adding combined weight";
  std::cout << " ..." << std::endl;

  enumForPUWeights whichSample = kFall10;
  if(!isData){
    if     (_fFile.Contains("S11"        )) whichSample = kSummer11;
    else if(_fFile.Contains("FAST"       )) whichSample = kFall10;
    else if(_fFile.Contains("F11_PU_Up"  )) whichSample = kFall11Plus05;
    else if(_fFile.Contains("F11_PU_Down")) whichSample = kFall11Minus05;
    else if(_fFile.Contains("F11"        )) whichSample = kFall11;
  }
  int whichPDFUncertainty = -1;
  bool upVariation = true; // value does not matter as long as _whichPDFUncertainty_==-1
  if(!isData){
    if(_fFile.Contains("PDF")){
      TString pdfUncert = _fFile;
      upVariation = pdfUncert.Contains("Up") ? true : false;
      pdfUncert = TString(pdfUncert(pdfUncert.Index("PDF",TString::kIgnoreCase)+3,2).Data());
      whichPDFUncertainty = pdfUncert.Atoi();
    }
    std::cout << "whichPDFUncertainty: " << whichPDFUncertainty << " (" << _fFile << ")" << std::endl;
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
  TBranch * br2 = tempTree->Branch("meanWMass", &meanWMass, "meanWMass/D");
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
    br2->Fill();
  }
  jets->Delete();
  delete[] w1Masses;
  delete[] w2Masses;
  std::cout << "Adapted  " << tempTree->GetName();
  if(!isData) std::cout << " and added  combined weight";
  std::cout << " ..." << std::endl;
}

TH2F* Analysis::GetH2Mass() {
  return _hMass;
}

TH2F* Analysis::GetH2MassError() {
  return _hMassError;
}

TH2F* Analysis::GetH2JES() {
  return _hJES;
}

TH2F* Analysis::GetH2JESError() {
  return _hJESError;
}

TH2F* Analysis::GetH2MassConstJES() {
  return _hMassConstJES;
}

TH2F* Analysis::GetH2MassConstJESError() {
  return _hMassConstJESError;
}

TH2F* Analysis::GetH2FSig() {
  return _hFSig;
}

TH2F* Analysis::GetH2FSigError() {
  return _hFSigError;
}

TH2F* Analysis::GetH2MassfSig() {
  return _hMassfSig;
}

TH2F* Analysis::GetH2MassfSigError() {
  return _hMassfSigError;
}

TH2F* Analysis::GetH2JESfSig() {
  return _hJESfSig;
}

TH2F* Analysis::GetH2JESfSigError() {
  return _hJESfSigError;
}

TString Analysis::GetIdentifier() {
  return _fIdentifier;
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
  double pt, eta, eff, effyScale_pt;

  std::vector<double> oneMinusBEffies(0) , oneMinusBEffies_scaled(0) ;
  std::vector<double> oneMinusBMistags(0), oneMinusBMistags_scaled(0);
  for(int i = 0; i < Njet; ++i){
    pt  = ((TLorentzVector*)jets->At(i))->Pt();
    eta = ((TLorentzVector*)jets->At(i))->Eta();

    if(pt > 670.)
      effyScale_pt    = 0.901615*((1.+(0.552628*670.))/(1.+(0.547195*670.)));
    if(pt < 30.)
      effyScale_pt    = 0.901615*((1.+(0.552628*30.))/(1.+(0.547195*30.)));
    else
      effyScale_pt    = 0.901615*((1.+(0.552628*pt))/(1.+(0.547195*pt)));

    if(pdgId[i] == 5 || pdgId[i] == -5){
      eff = _bTagEff->GetBinContent(_bTagEff->FindBin(pt,std::abs(eta)));
      oneMinusBEffies       .push_back(1.- eff);
      oneMinusBEffies_scaled.push_back(1.-(eff*effyScale_pt));
    }
    else if(pdgId[i] == 4 || pdgId[i] == -4){
      eff = _cTagEff->GetBinContent(_cTagEff->FindBin(pt,std::abs(eta)));
      oneMinusBMistags       .push_back(1.- eff);
      oneMinusBMistags_scaled.push_back(1.-(eff*effyScale_pt));
    }
    else{
      eff = _lTagEff->GetBinContent(_lTagEff->FindBin(pt,std::abs(eta)));
      oneMinusBMistags       .push_back(1.- eff);
      oneMinusBMistags_scaled.push_back(1.-(eff*(((0.948463+(0.00288102*pt))+(-7.98091e-06*(pt*pt)))+(5.50157e-09*(pt*(pt*pt)))) ));
    }
  }
  bTaggingEfficiency        = eventBTagProbability_(oneMinusBEffies       , oneMinusBMistags       );
  bTaggingEfficiency_scaled = eventBTagProbability_(oneMinusBEffies_scaled, oneMinusBMistags_scaled);

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

xml::XMLDocument* Analysis::_config(0);

void
Analysis::ReadConfigFromXMLFile(){
  if(!_config){
    _config = new xml::XMLDocument();
    TString xmlFilePath = "/afs/naf.desy.de/group/cms/scratch/eschliec/TopMass_hg_devel/Analyzer/Configuration_alljets.xml";
    int errorID = _config->LoadFile(xmlFilePath);
    if(errorID) {
      std::cout << "Parsing of XML file (" << xmlFilePath << ") failed with error " << errorID << "!" << std::endl;
      assert(!errorID);
    }
    xml::XMLNode *analysisConfiguration = 0;
    xml::XMLElement *configParameter = 0;

    analysisConfiguration = _config->FirstChild();
    while(!TString(analysisConfiguration->Value()).EqualTo("analysisConfig")){
      analysisConfiguration = analysisConfiguration->NextSibling();
    }
    if(analysisConfiguration->NoChildren()){
      std::cout << "No configuration contained in *" << analysisConfiguration->Value() << "* object in XMLFile:\n" << xmlFilePath << std::endl;
      assert(0);
    }
    do{
      if(!configParameter) configParameter = analysisConfiguration->FirstChildElement();
      else       configParameter = configParameter->NextSiblingElement();
      std::cout << "child: " << configParameter->Value() << std::endl;
      if(TString(configParameter->Value()).EqualTo("selection")) _selection = configParameter->GetText();
      if(TString(configParameter->Value()).EqualTo("treeVariables")){
        xml::XMLElement *variable = 0;
        do{
          if(!variable) variable = configParameter->FirstChildElement();
          else          variable = variable->NextSiblingElement();

          TString variableType = variable->Value();

          //std::cout << variableType << ": " << variable->GetText() << std::endl;

          /*
          if     (variableType.EqualTo("short" )) _shortVariables [variable->GetText()] = -1;
          else if(variableType.EqualTo("int"   )) _intVariables   [variable->GetText()] = -1;
          else if(variableType.EqualTo("float" )) _floatVariables [variable->GetText()] = -1.;
          else if(variableType.EqualTo("double")) _doubleVariables[variable->GetText()] = -1.;
          else if(variableType.EqualTo("uint"  )) _uintVariables  [variable->GetText()] =  0;
          else if(variableType.EqualTo("array")){
            xml::XMLElement *ushortArray = variable->FirstChildElement("ushort");
            xml::XMLElement *shortArray  = variable->FirstChildElement("short");
            xml::XMLElement *doubleArray = variable->FirstChildElement("double");
            int arrayLength = atoi(variable->FirstChildElement("length")->GetText());
            if     (ushortArray) _ushortArrayVariables[ushortArray->GetText()] = new unsigned short[arrayLength];
            else if(shortArray ) _shortArrayVariables [ shortArray->GetText()] = new short[arrayLength];
            else if(doubleArray) _doubleArrayVariables[doubleArray->GetText()] = new double[arrayLength];
          }
          else if(variableType.EqualTo("TClonesArray")){
            xml::XMLElement *lorentzVectorArray = variable->FirstChildElement("TLorentzVector");
            if(lorentzVectorArray) _TClonesArrayVariables[lorentzVectorArray->GetText()] = new TClonesArray("TLorentzVector");
          }
          */

          if(variableType.EqualTo("short") || variableType.EqualTo("int") ||
             variableType.EqualTo("float") || variableType.EqualTo("double")) _variables[variableType][variable->GetText()] = -1.;
          else if(variableType.EqualTo("uint")) _variables[variableType][variable->GetText()] =  0;
          else if(variableType.EqualTo("array")){
            xml::XMLElement *ushortArray = variable->FirstChildElement("ushort");
            xml::XMLElement *shortArray  = variable->FirstChildElement("short");
            xml::XMLElement *doubleArray = variable->FirstChildElement("double");
            TString length = variable->FirstChildElement("length")->GetText();
            int arrayLength = atoi(variable->FirstChildElement("maxlength")->GetText());
            if     (ushortArray) _variables[variableType+TString("_ushort_")+length][ushortArray->GetText()] = new unsigned short[arrayLength];
            else if(shortArray ) _variables[variableType+TString("_short_" )+length][ shortArray->GetText()] = new short[arrayLength];
            else if(doubleArray) _variables[variableType+TString("_double_")+length][doubleArray->GetText()] = new double[arrayLength];
          }
          else if(variableType.EqualTo("TClonesArray")){
            xml::XMLElement *lorentzVectorArray = variable->FirstChildElement("TLorentzVector");
            if(lorentzVectorArray) _variables[variableType+TString("_TLorentzVector")][lorentzVectorArray->GetText()] = new TClonesArray("TLorentzVector");
          }

        }
        while(variable != configParameter->LastChildElement());
      }
    }
    while(configParameter != analysisConfiguration->LastChild());

    std::cout << "_selection = " << _selection << std::endl;

    /*
    for(std::map<TString, short>::const_iterator it = _shortVariables.begin(); it != _shortVariables.end(); ++it){
      std::cout << "short: " << it->first << " = " << it->second << std::endl;
    }
    for(std::map<TString, int>::const_iterator it = _intVariables.begin(); it != _intVariables.end(); ++it){
      std::cout << "int: " << it->first << " = " << it->second << std::endl;
    }
    for(std::map<TString, float>::const_iterator it = _floatVariables.begin(); it != _floatVariables.end(); ++it){
      std::cout << "float: " << it->first << " = " << it->second << std::endl;
    }
    for(std::map<TString, double>::const_iterator it = _doubleVariables.begin(); it != _doubleVariables.end(); ++it){
      std::cout << "double: " << it->first << " = " << it->second << std::endl;
    }
    for(std::map<TString, unsigned int>::const_iterator it = _uintVariables.begin(); it != _uintVariables.end(); ++it){
      std::cout << "uint: " << it->first << " = " << it->second << std::endl;
    }
    for(std::map<TString, unsigned short*>::const_iterator it = _ushortArrayVariables.begin(); it != _ushortArrayVariables.end(); ++it){
      std::cout << "ushort: " << it->first << " <- Array" << std::endl;
    }
    for(std::map<TString, short*>::const_iterator it = _shortArrayVariables.begin(); it != _shortArrayVariables.end(); ++it){
      std::cout << "short: " << it->first << " <- Array" << std::endl;
    }
    for(std::map<TString, double*>::const_iterator it = _doubleArrayVariables.begin(); it != _doubleArrayVariables.end(); ++it){
      std::cout << "double: " << it->first << " <- Array" << std::endl;
    }
    for(std::map<TString, TClonesArray*>::const_iterator it = _TClonesArrayVariables.begin(); it != _TClonesArrayVariables.end(); ++it){
      std::cout << "TClonesArray: " << it->first << " <- Array"<< std::endl;
    }
    */

    /*
    TTree * testTree = new TTree("testTree","testTree");
    for(std::map<TString, std::map<TString, variableTypes> >::iterator it1 = _variables.begin(); it1 != _variables.end(); ++it1){
      for(std::map<TString, variableTypes>::iterator it2 = it1->second.begin(); it2 != it1->second.end(); ++it2){
        std::cout << it1->first <<  ": " << it2->first << " = " << it2->second << std::endl;
        TBranch* br = 0;
        if     (it1->first.EqualTo("double")){ br = testTree->Branch(it2->first, &it2->second, it2->first+TString("/D")); }
        else if(it1->first.EqualTo("float" )){ br = testTree->Branch(it2->first, &it2->second, it2->first+TString("/F")); }
        else if(it1->first.EqualTo("int"   )){ br = testTree->Branch(it2->first, &it2->second, it2->first+TString("/I")); }
        else if(it1->first.EqualTo("uint"  )){ br = testTree->Branch(it2->first, &it2->second, it2->first+TString("/i")); }
        else if(it1->first.EqualTo("short" )){ br = testTree->Branch(it2->first, &it2->second, it2->first+TString("/S")); }
        else if(it1->first.Contains("array")){
          TString length = it1->first(it1->first.Last('_')+1, it1->first.Length());
          TString branchName = it2->first+TString("[")+length+TString("]");
          std::cout << branchName << " " << length << std::endl;
          if(it1->first.Contains("ushort")){
            br = testTree->Branch(it2->first, &it2->second, branchName+TString("/s"));
          }
          else if(it1->first.Contains("short")){
            br = testTree->Branch(it2->first, &it2->second, branchName+TString("/S"));
          }
          else if(it1->first.Contains("double")){
            br = testTree->Branch(it2->first, &it2->second, branchName+TString("/D"));
          }
        }
        else if(it1->first.Contains("TClonesArray")){
          if(it1->first.Contains("TLorentzVector")){
            std::cout << it1->first << std::endl;
          }
        }
        if(br) std::cout << br->GetTitle() << std::endl;
      }
    }
*/

  }
}

void
Analysis::SetBranchStatuses(TTree* tree){
  tree->SetBranchStatus("*", 0);
  for(std::map<TString, std::map<TString, variableTypes> >::iterator it1 = _variables.begin(); it1 != _variables.end(); ++it1){
    for(std::map<TString, variableTypes>::iterator it2 = it1->second.begin(); it2 != it1->second.end(); ++it2){
      std::cout << it1->first <<  ": " << it2->first << " = " << it2->second << std::endl;
      tree->SetBranchStatus(it2->first, 1);
    }
  }
}

void
Analysis::SetBranchAddresses(TTree* tree){

  for(std::map<TString, std::map<TString, variableTypes> >::iterator it1 = _variables.begin(); it1 != _variables.end(); ++it1){
    for(std::map<TString, variableTypes>::iterator it2 = it1->second.begin(); it2 != it1->second.end(); ++it2){
      if     (it1->first.EqualTo("double")){ tree->SetBranchAddress(it2->first, &it2->second); }
      else if(it1->first.EqualTo("float" )){ tree->SetBranchAddress(it2->first, &it2->second); }
      else if(it1->first.EqualTo("int"   )){ tree->SetBranchAddress(it2->first, &it2->second); }
      else if(it1->first.EqualTo("uint"  )){ tree->SetBranchAddress(it2->first, &it2->second); }
      else if(it1->first.EqualTo("short" )){ tree->SetBranchAddress(it2->first, &it2->second); }
      else if(it1->first.Contains("array")){
        if(it1->first.Contains("ushort")){
          tree->SetBranchAddress(it2->first, boost::get<unsigned short*>(it2->second));
        }
        else if(it1->first.Contains("short")){
          tree->SetBranchAddress(it2->first, boost::get<short*>(it2->second));
        }
        else if(it1->first.Contains("double")){
          tree->SetBranchAddress(it2->first, boost::get<double*>(it2->second));
        }
      }
      else if(it1->first.Contains("TClonesArray")){
        if(it1->first.Contains("TLorentzVector")){
          tree->SetBranchAddress(it2->first, &it2->second);
        }
      }
    }
  }

}
