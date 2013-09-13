/*
 * RandomSubsetCreatorAllJets.cxx
 *
 *  Created on: Oct 23, 2012
 *      Author: eschliec
 */

#include "RandomSubsetCreatorAllJets.h"

#include "Analysis.h"
#include "ProgramOptionsReader.h"

#include <iostream>

#include "TChain.h"
#include "TEventList.h"
#include "TLorentzVector.h"
#include "TRandom3.h"
#include "TROOT.h"
#include "TSystem.h"

//#include "LHAPDF/LHAPDF.h"

typedef ProgramOptionsReader po;

RandomSubsetCreatorAllJets::RandomSubsetCreatorAllJets() :
selection_  (po::GetOption<std::string>("analysisConfig.selection" )), // filled from program options
samplePath_ (po::GetOption<std::string>("analysisConfig.samplePath")), // filled from program options
fIdentifier_(po::GetOption<std::string>("input")), // filled from program options
fFile_(samplePath_+fIdentifier_+TString(".root")),
tmpFile_(0), // has to survive until destructor is called
bTagEff_(0), // filled in FetchBTagEfficiencyHistograms
cTagEff_(0), // filled in FetchBTagEfficiencyHistograms
lTagEff_(0)  // filled in FetchBTagEfficiencyHistograms
{
  FetchBTagEfficiencyHistograms();
}

RandomSubsetCreatorAllJets::~RandomSubsetCreatorAllJets() {
  delete bTagEff_;
  delete cTagEff_;
  delete lTagEff_;
  tmpFile_->Close();
  delete tmpFile_;
}

void
RandomSubsetCreatorAllJets::FetchBTagEfficiencyHistograms() {
  TFile * bTagFile = TFile::Open(samplePath_+TString("bTagFile.root"));
  gROOT->cd();
  bTagEff_ = (TH2F*)bTagFile->Get("histb")->Clone();
  cTagEff_ = (TH2F*)bTagFile->Get("histc")->Clone();
  lTagEff_ = (TH2F*)bTagFile->Get("histl")->Clone();
  bTagFile->Close();

  delete bTagFile;
}

TTree*
RandomSubsetCreatorAllJets::CreateRandomSubset() {
  std::cout << "Create random subset..." << std::endl;

  TChain* tmpChainSig = new TChain("FullHadTreeWriter/tree");
  tmpChainSig->Add(fFile_);

  tmpChainSig->SetBranchStatus("*",0);
  tmpChainSig->SetBranchStatus("nCombos", 1);
  tmpChainSig->SetBranchStatus("comboTypes", 1);
  tmpChainSig->SetBranchStatus("fitAssigns", 1);
  tmpChainSig->SetBranchStatus("topMasses", 1);
  tmpChainSig->SetBranchStatus("topMass", 1);
  tmpChainSig->SetBranchStatus("w1Mass", 1);
  tmpChainSig->SetBranchStatus("w2Mass", 1);
  tmpChainSig->SetBranchStatus("probs", 1);
  tmpChainSig->SetBranchStatus("prob", 1);
  tmpChainSig->SetBranchStatus("dRbb", 1);

  tmpChainSig->SetBranchStatus("Njet", 1);
  tmpChainSig->SetBranchStatus("jets", 1);
  //_fChain->SetBranchStatus("bTag_CSV", 1);
  tmpChainSig->SetBranchStatus("partonFlavour", 1);

  tmpChainSig->SetBranchStatus("runNumber", 1);
  tmpChainSig->SetBranchStatus("luminosityBlockNumber", 1);
  tmpChainSig->SetBranchStatus("eventNumber", 1);

  tmpChainSig->SetBranchStatus("nPU", 1);
  tmpChainSig->SetBranchStatus("nPUTru", 1);

  tmpChainSig->SetBranchStatus("MCweight", 1);

  if(fFile_.Contains("PDF")){
    tmpChainSig->SetBranchStatus("id1", 1);
    tmpChainSig->SetBranchStatus("id2", 1);
    tmpChainSig->SetBranchStatus("x1", 1);
    tmpChainSig->SetBranchStatus("x2", 1);
    tmpChainSig->SetBranchStatus("Q", 1);
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
  short * comboTypes = new short[kMAXCombos];
  short * pdgId = new short[kMAXNJets];
  short * fitAssigns = new short[6];
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
  if(fIdentifier_.Contains("jet4_02")) triggerUncertainty += " & jets[3].Pt() > 62.";
  if(fIdentifier_.Contains("jet4_05")) triggerUncertainty += " & jets[3].Pt() > 65.";
  if(fIdentifier_.Contains("jet4_10")) triggerUncertainty += " & jets[3].Pt() > 70.";
  if(fIdentifier_.Contains("jet5_02")) triggerUncertainty += " & jets[4].Pt() > 52.";
  if(fIdentifier_.Contains("jet5_05")) triggerUncertainty += " & jets[4].Pt() > 55.";
  if(fIdentifier_.Contains("jet5_10")) triggerUncertainty += " & jets[4].Pt() > 60.";
  if(fIdentifier_.Contains("jet6_02")) triggerUncertainty += " & jets[5].Pt() > 42.";
  if(fIdentifier_.Contains("jet6_05")) triggerUncertainty += " & jets[5].Pt() > 45.";
  if(fIdentifier_.Contains("jet6_10")) triggerUncertainty += " & jets[5].Pt() > 50.";

  TTree *tmpTreeSig = 0, *tmpTreeBkg = 0, *returnedTree = 0;

  // check existence of a temp directory and create one if not available
  TString tmpFilePath(gSystem->Getenv("TMPDIR"));
  if(tmpFilePath.IsNull() || !tmpFilePath.Length()){
    tmpFilePath = gSystem->GetFromPipe("mktemp -d");
    gSystem->Setenv("TMPDIR", tmpFilePath);
  }
  std::cout << "Directory to be used for temporary files: " << tmpFilePath << std::endl;

  double lumi = po::GetOption<double>("lumi");
  if (lumi!=0) {
    tmpFile_ = TFile::Open(tmpFilePath+TString("/tempTree.root"), "READ");
    if(tmpFile_ && !tmpFile_->IsZombie()){
      tmpTreeSig = (TTree*)tmpFile_->Get(TString("fullTree_")+fIdentifier_);
      tmpTreeBkg = (TTree*)tmpFile_->Get("fullTreeBkg");
      std::cout << "Reading: " << tmpFile_->GetName() << std::endl;
    }
    else{
      tmpFile_ = TFile::Open(tmpFilePath+TString("/tempTree.root"), "RECREATE", "", 0);
      std::cout << "Creating: " << tmpFile_->GetName() << std::endl;
    }
    if(!tmpTreeSig || !tmpTreeBkg){
      tmpFile_->ReOpen("UPDATE");
      std::cout << "Updating: " << tmpFile_->GetName() << std::endl;

      if(!tmpTreeSig){
        tmpChainSig->Draw(">>selectedEvents",selection_+triggerUncertainty,"goff");
        TEventList* selectedEvents = (TEventList*)gDirectory->Get("selectedEvents");

        tmpChainSig->SetBranchAddress("nCombos", &nCombos);
        tmpChainSig->SetBranchAddress("comboTypes", comboTypes);
        tmpChainSig->SetBranchAddress("fitAssigns", fitAssigns);
        tmpChainSig->SetBranchAddress("topMasses", topMasses);
        tmpChainSig->SetBranchAddress("topMass", &topMass);
        tmpChainSig->SetBranchAddress("w1Mass", w1Masses);
        tmpChainSig->SetBranchAddress("w2Mass", w2Masses);
        tmpChainSig->SetBranchAddress("probs", probs);
        tmpChainSig->SetBranchAddress("prob", &prob);
        tmpChainSig->SetBranchAddress("dRbb", &dRbb);

        tmpChainSig->SetBranchAddress("Njet", &Njet);
        tmpChainSig->SetBranchAddress("jets", &jets);
        //fChain->SetBranchAddress("bTag_CSV", bTag);
        tmpChainSig->SetBranchAddress("partonFlavour", pdgId);

        tmpChainSig->SetBranchAddress("runNumber", &runNumber);
        tmpChainSig->SetBranchAddress("luminosityBlockNumber", &luminosityBlockNumber);
        tmpChainSig->SetBranchAddress("eventNumber", &eventNumber);

        tmpChainSig->SetBranchAddress("dRbb", &dRbb);

        tmpChainSig->SetBranchAddress("nPU", &nPU);
        tmpChainSig->SetBranchAddress("nPUTru", &nPUTru);

        tmpChainSig->SetBranchAddress("MCweight", &MCweight);

        if(fFile_.Contains("PDF")){
          tmpChainSig->SetBranchAddress("id1", &id1);
          tmpChainSig->SetBranchAddress("id2", &id2);
          tmpChainSig->SetBranchAddress("x1", &x1);
          tmpChainSig->SetBranchAddress("x2", &x2);
          tmpChainSig->SetBranchAddress("Q", &Q);
        }

        tmpTreeSig = tmpChainSig->CloneTree(0);
        tmpTreeSig->SetName(TString("fullTree_")+fIdentifier_);

        for(int idx = 0 , l = selectedEvents->GetN(); idx < l; ++idx){
          tmpChainSig->GetEntry(selectedEvents->GetEntry(idx));
          nCombos = 1;
          tmpTreeSig->Fill();
        }
        AddWeights(tmpTreeSig);
        tmpTreeSig->Write(0,TObject::kOverwrite);
        delete tmpTreeSig;
      }

      if(!tmpTreeBkg){
        TFile * fileBkg = TFile::Open(samplePath_+TString("QCDEstimationMix_2011_skimmed2.root"));
        if(TString(fileBkg->GetName()).Contains("_skimmed")) tmpTreeBkg = (TTree*)fileBkg->Get("tree");
        else  tmpTreeBkg = (TTree*)fileBkg->Get("analyzeFullHadEventMixer/tree");

        tmpTreeBkg->SetBranchStatus("*",0);
        tmpTreeBkg->SetBranchStatus("topMasses", 1);
        tmpTreeBkg->SetBranchStatus("topMass", 1);
        tmpTreeBkg->SetBranchStatus("w1Mass", 1);
        tmpTreeBkg->SetBranchStatus("w2Mass", 1);
        tmpTreeBkg->SetBranchStatus("probs", 1);
        tmpTreeBkg->SetBranchStatus("prob", 1);
        tmpTreeBkg->SetBranchStatus("dRbb", 1);

        tmpTreeBkg->SetBranchStatus("jets", 1);

        tmpTreeBkg->SetBranchStatus("runNumber", 1);
        tmpTreeBkg->SetBranchStatus("luminosityBlockNumber", 1);
        tmpTreeBkg->SetBranchStatus("eventNumber", 1);

        tmpFile_->cd();

        TString selBkg = selection_; selBkg.ReplaceAll("dRbb","dRbb[0]");
        TTree* tmpTreeBkgHelper = tmpTreeBkg->CopyTree(selBkg);
        tmpTreeBkgHelper->SetName("fullTreeBkg");
        AddWeights(tmpTreeBkgHelper,true);
        tmpTreeBkgHelper->Write(0,TObject::kOverwrite);
        delete tmpTreeBkgHelper;

        delete tmpTreeBkg;
        fileBkg->Close();
        delete fileBkg;
      }

      tmpFile_->Close();
      tmpFile_ = TFile::Open(tmpFilePath+TString("/tempTree.root"),"READ");
      tmpTreeSig = (TTree*)tmpFile_->Get(TString("fullTree_")+fIdentifier_);
      tmpTreeBkg = (TTree*)tmpFile_->Get("fullTreeBkg");
    }

    TRandom3* myRandom = new TRandom3(0);
    std::cout << "Random seed: " << myRandom->GetSeed() << std::endl;

    tmpTreeSig->SetCacheSize(10000000000);
    tmpTreeBkg->SetCacheSize(10000000000);

    tmpTreeSig->SetBranchStatus("*",0);
    tmpTreeSig->SetBranchStatus("CombinedWeight",1);
    tmpTreeSig->SetBranchStatus("meanWMass",1);
    tmpTreeSig->SetBranchStatus("topMasses",1);
    tmpTreeSig->SetBranchStatus("topMass",1);
    tmpTreeSig->SetBranchStatus("w1Mass",1);
    tmpTreeSig->SetBranchStatus("w2Mass",1);
    tmpTreeSig->SetBranchStatus("probs",1);
    tmpTreeSig->SetBranchStatus("prob",1);
    tmpTreeSig->SetBranchStatus("dRbb",1);
    tmpTreeSig->SetBranchStatus("jets",1);
    tmpTreeSig->SetBranchStatus("nCombos",1);
    tmpTreeSig->SetBranchStatus("comboTypes",1);
    tmpTreeSig->SetBranchStatus("fitAssigns",1);
    tmpTreeSig->SetBranchStatus("runNumber", 1);
    tmpTreeSig->SetBranchStatus("luminosityBlockNumber", 1);
    tmpTreeSig->SetBranchStatus("eventNumber", 1);

    tmpTreeSig->SetBranchAddress("CombinedWeight",&CombinedWeight);
    tmpTreeSig->SetBranchAddress("meanWMass",&meanWMass);
    tmpTreeSig->SetBranchAddress("topMasses",topMasses);
    tmpTreeSig->SetBranchAddress("topMass",&topMass);
    tmpTreeSig->SetBranchAddress("w1Mass",w1Masses);
    tmpTreeSig->SetBranchAddress("w2Mass",w2Masses);
    tmpTreeSig->SetBranchAddress("probs",probs);
    tmpTreeSig->SetBranchAddress("prob",&prob);
    tmpTreeSig->SetBranchAddress("dRbb",&dRbb);
    tmpTreeSig->GetBranch("jets")->SetAutoDelete(kFALSE);
    tmpTreeSig->SetBranchAddress("jets", &jets);
    tmpTreeSig->SetBranchAddress("nCombos",&nCombos);
    tmpTreeSig->SetBranchAddress("comboTypes",comboTypes);
    tmpTreeSig->SetBranchAddress("fitAssigns",fitAssigns);
    tmpTreeSig->SetBranchAddress("runNumber", &runNumber);
    tmpTreeSig->SetBranchAddress("luminosityBlockNumber", &luminosityBlockNumber);
    tmpTreeSig->SetBranchAddress("eventNumber", &eventNumber);

    tmpTreeBkg->SetBranchAddress("meanWMass",&meanWMass);
    tmpTreeBkg->SetBranchAddress("topMasses", topMasses);
    tmpTreeBkg->SetBranchAddress("topMass", &topMass);
    tmpTreeBkg->SetBranchAddress("w1Mass", w1Masses);
    tmpTreeBkg->SetBranchAddress("w2Mass", w2Masses);
    tmpTreeBkg->SetBranchAddress("probs", probs);
    tmpTreeBkg->SetBranchAddress("prob", &prob);
    tmpTreeBkg->SetBranchAddress("dRbb", &dRbb);
    tmpTreeBkg->GetBranch("jets")->SetAutoDelete(kFALSE);
    tmpTreeBkg->SetBranchAddress("jets", &jets);
    tmpTreeBkg->SetBranchAddress("nCombos",&nCombos);

    tmpTreeBkg->SetBranchAddress("runNumber", &runNumber);
    tmpTreeBkg->SetBranchAddress("luminosityBlockNumber", &luminosityBlockNumber);
    tmpTreeBkg->SetBranchAddress("eventNumber", &eventNumber);

    double maxMCWeight = tmpTreeSig->GetMaximum("CombinedWeight");

    if (maxMCWeight == -1) { std::cout << "Running over data?" << std::endl; }

    double fSig = po::GetOption<double>("fsig"); // 0.504; // 0.539; //
    if(fFile_.Contains("fSig_Up"))
      fSig += 0.10;
    else if(fFile_.Contains("fSig_Down"))
      fSig -= 0.10;
    else if(fFile_.Contains("BJES_Up"))
      fSig *= 1.0124;
    else if(fFile_.Contains("BJES_Down"))
      fSig *= 0.9895;
    else if(fFile_.Contains("JES_Up"))
      fSig *= 1.1042;
    else if(fFile_.Contains("JES_Down"))
      fSig *= 0.8967;
    else if(fFile_.Contains("JER_Up"))
      fSig *= 0.9816;
    else if(fFile_.Contains("JER_Down"))
      fSig *= 1.0192;
    else if(fFile_.Contains("BTAG_Up"))
      fSig *= 1.0509;
    else if(fFile_.Contains("BTAG_Down"))
      fSig *= 0.9565;
    else if(fFile_.Contains("Scale_Up"))
      fSig *= 0.8726;
    else if(fFile_.Contains("Scale_Down"))
      fSig *= 1.1184;
    else if(fFile_.Contains("Matching_Up"))
      fSig *= 0.9718;
    else if(fFile_.Contains("P11_NoCR"))
      fSig *= 1.0294;
    else if(fFile_.Contains("P11mpiHi"))
      fSig *= 1.0202;
    else if(fFile_.Contains("P11TeV"))
      fSig *= 1.0234;
    else if(fFile_.Contains("jet4_02"))
      fSig *= 0.9256;
    else if(fFile_.Contains("jet5_02"))
      fSig *= 0.9431;
    else if(fFile_.Contains("jet6_02"))
      fSig *= 0.9305;

    int permsMC  = tmpTreeSig->GetEntries();
    int permsBkg = tmpTreeBkg->GetEntries();
    int eventsPE = myRandom->Poisson(2767./3544.844*lumi); // add poisson
    int eventsDrawn = 0;

    tmpFile_->ReOpen("UPDATE");
    returnedTree = tmpTreeSig->CloneTree(0);
    returnedTree->SetName("tree");
    // remove unnecessary branches from fTree
    returnedTree->SetBranchStatus("CombinedWeight",0);
    returnedTree->SetBranchStatus("jets",0);
    //fTree->SetBranchStatus("runNumber", 0);
    //fTree->SetBranchStatus("luminosityBlockNumber", 0);
    //fTree->SetBranchStatus("eventNumber", 0);

    //double& CombinedWeight = GetVariable("CombinedWeight");

    int signalDrawn = 0, backgroundDrawn = 0;
    if(lumi>0) {
      //int drawCounter = 0;
      while (eventsDrawn < eventsPE) {
        double fSigRndm = myRandom->Uniform(0.,1.);
        if(fSigRndm < fSig){
          int oldEventsDrawn = eventsDrawn;
          do {
            int drawn = myRandom->Integer(permsMC);
            //++drawCounter;
            tmpTreeSig->GetEntry(drawn);
            //std::cout << CombinedWeight << " " << maxMCWeight << std::endl;
            if (CombinedWeight > myRandom->Uniform(0., maxMCWeight)) {
              // reduce size of tree before filling
              nCombos = 1;
              //if(fTree->Fill() == -1) std::cout << eventsDrawn << " " << meanWMass << " " << topMasses[0] << " " << topMass << " " << w1Masses[0] << " " << w2Masses[0] << " " << probs[0] << " " << prob << " " << dRbb << " " << nCombos << " " << comboTypes[0] << " " << runNumber << " " << luminosityBlockNumber << " " << eventNumber << " " << std::endl;
              returnedTree->Fill();

              ++eventsDrawn;
              ++signalDrawn;
            }
          }
          while(oldEventsDrawn == eventsDrawn);
        }
        else{
          int drawn = myRandom->Integer(permsBkg);
          //++drawCounter;
          tmpTreeBkg->GetEntry(drawn);
          // reduce size of tree before filling
          nCombos = 1;
          // set missing variable for background tree
          comboTypes[0] = 0;
          //if(fTree->Fill() == -1) std::cout << eventsDrawn << " " << meanWMass << " " << topMasses[0] << " " << topMass << " " << w1Masses[0] << " " << w2Masses[0] << " " << probs[0] << " " << prob << " " << dRbb << " " << nCombos << " " << comboTypes[0] << " " << runNumber << " " << luminosityBlockNumber << " " << eventNumber << " " << std::endl;
          returnedTree->Fill();
          ++eventsDrawn;
          ++backgroundDrawn;
        }
      }
      //std::cout << "DRAWCOUNTER: " << drawCounter << std::endl;
    }
    else {
      for(int ev = 0; ev < permsMC; ++ev){
        tmpTreeSig->GetEntry(ev);
        returnedTree->Fill();
        ++signalDrawn;
      }
      std::cout << "wanted / available events: " << permsMC*((1.-fSig)/fSig) << " / " << permsBkg << std::endl;
      for(int ev = 0; ev < permsMC*((1.-fSig)/fSig); ++ev){
        int drawn = myRandom->Integer(permsBkg);
        tmpTreeBkg->GetEntry(drawn);
        returnedTree->Fill();
        ++backgroundDrawn;
      }
    }
    returnedTree->Write(0,TObject::kOverwrite);
    std::cout << "Events drawn: " << eventsDrawn << " (sig: " << signalDrawn << ", bkg: " << backgroundDrawn << ") -> fSig: " << double(signalDrawn)/double(signalDrawn+backgroundDrawn) << " (default: " << fSig << ")" << std::endl;
    delete returnedTree;
    delete tmpTreeSig;
    delete tmpTreeBkg;
    tmpFile_->Close();
    tmpFile_ = TFile::Open(tmpFilePath+TString("/tempTree.root"));
    returnedTree = (TTree*)tmpFile_->Get("tree");
    delete myRandom;
  }
  else {
    TFile* myFile = TFile::Open(fFile_, "READ");
    TTree* myTree = (TTree*)myFile->Get("analyzeKinFit/eventTree");

    tmpFile_ = new TFile(tmpFilePath+TString("/tempTree.root"), "RECREATE");
    tmpFile_->cd();
    tmpTreeSig = tmpChainSig->CopyTree(selection_);
    AddWeights(tmpTreeSig,true);

    myTree = myTree->CopyTree("fitProb[0] > 0.09 && fitTop1[0].M() > 100.0 && fitTop2[0].M() < 550.0");
    myTree->SetName("eventTreeNew");

    tmpTreeSig->AddFriend("eventTreeNew");

    tmpFile_->Write(0,TObject::kOverwrite);
    delete tmpTreeSig;
    tmpFile_->Close();
    tmpFile_ = TFile::Open(tmpFilePath+TString("/tempTree.root"));
    returnedTree = (TTree*)tmpFile_->Get("tree");
  }
  delete tmpChainSig;
  delete[] topMasses;
  delete[] w1Masses;
  delete[] w2Masses;
  delete[] probs;
  delete[] comboTypes;
  delete[] fitAssigns;
  delete[] pdgId;
  jets->Delete();
  std::cout << "Created random subset." << std::endl;
  return returnedTree;
}

void
RandomSubsetCreatorAllJets::AddWeights(TTree* tempTree, bool isData) {
  std::cout << "Adapting " << tempTree->GetName();
  if(!isData) std::cout << " and adding combined weight";
  std::cout << " ..." << std::endl;

  enumForPUWeights whichSample = kFall10;
  if(!isData){
    if     (fFile_.Contains("S11"        )) whichSample = kSummer11;
    else if(fFile_.Contains("FAST"       )) whichSample = kFall10;
    else if(fFile_.Contains("F11_PU_Up"  )) whichSample = kFall11Plus05;
    else if(fFile_.Contains("F11_PU_Down")) whichSample = kFall11Minus05;
    else if(fFile_.Contains("F11"        )) whichSample = kFall11;
  }
  int whichPDFUncertainty = -1;
  bool upVariation = true; // value does not matter as long as _whichPDFUncertainty_==-1
  if(!isData){
    if(fFile_.Contains("PDF")){
      TString pdfUncert = fFile_;
      upVariation = pdfUncert.Contains("Up") ? true : false;
      pdfUncert = TString(pdfUncert(pdfUncert.Index("PDF",TString::kIgnoreCase)+3,2).Data());
      whichPDFUncertainty = pdfUncert.Atoi();
    }
    std::cout << "whichPDFUncertainty: " << whichPDFUncertainty << " (" << fFile_ << ")" << std::endl;
  }

  const int kMAXNJets(50);
  TString bTagAlgo = "CSV";
  int Njet;
  if(!isData) tempTree->SetBranchStatus ("Njet", 1);
  if(!isData) tempTree->SetBranchAddress("Njet", &Njet);
  //float bTag[kMAXNJets];
  //if(!isData) tempTree->SetBranchStatus (TString("bTag_")+bTagAlgo, 1);
  //if(!isData) tempTree->SetBranchAddress(TString("bTag_")+bTagAlgo, &bTag );
  short pdgId[kMAXNJets];
  if(!isData) tempTree->SetBranchStatus ("partonFlavour", 1);       // pdgId  *or*  partonFlavour
  if(!isData) tempTree->SetBranchAddress("partonFlavour", &pdgId); //  pdgId  *or*  partonFlavour
  TClonesArray * jets = new TClonesArray("TLorentzVector");
  tempTree->SetBranchStatus("jets", 1);
  tempTree->GetBranch("jets")->SetAutoDelete(kFALSE);
  tempTree->SetBranchAddress("jets", &jets);
  //unsigned int nCombos  = 0;
  //tempTree->SetBranchStatus ("nCombos", 1);
  //tempTree->SetBranchAddress("nCombos", &nCombos);
  int kMAXCombos = 12000;
  double * w1Masses  = new double[kMAXCombos];
  double * w2Masses  = new double[kMAXCombos];
  tempTree->SetBranchStatus("w1Mass", 1);
  tempTree->SetBranchStatus("w2Mass", 1);
  tempTree->SetBranchAddress("w1Mass",w1Masses);
  tempTree->SetBranchAddress("w2Mass",w2Masses);
  short * fitAssigns = new short[6];
  tempTree->SetBranchStatus("fitAssigns", 1);
  tempTree->SetBranchAddress("fitAssigns",fitAssigns);

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
    ////if(samplePath_.EndsWith("/20/")){
    //  TLorentzVector lQ    = *((TLorentzVector*)jets->At(fitAssigns[0]));
    //  TLorentzVector lQBar = *((TLorentzVector*)jets->At(fitAssigns[1]));
    //  TLorentzVector lP    = *((TLorentzVector*)jets->At(fitAssigns[3]));
    //  TLorentzVector lPBar = *((TLorentzVector*)jets->At(fitAssigns[4]));
    //  //double lQPt    = lQ   .Pt();
    //  //double lQBarPt = lQBar.Pt();
    //  //double lPPt    = lP   .Pt();
    //  //double lPBarPt = lPBar.Pt();
    //  //lQ    = lQ    * L7PartonCorrection(lQPt);
    //  //lQBar = lQBar * L7PartonCorrection(lQBarPt);
    //  //lP    = lP    * L7PartonCorrection(lPPt);
    //  //lPBar = lPBar * L7PartonCorrection(lPBarPt);
    //  double mW1 = (lQ + lQBar).M();
    //  double mW2 = (lP + lPBar).M();
    //  meanWMass = (mW1+mW2)/2.;
    ////}
    ////else{
      meanWMass = (w1Masses[0]+w2Masses[0])/2.0;
    ////}
    br2->Fill();
  }
  jets->Delete();
  delete[] w1Masses;
  delete[] w2Masses;
  std::cout << "Adapted  " << tempTree->GetName();
  if(!isData) std::cout << " and added  combined weight";
  std::cout << " ..." << std::endl;
}

// return the PU weights for the different samples
double
RandomSubsetCreatorAllJets::calcPUWeight_(enum enumForPUWeights sample, short nPU)
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

double
RandomSubsetCreatorAllJets::calcPDFWeight_(int whichPDFUncertainty, bool upVariation, double x1, int id1, float Q, double x2, int id2)
{
  if(whichPDFUncertainty < 0)
    return 1.;

  //const std::string NAME = "cteq6mE";
  //LHAPDF::initPDFSetM(1, NAME, LHAPDF::LHGRID);
  //LHAPDF::initPDFSetM(2, NAME, LHAPDF::LHGRID);
  //
  //LHAPDF::initPDFM(1, 0);
  //if(upVariation) LHAPDF::initPDFM(2, 2*whichPDFUncertainty-1);
  //else            LHAPDF::initPDFM(2, 2*whichPDFUncertainty);
  ////std::cout << whichPDFUncertainty << ": " << LHAPDF::xfxM(2, x1, Q, id1) << " " << LHAPDF::xfxM(2, x2, Q, id2) << " " << LHAPDF::xfxM(1, x1, Q, id1) << " " << LHAPDF::xfxM(1, x2, Q, id2) << std::endl;
  //return (LHAPDF::xfxM(2, x1, Q, id1)*LHAPDF::xfxM(2, x2, Q, id2)/(LHAPDF::xfxM(1, x1, Q, id1)*LHAPDF::xfxM(1, x2, Q, id2)));
  return 1.0;
}

double
RandomSubsetCreatorAllJets::calcBTagWeight_(int Njet, short * pdgId, TClonesArray * jets)
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
      eff = bTagEff_->GetBinContent(bTagEff_->FindBin(pt,std::abs(eta)));
      oneMinusBEffies       .push_back(1.- eff);
      oneMinusBEffies_scaled.push_back(1.-(eff*effyScale_pt));
    }
    else if(pdgId[i] == 4 || pdgId[i] == -4){
      eff = cTagEff_->GetBinContent(cTagEff_->FindBin(pt,std::abs(eta)));
      oneMinusBMistags       .push_back(1.- eff);
      oneMinusBMistags_scaled.push_back(1.-(eff*effyScale_pt));
    }
    else{
      eff = lTagEff_->GetBinContent(lTagEff_->FindBin(pt,std::abs(eta)));
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
RandomSubsetCreatorAllJets::eventBTagProbability_(std::vector<double> &oneMinusBEffies, std::vector<double> &oneMinusBMistags, bool verbose){
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

double RandomSubsetCreatorAllJets::L7PartonCorrection(double pt)
{
  double mW = 80.4;
  double p0 = 90.9183;
  double p1 = 4.21458e-2;
  double p2 = -1.62001;
  double xi = 1.01988005;

  return ((mW*xi)/(p0+p1*pt+p2*log(pt)));
}
