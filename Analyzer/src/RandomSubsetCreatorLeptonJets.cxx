#include "RandomSubsetCreatorLeptonJets.h"

#include "Analysis.h"
#include "ProgramOptionsReader.h"

#include <iostream>
#include <boost/program_options.hpp>
#include <boost/algorithm/string/split.hpp>
#include <boost/algorithm/string/classification.hpp>
#include <boost/progress.hpp>

#include "TChain.h"
#include "TEventList.h"
#include "TLorentzVector.h"
#include "TRandom3.h"
#include "TROOT.h"
#include "TSystem.h"

typedef ProgramOptionsReader po;

RandomSubsetCreatorLeptonJets::RandomSubsetCreatorLeptonJets() :
selection_  (po::GetOption<std::string>("analysisConfig.selection" )), // filled from program options
samplePath_ (po::GetOption<std::string>("analysisConfig.samplePath")), // filled from program options
fIdentifier_(po::GetOption<std::string>("input")), // filled from program options
fChannel_   (po::GetOption<std::string>("channel")),
fLumi  (po::GetOption<double>("lumi")),
fSig   (po::GetOption<double>("fsig")),
fBDisc (po::GetOption<double>("bdisc")),
fWeight(po::GetOption<std::string>("weight"))
{
  if (!strcmp(fChannel_, "muon") || !strcmp(fChannel_, "all")) {
    fTreeTTmu = PrepareTree(samplePath_+fIdentifier_+TString("_muon/analyzeTop.root"));
    fTreeWmu  = PrepareTree("/scratch/hh/lustre/cms/user/mseidel/Fall11_Wbb_muon/analyzeTop.root");
    fTreeSTmu = PrepareTree("/scratch/hh/lustre/cms/user/mseidel/Fall11_T_muon/analyzeTop.root");
  }

  if (!strcmp(fChannel_, "electron") || !strcmp(fChannel_, "all")) {
    fTreeTTe  = PrepareTree(samplePath_+fIdentifier_+TString("_electron/analyzeTop.root"));
    fTreeWe   = PrepareTree("/scratch/hh/lustre/cms/user/mseidel/Fall11_Wbb_electron/analyzeTop.root");
    fTreeSTe  = PrepareTree("/scratch/hh/lustre/cms/user/mseidel/Fall11_T_electron/analyzeTop.root");
  }
}

RandomSubsetCreatorLeptonJets::~RandomSubsetCreatorLeptonJets() {
}

TTree* RandomSubsetCreatorLeptonJets::CreateRandomSubset() {
  if (fLumi>0) {
    std::cout << "Create random subset..." << std::endl;

    time_t start, end;
    time(&start);
    time(&end);

    fTree = new TTree("fTree", "fTree");

    fTree->Branch("target", &target, "target/I");
    fTree->Branch("run", &run, "run/I");
    fTree->Branch("luminosityBlock", &luminosityBlock, "luminosityBlock/I");
    fTree->Branch("event", &event, "event/I");
    fTree->Branch("combi", &combi, "combi/I");
    fTree->Branch("nVertex", &nVertex, "nVertex/I");
    fTree->Branch("leptonId", &leptonId, "leptonId/I");

    fTree->Branch("hadTopMass", &hadTopMass, "hadTopMass/D");
    fTree->Branch("hadWRawMass", &hadWRawMass, "hadWRawMass/D");
    fTree->Branch("leptonPt", &leptonPt, "leptonPt/D");
    fTree->Branch("leptonC", &leptonC, "leptonC/D");
    fTree->Branch("hitFitProb", &hitFitProb, "hitFitProb/D");
    fTree->Branch("deltaThetaHadWHadB", &deltaThetaHadWHadB, "deltaThetaHadWHadB/D");
    fTree->Branch("deltaThetaHadQHadQBar", &deltaThetaHadQHadQBar, "deltaThetaHadQHadQBar/D");
    fTree->Branch("PUWeight", &PUWeight, "PUWeight/D");
    fTree->Branch("PUWeightUp", &PUWeightUp, "PUWeightUp/D");
    fTree->Branch("PUWeightDown", &PUWeightDown, "PUWeightDown/D");
    fTree->Branch("muWeight", &muWeight, "muWeight/D");
    fTree->Branch("bWeight", &bWeight, "bWeight/D");
    fTree->Branch("bWeight_bTagSFUp", &bWeight_bTagSFUp, "bWeight_bTagSFUp/D");
    fTree->Branch("bWeight_bTagSFDown", &bWeight_bTagSFDown, "bWeight_bTagSFDown/D");
    fTree->Branch("bWeight_misTagSFUp", &bWeight_misTagSFUp, "bWeight_misTagSFUp/D");
    fTree->Branch("bWeight_misTagSFDown", &bWeight_misTagSFDown, "bWeight_misTagSFDown/D");
    fTree->Branch("MCWeight", &MCWeight, "MCWeight/D");
    fTree->Branch("mcWeight", &mcWeight, "mcWeight/D");
    fTree->Branch("hadQBCSV", &hadQBCSV, "hadQBCSV/D");
    fTree->Branch("hadQBarBCSV", &hadQBarBCSV, "hadQBarBCSV/D");
    fTree->Branch("hadBBCSV", &hadBBCSV, "hadBBCSV/D");
    fTree->Branch("lepBBCSV", &lepBBCSV, "lepBBCSV/D");
    fTree->Branch("lepBBCSV", &lepBBCSV, "lepBBCSV/D");
    fTree->Branch("jetsPt", &jetsPt, "jetsPt[4]/D");
    fTree->Branch("pdfWeights", &pdfWeights, "pdfWeights[44]/D");

    TRandom3* random = new TRandom3(0);

    // DATA
    // eventTree->GetEntries("leptonPt>30 & bottomCSVJetMultiplicity > 1 & hitFitProb>0.2 & combi==0")
    // (Long64_t)5144
    double nEventsDataMuon      = 2906.;
    double nEventsDataElectron  = 2268.;

    int eventsPEMuon      = random->Poisson(nEventsDataMuon/5000.*fLumi);
    int eventsPEElectron  = random->Poisson(nEventsDataElectron/5000.*fLumi);

    if (!strcmp(fChannel_, "muon") || !strcmp(fChannel_, "all")) {
      DrawEvents(fTreeTTmu, eventsPEMuon*fSig);
      DrawEvents(fTreeWmu,  eventsPEMuon*(1.-fSig)*1./4.);
      DrawEvents(fTreeSTmu, eventsPEMuon*(1.-fSig)*3./4.);
    }

    if (!strcmp(fChannel_, "electron") || !strcmp(fChannel_, "all")) {
      DrawEvents(fTreeTTe, eventsPEElectron*fSig);
      DrawEvents(fTreeWe,  eventsPEElectron*(1.-fSig)*1./4.);
      DrawEvents(fTreeSTe, eventsPEElectron*(1.-fSig)*3./4.);
    }

    time(&end);
    std::cout << "Created random subset in " << difftime(end, start) << " seconds." << std::endl;
  }
  else {
    TList treeList;

    if (!strcmp(fChannel_, "muon") || !strcmp(fChannel_, "all")) treeList.Add(fTreeTTmu);
    if (!strcmp(fChannel_, "electron") || !strcmp(fChannel_, "all")) treeList.Add(fTreeTTe);

    fTree = TTree::MergeTrees(&treeList);
  }

  return fTree;
}


void RandomSubsetCreatorLeptonJets::DrawEvents(TTree* tempTree, double nEventsPE) {
  std::cout << "nEventsPE: " << nEventsPE << std::endl;

  mcWeight = 1;

  fTree->CopyAddresses(tempTree);

  TRandom3* random = new TRandom3(0);

  int permsMC  = tempTree->GetEntries("");

  std::vector<std::string> vWeight;
  boost::split( vWeight, fWeight, boost::is_any_of("-*"));

  double maxMCWeight = 1.; // calculate upper bound for combined MCWeight
  //double maxMCWeightB = 1.;

  for (unsigned int i = 0; i < vWeight.size(); ++i) {
    std::cout << vWeight[i] << std::endl;
    if (strncmp((TString) vWeight[i], "pdfWeights", 10)) maxMCWeight *= tempTree->GetMaximum((TString) vWeight[i]);
    else maxMCWeight *= tempTree->GetMaximum("pdfWeights");
  }

  std::cout << "maxMCWeight(" << fWeight << "):" << maxMCWeight  << std::endl;

  if (maxMCWeight ==  0) { std::cout << "Weight not active?" << std::endl; }
  if (maxMCWeight == -1) { std::cout << "Running over data?" << std::endl; }

  std::cout << "while (eventsDrawn < eventsPE)..." << std::endl;

  int eventsDrawn = 0;
  int nAttempts = 0;

  boost::progress_display progress((int)nEventsPE, std::cout);

  while (eventsDrawn < (int)nEventsPE) {
    int drawn = random->Integer(permsMC);
    tempTree->GetEntry(drawn);
    ++nAttempts;
    //if (combi != 0) tempTree->GetEntry(drawn-combi); // Biases towards events with more combis
    if (combi == 0) {
      double eventWeight = 1.;
      for (unsigned int i = 0; i < vWeight.size(); ++i) {
        if (!strcmp((TString) vWeight[i], "MCWeight")) eventWeight *= MCWeight;
        else if (!strcmp((TString) vWeight[i], "muWeight")) eventWeight *= muWeight;
        else if (!strcmp((TString) vWeight[i], "PUWeight")) eventWeight *= PUWeight;
        else if (!strcmp((TString) vWeight[i], "PUWeightUp")) eventWeight *= PUWeightUp;
        else if (!strcmp((TString) vWeight[i], "PUWeightDown")) eventWeight *= PUWeightDown;
        else if (!strcmp((TString) vWeight[i], "bWeight")) eventWeight *= bWeight;
        else if (!strcmp((TString) vWeight[i], "bWeight_bTagSFUp")) eventWeight *= bWeight_bTagSFUp;
        else if (!strcmp((TString) vWeight[i], "bWeight_bTagSFDown")) eventWeight *= bWeight_bTagSFDown;
        else if (!strcmp((TString) vWeight[i], "bWeight_misTagSFUp")) eventWeight *= bWeight_misTagSFUp;
        else if (!strcmp((TString) vWeight[i], "bWeight_misTagSFDown")) eventWeight *= bWeight_misTagSFDown;
        else if (!strncmp((TString) vWeight[i], "pdfWeights", 10)) {
          std::string sub = vWeight[i].substr(11);
          int pdfWeightN = atol(sub.c_str());
          //std::cout << "pdfWeightN: " << pdfWeightN << std::endl;
          eventWeight *= pdfWeights[pdfWeightN];
        }
      }

      if (eventWeight > random->Uniform(0, maxMCWeight)) {
        for (int iComb = 0; iComb < 4; iComb++) {
          tempTree->GetEntry(drawn + iComb);

          if ((iComb != 0 && combi == 0)) {
            break;
          }

          fTree->Fill();
        }

        //std::cout << "mcWeight: " << mcWeight << std::endl;

        if (mcWeight < 0) {
          eventsDrawn += -1;
          //progress    += -1;
        }
        else {
          ++eventsDrawn;
          ++progress;
        }
      }
    }
  }
  std::cout << eventsDrawn << " events drawn in " << nAttempts << " attempts." << std::endl;
}


TTree* RandomSubsetCreatorLeptonJets::PrepareTree(TString file) {
  TChain* chain = new TChain("analyzeHitFit/eventTree");
  chain->Add(file);

  chain->SetBranchStatus("*",0);
  chain->SetBranchStatus("target", 1);
  chain->SetBranchStatus("hadTopMass", 1);
  chain->SetBranchStatus("hadWRawMass", 1);
  chain->SetBranchStatus("leptonPt", 1);
  chain->SetBranchStatus("leptonC", 1);
  chain->SetBranchStatus("leptonId", 1);
  chain->SetBranchStatus("hitFitProb", 1);
  chain->SetBranchStatus("run", 1);
  chain->SetBranchStatus("luminosityBlock", 1);
  chain->SetBranchStatus("event", 1);
  chain->SetBranchStatus("combi", 1);
  chain->SetBranchStatus("deltaThetaHadWHadB", 1);
  chain->SetBranchStatus("deltaThetaHadQHadQBar", 1);
  chain->SetBranchStatus("PUWeight", 1);
  chain->SetBranchStatus("PUWeightUp", 1);
  chain->SetBranchStatus("PUWeightDown", 1);
  chain->SetBranchStatus("muWeight", 1);
  chain->SetBranchStatus("bWeight", 1);
  chain->SetBranchStatus("bWeight_bTagSFUp", 1);
  chain->SetBranchStatus("bWeight_bTagSFDown", 1);
  chain->SetBranchStatus("bWeight_misTagSFUp", 1);
  chain->SetBranchStatus("bWeight_misTagSFDown", 1);
  chain->SetBranchStatus("MCWeight", 1);
  chain->SetBranchStatus("mcWeight", 1);
  chain->SetBranchStatus("nVertex", 1);
  chain->SetBranchStatus("hadQBCSV", 1);
  chain->SetBranchStatus("hadQBarBCSV", 1);
  chain->SetBranchStatus("hadBBCSV", 1);
  chain->SetBranchStatus("lepBBCSV", 1);
  chain->SetBranchStatus("jetsPt", 1);
  chain->SetBranchStatus("pdfWeights", 1);

  TString selection(selection_);
  selection += " & hadQBCSV<";    selection += fBDisc;
  selection += " & hadQBarBCSV<"; selection += fBDisc;
  selection += " & hadBBCSV>";    selection += fBDisc;
  selection += " & lepBBCSV>";    selection += fBDisc;

  std::cout << file << ": " << chain->GetEntries(selection + " & combi==0") << " events" << std::endl;

  return chain->CopyTree(selection);
}
