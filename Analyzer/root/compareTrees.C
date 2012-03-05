#include "TTree.h"
#include "TFile.h"
#include "TRandom.h"
   
int run, event;
Float_t x,y,z;

void compareTrees() {
  // The two TTrees created above are compared.
  // The subset of entries in the small TTree must be identical
  // to the entries in the original TTree.

  TFile *f = new TFile("/scratch/hh/current/cms/user/mseidel/Run2011/analyzeTop_event2bW.root");
  TTree *T  = (TTree*)f->Get("eventTree");
  T->BuildIndex("run","event");
  TFile *ff = new TFile("/scratch/hh/current/cms/user/mseidel/Run2011_MuHad/analyzeTop_event2bW.root");
  TTree *TF = (TTree*)ff->Get("eventTree");
  TF->BuildIndex("run","event");
  
  int run, event;
  double leptonPt, hadBRawPt, hadQRawPt, hadQBarRawPt;
  T->SetBranchAddress("run", &run);
  T->SetBranchAddress("event", &event);
  T->SetBranchAddress("leptonPt", &leptonPt);
  T->SetBranchAddress("hadBRawPt", &hadBRawPt);
  T->SetBranchAddress("hadQRawPt", &hadQRawPt);
  T->SetBranchAddress("hadQBarRawPt", &hadQBarRawPt);
  
  int frun, fevent;
  double fleptonPt, fhadBRawPt, fhadQRawPt, fhadQBarRawPt;
  TF->SetBranchAddress("run", &frun);
  TF->SetBranchAddress("event", &fevent);
  TF->SetBranchAddress("leptonPt", &fleptonPt);
  TF->SetBranchAddress("hadBRawPt", &fhadBRawPt);
  TF->SetBranchAddress("hadQRawPt", &fhadQRawPt);
  TF->SetBranchAddress("hadQBarRawPt", &fhadQBarRawPt);

  printf("Events in SingleMu only\n");
  int nentries = T->GetEntries();
  int nok = 0;
  for (int i = 0; i < nentries; ++i) {
    T->GetEntry(i);
    if (TF->GetEntryWithIndex(run, event) < 0) {
      printf("i=%lld, run=%d, event=%d, leptonPt=%d, hadBPt=%d, hadQPt=%d, hadQBarPt=%d \n", i, run, event, leptonPt, hadBRawPt, hadQRawPt, hadQBarRawPt);
      ++nok;
    } 
  }
  printf("nok = %d, fentries=%lld\n", nok, nentries);
  //*/
  
  printf("Events in MuHad only\n");
  nentries = TF->GetEntries();
  nok = 0;
  for (int i = 0; i < nentries; i++) {
    TF->GetEntry(i);
    if (T->GetEntryWithIndex(frun, fevent) < 0) {
      printf("i=%lld, frun=%d, fevent=%d, leptonPt=%d, hadBPt=%d, hadQPt=%d, hadQBarPt=%d, \n", i, frun, fevent, fleptonPt, fhadBRawPt, fhadQRawPt, fhadQBarRawPt);
      ++nok;
    } 
  }
  printf("nok = %d, fentries=%lld\n", nok, nentries);

  delete f;
  delete ff;
}
