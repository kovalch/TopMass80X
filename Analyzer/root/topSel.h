//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Tue May  5 10:30:47 2015 by ROOT version 6.02/05
// from TTree eventTree/
// found on file: /nfs/dust/cms/user/mseidel/trees_2015/Summer12_TTJetsMS1725_1.00_muon/job_111_analyzeTop.root
//////////////////////////////////////////////////////////

#ifndef topSel_h
#define topSel_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <TSelector.h>

// Header file for the classes stored in the TTree if any.
#include "TopMass/TopEventTree/interface/TopEvent.h"
#include "TObject.h"
#include "TVector3.h"
#include "TLorentzVector.h"
#include "TopMass/TopEventTree/interface/JetEvent.h"
#include "TVector2.h"
#include "TopMass/TopEventTree/interface/WeightEvent.h"

#include "TH1.h"

class topSel : public TSelector {
public :
  TTree          *fChain;   //!pointer to the analyzed TTree or TChain

  TopEvent        *top;
  JetEvent        *jet;
  WeightEvent     *weight;
  // List of branches
  TBranch        *b_top;   //!
  TBranch        *b_jet;   //!
  TBranch        *b_weight;   //!

  topSel(TTree * /*tree*/ =0) : fChain(0) {
    top = new TopEvent();
    jet = new JetEvent();
    weight = new WeightEvent();
  }
  virtual ~topSel() { }
  virtual Int_t   Version() const { return 2; }
  virtual void    Begin(TTree *tree);
  virtual void    SlaveBegin(TTree *tree);
  virtual void    Init(TTree *tree);
  virtual Bool_t  Notify();
  virtual Bool_t  Process(Long64_t entry);
  virtual Int_t   GetEntry(Long64_t entry, Int_t getall = 0) { return fChain ? fChain->GetTree()->GetEntry(entry, getall) : 0; }
  virtual void    SetOption(const char *option) { fOption = option; }
  virtual void    SetObject(TObject *obj) { fObject = obj; }
  virtual void    SetInputList(TList *input) { fInput = input; }
  virtual TList  *GetOutputList() const { return fOutput; }
  virtual void    SlaveTerminate();
  virtual void    Terminate();

  ClassDef(topSel,0);


private:
  TH1D* hmass;
};

#endif

#ifdef topSel_cxx
void topSel::Init(TTree *tree)
{
  // The Init() function is called when the selector needs to initialize
  // a new tree or chain. Typically here the branch addresses and branch
  // pointers of the tree will be set.
  // It is normally not necessary to make changes to the generated
  // code, but the routine can be extended by the user if needed.
  // Init() will be called many times when running on PROOF
  // (once per file to be processed).

  // Set branch addresses and branch pointers
  if (!tree) return;
  fChain = tree;
  //fChain->SetMakeClass(1);

  fChain->SetBranchAddress("top.", &top, &b_top);
  fChain->SetBranchAddress("jet.", &jet, &b_jet);
  fChain->SetBranchAddress("weight.", &weight, &b_weight);
}

Bool_t topSel::Notify()
{
  // The Notify() function is called when a new file is opened. This
  // can be either for a new TTree in a TChain or when when a new TTree
  // is started when using PROOF. It is normally not necessary to make changes
  // to the generated code, but the routine can be extended by the
  // user if needed. The return value is currently not used.

  return kTRUE;
}

#endif // #ifdef topSel_cxx
