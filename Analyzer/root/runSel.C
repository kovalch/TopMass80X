{

  gSystem->AddIncludePath("-I../../TopEventTree/interface");
  gSystem->AddIncludePath("-I../../..");
  
  
  //gROOT->ProcessLine(".L ../../TopEventTree/src/BRegJetEvent.cc+");
  gROOT->ProcessLine(".L ../../TopEventTree/src/JetEvent.cc+");
  gROOT->ProcessLine(".L ../../TopEventTree/src/TopEvent.cc+");
  gROOT->ProcessLine(".L ../../TopEventTree/src/WeightEvent.cc+");

  TFile *_file0 = TFile::Open("/nfs/dust/cms/user/mseidel/trees_2015/Summer12_TTJetsMS1725_1.00_muon/job_111_analyzeTop.root");
  analyzeHitFit->cd();
  eventTree->Process("topSel.C+");


}

