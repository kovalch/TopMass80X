{

  gSystem->AddIncludePath("-I../../TopEventTree/interface");
  gSystem->AddIncludePath("-I../../..");
  
  
  gROOT->ProcessLine(".L ../../TopEventTree/src/BRegJetEvent.cc+");
  gROOT->ProcessLine(".L ../../TopEventTree/src/JetEvent.cc+");
  gROOT->ProcessLine(".L ../../TopEventTree/src/TopEvent.cc+");
  gROOT->ProcessLine(".L ../../TopEventTree/src/WeightEvent.cc+");
  

  gROOT->ProcessLine(".x readTopTree.C++");


}
