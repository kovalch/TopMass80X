{
  TFile *_file1 = TFile::Open("/nfs/dust/cms/user/eschliec/TopMass/2012/Skim_05B/Run2012_Mixing8_alljets.root");
  analyzeKinFit.cd();
  eventTree->Draw("top.fitTop1.M()>>hMix8(450,100,550)","weight.combinedWeight*(top.fitProb > 0.1 && (pow(top.fitB1.Eta()-top.fitB2.Eta(),2) + pow(TVector2::Phi_mpi_pi(top.fitB1.Phi()-top.fitB2.Phi()),2)) > 2.0*2.0 && max(max(max(max(max(top.recoJetIdxB1,top.recoJetIdxW1Prod1),top.recoJetIdxW1Prod2),top.recoJetIdxB2),top.recoJetIdxW2Prod1),top.recoJetIdxW2Prod2) < 6 && jet.alternativeJet[3].Pt() > 60 && jet.jet[3].Pt() > 60)", "norm");
  TFile *_file2 = TFile::Open("/nfs/dust/cms/user/eschliec/TopMass/2012/Skim_05B/Run2012_Background_alljets.root");
  analyzeKinFit.cd();
  eventTree->Draw("top.fitTop1.M()>>hBkg(450,100,550)","weight.combinedWeight*(top.fitProb > 0.1 && (pow(top.fitB1.Eta()-top.fitB2.Eta(),2) + pow(TVector2::Phi_mpi_pi(top.fitB1.Phi()-top.fitB2.Phi()),2)) > 2.0*2.0 && max(max(max(max(max(top.recoJetIdxB1,top.recoJetIdxW1Prod1),top.recoJetIdxW1Prod2),top.recoJetIdxB2),top.recoJetIdxW2Prod1),top.recoJetIdxW2Prod2) < 6 && jet.alternativeJet[3].Pt() > 60 && jet.jet[3].Pt() > 60)", "sames norm");
}
