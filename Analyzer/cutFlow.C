// module load pod
// pod-server start
// pod-submit  -n 40 -r ge
void cutFlow() {
  TString connect = gSystem.GetFromPipe("pod-info -c");
  TProof *p = TProof::Open(connect);
  TChain eventTree("analyzeKinFit/eventTree");
  eventTree.Add("/nfs/dust/cms/user/eschliec/GRID-CONTROL_JOBS/TopMassTreeWriter_06_Data01/MJP12B_Winter14_v1_data/*.root");
  eventTree.Add("/nfs/dust/cms/user/eschliec/GRID-CONTROL_JOBS/TopMassTreeWriter_06_Data01/MJP12C1_Winter14_v1_data/*.root");
  eventTree.Add("/nfs/dust/cms/user/eschliec/GRID-CONTROL_JOBS/TopMassTreeWriter_06_Data01/MJP12C2_Winter14_v1_data/*.root");
  eventTree.Add("/nfs/dust/cms/user/eschliec/GRID-CONTROL_JOBS/TopMassTreeWriter_06_Data01/MJP12D1_Winter14_v1_data/*.root");
  eventTree.Add("/nfs/dust/cms/user/eschliec/GRID-CONTROL_JOBS/TopMassTreeWriter_06_Data01/MJP12D2_Winter14_v1_data/*.root");
  eventTree.Add("/nfs/dust/cms/user/eschliec/GRID-CONTROL_JOBS/TopMassTreeWriter_06_Data01/MJP12D3_Winter14_v1_data/*.root");
  eventTree.SetProof(1);


 //six jets
 eventTree.Draw("1 >>hist(1,0,2)","jet.jet[3].Pt() > 60 && jet.jet[5].Pt() > 30 && jet.jet[4].Pt() > 40");
 std::cout << "6j selected events:" << (int)hist->GetBinContent(1) << '\n';
 
 //two b-jets in first 6 jets
 eventTree->Draw("1 >>hist(1,0,2)","jet.jet[3].Pt() > 60  && jet.jet[4].Pt() > 40 && jet.jet[5].Pt() > 30 &&  Min$(max(max(max(max(max(top.recoJetIdxB1,top.recoJetIdxW1Prod1),top.recoJetIdxW1Prod2),top.recoJetIdxB2),top.recoJetIdxW2Prod1),top.recoJetIdxW2Prod2)) == 5","");
 std::cout << "2b selected events:" << hist->GetBinContent(1) << '\n';

 eventTree->Draw("1 >>hist(1,0,2)","jet.jet[3].Pt() > 60  && jet.jet[4].Pt() > 40 && jet.jet[5].Pt() > 30 &&  (((MinIf$(max(max(max(max(max(top.recoJetIdxB1,top.recoJetIdxW1Prod1),top.recoJetIdxW1Prod2),top.recoJetIdxB2),top.recoJetIdxW2Prod1),top.recoJetIdxW2Prod2),top.fitProb > 0.10 && (pow(top.fitB1.Eta()-top.fitB2.Eta(),2) + pow(TVector2::Phi_mpi_pi(top.fitB1.Phi()-top.fitB2.Phi()),2)) > 2.0*2.0) == 5) &&Length$(top.fitProb)>1) || ((Min$(max(max(max(max(max(top.recoJetIdxB1,top.recoJetIdxW1Prod1),top.recoJetIdxW1Prod2),top.recoJetIdxB2),top.recoJetIdxW2Prod1),top.recoJetIdxW2Prod2)) == 5) && top.fitProb[0] > 0.10 && (pow(top.fitB1[0].Eta()-top.fitB2[0].Eta(),2) + pow(TVector2::Phi_mpi_pi(top.fitB1[0].Phi()-top.fitB2[0].Phi()),2)) > 2.0*2.0 && Length$(top.fitProb)==1))","");
 std::cout << "final:" << hist->GetBinContent(1) << '\n';

 TChain eventTreeMC("analyzeKinFit/eventTree");
 eventTreeMC.Add("/nfs/dust/cms/user/eschliec/GRID-CONTROL_JOBS/TopMassTreeWriter_05_MC01/Z2_S12_ABS_JES_100_172_5_MadSpin_sig/*.root");
 eventTreeMC.SetProof(1);
 
 //six jets
 eventTreeMC.Draw("1 >>hist(1,0,2)","(jet.jet[3].Pt() > 60 && jet.jet[5].Pt() > 30 && jet.jet[4].Pt() > 40) * weight.combinedWeight* (0.5*TMath::Erf((jet.jet[3].Pt()-45.8627)/18.2471)+0.5)");
 std::cout << "6j selected events:" << (int)hist->GetBinContent(1) << '\n';

 
 //two b-jets in first 6 jets
 eventTreeMC->Draw("1 >>hist(1,0,2)","weight.combinedWeight*(0.5*TMath::Erf((jet.jet[3].Pt()-45.8627)/18.2471)+0.5)*(jet.jet[3].Pt() > 60  && jet.jet[4].Pt() > 40 && jet.jet[5].Pt() > 30 &&  Min$(max(max(max(max(max(top.recoJetIdxB1,top.recoJetIdxW1Prod1),top.recoJetIdxW1Prod2),top.recoJetIdxB2),top.recoJetIdxW2Prod1),top.recoJetIdxW2Prod2)) == 5)","");
 std::cout << "2b selected events:" << hist->GetBinContent(1) << '\n';
 
}
