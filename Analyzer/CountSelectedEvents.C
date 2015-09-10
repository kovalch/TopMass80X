void CountSelectedEvents()
{
  analyzeKinFit->cd();
  eventTree->Draw("1 >>hist(1,0,2)","weight.combinedWeight*((top.combinationType[0]!=0)?(0.5*TMath::Erf((jet.jet[3].Pt()-45.8627)/18.2471)+0.5):1.0)*(top.fitProb[0] > 0.10 && (pow(top.fitB1[0].Eta()-top.fitB2[0].Eta(),2) + pow(TVector2::Phi_mpi_pi(top.fitB1[0].Phi()-top.fitB2[0].Phi()),2)) > 2.0*2.0 && max(max(max(max(max(top.recoJetIdxB1[0],top.recoJetIdxW1Prod1[0]),top.recoJetIdxW1Prod2[0]),top.recoJetIdxB2[0]),top.recoJetIdxW2Prod1[0]),top.recoJetIdxW2Prod2[0]) < 6 && jet.jet[3].Pt() > 60 && jet.jet[0].Pt() > 100)");
  std::cout << "selected events:" << hist->GetBinContent(1) << "  fsig:" << hist->GetBinContent(1)/(double)5560 << '\n';
}
