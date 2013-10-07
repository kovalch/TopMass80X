#include "TChain.h"
#include "TFile.h"
#include "TString.h"
#include "TTree.h"

#include <string>
#include <vector>

void SkimTree()
{
  std::vector<std::string> files = {//"MJP12*_v1_data.root"
     //  "Z2_S12_ABS_JES_096_161_5_sig.root","Z2_S12_ABS_JES_098_161_5_sig.root","Z2_S12_ABS_JES_100_161_5_sig.root","Z2_S12_ABS_JES_102_161_5_sig.root","Z2_S12_ABS_JES_104_161_5_sig.root"
     // ,"Z2_S12_ABS_JES_096_163_5_sig.root","Z2_S12_ABS_JES_098_163_5_sig.root","Z2_S12_ABS_JES_100_163_5_sig.root","Z2_S12_ABS_JES_102_163_5_sig.root","Z2_S12_ABS_JES_104_163_5_sig.root"
     // ,"Z2_S12_ABS_JES_096_166_5_sig.root","Z2_S12_ABS_JES_098_166_5_sig.root","Z2_S12_ABS_JES_100_166_5_sig.root","Z2_S12_ABS_JES_102_166_5_sig.root","Z2_S12_ABS_JES_104_166_5_sig.root"
     // ,"Z2_S12_ABS_JES_096_169_5_sig.root","Z2_S12_ABS_JES_098_169_5_sig.root","Z2_S12_ABS_JES_100_169_5_sig.root","Z2_S12_ABS_JES_102_169_5_sig.root","Z2_S12_ABS_JES_104_169_5_sig.root"
     // ,"Z2_S12_*ABS_JES_096_172_5_sig.root","Z2_S12_*ABS_JES_098_172_5_sig.root","Z2_S12_*ABS_JES_100_172_5_sig.root","Z2_S12_*ABS_JES_102_172_5_sig.root","Z2_S12_*ABS_JES_104_172_5_sig.root"
     // ,"Z2_S12_ABS_JES_096_175_5_sig.root","Z2_S12_ABS_JES_098_175_5_sig.root","Z2_S12_ABS_JES_100_175_5_sig.root","Z2_S12_ABS_JES_102_175_5_sig.root","Z2_S12_ABS_JES_104_175_5_sig.root"
     // ,"Z2_S12_ABS_JES_096_178_5_sig.root","Z2_S12_ABS_JES_098_178_5_sig.root","Z2_S12_ABS_JES_100_178_5_sig.root","Z2_S12_ABS_JES_102_178_5_sig.root","Z2_S12_ABS_JES_104_178_5_sig.root"
     // ,"Z2_S12_ABS_JES_096_181_5_sig.root","Z2_S12_ABS_JES_098_181_5_sig.root","Z2_S12_ABS_JES_100_181_5_sig.root","Z2_S12_ABS_JES_102_181_5_sig.root","Z2_S12_ABS_JES_104_181_5_sig.root"
     // ,"Z2_S12_ABS_JES_096_184_5_sig.root","Z2_S12_ABS_JES_098_184_5_sig.root","Z2_S12_ABS_JES_100_184_5_sig.root","Z2_S12_ABS_JES_102_184_5_sig.root","Z2_S12_ABS_JES_104_184_5_sig.root"
     // ,"Z2_S12_Matching_Up_sig.root","Z2_S12_Matching_Down_sig.root","Z2_S12_Scale_Up_sig.root","Z2_S12_Scale_Down_sig.root"
    "Z2_S12_MassiveBinDecay_ABS_JES_096_172_5_sig.root","Z2_S12_MassiveBinDecay_ABS_JES_098_172_5_sig.root"/*,"Z2_S12_MassiveBinDecay_ABS_JES_100_172_5_sig.root"*/,"Z2_S12_MassiveBinDecay_ABS_JES_102_172_5_sig.root","Z2_S12_MassiveBinDecay_ABS_JES_104_172_5_sig.root"
     // ,"Z2_S12_MCNLO*_sig.root","Z2_S12_POWHEG_sig.root","Z2_S12_POWHER_sig.root"
     // ,"Z2_S12_*JER_Up_sig.root","Z2_S12_*JER_Down_sig.root","Z2_S12_*BJES_Up_sig.root","Z2_S12_*BJES_Down_sig.root"
     // ,"Z2_S12_*CORFLAVQUARKJES_Up_sig.root","Z2_S12_*CORFLAVQUARKJES_Down_sig.root","Z2_S12_*CORTOTNOFLAVJES_Up_sig.root","Z2_S12_*CORTOTNOFLAVJES_Down_sig.root"
     // ,"Z2_S12_*P11mpiHi_sig.root","Z2_S12_*P11NoCR_sig.root","Z2_S12_*P11_sig.root","Z2_S12_*P11TeV_sig.root"
     // "QCDMixing_MJPS12*_v1_data.root"
  };

  for(std::string& fileName : files){
    std::cout << fileName << std::endl;
    std::string path = "/scratch/hh/dust/naf/cms/user/eschliec/TopMass/2012/";
    std::string oldPath = path+"/02/";
    std::string newPath = path+"/Skim_02/";
    std::string treeFolder = "analyzeKinFit";
    TChain* chain = new TChain((treeFolder+std::string("/eventTree")).c_str());
    //chain->Add((oldPath+std::string("Z2_S12_Had*_ABS_JES_100_172_5_sig.root")).c_str());
    //chain->Add((oldPath+std::string("Z2_S12_Semi_ABS_JES_100_172_5_sig.root")).c_str());
    //chain->Add((oldPath+std::string("Z2_S12_Lept_ABS_JES_100_172_5_sig.root")).c_str());
    chain->Add((oldPath+fileName).c_str());

    //TFile* newFile = TFile::Open((newPath+std::string("Z2_S12_ABS_JES_100_172_5_sig_5.root")).c_str(), "RECREATE");
    TString newFileName = fileName;
    newFileName.ReplaceAll("*","");
    TFile* newFile = TFile::Open((newPath+std::string(newFileName)).c_str(), "RECREATE");
    newFile->mkdir(treeFolder.c_str())->cd();
    std::string selection = "top.fitProb[0] > 0.1 && (pow(top.fitB1[0].Eta()-top.fitB2[0].Eta(),2) + pow(TVector2::Phi_mpi_pi(top.fitB1[0].Phi()-top.fitB2[0].Phi()),2)) > 1.5*1.5";
    //std::string selection = "top.fitProb > 0.01";
    chain->CopyTree(selection.c_str());

    newFile->Write();
    newFile->Close();
    delete chain;
  }
}
