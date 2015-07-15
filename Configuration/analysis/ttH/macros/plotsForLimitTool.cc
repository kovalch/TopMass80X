#include <TString.h>
#include <TFile.h>
#include <TTree.h>
#include <TBranch.h>
#include <TROOT.h>
#include <TFile.h>
#include <TH1D.h>
#include <TH1.h>
#include <TObjArray.h>

#include <iostream>
#include <fstream>
#include <cstdlib>
#include <string>
#include <sstream>
#include <typeinfo>


void plotsForLimitTool()
{
    TString sampleArray[19] = {"BTAGDISCR_BPURITY_DOWN", "BTAGDISCR_BSTAT2_DOWN", "BTAGDISCR_CERR2_DOWN", "BTAGDISCR_LSTAT1_DOWN", "JES_DOWN", "BTAGDISCR_BPURITY_UP", "BTAGDISCR_BSTAT2_UP", "BTAGDISCR_CERR2_UP", "BTAGDISCR_LSTAT1_UP", "JES_UP", "BTAGDISCR_BSTAT1_DOWN", "BTAGDISCR_CERR1_DOWN", "BTAGDISCR_LPURITY_DOWN", "BTAGDISCR_LSTAT2_DOWN", "Nominal", "BTAGDISCR_BSTAT1_UP",  "BTAGDISCR_CERR1_UP",  "BTAGDISCR_LPURITY_UP", "BTAGDISCR_LSTAT2_UP"};
    
    
    for (TString *sample = &sampleArray[0]; sample!=&sampleArray[19];++sample)
    {
        const TString* inpath = new TString("../Plots/"+*sample+"/combined/");
        std::cout<<"the sample is = "<<*sample<<std::endl;
        
        const TString* outpath = new TString("./");
        
        // Open the input files
        TFile* sampleFile = new TFile(inpath->Copy().Append("mvaEventA_mvaWeight_c1_step7_cate3_source.root"));
        if(!sampleFile) std::cout<<"\n\tWe didn't find the "+*sample+" input!!\n";
        
        if (*sample == "Nominal")
        {
            TH1D* histo_data = (TH1D*) sampleFile->Get("mvaEventA_mvaWeight_c1_step7_cate3_data");
            TH1D* histo_bkg = (TH1D*) sampleFile->Get("mvaEventA_mvaWeight_c1_step7_cate3_allmc");
            TH1D* histo_signalmc = (TH1D*) sampleFile->Get("mvaEventA_mvaWeight_c1_step7_cate3_signalmc");
            
            if (!histo_data) std::cout<<"\n\tWe didn't find the data histogram!!\n";
            if (!histo_bkg) std::cout<<"\n\tWe didn't find the bckg histogram!!\n";
            if (!histo_signalmc) std::cout<<"\n\tWe didn't find the signalmc histogram!!\n";
            
            // Access integral values to determine observation and rates
            std::cout<<"The integral for data is = "<<histo_data->Integral()<<std::endl;
            std::cout<<"The integral for bkgd is = "<<histo_bkg->Integral()<<std::endl;
            std::cout<<"The integral for signal is = "<<histo_signalmc->Integral()<<std::endl;
            
            TFile* outputFile(0);
            outputFile = new TFile("mvaEventA_mvaWeight_"+*sample+".root","RECREATE");
            
            // Write histograms to file
            histo_data->Write("data_obs");
            histo_bkg->Write("bkgd");
            histo_signalmc->Write("sig");
            outputFile->Write();
            outputFile->Close();
        }
        else
        {
            TH1D* histo_bkg = (TH1D*) sampleFile->Get("mvaEventA_mvaWeight_c1_step7_cate3_allmc");
            TH1D* histo_signalmc = (TH1D*) sampleFile->Get("mvaEventA_mvaWeight_c1_step7_cate3_signalmc");
            
            if (!histo_bkg) std::cout<<"\n\tWe didn't find the bckg histogram!!\n";
            if (!histo_signalmc) std::cout<<"\n\tWe didn't find the signalmc histogram!!\n";
            
            std::cout<<"The integral for"+*sample+"bkgd is = "<<histo_bkg->Integral()<<std::endl;
            std::cout<<"The integral for"+*sample+"signal is = "<<histo_signalmc->Integral()<<std::endl;
            
            TFile* outputFile(0);
            outputFile = new TFile("mvaEventA_mvaWeight_"+*sample+".root","RECREATE");
            
            if (*sample == "BTAGDISCR_BPURITY_DOWN") histo_bkg->Write("bkgd_BTAGDISCR_BPURITY_Down");
            if (*sample == "BTAGDISCR_BPURITY_UP") histo_bkg->Write("bkgd_BTAGDISCR_BPURITY_Up");
            if (*sample == "BTAGDISCR_BSTAT1_DOWN") histo_bkg->Write("bkgd_BTAGDISCR_BSTAT1_Down");
            if (*sample == "BTAGDISCR_BSTAT1_UP") histo_bkg->Write("bkgd_BTAGDISCR_BSTAT1_Up");
            if (*sample == "BTAGDISCR_BSTAT2_DOWN") histo_bkg->Write("bkgd_BTAGDISCR_BSTAT2_Down");
            if (*sample == "BTAGDISCR_BSTAT2_UP") histo_bkg->Write("bkgd_BTAGDISCR_BSTAT2_Up");
            if (*sample == "BTAGDISCR_CERR1_DOWN") histo_bkg->Write("bkgd_BTAGDISCR_CERR1_Down");
            if (*sample == "BTAGDISCR_CERR1_UP") histo_bkg->Write("bkgd_BTAGDISCR_CERR1_Up");
            if (*sample == "BTAGDISCR_CERR2_DOWN") histo_bkg->Write("bkgd_BTAGDISCR_CERR2_Down");
            if (*sample == "BTAGDISCR_CERR2_UP") histo_bkg->Write("bkgd_BTAGDISCR_CERR2_Up");
            if (*sample == "BTAGDISCR_LPURITY_DOWN") histo_bkg->Write("bkgd_BTAGDISCR_LPURITY_Down");
            if (*sample == "BTAGDISCR_LPURITY_UP") histo_bkg->Write("bkgd_BTAGDISCR_LPURITY_Up");
            if (*sample == "BTAGDISCR_LSTAT1_DOWN") histo_bkg->Write("bkgd_BTAGDISCR_LSTAT1_Down");
            if (*sample == "BTAGDISCR_LSTAT1_UP") histo_bkg->Write("bkgd_BTAGDISCR_LSTAT1_Up");
            if (*sample == "BTAGDISCR_LSTAT2_DOWN") histo_bkg->Write("bkgd_BTAGDISCR_LSTAT2_Down");
            if (*sample == "BTAGDISCR_LSTAT2_UP") histo_bkg->Write("bkgd_BTAGDISCR_LSTAT2_Up");
            if (*sample == "JES_DOWN") histo_bkg->Write("bkgd_JES_Down");
            if (*sample == "JES_UP") histo_bkg->Write("bkgd_JES_Up");
            
            if (*sample == "BTAGDISCR_BPURITY_DOWN") histo_signalmc->Write("sig_BTAGDISCR_BPURITY_Down");
            if (*sample == "BTAGDISCR_BPURITY_UP") histo_signalmc->Write("sig_BTAGDISCR_BPURITY_Up");
            if (*sample == "BTAGDISCR_BSTAT1_DOWN") histo_signalmc->Write("sig_BTAGDISCR_BSTAT1_Down");
            if (*sample == "BTAGDISCR_BSTAT1_UP") histo_signalmc->Write("sig_BTAGDISCR_BSTAT1_Up");
            if (*sample == "BTAGDISCR_BSTAT2_DOWN") histo_signalmc->Write("sig_BTAGDISCR_BSTAT2_Down");
            if (*sample == "BTAGDISCR_BSTAT2_UP") histo_signalmc->Write("sig_BTAGDISCR_BSTAT2_Up");
            if (*sample == "BTAGDISCR_CERR1_DOWN") histo_signalmc->Write("sig_BTAGDISCR_CERR1_Down");
            if (*sample == "BTAGDISCR_CERR1_UP") histo_signalmc->Write("sig_BTAGDISCR_CERR1_Up");
            if (*sample == "BTAGDISCR_CERR2_DOWN") histo_signalmc->Write("sig_BTAGDISCR_CERR2_Down");
            if (*sample == "BTAGDISCR_CERR2_UP") histo_signalmc->Write("sig_BTAGDISCR_CERR2_Up");
            if (*sample == "BTAGDISCR_LPURITY_DOWN") histo_signalmc->Write("sig_BTAGDISCR_LPURITY_Down");
            if (*sample == "BTAGDISCR_LPURITY_UP") histo_signalmc->Write("sig_BTAGDISCR_LPURITY_Up");
            if (*sample == "BTAGDISCR_LSTAT1_DOWN") histo_signalmc->Write("sig_BTAGDISCR_LSTAT1_Down");
            if (*sample == "BTAGDISCR_LSTAT1_UP") histo_signalmc->Write("sig_BTAGDISCR_LSTAT1_Up");
            if (*sample == "BTAGDISCR_LSTAT2_DOWN") histo_signalmc->Write("sig_BTAGDISCR_LSTAT2_Down");
            if (*sample == "BTAGDISCR_LSTAT2_UP") histo_signalmc->Write("sig_BTAGDISCR_LSTAT2_Up");
            if (*sample == "JES_DOWN") histo_signalmc->Write("sig_JES_Down");
            if (*sample == "JES_UP") histo_signalmc->Write("sig_JES_Up");
            
            outputFile->Write();
            outputFile->Close();
        }
    }
}
    
    
    
    
    
    
    
    
    
    
    
    
    
    