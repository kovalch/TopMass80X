
#include <fstream>
#include <iostream>
// #include <cstdio>
#include <sstream>
// #include <cmath>
// #include <iomanip>

#include <TCanvas.h>
#include <TLegend.h>
// #include <TExec.h>
#include <TStyle.h>
#include <TSystem.h>
// #include <TMath.h>
// #include <TROOT.h>
#include <THStack.h>
#include <TDirectory.h>
#include <TFile.h>
#include <TString.h>
#include <TTree.h>
#include <TH1.h>
#include <TH2.h>
#include <TGraph.h>

#include "utils.h"
#include "Plotter.h"
#include "ttbarUtils.h"
#include "Samples.h"
#include "Sample.h"
#include "../../common/include/utils.h"
#include "../../common/include/plotterUtils.h"
#include "UsefulTools.h"
#include "Output.h"

#include "TUnfold.h"
#include <TLine.h>
#include <TMath.h>
#include "TUnfoldBinning.h"
#include "TUnfoldDensity.h"





Plotter::Plotter(const Samples& samples, const double lumi, double topxsec):

samples_(samples),
lumi_(lumi),
topxsec_(topxsec),
nD_(-999),
v_plotName_(std::vector<TString>(0)),
v_plotTitle_(std::vector<TString>(0)),
v_plotUnits_(std::vector<TString>(0)),

plotsFolder_("Plots/"),

//Control plots//
v_cpNBins_(0),
v_R1_(std::vector<Double_t>(0)),
v_R2_(std::vector<Double_t>(0)),
vv_SampleHist_(std::vector<std::vector<TH1D* > >(0)),
v_histRecoAllBins_(0),
vvv_SampleHist_(std::vector<std::vector<std::vector<TH1D* > > >(0)),
vv_UnderflowSampleHist_(std::vector<std::vector<TH1D* > >(0)),
vv_OverflowSampleHist_(std::vector<std::vector<TH1D* > >(0)),

//Unfolding//
detectorBinning_(0),
generatorBinning_(0),
detectorDistribution_(0),
generatorDistribution_(0),
histMigration_(0),
histBgrUo_(0),
histBgr_(0),
histData_(0),
unfoldedData_(0),

//Closure test//
histSR_(0),
histUf_(0),
histSG_(0),
histRWSG_(0),

histEff_(0),
histGen_(0),
histPurity_(0),
histStability_(0),
histRecoGen_(0),
histPurityAllBins_(0),
histStabilityAllBins_(0),
histRecoGenAllBins_(0),
histEffAllBins_(0),
histGenAllBins_(0),

v_coarseBins_(std::vector<std::vector<Double_t> >(0)),
v_fineBins_(std::vector<std::vector<Double_t> >(0)),
v_uReco_(std::vector<int >(0)),
v_oReco_(std::vector<int >(0)),
v_uTrue_(std::vector<int >(0)),
v_oTrue_(std::vector<int >(0)),

// Tree branches //
entry_(-999),
entry0_(-999),
eventWeight_(-999), 
trueLevelWeight_(-999),
trueLevelWeight0_(-999),
branchVals_(std::vector<Float_t>(0)),
branchValsGen_(std::vector<Float_t>(0)),
branchValsGen0_(std::vector<Float_t>(0)),

histMCReco_(0)



{
    //FIXME: make it in better way
    v_BR_.push_back(0.01166);//ee
    v_BR_.push_back(0.02332);//emu
    v_BR_.push_back(0.01166);//mumu
    v_BR_.push_back(0.04666);//combined
     
     // switch on histogram errors
    TH1::SetDefaultSumw2();
     
}



void Plotter::setOptions(const std::vector<TString> v_plotName)
{
    
    //std::cout<<"[Plotter]:--- Beginning of plot options settings\n\n";
    
    //Clearing previous declarations  
    this->clearMemory();
    
    v_plotName_ = v_plotName;
    nD_ = (int)v_plotName_.size();
    TString nDflag;
    if(nD_==1) nDflag="1d";
    else if(nD_==2) nDflag="2d";
    else{
         std::cerr<<"ERROR in Plotter! nD_ is not equal to 1 or 2.\n...break\n"<<std::endl;
         exit(237);
     }

    branchVals_ = std::vector<Float_t>(nD_,-999);
    branchValsGen_ = std::vector<Float_t>(nD_,-999);
    branchValsGen0_ = std::vector<Float_t>(nD_,-999);
    
    // reading Control Cards
    for(auto plotName : v_plotName_){
        const std::string ccFileName(common::CMSSW_BASE() + "/src/TopAnalysis/Configuration/analysis/diLeptonic/ControlCards/" + plotName + ".txt");
        std::ifstream ccFileStream(ccFileName.data(), std::ifstream::in);
        if (!ccFileStream.good()) {
          std::cerr<<"Error in Plotter! Cannot find file with name: "<< ccFileName <<"\n...break\n"<<std::endl;
          exit(12);
        }

       // Loop over all lines in ccFile
       bool isInBlock = false;
       while(ccFileStream.good()){
        // Read cc-File
            std::string line;
            getline(ccFileStream, line);
            line.erase(0, line.find_first_not_of(" \t"));
            if (line.size() == 0 || line[0] == '#') continue;

            // Loop over words in cc-File line and fill vWord
             std::vector<TString> vWord;
             std::string word;
             for (std::stringstream ss(line); ss >> word; ){
                vWord.push_back(word);
             }
             
             if(vWord.at(0) == "cpnbins")
             {
                     v_cpNBins_.push_back((vWord.at(1)).Atoi());
                     v_R1_.push_back((vWord.at(2)).Atof());
                     v_R2_.push_back((vWord.at(3)).Atof());
             }
             if(vWord.at(0) == "title")v_plotTitle_.push_back(vWord.at(1));
             
             if(vWord.at(0) == "units" ){
                if(vWord.size()>1) v_plotUnits_.push_back(vWord.at(1)); 
                if(vWord.size()<2) v_plotUnits_.push_back("");
             }
             
             
             //Pass to dimensional block
             if( (int)vWord.size()==1 && vWord.at(0) == nDflag ){
                 if(isInBlock) isInBlock = false;
                 if(!isInBlock) isInBlock = true;
            }
             
             // Reading only from dimensional block we are interesting
             if(isInBlock){
                 std::vector<Double_t > binsVector;
                 
                if(vWord.at(0) == "coarse"){
                     for(auto word : vWord) binsVector.push_back(word.Atof());
                     binsVector.erase(binsVector.begin(), binsVector.begin() + 1);
                     v_coarseBins_.push_back(binsVector);
                }
                if(vWord.at(0) == "fine"){
                     for(auto word : vWord) binsVector.push_back(word.Atof());
                     binsVector.erase(binsVector.begin(), binsVector.begin() + 1);
                     v_fineBins_.push_back(binsVector);
                }
                if(vWord.at(0) == "uoTrue")
                {
                    v_uTrue_.push_back(vWord.at(1).Atoi());
                    v_oTrue_.push_back(vWord.at(2).Atoi());
                }
                if(vWord.at(0) == "uoReco")
                {
                    v_uReco_.push_back(vWord.at(1).Atoi());
                    v_oReco_.push_back(vWord.at(2).Atoi());
                }
                 
             }//read from dimensional block
             else continue;
             
       }//stream line
       
    }//dimension name

    //Set binning scheme
    //Definition of TUnfoldBinning no detector and generator level
    detectorBinning_ = new TUnfoldBinning("detector");
    detectorDistribution_ = detectorBinning_->AddBinning("detectordistribution");
    for(int i=0;i<nD_;++i)detectorDistribution_->AddAxis(v_plotName_.at(i),(int)v_fineBins_.at(i).size()-1,v_fineBins_.at(i).data(),v_uReco_.at(i), v_oReco_.at(i));
    
    generatorBinning_ = new TUnfoldBinning("generator");
    generatorDistribution_ = generatorBinning_->AddBinning("signal");
    for(int i=0;i<nD_;++i)generatorDistribution_->AddAxis("gen_"+v_plotName_.at(i),(int)v_coarseBins_.at(i).size()-1, v_coarseBins_.at(i).data() , v_uTrue_.at(i), v_oTrue_.at(i));
    
    //std::cout<<"[Plotter]:--- Finishing of plot options settings\n\n";
}



void Plotter::clearMemory()
{
    v_plotName_.clear();
    v_plotTitle_.clear();
    v_plotUnits_.clear();
    nD_ = -999;
    branchVals_.clear();
    branchValsGen_.clear();
    branchValsGen0_.clear();
    v_coarseBins_.clear();
    v_fineBins_.clear();
    v_uReco_.clear();
    v_oReco_.clear();
    v_uTrue_.clear();
    v_oTrue_.clear();
    v_cpNBins_.clear();
    v_R1_.clear();
    v_R2_.clear();
    
    if(detectorDistribution_)detectorDistribution_->Delete();
    if(detectorBinning_)detectorBinning_->Delete();
    if(generatorDistribution_)generatorDistribution_->Delete();
    if(generatorBinning_)generatorBinning_->Delete();
    
}

void Plotter::clearMemoryPerSystematicChannel()
{
    for(auto v_p : vv_SampleHist_)
     for(auto p : v_p)
         delete p;
    vv_SampleHist_.clear();
    
    for(auto v_p : vv_UnderflowSampleHist_)
     for(auto p : v_p)
         delete p;
    vv_UnderflowSampleHist_.clear();
    
    for(auto v_p : vv_OverflowSampleHist_)
     for(auto p : v_p)
         delete p;
    vv_OverflowSampleHist_.clear();
    
    for(auto vv_p : vvv_SampleHist_)
     for(auto v_p : vv_p)
      for(auto p : v_p)
         delete p;
    vvv_SampleHist_.clear();
    
    if(histMigration_)histMigration_->Delete();
    if(histBgrUo_)histBgrUo_->Delete();
    if(histBgr_)histBgr_->Delete();
    if(histData_)histData_->Delete();
    if(unfoldedData_)unfoldedData_->Delete();
    
    //Closurer test//
    if(histSR_)histSR_->Delete();
    if(histUf_)histUf_->Delete();
    if(histSG_)histSG_->Delete();
    if(histRWSG_)histRWSG_->Delete();
    
    for(auto p : v_histRecoAllBins_)
        delete p;
    v_histRecoAllBins_.clear();
    
    if(histMCReco_)histMCReco_->Delete();
    
    
    if(histEff_)histEff_->Delete();
    if(histGen_)histGen_->Delete();
    if(histPurity_)histPurity_->Delete();
    if(histStability_)histStability_->Delete();
    if(histRecoGen_)histRecoGen_->Delete();
    
    
    if(histEffAllBins_)histEffAllBins_->Delete();
    if(histGenAllBins_)histGenAllBins_->Delete();
    if(histPurityAllBins_)histPurityAllBins_->Delete();
    if(histStabilityAllBins_)histStabilityAllBins_->Delete();
    if(histRecoGenAllBins_)histRecoGenAllBins_->Delete();
    
}



void Plotter::prepareHistograms(const std::vector<Sample>& v_sample)
{
    for(int ind=0;ind<nD_;++ind){
        std::vector<TH1D* > v_SampleHist;
        std::vector<TH1D* > v_UnderflowSampleHist;
        std::vector<TH1D* > v_OverflowSampleHist;
        
        for(size_t iSample = 0; iSample < v_sample.size(); ++iSample){
            const auto& sample(v_sample.at(iSample));
            TH1D* sampleHist = new TH1D(v_plotName_.at(ind)+"_cp" + std::to_string(ind) + std::to_string((int)iSample)  ,v_plotTitle_.at(ind),v_cpNBins_.at(ind),v_R1_.at(ind),v_R2_.at(ind));//FIXME: remove this: + std::to_string(ind) + std::to_string((int)iSample)
            sampleHist->SetFillColor(sample.color());
            sampleHist->SetLineColor(sample.color());
            v_SampleHist.push_back(sampleHist);
            v_UnderflowSampleHist.push_back((TH1D*)(sampleHist->Clone(sampleHist->GetName() + TString("u"))));
            v_OverflowSampleHist.push_back((TH1D*)(sampleHist->Clone(sampleHist->GetName() + TString("o"))));
        }
        vv_SampleHist_.push_back(v_SampleHist);
        vv_UnderflowSampleHist_.push_back(v_UnderflowSampleHist);
        vv_OverflowSampleHist_.push_back(v_OverflowSampleHist);
    }

    for(size_t iSample = 0; iSample < v_sample.size(); ++iSample){
        const auto& sample(v_sample.at(iSample));
        TH1* sampleHist =generatorBinning_->CreateHistogram("histRecoAllBins"+(TString)std::to_string((int)iSample).data());
        sampleHist->SetFillColor(sample.color());
        sampleHist->SetLineColor(sample.color());
        v_histRecoAllBins_.push_back((TH1*)sampleHist->Clone());
        sampleHist->Delete();
    }
    
    //setinng vvv_SampleHist_
    for(int ind=0;ind<nD_;++ind){
        std::vector<std::vector<TH1D* > > vv_SampleHist;
        for(int jnd=0;jnd<nD_;++jnd){
            if(ind==jnd)continue;
            const auto& v_bins = v_coarseBins_.at(jnd);
            for(int ibin=0;ibin<(int)v_bins.size()-1;++ibin){
                std::vector<TH1D* > v_SampleHist;
                for(size_t iSample = 0; iSample < v_sample.size(); ++iSample){
                    const auto& sample(v_sample.at(iSample));
                    TH1D* sampleHist = new TH1D(v_plotName_.at(ind)+"_cpBin" + std::to_string(ind) + std::to_string(ibin) + std::to_string((int)iSample)  ,v_plotTitle_.at(ind),
                                                (int)v_coarseBins_.at(ind).size()-1,v_coarseBins_.at(ind).data());
                    sampleHist->SetFillColor(sample.color());
                    sampleHist->SetLineColor(sample.color());
                    v_SampleHist.push_back(sampleHist);
                }
                vv_SampleHist.push_back(v_SampleHist);
            }
        }
        vvv_SampleHist_.push_back(vv_SampleHist);
    }
    
    
    
    //Definition of unfolding histograms
    histMigration_ = TUnfoldBinning::CreateHistogramOfMigrations(generatorBinning_,detectorBinning_,"histMigration");
    histBgrUo_ = detectorBinning_->CreateHistogram("histBgrUo");
    histBgr_ = detectorBinning_->CreateHistogram("histBgr");
    histData_ = detectorBinning_->CreateHistogram("histDataR");
    
    histMCReco_ = detectorBinning_->CreateHistogram("histDataReco");
    
    
        histEff_ = generatorBinning_->CreateHistogram("histEff",kTRUE);
        histGen_ = generatorBinning_->CreateHistogram("histGen",kTRUE);
        histPurity_ = generatorBinning_->CreateHistogram("histPurity",kTRUE);
        histStability_ = generatorBinning_->CreateHistogram("histStability",kTRUE);
        histRecoGen_ = generatorBinning_->CreateHistogram("histRecoGen",kTRUE);
        
        histEffAllBins_ = generatorBinning_->CreateHistogram("histEffAllBins");
        histGenAllBins_ = generatorBinning_->CreateHistogram("histGenAllBins");
        histPurityAllBins_ = generatorBinning_->CreateHistogram("histPurityAllBins");
        histStabilityAllBins_ = generatorBinning_->CreateHistogram("histStabilityAllBins");
        histRecoGenAllBins_ = generatorBinning_->CreateHistogram("histRecoGenAllBins");
    
    //Closure test//
    histSR_ = detectorBinning_->CreateHistogram("histSR");
    histSG_ = generatorBinning_->CreateHistogram("histSG");
    histRWSG_ = generatorBinning_->CreateHistogram("histRWSG");
        
}



void Plotter::producePlots()
{
    //std::cout<<"[Plotter]:--- Beginning of plot production\n\n";
    
    // Access correction factors
    const SystematicChannelFactors globalWeights = this->scaleFactors("step8");
    
    // Loop over all channels and systematics and produce plots
    const SystematicChannelSamples& m_systematicChannelSample(samples_.getSystematicChannelSamples());
    for(const auto& systematicChannelSamples : m_systematicChannelSample){                     //systematic loop
        const Systematic::Systematic& systematic(systematicChannelSamples.first);
        for(const auto& channelSample : systematicChannelSamples.second){                     //channel loop
            this->clearMemoryPerSystematicChannel();
            const Channel::Channel& channel(channelSample.first);
            const std::vector<Sample>& v_sample(channelSample.second);
            const std::vector<double>& v_weight(globalWeights.at(systematic).at(channel));
            
            //Set plotsFolder_
            plotsFolder_ = common::assignFolder("Plots",channel,systematic,v_plotName_.at(0)+"-"+v_plotName_.at(1));
            this->prepareHistograms(v_sample);
            
            TString sampleInfoFolder = common::assignFolder("Plots",channel,systematic) + "/sampleInfo.txt";
            Output sampleInfo("sample");
            
            
            double nData=0,nBg=0,nRecoSig=0,nGenSig=0;
            
            std::cout << systematic.name() << std::endl;
            for(size_t iSample = 0; iSample < v_sample.size(); ++iSample){                  //sample loop
                const auto& sample(v_sample.at(iSample));
                
                if(sample.legendEntry() == "QCD Multijet") continue;
                if(sample.legendEntry() == "W+Jets") continue;
                
                TFile* dataFile=new TFile(sample.inputFile().ReplaceAll("selectionRoot","ddaInput"));
                TTree *dataTree0=(TTree *) dataFile->Get("ttBar_treeVariables_step0"); 
                TTree *dataTree8=(TTree *) dataFile->Get("ttBar_treeVariables_step8");
                if(!dataTree0 || !dataTree8) {
                     std::cout<<"could not read 'ttBar_treeVariables_step' tree from " << dataFile->GetName() << "\n";
                 }
                 
                // set branches
                dataTree0->ResetBranchAddresses();
                dataTree0->SetBranchAddress("entry",&entry0_);
                dataTree0->SetBranchAddress("trueLevelWeight",&trueLevelWeight0_);
                for(int i=0;i<nD_;++i)dataTree0->SetBranchAddress("gen_"+v_plotName_.at(i),&(branchValsGen0_.at(i)));
                
                dataTree8->ResetBranchAddresses();
                dataTree8->SetBranchAddress("entry",&entry_);
                dataTree8->SetBranchAddress("eventWeight",&eventWeight_);
                dataTree8->SetBranchAddress("trueLevelWeight",&trueLevelWeight_);
                dataTree8->SetBranchAddress("isTopGen",&isTopGen_);
                dataTree8->SetBranchAddress("isKinReco",&isKinReco_);
                for(int i=0;i<nD_;++i)dataTree8->SetBranchAddress(v_plotName_.at(i),&(branchVals_.at(i)));
                for(int i=0;i<nD_;++i)dataTree8->SetBranchAddress("gen_"+v_plotName_.at(i),&(branchValsGen_.at(i)));
                
                sampleInfo.add("entries", utils::numToString(dataTree8->GetEntriesFast()));
                sampleInfo.add("weight",utils::numToString(v_weight.at(iSample)));
                sampleInfo.add("file",sample.inputFile());
                sampleInfo.add("name",sample.legendEntry());
                
                
                //loop over tree8 events
                Int_t lastEvent=0;
                for(Int_t ievent=0;ievent<dataTree8->GetEntriesFast();ievent++){
                    if(dataTree8->GetEntry(ievent)<=0) break;

                    //Control plots//
                    for(int ind=0;ind<nD_;++ind){
                        vv_SampleHist_.at(ind).at(iSample)->Fill(branchVals_.at(ind),eventWeight_*v_weight.at(iSample));
                        for(int jnd=0;jnd<nD_;++jnd){
                            if(ind==jnd)continue;
                            const auto& v_bins = v_coarseBins_.at(jnd);
                            if(branchVals_.at(jnd)<v_bins.at(0))vv_UnderflowSampleHist_.at(ind).at(iSample)->Fill(branchVals_.at(ind),eventWeight_*v_weight.at(iSample));
                            if(branchVals_.at(jnd)>v_bins.at((int)v_bins.size()-1))vv_OverflowSampleHist_.at(ind).at(iSample)->Fill(branchVals_.at(ind),eventWeight_*v_weight.at(iSample));
                            for(int ibin=0;ibin<(int)v_bins.size()-1;++ibin){
                                   if(branchVals_.at(jnd) >= v_bins.at(ibin) && branchVals_.at(jnd) < v_bins.at(ibin+1)){
                                        vvv_SampleHist_.at(ind).at(ibin).at(iSample)->Fill(branchVals_.at(ind),eventWeight_*v_weight.at(iSample));
                                    }
                            }
                        }
                    }

                    v_histRecoAllBins_.at(iSample)->Fill(genBin_(branchVals_),eventWeight_*v_weight.at(iSample));
                    
                    // ... //
                    
                    if(isTopGen_){
                        for(Int_t ievent0=lastEvent;ievent0<dataTree0->GetEntriesFast();ievent0++) {
                            if(dataTree0->GetEntry(ievent0)<=0) break;
                                
                                if(nD_==1)histGen_->Fill(branchValsGen0_.at(0),trueLevelWeight0_*v_weight.at(iSample));
                                if(nD_==2)((TH2*)histGen_)->Fill(branchValsGen0_.at(0),branchValsGen0_.at(1),trueLevelWeight0_*v_weight.at(iSample));
                                
                                histSG_->Fill(genBin_(branchValsGen0_),trueLevelWeight0_*v_weight.at(iSample));
                                //histSG_->Fill(genBin_(branchValsGen0_),1*v_weight.at(iSample));
                                //Closure test
                                //histRWSG_->Fill(genBin_(branchValsGen0_),(rewTopPtUp(branchValsGen0_.at(1)))*v_weight.at(iSample));//pt
                                //histRWSG_->Fill(genBin_(branchValsGen0_),(rewTopPtDown(branchValsGen0_.at(1)))*v_weight.at(iSample));//pt down
                                //histRWSG_->Fill(genBin_(branchValsGen0_),(rewTopRapidityUp(branchValsGen0_.at(0)))*v_weight.at(iSample));//y up
                                //histRWSG_->Fill(genBin_(branchValsGen0_),(rewTopRapidityDown(branchValsGen0_.at(0)))*v_weight.at(iSample));//y down
                                //histRWSG_->Fill(genBin_(branchValsGen0_),1*v_weight.at(iSample));//1
                                
                                histGenAllBins_->Fill(genBin_(branchValsGen0_),trueLevelWeight0_*v_weight.at(iSample));
                                
                                if(entry0_ != entry_)
                                {
                                    histMigration_->Fill(genBin_(branchValsGen0_),0.,trueLevelWeight0_*v_weight.at(iSample));
                                    //histMigration_->Fill(genBin_(branchValsGen0_),0.,1*v_weight.at(iSample));
                                    nGenSig=nGenSig+trueLevelWeight0_*v_weight.at(iSample);
                                }
                                else if(entry0_ == entry_)
                                {
                                    lastEvent = ievent0+1;
                                    break;
                                }
                        }
                    }
                    
                    if(v_sample.at(iSample).sampleType() == Sample::data){
                        histData_->Fill(recoBin_(branchVals_),eventWeight_*v_weight.at(iSample)); 
                        nData=nData+eventWeight_*v_weight.at(iSample);
                    }
                    if(v_sample.at(iSample).sampleType() != Sample::data && isKinReco_ && !isTopGen_){
                        histBgr_->Fill(recoBin_(branchVals_),eventWeight_*v_weight.at(iSample));
                        nBg=nBg+eventWeight_*v_weight.at(iSample);
                    }
                    if(isKinReco_&&isTopGen_){
                       histMigration_->Fill(genBin_(branchValsGen_),recoBin_(branchVals_),eventWeight_*v_weight.at(iSample));
                       histMigration_->Fill(genBin_(branchValsGen_),0.,trueLevelWeight_*v_weight.at(iSample)-eventWeight_*v_weight.at(iSample));
                       nRecoSig=nRecoSig+eventWeight_*v_weight.at(iSample);
                       nGenSig=nGenSig+trueLevelWeight_*v_weight.at(iSample);
                       histSR_->Fill(recoBin_(branchVals_),eventWeight_*v_weight.at(iSample));
                       
                       //double weightCT = 1.;
                       //Closure test
                       //weightCT = (rewTopPtUp(branchValsGen_.at(1)));//pt up
                       //weightCT = (rewTopPtDown(branchValsGen_.at(1)));//pt down
                       //weightCT = (rewTopRapidityUp(branchValsGen_.at(0)));//y up
                       //weightCT = (rewTopRapidityDown(branchValsGen_.at(0)));//y down
                       
                       //histMigration_->Fill(genBin_(branchValsGen_),recoBin_(branchVals_),1*v_weight.at(iSample));
                       //histSR_->Fill(recoBin_(branchVals_),weightCT*v_weight.at(iSample));
                       
                        if(nD_ == 1)
                        {
                            if((branchValsGen_.at(0) <v_coarseBins_.at(0).at(0) && v_uTrue_.at(0)==0) 
                                || (branchValsGen_.at(0)>v_coarseBins_.at(0).back() && v_oTrue_.at(0)==0) )
                            {
                                histBgrUo_->Fill(recoBin_(branchVals_),eventWeight_*v_weight.at(iSample));
                                //histBgrUo_->Fill(recoBin_(branchVals_),weightCT*v_weight.at(iSample));
                            }
                        }
                        else if(nD_ == 2)
                        {
                            if( (branchValsGen_.at(0) <v_coarseBins_.at(0).at(0) && v_uTrue_.at(0)==0) 
                                || (branchValsGen_.at(0)>v_coarseBins_.at(0).back() && v_oTrue_.at(0)==0)
                                || (branchValsGen_.at(1) <v_coarseBins_.at(1).at(0) && v_uTrue_.at(1)==0) 
                                || (branchValsGen_.at(1)>v_coarseBins_.at(1).back() && v_oTrue_.at(1)==0) )
                            {
                                histBgrUo_->Fill(recoBin_(branchVals_),eventWeight_*v_weight.at(iSample));
                                //histBgrUo_->Fill(recoBin_(branchVals_),weightCT*v_weight.at(iSample));
                            }
                        }
                        
                       Int_t recBin_temp = -999;
                       Int_t genBin_temp = -999;
                       
                       if(nD_ == 1){
                            recBin_temp = generatorDistribution_->GetGlobalBinNumber(branchVals_.at(0));
                            genBin_temp = generatorDistribution_->GetGlobalBinNumber(branchValsGen_.at(0));
                            histEff_->Fill(branchVals_.at(0),eventWeight_*v_weight.at(iSample));
                            histPurity_->Fill(branchVals_.at(0),eventWeight_*v_weight.at(iSample));
                            histStability_->Fill(branchValsGen_.at(0),eventWeight_*v_weight.at(iSample));
                       }
                       else if(nD_ == 2){
                            recBin_temp = generatorDistribution_->GetGlobalBinNumber(branchVals_.at(0),branchVals_.at(1));
                            genBin_temp = generatorDistribution_->GetGlobalBinNumber(branchValsGen_.at(0),branchValsGen_.at(1));
                            ((TH2*)histEff_)->Fill(branchVals_.at(0),branchVals_.at(1),eventWeight_*v_weight.at(iSample));
                            ((TH2*)histPurity_)->Fill(branchVals_.at(0),branchVals_.at(1),eventWeight_*v_weight.at(iSample));
                            ((TH2*)histStability_)->Fill(branchValsGen_.at(0),branchValsGen_.at(1),eventWeight_*v_weight.at(iSample));
                       }
                       histEffAllBins_->Fill(recBin_temp,eventWeight_*v_weight.at(iSample));
                       histPurityAllBins_->Fill(recBin_temp,eventWeight_*v_weight.at(iSample));
                       histStabilityAllBins_->Fill(genBin_temp,eventWeight_*v_weight.at(iSample));
                       
                       if(recBin_temp==genBin_temp){
                           if(nD_ == 1)histRecoGen_->Fill(branchVals_.at(0),eventWeight_*v_weight.at(iSample));
                           if(nD_ == 2)((TH2*)histRecoGen_)->Fill(branchVals_.at(0),branchVals_.at(1),eventWeight_*v_weight.at(iSample));
                           histRecoGenAllBins_->Fill(recBin_temp,eventWeight_*v_weight.at(iSample));
                       }
                       
                    }
                    if(v_sample.at(iSample).sampleType() != Sample::data && isKinReco_ ){
                        histMCReco_->Fill(recoBin_(branchVals_),eventWeight_*v_weight.at(iSample)); 
                    }
                    
                    
                }
                
            
                dataTree0->Delete();
                dataTree8->Delete();
                dataFile->Close();
                dataFile->Delete();
            
                
                
            }//samples loop
            
            //Cross sections//
            // Full 
                topxsec_ = (nData - nBg)/(nRecoSig/nGenSig)/lumi_/v_BR_.at(channel);
                sampleInfo.addLine("fullTopXSec = "+utils::numToString(topxsec_));
            sampleInfo.save(sampleInfoFolder);
            
            //Control plots//
            for(int ind=0;ind<nD_;++ind){
               writePlotCP( v_sample ,vv_SampleHist_.at(ind),ind);
                    if(v_uTrue_.at(!ind)&&nD_==2)writePlotCP( v_sample ,vv_UnderflowSampleHist_.at(ind),ind,-1);
                    if(v_oTrue_.at(!ind)&&nD_==2)writePlotCP( v_sample ,vv_OverflowSampleHist_.at(ind),ind,-2);
                for(int jnd=0;jnd<nD_;++jnd){
                  if(ind==jnd)continue;
                   const auto& v_bins = v_coarseBins_.at(jnd);
                   for(int ibin=0;ibin<(int)v_bins.size()-1;++ibin){
                      writePlotCP( v_sample ,vvv_SampleHist_.at(ind).at(ibin),ind,ibin);
                 }
                }
            }
            if(nD_==2)writePlotCPAllBins(v_sample,v_histRecoAllBins_);
            // ... //
            
            
            //Unfolding//
            if(nD_==2)runUnfolding(histMigration_,histData_,histBgr_,histBgrUo_,histUf_);
            //runUnfolding(histMigration_,histSR_,0,histBgrUo_,histUf_);
            //writePlotCT(histSG_,histRWSG_,histUf_);
            // ... //
            
            //Cross sections//
            // Diff
                if(nD_==2)writePlotXSec(unfoldedData_,histGen_);
            //writePlotEPS();
                if(nD_==2)writePlotEPSAllBins();
            // ...
            
        }//channel loop
        
    }//systematics loop
    
    //std::cout<<"\n[Plotter]:=== Finishing of plot production\n\n";
}


int Plotter::genBin_(const std::vector<float>& val)
{
    int bin = -999;
    if(nD_ == 1) bin = generatorDistribution_->GetGlobalBinNumber(val.at(0));
    if(nD_ == 2) bin = generatorDistribution_->GetGlobalBinNumber(val.at(0),val.at(1));
    return bin;
}



int Plotter::recoBin_(const std::vector<float>& val)
{
    int bin = -999;
    if(nD_ == 1) bin = detectorDistribution_->GetGlobalBinNumber(val.at(0));
    if(nD_ == 2) bin = detectorDistribution_->GetGlobalBinNumber(val.at(0),val.at(1));
    return bin;
}


//Closure test//
void Plotter::writePlotCT(TH1* histSG,TH1* histRW,TH1* histUf)
{
    // Prepare canvas and legend
    TCanvas* canvas = utils::setCanvas();
    TLegend* legend = utils::setLegend();
    canvas->cd();
    
    histSG->SetLineColor(2);
    histSG->SetLineWidth(2);
    histSG->SetStats(0);
    legend->AddEntry(histSG,"gen","l");
    histSG->Draw("hist");
    
    histRW->SetLineColor(kMagenta);
    histRW->SetLineWidth(2);
    histRW->SetLineStyle(2);
    legend->AddEntry(histRW,"reweighted gen","l");
    histRW->Draw("hist same");
    
    histUf->SetMarkerStyle(20);
    histUf->Sumw2();
    legend->AddEntry(histUf,"unfolded","pe");
    histUf->Draw("e1 same");
    
    //Lines and bins//
    double yMax = ( histUf->GetMaximum() > histRW->GetMaximum() ? histUf->GetMaximum()*1.2 : histRW->GetMaximum()*1.2 );
    int Nbins1 = (int)(v_coarseBins_.at(1).size()-1)+(int)(v_oTrue_.at(1)+v_uTrue_.at(1));
    int Nbins0 = (int)(v_coarseBins_.at(0).size()-1)+(int)(v_oTrue_.at(0)+v_uTrue_.at(0));
    TH1D hist("hist","",Nbins1,0.5,Nbins0*Nbins1+0.5);
    hist.SetTitle(";"+v_plotTitle_.at(1)+" "+v_plotUnits_.at(1));
    hist.SetStats(0);
    hist.GetXaxis()->SetTickLength(0);
    hist.SetAxisRange(0,yMax ,"Y");
        for(int i=1;i<=Nbins1;++i)
        {
            TString binName = "";
            int isUnderflow = v_uTrue_.at(1);
            int isOverflow =  v_oTrue_.at(1);
            if(isUnderflow&&i==1) binName = "underflow";
            else if(isOverflow&&i==Nbins1){
                 binName = "overflow";
            }
            else{
                binName = utils::makeBinTitle("",v_coarseBins_.at(1).at(i-isUnderflow-1),v_coarseBins_.at(1).at(i-isUnderflow));
            }
            hist.GetXaxis()->SetBinLabel(i,binName.Data());
        }
        hist.SetTitle( utils::makeTitleBins(v_plotTitle_.at(0) +" "+ v_plotUnits_.at(0),v_coarseBins_.at(0),v_uTrue_.at(0),v_oTrue_.at(0) ) );
    hist.Draw();
 
    histSG->Draw("hist same");
    histRW->Draw("hist same");
    histUf->Draw("e1 same");
    for(int i=1;i<Nbins1;++i)
        {
            double xLine = i*(Nbins0) + 0.5;
            TLine ln(xLine,-0.03*yMax,xLine,yMax);//FIXME: y1 to be related to pad size... 
            ln.DrawClone("same");
        }
    legend->Draw("same");
    utils::drawRatio(histUf, histRW, NULL,0.5,2,(TH1D*)hist.Clone());

    writeCanvas(canvas,"CT");
}


double Plotter::rewTopPtUp(double pt)
{
    //func y = ax+b
    double e = 0.2;
    double b = 1 + e;
    double a = (1. - b)/100.;
    return a*pt + b;
}


double Plotter::rewTopPtDown(double pt)
{
    //func y = ax+b
    double e = 0.2;
    double b = 1 - e;
    double a = (1. - b)/100.;
    return a*pt + b;
}


double Plotter::rewTopRapidityUp(double y)
{
    //func y = axx+b
    double b = 1 + 0.2;
    double a = (1 - 0.1 - b)/1.2/1.2;
    return a*y*y+b;
}


double Plotter::rewTopRapidityDown(double y)
{
    //func y = axx+b
    double b =1 - 0.2;
    double a = (1 + 0.1 - b)/1.2/1.2;
    return a*y*y+b;
}


void Plotter::runUnfolding(const TH2* histMigration,const TH1* histInput,
                           const TH1* histBgr,const TH1* histBgrUo,
                           TH1*& histUf)
{
     // preserve the area
    TUnfold::EConstraint constraintMode= TUnfold::kEConstraintArea;
    //TUnfold::EConstraint constraintMode= TUnfold::kEConstraintNone;
    
    // basic choice of regularisation scheme:
    //    curvature (second derivative)
   TUnfold::ERegMode regMode = TUnfold::kRegModeCurvature;
//     TUnfold::ERegMode regMode = TUnfold::kRegModeDerivative;
   // TUnfold::ERegMode regMode = TUnfold::kRegModeSize;
   // TUnfold::ERegMode regMode = TUnfold::kRegModeNone;

    // density flags
    TUnfoldDensity::EDensityMode densityFlags = TUnfoldDensity::kDensityModeBinWidth;

    // detailed steering for regularisation
    const char *REGULARISATION_DISTRIBUTION=0;
    const char *REGULARISATION_AXISSTEERING=0;//""
    //detectorBinning_->PrintStream(std::cout);
    //generatorBinning_->PrintStream(std::cout);
    
    // set up matrix of migrations
    
    TUnfoldDensity unfold(histMigration,TUnfold::kHistMapOutputHoriz,
                        regMode,constraintMode,densityFlags,
                        generatorBinning_,detectorBinning_,
                        REGULARISATION_DISTRIBUTION,
                        REGULARISATION_AXISSTEERING);
    
    // define the input vector (the measured data distribution)

    // special test to use input covariance matrix
    TH2D *inputEmatrix = detectorBinning_->CreateErrorMatrixHistogram("input_covar",true);
    for(int i=1;i<=inputEmatrix->GetNbinsX();i++) {
       Double_t e=histInput->GetBinError(i);
       inputEmatrix->SetBinContent(i,i,e*e);
       // test: non-zero covariance where variance is zero
       //   if(e<=0.) inputEmatrix->SetBinContent(i,i+1,1.0);
    }
    unfold.SetInput(histInput,1.0,0.0,inputEmatrix);
    //unfold.SetInput(histInput,0.0,0.0,inputEmatrix);
    
    if(histBgr){unfold.SubtractBackground(histBgr,"bgrSum");}
    if(histBgrUo)unfold.SubtractBackground(histBgrUo,"bgrUo");
      // run the unfolding
  Int_t nScan=50;
  TSpline *logTauX,*logTauY;
  //TSpline *rhoLogTau=0;
  TGraph *lCurve;
  // this method scans the parameter tau and finds the kink in the L curve
  // finally, the unfolding is done for the best choice of tau
  Int_t iBest=unfold.ScanLcurve(nScan,0.,0.,&lCurve,&logTauX,&logTauY);

//     // save graphs with one point to visualize best choice of tau
//     Double_t t[1],x[1],y[1];
//     logTauX->GetKnot(iBest,t[0],x[0]);
//     logTauY->GetKnot(iBest,t[0],y[0]);
//     TGraph *bestLcurve=new TGraph(1,x,y);
//     TGraph *bestLogTauLogChi2=new TGraph(1,t,x);

      // create graphs with one point to visualize best choice of tau
      //Double_t t[1],rho[1];
       Double_t x[1],y[1];
      //rhoLogTau->GetKnot(iBest,t[0],rho[0]);
      lCurve->GetPoint(iBest,x[0],y[0]);
      //TGraph *bestRhoLogTau=new TGraph(1,t,rho);
      TGraph *bestLCurve=new TGraph(1,x,y);
      //Double_t *tAll=new Double_t[nScan],*rhoAll=new Double_t[nScan];
      //for(Int_t i=0;i<nScan;i++) {
        //rhoLogTau->GetKnot(i,tAll[i],rhoAll[i]);
      //}
      //TGraph *knots=new TGraph(nScan,tAll,rhoAll);
    
  //unfold.DoUnfold(0); // do unfold with 0 tau^2 strength 
    
  std::cout<<"chi**2="<<unfold.GetChi2A()<<"+"<<unfold.GetChi2L()
      <<" / "<<unfold.GetNdf()<<"\n";
    
    
    TCanvas* canvas = utils::setCanvas();
    TLegend* legend = utils::setLegend();
    
    
    //Tau parameter
     // rhoLogTau->Draw();
///    knots->Draw("*ap");
///    bestRhoLogTau->SetMarkerColor(kRed);
///    bestRhoLogTau->Draw("*ap,same");
///    writeCanvas(canvas,"Tau");
    
    //L-curve
    lCurve->Draw("AL");
    bestLCurve->SetMarkerColor(kRed);
    bestLCurve->Draw("*");
    writeCanvas(canvas,"Lcurve");
    
    
    const char* distributionName = 0;
    const char* projectionMode = 0;
    bool useAxisBinning = kTRUE;
    
    
    //Data and GetFoldedOutput
///    TH1* foldedOutput = unfold.GetFoldedOutput("foldedOutput", "folded output", distributionName, projectionMode , kFALSE, kTRUE);
///    foldedOutput->SetStats(0);
///    foldedOutput->SetLineColor(1);
///    foldedOutput->Draw();
///    histInput->DrawClone("pe,same");
///    histMCReco_->SetLineColor(2);
///    histMCReco_->Draw("samehist");
///    legend->AddEntry(foldedOutput,"folded","l");
///    legend->AddEntry(histInput,"reco data","pe");
///    legend->AddEntry(histMCReco_,"mc reco","l");
///    legend->Draw("same");
///    writeCanvas(canvas,foldedOutput->GetName());
///    legend->Clear();
    
    
    
    histUf = unfold.GetOutput("histUf","unfolding result",distributionName,projectionMode,kFALSE);
    
    unfoldedData_ = unfold.GetOutput("unfoldedData_","unfolding result",distributionName,projectionMode,useAxisBinning);
    //unfoldedData_ = (TH1D*)unfold.GetOutput("unfoldedData","unfolding result",distributionName,"*[UO]",useAxisBinning);
    //unfoldedData_ = (TH1D*)unfold.GetOutput("unfoldedData","unfolding result",distributionName,"*",useAxisBinning);
    //unfoldedData_->GetYaxis()->SetRangeUser(0,11500);
    unfoldedData_->SetStats(0);
    unfoldedData_->Draw("colz");
    writeCanvas(canvas,unfoldedData_->GetName());
    
    
    // Probability matrix //TH2*
    TH2* probabilityMatrix = unfold.GetProbabilityMatrix("probabilityMatrix","probability matrix",useAxisBinning);
    probabilityMatrix->SetStats(0);
    probabilityMatrix->Draw("colz");
    writeCanvas(canvas,probabilityMatrix->GetName());
    
    // rho ij total //TH2*
    TH2* rhoIJtotal = unfold.GetRhoIJtotal("rhoIJtotal","rho ij total matrix",distributionName,projectionMode,useAxisBinning);
    rhoIJtotal->SetStats(0);
    rhoIJtotal->Draw("colz");
    writeCanvas(canvas,rhoIJtotal->GetName());
    
    
    // Convenient, rho ij total //TH2*
    TH2* rhoIJtotalConv = unfold.GetRhoIJtotal("rhoIJtotalConv","convenient rho ij total matrix",distributionName,projectionMode,kTRUE);
    rhoIJtotalConv->SetStats(0);
    TH2* errorIJtotal = unfold.GetEmatrixTotal("errorIJtotal","error ij total matrix",distributionName,projectionMode,kTRUE);
    TH1* unfoldedData = unfold.GetOutput("unfoldedData","unfolding result",distributionName,projectionMode,kFALSE);
    for(Int_t i=1;i<=rhoIJtotalConv->GetNbinsX();i++) {
         Double_t e_i=errorIJtotal->GetBinContent(i,i);
         if(e_i>0.0) e_i=TMath::Sqrt(e_i);
         else e_i=0.0;
         Double_t data_i=unfoldedData->GetBinContent(i);
         if(data_i>0.0&&e_i<data_i)
           rhoIJtotalConv->SetBinContent(i,i,e_i/data_i);
         else rhoIJtotalConv->SetBinContent(i,i,1);
      }
    rhoIJtotalConv->Draw("colz");
    writeCanvas(canvas,rhoIJtotalConv->GetName());
    
    
    
    // Error matrices info
    TH1* t_Ematrix = generatorBinning_->CreateHistogram("t_Ematrix");
    TH1* i_Ematrix = generatorBinning_->CreateHistogram("i_Ematrix");
    TH1* d_Ematrix = generatorBinning_->CreateHistogram("d_Ematrix");
    TH1* a_Ematrix = generatorBinning_->CreateHistogram("a_Ematrix");
    TH1* b_Ematrix = generatorBinning_->CreateHistogram("b_Ematrix");
    
    TH2* t_Ematrix_2d = unfold.GetEmatrixTotal("t_Ematrix_2d", "error matrix (t)");
    TH2* i_Ematrix_2d = unfold.GetEmatrixInput("i_Ematrix_2d", "error matrix (i)");
    TH2* d_Ematrix_2d = unfold.GetEmatrixSysBackgroundUncorr("d_Ematrix_2d", "error matrix (d)");
    TH2* a_Ematrix_2d = unfold.GetEmatrixSysUncorr("a_Ematrix_2d", "error matrix (a)");
    //std::cout << " TH2* b_Ematrix_2d=generatorBinning_->CreateErrorMatrixHistogram(\"b_Ematrix_2d\",kTRUE);   " << std::endl;
    //TH2* b_Ematrix_2d=generatorBinning_->CreateErrorMatrixHistogram("b_Ematrix_2d",kTRUE);
    //std::cout << " ... " << std::endl;
    //for(Int_t i=1;i<=i_Ematrix->GetNbinsX();i++)b_Ematrix_2d->SetBinContent(i,i, 1000 );
    //unfold.GetEmatrixSysSource(b_Ematrix_2d,"b_Ematrix_2d");
    
    
    
    for(Int_t i=1;i<=i_Ematrix->GetNbinsX();i++) 
    {
        double data = unfoldedData->GetBinContent(i);
        if(data>=0){
            t_Ematrix->SetBinContent(i, TMath::Sqrt(t_Ematrix_2d->GetBinContent(i,i))/data );
            i_Ematrix->SetBinContent(i, TMath::Sqrt(i_Ematrix_2d->GetBinContent(i,i))/data );
            d_Ematrix->SetBinContent(i, TMath::Sqrt(d_Ematrix_2d->GetBinContent(i,i))/data );
            if(a_Ematrix_2d->GetBinContent(i,i)>=0){
                a_Ematrix->SetBinContent(i, TMath::Sqrt(a_Ematrix_2d->GetBinContent(i,i))/data ); }
            //b_Ematrix->SetBinContent(i, TMath::Sqrt(b_Ematrix_2d->GetBinContent(i,i))/data );
        }
        
    }

    t_Ematrix->SetStats(0);
    t_Ematrix->GetYaxis()->SetRangeUser(-0.01,1);
    
    t_Ematrix->SetMarkerStyle(20);
    i_Ematrix->SetMarkerStyle(21);
    d_Ematrix->SetMarkerStyle(22);
    a_Ematrix->SetMarkerStyle(23);
    b_Ematrix->SetMarkerStyle(24);
    
    t_Ematrix->SetMarkerColor(1);
    i_Ematrix->SetMarkerColor(2);
    d_Ematrix->SetMarkerColor(3);
    a_Ematrix->SetMarkerColor(4);
    b_Ematrix->SetMarkerColor(5);
    
    t_Ematrix->SetMarkerSize(0.6);
    i_Ematrix->SetMarkerSize(0.6);
    d_Ematrix->SetMarkerSize(0.6);
    a_Ematrix->SetMarkerSize(0.6);
    b_Ematrix->SetMarkerSize(0.6);
    
    t_Ematrix->Draw("pe");
    i_Ematrix->Draw("same,pe");
    d_Ematrix->Draw("same,pe");
    a_Ematrix->Draw("same,pe");
    b_Ematrix->Draw("same,pe");
    
    legend->SetX1(0.1);
    legend->SetX2(0.4);
    legend->Clear();
    legend->AddEntry(t_Ematrix,"Total","pe");
    legend->AddEntry(i_Ematrix,"Input","pe");
    legend->AddEntry(d_Ematrix,"SysBackgroundUncorr","pe");
    legend->AddEntry(a_Ematrix,"SysUncorr","pe");
    //legend->AddEntry(b_Ematrix,"SysSource","pe");
    
    legend->Draw("same");
    writeCanvas(canvas,"ematricesInfo");
    legend->Clear();
    
    t_Ematrix->Delete();
    i_Ematrix->Delete();
    d_Ematrix->Delete();
    a_Ematrix->Delete();
    b_Ematrix->Delete();
    t_Ematrix_2d->Delete();
    i_Ematrix_2d->Delete();
    d_Ematrix_2d->Delete();
    a_Ematrix_2d->Delete();
    //b_Ematrix_2d->Delete();
    // ...
    
    
    //Draw input correlation matrix
    histMigration_->Draw("colz");
    writeCanvas(canvas,histMigration_->GetName());
    
    
    //Delete objects
    delete canvas;
    delete legend;
    ///ematrixTotal->Delete();
    ///ematrixInput->Delete();
    rhoIJtotal->Delete();
    probabilityMatrix->Delete();
    lCurve->Delete();
    bestLCurve->Delete();
    
    rhoIJtotalConv->Delete();
    errorIJtotal->Delete();
    unfoldedData->Delete();
    
    
    
}



void Plotter::writeCanvas( TCanvas* canvas, const TString name )
{
        if(nD_ == 2)canvas->Print(plotsFolder_+ name + "_" + v_plotName_.at(0) + "_vs_" + v_plotName_.at(1)  + ".pdf");
        if(nD_ == 1)canvas->Print(plotsFolder_+ name + "_" + v_plotName_.at(0) + ".pdf");
    canvas->Clear();
}


void Plotter::writePlotCP(const std::vector<Sample>& v_sample ,const std::vector<TH1D* >& v_SampleHist,const int& ind, const int binNum)
{
    int jnd = int((ind==0)?1:0);
    
    TCanvas* canvas = utils::setCanvas();
    TLegend* legend = utils::setLegend();
    canvas->cd();
    
    THStack* stack = new THStack();
    TH1D* histData = (TH1D*)v_SampleHist.at(0)->Clone("name");
    legend->AddEntry(histData,v_sample.at(0).legendEntry(),"pe");
    
    for(size_t iSample = 1; iSample < v_sample.size(); ++iSample)
    {   
        TH1D* hist = v_SampleHist.at(iSample);
        if(v_sample.at(iSample).sampleType() != Sample::data) stack->Add(hist);
        else histData->Add(hist);
        
        if(v_sample.at(iSample-1).legendEntry()!=v_sample.at(iSample).legendEntry())
            legend->AddEntry(hist,v_sample.at(iSample).legendEntry(),"f");
        
    }
    
    //Title
    TString plotTitle="";
    if(binNum == -1 ) plotTitle = v_plotTitle_.at(jnd) + " underflow bin" ;
    if(binNum == -2 ) plotTitle = v_plotTitle_.at(jnd) + " overflow bin" ;
    if(binNum >=0) plotTitle = utils::makeBinTitle(v_plotTitle_.at(jnd),v_coarseBins_.at(jnd).at(binNum),v_coarseBins_.at(jnd).at(binNum+1));
    if(binNum == -1 || binNum == -2)plotTitle = plotTitle + ";" + v_plotTitle_.at(ind) + ", " + v_plotUnits_.at(ind) + ";";
    else plotTitle = plotTitle  + ";" + v_plotTitle_.at(ind) + " " + v_plotUnits_.at(ind) + ";";
    
    //Draw
    TH1* stacksum = common::summedStackHisto(stack);
    histData->SetAxisRange(0,( histData->GetMaximum() > stacksum->GetMaximum() ? histData->GetMaximum()*1.2 : stacksum->GetMaximum()*1.2 ),"Y");
    histData->SetStats(0);
    histData->SetTitle(plotTitle);
    histData->Draw("e1");
    stack->Draw("same HIST");
    histData->Draw("e1 same");
    legend->Draw("same");
    UsefulTools::DrawCMSLabels(lumi_,2);
    UsefulTools::DrawDecayChLabel("emu");
    utils::drawRatio(histData, stacksum, NULL,0.5,2);
    
    //Print
    if(nD_==2&&binNum!=-999){
        if(binNum>=0)canvas->Print(plotsFolder_+ "CP_" + v_plotName_.at(ind) + "_IN_" + v_plotName_.at(jnd)  + "_" + std::to_string(binNum) + ".pdf");
        if(binNum==-1)canvas->Print(plotsFolder_+ "CP_" + v_plotName_.at(ind) + "_IN_" + v_plotName_.at(jnd)  + "_" + "u" + ".pdf");
        if(binNum==-2)canvas->Print(plotsFolder_+ "CP_" + v_plotName_.at(ind) + "_IN_" + v_plotName_.at(jnd)  + "_" + "o" + ".pdf");
    }
    if(nD_==1||binNum==-999)canvas->Print(plotsFolder_+ "CP_" + v_plotName_.at(ind) + ".pdf");
    
    //Delete
    stack->Delete();
    histData->Delete();
    delete canvas;
    delete legend;
    stacksum->Delete();
    
}

void Plotter::writePlotCPAllBins(const std::vector<Sample>& v_sample ,const std::vector<TH1* >& v_SampleHist)
{
    // Prepare canvas and legend
    TCanvas* canvas = utils::setCanvas();
    TLegend* legend = utils::setLegend();
    canvas->cd();
    
    THStack* stack = new THStack();
    TH1D* histData = (TH1D*)v_SampleHist.at(0)->Clone("name");
    legend->AddEntry(histData,v_sample.at(0).legendEntry(),"pe");
    for(size_t iSample = 1; iSample < v_sample.size(); ++iSample)
    {
        TH1D* hist = (TH1D*)v_SampleHist.at(iSample);
        if(v_sample.at(iSample).sampleType() != Sample::data) stack->Add(hist);
        else histData->Add(hist);
        
        if(v_sample.at(iSample-1).legendEntry()!=v_sample.at(iSample).legendEntry())
            legend->AddEntry(hist,v_sample.at(iSample).legendEntry(),"f");
        
    }
    TH1* stacksum = common::summedStackHisto(stack);
    double yMax = ( histData->GetMaximum() > stacksum->GetMaximum() ? histData->GetMaximum()*1.2 : stacksum->GetMaximum()*1.2 );
    
    int Nbins1 = (int)(v_coarseBins_.at(1).size()-1)+(int)(v_oTrue_.at(1)+v_uTrue_.at(1));
    int Nbins0 = (int)(v_coarseBins_.at(0).size()-1)+(int)(v_oTrue_.at(0)+v_uTrue_.at(0));
    
    TH1D hist("hist","",Nbins1,0.5,Nbins0*Nbins1+0.5);
    hist.SetTitle(";"+v_plotTitle_.at(1)+" "+v_plotUnits_.at(1));
    hist.SetStats(0);
    hist.GetXaxis()->SetTickLength(0);
    hist.SetAxisRange(0,yMax ,"Y");
    
    
        for(int i=1;i<=Nbins1;++i)
        {
            TString binName = "";
            int isUnderflow = v_uTrue_.at(1);
            int isOverflow =  v_oTrue_.at(1);
            if(isUnderflow&&i==1) binName = "underflow";
            else if(isOverflow&&i==Nbins1){
                 binName = "overflow";
            }
            else{
                binName = utils::makeBinTitle("",v_coarseBins_.at(1).at(i-isUnderflow-1),v_coarseBins_.at(1).at(i-isUnderflow));
            }
            hist.GetXaxis()->SetBinLabel(i,binName.Data());
        }
        hist.SetTitle( utils::makeTitleBins(v_plotTitle_.at(0) +" "+ v_plotUnits_.at(0),v_coarseBins_.at(0),v_uTrue_.at(0),v_oTrue_.at(0) ) );
    hist.Draw();
    
    histData->SetAxisRange(0,( histData->GetMaximum() > stacksum->GetMaximum() ? histData->GetMaximum()*1.2 : stacksum->GetMaximum()*1.2 ),"Y");
    histData->SetStats(0);
    stack->Draw("same HIST");
    histData->Draw("e1 same");
    for(int i=1;i<Nbins1;++i)
        {
            double xLine = i*(Nbins0) + 0.5;
            TLine ln(xLine,-0.03*yMax,xLine,yMax);//FIXME: y1 to be related to pad size... 
            ln.DrawClone("same");

        }
    legend->Draw("same");
    utils::drawRatio(histData, stacksum, NULL,0.5,2,(TH1D*)hist.Clone());
    
    canvas->Print(plotsFolder_+ "CP_AllBins_" + v_plotName_.at(0) + "_vs_"+ v_plotName_.at(1) + ".pdf");
    
    //Delete objects
    delete canvas;
    delete legend;
    stacksum->Delete();
}

void Plotter::writePlotEPSAllBins()
{
    // Prepare canvas and legend
    TCanvas* canvas = utils::setCanvas();
    TLegend* legend = utils::setLegend();
    legend->SetX1(0.1);
    legend->SetX2(0.4);
    legend->SetY1(0.7);
    legend->SetY2(0.85);
    
    canvas->cd();
    
    histEffAllBins_->Divide(histEffAllBins_,histGenAllBins_,1,1,"B");
    histEffAllBins_->SetStats(0);
    histEffAllBins_->SetMarkerStyle(20);
    histEffAllBins_->SetMarkerColor(3);
    legend->AddEntry(histEffAllBins_,"Efficiency","pe");
    
    histPurityAllBins_->Divide(histRecoGenAllBins_,histPurityAllBins_,1,1,"B");
    histPurityAllBins_->SetStats(0);
    histPurityAllBins_->SetMarkerStyle(23);
    histPurityAllBins_->SetMarkerColor(4);
    legend->AddEntry(histPurityAllBins_,"Purity","pe");
    
    histStabilityAllBins_->Divide(histRecoGenAllBins_,histStabilityAllBins_,1,1,"B");
    histStabilityAllBins_->SetStats(0);
    histStabilityAllBins_->SetMarkerStyle(22);
    histStabilityAllBins_->SetMarkerColor(2);
    legend->AddEntry(histStabilityAllBins_,"Stability","pe");
    
    
    double yMax = histPurityAllBins_->GetMaximum() > histStabilityAllBins_->GetMaximum() ? histPurityAllBins_->GetMaximum()*1.15 : histStabilityAllBins_->GetMaximum()*1.15;
    histPurityAllBins_->SetAxisRange(0,yMax ,"Y");
    
    int Nbins1 = (int)(v_coarseBins_.at(1).size()-1)+(int)(v_oTrue_.at(1)+v_uTrue_.at(1));
    int Nbins0 = (int)(v_coarseBins_.at(0).size()-1)+(int)(v_oTrue_.at(0)+v_uTrue_.at(0));
    
    
    
    
    TH1D hist("hist","",Nbins1,0.5,Nbins0*Nbins1+0.5);
    hist.SetTitle(";"+v_plotTitle_.at(1)+" "+v_plotUnits_.at(1));
    hist.SetStats(0);
    hist.GetXaxis()->SetTickLength(0);
    hist.SetAxisRange(0,yMax ,"Y");
    hist.SetAxisRange(0,1 ,"Y");
    
    
    if(nD_ == 2){
        for(int i=1;i<=Nbins1;++i)
        {
            TString binName = "";
            int isUnderflow = v_uTrue_.at(1);
            int isOverflow =  v_oTrue_.at(1);
            if(isUnderflow&&i==1) binName = "underflow";
            else if(isOverflow&&i==Nbins1){
                 binName = "overflow";
            }
            else{
                binName = utils::makeBinTitle("",v_coarseBins_.at(1).at(i-isUnderflow-1),v_coarseBins_.at(1).at(i-isUnderflow));
            }
            
            hist.GetXaxis()->SetBinLabel(i,binName.Data());
        }
        hist.SetTitle( utils::makeTitleBins(v_plotTitle_.at(0) +" "+ v_plotUnits_.at(0),v_coarseBins_.at(0),v_uTrue_.at(0),v_oTrue_.at(0) ) );
        
        
        
    hist.Draw();
    
    
        for(int i=1;i<Nbins1;++i)
        {
            double xLine = i*(Nbins0) + 0.5;
            TLine ln(xLine,-0.03,xLine,yMax);//FIXME: y1 to be related to pad size... 
            ln.DrawClone("same");

        }
    }
    
    
    histPurityAllBins_->Draw("e same");
    histStabilityAllBins_->Draw("e same");
    histEffAllBins_->Draw("e same");
    legend->Draw("same");
    if(nD_ == 1)canvas->Print(plotsFolder_+ "EPS_AllBins_" + v_plotName_.at(0) + ".pdf");
    if(nD_ == 2)canvas->Print(plotsFolder_+ "EPS_AllBins_" + v_plotName_.at(0) + "_vs_"+ v_plotName_.at(1) + ".pdf");
    
    //Delete objects
    delete canvas;
    delete legend;
    
}



void Plotter::writePlotEPS()
{
    // Prepare canvas and legend
    TCanvas* canvas = utils::setCanvas();
    TLegend* legend = utils::setLegend();
    canvas->cd();
    
    histPurity_->Divide(histRecoGen_,histPurity_,1,1,"B");
    histPurity_->SetMarkerStyle(23);
    histPurity_->SetMarkerColor(4);
    
    histStability_->Divide(histRecoGen_,histStability_,1,1,"B");
    histStability_->SetMarkerStyle(22);
    histStability_->SetMarkerColor(2);
    
    legend->AddEntry(histPurity_,"Purity","pe");
    legend->AddEntry(histStability_,"Stability","pe");
    
    
    if(nD_==1){
          
          histPurity_->SetAxisRange(0,histPurity_->GetMaximum() > histStability_->GetMaximum() ? histPurity_->GetMaximum()*1.15 : histStability_->GetMaximum()*1.15 ,"Y");
          histPurity_->Draw("e");
          histStability_->Draw("e same");
          //legend->Draw("same");
          canvas->Print(plotsFolder_+ "EPS_" + v_plotName_.at(0) + ".pdf");
    }
    
        if(nD_==2){
        TH2* h2dPurity = (TH2*)histPurity_->Clone();
        TH2* h2dStability = (TH2*)histStability_->Clone();
        
        for(int iy=-1;iy<=(int)v_coarseBins_.at(1).size()-1;iy++)
        {
            TString plotTitle = "";
            
            if(iy==-1)
            {
                if(v_uTrue_.at(1))plotTitle = v_plotTitle_.at(1) + " underflow bin";
                if(!v_uTrue_.at(1))continue;
            }
            else if(iy == (int)v_coarseBins_.at(1).size()-1)
            {
                if(v_oTrue_.at(1))plotTitle = v_plotTitle_.at(1) + " overflow bin";
                if(!v_oTrue_.at(1))continue;
            }
            else plotTitle = utils::makeBinTitle(v_plotTitle_.at(1),v_coarseBins_.at(1).at(iy),v_coarseBins_.at(1).at(iy+1));
            plotTitle = plotTitle + ";" + v_plotTitle_.at(0) + ", " + v_plotUnits_.at(0);
            
            TH1* h_tempP = (TH1*)(h2dPurity->ProjectionX("tempP_iy",iy+1,iy+1,"e"));
            h_tempP->SetStats(0);
            h_tempP->SetTitle(plotTitle);
            
            TH1* h_tempS = (TH1*)(h2dStability->ProjectionX("tempS_iy",iy+1,iy+1,"e"));
            h_tempS->SetTitle(plotTitle);
            
            h_tempP->SetAxisRange(0,h_tempP->GetMaximum() > h_tempS->GetMaximum() ? h_tempP->GetMaximum()*1.15 : h_tempS->GetMaximum()*1.15 ,"Y");
            h_tempP->Draw("e");
            h_tempS->Draw("e same");
            //legend->Draw("same");
            
            if(iy==-1)canvas->Print(plotsFolder_+ "EPS_" + v_plotName_.at(0) + "_IN_"+ v_plotName_.at(1)+ "_" + "u" + ".pdf");
            else if(iy == (int)v_coarseBins_.at(1).size()-1)canvas->Print(plotsFolder_+ "EPS_" + v_plotName_.at(0) + "_IN_"+ v_plotName_.at(1)+ "_" + "o" + ".pdf");
            else canvas->Print(plotsFolder_+ "EPS_" + v_plotName_.at(0) + "_IN_"+ v_plotName_.at(1)+ "_" + std::to_string(iy) + ".pdf");
            canvas->Clear();
            
            h_tempP->Delete();
            h_tempS->Delete();
        }
        

        for(int ix=-1;ix<=(int)v_coarseBins_.at(0).size()-1;ix++)
        {
            TString plotTitle = "";
            
            if(ix==-1)
            {
                if(v_uTrue_.at(0))plotTitle = v_plotTitle_.at(0) + " underflow bin";
                if(!v_uTrue_.at(0))continue;
            }
            else if(ix == (int)v_coarseBins_.at(0).size()-1)
            {
                if(v_oTrue_.at(0))plotTitle = v_plotTitle_.at(0) + " overflow bin";
                if(!v_oTrue_.at(0))continue;
            }
            else plotTitle = utils::makeBinTitle(v_plotTitle_.at(0),v_coarseBins_.at(0).at(ix),v_coarseBins_.at(0).at(ix+1));
            plotTitle = plotTitle + ";" + v_plotTitle_.at(1) + ", " + v_plotUnits_.at(1);
            
            TH1* h_tempP = (TH1*)(h2dPurity->ProjectionY("tempP_ix",ix+1,ix+1,"e"));
            h_tempP->SetStats(0);
            h_tempP->SetTitle(plotTitle);
            
            TH1* h_tempS = (TH1*)(h2dStability->ProjectionY("tempS_ix",ix+1,ix+1,"e"));
            h_tempS->SetTitle(plotTitle);
            
            h_tempP->SetAxisRange(0,h_tempP->GetMaximum() > h_tempS->GetMaximum() ? h_tempP->GetMaximum()*1.15 : h_tempS->GetMaximum()*1.15 ,"Y");
            h_tempP->Draw("e");
            h_tempS->Draw("e same");
            //legend->Draw("same");
            
           if(ix==-1)canvas->Print(plotsFolder_+ "EPS_" + v_plotName_.at(1) + "_IN_"+ v_plotName_.at(0)+ "_" + "u" + ".pdf");
           else if(ix == (int)v_coarseBins_.at(0).size()-1)canvas->Print(plotsFolder_+ "EPS_" + v_plotName_.at(1) + "_IN_"+ v_plotName_.at(0)+ "_" + "o" + ".pdf");
           else canvas->Print(plotsFolder_+ "EPS_" + v_plotName_.at(1) + "_IN_"+ v_plotName_.at(0)+ "_" + std::to_string(ix) + ".pdf");
           canvas->Clear();
                
            h_tempP->Delete();
            h_tempS->Delete();
            
        }

        h2dPurity->Delete();
        h2dStability->Delete();
      
     }


    //Delete objects
    delete canvas;
    delete legend;
    
}



void Plotter::writePlotXSec(const TH1* hData,const TH1* hMC)
{
    // Prepare canvas and legend
    TCanvas* canvas = utils::setCanvas();
    TLegend* legend = utils::setLegend();
    canvas->cd();
    
         TH1* histData = (TH1*)hData->Clone("data");
         histData->Scale(1.0/lumi_);
         histData->Scale(1.0/topxsec_);
         histData->Scale(1.0/v_BR_.at(Channel::emu));//FIXME: for different channels
         histData->SetStats(0);

         TH1* histMC   = (TH1*)hMC->Clone("mc");
         histMC->Scale(1.0/lumi_);
         histMC->Scale(1.0/topxsec_);
         histMC->Scale(1.0/v_BR_.at(Channel::emu));//FIXME: for different channels
         histMC->SetLineColor(2);
         

    if(nD_==1){
          histData->Scale(1.,"width");
          histMC->Scale(1.,"width");
          
          histData->SetAxisRange(0,( histData->GetMaximum() > histMC->GetMaximum() ? histData->GetMaximum()*1.15 : histMC->GetMaximum()*1.15 ),"Y");
          histData->Draw("e");
          histMC->Draw("hist same");
          histData->Draw("e,same");
          utils::drawRatio(histData, histMC, NULL,0.5,2);
          canvas->Print(plotsFolder_+ "xSec_" + v_plotName_.at(0) + ".pdf");
    }
    
    if(nD_==2){
        TH2* h2dData = (TH2*)histData->Clone();
        TH2* h2dMC = (TH2*)histMC->Clone();
        for(int iy=-1;iy<=(int)v_coarseBins_.at(1).size()-1;iy++)
        {
            TString plotTitle = "";
            
            if(iy==-1)
            {
                if(v_uTrue_.at(1))plotTitle = v_plotTitle_.at(1) + " underflow bin";
                if(!v_uTrue_.at(1))continue;
            }
            else if(iy == (int)v_coarseBins_.at(1).size()-1)
            {
                if(v_oTrue_.at(1))plotTitle = v_plotTitle_.at(1) + " overflow bin";
                if(!v_oTrue_.at(1))continue;
            }
            else plotTitle = utils::makeBinTitle(v_plotTitle_.at(1),v_coarseBins_.at(1).at(iy),v_coarseBins_.at(1).at(iy+1));
            
            plotTitle = plotTitle + ";" + v_plotTitle_.at(0) + ", " + v_plotUnits_.at(0);
            
            TH1* h_temp = (TH1*)(h2dData->ProjectionX("tempData_iy",iy+1,iy+1,"e"));
            h_temp->Scale(1,"width");
            TH1* h_temp_mc =  (TH1*)(h2dMC->ProjectionX("tempMC_iy",iy+1,iy+1,"e"));
            h_temp_mc->Scale(1,"width");
            h_temp->SetStats(0);
            h_temp->SetTitle(plotTitle);
            
                h_temp->SetAxisRange(0,( h_temp->GetMaximum() > h_temp_mc->GetMaximum() ? h_temp->GetMaximum()*1.15 : h_temp_mc->GetMaximum()*1.15 ),"Y");
                h_temp->Draw("e");
                h_temp_mc->Draw("hist same");
                h_temp->Draw("e,same");
                utils::drawRatio(h_temp, h_temp_mc, NULL,0.5,2);
                
                if(iy==-1)canvas->Print(plotsFolder_+ "xSec_" + v_plotName_.at(0) + "_IN_"+ v_plotName_.at(1)+ "_" + "u" + ".pdf");
                else if(iy == (int)v_coarseBins_.at(1).size()-1)canvas->Print(plotsFolder_+ "xSec_" + v_plotName_.at(0) + "_IN_"+ v_plotName_.at(1)+ "_" + "o" + ".pdf");
                else{
                    
                    Output xsecInfo("xsec");
                    std::vector<double> v_bin = {v_coarseBins_.at(1).at(iy),v_coarseBins_.at(1).at(iy+1)};
                        xsecInfo.add("bin",v_bin);
                        xsecInfo.add("bins",v_coarseBins_.at(0));
                    std::vector<double> v_topxsec = {topxsec_};
                        xsecInfo.add("inc.xsec",v_topxsec);
                    std::vector<double> v_xsec;
                    std::vector<double> v_stat;
                    std::vector<double> v_xsec_mc;
                    for(int i=1;i<=h_temp->GetNbinsX();i++){
                        v_xsec.push_back(h_temp->GetBinContent(i));
                        v_stat.push_back(h_temp->GetBinError(i));
                        v_xsec_mc.push_back(h_temp_mc->GetBinContent(i));
                    }
                    xsecInfo.add("stat",v_stat);
                    xsecInfo.add("xsec",v_xsec);
                    xsecInfo.add("mc",v_xsec_mc);
                    xsecInfo.save(plotsFolder_+ "xSec_" + v_plotName_.at(0) + "_IN_"+ v_plotName_.at(1)+ "_" + std::to_string(iy) + ".txt");
                    
                    canvas->Print(plotsFolder_+ "xSec_" + v_plotName_.at(0) + "_IN_"+ v_plotName_.at(1)+ "_" + std::to_string(iy) + ".pdf");
                }
                canvas->Clear();
                h_temp->Delete();
                h_temp_mc->Delete();
        }
        

        for(int ix=-1;ix<=(int)v_coarseBins_.at(0).size()-1;ix++)
        {
            TString plotTitle = "";
            
            if(ix==-1)
            {
                if(v_uTrue_.at(0))plotTitle = v_plotTitle_.at(0) + " underflow bin";
                if(!v_uTrue_.at(0))continue;
            }
            else if(ix == (int)v_coarseBins_.at(0).size()-1)
            {
                if(v_oTrue_.at(0))plotTitle = v_plotTitle_.at(0) + " overflow bin";
                if(!v_oTrue_.at(0))continue;
            }
            else plotTitle = utils::makeBinTitle(v_plotTitle_.at(0),v_coarseBins_.at(0).at(ix),v_coarseBins_.at(0).at(ix+1));
            
            plotTitle = plotTitle + ";" + v_plotTitle_.at(1) + ", " + v_plotUnits_.at(1);
            
            TH1* h_temp = (TH1*)(h2dData->ProjectionY("tempData_ix",ix+1,ix+1,"e"));
            h_temp->Scale(1,"width");
            TH1* h_temp_mc =  (TH1*)(h2dMC->ProjectionY("tempMC_ix",ix+1,ix+1,"e"));
            h_temp_mc->Scale(1,"width");
            h_temp->SetStats(0);
            h_temp->SetTitle(plotTitle);
            
                h_temp->SetAxisRange(0,( h_temp->GetMaximum() > h_temp_mc->GetMaximum() ? h_temp->GetMaximum()*1.15 : h_temp_mc->GetMaximum()*1.15 ),"Y");
                h_temp->Draw("e");
                h_temp_mc->Draw("hist same");
                h_temp->Draw("e,same");
//                 legend->AddEntry(h_temp,"Data","p");
//                 legend->AddEntry(h_temp_mc,"MadGraph+Pythia","l");
//                 legend->DrawClone("same");
                utils::drawRatio(h_temp, h_temp_mc, NULL,0.5,2);
                
                if(ix==-1)canvas->Print(plotsFolder_+ "xSec_" + v_plotName_.at(1) + "_IN_"+ v_plotName_.at(0)+ "_" + "u" + ".pdf");
                else if(ix == (int)v_coarseBins_.at(0).size()-1)canvas->Print(plotsFolder_+ "xSec_" + v_plotName_.at(1) + "_IN_"+ v_plotName_.at(0)+ "_" + "o" + ".pdf");
                else{
                    
                    Output xsecInfo("xsec");
                    std::vector<double> v_bin = {v_coarseBins_.at(0).at(ix),v_coarseBins_.at(0).at(ix+1)};
                        xsecInfo.add("bin",v_bin);
                        xsecInfo.add("bins",v_coarseBins_.at(1));
                    std::vector<double> v_topxsec = {topxsec_};
                        xsecInfo.add("inc.xsec",v_topxsec);
                    std::vector<double> v_xsec;
                    std::vector<double> v_stat;
                    std::vector<double> v_xsec_mc;
                    for(int i=1;i<=h_temp->GetNbinsX();i++){
                        v_xsec.push_back(h_temp->GetBinContent(i));
                        v_stat.push_back(h_temp->GetBinError(i));
                        v_xsec_mc.push_back(h_temp_mc->GetBinContent(i));
                    }
                    xsecInfo.add("stat",v_stat);
                    xsecInfo.add("xsec",v_xsec);
                    xsecInfo.add("mc",v_xsec_mc);
                    xsecInfo.save(plotsFolder_+ "xSec_" + v_plotName_.at(1) + "_IN_"+ v_plotName_.at(0)+ "_" + std::to_string(ix) + ".txt");
                    
                    canvas->Print(plotsFolder_+ "xSec_" + v_plotName_.at(1) + "_IN_"+ v_plotName_.at(0)+ "_" + std::to_string(ix) + ".pdf");
                }
                canvas->Clear();
                h_temp->Delete();
                h_temp_mc->Delete();
        }

      h2dData->Delete();
      h2dMC->Delete();
     }

          histMC->Delete();
          histData->Delete();

    //Delete objects
    delete canvas;
    delete legend;
}



const SystematicChannelFactors& Plotter::scaleFactors(const TString& name)
{
    const TString fullStepname = ttbar::extractSelectionStepAndJetCategory(name);
    
    // Check if map contains already scale factors for this step, else access and fill them
    if(m_stepFactors_.find(fullStepname) == m_stepFactors_.end()){
        m_stepFactors_[fullStepname] = samples_.globalWeights(fullStepname).first;
    }
    
    return m_stepFactors_.at(fullStepname);
}

