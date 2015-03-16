
#include <fstream>
#include <iostream>
#include <sstream>

#include <TCanvas.h>
#include <TLegend.h>
#include <TStyle.h>
#include <TSystem.h>
#include <THStack.h>
#include <TDirectory.h>
#include <TFile.h>
#include <TString.h>
#include <TTree.h>
#include <TH1.h>
#include <TH2.h>
#include <TGraph.h>
#include <TLine.h>
#include <TMath.h>

#include "utils.h"
#include "PlotterBTop.h"
#include "ttbarUtils.h"
#include "Samples.h"
#include "Sample.h"
#include "../../common/include/utils.h"
#include "../../common/include/plotterUtils.h"
#include "UsefulTools.h"
#include "Output.h"

#include "TUnfold.h"
#include "TUnfoldBinning.h"
#include "TUnfoldDensity.h"





PlotterBTop::PlotterBTop(const Samples& samples, const double lumi, const double topxsec, const TString& specname):

samples_(samples),
lumi_(lumi),
topxsec_(topxsec),
specname_(specname),
plotName_(""),
plotTitle_(""),
plotUnits_(""),
plotsFolder_("Plots/"),

//Control plots//
cpNBins_(-999),
R1_(-999),
R2_(-999),
v_SampleHist_(std::vector<TH1D* >(0)),
v_UnderflowSampleHist_(std::vector<TH1D* >(0)),
v_OverflowSampleHist_(std::vector<TH1D* >(0)),
hRecoGen_(0),
recogenBins_(std::vector<Double_t>(0)),


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
// histRWSG_(0),
// 
histGen_(0),
histPurity_(0),
histStability_(0),
histRecoGen_(0),
histPurityAllBins_(0),
histStabilityAllBins_(0),
histRecoGenAllBins_(0),
histEffAllBins_(0),
histGenAllBins_(0),

coarseBins_(std::vector<Double_t>(0)),
fineBins_(std::vector<Double_t>(0)),
uReco_(-999),
oReco_(-999),
uTrue_(-999),
oTrue_(-999),

// Tree branches //
entry_(-999),
entry0_(-999),
eventWeight_(-999),
trueLevelWeight_(-999),
trueLevelWeight0_(-999),
top_pt_(-999),
topbar_pt_(-999),
branchVal_(-999),
branchValGen_(-999),
branchValGen0_(-999)/*,


histMCReco_(0)*/



{
    //FIXME: make it in better way
    v_BR_.push_back(0.01166);//ee
    v_BR_.push_back(0.02332);//emu
    v_BR_.push_back(0.01166);//mumu
    v_BR_.push_back(0.04666);//combined
     
     // switch on histogram errors
    TH1::SetDefaultSumw2();
     
}



void PlotterBTop::setOptions(const std::vector<TString> v_plotName)
{
    
    //std::cout<<"[PlotterBTop]:--- Beginning of plot options settings\n\n";
    
    //Clearing previous declarations  
    this->clearMemory();
    
    plotName_ = v_plotName.at(0);
    TString nDflag;
    nDflag="1d";

    branchVal_ = -999;
    branchValGen_ = -999;
    branchValGen0_ = -999;
    
    // reading Control Card
        const std::string ccFileName(common::CMSSW_BASE() + "/src/TopAnalysis/Configuration/analysis/diLeptonic/ControlCardsBTop/" + plotName_ + ".txt");
        std::ifstream ccFileStream(ccFileName.data(), std::ifstream::in);
        if (!ccFileStream.good()) {
          std::cerr<<"Error in PlotterBTop! Cannot find file with name: "<< ccFileName <<"\n...break\n"<<std::endl;
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
                     cpNBins_ = (vWord.at(1)).Atoi();
                     R1_ = (vWord.at(2)).Atof();
                     R2_ = (vWord.at(3)).Atof();
             }
             if(vWord.at(0) == "title")plotTitle_ = vWord.at(1);
             
             if(vWord.at(0) == "units" ){
                if(vWord.size()>1) plotUnits_ = vWord.at(1); 
                if(vWord.size()<2) plotUnits_ = "";
             }
             
             if(vWord.at(0) == "recogen" ){
                for(auto word : vWord)recogenBins_.push_back(word.Atof());
                recogenBins_.erase(recogenBins_.begin(),recogenBins_.begin()+1);
             }
             
             //Pass to dimensional block
             if( (int)vWord.size()==1 && vWord.at(0) == nDflag ){
                 if(isInBlock) isInBlock = false;
                 if(!isInBlock) isInBlock = true;
            }
             
             // Reading only from dimensional block we are interesting
             if(isInBlock){

                if(vWord.at(0) == "coarse"){
                     for(auto word : vWord)coarseBins_.push_back(word.Atof());
                     coarseBins_.erase(coarseBins_.begin(),coarseBins_.begin()+1);
                }
                if(vWord.at(0) == "fine"){
                     for(auto word : vWord)fineBins_.push_back(word.Atof());
                     fineBins_.erase(fineBins_.begin(),fineBins_.begin()+1);
                }
                if(vWord.at(0) == "uoTrue")
                {
                    uTrue_ = vWord.at(1).Atoi();
                    oTrue_ = vWord.at(2).Atoi();
                }
                if(vWord.at(0) == "uoReco")
                {
                    uReco_ = vWord.at(1).Atoi();
                    oReco_ = vWord.at(2).Atoi();
                }
                 
             }//read from dimensional block
             else continue;
             
       }//stream line

    //Set binning scheme
    //Definition of TUnfoldBinning no detector and generator level
    detectorBinning_ = new TUnfoldBinning("detector");
    detectorDistribution_ = detectorBinning_->AddBinning("detectordistribution");
    detectorDistribution_->AddAxis(plotName_,(int)fineBins_.size()-1,fineBins_.data(),uReco_, oReco_);
    
    
    generatorBinning_ = new TUnfoldBinning("generator");
    generatorDistribution_ = generatorBinning_->AddBinning("signal");
    generatorDistribution_->AddAxis("gen_"+plotName_,(int)coarseBins_.size()-1, coarseBins_.data() , uTrue_, oTrue_);
    
    //std::cout<<"[PlotterBTop]:--- Finishing of plot options settings\n\n";
}



void PlotterBTop::clearMemory()
{
    plotName_.Clear();
    plotTitle_.Clear();
    plotUnits_.Clear();
    branchVal_ = -999;
    branchValGen_ = -999;
    branchValGen0_ = -999;
    coarseBins_.clear();
    fineBins_.clear();
    recogenBins_.clear();
    uReco_ = -999;
    oReco_ = -999;
    uTrue_ = -999;
    oTrue_ = -999;
    cpNBins_ = -999;
    R1_ = -999;
    R2_ = -999;
    
    if(detectorDistribution_)detectorDistribution_->Delete();
    if(detectorBinning_)detectorBinning_->Delete();
    if(generatorDistribution_)generatorDistribution_->Delete();
    if(generatorBinning_)generatorBinning_->Delete();
    
}

void PlotterBTop::clearMemoryPerSystematicChannel()
{
     for(auto p : v_SampleHist_)
         delete p;
    v_SampleHist_.clear();
    
    for(auto p : v_UnderflowSampleHist_)
         delete p;
    v_UnderflowSampleHist_.clear();
    
    for(auto p : v_OverflowSampleHist_)
         delete p;
    v_OverflowSampleHist_.clear();
    
    if(hRecoGen_)hRecoGen_->Delete();
    
    if(histMigration_)histMigration_->Delete();
    if(histBgrUo_)histBgrUo_->Delete();
    if(histBgr_)histBgr_->Delete();
    if(histData_)histData_->Delete();
    if(unfoldedData_)unfoldedData_->Delete();
    
    //Closurer test//
    if(histSR_)histSR_->Delete();
    if(histUf_)histUf_->Delete();
    if(histSG_)histSG_->Delete();
//     if(histRWSG_)histRWSG_->Delete();
    
//     if(histMCReco_)histMCReco_->Delete();
    
    
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



void PlotterBTop::prepareHistograms(const std::vector<Sample>& v_sample)
{
        
        for(size_t iSample = 0; iSample < v_sample.size(); ++iSample){
            const auto& sample(v_sample.at(iSample));
            TH1D* sampleHist = new TH1D(plotName_+"_cp" + std::to_string((int)iSample)  ,plotTitle_,cpNBins_,R1_,R2_);//FIXME: remuve this: + std::to_string(ind) + std::to_string((int)iSample)
            sampleHist->SetFillColor(sample.color());
            sampleHist->SetLineColor(sample.color());
            v_SampleHist_.push_back(sampleHist);
            v_UnderflowSampleHist_.push_back((TH1D*)(sampleHist->Clone(sampleHist->GetName() + TString("u"))));
            v_OverflowSampleHist_.push_back((TH1D*)(sampleHist->Clone(sampleHist->GetName() + TString("o"))));
        }
        
        hRecoGen_ = new TH2D(plotName_+"_recogen",specname_+";reco "+plotTitle_+plotUnits_+";gen "+plotTitle_+plotUnits_,
                             recogenBins_.at(0),recogenBins_.at(1),recogenBins_.at(2),recogenBins_.at(3),recogenBins_.at(4),recogenBins_.at(5));
    
    //Definition of unfolding histograms
    histMigration_ = TUnfoldBinning::CreateHistogramOfMigrations(generatorBinning_,detectorBinning_,"histMigration");
    histBgrUo_ = detectorBinning_->CreateHistogram("histBgrUo");
    histBgr_ = detectorBinning_->CreateHistogram("histBgr");
    histData_ = detectorBinning_->CreateHistogram("histDataR");
    
//     histMCReco_ = detectorBinning_->CreateHistogram("histDataReco");
//     
//     
        histGen_ = generatorBinning_->CreateHistogram("histGen",kTRUE);
        histPurity_ = generatorBinning_->CreateHistogram("histPurity",kTRUE);
        histStability_ = generatorBinning_->CreateHistogram("histStability",kTRUE);
        histRecoGen_ = generatorBinning_->CreateHistogram("histRecoGen",kTRUE);
        
        histEffAllBins_ = generatorBinning_->CreateHistogram("histEffAllBins");
        histGenAllBins_ = generatorBinning_->CreateHistogram("histGenAllBins");
        histPurityAllBins_ = generatorBinning_->CreateHistogram("histPurityAllBins");
        histStabilityAllBins_ = generatorBinning_->CreateHistogram("histStabilityAllBins");
        histRecoGenAllBins_ = generatorBinning_->CreateHistogram("histRecoGenAllBins");
    
//     //Closure test//
    histSR_ = detectorBinning_->CreateHistogram("histSR");
    histSG_ = generatorBinning_->CreateHistogram("histSG");
//     histRWSG_ = generatorBinning_->CreateHistogram("histRWSG");
        
}



void PlotterBTop::producePlots()
{
    //std::cout<<"[PlotterBTop]:--- Beginning of plot production\n\n";
    
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
            plotsFolder_ = common::assignFolder("Plots",channel,systematic,plotName_);
            this->prepareHistograms(v_sample);
            
            TString sampleInfoFolder = common::assignFolder("Plots",channel,systematic) + "/sampleInfo.txt";
            Output sampleInfo("sample");
            
            for(size_t iSample = 0; iSample < v_sample.size(); ++iSample){                  //sample loop
                const auto& sample(v_sample.at(iSample));
                
                if(sample.legendEntry() == "QCD Multijet") continue;
                if(sample.legendEntry() == "W+Jets") continue;
                
                TFile* dataFile=new TFile(sample.inputFile().ReplaceAll("selectionRoot","btopInput"));
                //TTree *dataTree0=(TTree *) dataFile->Get("bTop_treeVariables_step0"); 
                TTree *dataTree8=(TTree *) dataFile->Get("bTop_treeVariables_step8");
                if(/*!dataTree0 ||*/ !dataTree8) {
                     std::cout<<"could not read 'bTop_treeVariables_step' tree from " << dataFile->GetName() << "\n";
                 }
                 
                // set branches
//                 dataTree0->ResetBranchAddresses();
//                 dataTree0->SetBranchAddress("entry",&entry0_);
//                 dataTree0->SetBranchAddress("trueLevelWeight",&trueLevelWeight0_);
//                 dataTree0->SetBranchAddress("gen_"+plotName_,&branchValGen0_);
                
                dataTree8->ResetBranchAddresses();
                dataTree8->SetBranchAddress("entry",&entry_);
                dataTree8->SetBranchAddress("eventWeight",&eventWeight_);
                dataTree8->SetBranchAddress("trueLevelWeight",&trueLevelWeight_);
                dataTree8->SetBranchAddress("isTopGen",&isTopGen_);
                dataTree8->SetBranchAddress("isKinReco",&isKinReco_);
                dataTree8->SetBranchAddress(plotName_,&branchVal_);
                dataTree8->SetBranchAddress("gen_"+plotName_,&branchValGen_);
                dataTree8->SetBranchAddress("top_pt",&top_pt_);
                dataTree8->SetBranchAddress("topbar_pt",&topbar_pt_);
                
                sampleInfo.add("entries", utils::numToString(dataTree8->GetEntriesFast()));
                sampleInfo.add("weight",utils::numToString(v_weight.at(iSample)));
                sampleInfo.add("file",sample.inputFile());
                sampleInfo.add("name",sample.legendEntry());
                
                //loop over tree8 events
                Int_t lastEvent=0;
                for(Int_t ievent=0;ievent<dataTree8->GetEntriesFast();ievent++){
                    if(dataTree8->GetEntry(ievent)<=0) break;

                    if(plotName_=="top_pt")branchVal_=top_pt_;
                    if(top_pt_ < 300 || topbar_pt_ < 300)continue;
                    //Control plots//
                        v_SampleHist_.at(iSample)->Fill(branchVal_,eventWeight_*v_weight.at(iSample));
                        v_UnderflowSampleHist_.at(iSample)->Fill(branchVal_,eventWeight_*v_weight.at(iSample));
                        v_OverflowSampleHist_.at(iSample)->Fill(branchVal_,eventWeight_*v_weight.at(iSample));
                        if(isTopGen_&&plotName_!="mlblbmet")hRecoGen_->Fill(branchVal_,branchValGen_,eventWeight_*v_weight.at(iSample));
                    // ... //
                    
                    if(isTopGen_){
//                         for(Int_t ievent0=lastEvent;ievent0<dataTree0->GetEntriesFast();ievent0++) {
//                             if(dataTree0->GetEntry(ievent0)<=0) break;
//                                 
//                                 histGen_->Fill(branchValGen0_,trueLevelWeight0_*v_weight.at(iSample));
//                                 
//                                 histSG_->Fill(genBin_(branchValGen0_),trueLevelWeight0_*v_weight.at(iSample));
//                                 //histSG_->Fill(genBin_(branchValGen0_),1*v_weight.at(iSample));
//                                 //Closure test
//                                 //histRWSG_->Fill(genBin_(branchValsGen0_),(rewTopPtUp(branchValsGen0_.at(1)))*v_weight.at(iSample));//pt
//                                 //histRWSG_->Fill(genBin_(branchValsGen0_),(rewTopPtDown(branchValsGen0_.at(1)))*v_weight.at(iSample));//pt down
//                                 //histRWSG_->Fill(genBin_(branchValsGen0_),(rewTopRapidityUp(branchValsGen0_.at(0)))*v_weight.at(iSample));//y up
//                                 //histRWSG_->Fill(genBin_(branchValsGen0_),(rewTopRapidityDown(branchValsGen0_.at(0)))*v_weight.at(iSample));//y down
//                                 //histRWSG_->Fill(genBin_(branchValsGen0_),1*v_weight.at(iSample));//1
//                                 
//                                 histGenAllBins_->Fill(genBin_(branchValGen0_),trueLevelWeight0_*v_weight.at(iSample));
//                                 
//                                 if(entry0_ != entry_)
//                                 {
//                                     histMigration_->Fill(genBin_(branchValGen0_),0.,trueLevelWeight0_*v_weight.at(iSample));
//                                     //histMigration_->Fill(genBin_(branchValsGen0_),0.,1*v_weight.at(iSample));
//                                 }
//                                 else if(entry0_ == entry_)
//                                 {
//                                     lastEvent = ievent0+1;
//                                     break;
//                                 }
//                         }
                    }
//                     
//                     if(v_sample.at(iSample).sampleType() == Sample::data){
//                         histData_->Fill(recoBin_(branchVal_),eventWeight_*v_weight.at(iSample)); 
//                     }
//                     if(v_sample.at(iSample).sampleType() != Sample::data && isKinReco_ && !isTopGen_){
//                         histBgr_->Fill(recoBin_(branchVal_),eventWeight_*v_weight.at(iSample)); 
//                     }
//                     if(isKinReco_&&isTopGen_){
//                        histMigration_->Fill(genBin_(branchValGen_),recoBin_(branchVal_),eventWeight_*v_weight.at(iSample));
//                        histMigration_->Fill(genBin_(branchValGen_),0.,trueLevelWeight_*v_weight.at(iSample)-eventWeight_*v_weight.at(iSample));
//                        histSR_->Fill(recoBin_(branchVal_),eventWeight_*v_weight.at(iSample));
//                        
//                        //double weightCT = 1.;
//                        //Closure test
//                        //weightCT = (rewTopPtUp(branchValsGen_.at(1)));//pt up
//                        //weightCT = (rewTopPtDown(branchValsGen_.at(1)));//pt down
//                        //weightCT = (rewTopRapidityUp(branchValsGen_.at(0)));//y up
//                        //weightCT = (rewTopRapidityDown(branchValsGen_.at(0)));//y down
//                        
//                        //histMigration_->Fill(genBin_(branchValGen_),recoBin_(branchVals_),1*v_weight.at(iSample));
//                        //histSR_->Fill(recoBin_(branchVals_),weightCT*v_weight.at(iSample));
//                        
//                             if((branchValGen_ <coarseBins_.at(0) && uTrue_==0)  || (branchValGen_>coarseBins_.back() && oTrue_==0) )
//                             {
//                                 histBgrUo_->Fill(recoBin_(branchVal_),eventWeight_*v_weight.at(iSample));
//                                 //histBgrUo_->Fill(recoBin_(branchVal_),weightCT*v_weight.at(iSample));
//                             }
//                         
//                        Int_t recBin_temp = -999;
//                        Int_t genBin_temp = -999;
//                        
//                        recBin_temp = generatorDistribution_->GetGlobalBinNumber(branchVal_);
//                        genBin_temp = generatorDistribution_->GetGlobalBinNumber(branchValGen_);
//                        histPurity_->Fill(branchVal_,eventWeight_*v_weight.at(iSample));
//                        histStability_->Fill(branchValGen_,eventWeight_*v_weight.at(iSample));
// 
//                        histEffAllBins_->Fill(recBin_temp,eventWeight_*v_weight.at(iSample));
//                        histPurityAllBins_->Fill(recBin_temp,eventWeight_*v_weight.at(iSample));
//                        histStabilityAllBins_->Fill(genBin_temp,eventWeight_*v_weight.at(iSample));
//                        
//                        if(recBin_temp==genBin_temp){
//                            histRecoGen_->Fill(branchVal_,eventWeight_*v_weight.at(iSample));
//                            histRecoGenAllBins_->Fill(recBin_temp,eventWeight_*v_weight.at(iSample));
//                        }
//                        
//                     }
//                     if(v_sample.at(iSample).sampleType() != Sample::data && isKinReco_ ){
//                         //histMCReco_->Fill(recoBin_(branchVals_),eventWeight_*v_weight.at(iSample)); 
//                     }
                    
                    
                } // dataTree8 - loop
                
            
//                 dataTree0->Delete();
                dataTree8->Delete();
                dataFile->Delete();
            
                
                
            }//samples loop
            
            
            sampleInfo.save(sampleInfoFolder);
            
            //Control plots//
               writePlotCP(v_sample ,v_SampleHist_);
               
            // ... //
            
            
            //Unfolding//
            //runUnfolding(histMigration_,histData_,histBgr_,histBgrUo_,histUf_);
            //std::cout << "N entries in histBgrUo_ and histData_: " << histBgrUo_->Integral() << " " << histData_->Integral() << std::endl;
            //runUnfolding(histMigration_,histSR_,0,histBgrUo_,histUf_);
            //writePlotCT(histSG_,histRWSG_,histUf_);
            // ... //
            
            //Cross sections//
            //writePlotXSec(unfoldedData_,histGen_);
            //writePlotEPS();
            //writePlotEPSAllBins();
            
        }//channel loop
        
    }//systematics loop
    
    //std::cout<<"\n[PlotterBTop]:=== Finishing of plot production\n\n";
}


int PlotterBTop::genBin_(const float& val)
{
    int bin = -999;
    bin = generatorDistribution_->GetGlobalBinNumber(val);
    return bin;
}



int PlotterBTop::recoBin_(const float& val)
{
    int bin = -999;
    bin = detectorDistribution_->GetGlobalBinNumber(val);
    return bin;
}


//Closure test//
// void PlotterBTop::writePlotCT(TH1* histSG,TH1* histRW,TH1* histUf)
// {
//     // Prepare canvas and legend
//     TCanvas* canvas = utils::setCanvas();
//     TLegend* legend = this->setLegend();
//     canvas->cd();
//     
//     histSG->SetLineColor(2);
//     histSG->SetLineWidth(2);
//     histSG->SetStats(0);
//     legend->AddEntry(histSG,"gen","l");
//     histSG->Draw("hist");
//     
//     histRW->SetLineColor(kMagenta);
//     histRW->SetLineWidth(2);
//     histRW->SetLineStyle(2);
//     legend->AddEntry(histRW,"reweighted gen","l");
//     histRW->Draw("hist same");
//     
//     histUf->SetMarkerStyle(20);
//     histUf->Sumw2();
//     legend->AddEntry(histUf,"unfolded","pe");
//     histUf->Draw("e1 same");
//     
//     //Lines and bins//
//     double yMax = ( histUf->GetMaximum() > histRW->GetMaximum() ? histUf->GetMaximum()*1.2 : histRW->GetMaximum()*1.2 );
//     int Nbins1 = (int)(v_coarseBins_.at(1).size()-1)+(int)(v_oTrue_.at(1)+v_uTrue_.at(1));
//     int Nbins0 = (int)(v_coarseBins_.at(0).size()-1)+(int)(v_oTrue_.at(0)+v_uTrue_.at(0));
//     TH1D hist("hist","",Nbins1,0.5,Nbins0*Nbins1+0.5);
//     hist.SetTitle(";"+v_plotTitle_.at(1)+" "+v_plotUnits_.at(1));
//     hist.SetStats(0);
//     hist.GetXaxis()->SetTickLength(0);
//     hist.SetAxisRange(0,yMax ,"Y");
//         for(int i=1;i<=Nbins1;++i)
//         {
//             TString binName = "";
//             int isUnderflow = v_uTrue_.at(1);
//             int isOverflow =  v_oTrue_.at(1);
//             if(isUnderflow&&i==1) binName = "underflow";
//             else if(isOverflow&&i==Nbins1){
//                  binName = "overflow";
//             }
//             else{
//                 binName = utils::makeBinTitle("",v_coarseBins_.at(1).at(i-isUnderflow-1),v_coarseBins_.at(1).at(i-isUnderflow));
//             }
//             hist.GetXaxis()->SetBinLabel(i,binName.Data());
//         }
//         hist.SetTitle( utils::makeTitleBins(v_plotTitle_.at(0) +" "+ v_plotUnits_.at(0),v_coarseBins_.at(0),v_uTrue_.at(0),v_oTrue_.at(0) ) );
//     hist.Draw();
//  
//     histSG->Draw("hist same");
//     histRW->Draw("hist same");
//     histUf->Draw("e1 same");
//     for(int i=1;i<Nbins1;++i)
//         {
//             double xLine = i*(Nbins0) + 0.5;
//             TLine ln(xLine,-0.03*yMax,xLine,yMax);//FIXME: y1 to be related to pad size... 
//             ln.DrawClone("same");
//         }
//     legend->Draw("same");
//     utils::drawRatio(histUf, histRW, NULL,0.5,2,(TH1D*)hist.Clone());
// 
//     writeCanvas(canvas,"CT");
// }


// double PlotterBTop::rewTopPtUp(double pt)
// {
//     //func y = ax+b
//     double e = 0.2;
//     double b = 1 + e;
//     double a = (1. - b)/100.;
//     return a*pt + b;
// }
// 
// 
// double PlotterBTop::rewTopPtDown(double pt)
// {
//     //func y = ax+b
//     double e = 0.2;
//     double b = 1 - e;
//     double a = (1. - b)/100.;
//     return a*pt + b;
// }
// 
// 
// double PlotterBTop::rewTopRapidityUp(double y)
// {
//     //func y = axx+b
//     double b = 1 + 0.2;
//     double a = (1 - 0.1 - b)/1.2/1.2;
//     return a*y*y+b;
// }
// 
// 
// double PlotterBTop::rewTopRapidityDown(double y)
// {
//     //func y = axx+b
//     double b =1 - 0.2;
//     double a = (1 + 0.1 - b)/1.2/1.2;
//     return a*y*y+b;
// }


void PlotterBTop::runUnfolding(const TH2* histMigration,const TH1* histInput,
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
    TLegend* legend = this->setLegend();
    
    
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
    
    unfoldedData_ = unfold.GetOutput("unfoldedData","unfolding result",distributionName,projectionMode,useAxisBinning);
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



void PlotterBTop::writeCanvas( TCanvas* canvas, const TString name )
{
    canvas->Print(plotsFolder_+ name + "_" + plotName_ + ".pdf");
    canvas->Clear();
}


void PlotterBTop::writePlotCP(const std::vector<Sample>& v_sample ,const std::vector<TH1D* >& v_SampleHist)
{
    TCanvas* canvas = utils::setCanvas();
    TLegend* legend = this->setLegend();
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
    TString plotTitle = specname_ + ";" + plotTitle_ + " " + plotUnits_ + ";";
    
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
    canvas->Print(plotsFolder_+ "CP_" + plotName_ + ".pdf");
    canvas->Print(plotsFolder_+ "CP_" + plotName_ + ".png");
    canvas->Print(plotsFolder_+ "CP_" + plotName_ + ".root");
    
    //Draw recoge
    delete canvas;
    canvas = utils::setCanvas();
    hRecoGen_->Draw("colz");
    canvas->Print(plotsFolder_ + "RecoGen_" + plotName_ + ".pdf");
    canvas->Print(plotsFolder_ + "RecoGen_" + plotName_ + ".png");
    canvas->Print(plotsFolder_ + "RecoGen_" + plotName_ + ".root");
    
    
    
    
    //Delete
    stack->Delete();
    histData->Delete();
    delete canvas;
    delete legend;
    stacksum->Delete();
    
}


void PlotterBTop::writePlotEPSAllBins()
{
    // Prepare canvas and legend
    TCanvas* canvas = utils::setCanvas();
    TLegend* legend = this->setLegend();
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
    
    histPurityAllBins_->Draw("e same");
    histStabilityAllBins_->Draw("e same");
    histEffAllBins_->Draw("e same");
    legend->Draw("same");
    canvas->Print(plotsFolder_+ "EPS_AllBins_" + plotName_ + ".pdf");
    
    //Delete objects
    delete canvas;
    delete legend;
    
}



void PlotterBTop::writePlotEPS()
{
    // Prepare canvas and legend
    TCanvas* canvas = utils::setCanvas();
    TLegend* legend = this->setLegend();
    canvas->cd();
    
    histPurity_->Divide(histRecoGen_,histPurity_,1,1,"B");
    histPurity_->SetMarkerStyle(23);
    histPurity_->SetMarkerColor(4);
    
    histStability_->Divide(histRecoGen_,histStability_,1,1,"B");
    histStability_->SetMarkerStyle(22);
    histStability_->SetMarkerColor(2);
    
    legend->AddEntry(histPurity_,"Purity","pe");
    legend->AddEntry(histStability_,"Stability","pe");
    
    histPurity_->SetAxisRange(0,histPurity_->GetMaximum() > histStability_->GetMaximum() ? histPurity_->GetMaximum()*1.15 : histStability_->GetMaximum()*1.15 ,"Y");
    histPurity_->Draw("e");
    histStability_->Draw("e same");
    //legend->Draw("same");
    canvas->Print(plotsFolder_+ "EPS_" + plotName_ + ".pdf");
    
    //Delete objects
    delete canvas;
    delete legend;
    
}



void PlotterBTop::writePlotXSec(const TH1* hData,const TH1* hMC)
{
    // Prepare canvas and legend
    TCanvas* canvas = utils::setCanvas();
    TLegend* legend = this->setLegend();
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
        
        histData->Scale(1.,"width");
        histMC->Scale(1.,"width");
        
        histData->SetAxisRange(0,( histData->GetMaximum() > histMC->GetMaximum() ? histData->GetMaximum()*1.15 : histMC->GetMaximum()*1.15 ),"Y");
        histData->Draw("e");
        histMC->Draw("hist same");
        histData->Draw("e,same");
        utils::drawRatio(histData, histMC, NULL,0.5,2);
        canvas->Print(plotsFolder_+ "xSec_" + plotName_ + ".pdf");
        
        histMC->Delete();
        histData->Delete();

    //Delete objects
    delete canvas;
    delete legend;
}



TLegend* PlotterBTop::setLegend()
{
    TLegend* legend = new TLegend(0.73,0.60,0.90,0.85);
    legend->SetFillStyle(0);
    legend->SetBorderSize(0);
    legend->SetX1NDC(1.0 - gStyle->GetPadRightMargin() - gStyle->GetTickLength() - 0.25);
    legend->SetY1NDC(1.0 - gStyle->GetPadTopMargin()  - gStyle->GetTickLength() - 0.05 - legend->GetBorderSize()*0.04);
    legend->SetX2NDC(1.0 - gStyle->GetPadRightMargin() - gStyle->GetTickLength());
    legend->SetY2NDC(1.0 - gStyle->GetPadTopMargin()  - gStyle->GetTickLength());
    legend->Clear();
    return legend;
}



const SystematicChannelFactors& PlotterBTop::scaleFactors(const TString& name)
{
    const TString fullStepname = ttbar::extractSelectionStepAndJetCategory(name);
    
    // Check if map contains already scale factors for this step, else access and fill them
    if(m_stepFactors_.find(fullStepname) == m_stepFactors_.end()){
        m_stepFactors_[fullStepname] = samples_.globalWeights(fullStepname).first;
    }
    
    return m_stepFactors_.at(fullStepname);
}

