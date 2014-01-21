#include <fstream>
#include <iostream>
#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <algorithm>
#include <sstream>

#include <TSystem.h>
#include <TFile.h>
#include <TH1.h>
#include <TH1D.h>
#include <TH2.h>
#include <TSelectorList.h>
#include <TString.h>
#include <TObjArray.h>
#include <TRandom3.h>
#include <Rtypes.h>
#include <TMath.h>

#include "ScaleFactors.h"
#include "classes.h"
#include "utils.h"
#include "TopAnalysis/ZTopUtils/ext/interface/JetCorrectorParameters.h"
#include "TopAnalysis/ZTopUtils/ext/interface/JetCorrectionUncertainty.h"









std::string common::assignFolder(const char* baseDir, const TString& channel, const TString& systematic)
{
    std::string path("");
    
    // Create all subdirectories contained in baseDir
    TObjArray* a_subDir = TString(baseDir).Tokenize("/");
    for(Int_t iSubDir = 0; iSubDir < a_subDir->GetEntriesFast(); ++iSubDir){
        const TString& subDir = a_subDir->At(iSubDir)->GetName();
        path.append(subDir);
        path.append("/");
        gSystem->MakeDirectory(path.c_str());
    }
    
    // Create subdirectories for systematic and channel
    path.append(systematic);
    path.append("/");
    gSystem->MakeDirectory(path.c_str());
    path.append(channel);
    path.append("/");
    gSystem->MakeDirectory(path.c_str());
    
    return path;
}



std::string common::accessFolder(const char* baseDir, const TString& channel,
                         const TString& systematic, const bool allowNonexisting)
{
    // Build directory path
    std::string path(baseDir);
    path.append("/");
    path.append(systematic);
    path.append("/");
    path.append(channel);
    path.append("/");
    
    // Check if directory really exists
    if(!gSystem->OpenDirectory(path.c_str())){
        if(allowNonexisting){
            // It is allowed to request a folder which does not exist, so return empty string silently
            return "";
        }
        else{
            std::cerr<<"ERROR! Request to access directory is not possible, because it does not exist. Directory name: "<<path
                     <<"\n...break\n"<<std::endl;
            exit(237);
        }
    }
    
    return path;
}







double ScaleFactorHelpers::get2DSF(TH2* histo, const double x, const double y)
{
    int xbin, ybin, dummy;
    histo->GetBinXYZ(histo->FindBin(x, y), xbin, ybin, dummy);
    //overflow to last bin
    xbin = std::min(xbin, histo->GetNbinsX());
    ybin = std::min(ybin, histo->GetNbinsY());
    return histo->GetBinContent(xbin, ybin);
}



double ScaleFactorHelpers::median(TH1* h1)
{ 
   int nBin = h1->GetXaxis()->GetNbins();
   std::vector<double> x(nBin);
   h1->GetXaxis()->GetCenter(&x[0]);
   TH1D* h1D = dynamic_cast<TH1D*>(h1);
   if(!h1D){
       std::cerr << "Median needs a TH1D!\n";
       exit(7);
   }
   const double* y = h1D->GetArray(); 
   // exclude underflow/overflows from bin content array y
   return TMath::Median(nBin, &x[0], &y[1]); 
}













// --------------------------- Methods for LeptonScaleFactors ---------------------------------------------



LeptonScaleFactors::LeptonScaleFactors(const char* electronSFInputFileName,
                                       const char* muonSFInputFileName,
                                       const LeptonScaleFactors::Systematic& systematic):
h2_ElectronIDSFpteta(0),
h2_MuonIDSFpteta(0)
{
    std::cout<<"--- Beginning preparation of lepton scale factors\n";
    
    for(const auto& lepton : {electron, muon}){
        std::string inputFileName(common::DATA_PATH_COMMON());
        inputFileName.append("/");
        std::string histogramName;
        if(lepton == electron){
            inputFileName.append(electronSFInputFileName);
            histogramName = "ElectronIdIsoSF";
            h2_ElectronIDSFpteta = this->prepareLeptonIDSF(inputFileName, histogramName, systematic);
        }
        else if(lepton == muon){
            inputFileName.append(muonSFInputFileName);
            histogramName = "MuonIdIsoSF";
            h2_MuonIDSFpteta = this->prepareLeptonIDSF(inputFileName, histogramName, systematic);
        }
    }
    
    const bool electronInputFound(h2_ElectronIDSFpteta);
    const bool muonInputFound(h2_MuonIDSFpteta);
    if(electronInputFound || muonInputFound){
        std::cout<<"Found lepton Id/Iso scale factors for: ";
        if(electronInputFound) std::cout<<"electron , ";
        if(muonInputFound) std::cout<<"muon , ";
        std::cout<<std::endl;
    }
    if(!electronInputFound || !muonInputFound){
        std::cout<<"Could NOT find lepton Id/Iso scale factors for: ";
        if(!electronInputFound) std::cout<<"electron , ";
        if(!muonInputFound) std::cout<<"muon , ";
        std::cout<<std::endl;
    }
    std::cout<<"=== Finishing preparation of lepton scale factors\n\n";
}



TH2* LeptonScaleFactors::prepareLeptonIDSF(const std::string& inputFileName,
                                           const std::string& histogramName,
                                           const LeptonScaleFactors::Systematic& systematic)const
{
    // Access file containing scale factors
    TFile scaleFactorFile(inputFileName.c_str());
    if (scaleFactorFile.IsZombie()){
        std::cout<<"File containing lepton Id/Iso scale factors not found: "<<inputFileName
                 <<"\nAssuming ScaleFactor = 1.\n";
        return 0;
    }
    
    // Access histogram containing scale factors
    TH2* h_scaleFactorPtEta(0);
    h_scaleFactorPtEta = dynamic_cast<TH2*>(scaleFactorFile.Get(histogramName.c_str()));
    if(!h_scaleFactorPtEta){
        std::cout<<"TH2 for lepton Id/Iso scale factors not found: "<<histogramName
                 <<"\nAssuming ScaleFactor = 1.\n";
        return 0;
    }
    
    // Apply systematic variations
    if(systematic != nominal){
        double factor(0.);
        if(systematic == vary_up) factor = 1.;
        else if(systematic == vary_down) factor = -1.;
        else{
            std::cerr<<"ERROR in LeptonScaleFactors! Systematic with undefined behaviour requested\n...break\n"<<std::endl;
            exit(21);
        }
        
        for (int i = 1; i <= h_scaleFactorPtEta->GetNbinsX(); ++i) {
            for (int j = 1; j <= h_scaleFactorPtEta->GetNbinsY(); ++j) {
                h_scaleFactorPtEta->SetBinContent(i, j,
                    h_scaleFactorPtEta->GetBinContent(i,j) + factor*h_scaleFactorPtEta->GetBinError(i,j));
            }
        }
    }
    
    // Store histogram in memory and close file
    h_scaleFactorPtEta->SetDirectory(0);
    scaleFactorFile.Close();
    
    return h_scaleFactorPtEta;
}



double LeptonScaleFactors::getLeptonIDSF(const int leadingLeptonIndex, const int nLeadingLeptonIndex,
                                         const VLV& leptons, const std::vector<int>& lepPdgIds)const
{
    if(!h2_ElectronIDSFpteta || !h2_MuonIDSFpteta) return 1.;
    
    const LV& leadingLepton(leptons.at(leadingLeptonIndex));
    const LV& nLeadingLepton(leptons.at(nLeadingLeptonIndex));
    const int leadingPdgId(std::abs(lepPdgIds.at(leadingLeptonIndex)));
    const int nLeadingPdgId(std::abs(lepPdgIds.at(nLeadingLeptonIndex)));
    
    if(leadingPdgId==11 && nLeadingPdgId==11)
        return ScaleFactorHelpers::get2DSF(h2_ElectronIDSFpteta, leadingLepton.Eta(), leadingLepton.pt()) *
               ScaleFactorHelpers::get2DSF(h2_ElectronIDSFpteta, nLeadingLepton.Eta(), nLeadingLepton.pt());
    if(leadingPdgId==13 && nLeadingPdgId==13)
        return ScaleFactorHelpers::get2DSF(h2_MuonIDSFpteta, leadingLepton.Eta(), leadingLepton.pt()) *
               ScaleFactorHelpers::get2DSF(h2_MuonIDSFpteta, nLeadingLepton.Eta(), nLeadingLepton.pt());
    if(leadingPdgId==13 && nLeadingPdgId==11)
        return ScaleFactorHelpers::get2DSF(h2_MuonIDSFpteta, leadingLepton.Eta(), leadingLepton.pt()) *
               ScaleFactorHelpers::get2DSF(h2_ElectronIDSFpteta, nLeadingLepton.Eta(), nLeadingLepton.pt());
    if(leadingPdgId==11 && nLeadingPdgId==13)
        return ScaleFactorHelpers::get2DSF(h2_ElectronIDSFpteta, leadingLepton.Eta(), leadingLepton.pt()) *
               ScaleFactorHelpers::get2DSF(h2_MuonIDSFpteta, nLeadingLepton.Eta(), nLeadingLepton.pt());
    std::cout<<"WARNING in method getLeptonIDSF! LeptonPdgIds are not as expected (pdgId1, pdgId2): "
             <<leadingPdgId<<" , "<<nLeadingPdgId<<"\n...will return scale factor = 1.\n";
    return 1.;
}



double LeptonScaleFactors::scaleFactorAllLeptons(const std::vector<int>& allLeptonIndices,
                                                 const VLV& leptons, const std::vector<int>& lepPdgIds)const
{
    if(!h2_ElectronIDSFpteta || !h2_MuonIDSFpteta) return 1.;
    
    double result(1.);
    
    for(const int index : allLeptonIndices){
        const int absPdgId(std::abs(lepPdgIds.at(index)));
        TH2* histo(0);
        if(absPdgId==11) histo = h2_ElectronIDSFpteta;
        else if(absPdgId==13) histo = h2_MuonIDSFpteta;
        else{
            std::cout<<"WARNING in method scaleFactorAllLeptons! LeptonPdgId is not as expected: "
                     <<absPdgId<<"\n...will use scale factor = 1.\n";
             continue;
        }
        
        result *= ScaleFactorHelpers::get2DSF(histo, leptons.at(index).eta(), leptons.at(index).pt());
    }
    
    return result;
}













// --------------------------- Methods for TriggerScaleFactors ---------------------------------------------



TriggerScaleFactors::TriggerScaleFactors(const char* inputFileSuffix,
                                         const std::vector<std::string>& channels,
                                         const Systematic& systematic):
h2_eeTrigSFeta(0),
h2_emuTrigSFeta(0),
h2_mumuTrigSFeta(0)
{
    std::cout<<"--- Beginning preparation of trigger scale factors\n";
    
    bool someChannelFound(false);
    bool someChannelNotFound(false);
    std::stringstream channelsFound;
    std::stringstream channelsNotFound;
    
    std::vector<TString> triggerInputFileNames;
    for(const auto& channel : channels){
        TString fullName(common::DATA_PATH_COMMON());
        fullName.Append("/triggerSummary_");
        fullName.Append(channel);
        fullName.Append(inputFileSuffix);
        triggerInputFileNames.push_back(fullName);
        TH2* h2_triggerSF = this->prepareTriggerSF(fullName, systematic);
        if(channel == "ee") h2_eeTrigSFeta = h2_triggerSF;
        else if(channel == "emu") h2_emuTrigSFeta = h2_triggerSF;
        else if(channel == "mumu") h2_mumuTrigSFeta = h2_triggerSF;
        else{
            std::cerr<<"ERROR in TriggerScaleFactors! Invalid channel requested: "<<channel<<"\n...break\n"<<std::endl;
            exit(23);
        }
        if(h2_triggerSF){
            someChannelFound = true;
            channelsFound<<channel<<" , ";
        }
        else{
            someChannelNotFound = true;
            channelsNotFound<<channel<<" , ";
        }
    }
    
    if(someChannelFound) std::cout<<"Found trigger scale factors for requested channels: "<<channelsFound.str()<<"\n";
    if(someChannelNotFound)  std::cout<<"Could NOT find trigger scale factors for requested channels: "<<channelsNotFound.str()<<"\n";
    std::cout<<"=== Finishing preparation of trigger scale factors\n\n";
}



TH2* TriggerScaleFactors::prepareTriggerSF(const TString& fileName, const Systematic& systematic)const
{
    TFile trigEfficiencies(fileName);
    if(trigEfficiencies.IsZombie()){
        std::cout << "Trigger efficiencies not found. Assuming ScaleFactor = 1.\n";
        std::cout << "Currently triggerEfficieny files can be found in diLeptonic/data folder\n\n";
        return 0;
    }
    
    //Right now pT efficiency flat ==> Not used
    TH2* h_TrigSFeta(0);
    h_TrigSFeta = dynamic_cast<TH2*>(trigEfficiencies.Get("scalefactor eta2d with syst"));
    if(!h_TrigSFeta){
        std::cout<<"TH2 >>TH scalefactor eta<< is not in the file "<<trigEfficiencies.GetName()<<"\n";
        return 0;
    }
    
    if(systematic != nominal){
        double factor(0.);
        if(systematic == vary_up) factor = 1.;
        else if(systematic == vary_down) factor = -1.;
        else{
            std::cerr<<"ERROR in TriggerScaleFactors! Systematic with undefined behaviour requested\n...break\n"<<std::endl;
            exit(24);
        }
        
        for (int i = 1; i <= h_TrigSFeta->GetNbinsX(); ++i) {
            for (int j = 1; j <= h_TrigSFeta->GetNbinsY(); ++j) {
                h_TrigSFeta->SetBinContent(i, j,
                    h_TrigSFeta->GetBinContent(i,j) + factor*h_TrigSFeta->GetBinError(i,j));
            }
        }
    }
    
    h_TrigSFeta->SetDirectory(0);
    trigEfficiencies.Close();
    
    return h_TrigSFeta;
}



double TriggerScaleFactors::getTriggerSF(const int leptonXIndex, const int leptonYIndex,
                                         const VLV& leptons, const TString& channel)const
{
    TH2* h_TrigSFeta(0);
    if(channel == "ee") h_TrigSFeta = h2_eeTrigSFeta;
    else if(channel == "emu") h_TrigSFeta = h2_emuTrigSFeta;
    else if(channel == "mumu") h_TrigSFeta = h2_mumuTrigSFeta;
    else{
        std::cerr<<"ERROR in TriggerScaleFactors! Invalid channel requested: "<<channel<<"\n...break\n"<<std::endl;
        exit(25);
    }
    if (!h_TrigSFeta) return 1.;
    
    //For 'ee' and 'mumu' channels Xaxis of the 2D plots is the highest pT lepton
    // for the 'emu' channel Xaxis is the electron and Y axis muon
    const LV& leptonX(leptons.at(leptonXIndex));
    const LV& leptonY(leptons.at(leptonYIndex));
    return ScaleFactorHelpers::get2DSF(h_TrigSFeta, std::fabs(leptonX.eta()), std::fabs(leptonY.eta()));
}















// --------------------------- Methods for BtagScaleFactors ---------------------------------------------




BtagScaleFactors::BtagScaleFactors(const char* btagEfficiencyInputDir,
                             const char* btagEfficiencyOutputDir,
                             const std::vector<std::string>& channels,
                             TString systematic):
inputDirName_(btagEfficiencyInputDir),
outputDirName_(btagEfficiencyOutputDir),
fileName_("ttbarsignalplustau.root")
{
    std::cout<<"--- Beginning preparation of b-tagging scale factors\n";
    if (systematic == "" || systematic.Contains("PDF") || systematic.Contains("closure")) systematic = "Nominal";
    // Check if all relevant input files are available
    bool allInputFilesAvailable(true);
    for(const auto& channel : channels){
        std::string bTagInputFile = common::accessFolder(inputDirName_.c_str(),channel, systematic, true).append(channel).append("_").append(fileName_);
        ifstream inputFileStream(bTagInputFile);
        // Setting the file and sample name for each channel in the map if the file exists
        if(inputFileStream.is_open() && bTagInputFile.length() > fileName_.length() + channel.length() + 1) {
            channelFileNames_[channel] = bTagInputFile;
            
            std::string sampleName(bTagInputFile);
            sampleName.erase(0,inputDirName_.length());
            while(sampleName.at(0) == '/') sampleName.erase(0,1);
            channelSampleNames_[channel] = sampleName;
        } else {
            std::cout<< "******************************************************\n"
                     << "Btag efficiency file [" << bTagInputFile << "] doesn't exist.\n"
                     << "RUNNING WITHOUT BTAGSF!!!\n"
                     << "To create the file, run (for each systematic 'SYST'):\n"
                     << "\t> ./build/load_Analysis -f ttbarsignalplustau.root -c emu -s SYST\n"
                     << "\t> ./build/load_Analysis -f ttbarsignalplustau.root -c ee -s SYST\n"
                     << "\t> ./build/load_Analysis -f ttbarsignalplustau.root -c mumu -s SYST\n"
                     << " and move the selectionRoot/BTagEff directory to the cwd:\n"
                     << "\t> mv selectionRoot/BTagEff .\n"
                     << "This error is NOT fatal, not applying any b-tag rescaling\n"
                     << "*******************************************************\n";
            allInputFilesAvailable = false;
            break;
        }   // If file couldn't be opened for reading
    }   // End of loop over channels
    
    if(!allInputFilesAvailable){
        std::cout<<"Not all input files for b-tagging efficiencies available\n"
                 <<"\t-->  Efficiencies will not be used, but produced in Analysis\n";
        setMakeEff(true);
        // Resetting the root file names for storing btagging efficiencies for each channel in the map
        for(const auto& channel : channels){
            std::string bTagOutputFile = common::assignFolder(outputDirName_.c_str(), channel, systematic).append(channel).append("_").append(fileName_);
            channelFileNames_[channel] = bTagOutputFile;
            
            std::string sampleName(bTagOutputFile);
            sampleName.erase(0,outputDirName_.length());
            while(sampleName.at(0) == '/') sampleName.erase(0,1);
            channelSampleNames_[channel] = sampleName;
        }
    } else setMakeEff(false);

    // Set systematic if it is an allowed one for btag efficiencies, else set to nominal
    if(systematic == "BTAG_UP") setSystematic(systematics::heavyup);
    else if(systematic == "BTAG_DOWN") setSystematic(systematics::heavydown);
    else if(systematic == "BTAG_PT_UP") setSystematic(systematics::heavyuppt);
    else if(systematic == "BTAG_PT_DOWN") setSystematic(systematics::heavydownpt);
    else if(systematic == "BTAG_ETA_UP") setSystematic(systematics::heavyupeta);
    else if(systematic == "BTAG_ETA_DOWN") setSystematic(systematics::heavydowneta);
    else if(systematic == "BTAG_LJET_UP") setSystematic(systematics::lightup);
    else if(systematic == "BTAG_LJET_DOWN") setSystematic(systematics::lightdown);
    else if(systematic == "BTAG_LJET_PT_UP") setSystematic(systematics::lightuppt);
    else if(systematic == "BTAG_LJET_PT_DOWN") setSystematic(systematics::lightdownpt);
    else if(systematic == "BTAG_LJET_ETA_UP") setSystematic(systematics::lightupeta);
    else if(systematic == "BTAG_LJET_ETA_DOWN") setSystematic(systematics::lightdowneta);
    else setSystematic(systematics::nominal);
    
    std::cout<<"=== Finishing preparation of b-tagging scale factors\n\n";
}



bool BtagScaleFactors::makeEfficiencies()
{
    return getMakeEff();
}



void BtagScaleFactors::prepareBTags(TSelectorList* output, const std::string& channel)
{
    std::string inputFileName = channelFileNames_.at(channel);
    std::string sampleName = channelSampleNames_.at(channel);

    // Set pointer to output, so that histograms are owned by it
    selectorList_ = output;

    // Stopping if no efficiency histograms need to be read (production mode)
    if(getMakeEff()) {
        setSampleName(sampleName);
        return;
    }


    std::vector<TH2D> histos;
    std::vector<TH2D> effhistos;
    std::vector<float> medians;

    // Load per-jet efficiencies file
    TFile file(inputFileName.c_str(),"READ");

    // Disabling referencing histograms to gDirectory to prevent crashing when closing the root file
    TH1::AddDirectory(false);
    TH2D* tempHisto;
    // Reading the tag and efficiency histograms from the file and adding them to the maps in proper order
    for(int id = 0; ; ++id) {
        std::string histoName = histoNameAtId(id, histoTypes::tag);
        if(histoName != "") {
            file.GetObject(histoName.c_str(), tempHisto);
            if(!tempHisto) 
                throw std::runtime_error("BtagScaleFactors::prepareBTags Couldn't find all [tag] histograms in the input root file");
            histos.push_back(*tempHisto);
            tempHisto = 0;
        }
        
        std::string effHistoName = histoNameAtId(id, histoTypes::eff);
        if(effHistoName != "") {
            file.GetObject(effHistoName.c_str(), tempHisto);
            if(!tempHisto) 
                throw std::runtime_error("BtagScaleFactors::prepareBTags Couldn't find all [eff] histograms in the input root file");
            effhistos.push_back(*tempHisto);
        }
        
        if(histoName == "" && effHistoName == "") break;
    }
    // Reading and extracting madian values from the histogram
    TH1* medianHisto;
    file.GetObject("medians", medianHisto);
    if(!medianHisto) throw std::runtime_error("BtagScaleFactors::prepareBTags Couldn't find [medians] histogram in the input root file");
    for(int i = 0; i<medianHisto->GetNbinsX(); ++i) medians.push_back((float)medianHisto->GetBinContent(i+1));

    histos_[sampleName] = histos;
    effhistos_[sampleName] = effhistos;
    medianMap_[sampleName] = medians;

    file.Close();

    if(setSampleName(sampleName) < 0) throw std::runtime_error("BtagScaleFactors::prepareBTags Tried to set a non-existing sampleName");

}




void BtagScaleFactors::indexOfBtags(std::vector<int>& bjetIndices,
                                 const std::vector<int>& jetIndices,
                                 const VLV& jets,
                                 const std::vector<int>& jetPartonFlavours,
                                 const std::vector<double>& btagDiscriminants)const
{
    bjetIndices.clear();

    std::vector<int> tagged_indices;
    for(const int index : jetIndices){
        //Skip jets where there is no partonFlavour
        if(jetPartonFlavours.at(index) == 0) continue;
        LV jet = jets.at(index);

        // Preparing a seed for the random jet retagging
        const unsigned int seed = std::abs( static_cast<int>( 1.e6*sin( 1.e6*jet.Phi() ) ) );
        bool isTagged = jetIsTagged( jet.pt(), std::fabs(jet.eta()), jetPartonFlavours.at(index), btagDiscriminants.at(index), seed );
        if(isTagged) tagged_indices.push_back(index);
    }
    bjetIndices = tagged_indices;
}




double BtagScaleFactors::calculateBtagSF(const std::vector<int>& jetIndices,
                                      const VLV& jets,
                                      const std::vector<int>& jetPartonFlavours)
{
    resetCounter();
    for(const int index : jetIndices){
        countJet(jets.at(index).pt(), std::fabs(jets.at(index).eta()), jetPartonFlavours.at(index));
    }

    // per-event SF calculation
    double scale_factor = getEventSF();

    if(std::fabs(scale_factor-1.) > 0.05) scale_factor = 1.;

    return scale_factor;
}




void BtagScaleFactors::fillBtagHistograms(const std::vector<int>& jetIndices,
                                          const std::vector<double>& bTagDiscriminant,
                                          const VLV& jets,
                                          const std::vector<int>& jetPartonFlavours,
                                          const double weight)
{
    for(const int index : jetIndices){
        fillEff(jets.at(index).pt(), std::fabs(jets.at(index).eta()), std::abs(jetPartonFlavours.at(index)), bTagDiscriminant.at(index), weight);
    }
}



void BtagScaleFactors::produceBtagEfficiencies(const std::string& channel)
{
    std::string outputFileName = channelFileNames_.at(channel);
    std::string sampleName = channelSampleNames_.at(channel);

    TFile file(outputFileName.c_str(),"RECREATE");

    // Creating the histograms
    makeEffs();

    // Writing each histogram to file
    std::vector<TH2D> histos = histos_.at(sampleName);
    for(TH2D histo : histos) {
        histo.Write();
    }

    // Writing each efficiency hitogram to file
    std::vector<TH2D> effhistos = effhistos_.at(sampleName);
    for(TH2D histo : effhistos) {
        histo.Write();
    }
    
    // Writing medians to the histogram and then to the file
    size_t nMedians = (size_t)medians::length_median;
    TH1D histo("medians","Medians;Property id;Median value",nMedians, 0, nMedians);
    
    if(nMedians != medianMap_.at(sampleName).size()) {
        throw std::range_error("BtagScaleFactors::produceBtagEfficiencies Numbers of stored and designed median values differ");
    }
    
    for(size_t i = 0; i<nMedians; ++i) {
        histo.SetBinContent(i+1, medianMap_.at(sampleName).at(i));
    }
    histo.Write();

    
    file.Close();


    std::cout<<"Done with production of b-tag efficiency file: "<< outputFileName
             <<"\n\n"<<std::endl;
}








// --------------------------- Methods for JetEnergyResolutionScaleFactors ---------------------------------------------




JetEnergyResolutionScaleFactors::JetEnergyResolutionScaleFactors(const JetEnergyResolutionScaleFactors::Systematic& systematic)
{
    std::cout<<"--- Beginning preparation of JER scale factors\n";
    
    // Hardcoded eta ranges
    v_etaRange_ = {0.5, 1.1, 1.7, 2.3, 10};
    
    // Hardcoded scale factors for eta ranges, nominal is = {1.052, 1.057, 1.096, 1.134, 1.288};
    if(systematic == vary_up){
        v_etaScaleFactor_ = {1.115, 1.114, 1.161, 1.228, 1.488};
        std::cout<<"Apply systematic variation: up\n";
    }
    else if(systematic == vary_down){
        v_etaScaleFactor_ = {0.990, 1.001, 1.032, 1.042, 1.089};
        std::cout<<"Apply systematic variation: down\n";
    }
    else{
        std::cerr<<"Error in constructor of JetEnergyResolutionScaleFactors! Systematic variation not allowed\n...break\n"<<std::endl;
        exit(621);
    }
    
    // Check correct size
    if(v_etaRange_.size() != v_etaScaleFactor_.size()){
        std::cerr<<"Error in constructor of JetEnergyResolutionScaleFactors! "
                 <<"Different number of intervals in eta and corresponding scale factors (intervals, scale factors): "
                 <<v_etaRange_.size()<<" , "<<v_etaScaleFactor_.size()<<"\n...break\n"<<std::endl;
        exit(622);
    }
    
    std::cout<<"=== Finishing preparation of JER scale factors\n\n";
}



void JetEnergyResolutionScaleFactors::applySystematic(VLV* jets, VLV* jetsForMET, LV* met,
                                                      const std::vector<double>* jetJERSF, const std::vector<double>* jetForMETJERSF,
                                                      const VLV* associatedGenJet, const VLV* associatedGenJetForMET)const
{
    // This first loop will correct the jet collection that is used for jet selections
    for(size_t iJet = 0; iJet < jets->size(); ++iJet){
        size_t jetEtaBin = 0;

        for(size_t iBin = 0; iBin < v_etaRange_.size(); ++iBin){
            if(std::fabs(jets->at(iJet).eta()) < v_etaRange_.at(iBin)){
                jetEtaBin = iBin;
                break;
            }
        }

        if(jetJERSF->at(iJet) != 0.){
            jets->at(iJet) *= 1./jetJERSF->at(iJet);

            // FIXME: should this factor really be =0. in case no associatedGenJet is found ?
            double factor = 0.;

            if(associatedGenJet->at(iJet).pt() != 0.) factor = 1. + (v_etaScaleFactor_.at(jetEtaBin) - 1.)*(1. - (associatedGenJet->at(iJet).pt()/jets->at(iJet).pt()));
            if(jetJERSF->at(iJet) == 1.) factor = 1.;

            jets->at(iJet) *= factor;
        }
    }

    // This second loop will correct the jet collection that is used to modify the MET
    double JEC_dpX =0.;
    double JEC_dpY =0.;
    for(size_t iJet = 0; iJet < jetsForMET->size(); ++iJet){
        size_t jetEtaBin = 0;
        for(size_t iBin = 0; iBin < v_etaRange_.size(); ++iBin){
            if(std::fabs(jetsForMET->at(iJet).eta()) < v_etaRange_.at(iBin)){
                jetEtaBin = iBin;
                break;
            }
        }

        if(jetForMETJERSF->at(iJet) != 0.){
            const double dpX = jetsForMET->at(iJet).px();
            const double dpY = jetsForMET->at(iJet).py();
            jetsForMET->at(iJet) *= 1./jetForMETJERSF->at(iJet);

            // FIXME: should this factor really be =0. in case no associatedGenJet is found ?
            double factor = 0.;
            if(associatedGenJetForMET->at(iJet).pt() != 0.) factor = 1. + (v_etaScaleFactor_.at(jetEtaBin) - 1.)*(1. - (associatedGenJetForMET->at(iJet).pt()/jetsForMET->at(iJet).pt()));
            if(jetForMETJERSF->at(iJet) == 1.) factor = 1.;

            jetsForMET->at(iJet) *= factor;
            JEC_dpX += jetsForMET->at(iJet).px() - dpX;
            JEC_dpY += jetsForMET->at(iJet).py() - dpY;
        }
    }

    // Adjust the MET
    const double scaledMETPx = met->px() - JEC_dpX;
    const double scaledMETPy = met->py() - JEC_dpY;
    met->SetPt(std::sqrt(scaledMETPx*scaledMETPx + scaledMETPy*scaledMETPy));
}








// --------------------------- Methods for JetEnergyScaleScaleFactors ---------------------------------------------




JetEnergyScaleScaleFactors::JetEnergyScaleScaleFactors(const char* jesUncertaintySourceFile,
                                                       const JetEnergyScaleScaleFactors::Systematic& systematic):
jetCorrectionUncertainty_(0),
varyUp_(false)
{
    std::cout<<"--- Beginning preparation of JES scale factors\n";
    
    std::string inputFileName(common::DATA_PATH_COMMON());
    inputFileName.append("/");
    inputFileName.append(jesUncertaintySourceFile);
    jetCorrectionUncertainty_ = new ztop::JetCorrectionUncertainty(ztop::JetCorrectorParameters(inputFileName.data(), "Total"));
    
    if(systematic == vary_up){
        varyUp_ = true;
        std::cout<<"Apply systematic variation: up\n";
    }
    else if(systematic == vary_down){
        varyUp_ = false;
        std::cout<<"Apply systematic variation: down\n";
    }
    else{
        std::cerr<<"Error in constructor of JetEnergyScaleScaleFactors! Systematic variation not allowed\n...break\n"<<std::endl;
        exit(624);
    }
    
    std::cout<<"=== Finishing preparation of JES scale factors\n\n";
}



JetEnergyScaleScaleFactors::~JetEnergyScaleScaleFactors()
{
    delete jetCorrectionUncertainty_;
}



void JetEnergyScaleScaleFactors::applySystematic(VLV* jets, VLV* jetsForMET, LV* met)const
{
    // This first loop will correct the jet collection that is used for jet selections
    for(size_t iJet = 0; iJet < jets->size(); ++iJet){
        jetCorrectionUncertainty_->setJetPt(jets->at(iJet).pt());
        jetCorrectionUncertainty_->setJetEta(jets->at(iJet).eta());
        const double dunc = jetCorrectionUncertainty_->getUncertainty(true);

        if(varyUp_) jets->at(iJet) *= 1. + dunc;
        else jets->at(iJet) *= 1. - dunc;
    }
    
    // This second loop will correct the jet collection that is used for modifying MET
    double JEC_dpX =0.;
    double JEC_dpY =0.;
    for(size_t iJet = 0; iJet < jetsForMET->size(); ++iJet){
        
        const double dpX = jetsForMET->at(iJet).px();
        const double dpY = jetsForMET->at(iJet).py();

        jetCorrectionUncertainty_->setJetPt(jetsForMET->at(iJet).pt());
        jetCorrectionUncertainty_->setJetEta(jetsForMET->at(iJet).eta());
        const double dunc = jetCorrectionUncertainty_->getUncertainty(true);

        if(varyUp_) jetsForMET->at(iJet) *= 1. + dunc;
        else jetsForMET->at(iJet) *= 1. - dunc;

        JEC_dpX += jetsForMET->at(iJet).px() - dpX;
        JEC_dpY += jetsForMET->at(iJet).py() - dpY;
    }

    // Adjust the MET
    const double scaledMETPx = met->px() - JEC_dpX;
    const double scaledMETPy = met->py() - JEC_dpY;
    met->SetPt(std::sqrt(scaledMETPx*scaledMETPx + scaledMETPy*scaledMETPy));
}






