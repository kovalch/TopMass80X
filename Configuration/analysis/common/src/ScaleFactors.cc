#include <fstream>
#include <iostream>
#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <algorithm>
#include <sstream>

#include <TFile.h>
#include <TH1.h>
#include <TH1D.h>
#include <TH2.h>
#include <TSelectorList.h>
#include <TString.h>
#include <TRandom3.h>
#include <Rtypes.h>
#include <TMath.h>

#include "ScaleFactors.h"
#include "sampleHelpers.h"
#include "classes.h"
#include "utils.h"
#include "TopAnalysis/ZTopUtils/interface/PUReweighter.h"
#include "TopAnalysis/ZTopUtils/ext/interface/JetCorrectorParameters.h"
#include "TopAnalysis/ZTopUtils/ext/interface/JetCorrectionUncertainty.h"









double ScaleFactorHelpers::get2DSF(const TH2* const histo, const double x, const double y)
{
    int xbin, ybin, dummy;
    histo->GetBinXYZ(histo->FindFixBin(x, y), xbin, ybin, dummy);
    //overflow to last bin
    xbin = std::min(xbin, histo->GetNbinsX());
    ybin = std::min(ybin, histo->GetNbinsY());
    return histo->GetBinContent(xbin, ybin);
}








// --------------------------- Methods for PileupScaleFactors ---------------------------------------------


PileupScaleFactors::PileupScaleFactors(const std::string& inputFilename,
                                       const std::string& mcEra, const std::string& pileupScenario,
                                       const Systematic::Systematic& systematic):
puReweighter_(new ztop::PUReweighter())
{
    std::cout<<"--- Beginning preparation of pileup reweighter\n";
    
    if(mcEra == "Summer12"){
        if(pileupScenario == "S10"){
            puReweighter_->setMCDistrSum12("S10");
        }
        else{
            std::cerr<<"ERROR in constructor of PileupScaleFactors! Given pileup tag is not defined (MC era, pileup tag): "
                     <<mcEra<<" , "<<pileupScenario<<"\n...break\n"<<std::endl;
            exit(33);
        }
    }
    else{
        std::cerr<<"ERROR in constructor of PileupScaleFactors! Given MC era is not defined: "
                 <<mcEra<<"\n...break\n"<<std::endl;
        exit(33);
    }
    
    TString pileupInput(common::DATA_PATH_COMMON());
    pileupInput.Append("/").Append(inputFilename);
    if(systematic.type() == Systematic::pu){
        if(systematic.variation() == Systematic::up) pileupInput.ReplaceAll(".root", "_sysUp.root");
        else if(systematic.variation() == Systematic::down) pileupInput.ReplaceAll(".root", "_sysDown.root");
    }
    std::cout<<"Using PU input file:\n"<<pileupInput<<std::endl;
    puReweighter_->setDataTruePUInput(pileupInput.Data());
    
    std::cout<<"=== Finishing preparation of pileup reweighter\n\n";
}



double PileupScaleFactors::getSF(const size_t trueVertexMultiplicity)const
{
    return puReweighter_->getPUweight(trueVertexMultiplicity);
}







// --------------------------- Methods for LeptonScaleFactors ---------------------------------------------



LeptonScaleFactors::LeptonScaleFactors(const char* electronSFInputFileName,
                                       const char* muonSFInputFileName,
                                       const Systematic::Systematic& systematic):
h2_electronSFpteta_(0),
h2_muonSFpteta_(0)
{
    std::cout<<"--- Beginning preparation of lepton scale factors\n";
    
    SystematicInternal systematicInternal(nominal);
    if(systematic.type() == Systematic::lept){
        if(systematic.variation() == Systematic::up) systematicInternal = vary_up;
        else if(systematic.variation() == Systematic::down) systematicInternal = vary_down;
    }
    
    for(const auto& lepton : {electron, muon}){
        std::string inputFileName(common::DATA_PATH_COMMON());
        inputFileName.append("/");
        std::string histogramName;
        if(lepton == electron){
            inputFileName.append(electronSFInputFileName);
            histogramName = "ElectronIdIsoSF";
            h2_electronSFpteta_ = this->prepareSF(inputFileName, histogramName, systematicInternal);
        }
        else if(lepton == muon){
            inputFileName.append(muonSFInputFileName);
            histogramName = "MuonIdIsoSF";
            h2_muonSFpteta_ = this->prepareSF(inputFileName, histogramName, systematicInternal);
        }
    }
    
    const bool electronInputFound(h2_electronSFpteta_);
    const bool muonInputFound(h2_muonSFpteta_);
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



const TH2* LeptonScaleFactors::prepareSF(const std::string& inputFileName,
                                         const std::string& histogramName,
                                         const LeptonScaleFactors::SystematicInternal& systematic)const
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



double LeptonScaleFactors::getSFDilepton(const int leadingLeptonIndex, const int nLeadingLeptonIndex,
                                         const VLV& leptons, const std::vector<int>& lepPdgIds)const
{
    if(!h2_electronSFpteta_ || !h2_muonSFpteta_) return 1.;

    const LV& leadingLepton(leptons.at(leadingLeptonIndex));
    const LV& nLeadingLepton(leptons.at(nLeadingLeptonIndex));
    const int leadingPdgId(std::abs(lepPdgIds.at(leadingLeptonIndex)));
    const int nLeadingPdgId(std::abs(lepPdgIds.at(nLeadingLeptonIndex)));

    if(leadingPdgId==11 && nLeadingPdgId==11)
        return ScaleFactorHelpers::get2DSF(h2_electronSFpteta_, leadingLepton.Eta(), leadingLepton.pt()) *
                ScaleFactorHelpers::get2DSF(h2_electronSFpteta_, nLeadingLepton.Eta(), nLeadingLepton.pt());
    if(leadingPdgId==13 && nLeadingPdgId==13)
        return ScaleFactorHelpers::get2DSF(h2_muonSFpteta_, leadingLepton.Eta(), leadingLepton.pt()) *
                ScaleFactorHelpers::get2DSF(h2_muonSFpteta_, nLeadingLepton.Eta(), nLeadingLepton.pt());
    if(leadingPdgId==13 && nLeadingPdgId==11)
        return ScaleFactorHelpers::get2DSF(h2_muonSFpteta_, leadingLepton.Eta(), leadingLepton.pt()) *
                ScaleFactorHelpers::get2DSF(h2_electronSFpteta_, nLeadingLepton.Eta(), nLeadingLepton.pt());
    if(leadingPdgId==11 && nLeadingPdgId==13)
        return ScaleFactorHelpers::get2DSF(h2_electronSFpteta_, leadingLepton.Eta(), leadingLepton.pt()) *
                ScaleFactorHelpers::get2DSF(h2_muonSFpteta_, nLeadingLepton.Eta(), nLeadingLepton.pt());
    std::cout<<"WARNING in method getLeptonIDSF! LeptonPdgIds are not as expected (pdgId1, pdgId2): "
             <<leadingPdgId<<" , "<<nLeadingPdgId<<"\n...will return scale factor = 1.\n";
    return 1.;
}



double LeptonScaleFactors::getSFAllLeptons(const std::vector<int>& allLeptonIndices,
                                           const VLV& leptons, const std::vector<int>& lepPdgIds)const
{
    if(!h2_electronSFpteta_ || !h2_muonSFpteta_) return 1.;

    double result(1.);

    for(const int index : allLeptonIndices){
        const int absPdgId(std::abs(lepPdgIds.at(index)));
        const TH2* histo(0);
        if(absPdgId==11) histo = h2_electronSFpteta_;
        else if(absPdgId==13) histo = h2_muonSFpteta_;
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
                                         const std::vector<Channel::Channel>& channels,
                                         const Systematic::Systematic& systematic):
h2_eeSFeta_(0),
h2_emuSFeta_(0),
h2_mumuSFeta_(0)
{
    std::cout<<"--- Beginning preparation of trigger scale factors\n";
    
    SystematicInternal systematicInternal(nominal);
    if(systematic.type() == Systematic::trig){
        if(systematic.variation() == Systematic::up) systematicInternal = vary_up;
        else if(systematic.variation() == Systematic::down) systematicInternal = vary_down;
    }
    
    bool someChannelFound(false);
    bool someChannelNotFound(false);
    std::stringstream channelsFound;
    std::stringstream channelsNotFound;
    
    std::vector<TString> triggerInputFileNames;
    for(const auto& channel : channels){
        const TString channelName = Channel::convert(channel);
        TString fullName(common::DATA_PATH_COMMON());
        fullName.Append("/triggerSummary_");
        fullName.Append(channelName);
        fullName.Append(inputFileSuffix);
        triggerInputFileNames.push_back(fullName);
        const TH2* const h2_triggerSF = this->prepareSF(fullName, systematicInternal);
        if(channel == Channel::ee) h2_eeSFeta_ = h2_triggerSF;
        else if(channel == Channel::emu) h2_emuSFeta_ = h2_triggerSF;
        else if(channel == Channel::mumu) h2_mumuSFeta_ = h2_triggerSF;
        else{
            std::cerr<<"ERROR in TriggerScaleFactors! Invalid channel requested: "
                     <<channelName<<"\n...break\n"<<std::endl;
            exit(23);
        }
        if(h2_triggerSF){
            someChannelFound = true;
            channelsFound<<channelName<<" , ";
        }
        else{
            someChannelNotFound = true;
            channelsNotFound<<channelName<<" , ";
        }
    }
    
    if(someChannelFound) std::cout<<"Found trigger scale factors for requested channels: "<<channelsFound.str()<<"\n";
    if(someChannelNotFound)  std::cout<<"Could NOT find trigger scale factors for requested channels: "<<channelsNotFound.str()<<"\n";
    
    std::cout<<"=== Finishing preparation of trigger scale factors\n\n";
}



const TH2* TriggerScaleFactors::prepareSF(const TString& fileName, const SystematicInternal& systematic)const
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

        for(int i = 1; i <= h_TrigSFeta->GetNbinsX(); ++i){
            for(int j = 1; j <= h_TrigSFeta->GetNbinsY(); ++j){
                h_TrigSFeta->SetBinContent(i, j, h_TrigSFeta->GetBinContent(i,j) + factor*h_TrigSFeta->GetBinError(i,j));
            }
        }
    }
    
    h_TrigSFeta->SetDirectory(0);
    trigEfficiencies.Close();
    
    return h_TrigSFeta;
}



double TriggerScaleFactors::getSF(const int leptonXIndex, const int leptonYIndex,
                                  const VLV& leptons, const Channel::Channel& channel)const
{
    const TH2* h_TrigSFeta(0);
    if(channel == Channel::ee) h_TrigSFeta = h2_eeSFeta_;
    else if(channel == Channel::emu) h_TrigSFeta = h2_emuSFeta_;
    else if(channel == Channel::mumu) h_TrigSFeta = h2_mumuSFeta_;
    else{
        std::cerr<<"ERROR in TriggerScaleFactors! Invalid channel requested: "
                 <<Channel::convert(channel)<<"\n...break\n"<<std::endl;
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
                                   const std::vector<Channel::Channel>& channels,
                                   const Systematic::Systematic& systematic,
                                   const CorrectionMode& correctionMode):
bTagBase(),
selectorList_(0),
channel_(Channel::undefined),
correctionMode_(correctionMode)
{
    std::cout<<"--- Beginning preparation of b-tagging scale factors\n";
    
    // Set systematic if it is an allowed one for btag efficiencies, else set to nominal
    ztop::bTagBase::systematics systematicInternal(nominal);
    if(systematic.type() == Systematic::btag){
        if(systematic.variation() == Systematic::up) systematicInternal = heavyup;
        else if(systematic.variation() == Systematic::down) systematicInternal = heavydown;
    }
    if(systematic.type() == Systematic::btagPt){
        if(systematic.variation() == Systematic::up) systematicInternal = heavyuppt;
        else if(systematic.variation() == Systematic::down) systematicInternal = heavydownpt;
    }
    if(systematic.type() == Systematic::btagEta){
        if(systematic.variation() == Systematic::up) systematicInternal = heavyupeta;
        else if(systematic.variation() == Systematic::down) systematicInternal = heavydowneta;
    }
    if(systematic.type() == Systematic::btagLjet){
        if(systematic.variation() == Systematic::up) systematicInternal = lightup;
        else if(systematic.variation() == Systematic::down) systematicInternal = lightdown;
    }
    if(systematic.type() == Systematic::btagLjetPt){
        if(systematic.variation() == Systematic::up) systematicInternal = lightuppt;
        else if(systematic.variation() == Systematic::down) systematicInternal = lightdownpt;
    }
    if(systematic.type() == Systematic::btagLjetEta){
        if(systematic.variation() == Systematic::up) systematicInternal = lightupeta;
        else if(systematic.variation() == Systematic::down) systematicInternal = lightdowneta;
    }
    this->setSystematic(systematicInternal);
    
    // FIXME: need to adapt behaviour for requested correction mode
    //if(correctionMode_ == none) ;
    //else if(correctionMode_==greaterEqualOneTagReweight || correctionMode_==randomNumberRetag) ;
    //else if(correctionMode_==discriminatorReweight) ;
    //else ;
    
    // Hardcoded filename, since btag efficiencies are produced and used always from this one
    const TString filename("ttbarsignalplustau.root"); 
    
    // Set the sample name for each channel
    for(const auto& channel : channels) m_channelSamplename_[channel] = Channel::convert(channel)+"_"+filename;
    
    // Checking whether files with btag efficiencies exist for all channels
    bool allInputFilesAvailable(true);
    for(const auto& channel : channels){
        TString btagInputFile = common::accessFolder(btagEfficiencyInputDir, channel, systematic, true);
        if(btagInputFile == "") allInputFilesAvailable = false;
        btagInputFile.Append(m_channelSamplename_.at(channel));
        
        ifstream inputFileStream;
        if(allInputFilesAvailable) inputFileStream.open(btagInputFile);
        if(inputFileStream.is_open()){
            // Setting the file and sample name for each channel in the map if the file exists
            m_channelFilename_[channel] = btagInputFile;
            inputFileStream.close();
        }
        else{
            std::cout<< "******************************************************\n"
                    << "Btag efficiency file [" << btagInputFile << "] doesn't exist.\n"
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
        }
    }

    if(!allInputFilesAvailable){
        std::cout<<"Not all input files for b-tagging efficiencies available\n"
                 <<"\t-->  Efficiencies will not be used, but produced in Analysis\n";
        this->setMakeEff(true);
        // Setting the root file names for storing btagging efficiencies for each channel in the map
        for(const auto& channel : channels){
            const TString btagOutputFile = common::assignFolder(btagEfficiencyOutputDir, channel, systematic).Append(m_channelSamplename_.at(channel));
            m_channelFilename_[channel] = btagOutputFile;
        }
    }
    else{
        std::cout<<"Found all input files for b-tagging efficiencies\n";
        this->setMakeEff(false);
    }

    std::cout<<"=== Finishing preparation of b-tagging scale factors\n\n";
}



void BtagScaleFactors::prepareSF(const Channel::Channel& channel)
{
    channel_ = channel;
    
    // Stopping if no efficiency histograms need to be read (production mode)
    if(this->makeEfficiencies()) return;
    
    std::vector<TH2D> histos;
    std::vector<TH2D> effhistos;
    std::vector<float> medians;

    // Load per-jet efficiencies file
    const TString& inputFileName = m_channelFilename_.at(channel_);
    TFile file(inputFileName,"READ");

    // Disabling referencing histograms to gDirectory to prevent crashing when closing the root file
    TH1::AddDirectory(false);
    TH2D* tempHisto(0);

    //read the histograms in the right order from file. this is sufficient, since getJetHistoOrderedNames() and
    //getEffHistoOrderedNames() take care of the right index <-> name association
    //get filled histograms
    for(size_t i = 0; i < this->getJetHistoOrderedNames().size(); ++i){
        file.GetObject(this->getJetHistoOrderedNames().at(i).c_str(), tempHisto);
        if(!tempHisto)
            throw std::runtime_error("BtagScaleFactors::prepareBtagSF Couldn't find all histograms in the input root file");
        histos.push_back(*tempHisto);
        tempHisto = 0;
    }
    
    //get efficiency histograms
    for(size_t i = 0; i < this->getEffHistoOrderedNames().size(); ++i){
        file.GetObject(this->getEffHistoOrderedNames().at(i).c_str(), tempHisto);
        if(!tempHisto)
            throw std::runtime_error("BtagScaleFactors::prepareBtagSF Couldn't find all histograms in the input root file");
        effhistos.push_back(*tempHisto);
        tempHisto = 0;
    }
    
    // Reading the tag and efficiency histograms from the file and adding them to the maps in proper order
    /* nazars impl  for(int id = 0; ; ++id) {
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
     */
    
    // Reading and extracting madian values from the histogram
    TH1* medianHisto(0);
    file.GetObject("medians", medianHisto);
    if(!medianHisto) throw std::runtime_error("BtagScaleFactors::prepareBTags Couldn't find [medians] histogram in the input root file");
    for(int i = 0; i < medianHisto->GetNbinsX(); ++i) medians.push_back((float)medianHisto->GetBinContent(i+1));

    const TString& sampleName = m_channelSamplename_.at(channel_);
    histos_[sampleName.Data()] = histos;
    effhistos_[sampleName.Data()] = effhistos;
    medianMap_[sampleName.Data()] = medians;

    file.Close();

    if(this->setSampleName(sampleName.Data()) < 0)
        throw std::runtime_error("BtagScaleFactors::prepareBTags Tried to set a non-existing sampleName");
}



bool BtagScaleFactors::makeEfficiencies()
{
    return this->getMakeEff();
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
        const LV& jet = jets.at(index);

        // Preparing a seed for the random jet retagging
        const unsigned int seed = std::abs( static_cast<int>( 1.e6*sin( 1.e6*jet.Phi() ) ) );
        const bool isTagged = this->jetIsTagged( jet.pt(), std::fabs(jet.eta()), jetPartonFlavours.at(index), btagDiscriminants.at(index), seed );
        if(isTagged) tagged_indices.push_back(index);
    }
    bjetIndices = tagged_indices;
}



double BtagScaleFactors::getSFGreaterEqualOneTag(const std::vector<int>& jetIndices,
                                                 const VLV& jets,
                                                 const std::vector<int>& jetPartonFlavours)
{
    this->resetCounter();
    for(const int index : jetIndices){
        this->countJet(jets.at(index).pt(), std::abs(jets.at(index).eta()), jetPartonFlavours.at(index));
    }

    // per-event SF calculation
    //double scale_factor = getEventSF();

    //this is plain wrong and especially leads to underestimated uncertainties!
    //if(std::abs(scale_factor-1.) > 0.05) scale_factor = 1.;

    return this->getEventSF();
}



void BtagScaleFactors::bookHistograms(TSelectorList* output)
{
    // Set pointer to output, so that selectorList_histograms are owned by it
    // (not used in current implementation)
    selectorList_ = output;

    const TString& sampleName = m_channelSamplename_.at(channel_);
    this->setSampleName(sampleName.Data());
}



void BtagScaleFactors::fillHistograms(const std::vector<int>& jetIndices,
                                      const std::vector<double>& bTagDiscriminant,
                                      const VLV& jets,
                                      const std::vector<int>& jetPartonFlavours,
                                      const double& weight)
{
    for(const int index : jetIndices){
        this->fillEff(jets.at(index).pt(), std::abs(jets.at(index).eta()),
            std::abs(jetPartonFlavours.at(index)), bTagDiscriminant.at(index), weight);
    }
}



void BtagScaleFactors::produceEfficiencies()
{
    const TString& outputFileName = m_channelFilename_.at(channel_);
    const TString& sampleName = m_channelSamplename_.at(channel_);

    TFile file(outputFileName, "RECREATE");

    // Creating the histograms
    this->makeEffs();

    // Writing each histogram to file
    const std::vector<TH2D>& histos = histos_.at(sampleName.Data());
    for(const TH2D& histo : histos) {
        histo.Write();
    }

    // Writing each efficiency hitogram to file
    const std::vector<TH2D>& effhistos = effhistos_.at(sampleName.Data());
    for(const TH2D& histo : effhistos){
        histo.Write();
    }

    // Writing medians to the histogram and then to the file
    const size_t nMedians = (size_t)medians::length_median;
    TH1D histo("medians", "Medians;Property id;Median value", nMedians, 0, nMedians);

    if(nMedians != medianMap_.at(sampleName.Data()).size()){
        throw std::range_error("BtagScaleFactors::produceBtagEfficiencies Numbers of stored and designed median values differ");
    }

    for(size_t i = 0; i<nMedians; ++i){
        histo.SetBinContent(i+1, medianMap_.at(sampleName.Data()).at(i));
    }
    histo.Write();
    
    file.Close();
    
    std::cout<<"Done with production of b-tag efficiency file: "<<outputFileName<<"\n\n"<<std::endl;
}








// --------------------------- Methods for JetEnergyResolutionScaleFactors ---------------------------------------------




JetEnergyResolutionScaleFactors::JetEnergyResolutionScaleFactors(const Systematic::Systematic& systematic)
{
    std::cout<<"--- Beginning preparation of JER scale factors\n";
    
    SystematicInternal systematicInternal(undefined);
    if(systematic.type() == Systematic::jer){
        if(systematic.variation() == Systematic::up) systematicInternal = vary_up;
        else if(systematic.variation() == Systematic::down) systematicInternal = vary_down;
    }
    else{
        std::cerr<<"ERROR in constructor of JetEnergyResolutionScaleFactors! Systematic is invalid: "
                 <<Systematic::convertType(systematic.type())<<"\n...break\n"<<std::endl;
        exit(98);
    }
    
    // Hardcoded eta ranges
    v_etaRange_ = {0.5, 1.1, 1.7, 2.3, 10.};
    
    // Hardcoded scale factors for eta ranges, nominal is = {1.052, 1.057, 1.096, 1.134, 1.288};
    if(systematicInternal == vary_up){
        v_etaScaleFactor_ = {1.115, 1.114, 1.161, 1.228, 1.488};
        std::cout<<"Apply systematic variation: up\n";
    }
    else if(systematicInternal == vary_down){
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
                                                       const Systematic::Systematic& systematic):
jetCorrectionUncertainty_(0),
varyUp_(false)
{
    std::cout<<"--- Beginning preparation of JES scale factors\n";
    
    SystematicInternal systematicInternal(undefined);
    if(systematic.type() == Systematic::jes){
        if(systematic.variation() == Systematic::up) systematicInternal = vary_up;
        else if(systematic.variation() == Systematic::down) systematicInternal = vary_down;
    }
    else{
        std::cerr<<"ERROR in constructor of JetEnergyScaleScaleFactors! Systematic is invalid: "
                 <<Systematic::convertType(systematic.type())<<"\n...break\n"<<std::endl;
        exit(98);
    }
    
    std::string inputFileName(common::DATA_PATH_COMMON());
    inputFileName.append("/");
    inputFileName.append(jesUncertaintySourceFile);
    jetCorrectionUncertainty_ = new ztop::JetCorrectionUncertainty(ztop::JetCorrectorParameters(inputFileName.data(), "Total"));
    
    if(systematicInternal == vary_up){
        varyUp_ = true;
        std::cout<<"Apply systematic variation: up\n";
    }
    else if(systematicInternal == vary_down){
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






