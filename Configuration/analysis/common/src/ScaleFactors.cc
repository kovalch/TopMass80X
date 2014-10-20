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
h2_mumuSFeta_(0),
h2_channelSF_(0)
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



void TriggerScaleFactors::prepareChannel(const Channel::Channel& channel)
{
    if(channel == Channel::ee) h2_channelSF_ = h2_eeSFeta_;
    else if(channel == Channel::emu) h2_channelSF_ = h2_emuSFeta_;
    else if(channel == Channel::mumu) h2_channelSF_ = h2_mumuSFeta_;
    else{
        std::cerr<<"ERROR in TriggerScaleFactors::prepareChannel()! Invalid channel requested: "
                 <<Channel::convert(channel)<<"\n...break\n"<<std::endl;
        exit(25);
    }
}



double TriggerScaleFactors::getSF(const int leptonXIndex, const int leptonYIndex,
                                  const VLV& leptons)const
{
    if (!h2_channelSF_) return 1.;

    //For 'ee' and 'mumu' channels Xaxis of the 2D plots is the highest pT lepton
    // for the 'emu' channel Xaxis is the electron and Y axis muon
    const LV& leptonX(leptons.at(leptonXIndex));
    const LV& leptonY(leptons.at(leptonYIndex));
    return ScaleFactorHelpers::get2DSF(h2_channelSF_, std::fabs(leptonX.eta()), std::fabs(leptonY.eta()));
}















// --------------------------- Methods for BtagScaleFactors ---------------------------------------------




BtagScaleFactors::BtagScaleFactors(const char* efficiencyInputDir,
                                   const char* efficiencyOutputDir,
                                   const char* inputFileHeavyFlavour,
                                   const char* inputFileLightFlavour,
                                   const std::vector<Channel::Channel>& channels,
                                   const Systematic::Systematic& systematic,
                                   const CorrectionMode& correctionMode):
bTagBase(),
selectorList_(0),
channel_(Channel::undefined),
correctionMode_(correctionMode)
{
    std::cout<<"--- Beginning preparation of b-tagging scale factors\n";
    
    // Set systematic if it is an allowed one for b-tag efficiencies, else set to nominal
    this->btagSystematic(systematic);
    
    // Adapt behaviour for requested correction mode
    this->setMakeEff(false);
    if(correctionMode_ == none){
        std::cout<<"Correction mode: none\n";
        std::cout<<"No corrections will be applied\n";
    }
    else if(correctionMode_==greaterEqualOneTagReweight || correctionMode_==randomNumberRetag){
        if(correctionMode_ == greaterEqualOneTagReweight)
            std::cout<<"Correction mode: Event scale factor correcting (mis-)tagging efficiency for >=1 b-tagged jet\n";
        else if(correctionMode_ == randomNumberRetag)
            std::cout<<"Correction mode: Random number based (un-)tagging of b-/c-/l- jets correcting b-tag multiplicity\n";
        std::cout<<"Requires b-/c-/l-tagging efficiency histograms from MC\n";
        this->prepareEfficiencies(efficiencyInputDir, efficiencyOutputDir, channels, systematic);
    }
    else if(correctionMode_ == discriminatorReweight){
        std::cout<<"Correction mode: Event scale factor correcting the b-tag discriminator distribution\n";
        std::cout<<"Requires files containing the official flavour-specific per-jet correction factors\n";
        this->prepareDiscriminatorReweighting(inputFileHeavyFlavour, inputFileLightFlavour);
    }
    else{
        std::cerr<<"ERROR in constructor of BtagScaleFactors! Requested a correction mode which is not implemented\n...break\n"<<std::endl;
        exit(90);
    }

    std::cout<<"=== Finishing preparation of b-tagging scale factors\n\n";
}



void BtagScaleFactors::btagSystematic(const Systematic::Systematic& systematic)
{
    // Check if given systematic should be processed, i.e. whether it is valid for the given correction mode, else stop
    if(std::find(Systematic::btagTypes.begin(), Systematic::btagTypes.end(), systematic.type()) != Systematic::btagTypes.end()){
        std::ostringstream stream;
        stream<<"Correction mode of b-tag not compatible with specified systematic: "<<systematic.name()
              <<"\nStop running of analysis --- this is NOT an error, just avoiding double running\n";
        if(correctionMode_ == none){
            std::cout<<stream.str()<<std::endl;
            exit(111);
        }
        else if(correctionMode_==greaterEqualOneTagReweight || correctionMode_==randomNumberRetag){
            const std::vector<Systematic::Type>& efficiencyTypes = Systematic::btagEfficiencyCorrectionTypes;
            if(std::find(efficiencyTypes.begin(), efficiencyTypes.end(), systematic.type()) == efficiencyTypes.end()){
                std::cout<<stream.str()<<std::endl;
                exit(111);
            }
        }
        else if(correctionMode_ == discriminatorReweight){
            const std::vector<Systematic::Type>& discriminatorTypes = Systematic::btagDiscriminatorReweightTypes;
            if(std::find(discriminatorTypes.begin(), discriminatorTypes.end(), systematic.type()) == discriminatorTypes.end()){
                std::cout<<stream.str()<<std::endl;
                exit(111);
            }
        }
    }
    
    // Assign the internal systematic
    ztop::bTagBase::systematics systematicInternal(nominal);
    if(systematic.type() == Systematic::btag){
        if(systematic.variation() == Systematic::up) systematicInternal = heavyup;
        else if(systematic.variation() == Systematic::down) systematicInternal = heavydown;
    }
    else if(systematic.type() == Systematic::btagPt){
        if(systematic.variation() == Systematic::up) systematicInternal = heavyuppt;
        else if(systematic.variation() == Systematic::down) systematicInternal = heavydownpt;
    }
    else if(systematic.type() == Systematic::btagEta){
        if(systematic.variation() == Systematic::up) systematicInternal = heavyupeta;
        else if(systematic.variation() == Systematic::down) systematicInternal = heavydowneta;
    }
    else if(systematic.type() == Systematic::btagLjet){
        if(systematic.variation() == Systematic::up) systematicInternal = lightup;
        else if(systematic.variation() == Systematic::down) systematicInternal = lightdown;
    }
    else if(systematic.type() == Systematic::btagLjetPt){
        if(systematic.variation() == Systematic::up) systematicInternal = lightuppt;
        else if(systematic.variation() == Systematic::down) systematicInternal = lightdownpt;
    }
    else if(systematic.type() == Systematic::btagLjetEta){
        if(systematic.variation() == Systematic::up) systematicInternal = lightupeta;
        else if(systematic.variation() == Systematic::down) systematicInternal = lightdowneta;
    }
    else if(systematic.type() == Systematic::btagDiscrBstat1){
        if(systematic.variation() == Systematic::up) systematicInternal = hfstat1up;
        else if(systematic.variation() == Systematic::down) systematicInternal = hfstat1down;
    }
    else if(systematic.type() == Systematic::btagDiscrBstat2){
        if(systematic.variation() == Systematic::up) systematicInternal = hfstat2up;
        else if(systematic.variation() == Systematic::down) systematicInternal = hfstat2down;
    }
    else if(systematic.type() == Systematic::btagDiscrLstat1){
        if(systematic.variation() == Systematic::up) systematicInternal = lfstat1up;
        else if(systematic.variation() == Systematic::down) systematicInternal = lfstat1down;
    }
    else if(systematic.type() == Systematic::btagDiscrLstat2){
        if(systematic.variation() == Systematic::up) systematicInternal = lfstat2up;
        else if(systematic.variation() == Systematic::down) systematicInternal = lfstat2down;
    }
    else if(systematic.type() == Systematic::btagDiscrPurity){
        if(systematic.variation() == Systematic::up) systematicInternal = purityup;
        else if(systematic.variation() == Systematic::down) systematicInternal = puritydown;
    }
    else if(systematic.type() == Systematic::jes){
        // This variation covers the with JES correlated uncertainty, only needed for discriminator reweighting
        if(correctionMode_ == discriminatorReweight){
            if(systematic.variation() == Systematic::up) systematicInternal = jesup;
            else if(systematic.variation() == Systematic::down) systematicInternal = jesdown;
        }
    }
    this->setSystematic(systematicInternal);
}



void BtagScaleFactors::prepareChannel(const Channel::Channel& channel)
{
    channel_ = channel;
    
    if(correctionMode_==randomNumberRetag || correctionMode_==greaterEqualOneTagReweight) this->loadEfficiencies();
}



void BtagScaleFactors::algorithmAndWorkingPoint(const Btag::Algorithm& algorithm,
                                                const Btag::WorkingPoint& workingPoint)
{
    if(algorithm == Btag::csv){
        if(workingPoint == Btag::L) this->setWorkingPoint(csvl_wp);
        else if(workingPoint == Btag::M) this->setWorkingPoint(csvm_wp);
        else if(workingPoint == Btag::T) this->setWorkingPoint(csvt_wp);
        else{
            std::cerr<<"ERROR in BtagScaleFactors::algorithmAndWorkingPoint()! Working point is not implemented for algorithm 'csv'\n...break\n"<<std::endl;
            exit(232);
        }
    }
    else{
        std::cerr<<"ERROR in BtagScaleFactors::algorithmAndWorkingPoint()! Algorithm is not implemented\n...break\n"<<std::endl;
        exit(232);
    }
}



void BtagScaleFactors::indexOfBtags(std::vector<int>& bjetIndices,
                                    const std::vector<int>& jetIndices,
                                    const VLV& jets,
                                    const std::vector<int>& jetPartonFlavours,
                                    const std::vector<double>& btagDiscriminants)const
{
    // In case efficiencies are produced, the input needs to be without any retagging
    if(this->makeEfficiencies()) return;
    
    if(correctionMode_ != randomNumberRetag) return;
    
    bjetIndices.clear();
    for(const int index : jetIndices){
        // Skip jets where there is no partonFlavour
        if(jetPartonFlavours.at(index) == 0) continue;
        const LV& jet = jets.at(index);

        // Preparing a seed for the random jet retagging
        const unsigned int seed = std::abs( static_cast<int>( 1.e6*sin( 1.e6*jet.Phi() ) ) );
        const bool isTagged = this->jetIsTagged( jet.pt(), std::fabs(jet.eta()), jetPartonFlavours.at(index), btagDiscriminants.at(index), seed );
        if(isTagged) bjetIndices.push_back(index);
    }
}



double BtagScaleFactors::getSF(const std::vector<int>& jetIndices,
                               const VLV& jets,
                               const std::vector<int>& jetPartonFlavours,
                               const std::vector<double>& btagDiscriminators)
{
    // In case efficiencies are produced, the input needs to be without scale factors
    if(this->makeEfficiencies()) return 1.;
    
    if(correctionMode_ == greaterEqualOneTagReweight)
        return this->scaleFactorGreaterEqualOneTag(jetIndices, jets, jetPartonFlavours);
    else if(correctionMode_ == discriminatorReweight)
        return this->scaleFactorDiscriminatorReweight(jetIndices, jets, jetPartonFlavours, btagDiscriminators);
    else return 1.;
}



double BtagScaleFactors::scaleFactorDiscriminatorReweight(const std::vector<int>& jetIndices,
                                                          const VLV& jets,
                                                          const std::vector<int>& jetPartonFlavours,
                                                          const std::vector<double>& btagDiscriminators)const
{
    double result(1.);
    
    // Loop over all jets and multiply the individual weights
    for(const int index : jetIndices){
        const LV& jet = jets.at(index);
        const int& jetPartonFlavour = jetPartonFlavours.at(index);
        const double& btagDiscriminator = btagDiscriminators.at(index);
        const double weight = this->getJetDiscrShapeWeight(jet.pt(), std::abs(jet.eta()), jetPartonFlavour, btagDiscriminator);
        result *= weight;
    }
    
    return result;
}



double BtagScaleFactors::scaleFactorGreaterEqualOneTag(const std::vector<int>& jetIndices,
                                                       const VLV& jets,
                                                       const std::vector<int>& jetPartonFlavours)
{
    this->resetCounter();
    for(const int index : jetIndices){
        this->countJet(jets.at(index).pt(), std::abs(jets.at(index).eta()), jetPartonFlavours.at(index));
    }

    // per-event SF calculation
    //this is plain wrong and especially leads to underestimated uncertainties!
    //double scale_factor = this->getEventSF();
    //if(std::abs(scale_factor-1.) > 0.05) scale_factor = 1.;

    return this->getEventSF();
}



bool BtagScaleFactors::makeEfficiencies()const
{
    return this->getMakeEff();
}



void BtagScaleFactors::prepareEfficiencies(const char* efficiencyInputDir,
                                           const char* efficiencyOutputDir,
                                           const std::vector<Channel::Channel>& channels,
                                           const Systematic::Systematic& systematic)
{
    // Hardcoded filename, since btag efficiencies are produced and used always from this one
    const TString filename("ttbarsignalplustau.root"); 
    
    // Set the sample name for each channel
    for(const auto& channel : channels) m_channelSamplename_[channel] = Channel::convert(channel)+"_"+filename;
    
    // Checking whether files with btag efficiencies exist for all channels
    bool allInputFilesAvailable(true);
    for(const auto& channel : channels){
        TString btagInputFile = common::accessFolder(efficiencyInputDir, channel, systematic, true);
        if(btagInputFile == "") allInputFilesAvailable = false;
        btagInputFile.Append(filename);
        
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
                    << "To create the file, run (for each systematic 'SYST'):\n"
                    << "\t> ./install/bin/load_Analysis -f ttbarsignalplustau.root -c emu -s SYST\n"
                    << "\t> ./install/bin/load_Analysis -f ttbarsignalplustau.root -c ee -s SYST\n"
                    << "\t> ./install/bin/load_Analysis -f ttbarsignalplustau.root -c mumu -s SYST\n"
                    << " and move the selectionRoot/BTagEff directory to the cwd:\n"
                    << "\t> mv selectionRoot/BTagEff .\n"
                    << "WARNING: Analysis results are not produced\n"
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
            const TString btagOutputFile = common::assignFolder(efficiencyOutputDir, channel, systematic).Append(filename);
            m_channelFilename_[channel] = btagOutputFile;
        }
    }
    else{
        std::cout<<"Found all input files for b-tagging efficiencies\n";
    }
}



void BtagScaleFactors::loadEfficiencies()
{
    // Stop if no efficiency histograms need to be read (production mode)
    if(this->makeEfficiencies()) return;
    
    // Load per-jet efficiencies file
    const TString& inputFileName = m_channelFilename_.at(channel_);
    TFile file(inputFileName,"READ");

    // Disabling referencing histograms to gDirectory to prevent crashing when closing the root file
    TH1::AddDirectory(false);
    TH2D* tempHisto(0);

    // Read the histograms for jets and b-jets of different flavour in the right order from file
    // getJetHistoOrderedNames() takes care of the right index <-> name association
    std::vector<TH2D> histos;
    for(size_t iHisto = 0; iHisto < this->getJetHistoOrderedNames().size(); ++iHisto){
        file.GetObject(this->getJetHistoOrderedNames().at(iHisto).c_str(), tempHisto);
        if(!tempHisto) throw std::runtime_error("ERROR in BtagScaleFactors::prepareBtagSF()! Couldn't find all histograms in the input root file");
        histos.push_back(*tempHisto);
        tempHisto = 0;
    }
    
    // Read efficiency histograms in the right order from file
    // getEffHistoOrderedNames() takes care of the right index <-> name association
    std::vector<TH2D> effhistos;
    for(size_t iHisto = 0; iHisto < this->getEffHistoOrderedNames().size(); ++iHisto){
        file.GetObject(this->getEffHistoOrderedNames().at(iHisto).c_str(), tempHisto);
        if(!tempHisto) throw std::runtime_error("BtagScaleFactors::prepareBtagSF Couldn't find all histograms in the input root file");
        effhistos.push_back(*tempHisto);
        tempHisto = 0;
    }
    
    // Extract median values from histogram
    std::vector<float> medians;
    TH1* medianHisto(0);
    file.GetObject("medians", medianHisto);
    if(!medianHisto) throw std::runtime_error("BtagScaleFactors::prepareBTags Couldn't find [medians] histogram in the input root file");
    for(int iBin = 0; iBin < medianHisto->GetNbinsX(); ++iBin) medians.push_back((float)medianHisto->GetBinContent(iBin+1));
    
    // Set data members for using efficiency in correction tools
    const std::string& sampleName = m_channelSamplename_.at(channel_);
    histos_[sampleName] = histos;
    effhistos_[sampleName] = effhistos;
    medianMap_[sampleName] = medians;
    if(this->setSampleName(sampleName) < 0) throw std::runtime_error("ERROR in BtagScaleFactors::prepareSF()! Tried to set a non-existing sampleName");
    
    // Cleanup
    file.Close();
}



void BtagScaleFactors::bookEfficiencyHistograms(TSelectorList* output)
{
    // If no efficiencies should be produced, return
    if(!this->makeEfficiencies()) return;
    
    // Set pointer to output, so that selectorList_ histograms are owned by it
    // (not used in current implementation)
    selectorList_ = output;
    
    // Book histograms
    this->setSampleName(m_channelSamplename_.at(channel_));
}



void BtagScaleFactors::fillEfficiencyHistograms(const std::vector<int>& jetIndices,
                                                const std::vector<double>& btagDiscriminators,
                                                const VLV& jets,
                                                const std::vector<int>& jetPartonFlavours,
                                                const double& weight)
{
    // If no efficiencies should be produced, return
    if(!this->makeEfficiencies()) return;
    
    for(const int index : jetIndices){
        this->fillEff(jets.at(index).pt(), std::abs(jets.at(index).eta()),
                      std::abs(jetPartonFlavours.at(index)), btagDiscriminators.at(index), weight);
    }
}



void BtagScaleFactors::produceEfficiencies()
{
    // If no efficiencies should be produced, return
    if(!this->makeEfficiencies()) return;
    
    const TString& outputFileName = m_channelFilename_.at(channel_);
    TFile file(outputFileName, "RECREATE");

    // Creating the histograms
    this->makeEffs();

    const std::string& sampleName = m_channelSamplename_.at(channel_);

    // Writing each jet and b-jet histogram to file
    const std::vector<TH2D>& histos = histos_.at(sampleName);
    for(const TH2D& histo : histos) histo.Write();

    // Writing each efficiency hitogram to file
    const std::vector<TH2D>& effhistos = effhistos_.at(sampleName);
    for(const TH2D& histo : effhistos) histo.Write();

    // Writing medians to the histogram and then to the file
    const size_t nMedians = (size_t)medians::length_median;
    TH1D histo("medians", "Medians;Property id;Median value", nMedians, 0, nMedians);
    if(nMedians != medianMap_.at(sampleName).size())
        throw std::range_error("ERROR in BtagScaleFactors::produceBtagEfficiencies()! Numbers of stored and designed median values differ");
    for(size_t i = 0; i < nMedians; ++i) histo.SetBinContent(i+1, medianMap_.at(sampleName).at(i));
    
    // Cleanup
    histo.Write();
    file.Close();
    
    std::cout<<"Done with production of b-tag efficiency file: "<<outputFileName<<"\n\n"<<std::endl;
}



void BtagScaleFactors::prepareDiscriminatorReweighting(const char* inputFileHeavyFlavour, const char* inputFileLightFlavour)
{
    TString inputFolder = common::DATA_PATH_COMMON();
    inputFolder.Append("/");
    this->readShapeReweightingFiles(inputFolder+inputFileHeavyFlavour, inputFolder+inputFileLightFlavour);
    std::cout<<"Found all input files for scale factors\n";
}








// --------------------------- Methods for JetEnergyResolutionScaleFactors ---------------------------------------------




JetEnergyResolutionScaleFactors::JetEnergyResolutionScaleFactors(const char* scaleFactorSource,
                                                                 const Systematic::Systematic& systematic)
{
    std::cout<<"--- Beginning preparation of JER scale factors\n";
    
    // Set internal systematic
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
    
    // Sanity check
    std::cout<<"Using scale factor source: "<<scaleFactorSource<<"\n";
    if(systematicInternal == vary_up){
        std::cout<<"Apply systematic variation: up\n";
    }
    else if(systematicInternal == vary_down){
        std::cout<<"Apply systematic variation: down\n";
    }
    else{
        std::cerr<<"Error in constructor of JetEnergyResolutionScaleFactors! Systematic variation not allowed\n...break\n"<<std::endl;
        exit(621);
    }
    
    // Set scale factors corresponding to specified source
    this->scaleFactors(scaleFactorSource, systematicInternal);
    
    
    // Check correct size
    if(v_etaRange_.size() != v_etaScaleFactor_.size()){
        std::cerr<<"Error in constructor of JetEnergyResolutionScaleFactors! "
                 <<"Different number of intervals in eta and corresponding scale factors (intervals, scale factors): "
                 <<v_etaRange_.size()<<" , "<<v_etaScaleFactor_.size()<<"\n...break\n"<<std::endl;
        exit(622);
    }
    
    std::cout<<"=== Finishing preparation of JER scale factors\n\n";
}



void JetEnergyResolutionScaleFactors::scaleFactors(const std::string& scaleFactorSource,
                                                   const SystematicInternal& systematicInternal)
{
    // Scale factor values for all sources
    if(scaleFactorSource == "jer2011"){
        // Eta ranges
        v_etaRange_ = {0.5, 1.1, 1.7, 2.3, 10.};
        // Scale factors for eta ranges, nominal is = {1.052, 1.057, 1.096, 1.134, 1.288};
        if(systematicInternal == vary_up)
            v_etaScaleFactor_ = {1.115, 1.114, 1.161, 1.228, 1.488};
        else if(systematicInternal == vary_down)
            v_etaScaleFactor_ = {0.990, 1.001, 1.032, 1.042, 1.089};
    }
    else if(scaleFactorSource == "jer2012"){
        // Eta ranges
        v_etaRange_ = {0.5, 1.1, 1.7, 2.3, 2.8, 3.2, 5.0};
        // Scale factors for eta ranges, nominal is = {1.079, 1.099, 1.121, 1.208, 1.254, 1.395, 1.056};
        if(systematicInternal == vary_up)
            v_etaScaleFactor_ = {1.105, 1.127, 1.150, 1.254, 1.316, 1.458, 1.247};
        else if(systematicInternal == vary_down)
            v_etaScaleFactor_ = {1.053, 1.071, 1.092, 1.162, 1.192, 1.332, 0.865};
    }
    else{
        std::cerr<<"ERROR in JetEnergyResolutionScaleFactors::scaleFactors()! Name of scale factor source not defined: "
                 <<scaleFactorSource<<"\n...break\n"<<std::endl;
        exit(883);
    }
}



void JetEnergyResolutionScaleFactors::applyJetSystematic(VLV* const v_jet,
                                                         const std::vector<double>* const v_jetJerSF,
                                                         const VLV* const v_associatedGenJet)const
{
    // This loop corrects the jet collection used for jet selections
    for(size_t iJet = 0; iJet < v_jet->size(); ++iJet){
        LV& jet = v_jet->at(iJet);
        const double& jetJerSF = v_jetJerSF->at(iJet);
        const LV& associatedGenJet = v_associatedGenJet->at(iJet);
        
        this->scaleJet(jet, jetJerSF, associatedGenJet);
    }
}



void JetEnergyResolutionScaleFactors::applyMetSystematic(VLV* const v_jetForMet, LV* const met,
                                                         const std::vector<double>* const v_jetForMetJerSF,
                                                         const VLV* const v_associatedGenJetForMet)const
{
    // This loop corrects the jet collection used to modify the MET
    double deltaPx = 0.;
    double deltaPy = 0.;
    for(size_t iJet = 0; iJet < v_jetForMet->size(); ++iJet){
        LV& jet = v_jetForMet->at(iJet);
        const double& jetJerSF = v_jetForMetJerSF->at(iJet);
        const LV& associatedGenJet = v_associatedGenJetForMet->at(iJet);
        
        const double storedPx = jet.px();
        const double storedPy = jet.py();
        if(this->scaleJet(jet, jetJerSF, associatedGenJet)){
            deltaPx += jet.px() - storedPx;
            deltaPy += jet.py() - storedPy;
        }
    }
    
    // Adjust the MET
    const double scaledMetPx = met->px() - deltaPx;
    const double scaledMetPy = met->py() - deltaPy;
    met->SetPt(std::sqrt(scaledMetPx*scaledMetPx + scaledMetPy*scaledMetPy));
}



bool JetEnergyResolutionScaleFactors::scaleJet(LV& jet, const double& jetJerSF, const LV& associatedGenJet)const
{
    // Avoid division by 0, and do not scale if original scale factor is 1.
    if(jetJerSF==0. || jetJerSF==1.) return false;
    
    // Do not scale if no associated genJet is found
    if(associatedGenJet.pt() == 0.) return false;
    
    // Do not scale if jet is outside range defined for scale factors (i.e. in a range where no scale factors are derived)
    const int jetEtaBin = this->jetEtaBin(jet);
    if(jetEtaBin == -1) return false;
    
    // Scale the jet
    jet *= 1./jetJerSF;
    const double factor = 1. + (v_etaScaleFactor_.at(jetEtaBin) - 1.)*(1. - (associatedGenJet.pt()/jet.pt()));
    jet *= factor;
    return true;
}



int JetEnergyResolutionScaleFactors::jetEtaBin(const LV& jet)const
{
    int result = -1;
    
    const double eta = jet.eta();
    for(size_t iBin = 0; iBin < v_etaRange_.size(); ++iBin){
        if(std::fabs(eta) < v_etaRange_.at(iBin)){
            result = iBin;
            break;
        }
    }
    
    return result;
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



void JetEnergyScaleScaleFactors::applyJetSystematic(VLV* v_jet)const
{
    // This loop corrects the jet collection used for jet selections
    for(size_t iJet = 0; iJet < v_jet->size(); ++iJet){
        LV& jet = v_jet->at(iJet);
        
        this->scaleJet(jet);
    }
}



void JetEnergyScaleScaleFactors::applyMetSystematic(VLV* v_jetForMet, LV* met)const
{
    // This loop corrects the jet collection used to modify the MET
    double deltaPx = 0.;
    double deltaPy = 0.;
    for(size_t iJet = 0; iJet < v_jetForMet->size(); ++iJet){
        LV& jet = v_jetForMet->at(iJet);
        
        const double storedPx = jet.px();
        const double storedPy = jet.py();
        this->scaleJet(jet);
        deltaPx += jet.px() - storedPx;
        deltaPy += jet.py() - storedPy;
    }
    
    // Adjust the MET
    const double scaledMetPx = met->px() - deltaPx;
    const double scaledMetPy = met->py() - deltaPy;
    met->SetPt(std::sqrt(scaledMetPx*scaledMetPx + scaledMetPy*scaledMetPy));
}



void JetEnergyScaleScaleFactors::scaleJet(LV& jet)const
{
    jetCorrectionUncertainty_->setJetPt(jet.pt());
    jetCorrectionUncertainty_->setJetEta(jet.eta());
    const double uncertainty = jetCorrectionUncertainty_->getUncertainty(true);

    if(varyUp_) jet *= 1. + uncertainty;
    else jet *= 1. - uncertainty;
}






// --------------------------- Methods for TopPtScaleFactors ---------------------------------------------


TopPtScaleFactors::TopPtScaleFactors(const Systematic::Systematic& systematic):
reweightValue_(0)
{
    std::cout<<"--- Beginning preparation of top-pt scale factors\n";
    
    // Set up internal systematic
    SystematicInternal systematicInternal(nominal);
    if(systematic.type() == Systematic::topPt){
        if(systematic.variation() == Systematic::up){
            systematicInternal = vary_up;
            std::cout<<"Using systematic variation up\n";
        }
        else if(systematic.variation() == Systematic::down){
            systematicInternal = vary_down;
            std::cout<<"Using systematic variation down\n";
        }
    }
    
    // Set up reweighting function
    this->reweightFunction(systematicInternal);
    
    std::cout<<"=== Finishing preparation of top-pt scale factors\n\n";
}


double TopPtScaleFactors::getSF(const double& topPt, const double& antiTopPt)const
{
    return std::sqrt(reweightValue_(topPt) * reweightValue_(antiTopPt));
}



void TopPtScaleFactors::reweightFunction(const SystematicInternal& systematic)
{
    reweightValue_ = [systematic](const double& pt) -> double {
        
        // Fit parameters f0 and f1 of function exp(f0 - f1*pt)
        // only dilepton
        constexpr double p0 = 0.128;
        constexpr double p1 = 0.00121;
        // only semileptons
        //constexpr double p0 = 0.130;
        //constexpr double p1 = 0.00116;
        // dilepton + semileptons
        //constexpr double p0 = 0.130;
        //constexpr double p1 = 0.00118;
        
        if(systematic == nominal){
            return std::exp(p0-p1*pt);
        }
        else if(systematic == vary_up){
            // Systematic uncertainty such, that compatible with no correction (scale factor =1.)
            return 1.;
        }
        else if(systematic == vary_down){
            // Systematic uncertainty such, that twice the correction is assigned
            return 2.*std::exp(p0-p1*pt) - 1.;
        }
        else{
            std::cerr<<"ERROR in TopPtScaleFactors::reweightFunction()! Systematic is invalid\n...break\n"<<std::endl;
            exit(5);
        }
    };
}







