#include <iostream>
#include <fstream>
#include <cstdlib>
#include <vector>
#include <string>
#include <algorithm>

#include <TString.h>
#include <TFile.h>
#include <TTree.h>
#include <TProof.h>
#include <TSelector.h>
#include <TObjString.h>
#include <TChain.h>
#include <TH1.h>

#include "HiggsAnalysis.h"
#include "analysisHelpers.h"
#include "JetCategories.h"
#include "MvaTreeHandler.h"
#include "MvaReader.h"

#include "MvaValidation.h"
#include "DijetAnalyzer.h"
#include "AnalysisHistograms.h"
#include "BasicAnalyzer.h"
#include "JetMatchAnalyzer.h"
#include "JetChargeAnalyzer.h"
#include "Playground.h"
#include "../../common/include/sampleHelpers.h"
#include "../../common/include/utils.h"
#include "../../common/include/CommandLineParameters.h"
#include "../../common/include/KinematicReconstruction.h"
#include "../../common/include/ScaleFactors.h"
#include "../../common/include/BTagUtils.h"
#include "TopAnalysis/ZTopUtils/interface/PUReweighter.h"




/// Set pileup distribution file corresponding to data sample in use
/// The file ending is automatically adjusted for different systematics
//constexpr const char* PileupInputFILE = "Data_PUDist_19624pb.root";
//constexpr const char* PileupInputFILE = "Data_PUDist_19789pb.root";
constexpr const char* PileupInputFILE = "Data_PUDist_Full2012ReReco_FinalRecommendation.root";



/// Input file for electron ID scale factor
//constexpr const char* ElectronSFInputFILE = "ElectronSFtop12028.root";
//constexpr const char* ElectronSFInputFILE = "ElectronSFtop12028_19fb.root";
constexpr const char* ElectronSFInputFILE = "ElectronSF_198fbReReco.root";

/// Input file for muon ID scale factor
//constexpr const char* MuonSFInputFILE = "MuonSFtop12028.root";
//constexpr const char* MuonSFInputFILE = "MuonSFtop12028_19fb.root";
constexpr const char* MuonSFInputFILE = "MuonSF_198fbReReco.root";



/// File ending of dilepton trigger scale factors input file
//constexpr const char* TriggerSFInputSUFFIX = ".root";
//constexpr const char* TriggerSFInputSUFFIX = "_19fb.root";
constexpr const char* TriggerSFInputSUFFIX = "_rereco198fb.root";



/// File containing the uncertainties associated to JES
//constexpr const char* JesUncertaintySourceFILE = "Fall12_V7_DATA_UncertaintySources_AK5PFchs.txt";
// constexpr const char* JesUncertaintySourceFILE = "Summer13_V1_DATA_UncertaintySources_AK5PFchs.txt";
constexpr const char* JesUncertaintySourceFILE = "Summer13_V4_DATA_UncertaintySources_AK5PFchs.txt";




/// Folder where to find the b-/c-/l-tagging efficiencies
constexpr const char* BtagEfficiencyInputDIR = "BTagEff";

/// Folder for b-tag efficiency file storage (in case efficiencies are produced)
constexpr const char* BtagEfficiencyOutputDIR = "selectionRoot/BTagEff";



/// Folder for storage of MVA input TTree
constexpr const char* MvaInputDIR = "mvaInput";


/// Histogram containing the 2D distribution of MVA weights (needs to fit with the two weights also specified here)
constexpr const char* Mva2dWeightsFILE = "mvaOutput/Nominal/combined/weights/weights2d.root";


/// Folder for basic analysis output
constexpr const char* AnalysisOutputDIR = "selectionRoot";









void load_HiggsAnalysis(const TString& validFilenamePattern,
                        const Channel::Channel& channel,
                        const Systematic::Systematic& systematic,
                        const int dy,
                        const int jetCategoriesId,
                        const std::vector<AnalysisMode::AnalysisMode>& v_analysisMode)
{
    std::cout<<std::endl;

    // Set up the channels to run over
    std::vector<Channel::Channel> channels;
    if(channel != Channel::undefined){
        channels.push_back(channel);
    }
    else{
        channels = Channel::realChannels;
    }

    // Set up kinematic reconstruction
    KinematicReconstruction* kinematicReconstruction(0);
    kinematicReconstruction = new KinematicReconstruction();

    // Set up pileup reweighter
    std::cout<<"--- Beginning preparation of pileup reweighter\n";
    ztop::PUReweighter* puReweighter = new ztop::PUReweighter();
    puReweighter->setMCDistrSum12("S10");
    TString pileupInput(common::DATA_PATH_COMMON());
    pileupInput.Append("/").Append(PileupInputFILE);
    if(systematic == Systematic::pu_up) pileupInput.ReplaceAll(".root", "_sysUp.root");
    else if(systematic == Systematic::pu_down) pileupInput.ReplaceAll(".root", "_sysDown.root");
    std::cout<<"Using PU input file:\n"<<pileupInput<<std::endl;
    puReweighter->setDataTruePUInput(pileupInput.Data());
    std::cout<<"=== Finishing preparation of pileup reweighter\n\n";

    // Set up lepton efficiency scale factors
    LeptonScaleFactors::Systematic leptonSFSystematic(LeptonScaleFactors::nominal);
    if(systematic == Systematic::lept_up) leptonSFSystematic = LeptonScaleFactors::vary_up;
    else if(systematic == Systematic::lept_down) leptonSFSystematic = LeptonScaleFactors::vary_down;
    const LeptonScaleFactors leptonScaleFactors(ElectronSFInputFILE, MuonSFInputFILE, leptonSFSystematic);

    // Set up trigger efficiency scale factors (do it for all channels)
    TriggerScaleFactors::Systematic triggerSFSystematic(TriggerScaleFactors::nominal);
    if(systematic == Systematic::trig_up) triggerSFSystematic = TriggerScaleFactors::vary_up;
    else if(systematic == Systematic::trig_down) triggerSFSystematic = TriggerScaleFactors::vary_down;
    const TriggerScaleFactors triggerScaleFactors(TriggerSFInputSUFFIX,
                                                  Channel::convertChannels(Channel::realChannels),
                                                  triggerSFSystematic);
    
    BtagScaleFactors *btagScaleFactors = 0;
    BTagSFGeneric *bTagSFGeneric = 0;

    // Setting the flag to choose which kind of btag scale factors to use
    bool useGenericBTagSF = false;

    if(useGenericBTagSF){
        // Set up generic btag efficiency scale factors
        bTagSFGeneric = new BTagSFGeneric(BtagEfficiencyInputDIR,
                                          BtagEfficiencyOutputDIR,
                                          Channel::convertChannels(Channel::realChannels),
                                          Systematic::convertSystematic(systematic));
    }
    else{
        // Set up btag efficiency scale factors (do it for all channels)
        btagScaleFactors = new BtagScaleFactors(BtagEfficiencyInputDIR,
                                                BtagEfficiencyOutputDIR,
                                                Channel::convertChannels(Channel::realChannels),
                                                Systematic::convertSystematic(systematic));
    }
    
    // Set up JER systematic scale factors
    JetEnergyResolutionScaleFactors* jetEnergyResolutionScaleFactors(0);
    if(systematic==Systematic::jer_up || systematic==Systematic::jer_down){
        JetEnergyResolutionScaleFactors::Systematic jerSystematic(JetEnergyResolutionScaleFactors::vary_up);
        if(systematic == Systematic::jer_down) jerSystematic = JetEnergyResolutionScaleFactors::vary_down;
        jetEnergyResolutionScaleFactors = new JetEnergyResolutionScaleFactors(jerSystematic);
    }
    
    // Set up JES systematic scale factors
    JetEnergyScaleScaleFactors* jetEnergyScaleScaleFactors(0);
    if(systematic==Systematic::jes_up || systematic==Systematic::jes_down){
        JetEnergyScaleScaleFactors::Systematic jesSystematic(JetEnergyScaleScaleFactors::vary_up);
        if(systematic == Systematic::jes_down) jesSystematic = JetEnergyScaleScaleFactors::vary_down;
        jetEnergyScaleScaleFactors = new JetEnergyScaleScaleFactors(JesUncertaintySourceFILE, jesSystematic);
    }
    

    // Set up jet categories
    JetCategories* jetCategories(0);
    if(jetCategoriesId == 0) jetCategories = new JetCategories(2, 4, 1, 3, true, true); // Default categories
    else if(jetCategoriesId == 1) jetCategories = new JetCategories(2, 5, 0, 5, true, true); // Overview categories
    else if(jetCategoriesId == 2) jetCategories = new JetCategories(4, 4, 1, 3, true, true); // 4-jet categories of default categories
    else if(jetCategoriesId == 3) jetCategories = new JetCategories(4, 4, 2, 4, true, true); // Nazar's categories
    if(!jetCategories){
        std::cerr<<"Error in load_Analysis! No jet categories defined\n...break\n"<<std::endl;
        exit(832);
    }

    // Vector for setting up all analysers
    std::vector<AnalysisHistogramsBase*> v_analysisHistograms;

    // Set up event yield histograms
    EventYieldHistograms* eventYieldHistograms(0);
    eventYieldHistograms = new EventYieldHistograms({"1", "2", "3", "4", "5", "6", "7"}, {"7"}, jetCategories);
    v_analysisHistograms.push_back(eventYieldHistograms);

    // Set up Drell-Yan scaling histograms
    DyScalingHistograms* dyScalingHistograms(0);
    dyScalingHistograms = new DyScalingHistograms({"4", "5", "6", "7"}, "5");
    v_analysisHistograms.push_back(dyScalingHistograms);

    // Set up basic histograms
    BasicHistograms* basicHistograms(0);
    if(std::find(v_analysisMode.begin(), v_analysisMode.end(), AnalysisMode::cp) != v_analysisMode.end()){
        basicHistograms = new BasicHistograms({"1", "2", "3", "4", "5", "6", "7"}, {"7"}, jetCategories);
        v_analysisHistograms.push_back(basicHistograms);
    }

    // Set up jet charge analyzer
    JetChargeAnalyzer* jetChargeAnalyzer(0);
    if(std::find(v_analysisMode.begin(), v_analysisMode.end(), AnalysisMode::charge) != v_analysisMode.end()){
        jetChargeAnalyzer = new JetChargeAnalyzer({"7"}, {"7"}, jetCategories);
        v_analysisHistograms.push_back(jetChargeAnalyzer);
    }

    // Set up jet match analyzer
    JetMatchAnalyzer* jetMatchAnalyzer(0);
    if(std::find(v_analysisMode.begin(), v_analysisMode.end(), AnalysisMode::match) != v_analysisMode.end()){
        jetMatchAnalyzer = new JetMatchAnalyzer({"7"}, {"7"}, jetCategories);
        v_analysisHistograms.push_back(jetMatchAnalyzer);
    }

    // Set up playground
    Playground* playground(0);
    if(std::find(v_analysisMode.begin(), v_analysisMode.end(), AnalysisMode::playg) != v_analysisMode.end()){
        //playground = new Playground({"1", "2", "3", "4", "5", "6", "7"},{"6", "7"}, jetCategories);
        playground = new Playground({"1", "2", "3", "4", "5", "6", "7"});
        v_analysisHistograms.push_back(playground);
    }

    // Set up DijetAnalyzer
    DijetAnalyzer* dijetAnalyzer(0);
    if(std::find(v_analysisMode.begin(), v_analysisMode.end(), AnalysisMode::dijet) != v_analysisMode.end()){
        dijetAnalyzer = new DijetAnalyzer(Mva2dWeightsFILE, "", "", {}, {"7"}, jetCategories);
        v_analysisHistograms.push_back(dijetAnalyzer);
    }

    // Set up MVA validation, including reading in MVA weights in case they exist
    MvaValidation* mvaValidation(0);
    if(std::find(v_analysisMode.begin(), v_analysisMode.end(), AnalysisMode::mvaA) != v_analysisMode.end()){
        mvaValidation = new MvaValidation(Mva2dWeightsFILE, {"7"}, {"7"}, jetCategories);
        v_analysisHistograms.push_back(mvaValidation);
    }

    // Set up production of MVA input tree
    MvaTreeHandler* mvaTreeHandler(0);
    if(std::find(v_analysisMode.begin(), v_analysisMode.end(), AnalysisMode::mvaP) != v_analysisMode.end()){
        mvaTreeHandler = new MvaTreeHandler(MvaInputDIR, {"7"}, {"7"}, jetCategories);
    }

    // Set up the analysis
    HiggsAnalysis* selector = new HiggsAnalysis();
    selector->SetAnalysisOutputBase(AnalysisOutputDIR);
    selector->SetKinematicReconstruction(kinematicReconstruction);
    selector->SetPUReweighter(puReweighter);
    selector->SetLeptonScaleFactors(leptonScaleFactors);
    selector->SetTriggerScaleFactors(triggerScaleFactors);
    if(btagScaleFactors) selector->SetBtagScaleFactors(*btagScaleFactors);
    else if(bTagSFGeneric) selector->SetBtagScaleFactors(*bTagSFGeneric);
    selector->SetJetEnergyResolutionScaleFactors(jetEnergyResolutionScaleFactors);
    selector->SetJetEnergyScaleScaleFactors(jetEnergyScaleScaleFactors);
    
    selector->SetMvaInputProduction(mvaTreeHandler);
    selector->SetAllAnalysisHistograms(v_analysisHistograms);

    // Access selectionList containing all input sample nTuples
    ifstream infile("selectionList.txt");
    if(!infile.good()){
        std::cerr<<"Error! Please check the selectionList.txt file!\n"<<std::endl;
        exit(773);
    }

    // Loop over all input files
    int filecounter = 0;
    while(!infile.eof()){

        // Access nTuple input filename from selectionList
        TString filename;
        infile>>filename;
        if(filename=="" || filename[0]=='#') continue; //empty or commented line? --> skip
        if(!filename.Contains(validFilenamePattern)) continue;
        std::cout<<std::endl;
        std::cout<<"PROCESSING File "<<++filecounter<<" ("<<filename<<") from selectionList.txt"<<std::endl;
        std::cout<<std::endl;

        // Open nTuple file
        TFile file(filename);
        if(file.IsZombie()){
            std::cerr<<"ERROR! Cannot open nTuple file with name: "<<filename<<std::endl;
            exit(853);
        }

        // Check whether nTuple can be found
        TTree* tree = dynamic_cast<TTree*>(file.Get("writeNTuple/NTuple"));
        if(!tree){
            std::cerr<<"ERROR! TTree (=nTuple) not found in file!\n"<<std::endl;
            exit(854);
        }

        // Access information about samples, stored in the nTuples
        TObjString* channel_from_file = dynamic_cast<TObjString*>(file.Get("writeNTuple/channelName"));
        TObjString* systematics_from_file = dynamic_cast<TObjString*>(file.Get("writeNTuple/systematicsName"));
        TObjString* samplename = dynamic_cast<TObjString*>(file.Get("writeNTuple/sampleName"));
        TObjString* o_isSignal = dynamic_cast<TObjString*>(file.Get("writeNTuple/isSignal"));
        TObjString* o_isHiggsSignal = dynamic_cast<TObjString*>(file.Get("writeNTuple/isHiggsSignal"));
        TObjString* o_isMC = dynamic_cast<TObjString*>(file.Get("writeNTuple/isMC"));
        TH1* weightedEvents = dynamic_cast<TH1*>(file.Get("EventsBeforeSelection/weightedEvents"));
        if(!channel_from_file || !systematics_from_file || !o_isSignal || !o_isMC || !samplename){
            std::cout<<"Error: info about sample missing!"<<std::endl;
            return;
        }

        // Configure information about samples
        const bool isTopSignal = o_isSignal->GetString() == "1";
        const bool isMC = o_isMC->GetString() == "1";
        const bool isHiggsSignal(o_isHiggsSignal && o_isHiggsSignal->GetString()=="1");
        const bool isHiggsInclusive(isHiggsSignal && samplename->GetString()=="ttbarhiggsinclusive");
        const bool isTtbarV(samplename->GetString()=="ttbarw" || samplename->GetString()=="ttbarz");

        if(!isMC && systematic!=Systematic::undefined){
            std::cout<<"Sample is DATA, so not running again for systematic variation\n";
            continue;
        }

        // Is the channel given in the file?
        // FIXME: In case of data, use the channel as specified in the file, or what is this for?
        if(channel_from_file->GetString() != ""){
            channels.clear();
            channels.push_back(Channel::convertChannel(static_cast<std::string>(channel_from_file->GetString())));
        }

        // Loop over channels and run selector
        for(const auto& selectedChannel : channels){

            // Set output file name
            TString outputfilename(filename);
            if(outputfilename.Contains('/')){
                Ssiz_t last = outputfilename.Last('/');
                outputfilename = outputfilename.Data() + last + 1;
            }
            const TString channelName = Channel::convertChannel(selectedChannel);
            if(!outputfilename.BeginsWith(channelName + "_")) outputfilename.Prepend(channelName + "_");
            //outputfile is now channel_filename.root
            if(dy){
                if(outputfilename.First("_dy") == kNPOS){
                    std::cerr << "DY variations must be run on DY samples!\n";
                    std::cerr << outputfilename << " must start with 'channel_dy'\n";
                    exit(1);
                }
                const Channel::Channel dyChannel = dy == 11 ? Channel::ee : dy == 13 ? Channel::mumu : Channel::tautau;
                outputfilename.ReplaceAll("_dy", TString("_dy").Append(Channel::convertChannel(dyChannel)));
            }
            if(isHiggsInclusive)outputfilename.ReplaceAll("inclusive", "inclusiveOther");
            if(isTopSignal)outputfilename.ReplaceAll("signalplustau", "signalPlusOther");

            // Configure selector
            selector->SetChannel(channelName);
            selector->SetTopSignal(isTopSignal);
            selector->SetHiggsSignal(isHiggsSignal);
            selector->SetMC(isMC);
            selector->SetTrueLevelDYChannel(dy);
            if(systematic == Systematic::undefined){
                selector->SetSystematic(systematics_from_file->GetString());
            }
            else{
                selector->SetSystematic(Systematic::convertSystematic(systematic));
            }
            selector->SetWeightedEvents(weightedEvents);
            // FIXME: correction for MadGraph W decay branching fractions are not correctly applied
            // Recently it is done for W from ttbar decays, set via SetSamplename
            // Needs to be changed: for ttbarW, also correction for 3rd W needs to be applied, for ttbarhiggs corrections for 2 or 4 Ws needed, depending on Higgs decay (H->WW?)
            // and what about Wlnu sample, or possible others ?
            selector->SetSamplename(samplename->GetString(), systematics_from_file->GetString());
            selector->SetOutputfilename(outputfilename);
            selector->SetRunViaTau(0);
            selector->SetRunWithTtbb(0);
            selector->SetHiggsInclusiveSample(isHiggsInclusive);
            selector->SetHiggsInclusiveSeparation(false);

            // Set up nTuple chain and run selector
            TChain chain("writeNTuple/NTuple");
            chain.Add(filename);
            // chain.SetProof(); //will work from 5.34 onwards
            chain.Process(selector);

            // For splitting of dileptonic ttbar in tt+bb and tt+other
            if(isTopSignal && !isHiggsSignal && !isTtbarV){
                selector->SetRunWithTtbb(1);
                outputfilename.ReplaceAll("Other", "Bbbar");
                selector->SetOutputfilename(outputfilename);
                chain.Process(selector);
            }

            // For splitting of ttH inclusive decay in H->bb and other decays
            if(isHiggsInclusive){
                selector->SetHiggsInclusiveSeparation(true);
                outputfilename.ReplaceAll("Other", "Bbbar");
                selector->SetOutputfilename(outputfilename);
                chain.Process(selector);
            }
        }
        file.Close();
    }
}



int main(int argc, char** argv) {
    CLParameter<std::string> opt_filenamePattern("f", "Restrict to filename pattern, e.g. ttbar", false, 1, 1);
    CLParameter<std::string> opt_channel("c", "Specify a certain channel (ee, emu, mumu). No channel specified = run on all channels", false, 1, 1,
            common::makeStringCheck(Channel::convertChannels(Channel::allowedChannelsAnalysis)));
    CLParameter<std::string> opt_systematic("s", "Run with a systematic that runs on the nominal ntuples, e.g. 'PU_UP'", false, 1, 1,
            common::makeStringCheck(Systematic::convertSystematics(Systematic::allowedSystematicsHiggsAnalysis)));
    CLParameter<int> opt_dy("d", "Drell-Yan mode (11 for ee, 13 for mumu, 15 for tautau)", false, 1, 1,
            [](int dy){return dy == 11 || dy == 13 || dy == 15;});
    CLParameter<int> opt_jetCategoriesId("j", "ID for jet categories (# jets, # b-jets). If not specified, use default categories (=0)", false, 1, 1,
            [](int id){return id>=0 && id<=3;});
    CLParameter<std::string> opt_mode("m", "Mode of analysis: control plots (cp), Produce MVA input (mvaP), Apply MVA weights (mvaA), dijet analyser (dijet), playground (playg), charge analyser (charge), jet match analyser (match). Default is cp", false, 1, 100,
            common::makeStringCheck(AnalysisMode::convertAnalysisModes(AnalysisMode::allowedAnalysisModes)));
    CLAnalyser::interpretGlobal(argc, argv);

    // Set up a pattern for selecting only files from selectionList containing that pattern in filename
    const TString validFilenamePattern = opt_filenamePattern.isSet() ? opt_filenamePattern[0] : "";

    // Set up channel
    Channel::Channel channel(Channel::undefined);
    if(opt_channel.isSet()) channel = Channel::convertChannel(opt_channel[0]);

    // Set up systematic
    Systematic::Systematic systematic(Systematic::undefined);
    if(opt_systematic.isSet()) systematic = Systematic::convertSystematic(opt_systematic[0]);

    // Set up Drell-Yan mode
    const int dy = opt_dy.isSet() ? opt_dy[0] : 0;

    // Set up jet categories
    const int jetCategoriesId = opt_jetCategoriesId.isSet() ? opt_jetCategoriesId[0] : 0;

    // Set up analysis mode
    std::vector<AnalysisMode::AnalysisMode> v_analysisMode({AnalysisMode::cp});
    if(opt_mode.isSet()) v_analysisMode = AnalysisMode::convertAnalysisModes(opt_mode.getArguments());
    std::cout<<"\nRunning the following analysis modes:\n";
    for(const auto& analysisMode : v_analysisMode) std::cout<<AnalysisMode::convertAnalysisMode(analysisMode)<<" , ";
    std::cout<<"\n\n";

    // Start plotting
    //TProof* p = TProof::Open(""); // not before ROOT 5.34
    load_HiggsAnalysis(validFilenamePattern, channel, systematic, dy, jetCategoriesId, v_analysisMode);
    //delete p;
}



