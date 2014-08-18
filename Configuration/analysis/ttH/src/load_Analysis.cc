#include <iostream>
#include <fstream>
#include <cstdlib>
#include <vector>
#include <string>

#include <TString.h>
#include <TFile.h>
#include <TTree.h>
#include <TProof.h>
#include <TSelector.h>
#include <TObjString.h>
#include <TChain.h>
#include <TH1.h>
#include <Rtypes.h>

#include "HiggsAnalysis.h"
#include "analysisHelpers.h"
#include "JetCategories.h"
#include "MvaTreeHandlerBase.h"
#include "MvaTreeHandlerTopJets.h"
#include "MvaTreeHandlerEventClassification.h"
#include "AnalyzerBase.h"
#include "AnalyzerMvaTopJets.h"
#include "AnalyzerDijet.h"
#include "AnalyzerControlPlots.h"
#include "AnalyzerJetMatch.h"
#include "AnalyzerJetCharge.h"
#include "AnalyzerPlayground.h"
#include "AnalyzerEventWeight.h"
#include "AnalyzerGenEvent.h"
#include "../../common/include/sampleHelpers.h"
#include "../../common/include/utils.h"
#include "../../common/include/CommandLineParameters.h"
#include "../../common/include/KinematicReconstruction.h"
#include "../../common/include/ScaleFactors.h"
#include "../../common/include/Correctors.h"





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



/// Source for the uncertainties associated to JER
constexpr const char* JerUncertaintySourceNAME = "jer2011";
//constexpr const char* JerUncertaintySourceNAME = "jer2012";



/// File containing the uncertainties associated to JES
//constexpr const char* JesUncertaintySourceFILE = "Fall12_V7_DATA_UncertaintySources_AK5PFchs.txt";
//constexpr const char* JesUncertaintySourceFILE = "Summer13_V1_DATA_UncertaintySources_AK5PFchs.txt";
constexpr const char* JesUncertaintySourceFILE = "Summer13_V4_DATA_UncertaintySources_AK5PFchs.txt";
//constexpr const char* JesUncertaintySourceFILE = "Summer13_V5_DATA_UncertaintySources_AK5PFchs.txt";



/// The correction mode for the b-tagging
//constexpr BtagScaleFactors::CorrectionMode BtagCorrectionMODE = BtagScaleFactors::none;
//constexpr BtagScaleFactors::CorrectionMode BtagCorrectionMODE = BtagScaleFactors::randomNumberRetag;
constexpr BtagScaleFactors::CorrectionMode BtagCorrectionMODE = BtagScaleFactors::discriminatorReweight;

/// File for the official heavy flavour scale factors for b-tag discriminator reweighting
constexpr const char* BtagHeavyFlavourFILE = "csv_rwt_hf_20pt_7_2_14.root";

/// File for the official light flavour scale factors for b-tag discriminator reweighting
constexpr const char* BtagLightFlavourFILE = "csv_rwt_lf_20pt_7_2_14.root";



/// File containing the fits for the MVA MET recoil corrections in data
constexpr const char* MvaMetRecoilDataFILE = "METrecoil_Fits_DataSummer2013.root";

/// File containing the fits for the MVA MET recoil corrections in MC
constexpr const char* MvaMetRecoilMcFILE = "METrecoil_Fits_MCSummer2013.root";



/// Histogram containing the 2D distributions of MVA weights
constexpr const char* Mva2dWeightsFILE = "mvaOutput/Nominal/combined/weights/weights2d.root";









void load_Analysis(const TString& validFilenamePattern,
                   const int part,
                   const Channel::Channel& channel,
                   const Systematic::Systematic& systematic,
                   const int jetCategoriesId,
                   const std::vector<AnalysisMode::AnalysisMode>& v_analysisMode,
                   const Long64_t& maxEvents,
                   const Long64_t& skipEvents)
{
    std::cout<<std::endl;
    
    // Set up the channels to run over
    std::vector<Channel::Channel> channels;
    if(channel != Channel::undefined) channels.push_back(channel);
    else channels = Channel::realChannels;
    
    // Set up kinematic reconstruction
    KinematicReconstruction* kinematicReconstruction(0);
    //kinematicReconstruction = new KinematicReconstruction();
    
    // Set up kinematic reconstruction scale factors (null-pointer means no application)
    KinematicReconstructionScaleFactors* kinematicReconstructionScaleFactors(0);
    //kinematicReconstructionScaleFactors = new KinematicReconstructionScaleFactors(channels, systematic);
    if(systematic.type()==Systematic::kin && !kinematicReconstructionScaleFactors){
        std::cout<<"Systematic for kinematic reconstruction requested, but scale factors not applied"
                 <<"\nStop running of analysis --- this is NOT an error, just avoiding double running\n"<<std::endl;
        exit(1);
    }
    
    // Set up pileup reweighter
    const PileupScaleFactors* const pileupScaleFactors = new PileupScaleFactors(PileupInputFILE, "Summer12", "S10", systematic);
    
    // Set up lepton efficiency scale factors
    const LeptonScaleFactors leptonScaleFactors(ElectronSFInputFILE, MuonSFInputFILE, systematic);
    
    // Set up trigger efficiency scale factors
    TriggerScaleFactors triggerScaleFactors(TriggerSFInputSUFFIX, channels, systematic);
    
    // Set up JER systematic scale factors (null-pointer means no application)
    const JetEnergyResolutionScaleFactors* jetEnergyResolutionScaleFactors(0);
    if(systematic.type() == Systematic::jer) jetEnergyResolutionScaleFactors = new JetEnergyResolutionScaleFactors(JerUncertaintySourceNAME, systematic);
    
    // Set up JES systematic scale factors (null-pointer means no application)
    const JetEnergyScaleScaleFactors* jetEnergyScaleScaleFactors(0);
    if(systematic.type() == Systematic::jes) jetEnergyScaleScaleFactors = new JetEnergyScaleScaleFactors(JesUncertaintySourceFILE, systematic);
    
    // Set up top-pt reweighting scale factors (null-pointer means no application)
    const TopPtScaleFactors* topPtScaleFactors(0);
    //topPtScaleFactors = new TopPtScaleFactors(systematic);
    if(systematic.type()==Systematic::topPt && !topPtScaleFactors){
        std::cout<<"Systematic for top-pt reweighting requested, but scale factors not applied"
                 <<"\nStop running of analysis --- this is NOT an error, just avoiding double running\n"<<std::endl;
        exit(1);
    }
    
    // Set up recoil corrector for MVA MET in Drell-Yan samples (initialised only in first Drell-Yan sample)
    const MetRecoilCorrector* metRecoilCorrector(0);
    
    
    // Set up jet categories
    JetCategories* jetCategories(0);
    if(jetCategoriesId == 0) jetCategories = new JetCategories(2, 4, 1, 3, true, true); // Default categories
    else if(jetCategoriesId == 1) jetCategories = new JetCategories(2, 5, 0, 5, true, true); // Overview categories
    else if(jetCategoriesId == 2) jetCategories = new JetCategories(4, 4, 1, 3, true, true); // 4-jet categories of default categories
    else if(jetCategoriesId == 3) jetCategories = new JetCategories(4, 4, 2, 4, true, true); // Nazar's categories
    else if(jetCategoriesId == 4) jetCategories = new JetCategories(4, 4, 4, 4, true, false); // 4 exclusive b-tags for tt+bb diff. xsection (dijet)
    else if(jetCategoriesId == 5) jetCategories = new JetCategories(3, 3, 3, 3, true, true);  // 3 inclusive b-tags for tt+bb diff. xsection (jet)
    if(!jetCategories){
        std::cerr<<"Error in load_Analysis! No jet categories defined\n...break\n"<<std::endl;
        exit(832);
    }
    
    
    // Vector for setting up all analysers
    std::vector<AnalyzerBase*> v_analyzer;
    
    // Set up event yield histograms
    AnalyzerEventYields* analyzerEventYields(0);
    analyzerEventYields = new AnalyzerEventYields({"0a", "0b", "1", "2", "3", "4", "5", "6", "7"}, {"7"}, jetCategories);
    v_analyzer.push_back(analyzerEventYields);
    
    // Set up Drell-Yan scaling histograms
    AnalyzerDyScaling* analyzerDyScaling(0);
    analyzerDyScaling = new AnalyzerDyScaling({"4", "5", "6", "7"}, "5");
    v_analyzer.push_back(analyzerDyScaling);
    
    // Set up Heavy-Flavour fraction scaling histograms
    AnalyzerHfFracScaling* analyzerHfFracScaling(0);
    analyzerHfFracScaling = new AnalyzerHfFracScaling({"5", "6", "7"});
    v_analyzer.push_back(analyzerHfFracScaling);
    
    // Set up basic histograms
    AnalyzerControlPlots* analyzerControlPlots(0);
    if(std::find(v_analysisMode.begin(), v_analysisMode.end(), AnalysisMode::cp) != v_analysisMode.end()){
        analyzerControlPlots = new AnalyzerControlPlots({"1", "2", "3", "4", "5", "6", "7"}, {"7"}, jetCategories);
        v_analyzer.push_back(analyzerControlPlots);
    }
    
    // Set up jet charge analyzer
    AnalyzerJetCharge* analyzerJetCharge(0);
    if(std::find(v_analysisMode.begin(), v_analysisMode.end(), AnalysisMode::charge) != v_analysisMode.end()){
        analyzerJetCharge = new AnalyzerJetCharge({"7"});
        v_analyzer.push_back(analyzerJetCharge);
    }
    
    // Set up jet match analyzer
    AnalyzerJetMatch* analyzerJetMatch(0);
    if(std::find(v_analysisMode.begin(), v_analysisMode.end(), AnalysisMode::match) != v_analysisMode.end()){
        analyzerJetMatch = new AnalyzerJetMatch({"7"}, {"7"}, jetCategories);
        v_analyzer.push_back(analyzerJetMatch);
    }
    
    // Set up playground
    AnalyzerPlayground* anayzerPlayground(0);
    if(std::find(v_analysisMode.begin(), v_analysisMode.end(), AnalysisMode::playg) != v_analysisMode.end()){
        //analyzerPlayground = new AnalyzerPlayground({"1", "2", "3", "4", "5", "6", "7"},{"6", "7"}, jetCategories);
        anayzerPlayground = new AnalyzerPlayground({"1", "2", "3", "4", "5", "6", "7"});
        v_analyzer.push_back(anayzerPlayground);
    }
    
    // Set up DijetAnalyzer
    AnalyzerDijet* analyzerDijet(0);
    if(std::find(v_analysisMode.begin(), v_analysisMode.end(), AnalysisMode::dijet) != v_analysisMode.end()){
        analyzerDijet = new AnalyzerDijet(Mva2dWeightsFILE, "correct_step7_cate0_cate1_cate2_c1", "", {}, {"7"}, jetCategories, false, true);
        v_analyzer.push_back(analyzerDijet);
    }
    
    // Set up event weight analyzer
    AnalyzerEventWeight* analyzerEventWeight(0);
    if(std::find(v_analysisMode.begin(), v_analysisMode.end(), AnalysisMode::weight) != v_analysisMode.end()){
        analyzerEventWeight = new AnalyzerEventWeight({"0b", "1", "2", "3", "4", "5", "6", "7"}, {"7"}, jetCategories);
        v_analyzer.push_back(analyzerEventWeight);
    }
    
    // Set up gen event analyzer
    AnalyzerGenEvent* analyzerGenEvent(0);
    if(std::find(v_analysisMode.begin(), v_analysisMode.end(), AnalysisMode::genEvent) != v_analysisMode.end()){
        analyzerGenEvent = new AnalyzerGenEvent({"7"}, {"7"}, jetCategories);
        v_analyzer.push_back(analyzerGenEvent);
    }
    
    // Set up MVA validation for top system jet assignment
    AnalyzerMvaTopJets* analyzerMvaTopJets(0);
    if(std::find(v_analysisMode.begin(), v_analysisMode.end(), AnalysisMode::mvaTopA) != v_analysisMode.end()){
        analyzerMvaTopJets = new AnalyzerMvaTopJets(Mva2dWeightsFILE, {"7"}, {"7"}, jetCategories);
        v_analyzer.push_back(analyzerMvaTopJets);
    }
    
    
    // Vector setting up all tree handlers for MVA input variables
    std::vector<MvaTreeHandlerBase*> v_mvaTreeHandler;
    
    // Set up production of MVA input tree for top system jet assignment
    MvaTreeHandlerTopJets* mvaTreeHandlerTopJets(0);
    if(std::find(v_analysisMode.begin(), v_analysisMode.end(), AnalysisMode::mvaTopP) != v_analysisMode.end()){
        mvaTreeHandlerTopJets = new MvaTreeHandlerTopJets("mvaInput", {"7"}, {"7"}, jetCategories);
        v_mvaTreeHandler.push_back(mvaTreeHandlerTopJets);
    }
    
    // Set up production of MVA input tree for event classification
    MvaTreeHandlerEventClassification* mvaTreeHandlerEventClassification(0);
    if(std::find(v_analysisMode.begin(), v_analysisMode.end(), AnalysisMode::mvaEventP) != v_analysisMode.end()){
        mvaTreeHandlerEventClassification = new MvaTreeHandlerEventClassification("mvaInput", {"7"}, {"7"}, jetCategories);
        v_mvaTreeHandler.push_back(mvaTreeHandlerEventClassification);
    }
    
    
    // Set up the analysis
    HiggsAnalysis* selector = new HiggsAnalysis();
    selector->SetAnalysisOutputBase("selectionRoot");
    selector->SetKinematicReconstruction(kinematicReconstruction, kinematicReconstructionScaleFactors);
    selector->SetPileupScaleFactors(pileupScaleFactors);
    selector->SetLeptonScaleFactors(leptonScaleFactors);
    selector->SetTriggerScaleFactors(triggerScaleFactors);
    selector->SetJetEnergyResolutionScaleFactors(jetEnergyResolutionScaleFactors);
    selector->SetJetEnergyScaleScaleFactors(jetEnergyScaleScaleFactors);
    selector->SetTopPtScaleFactors(topPtScaleFactors);
    selector->SetAllAnalyzers(v_analyzer);
    selector->SetAllTreeHandlers(v_mvaTreeHandler);
    
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
        std::cout<<"\nPROCESSING File "<<++filecounter<<" ("<<filename<<") from selectionList.txt\n"<<std::endl;
        
        // Access the basic filename by stripping off the folders
        TString filenameBase(filename);
        if(filenameBase.Contains('/')){
            Ssiz_t last = filenameBase.Last('/');
            filenameBase = filenameBase.Data() + last + 1;
        }
        
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
            std::cerr<<"Error: info about sample missing!"<<std::endl; 
            exit(855);
        }
        const Channel::Channel channelFromFile = Channel::convert(channel_from_file->GetString());
        const Systematic::Systematic systematicFromFile = Systematic::Systematic(systematics_from_file->GetString());
        
        // Configure information about samples
        const bool isTopSignal = o_isSignal->GetString() == "1";
        const bool isMC = o_isMC->GetString() == "1";
        const bool isHiggsSignal(o_isHiggsSignal && o_isHiggsSignal->GetString()=="1");
        const bool isHiggsInclusive(isHiggsSignal && samplename->GetString()=="ttbarhiggsinclusive");
        const bool isTtbarV(samplename->GetString()=="ttbarw" || samplename->GetString()=="ttbarz");
        const bool isDrellYan(samplename->GetString()=="dy1050" || samplename->GetString()=="dy50inf");
        
        // Checks avoiding running on ill-defined configurations
        if(!isMC && systematic.type()!=Systematic::undefinedType){
            std::cout<<"Sample is DATA, so not running again for systematic variation\n";
            continue;
        }
        
        // Is the channel given in the file? This is true only for data which is preselected due to trigger,
        // and guarantees that only the proper channel is processed
        if(channelFromFile != Channel::undefined){
            channels.clear();
            channels.push_back(channelFromFile);
        }
        
        // If no systematic is specified, read it from the file and use this (used for systematic variations of signal samples)
        const Systematic::Systematic selectedSystematic = systematic.type()==Systematic::undefinedType ? systematicFromFile : systematic;
        
        // If for specific systematic variations the nominal btagging efficiencies should be used
        const Systematic::Systematic systematicForBtagEfficiencies = selectedSystematic;
        
        // Set up btag efficiency scale factors
        // This has to be done only after potentially setting systematic from file, since it is varied with signal systematics
        BtagScaleFactors btagScaleFactors("BTagEff", "selectionRoot/BTagEff", BtagHeavyFlavourFILE, BtagLightFlavourFILE,
                                          channels, systematicForBtagEfficiencies, BtagCorrectionMODE);
        
        // Configure selector
        selector->SetTopSignal(isTopSignal);
        selector->SetHiggsSignal(isHiggsSignal);
        selector->SetMC(isMC);
        selector->SetDrellYan(isDrellYan);
        selector->SetWeightedEvents(weightedEvents);
        // FIXME: correction for MadGraph W decay branching fractions are not correctly applied
        // Recently it is done for W from ttbar decays, set via SetGeneratorBools
        // Needs to be changed: for ttbarW, also correction for 3rd W needs to be applied, for ttbarhiggs corrections for 2 or 4 Ws needed, depending on Higgs decay (H->WW?)
        // and what about Wlnu sample, or possible others ?
        selector->SetSamplename(samplename->GetString());
        selector->SetGeneratorBools(samplename->GetString(), selectedSystematic);
        selector->SetSystematic(selectedSystematic);
        selector->SetBtagScaleFactors(btagScaleFactors);
        
        // Loop over channels and run selector
        for(const auto& selectedChannel : channels){
            
            // Set the channel
            const TString channelName = Channel::convert(selectedChannel);
            const TString outputfilename = filenameBase.BeginsWith(channelName+"_") ? filenameBase : channelName+"_"+filenameBase;
            triggerScaleFactors.prepareChannel(selectedChannel);
            btagScaleFactors.prepareChannel(selectedChannel);
            if(kinematicReconstructionScaleFactors) kinematicReconstructionScaleFactors->prepareChannel(selectedChannel);
            selector->SetChannel(selectedChannel);
            
            // Set up nTuple chain
            TChain chain("writeNTuple/NTuple");
            chain.Add(filename);
            // chain.SetProof(); //will work from 5.34 onwards
            
            // Split specific samples into subsamples and run the selector
            if(isDrellYan){ // For splitting of Drell-Yan sample in decay modes ee, mumu, tautau
                if(selector->useMvaMet() && !metRecoilCorrector){
                    metRecoilCorrector = new MetRecoilCorrector(MvaMetRecoilDataFILE, MvaMetRecoilMcFILE);
                    selector->SetMetRecoilCorrector(metRecoilCorrector);
                }
                if(part==0 || part==-1){ // output is DY->ee
                    selector->SetTrueLevelDYChannel(11);
                    const Channel::Channel dyChannel = Channel::ee;
                    TString modifiedOutputfilename(outputfilename);
                    modifiedOutputfilename.ReplaceAll("_dy", TString("_dy").Append(Channel::convert(dyChannel)));
                    selector->SetOutputfilename(modifiedOutputfilename);
                    chain.Process(selector, "", maxEvents, skipEvents);
                }
                if(part==1 || part==-1){ // output is DY->mumu
                    selector->SetTrueLevelDYChannel(13);
                    const Channel::Channel dyChannel = Channel::mumu;
                    TString modifiedOutputfilename(outputfilename);
                    modifiedOutputfilename.ReplaceAll("_dy", TString("_dy").Append(Channel::convert(dyChannel)));
                    selector->SetOutputfilename(modifiedOutputfilename);
                    chain.Process(selector, "", maxEvents, skipEvents);
                }
                if(part==2 || part==-1){ // output is DY->tautau
                    selector->SetTrueLevelDYChannel(15);
                    const Channel::Channel dyChannel = Channel::tautau;
                    TString modifiedOutputfilename(outputfilename);
                    modifiedOutputfilename.ReplaceAll("_dy", TString("_dy").Append(Channel::convert(dyChannel)));
                    selector->SetOutputfilename(modifiedOutputfilename);
                    chain.Process(selector, "", maxEvents, skipEvents);
                }
                if(part >= 3){
                    std::cerr<<"ERROR in load_Analysis! Specified part for Drell-Yan separation is not allowed (sample, part): "
                             <<outputfilename<<" , "<<part<<"\n...break\n"<<std::endl;
                    exit(647);
                }
                // Reset the selection
                selector->SetTrueLevelDYChannel(0);
            }
            else if(isTopSignal && !isHiggsSignal && !isTtbarV){ // For splitting of ttbar in production modes associated with heavy flavours
                //selector->SetRunViaTau(0); // This could be used for splitting of dileptonic ttbar in component with intermediate taus and without
                if(part==0 || part==-1){ // output is ttbar+other
                    selector->SetAdditionalBjetMode(0);
                    TString modifiedOutputfilename(outputfilename);
                    modifiedOutputfilename.ReplaceAll("signalplustau", "signalPlusOther");
                    selector->SetOutputfilename(modifiedOutputfilename);
                    selector->SetSampleForBtagEfficiencies(true);
                    chain.Process(selector, "", maxEvents, skipEvents);
                    selector->SetSampleForBtagEfficiencies(false);
                }
                if(part==1 || part==-1){ // output is ttbar+b with >=2 b-hadrons in a single jet
                    selector->SetAdditionalBjetMode(1);
                    TString modifiedOutputfilename(outputfilename);
                    modifiedOutputfilename.ReplaceAll("signalplustau", "signalPlus2B");
                    selector->SetOutputfilename(modifiedOutputfilename);
                    chain.Process(selector, "", maxEvents, skipEvents);
                }
                if(part==2 || part==-1){ // output is ttbar+b
                    selector->SetAdditionalBjetMode(2);
                    TString modifiedOutputfilename(outputfilename);
                    modifiedOutputfilename.ReplaceAll("signalplustau", "signalPlusB");
                    selector->SetOutputfilename(modifiedOutputfilename);
                    chain.Process(selector, "", maxEvents, skipEvents);
                }
                if(part==3 || part==-1){ // output is ttbar+bbbar
                    selector->SetAdditionalBjetMode(3);
                    TString modifiedOutputfilename(outputfilename);
                    modifiedOutputfilename.ReplaceAll("signalplustau", "signalPlusBbbar");
                    selector->SetOutputfilename(modifiedOutputfilename);
                    chain.Process(selector, "", maxEvents, skipEvents);
                }
                if(part==4 || part==-1){ // output is ttbar+ccbar
                    selector->SetAdditionalBjetMode(4);
                    TString modifiedOutputfilename(outputfilename);
                    modifiedOutputfilename.ReplaceAll("signalplustau", "signalPlusCcbar");
                    selector->SetOutputfilename(modifiedOutputfilename);
                    chain.Process(selector, "", maxEvents, skipEvents);
                }
                if(part >= 5){
                    std::cerr<<"ERROR in load_Analysis! Specified part for ttbar+HF separation is not allowed (sample, part): "
                             <<outputfilename<<" , "<<part<<"\n...break\n"<<std::endl;
                    exit(647);
                }
                // Reset the selection
                selector->SetAdditionalBjetMode(-999);
            }
            else if(isHiggsInclusive){ // For splitting of ttH inclusive decay in H->bb and other decays
                if(part==0 || part==-1){ // output is H->other
                    selector->SetInclusiveHiggsDecayMode(0);
                    TString modifiedOutputfilename(outputfilename);
                    modifiedOutputfilename.ReplaceAll("inclusive", "inclusiveOther");
                    selector->SetOutputfilename(modifiedOutputfilename);
                    chain.Process(selector, "", maxEvents, skipEvents);
                }
                if(part==1 || part==-1){ // output is H->bb
                    selector->SetInclusiveHiggsDecayMode(5);
                    TString modifiedOutputfilename(outputfilename);
                    modifiedOutputfilename.ReplaceAll("inclusive", "inclusiveBbbar");
                    selector->SetOutputfilename(modifiedOutputfilename);
                    chain.Process(selector, "", maxEvents, skipEvents);
                }
                if(part >= 2){
                    std::cerr<<"ERROR in load_Analysis! Specified part for Higgs inclusive separation is not allowed (sample, part): "
                             <<outputfilename<<" , "<<part<<"\n...break\n"<<std::endl;
                    exit(647);
                }
                // Reset the selection
                selector->SetInclusiveHiggsDecayMode(-999);
            }
            else{ // All other samples which are not split in subsamples
                if(part >= 0){
                    std::cerr<<"ERROR in load_Analysis! Specified part for sample separation is not allowed, this sample cannot be split (sample, part): "
                             <<outputfilename<<" , "<<part<<"\n...break\n"<<std::endl;
                    exit(647);
                }
                selector->SetOutputfilename(outputfilename);
                chain.Process(selector, "", maxEvents, skipEvents);
            }
        }
        file.Close();
    }
}



/// All systematics allowed for analysis step
/// Only systematics which run on the nominal ntuples, e.g. pileup reweighting
namespace Systematic{
    const std::vector<Type> allowedSystematics = {
        nominal,
        pu, lept, trig,
        jer, jes,
        btag, btagLjet,
        btagDiscrBstat1, btagDiscrBstat2,
        btagDiscrLstat1, btagDiscrLstat2,
        btagDiscrCerr1, btagDiscrCerr2,
        kin,
        topPt,
    };
}



int main(int argc, char** argv)
{
    CLParameter<std::string> opt_filenamePattern("f", "Restrict to filename pattern, e.g. ttbar", false, 1, 1);
    CLParameter<int> opt_partToRun("p", "Specify a part to be run for samples which are split in subsamples (via ID). Default is no splitting", false, 1, 1,
            [](int part){return part>=0 && part<5;});
    CLParameter<std::string> opt_channel("c", "Specify a certain channel (ee, emu, mumu). No channel specified = run on all channels", false, 1, 1,
            common::makeStringCheck(Channel::convert(Channel::allowedChannelsAnalysis)));
    CLParameter<std::string> opt_systematic("s", "Run with a systematic that runs on the nominal ntuples, e.g. 'PU_UP'", false, 1, 1,
            common::makeStringCheckBegin(Systematic::convertType(Systematic::allowedSystematics)));
    CLParameter<int> opt_jetCategoriesId("j", "ID for jet categories (# jets, # b-jets). If not specified, use default categories (=0)", false, 1, 1,
            [](int id){return id>=0 && id<=5;});
    CLParameter<std::string> opt_mode("m", "Mode of analysis: control plots (cp), "
                                           "dijet analyser (dijet), jet charge analyser (charge), jet match analyser (match), playground (playg), "
                                           "event weight analyser (weight), gen event analyser(genEvent), "
                                           "Produce MVA input or Apply MVA weights for top jets (mvaTopP/mvaTopA), "
                                           "Produce MVA input or Apply MVA weights for event classification (mvaEventP/mvaEventA), "
                                           "Produce MVA input or Apply MVA weights for jet charge (mvaChargeP/mvaChargeA). "
                                           "Default is cp", false, 1, 100,
            common::makeStringCheck(AnalysisMode::convert(AnalysisMode::allowedAnalysisModes)));
    CLParameter<Long64_t> opt_maxEvents("maxEvents", "Maximum number of events to process", false, 1, 1,
            [](const Long64_t mE){return mE > 0;});
    CLParameter<Long64_t> opt_skipEvents("skipEvents", "Number of events to be skipped", false, 1, 1,
            [](const Long64_t sE){return sE > 0;});
    CLAnalyser::interpretGlobal(argc, argv);
    
    // Set up a pattern for selecting only files from selectionList containing that pattern in filename
    const TString validFilenamePattern = opt_filenamePattern.isSet() ? opt_filenamePattern[0] : "";
    
    // Set up part to be run for splitted samples
    const int part = opt_partToRun.isSet() ? opt_partToRun[0] : -1;
    
    // Set up channel
    Channel::Channel channel(Channel::undefined);
    if(opt_channel.isSet()) channel = Channel::convert(opt_channel[0]);
    
    // Set up systematic
    Systematic::Systematic systematic(Systematic::undefinedSystematic());
    if(opt_systematic.isSet()) systematic = Systematic::Systematic(opt_systematic[0]);
    
    // Set up jet categories
    const int jetCategoriesId = opt_jetCategoriesId.isSet() ? opt_jetCategoriesId[0] : 0;
    
    // Set up analysis mode
    std::vector<AnalysisMode::AnalysisMode> v_analysisMode({AnalysisMode::cp});
    if(opt_mode.isSet()) v_analysisMode = AnalysisMode::convert(opt_mode.getArguments());
    std::cout<<"\nRunning the following analysis modes:\n";
    for(const auto& analysisMode : v_analysisMode) std::cout<<AnalysisMode::convert(analysisMode)<<" , ";
    std::cout<<"\n\n";
    
    // Set up maximum number of events to process
    const Long64_t bigNumber(TChain::kBigNumber);
    const Long64_t maxEvents = opt_maxEvents.isSet() ? opt_maxEvents[0] : bigNumber;
    
    // Set up number of events to be skipped
    const Long64_t skipEvents = opt_skipEvents.isSet() ? opt_skipEvents[0] : 0;
    
    // Start plotting
    //TProof* p = TProof::Open(""); // not before ROOT 5.34
    load_Analysis(validFilenamePattern, part, channel, systematic, jetCategoriesId, v_analysisMode, maxEvents, skipEvents);
    //delete p;
}



