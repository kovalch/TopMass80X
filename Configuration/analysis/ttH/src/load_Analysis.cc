#include <iostream>
#include <fstream>
#include <cstdlib>
#include <vector>
#include <set>
#include <string>
#include <sstream>

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
#include "AnalysisConfig.h"
#include "analysisHelpers.h"
#include "JetCategories.h"
#include "JetCharge.h"
#include "MvaTreeHandlerBase.h"
#include "MvaTreeHandlerTopJets.h"
#include "MvaTreeHandlerEventClassification.h"
#include "AnalyzerBase.h"
#include "AnalyzerMvaTopJets.h"
#include "AnalyzerMvaEventClassification.h"
#include "AnalyzerDijet.h"
#include "AnalyzerControlPlots.h"
#include "AnalyzerJetMatch.h"
#include "AnalyzerJetCharge.h"
#include "AnalyzerJetProperties.h"
#include "AnalyzerPlayground.h"
#include "AnalyzerEventWeight.h"
#include "AnalyzerGenEvent.h"
#include "AnalyzerKinematicReconstruction.h"
#include "../../common/include/sampleHelpers.h"
#include "../../common/include/utils.h"
#include "../../common/include/CommandLineParameters.h"
#include "../../common/include/KinematicReconstruction.h"
#include "../../common/include/ScaleFactors.h"
#include "../../common/include/Correctors.h"





void load_Analysis(const TString& validFilenamePattern,
                   const int part,
                   const Channel::Channel& channel,
                   const Systematic::Systematic& systematic,
                   const int systematicId,
                   const int jetCategoriesId,
                   const std::vector<AnalysisMode::AnalysisMode>& v_analysisMode,
                   const std::string& reweightingName,
                   const double& reweightingSlope,
                   const Long64_t& maxEvents,
                   const Long64_t& skipEvents)
{
    // Read analysis config from text file
    const AnalysisConfig analysisConfig;
    const AnalysisConfig::Corrections& corrections = analysisConfig.corrections();
    //analysisConfig.print();
    
    // Set up the channels to run over
    std::vector<Channel::Channel> channels;
    if(channel != Channel::undefined) channels.push_back(channel);
    else channels = Channel::realChannels;
    
    // Set up kinematic reconstruction
    const KinematicReconstruction* kinematicReconstruction(0);
    //kinematicReconstruction = new KinematicReconstruction(1, true);
    
    // Set up kinematic reconstruction scale factors (null-pointer means no application)
    KinematicReconstructionScaleFactors* kinematicReconstructionScaleFactors(0);
    //kinematicReconstructionScaleFactors = new KinematicReconstructionScaleFactors(channels, systematic);
    if(systematic.type()==Systematic::kin && !kinematicReconstructionScaleFactors){
        std::cout<<"Systematic for kinematic reconstruction requested, but scale factors not applied"
                 <<"\nStop running of analysis --- this is NOT an error, just avoiding double running\n"<<std::endl;
        exit(1);
    }
    
    // Set up pileup reweighter
    const PileupScaleFactors* const pileupScaleFactors = new PileupScaleFactors(corrections.pileupInputFile_, corrections.pileupMcEra_,
                                                                                corrections.pileupScenario_, systematic);
    
    // Set up lepton efficiency scale factors
    const LeptonScaleFactors leptonScaleFactors(corrections.electronSFInputFile_, corrections.muonSFInputFile_, systematic);
    
    // Set up trigger efficiency scale factors
    TriggerScaleFactors triggerScaleFactors(corrections.triggerSFInputSuffix_, channels, systematic);
    
    // Set up JER systematic scale factors (null-pointer means no application)
    const JetEnergyResolutionScaleFactors* jetEnergyResolutionScaleFactors(0);
    if(systematic.type() == Systematic::jer) jetEnergyResolutionScaleFactors = new JetEnergyResolutionScaleFactors(corrections.jerUncertaintySourceName_, systematic);
    
    // Set up JES systematic scale factors (null-pointer means no application)
    const JetEnergyScaleScaleFactors* jetEnergyScaleScaleFactors(0);
    if(systematic.type() == Systematic::jes) jetEnergyScaleScaleFactors = new JetEnergyScaleScaleFactors(corrections.jesUncertaintySourceFile_, systematic);
    
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
    
    // Set up jet charge
    JetCharge* jetCharge(0);
    jetCharge = new JetCharge(false, false);
    
    // Set up jet categories
    JetCategories* jetCategories(0);
    if(jetCategoriesId == 0) jetCategories = new JetCategories(2, 4, 1, 3, true, true); // Default categories
    else if(jetCategoriesId == 1) jetCategories = new JetCategories(2, 5, 0, 5, true, true); // Overview categories
    else if(jetCategoriesId == 2) jetCategories = new JetCategories(4, 4, 1, 3, true, true); // 4-jet categories of default categories
    else if(jetCategoriesId == 3) jetCategories = new JetCategories(4, 4, 2, 4, true, true); // Nazar's categories
    else if(jetCategoriesId == 4){ // For tt+bb and tt+b diff. xsection
        jetCategories = new JetCategories();
        jetCategories->addCategory(3, 3, true, false);
        jetCategories->addCategory(4, 4, true, false);
    }
    else if(jetCategoriesId == 5){ // For ttH baseline analysis
        jetCategories = new JetCategories();
        jetCategories->addCategory(3, 2, false, false);
        jetCategories->addCategory(4, 2, true, false);
        jetCategories->addCategory(3, 3, true, false);
        jetCategories->addCategory(4, 4, true, true);
    }
    if(!jetCategories){
        std::cerr<<"ERROR in load_Analysis()! No jet categories defined\n...break\n"<<std::endl;
        exit(832);
    }
    
    
    // Bools to define for which analysers and tree handlers genObjects are needed in first step
    bool genStudiesTtbb(false);
    bool genStudiesTth(false);
    
    
    // Vector for setting up all analysers
    std::vector<AnalyzerBase*> v_analyzer;
    
    // Set up event yield histograms
    AnalyzerEventYields* analyzerEventYields(0);
    analyzerEventYields = new AnalyzerEventYields({"0a", "0b", "0", "1a", "1", "2", "3", "4", "5", "6", "7", "4zWindow", "5zWindow", "6zWindow", "7zWindow"}, {"7"}, jetCategories);
    v_analyzer.push_back(analyzerEventYields);
    
    // Set up Drell-Yan scaling histograms
    AnalyzerDyScaling* analyzerDyScaling(0);
    analyzerDyScaling = new AnalyzerDyScaling({"4", "5", "6", "7"}, "5");
    v_analyzer.push_back(analyzerDyScaling);
    
    // Set up Heavy-Flavour fraction scaling histograms
    AnalyzerHfFracScaling* analyzerHfFracScaling(0);
    analyzerHfFracScaling = new AnalyzerHfFracScaling({"7"}, {"7"}, jetCategories);
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
    
    // Set up jet properties analyzer
    AnalyzerJetProperties* analyzerJetProperties(0);
    if(std::find(v_analysisMode.begin(), v_analysisMode.end(), AnalysisMode::jetProp) != v_analysisMode.end()){
        analyzerJetProperties = new AnalyzerJetProperties({"3", "4", "5", "6", "7"}, {"7"}, jetCategories);
        v_analyzer.push_back(analyzerJetProperties);
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
        analyzerDijet = new AnalyzerDijet("mvaOutput/Nominal/combined/weights/weights2d.root", "correct_step7_cate0_cate1_cate2_d144", "", {"0"}, {"7"}, jetCategories, false, true);
        v_analyzer.push_back(analyzerDijet);
        genStudiesTtbb = true;
    }
    
    // Set up event weight analyzer
    AnalyzerEventWeight* analyzerEventWeight(0);
    if(std::find(v_analysisMode.begin(), v_analysisMode.end(), AnalysisMode::weight) != v_analysisMode.end()){
        analyzerEventWeight = new AnalyzerEventWeight({"0", "1", "2", "3", "4", "5", "6", "7"}, {"7"}, jetCategories);
        v_analyzer.push_back(analyzerEventWeight);
    }
    
    // Set up gen event analyzer
    AnalyzerGenEvent* analyzerGenEvent(0);
    if(std::find(v_analysisMode.begin(), v_analysisMode.end(), AnalysisMode::genEvent) != v_analysisMode.end()){
        analyzerGenEvent = new AnalyzerGenEvent({"7"}, {"7"}, jetCategories);
        v_analyzer.push_back(analyzerGenEvent);
    }
    
    // Set up kinematic reconstruction analyzer
    AnalyzerKinematicReconstruction* analyzerKinematicReconstruction(0);
    if(std::find(v_analysisMode.begin(), v_analysisMode.end(), AnalysisMode::kinReco) != v_analysisMode.end()){
        analyzerKinematicReconstruction = new AnalyzerKinematicReconstruction({"7"}, {"7"}, jetCategories);
        v_analyzer.push_back(analyzerKinematicReconstruction);
    }
    
    // Set up MVA validation for top system jet assignment
    AnalyzerMvaTopJets* analyzerMvaTopJets(0);
    if(std::find(v_analysisMode.begin(), v_analysisMode.end(), AnalysisMode::mvaTopA) != v_analysisMode.end()){
        analyzerMvaTopJets = new AnalyzerMvaTopJets("mvaOutput/Nominal/combined/weights/weights2d.root", {"7"}, {"7"}, jetCategories);
        v_analyzer.push_back(analyzerMvaTopJets);
    }
    
    // Set up MVA validation for event classification
    AnalyzerMvaEventClassification* analyzerMvaEventClassification(0);
    if(std::find(v_analysisMode.begin(), v_analysisMode.end(), AnalysisMode::mvaEventA) != v_analysisMode.end()){
        analyzerMvaEventClassification = new AnalyzerMvaEventClassification({"7"}, {"7"}, jetCategories);
        v_analyzer.push_back(analyzerMvaEventClassification);
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
        mvaTreeHandlerEventClassification = new MvaTreeHandlerEventClassification("mvaInput_mvaEvent", {"7"}, {"7"}, jetCategories);
        v_mvaTreeHandler.push_back(mvaTreeHandlerEventClassification);
    }
    
    
    // Set up the analysis
    HiggsAnalysis* selector = new HiggsAnalysis(analysisConfig);
    selector->SetAnalysisOutputBase("selectionRoot");
    selector->SetKinematicReconstruction(kinematicReconstruction, kinematicReconstructionScaleFactors);
    selector->SetJetCharge(jetCharge);
    selector->SetPileupScaleFactors(pileupScaleFactors);
    selector->SetLeptonScaleFactors(leptonScaleFactors);
    selector->SetTriggerScaleFactors(triggerScaleFactors);
    selector->SetJetEnergyResolutionScaleFactors(jetEnergyResolutionScaleFactors);
    selector->SetJetEnergyScaleScaleFactors(jetEnergyScaleScaleFactors);
    selector->SetTopPtScaleFactors(topPtScaleFactors);
    selector->SetGenStudies(genStudiesTtbb, genStudiesTth);
    selector->SetAllAnalyzers(v_analyzer);
    selector->SetAllTreeHandlers(v_mvaTreeHandler);
    
    // Access selectionList containing all input sample nTuples
    ifstream infile("selectionList.txt");
    if(!infile.good()){
        std::cerr<<"ERROR in load_Analysis()! Cannot find selectionList.txt file\n...break\n"<<std::endl;
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
        
        // In case of reweighting, modify basic filename to contain reweighting name and slope
        if(reweightingName != ""){
            genStudiesTtbb = true;
            std::stringstream sstream;
            if(reweightingName == "nominal") sstream<<"_reweighted_"<<reweightingName<<".root";
            else sstream<<"_reweighted_"<<reweightingName<<"_"<<reweightingSlope<<".root";
            filenameBase.ReplaceAll(".root", sstream.str());
        }
        
        // For ttbar samples, adjust basic filename to represent physics process (signal --> Dilepton, bg --> Notdilepton)
        // In case of systematic ttbar sample, remove from filename part identifying the systematic (will be identified by output folder)
        if(filenameBase.BeginsWith("ttbarsignal")){
            if(filenameBase.BeginsWith("ttbarsignalplustau")) filenameBase = "ttbarsignalplustau.root";
            filenameBase.ReplaceAll("signal", "Dilepton");
        }
        else if(filenameBase.BeginsWith("ttbarbg")){
            if(filenameBase.BeginsWith("ttbarbg")) filenameBase = "ttbarbg.root";
            filenameBase.ReplaceAll("bg", "Notdilepton");
        }
        
        // Open nTuple file
        TFile file(filename);
        if(file.IsZombie()){
            std::cerr<<"ERROR in load_Analysis()! Cannot open nTuple file with name: "<<filename<<"\n...break\n"<<std::endl;
            exit(853);
        }
        
        // Check whether nTuple can be found
        TTree* tree = dynamic_cast<TTree*>(file.Get("writeNTuple/NTuple"));
        if(!tree){
            std::cerr<<"ERROR in load_Analysis()! TTree (=nTuple) not found in file\n...break\n"<<std::endl;
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
            std::cerr<<"ERROR in load_Analysis()! Info about sample missing\n...break\n"<<std::endl; 
            exit(855);
        }
        const Channel::Channel channelFromFile = Channel::convert(channel_from_file->GetString());
        const Systematic::Systematic systematicFromFile = Systematic::Systematic(systematics_from_file->GetString());
        
        // Configure information about samples
        const bool isTopSignal = o_isSignal->GetString() == "1";
        const bool isMC = o_isMC->GetString() == "1";
        const bool isHiggsSignal(o_isHiggsSignal && o_isHiggsSignal->GetString()=="1");
        const bool isHiggsInclusive(isHiggsSignal && samplename->GetString()=="ttbarhiggsinclusive");
        const bool isTtbarNotdilepton(samplename->GetString() == "ttbarbg");
        const bool isTtbarV(samplename->GetString()=="ttbarw" || samplename->GetString()=="ttbarz");
        const bool isDrellYan(samplename->GetString()=="dy1050" || samplename->GetString()=="dy50inf");
        
        // Checks avoiding running on ill-defined configurations
        if(!isMC && systematic.type()!=Systematic::undefinedType){
            std::cout<<"Sample is DATA, so not running systematic variation\n";
            continue;
        }
        if(systematic.type() == Systematic::pdf){
            if(!isTopSignal || isHiggsSignal || isTtbarV || systematicFromFile.type() != Systematic::nominal){
                std::cout<<"Sample is not ttbar dilepton or not nominal, so not running PDF variation\n";
                continue;
            }
            TH1* pdfWeights = dynamic_cast<TH1*>(file.Get("EventsBeforeSelection/pdfEventWeights"));
            if(!pdfWeights){
                std::cerr<<"ERROR in load_Analysis()! Cannot find histogram pdfEventWeights\n...break\n"<<std::endl;
                exit(831);
            }
            if(systematicId>=pdfWeights->GetNbinsX() || pdfWeights->GetBinContent(systematicId+1)<=0.){
                std::cout<<"ID specified for PDF variation is above number of variations stored in sample, so not running\n";
                continue;
            }
        }
        if(reweightingName != ""){
            if(!isTopSignal || isHiggsSignal || isTtbarV || systematicFromFile.type() != Systematic::nominal){
                std::cout<<"Sample is not ttbar dilepton or not nominal, so not running reweighting\n";
                continue;
            }
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
        const Systematic::Systematic systematicForBtagEfficiencies = (selectedSystematic.type()==Systematic::pdf || reweightingName!="") ? Systematic::nominalSystematic() : selectedSystematic;
        
        // Set up btag efficiency scale factors
        // This has to be done only after potentially setting systematic from file, since it is varied with signal systematics
        BtagScaleFactors btagScaleFactors("BTagEff", "selectionRoot/BTagEff",
                                          corrections.btagHeavyFlavourFile_, corrections.btagLightFlavourFile_,
                                          channels, systematicForBtagEfficiencies, corrections.btagCorrectionMode_);
        
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
        selector->SetReweightingName(reweightingName);
        selector->SetReweightingSlope(reweightingSlope);
        
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
                if(analysisConfig.selections().mvaMet_ && !metRecoilCorrector){
                    // Initialise recoil corrector for MVA MET in Drell-Yan samples (null-pointer means no application)
                    metRecoilCorrector = new MetRecoilCorrector(corrections.mvaMetRecoilDataFile_, corrections.mvaMetRecoilMcFile_);
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
                    std::cerr<<"ERROR in load_Analysis()! Specified part for Drell-Yan separation is not allowed (sample, part): "
                             <<outputfilename<<" , "<<part<<"\n...break\n"<<std::endl;
                    exit(647);
                }
                // Reset the selection
                selector->SetTrueLevelDYChannel(0);
            }
            else if(isTopSignal && !isHiggsSignal && !isTtbarV){ // For splitting of ttbar dilepton in production modes associated with heavy flavours, and decays via taus
                if(selectedSystematic.type() == Systematic::pdf) selector->SetPdfVariation(systematicId);
                const std::set<int> allowedPartIds({0, 101, 201, 102, 202, 103, 203, 4, -1});
                if(part==0 || part==-1){ // output is tt+other: both leptons from W->e/mu or W->tau->e/mu
                    selector->SetAdditionalBjetMode(0);
                    TString modifiedOutputfilename(outputfilename);
                    modifiedOutputfilename.ReplaceAll("plustau", "PlustauOther");
                    selector->SetOutputfilename(modifiedOutputfilename);
                    selector->SetSampleForBtagEfficiencies(true);
                    chain.Process(selector, "", maxEvents, skipEvents);
                    selector->SetSampleForBtagEfficiencies(false);
                }
                if(part==101 || part==-1){ // output is tt+b: both leptons from W->e/mu
                    selector->SetAdditionalBjetMode(101);
                    TString modifiedOutputfilename(outputfilename);
                    modifiedOutputfilename.ReplaceAll("plustau", "NotauB");
                    selector->SetOutputfilename(modifiedOutputfilename);
                    chain.Process(selector, "", maxEvents, skipEvents);
                }
                if(part==201 || part==-1){ // output is tt+b: 1+ lepton from W->tau->e/mu
                    selector->SetAdditionalBjetMode(201);
                    TString modifiedOutputfilename(outputfilename);
                    modifiedOutputfilename.ReplaceAll("plustau", "OnlytauB");
                    selector->SetOutputfilename(modifiedOutputfilename);
                    chain.Process(selector, "", maxEvents, skipEvents);
                }
                if(part==102 || part==-1){ // output is tt+2b with 2+ b-hadrons/jet: both leptons from W->e/mu
                    selector->SetAdditionalBjetMode(102);
                    TString modifiedOutputfilename(outputfilename);
                    modifiedOutputfilename.ReplaceAll("plustau", "Notau2b");
                    selector->SetOutputfilename(modifiedOutputfilename);
                    chain.Process(selector, "", maxEvents, skipEvents);
                }
                if(part==202 || part==-1){ // output is tt+2b with 2+ b-hadrons/jet: 1+ lepton from W->tau->e/mu
                    selector->SetAdditionalBjetMode(202);
                    TString modifiedOutputfilename(outputfilename);
                    modifiedOutputfilename.ReplaceAll("plustau", "Onlytau2b");
                    selector->SetOutputfilename(modifiedOutputfilename);
                    chain.Process(selector, "", maxEvents, skipEvents);
                }
                if(part==103 || part==-1){ // output is tt+bb: both leptons from W->e/mu
                    selector->SetAdditionalBjetMode(103);
                    TString modifiedOutputfilename(outputfilename);
                    modifiedOutputfilename.ReplaceAll("plustau", "NotauBbbar");
                    selector->SetOutputfilename(modifiedOutputfilename);
                    chain.Process(selector, "", maxEvents, skipEvents);
                }
                if(part==203 || part==-1){ // output is tt+bb: 1+ lepton from W->tau->e/mu
                    selector->SetAdditionalBjetMode(203);
                    TString modifiedOutputfilename(outputfilename);
                    modifiedOutputfilename.ReplaceAll("plustau", "OnlytauBbbar");
                    selector->SetOutputfilename(modifiedOutputfilename);
                    chain.Process(selector, "", maxEvents, skipEvents);
                }
                if(part==4 || part==-1){ // output is tt+cc: both leptons from W->e/mu or W->tau->e/mu
                    selector->SetAdditionalBjetMode(4);
                    TString modifiedOutputfilename(outputfilename);
                    modifiedOutputfilename.ReplaceAll("plustau", "PlustauCcbar");
                    selector->SetOutputfilename(modifiedOutputfilename);
                    chain.Process(selector, "", maxEvents, skipEvents);
                }
                if(allowedPartIds.count(part) < 1){
                    std::cerr<<"ERROR in load_Analysis()! Specified part for ttbar+HF separation is not allowed (sample, part): "
                             <<outputfilename<<" , "<<part<<"\n...break\n"<<std::endl;
                    std::cerr << "Allowed parts: ";
                    for(const int partId : allowedPartIds) std::cerr << " " << partId;
                    std::cerr << std::endl;
                    exit(647);
                }
                // Reset the selection
                selector->SetPdfVariation(-1);
                selector->SetAdditionalBjetMode(-999);
            }
            else if(isTtbarNotdilepton && analysisConfig.general().era_!=Era::run1_8tev){ // For splitting of ttbar non-dilepton in production modes associated with heavy flavours
                // FIXME: Check for era is workaround as long as 8tev ntuples do not have tt+xx ID
                if(part==0 || part==-1){ // output is tt+other
                    selector->SetAdditionalBjetMode(0);
                    TString modifiedOutputfilename(outputfilename);
                    modifiedOutputfilename.ReplaceAll("Notdilepton", "NotdileptonOther");
                    selector->SetOutputfilename(modifiedOutputfilename);
                    selector->SetSampleForBtagEfficiencies(true);
                    chain.Process(selector, "", maxEvents, skipEvents);
                    selector->SetSampleForBtagEfficiencies(false);
                }
                if(part==1 || part==-1){ // output is tt+b
                    selector->SetAdditionalBjetMode(1);
                    TString modifiedOutputfilename(outputfilename);
                    modifiedOutputfilename.ReplaceAll("Notdilepton", "NotdileptonB");
                    selector->SetOutputfilename(modifiedOutputfilename);
                    chain.Process(selector, "", maxEvents, skipEvents);
                }
                if(part==2 || part==-1){ // output is tt+2b with 2+ b-hadrons/jet
                    selector->SetAdditionalBjetMode(2);
                    TString modifiedOutputfilename(outputfilename);
                    modifiedOutputfilename.ReplaceAll("Notdilepton", "Notdilepton2b");
                    selector->SetOutputfilename(modifiedOutputfilename);
                    chain.Process(selector, "", maxEvents, skipEvents);
                }
                if(part==3 || part==-1){ // output is tt+bb
                    selector->SetAdditionalBjetMode(3);
                    TString modifiedOutputfilename(outputfilename);
                    modifiedOutputfilename.ReplaceAll("Notdilepton", "NotdileptonBbbar");
                    selector->SetOutputfilename(modifiedOutputfilename);
                    chain.Process(selector, "", maxEvents, skipEvents);
                }
                if(part==4 || part==-1){ // output is tt+cc
                    selector->SetAdditionalBjetMode(4);
                    TString modifiedOutputfilename(outputfilename);
                    modifiedOutputfilename.ReplaceAll("Notdilepton", "NotdileptonCcbar");
                    selector->SetOutputfilename(modifiedOutputfilename);
                    chain.Process(selector, "", maxEvents, skipEvents);
                }
                if(part >= 5){
                    std::cerr<<"ERROR in load_Analysis()! Specified part for ttbar non-dilepton separation is not allowed (sample, part): "
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
                    std::cerr<<"ERROR in load_Analysis()! Specified part for Higgs inclusive separation is not allowed (sample, part): "
                             <<outputfilename<<" , "<<part<<"\n...break\n"<<std::endl;
                    exit(647);
                }
                // Reset the selection
                selector->SetInclusiveHiggsDecayMode(-999);
            }
            else{ // All other samples which are not split in subsamples
                if(part >= 0){
                    std::cerr<<"ERROR in load_Analysis()! Specified part for sample separation is not allowed, this sample cannot be split (sample, part): "
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
        btagDiscrBpurity, btagDiscrLpurity,
        btagDiscrCerr1, btagDiscrCerr2,
        kin,
        topPt,
        pdf
    };
}



int main(int argc, char** argv)
{
    CLParameter<std::string> opt_filenamePattern("f", "Restrict to filename pattern, e.g. ttbar", false, 1, 1);
    CLParameter<int> opt_partToRun("p", "Specify a part to be run for samples which are split in subsamples (via ID). Default is no splitting", false, 1, 1,
            [](int part){return part%100>=0 && part%100<5;});
    CLParameter<std::string> opt_channel("c", "Specify a certain channel (ee, emu, mumu). No channel specified = run on all channels", false, 1, 1,
            common::makeStringCheck(Channel::convert(Channel::allowedChannelsAnalysis)));
    CLParameter<std::string> opt_systematic("s", "Run with a systematic that runs on the nominal ntuples, e.g. 'PU_UP'", false, 1, 1,
            common::makeStringCheckBegin(Systematic::convertType(Systematic::allowedSystematics)));
    CLParameter<int> opt_systematicId("sid", "ID of systematic variation for systematics requiring one, e.g. 'PDF'", false, 1, 1,
            [](int id){return id>=0;});
    CLParameter<int> opt_jetCategoriesId("j", "ID for jet categories (# jets, # b-jets). If not specified, use default categories (=0)", false, 1, 1,
            [](int id){return id>=0 && id<=5;});
    CLParameter<std::string> opt_mode("m", "Mode of analysis: control plots (cp), "
                                           "dijet analyser (dijet), jet charge analyser (charge), jet match analyser (match), jet proerties analyser (jetProp), "
                                           "playground (playg), "
                                           "event weight analyser (weight), gen event analyser(genEvent), "
                                           "kinematic reconstruction analyser(kinReco), "
                                           "Produce MVA input or Apply MVA weights for top jets (mvaTopP/mvaTopA), "
                                           "Produce MVA input or Apply MVA weights for event classification (mvaEventP/mvaEventA), "
                                           "Produce MVA input or Apply MVA weights for jet charge (mvaChargeP/mvaChargeA). "
                                           "Default is cp", false, 1, 100,
            common::makeStringCheck(AnalysisMode::convert(AnalysisMode::allowedAnalysisModes)));
    CLParameter<std::string> opt_reweightName("reweightName", "Name of the event reweighting to be applied (nominal, 1st_add_bjet_pt, 1st_add_bjet_eta, 2nd_add_bjet_pt, 2nd_add_bjet_eta, add_bjet_dR, add_bjet_Mjj). Default: none", false, 1, 1);
    CLParameter<double> opt_reweightSlope("reweightSlope", "Slope of the event reweighting to be applied. Default: 0.0 Has no effect with nominal reweighting", false, 1, 1);
    CLParameter<Long64_t> opt_maxEvents("maxEvents", "Maximum number of events to process", false, 1, 1,
            [](const Long64_t mE){return mE > 0;});
    CLParameter<Long64_t> opt_skipEvents("skipEvents", "Number of events to be skipped", false, 1, 1,
            [](const Long64_t sE){return sE > 0;});
    CLAnalyser::interpretGlobal(argc, argv);
    
    std::cout<<"\n"<<"--- Beginning setting up command line steering parameters\n";
    
    // Set up a pattern for selecting only files from selectionList containing that pattern in filename
    TString validFilenamePattern("");
    if(opt_filenamePattern.isSet()){
        validFilenamePattern = opt_filenamePattern[0];
        std::cout<<"Using file pattern: "<<validFilenamePattern<<"\n";
    }
    
    // Set up part to be run for splitted samples
    int part(-1);
    if(opt_partToRun.isSet()){
        part = opt_partToRun[0];
        std::cout<<"Run part of sample with ID: "<<part<<"\n";
    }
    
    // Set up channel
    Channel::Channel channel(Channel::undefined);
    if(opt_channel.isSet()){
        channel = Channel::convert(opt_channel[0]);
        std::cout<<"Decay channel: "<<Channel::convert(channel)<<"\n";
    }
    
    // Set up systematic
    Systematic::Systematic systematic(Systematic::undefinedSystematic());
    if(opt_systematic.isSet()){
        systematic = Systematic::Systematic(opt_systematic[0]);
        std::cout<<"Systematic variation: "<<systematic.name()<<"\n";
    }
    
    // Set up systematic ID for systematics requiring one (e.g. PDF)
    int systematicId(-1);
    if(opt_systematicId.isSet()){
        systematicId = opt_systematicId[0];
        if(systematic.type() == Systematic::pdf){
            const Systematic::Variation variation = !systematicId ? Systematic::central : systematicId%2 ? Systematic::up : Systematic::down;
            const int variationNumber = (systematicId+1)/2;
            systematic = Systematic::Systematic(Systematic::pdf, variation, variationNumber);
        }
        else{
            std::cerr<<"ERROR in load_Analysis executable! Systematic is not of type requiring systematic ID, but one specified: "
                     <<systematicId<<"\n...break\n"<<std::endl;
            exit(1);
        }
        std::cout<<"Systematic variation constructed from ID (ID, name): "<<systematicId<<" , "<<systematic.name()<<"\n";
    }
    else{
        if(systematic.type() == Systematic::pdf){
            std::cerr<<"ERROR in load_Analysis executable! Systematic requires systematic ID, but none specified\n...break\n"<<std::endl;
            exit(1);
        }
    }
    
    // Set up jet categories
    const int jetCategoriesId = opt_jetCategoriesId.isSet() ? opt_jetCategoriesId[0] : 0;
    std::cout<<"Running on jet categories with ID: "<<jetCategoriesId<<"\n";
    
    // Set up analysis mode
    std::vector<AnalysisMode::AnalysisMode> v_analysisMode({AnalysisMode::cp});
    if(opt_mode.isSet()) v_analysisMode = AnalysisMode::convert(opt_mode.getArguments());
    std::cout<<"Running the following analysis modes:\n";
    for(const auto& analysisMode : v_analysisMode) std::cout<<AnalysisMode::convert(analysisMode)<<" , ";
    std::cout<<"\n";
    
    // Set up name and slope of event reweighting
    if((opt_reweightName.isSet() && !opt_reweightSlope.isSet()) || (!opt_reweightName.isSet() && opt_reweightSlope.isSet())){
        
    }
    std::string reweightingName(""); 
    if(opt_reweightName.isSet()){
        reweightingName = opt_reweightName[0];
        if(!opt_reweightSlope.isSet()){
            std::cerr<<"ERROR in load_Analysis executable! Reweighting name specified, but no slope: "
                     <<reweightingName<<"\n...break\n"<<std::endl;
            exit(1);
        }
    }
    double reweightingSlope(0.);
    if(opt_reweightSlope.isSet()){
        reweightingSlope = opt_reweightSlope[0];
        if(!opt_reweightName.isSet()){
            std::cerr<<"ERROR in load_Analysis executable! Reweighting slope specified, but no name: "
                     <<reweightingSlope<<"\n...break\n"<<std::endl;
            exit(1);
        }
        std::cout<<"Apply reweighting (name, slope): "<<reweightingName<<" , "<<reweightingSlope<<"\n";
    }
    
    
    // Set up maximum number of events to process
    const Long64_t bigNumber(TChain::kBigNumber);
    Long64_t maxEvents(bigNumber);
    if(opt_maxEvents.isSet()){
        maxEvents = opt_maxEvents[0];
        std::cout<<"Number of events to process: "<<maxEvents<<"\n";
    }
    
    // Set up number of events to be skipped
    Long64_t skipEvents(0);
    if(opt_skipEvents.isSet()){
        skipEvents = opt_skipEvents[0];
        std::cout<<"Number of events to skip: "<<skipEvents<<"\n";
    }
    
    std::cout<<"=== Finishing setting up command line steering parameters\n\n";
    
    // Start analysis
    //TProof* p = TProof::Open(""); // not before ROOT 5.34
    load_Analysis(validFilenamePattern, part, channel, systematic, systematicId, jetCategoriesId,
                  v_analysisMode, reweightingName, reweightingSlope, maxEvents, skipEvents);
    //delete p;
}



