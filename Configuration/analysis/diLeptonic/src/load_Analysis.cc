#include <iostream>
#include <fstream>
#include <cmath>
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

#include "TopAnalysis.h"
#include "analysisHelpers.h"
#include "AnalyzerBaseClass.h"
#include "AnalyzerControlPlots.h"
#include "AnalyzerDoubleDiffXS.h"
#include "AnalyzerKinReco.h"
#include "TreeHandlerBase.h"
#include "TreeHandlerTTBar.h"
#include "TreeHandlerBoostedTop.h"
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
//constexpr BtagScaleFactors::CorrectionMode BtagCorrectionMODE = BtagScaleFactors::greaterEqualOneTagReweight;
constexpr BtagScaleFactors::CorrectionMode BtagCorrectionMODE = BtagScaleFactors::randomNumberRetag;
//constexpr BtagScaleFactors::CorrectionMode BtagCorrectionMODE = BtagScaleFactors::discriminatorReweight;

/// File for the official heavy flavour scale factors for b-tag discriminator reweighting
constexpr const char* BtagHeavyFlavourFILE = "csv_rwt_hf_20pt_7_2_14.root";

/// File for the official light flavour scale factors for b-tag discriminator reweighting
constexpr const char* BtagLightFlavourFILE = "csv_rwt_lf_20pt_7_2_14.root";



/// File containing the fits for the MVA MET recoil corrections in data
constexpr const char* MvaMetRecoilDataFILE = "METrecoil_Fits_DataSummer2013.root";

/// File containing the fits for the MVA MET recoil corrections in MC
constexpr const char* MvaMetRecoilMcFILE = "METrecoil_Fits_MCSummer2013.root";








void load_Analysis(const TString& validFilenamePattern, 
                   const Channel::Channel& channel, 
                   const Systematic::Systematic& systematic,
                   const std::vector<AnalysisMode::AnalysisMode>& v_analysisMode,
                   const int specific_PDF,
                   const int dy,
                   const TString& closure,
                   const double& slope,
                   const Long64_t& maxEvents,
                   const Long64_t& skipEvents
                  )
{   
    std::cout<<std::endl;
    
    // Set up the channels to run over
    std::vector<Channel::Channel> channels;
    if(channel != Channel::undefined) channels.push_back(channel);
    else channels = Channel::realChannels;
    
    // Set up kinematic reconstruction
    const KinematicReconstruction* kinematicReconstruction(0);
    kinematicReconstruction = new KinematicReconstruction(1, true);
    
    // Set up kinematic reconstruction scale factors (null-pointer means no application)
    KinematicReconstructionScaleFactors* kinematicReconstructionScaleFactors(0);
    kinematicReconstructionScaleFactors = new KinematicReconstructionScaleFactors(channels, systematic);
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
    
    
    // Vector for setting up all analysers
    std::vector<AnalyzerBaseClass*> v_analyzer;
    
    // Set up event yield histograms
    AnalyzerEventYields* analyzerEventYields(0);
    analyzerEventYields = new AnalyzerEventYields({"1", "2", "3", "4", "5", "6", "7", "8"});
    v_analyzer.push_back(analyzerEventYields);
    
    // Set up Drell-Yan scaling histograms
    AnalyzerDyScaling* analyzerDyScaling(0);
    analyzerDyScaling = new AnalyzerDyScaling({"4", "5", "6", "7", "8"}, "5");
    v_analyzer.push_back(analyzerDyScaling);
    
    // Set up basic histograms
    AnalyzerControlPlots* analyzerControlPlots(0);
    if(std::find(v_analysisMode.begin(), v_analysisMode.end(), AnalysisMode::cp) != v_analysisMode.end()){
        analyzerControlPlots = new AnalyzerControlPlots({"1", "2", "3", "4", "5", "6", "7", "8"});
        v_analyzer.push_back(analyzerControlPlots);
    }
    // Set up dda histograms
    AnalyzerDoubleDiffXS* analyzerDoubleDiffXS(0);
    if(std::find(v_analysisMode.begin(), v_analysisMode.end(), AnalysisMode::dda) != v_analysisMode.end()){
        analyzerDoubleDiffXS = new AnalyzerDoubleDiffXS({"0","8"});
        v_analyzer.push_back(analyzerDoubleDiffXS);
    }
    
    // Set up KinReco histograms
    AnalyzerKinReco* analyzerKinReco(0);
    if(std::find(v_analysisMode.begin(), v_analysisMode.end(), AnalysisMode::kinReco) != v_analysisMode.end()){
        analyzerKinReco = new AnalyzerKinReco({"7","8"});
        v_analyzer.push_back(analyzerKinReco);
    }
    
    
    // Vector setting up all tree handlers for MVA input variables
    std::vector<TreeHandlerBase*> v_treeHandler;
    
    // Set up production of TTree for ttbar analysis
    TreeHandlerTTBar* treeHandlerTTBar(0);
    if(std::find(v_analysisMode.begin(), v_analysisMode.end(), AnalysisMode::ddaTree) != v_analysisMode.end()){
        treeHandlerTTBar = new TreeHandlerTTBar("ddaInput", {"0","8"});
        v_treeHandler.push_back(treeHandlerTTBar);
    }
    
    // Set up production of TTree for Boosted Top analysis
    TreeHandlerBoostedTop* treeHandlerBoostedTop(0);
    if(std::find(v_analysisMode.begin(), v_analysisMode.end(), AnalysisMode::btopTree) != v_analysisMode.end()){
        treeHandlerBoostedTop = new TreeHandlerBoostedTop("btopInput", {"0","8"});
        v_treeHandler.push_back(treeHandlerBoostedTop);
    }
    
    // Set up the analysis
    TopAnalysis* selector = new TopAnalysis();
    selector->SetAnalysisOutputBase("selectionRoot");
    selector->SetKinematicReconstruction(kinematicReconstruction, kinematicReconstructionScaleFactors);
    selector->SetPileupScaleFactors(pileupScaleFactors);
    selector->SetLeptonScaleFactors(leptonScaleFactors);
    selector->SetTriggerScaleFactors(triggerScaleFactors);
    selector->SetJetEnergyResolutionScaleFactors(jetEnergyResolutionScaleFactors);
    selector->SetJetEnergyScaleScaleFactors(jetEnergyScaleScaleFactors);
    selector->SetTopPtScaleFactors(topPtScaleFactors);
    selector->SetAllAnalyzers(v_analyzer);
    selector->SetAllTreeHandlers(v_treeHandler);
    
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
        
        // Checks avoiding running on ill-defined configurations
        if(!isMC && systematic.type()!=Systematic::undefinedType){
            std::cout<<"Sample is DATA, so not running again for systematic variation\n";
            continue;
        }
        if (systematic.type()==Systematic::pdf && (!isTopSignal || !(systematicFromFile.type()==Systematic::nominal))) {
            std::cout << "Skipping file: is not signal or not nominal -- and running PDFs\n";
            continue;
        }
        
        // Is the channel given in the file? This is true only for data which is preselected due to trigger,
        // and guarantees that only the proper channel is processed
        if(channelFromFile != Channel::undefined){
            channels.clear();
            channels.push_back(channelFromFile);
        }
        
        // If no systematic is specified, read it from the file and use this (used for systematic variations of signal samples, and for nominal)
        const Systematic::Systematic selectedSystematic = systematic.type()==Systematic::undefinedType ? systematicFromFile : systematic;
        
        // If for specific systematic variations the nominal btagging efficiencies should be used
        const Systematic::Systematic systematicForBtagEfficiencies = (selectedSystematic.type()==Systematic::pdf || selectedSystematic.type()==Systematic::closure) ? Systematic::nominalSystematic() : selectedSystematic;
        
        // Set up btag efficiency scale factors
        // This has to be done only after potentially setting systematic from file, since it is varied with signal systematics
        BtagScaleFactors btagScaleFactors("BTagEff", "selectionRoot/BTagEff", BtagHeavyFlavourFILE, BtagLightFlavourFILE,
                                          channels, systematicForBtagEfficiencies, BtagCorrectionMODE);
        
        // Configure selector
        selector->SetTopSignal(isTopSignal);
        selector->SetHiggsSignal(isHiggsSignal);
        selector->SetMC(isMC);
        selector->SetWeightedEvents(weightedEvents);
        selector->SetSamplename(samplename->GetString());
        selector->SetGeneratorBools(samplename->GetString(), selectedSystematic);
        selector->SetSystematic(selectedSystematic);
        selector->SetBtagScaleFactors(btagScaleFactors);
        
        // Loop over channels and run selector
        for(const auto& selectedChannel : channels){
            
            // Set the channel
            const TString channelName = Channel::convert(selectedChannel);
            TString outputfilename = filenameBase.BeginsWith(channelName+"_") ? filenameBase : channelName+"_"+filenameBase;
            triggerScaleFactors.prepareChannel(selectedChannel);
            btagScaleFactors.prepareChannel(selectedChannel);
            if(kinematicReconstructionScaleFactors) kinematicReconstructionScaleFactors->prepareChannel(selectedChannel);
            selector->SetChannel(selectedChannel);
            
            // Set up nTuple chain
            TChain chain("writeNTuple/NTuple");
            chain.Add(filename);
            // chain.SetProof(); //will work from 5.34 onwards
            
            // Split Drell-Yan sample in decay modes ee, mumu, tautau
            selector->SetTrueLevelDYChannel(dy);
            if(dy){
                if(outputfilename.First("_dy") == kNPOS){ 
                    std::cerr << "DY variations must be run on DY samples!\n";
                    std::cerr << outputfilename << " must start with 'channel_dy'\n";
                    exit(1);
                }
                
                if(selector->useMvaMet() && !metRecoilCorrector){
                    // Initialise recoil corrector for MVA MET in Drell-Yan samples (null-pointer means no application)
                    metRecoilCorrector = new MetRecoilCorrector(MvaMetRecoilDataFILE, MvaMetRecoilMcFILE);
                    selector->SetMetRecoilCorrector(metRecoilCorrector);
                }
                selector->SetDrellYan(true);
                const Channel::Channel dyChannel = dy == 11 ? Channel::ee : dy == 13 ? Channel::mumu : Channel::tautau;
                outputfilename.ReplaceAll("_dy", TString("_dy").Append(Channel::convert(dyChannel)));
            }
            else{
                selector->SetDrellYan(false);
            }
            
            // Run the selector
            selector->SetRunViaTau(0);
            selector->SetOutputfilename(outputfilename);
            selector->SetClosureTest(closure, slope);
            if(isTopSignal) selector->SetSampleForBtagEfficiencies(true);
            chain.Process(selector, "", maxEvents, skipEvents);
            selector->SetSampleForBtagEfficiencies(false);
            
            // For running on PDF systematics
            if(selectedSystematic.type() == Systematic::pdf){
                TH1* pdfWeights = dynamic_cast<TH1*>(file.Get("EventsBeforeSelection/pdfEventWeights"));
                if(!pdfWeights){
                    std::cerr << "Error: pdfEventWeights histo missing!\n";
                    exit(831);
                }
                for(int pdf_no = 0; pdfWeights->GetBinContent(pdf_no+1) > 0; ++pdf_no){
                    if(specific_PDF >=0 && pdf_no != specific_PDF) continue;
                    //weightedEvents->SetBinContent(1, pdfWeights->GetBinContent(pdf_no+1));
                    selector->SetWeightedEvents(weightedEvents);
                    selector->SetPDF(pdf_no);
                    chain.Process(selector, "", maxEvents, skipEvents);
                }
                continue;
            }
            
            // For splitting of dileptonic ttbar in component with intermediate taus and without
            if(isTopSignal && closure == ""){
                selector->SetRunViaTau(1);
                outputfilename.ReplaceAll("signalplustau", "bgviatau");
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
        btag, btagPt, btagEta,
        btagLjet, btagLjetPt, btagLjetEta,
        kin,
        //topPt,
        pdf, closure,
    };
}



int main(int argc, char** argv) {
    CLParameter<std::string> opt_f("f", "Restrict to filename pattern, e.g. ttbar", false, 1, 1);
    CLParameter<std::string> opt_c("c", "Specify a certain channel (ee, emu, mumu). No channel specified = run on all channels", false, 1, 1,
            common::makeStringCheck(Channel::convert(Channel::allowedChannelsAnalysis)));
    CLParameter<std::string> opt_s("s", "Run with a systematic that runs on the nominal ntuples, e.g. 'PDF', 'PU_UP' or 'TRIG_DOWN'", false, 1, 1,
            common::makeStringCheckBegin(Systematic::convertType(Systematic::allowedSystematics)));
    CLParameter<std::string> opt_mode("m", "Mode of analysis: control plots (cp), "
                                           "double differential analysis (dda), kin. reco. efficiency plots (kinReco), tree for 2d unfolding (ddaTree), "
                                           "Default is cp, dda, kinReco", false, 1, 100,
            common::makeStringCheck(AnalysisMode::convert(AnalysisMode::allowedAnalysisModes)));
    CLParameter<int> opt_pdfno("pdf", "Run a certain PDF systematic only, sets -s PDF. Use e.g. --pdf n, where n=0 is central, 1=variation 1 up, 2=1down, 3=2up, 4=2down, ...", false, 1, 1);
    CLParameter<int> opt_dy("d", "Drell-Yan mode (11 for ee, 13 for mumu, 15 for tautau)", false, 1, 1,
            [](int dy){return dy == 11 || dy == 13 || dy == 15;});
    CLParameter<std::string> opt_closure("closure", "Enable the closure test. Valid: pttop|ytop|nominal", false, 1, 1,
            [](const std::string &c){return c == "pttop" || c == "ytop" || c == "nominal";});
    CLParameter<double> opt_closureSlope("slope", "Slope for closure test, use -0.01 to 0.01 for pt and -0.4 to 0.4 for ytop", false, 1, 1,
            [](double s){return std::abs(s) < 1;});
    CLParameter<Long64_t> opt_maxEvents("maxEvents", "Maximum number of events to process", false, 1, 1,
            [](const Long64_t mE){return mE > 0;});
    CLParameter<Long64_t> opt_skipEvents("skipEvents", "Number of events to be skipped", false, 1, 1,
            [](const Long64_t sE){return sE > 0;});
    CLAnalyser::interpretGlobal(argc, argv);

    TString validFilenamePattern = opt_f.isSet() ? opt_f[0] : "";
    
    // Set up channel
    Channel::Channel channel(Channel::undefined);
    if(opt_c.isSet()) channel = Channel::convert(opt_c[0]);
    
    // Set up systematic
    Systematic::Systematic systematic(Systematic::undefinedSystematic());
    if(opt_s.isSet()) systematic = Systematic::Systematic(opt_s[0]);
    
    // Set up analysis mode
    std::vector<AnalysisMode::AnalysisMode> v_analysisMode({AnalysisMode::cp});
    
    if(opt_mode.isSet()) v_analysisMode = AnalysisMode::convert(opt_mode.getArguments());
    std::cout<<"\nRunning the following analysis modes:\n";
    for(const auto& analysisMode : v_analysisMode) std::cout<<AnalysisMode::convert(analysisMode)<<" , ";
    std::cout<<"\n\n";
    
    int dy = opt_dy.isSet() ? opt_dy[0] : 0;
    TString closure = opt_closure.isSet() ? opt_closure[0] : "";
    double slope = 0;
    if (closure != "" && closure != "nominal") {
        if (!opt_closureSlope.isSet()) {
            std::cerr << "closure test: need slope!\n"; exit(1);
        } else {
            slope = opt_closureSlope[0];
        }
    }
    int pdf_no = -1;
    if (opt_pdfno.isSet()) {
        pdf_no = opt_pdfno[0];
        if(pdf_no < 0){
            std::cerr<<"ERROR! PDF number is negative: "<<pdf_no<<"\n...break\n"<<std::endl;
            std::exit(1);
        }
        if (!(systematic.type() == Systematic::undefinedType || systematic.type() == Systematic::pdf)) {
            std::cerr << "Insonsistent systematic parameter: " << systematic.name() << " cannot be used with PDF systematic!\n";
            std::exit(1);
        }
        const Systematic::Variation variation = !pdf_no ? Systematic::central : pdf_no%2 ? Systematic::up : Systematic::down;
        const int variationNumber = (pdf_no+1)/2;
        systematic = Systematic::Systematic(Systematic::pdf, variation, variationNumber);
        std::cout << "Running PDF variation: " << systematic.name() << "\n";
    }
    
    // Set up maximum number of events to process
    const Long64_t bigNumber(TChain::kBigNumber);
    const Long64_t maxEvents = opt_maxEvents.isSet() ? opt_maxEvents[0] : bigNumber;

    // Set up number of events to be skipped
    const Long64_t skipEvents = opt_skipEvents.isSet() ? opt_skipEvents[0] : 0;
    
//     TProof* p = TProof::Open(""); // not before ROOT 5.34
    load_Analysis(validFilenamePattern, channel, systematic, v_analysisMode, pdf_no, dy, closure, slope, maxEvents, skipEvents);
//     delete p;
}
