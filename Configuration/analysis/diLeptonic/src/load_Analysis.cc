#include <iostream>
#include <fstream>
#include <cmath>

#include <TString.h>
#include <TFile.h>
#include <TTree.h>
#include <TProof.h>
#include <TSelector.h>
#include <TObjString.h>
#include <TChain.h>
#include <TH1.h>

#include "TopAnalysis.h"
#include "../../common/include/CommandLineParameters.h"
#include "../../common/include/utils.h"
#include "../../common/include/ScaleFactors.h"
#include "../../common/include/sampleHelpers.h"
#include "../../common/include/KinematicReconstruction.h"
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
constexpr const char* JesUncertaintySourceFILE = "Summer13_V4_DATA_UncertaintySources_AK5PFchs.txt";



/// Folder where to find the b-/c-/l-tagging efficiencies
constexpr const char* BtagEfficiencyInputDIR = "BTagEff";

/// Folder for b-tag efficiency file storage (in case efficiencies are produced)
constexpr const char* BtagEfficiencyOutputDIR = "selectionRoot/BTagEff";



/// Folder for basic analysis output
constexpr const char* AnalysisOutputDIR = "selectionRoot";








const TString pdfDirName(int pdf_no) {
    TString result("PDF_");
    if (pdf_no == 0) result += "CENTRAL";
    else result += TString::Format("%d", (pdf_no+1)/2) + "_" + (pdf_no % 2 ? "UP" : "DOWN");
    return result;
}



void load_Analysis(TString validFilenamePattern, 
                   TString channel, 
                   TString systematic,
                   int specific_PDF,
                   int dy,
                   TString closure,
                   double slope
                  )
{   
    // Set up the channels to run over
    std::vector<std::string> channels;
    if(channel != ""){
        channels.push_back(static_cast<std::string>(channel));
    }
    else{
        channels = {"ee", "emu", "mumu"};
    }
    
    // Set up kinematic reconstruction
    std::cout<<std::endl;
    KinematicReconstruction* kinematicReconstruction(0);
    kinematicReconstruction = new KinematicReconstruction();
    
    // Set up pileup reweighter
    std::cout<<"--- Beginning preparation of pileup reweighter\n";
    ztop::PUReweighter* puReweighter = new ztop::PUReweighter();
    puReweighter->setMCDistrSum12("S10");
    TString pileupInput(common::DATA_PATH_COMMON());
    pileupInput.Append("/").Append(PileupInputFILE);
    if(systematic == "PU_UP") pileupInput.ReplaceAll(".root", "_sysUp.root");
    else if(systematic == "PU_DOWN") pileupInput.ReplaceAll(".root", "_sysDown.root");
    std::cout<<"Using PU input file:\n"<<pileupInput<<std::endl;
    puReweighter->setDataTruePUInput(pileupInput.Data());
    std::cout<<"=== Finishing preparation of pileup reweighter\n\n";
    
    // Set up lepton efficiency scale factors
    LeptonScaleFactors::Systematic leptonSFSystematic(LeptonScaleFactors::nominal);
    if(systematic == "LEPT_UP") leptonSFSystematic = LeptonScaleFactors::vary_up;
    else if(systematic == "LEPT_DOWN") leptonSFSystematic = LeptonScaleFactors::vary_down;
    const LeptonScaleFactors leptonScaleFactors(ElectronSFInputFILE, MuonSFInputFILE, leptonSFSystematic);
    
    // Set up trigger efficiency scale factors (do it for all channels)
    TriggerScaleFactors::Systematic triggerSFSystematic(TriggerScaleFactors::nominal);
    if(systematic == "TRIG_UP") triggerSFSystematic = TriggerScaleFactors::vary_up;
    else if(systematic == "TRIG_DOWN") triggerSFSystematic = TriggerScaleFactors::vary_down;
    const TriggerScaleFactors triggerScaleFactors(TriggerSFInputSUFFIX,
                                                  {"ee", "emu", "mumu"},
                                                  triggerSFSystematic);
    
    // Set up btag efficiency scale factors (do it for all channels)
    BtagScaleFactors btagScaleFactors(BtagEfficiencyInputDIR,
                                      BtagEfficiencyOutputDIR,
                                      channels,
                                      systematic);
    
    // Set up JER systematic scale factors
    JetEnergyResolutionScaleFactors* jetEnergyResolutionScaleFactors(0);
    if(systematic=="JER_UP" || systematic=="JER_DOWN"){
        JetEnergyResolutionScaleFactors::Systematic jerSystematic(JetEnergyResolutionScaleFactors::vary_up);
        if(systematic == "JER_DOWN") jerSystematic = JetEnergyResolutionScaleFactors::vary_down;
        jetEnergyResolutionScaleFactors = new JetEnergyResolutionScaleFactors(jerSystematic);
    }
    
    // Set up JES systematic scale factors
    JetEnergyScaleScaleFactors* jetEnergyScaleScaleFactors(0);
    if(systematic=="JES_UP" || systematic=="JES_DOWN"){
        JetEnergyScaleScaleFactors::Systematic jesSystematic(JetEnergyScaleScaleFactors::vary_up);
        if(systematic == "JES_DOWN") jesSystematic = JetEnergyScaleScaleFactors::vary_down;
        jetEnergyScaleScaleFactors = new JetEnergyScaleScaleFactors(JesUncertaintySourceFILE, jesSystematic);
    }
    
    
    // Set up the analysis
    TopAnalysis *selector = new TopAnalysis();
    selector->SetAnalysisOutputBase(AnalysisOutputDIR);
    selector->SetKinematicReconstruction(kinematicReconstruction);
    selector->SetPUReweighter(puReweighter);
    selector->SetLeptonScaleFactors(leptonScaleFactors);
    selector->SetTriggerScaleFactors(triggerScaleFactors);
    selector->SetBtagScaleFactors(btagScaleFactors);
    selector->SetJetEnergyResolutionScaleFactors(jetEnergyResolutionScaleFactors);
    selector->SetJetEnergyScaleScaleFactors(jetEnergyScaleScaleFactors);
    
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
        
        if(!isMC && systematic!=""){
            std::cout<<"Sample is DATA, so not running again for systematic variation\n";
            continue;
        }
        if (systematic=="PDF" && (!isTopSignal || !(systematics_from_file->GetString()=="Nominal"))) {
            std::cout << "Skipping file: is not signal or not nominal -- and running PDFs\n";
            continue;
        }
        
        // Is the channel given in the file?
        if(channel_from_file->GetString() != ""){
            channels.clear();
            channels.push_back(static_cast<std::string>(channel_from_file->GetString()));
        }
        
        // Loop over channels and run selector
        for(const auto& selectedChannel : channels){
            
            // Set output file name
            TString outputfilename(filename);
            if(outputfilename.Contains('/')){
                Ssiz_t last = outputfilename.Last('/');
                outputfilename = outputfilename.Data() + last + 1;
            }
            const TString channelName = selectedChannel;
            if(!outputfilename.BeginsWith(channelName + "_")) outputfilename.Prepend(channelName + "_");
            //outputfile is now channel_filename.root
            if(dy){
                if(outputfilename.First("_dy") == kNPOS){ 
                    std::cerr << "DY variations must be run on DY samples!\n";
                    std::cerr << outputfilename << " must start with 'channel_dy'\n";
                    exit(1);
                }
                const TString dyChannel = dy == 11 ? "ee" : dy == 13 ? "mumu" : "tautau";
                outputfilename.ReplaceAll("_dy", TString("_dy").Append(dyChannel));
            }
            
            // Configure selector
            selector->SetChannel(channelName);
            selector->SetTopSignal(isTopSignal);
            selector->SetHiggsSignal(isHiggsSignal);
            selector->SetMC(isMC);
            selector->SetTrueLevelDYChannel(dy);
            if(systematic == ""){
                selector->SetSystematic(systematics_from_file->GetString());
            }
            else{
                selector->SetSystematic(systematic);
            }
            selector->SetWeightedEvents(weightedEvents);
            selector->SetSamplename(samplename->GetString());
            selector->SetGeneratorBools(samplename->GetString(), systematics_from_file->GetString());
            selector->SetOutputfilename(outputfilename);
            selector->SetRunViaTau(0);
            selector->SetClosureTest(closure, slope);
            
            // Set up nTuple chain and run selector
            TChain chain("writeNTuple/NTuple");
            chain.Add(filename);
            // chain.SetProof(); //will work from 5.34 onwards
            chain.Process(selector);
            
            // For running on PDF systematics
            if(systematic == "PDF"){
                TH1* pdfWeights = dynamic_cast<TH1*>(file.Get("EventsBeforeSelection/pdfEventWeights"));
                if(!pdfWeights){
                    std::cerr << "Error: pdfEventWeights histo missing!\n";
                    exit(831);
                }
                for(int pdf_no = 0; pdfWeights->GetBinContent(pdf_no+1) > 0; ++pdf_no){
                    if(specific_PDF >=0 && pdf_no != specific_PDF) continue;
                    selector->SetSystematic(pdfDirName(pdf_no));
                    //weightedEvents->SetBinContent(1, pdfWeights->GetBinContent(pdf_no+1));
                    selector->SetWeightedEvents(weightedEvents);
                    selector->SetPDF(pdf_no);
                    chain.Process(selector);
                }
                continue;
            }
            
            // For splitting of dileptonic ttbar in component with intermediate taus and without
            if(isTopSignal && closure == ""){
                selector->SetRunViaTau(1);
                outputfilename.ReplaceAll("signalplustau", "bgviatau");
                selector->SetOutputfilename(outputfilename);
                chain.Process(selector);
            }
        }
        file.Close();
    }
}

int main(int argc, char** argv) {
    CLParameter<std::string> opt_f("f", "Restrict to filename pattern, e.g. ttbar", false, 1, 1);
    CLParameter<std::string> opt_s("s", "Run with a systematic that runs on the nominal ntuples, e.g. 'PDF', 'PU_UP' or 'TRIG_DOWN'", false, 1, 1);
    CLParameter<int> opt_pdfno("pdf", "Run a certain PDF systematic only, sets -s PDF. Use e.g. --pdf n, where n=0 is central, 1=variation 1 up, 2=1down, 3=2up, 4=2down, ...", false, 1, 1);
    CLParameter<std::string> opt_c("c", "Specify a certain channel (ee, emu, mumu). No channel specified = run on all channels", false, 1, 1,
            [](const std::string &ch){return ch == "" || ch == "ee" || ch == "emu" || ch == "mumu";});
    CLParameter<int> opt_dy("d", "Drell-Yan mode (11 for ee, 13 for mumu, 15 for tautau)", false, 1, 1,
            [](int dy){return dy == 11 || dy == 13 || dy == 15;});
    CLParameter<std::string> opt_closure("closure", "Enable the closure test. Valid: pttop|ytop|nominal", false, 1, 1,
            [](const std::string &c){return c == "pttop" || c == "ytop" || c == "nominal";});
    CLParameter<double> opt_closureSlope("slope", "Slope for closure test, use -0.01 to 0.01 for pt and -0.4 to 0.4 for ytop", false, 1, 1,
            [](double s){return std::abs(s) < 1;});

    CLAnalyser::interpretGlobal(argc, argv);

    TString validFilenamePattern = opt_f.isSet() ? opt_f[0] : "";
    TString syst = opt_s.isSet() ? opt_s[0] : "";
    TString channel = opt_c.isSet() ? opt_c[0] : "";
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
        if (pdf_no >= 0)
            std::cout << "Running PDF variation: " << pdfDirName(pdf_no) << "\n";
        if (!(syst == "" || syst == "PDF")) {
            std::cout << "Insonsistent systematic parameter: " << syst << " cannot be used with PDF systematic!\n";
            std::exit(1);
        }
        syst = "PDF";
    }
    
//     TProof* p = TProof::Open(""); // not before ROOT 5.34
    load_Analysis(validFilenamePattern, channel, syst, pdf_no, dy, closure, slope);
//     delete p;
}
