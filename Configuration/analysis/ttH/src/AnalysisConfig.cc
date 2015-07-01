#include <iostream>

#include "AnalysisConfig.h"
#include "../../common/include/TextfileReader.h"
#include "../../common/include/sampleHelpers.h"





AnalysisConfig::General::General():
era_(Era::undefined),
luminosity_(-999.),
luminosityUncertainty_(-999.)
{}



void AnalysisConfig::General::print()const
{
    std::cout<<"Beginning printing AnalysisConfig general\n"
             <<"\tera: "<<Era::convert(era_)<<"\n"
             <<"\tluminosity: "<<luminosity_<<"\n"
             <<"\tluminosity uncertainty: "<<luminosityUncertainty_<<"\n"
             <<"Finishing printing AnalysisConfig general\n\n";
}



AnalysisConfig::Corrections::Corrections():
pileupInputFile_(""),
triggerSFInputSuffix_(""),
electronSFInputFile_(""),
muonSFInputFile_(""),
jerUncertaintySourceName_(""),
jesUncertaintySourceFile_(""),
btagCorrectionMode_(Btag::undefinedCorrectionMode),
btagHeavyFlavourFile_(""),
btagLightFlavourFile_(""),
mvaMetRecoilDataFile_(""),
mvaMetRecoilMcFile_("")
{}



void AnalysisConfig::Corrections::print()const
{
    std::cout<<"Beginning printing AnalysisConfig corrections\n"
             <<"\tpileup file: "<<pileupInputFile_<<"\n"
             <<"\ttrigger SF suffix: "<<triggerSFInputSuffix_<<"\n"
             <<"\telectron SF file: "<<electronSFInputFile_<<"\n"
             <<"\tmuon SF file: "<<muonSFInputFile_<<"\n"
             <<"\tJER name: "<<jerUncertaintySourceName_<<"\n"
             <<"\tJES file: "<<jesUncertaintySourceFile_<<"\n"
             <<"\tb-tag correction mode: "<<Btag::convertCorrectionMode(btagCorrectionMode_)<<"\n"
             <<"\tb-tag HF file: "<<btagHeavyFlavourFile_<<"\n"
             <<"\tb-tag LF file: "<<btagLightFlavourFile_<<"\n"
             <<"\tMVA MET recoil data file: "<<mvaMetRecoilDataFile_<<"\n"
             <<"\tMVA MET recoil MC file: "<<mvaMetRecoilMcFile_<<"\n"
             <<"Finishing printing AnalysisConfig corrections\n\n";
}



AnalysisConfig::Selections::Selections():
leptonEtaCut_(-999.),
leptonPtCut_(-999.),
jetEtaCut_(-999.),
jetPtCut_(-999.),
deltaRLeptonJetCut_(-999.),
lead2JetPtCut_(-999.),
btagAlgorithm_(Btag::undefinedAlgorithm),
btagWorkingPoint_(Btag::undefinedWP),
mvaMet_(false),
metCut_(-999.),
genJetEtaCut_(-999.),
genJetPtCut_(-999.),
genDeltaRLeptonJetCut_(-999.)
{}



void AnalysisConfig::Selections::print()const
{
    std::cout<<"Beginning printing AnalysisConfig selections\n"
             <<"\tlepton eta: "<<leptonEtaCut_<<"\n"
             <<"\tlepton pt: "<<leptonPtCut_<<"\n"
             <<"\tjet eta: "<<jetEtaCut_<<"\n"
             <<"\tjet pt: "<<jetPtCut_<<"\n"
             <<"\tdeltaR(lepton, jet): "<<deltaRLeptonJetCut_<<"\n"
             <<"\tleading 2 jet pt: "<<lead2JetPtCut_<<"\n"
             <<"\tb-tag algorithm: "<<Btag::convertAlgorithm(btagAlgorithm_)<<"\n"
             <<"\tb-tag working point: "<<Btag::convertWorkingPoint(btagWorkingPoint_)<<"\n"
             <<"\tuse MVA MET: "<<mvaMet_<<"\n"
             <<"\tMET et: "<<metCut_<<"\n"
             <<"\tgenJet eta: "<<genJetEtaCut_<<"\n"
             <<"\tgenJet pt: "<<genJetPtCut_<<"\n"
             <<"\tdeltaR(genJet, genLepton): "<<genDeltaRLeptonJetCut_<<"\n"
             <<"Finishing printing AnalysisConfig selections\n\n";
}



AnalysisConfig::AnalysisConfig(const std::string& configfilename)
{
    std::cout<<"--- Beginning setting up analysis config\n";
    
    TextfileReader textfileReader;
    
    // Read general info
    textfileReader.setStartMarker("[ general ]");
    textfileReader.setEndMarker("[ end - general ]");
    textfileReader.readFile(configfilename);
    general_.era_ = Era::convert(textfileReader.getValue<std::string>("era"));
    general_.luminosity_ = textfileReader.getValue<double>("luminosity");
    general_.luminosityUncertainty_ = textfileReader.getValue<double>("luminosityUncertainty");
    textfileReader.clear();
    
    // Read corrections
    textfileReader.setStartMarker("[ corrections ]");
    textfileReader.setEndMarker("[ end - corrections ]");
    textfileReader.readFile(configfilename);
    textfileReader.setRequireValues(false);
    corrections_.pileupInputFile_ = textfileReader.getValue<std::string>("pileupInputFile", "");
    corrections_.triggerSFInputSuffix_ = textfileReader.getValue<std::string>("triggerSFInputSuffix", "");
    corrections_.electronSFInputFile_ = textfileReader.getValue<std::string>("electronSFInputFile", "");
    corrections_.muonSFInputFile_ = textfileReader.getValue<std::string>("muonSFInputFile", "");
    corrections_.jerUncertaintySourceName_ = textfileReader.getValue<std::string>("jerUncertaintySourceName", "");
    corrections_.jesUncertaintySourceFile_ = textfileReader.getValue<std::string>("jesUncertaintySourceFile", "");
    corrections_.btagCorrectionMode_ = Btag::convertCorrectionMode(textfileReader.getValue<std::string>("btagCorrectionMode", "noCorrection"));
    corrections_.btagHeavyFlavourFile_ = textfileReader.getValue<std::string>("btagHeavyFlavourFile", "");
    corrections_.btagLightFlavourFile_ = textfileReader.getValue<std::string>("btagLightFlavourFile", "");
    corrections_.mvaMetRecoilDataFile_ = textfileReader.getValue<std::string>("mvaMetRecoilDataFile", "");
    corrections_.mvaMetRecoilMcFile_ = textfileReader.getValue<std::string>("mvaMetRecoilMcFile", "");
    textfileReader.setRequireValues(true);
    textfileReader.clear();
    
    // Read object selection
    textfileReader.setStartMarker("[ objectSelection ]");
    textfileReader.setEndMarker("[ end - objectSelection ]");
    textfileReader.readFile(configfilename);
    selections_.leptonEtaCut_ = textfileReader.getValue<double>("leptonEtaCut");
    selections_.leptonPtCut_ = textfileReader.getValue<double>("leptonPtCut");
    selections_.jetEtaCut_ = textfileReader.getValue<double>("jetEtaCut");
    selections_.jetPtCut_ = textfileReader.getValue<double>("jetPtCut");
    selections_.deltaRLeptonJetCut_ = textfileReader.getValue<double>("deltaRLeptonJetCut");
    selections_.lead2JetPtCut_ = textfileReader.getValue<double>("lead2JetPtCut");
    selections_.btagAlgorithm_ = Btag::convertAlgorithm(textfileReader.getValue<std::string>("btagAlgorithm"));
    selections_.btagWorkingPoint_ = Btag::convertWorkingPoint(textfileReader.getValue<std::string>("btagWorkingPoint"));
    selections_.mvaMet_ = textfileReader.getValue<bool>("mvaMet");
    selections_.metCut_ = textfileReader.getValue<double>("metCut");
    selections_.genJetEtaCut_ = textfileReader.getValue<double>("genJetEtaCut");
    selections_.genJetPtCut_ = textfileReader.getValue<double>("genJetPtCut");
    selections_.genDeltaRLeptonJetCut_ = textfileReader.getValue<double>("genDeltaRLeptonJetCut");
    textfileReader.clear();
    
    std::cout<<"=== Finishing setting up analysis config\n\n";
}



void AnalysisConfig::print()const
{
    general_.print();
    corrections_.print();
    selections_.print();
}




