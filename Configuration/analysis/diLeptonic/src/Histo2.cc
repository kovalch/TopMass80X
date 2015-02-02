#include <vector>
#include <algorithm> 
#include <fstream>

#include <TString.h>

#include "Samples.h"
#include "plotterHelpers.h"
#include "GlobalScaleFactors.h"
#include "Plotter.h"
#include "FinalPlot.h"
#include "PlotterBTop.h"
#include <ttbarUtils.h>
#include "../../common/include/utils.h"
#include "../../common/include/sampleHelpers.h"
#include "../../common/include/CommandLineParameters.h"
#include "utils.h"




/// Set data luminosity in pb-1
//constexpr double Luminosity = 19624.8;
//constexpr double Luminosity = 19789.;
constexpr double Luminosity = 19712.;

///Inclusive top xsec in 
constexpr double topxsec = 245.102;

std::vector<std::vector<TString>> vv_plotName;



/// All systematics allowed for plotting
namespace Systematic{
    const std::vector<Type> allowedSystematics = {
        nominal, all,
        pu, lept, trig,
        jer, jes,
        btag, btagPt, btagEta,
        btagLjet, btagLjetPt, btagLjetEta,
        kin,
    };
}



int main(int argc, char** argv){
    
    // Get and check configuration parameters
     CLParameter<std::string> opt_analysis("a", "Specify analysis, valid: 2d - for double differential , btop - for boosted top.",true, 1,1,
                                          common::makeStringCheck( {"2d","btop"} ));
     CLParameter<std::string> opt_plotter("p", "Specify plotting stage, valid: p - for Plotter , fp - for FinalPlot, without option - for both.",false, 0,1,
                                          common::makeStringCheck( {"p","fp"} ));
     CLParameter<std::string> opt_channel("c", "Specify channel(s), valid: emu, ee, mumu, combined. Default: all channels", false, 1, 4,
                                          common::makeStringCheck(Channel::convert(Channel::realChannels)));
     CLParameter<std::string> opt_systematic("s", "Systematic variation - default is Nominal, use 'all' for all", false, 1, 100,
                                             common::makeStringCheckBegin(Systematic::convertType(Systematic::allowedSystematics)));
     CLParameter<std::string> opt_globalCorrection("g", "Specify global correction, valid: empty argument for none, Drell-Yan (dy), tt+HF (ttbb). Default: dy", false, 0, 2,
                                                   common::makeStringCheck(GlobalCorrection::convert(GlobalCorrection::allowedGlobalCorrections)));
     CLAnalyser::interpretGlobal(argc, argv);
    
     styleUtils::setHHStyle(*gStyle);

     
     //Set up plotter
     bool runPlotter=true;
     bool runFinalPlot=true;
     if(opt_plotter.isSet()){
         if(opt_plotter.getArguments().at(0)=="p")runFinalPlot=false;
         if(opt_plotter.getArguments().at(0)=="fp")runPlotter=false;
     }
     
    // Set up channels
    //std::vector<Channel::Channel> v_channel(Channel::realChannels);
    std::vector<Channel::Channel> v_channel(1,Channel::emu);//emu
    if(opt_channel.isSet()) v_channel = Channel::convert(opt_channel.getArguments());
    std::cout << "Processing channels: ";
    for(auto channel : v_channel) std::cout << Channel::convert(channel) << " ";
    std::cout << "\n\n";
    
    // Set up systematics
    std::vector<Systematic::Systematic> v_systematic = Systematic::allowedSystematicsAnalysis(Systematic::allowedSystematics);
    if(opt_systematic.isSet() && opt_systematic[0]!=Systematic::convertType(Systematic::all)) v_systematic = Systematic::setSystematics(opt_systematic.getArguments());
    else if(opt_systematic.isSet() && opt_systematic[0]==Systematic::convertType(Systematic::all)); // do nothing
    else{v_systematic.clear(); v_systematic.push_back(Systematic::nominalSystematic());}
    std::cout << "Processing systematics (use >>-s all<< to process all known systematics): "; 
    for(auto systematic : v_systematic) std::cout << systematic.name() << " ";
    std::cout << "\n\n";
    
    Samples SAMPLES;
    
    if(runPlotter){
        
        // Set up global corrections
        std::vector<GlobalCorrection::GlobalCorrection> v_globalCorrection({GlobalCorrection::dy});
        if(opt_globalCorrection.isSet()) v_globalCorrection = GlobalCorrection::convert(opt_globalCorrection.getArguments());
        std::cout << "\n";
        std::cout << "Using global corrections: ";
        for(auto globalCorrection : v_globalCorrection) std::cout << GlobalCorrection::convert(globalCorrection) << " ";
        std::cout << "\n\n";
        
        // Set up scale factors
        const bool dyCorrection = std::find(v_globalCorrection.begin(), v_globalCorrection.end(), GlobalCorrection::dy) != v_globalCorrection.end();
        const bool ttbbCorrection = std::find(v_globalCorrection.begin(), v_globalCorrection.end(), GlobalCorrection::ttbb) != v_globalCorrection.end();
        const GlobalScaleFactors* globalScaleFactors = new GlobalScaleFactors(v_channel, v_systematic, Luminosity, dyCorrection, ttbbCorrection);
    
        // Access all samples   
        const Samples samples("FileLists", v_channel, v_systematic, globalScaleFactors); // "FileLists" is a folder in diLeptonic, to create this folder : 
        SAMPLES = samples;
    }
    
    if(opt_analysis.getArguments().at(0)=="2d")
    {
        std::cout << "Running 2d analysis ..."<<  std::endl;
        
        const std::string nameListFile(common::CMSSW_BASE() + "/src/TopAnalysis/Configuration/analysis/diLeptonic/" + "NameList");
        ttbar::setPlotNames(nameListFile,vv_plotName);

        if(runPlotter){
            std::cout<<"--- Beginning with the Plotter\n\n";
            Plotter generalPlot(SAMPLES,Luminosity,topxsec);
            for(auto v_plotName : vv_plotName){
                generalPlot.setOptions(v_plotName);
                generalPlot.producePlots();
            }
            std::cout<<"\n=== Finishing with the plotting\n\n";
        }
        if(runFinalPlot){
            std::cout<<"--- Beginning with the FinalPlot plotting\n\n";
            FinalPlot finalPlot(v_channel,v_systematic,Luminosity,topxsec);
            for(auto v_plotName : vv_plotName){
                finalPlot.setOptions(v_plotName);
                finalPlot.producePlots(".pdf");
            }
            std::cout<<"\n=== Finishing with the FinalPlot plotting\n\n";
        }
        
    }
    else if(opt_analysis.getArguments().at(0)=="btop"){
        std::cout << "Running btop analysis ..."<<  std::endl;
        
        const std::string nameListFile(common::CMSSW_BASE() + "/src/TopAnalysis/Configuration/analysis/diLeptonic/" + "NameListBTop");
        ttbar::setPlotNames(nameListFile,vv_plotName);
        
        if(runPlotter){
            std::cout<<"--- Beginning with the Plotter\n\n";
            PlotterBTop generalPlot(SAMPLES,Luminosity,topxsec);
            for(auto v_plotName : vv_plotName){
                std::cout << v_plotName.at(0) << std::endl;
                generalPlot.setOptions(v_plotName);
                generalPlot.producePlots();
            }
            std::cout<<"\n=== Finishing with the plotting\n\n";
        }
        if(runFinalPlot){
            std::cout<<"--- Beginning with the FinalPlot plotting\n\n";
            
            std::cout<<"\n=== Finishing with the FinalPlot plotting\n\n";
        }
    }
    
    
    
    
    

}





