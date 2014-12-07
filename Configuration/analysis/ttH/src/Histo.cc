#include <iostream>
#include <cstdlib>
#include <vector>
#include <string>
#include <algorithm>

#include <TString.h>

#include "Samples.h"
#include "GlobalScaleFactors.h"
#include "EventYields.h"
#include "plotterHelpers.h"
#include "Plotter.h"
#include "HistoListReader.h"
#include "higgsUtils.h"
#include "../../common/include/utils.h"
#include "../../common/include/sampleHelpers.h"
#include "../../common/include/CommandLineParameters.h"





/// Set data luminosity in pb-1
//constexpr double Luminosity = 19624.8;
//constexpr double Luminosity = 19789.;
constexpr double Luminosity = 19712.;

/// Set relative luminosity uncertainty
constexpr double LuminosityUNCERTAINTY = 0.026;





void Histo(const std::vector<std::string>& v_plot, 
           const std::vector<Channel::Channel>& v_channel,
           const std::vector<Systematic::Systematic>& v_systematic,
           const std::vector<GlobalCorrection::GlobalCorrection> v_globalCorrection,
           const DrawMode::DrawMode& drawMode)
{
    // Set up scale factors
    const bool dyCorrection = std::find(v_globalCorrection.begin(), v_globalCorrection.end(), GlobalCorrection::dy) != v_globalCorrection.end();
    const bool ttbbCorrection = std::find(v_globalCorrection.begin(), v_globalCorrection.end(), GlobalCorrection::ttbb) != v_globalCorrection.end();
    const GlobalScaleFactors* globalScaleFactors = new GlobalScaleFactors(v_channel, v_systematic, Luminosity, LuminosityUNCERTAINTY, dyCorrection, ttbbCorrection);
    
    // Access all samples
    const Samples samples("FileLists_plot", v_channel, v_systematic, globalScaleFactors);
    
    // Produce event yields
    const EventYields eventYields("EventYields", samples);
    
    // Create Plotter
    Plotter plotter("Plots", samples, drawMode);
    
    // Access the histoList specifying printing parameters of histograms
    const std::string histoListFile(tth::DATA_PATH_TTH() + "/" + "HistoList_control");
    const HistoListReader histoList(histoListFile.data());
    if(histoList.isZombie()){
        std::cerr<<"Error in Histo! Cannot find HistoList with name: "<<histoListFile<<"\n...break\n"<<std::endl;
        exit(12);
    }
    
    // Loop over all histograms in histoList and print them
    std::cout<<"--- Beginning with the plotting\n\n";
    for(auto it = histoList.begin(); it != histoList.end(); ++it){
        
        // Access plot properties from histoList and check whether histogram name contains name pattern
        const PlotProperties& plotProperties = it->second;
        std::cout << "\nchecking " << plotProperties.name << std::endl;
        bool found = false;
        for(const auto& plot : v_plot){
            if(plot.size() && plot[0] == '+'){
                if(plotProperties.name.CompareTo(&plot[1], TString::kIgnoreCase) == 0){
                    found = true;
                    break;
                }
            }
            else if(plotProperties.name.Contains(plot, TString::kIgnoreCase)){
                found = true;
                break;
            }
        }
        if(!found){
            std::cout<<"... no histograms found, continue with next\n";
            continue;
        }
        
        // Set plot properties
        plotter.setOptions(plotProperties.name, plotProperties.specialComment, plotProperties.ytitle, plotProperties.xtitle, 
                               plotProperties.rebin, plotProperties.do_dyscale, plotProperties.logX, plotProperties.logY, 
                               plotProperties.ymin, plotProperties.ymax, plotProperties.xmin, plotProperties.xmax,
                               plotProperties.bins, plotProperties.xbinbounds, plotProperties.bincenters);
        
        // Loop over all systematics and all channels and write histograms
        plotter.producePlots();
    }
    std::cout<<"\n=== Finishing with the plotting\n\n";
}



/// All systematics allowed for plotting
namespace Systematic{
    const std::vector<Type> allowedSystematics = {
        nominal, all, allAvailable,
        pu, lept, trig,
        jer, jes,
        btag, 
        btagPt, btagEta,
        btagLjet, 
        btagLjetPt, btagLjetEta,
        btagDiscrBstat1, btagDiscrBstat2,
        btagDiscrLstat1, btagDiscrLstat2,
        btagDiscrPurity,
        kin,
        lumi,
        xsec_tt2b, xsec_ttcc,
        topPt,
        mass, match, scale,
        powheg, powhegHerwig, mcatnlo, perugia11, perugia11NoCR,
        pdf
    };
}



int main(int argc, char** argv){
    
    // Get and check configuration parameters
    CLParameter<std::string> opt_plot("p", "Name (pattern) of plot; multiple patterns possible; use '+Name' to match name exactly", false, 1, 100);
    CLParameter<std::string> opt_channel("c", "Specify channel(s), valid: emu, ee, mumu, combined. Default: all channels", false, 1, 4,
        common::makeStringCheck(Channel::convert(Channel::allowedChannelsPlotting)));
    CLParameter<std::string> opt_systematic("s", "Systematic variation - default is Nominal, use 'all' for all, use 'allAvailable' for all available", false, 1, 100,
        common::makeStringCheckBegin(Systematic::convertType(Systematic::allowedSystematics)));
    CLParameter<std::string> opt_globalCorrection("g", "Specify global correction, valid: empty argument for none, Drell-Yan (dy), tt+HF (ttbb). Default: dy", false, 0, 2,
        common::makeStringCheck(GlobalCorrection::convert(GlobalCorrection::allowedGlobalCorrections)));
    CLParameter<std::string> opt_drawMode("d", "Specify draw mode of Higgs sample, valid: stacked, overlaid, scaledoverlaid. Default: scaledoverlaid", false, 1, 1,
        common::makeStringCheck(DrawMode::convert(DrawMode::allowedDrawModes)));
    CLAnalyser::interpretGlobal(argc, argv);
    
    // Set up plots
    std::vector<std::string> v_plot {""};
    if(opt_plot.isSet()){
        v_plot = opt_plot.getArguments();
        std::cout<< "Processing only histograms containing in name: ";
        for(auto plot : v_plot)std::cout<< plot << " ";
        std::cout << "\n\n";
    }
    
    // Set up channels
    std::vector<Channel::Channel> v_channel(Channel::allowedChannelsPlotting);
    if(opt_channel.isSet()) v_channel = Channel::convert(opt_channel.getArguments());
    std::cout << "Processing channels: ";
    for(auto channel : v_channel) std::cout << Channel::convert(channel) << " ";
    std::cout << "\n\n";
    
    // Set up systematics
    std::vector<Systematic::Systematic> v_systematic = Systematic::allowedSystematicsAnalysis(Systematic::allowedSystematics);
    if(opt_systematic.isSet()){
        if(opt_systematic[0] == Systematic::convertType(Systematic::all))
            ; // do nothing
        else if(opt_systematic[0] == Systematic::convertType(Systematic::allAvailable))
            v_systematic = common::findSystematicsFromFilelists("FileLists_plot", v_channel, v_systematic);
        else
            v_systematic = Systematic::setSystematics(opt_systematic.getArguments());
    }
    else{
        v_systematic.clear();
        v_systematic.push_back(Systematic::nominalSystematic());
    }
    std::cout << "Processing systematics: "; 
    for(auto systematic : v_systematic) std::cout << systematic.name() << " ";
    std::cout << "\n\n";
    
    // Set up global corrections
    std::vector<GlobalCorrection::GlobalCorrection> v_globalCorrection({GlobalCorrection::dy});
    if(opt_globalCorrection.isSet()) v_globalCorrection = GlobalCorrection::convert(opt_globalCorrection.getArguments());
    std::cout << "\n";
    std::cout << "Using global corrections: ";
    for(auto globalCorrection : v_globalCorrection) std::cout << GlobalCorrection::convert(globalCorrection) << " ";
    std::cout << "\n\n";
    
    // Set up draw mode (default is scaledOverlaid)
    DrawMode::DrawMode drawMode(DrawMode::undefined);
    if(opt_drawMode.isSet()) drawMode = DrawMode::convert(opt_drawMode[0]);
    else drawMode = DrawMode::scaledoverlaid;
    std::cout << "\n";
    std::cout << "Draw mode of Higgs sample: "<<DrawMode::convert(drawMode);
    std::cout << "\n\n";
    
    // Start analysis
    Histo(v_plot, v_channel, v_systematic, v_globalCorrection, drawMode);
}



