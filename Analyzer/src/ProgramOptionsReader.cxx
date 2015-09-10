/*
 * ProgramOptionsReader.cxx
 *
 *  Created on: Oct 29, 2012
 *      Author: eschliec
 */

#include "ProgramOptionsReader.h"

#include <fstream>

ProgramOptionsReader::ProgramOptionsReader(int ac, char** av) {
  ReadProgramOptions(ac, av);
}

ProgramOptionsReader::~ProgramOptionsReader() {
  delete programOptions_;
}

boost::program_options::variables_map* ProgramOptionsReader::programOptions_(0);

void
ProgramOptionsReader::ReadProgramOptions(int ac, char** av) {
  if(!programOptions_){
    programOptions_ = new boost::program_options::variables_map();
    // Declare the supported options
    boost::program_options::options_description desc("Allowed options", 150);
    desc.add_options()
        ("help,h", "produce help message")
        ("test", boost::program_options::value<bool>()->default_value(false),
            "Speed-ups for code testing purposes (restricts sample size)\n"
        )
        ("seed", boost::program_options::value<int>()->default_value(-1), "Random seed for pseudo-experiments")
        ("outPath,o", boost::program_options::value<std::string>()->default_value(""))
        ("method,m", boost::program_options::value<std::string>()->default_value("Ideogram"),
            "Top mass measurement method\n"
            "  GenMatch: \tGaussian fit of correct permutations (MC only)\n"
            "  Ideogram: \tJoint likelihood fit of mt and JES given data sample\n"
            "  IdeogramNew: \tNew interface\n"
            "  IdeogramMin: \tMinimizer function\n"
            "  RooFit:   \tTemplate fit of mt (JES, fSig) given data sample using RooFit"
        )
        ("minPlot", boost::program_options::value<bool>()->default_value(false),
            "Plot IdeogramMin result\n"
        )
        ("temPlot", boost::program_options::value<bool>()->default_value(false),
            "Plot IdeogramMin templates\n"
        )
        ("constrainJSF", boost::program_options::value<bool>()->default_value(true),
            "Add Gaussian constraint to JSF\n"
        )
        ("constrainJSFweight", boost::program_options::value<double>()->default_value(0.5),
            "Weight of Gaussian constraint to JSF\n"
        )
        ("constrainJSFsigma", boost::program_options::value<double>()->default_value(-1),
            "Sigma value of Gaussian constraint to JSF\n"
        )
        ("task,t", boost::program_options::value<std::string>()->default_value("sm"),
            "Task to be done\n"
            "  sm: \tSingle measurement based on input parameters\n"
            "  pe: \tPerform pseudo-experiments, additional settings necessary\n"
            "  hc: \tPerform any hardcoded tasks in TopMass constructor\n"
            "  diff: \tDifferential top mass, specify binning option"
        )
        ("input,i", boost::program_options::value<std::string>()->default_value("MJP12*_v1_data"),
            "Identifier of input file to be analyzed: Z2_F11_172_5_sig"
        )
        ("channel,c", boost::program_options::value<std::string>()->default_value("alljets"),
          "Channel\n"
          "  electron: \t\n"
          "  muon: \t\n"
          "  lepton: \telectron+muon\n"
          "  alljets: \talljets"
          "  hamburg: \telectron+muon+alljets"
        )
        ("topBranchName", boost::program_options::value<std::string>()->default_value("top."),
          "Top branch name in eventTree\n"
          "  top: \t(default value)\n"
          "  bRegTop: \tchoose branch containing fit with b-regression applied\n"
        )
        ("topTreeName", boost::program_options::value<std::string>()->default_value("analyzeHitFit/eventTre"),
          "Top tree name in file\n"
          "  analyzeHitFit/eventTree: \tused for l+jets (default value)\n"
          "  analyzeKinFit/eventTree: \tused for all-jets\n"
        )
        ("binning,b", boost::program_options::value<std::string>()->default_value("deltaThetaHadWHadB"),
          "Phasespace binning\n"
          "  deltaThetaHadWHadB\n"
          "  hadTopPt\n"
          "  hadBEta"
        )
        ("weight,w", boost::program_options::value<std::string>()->default_value("muWeight*bWeight*PUWeight"),
                "Event weight used in pseudo-experiments")
        ("mass,M", boost::program_options::value<double>()->default_value(172.5), "Input top mass for pseudo-experiments")
        ("jes,J", boost::program_options::value<double>()->default_value(1.0), "Input JES for pseudo-experiments")
        ("bdisc,B", boost::program_options::value<double>()->default_value(0.679), "Threshold for b-jets")
        ("fsig,f", boost::program_options::value<double>()->default_value(0.767995487), "Input signal fraction for pseudo-experiments")
        ("fsigLept", boost::program_options::value<double>()->default_value(0.), "Input signal fraction for pseudo-experiments for LeptonJets channel")
        ("fsigJets", boost::program_options::value<double>()->default_value(0.), "Input signal fraction for pseudo-experiments for AllJets channel")
        ("lumi,L", boost::program_options::value<double>()->default_value(0.0), "Luminosity for each pseudo-experiment")
        ("lumiLept", boost::program_options::value<double>()->default_value(0.0), "Luminosity for each pseudo-experiment in Lepton+Jets")
        ("lumiJets", boost::program_options::value<double>()->default_value(0.0), "Luminosity for each pseudo-experiment in Alljets")
        ("number,N", boost::program_options::value<int>()->default_value(10000), "Number of pseudo-experiments per job")
        ("walltime,W", boost::program_options::value<double>()->default_value(10), "set walltime limit for pseudo-experiments in minutes")
        ("shape,S", boost::program_options::value<double>()->default_value(0.0), "Background shape scaling factor for top distribution")
        ("shape2,T", boost::program_options::value<double>()->default_value(0.0), "Background shape scaling factor for w distribution")
        ("permu,P", boost::program_options::value<double>()->default_value(0.0), "Change permutation fractions by: fUN-P, fWP+0.5P, fCP+0.5P")
        ("fastsim,F", boost::program_options::value<int>()->default_value(0), "use additional calibration for FastSim")
        ("preliminary,p", boost::program_options::value<int>()->default_value(1), "use \"Preliminary\" label for plots")
        ("cmsenergy,e", boost::program_options::value<int>()->default_value(7), "cms energy to be used (for example for plots)")
        ("pullWidth", boost::program_options::value<double>()->default_value(1.0), "pull width correction factor")
        ("muo_pullWidth", boost::program_options::value<double>()->default_value(1.0), "pull width correction factor (mu+jets)")
        ("ele_pullWidth", boost::program_options::value<double>()->default_value(1.0), "pull width correction factor (e+jets)")
        ("analysisConfig.selection", boost::program_options::value<std::string>())
        ("analysisConfig.selectionLept", boost::program_options::value<std::string>()->default_value(""))
        ("analysisConfig.selectionJets", boost::program_options::value<std::string>()->default_value(""))
        ("analysisConfig.selectionCP", boost::program_options::value<std::string>())
        ("analysisConfig.selectionWP", boost::program_options::value<std::string>())
        ("analysisConfig.selectionUN", boost::program_options::value<std::string>())
        ("analysisConfig.samplePath", boost::program_options::value<std::string>())
        ("analysisConfig.samplePathLept", boost::program_options::value<std::string>()->default_value(""))
        ("analysisConfig.samplePathJets", boost::program_options::value<std::string>()->default_value(""))
        ("analysisConfig.var1", boost::program_options::value<std::string>())
        ("analysisConfig.var2", boost::program_options::value<std::string>())
        ("analysisConfig.var3", boost::program_options::value<std::string>())
        ("analysisConfig.var4", boost::program_options::value<std::string>())
        ("analysisConfig.maxPermutations", boost::program_options::value<int>()->default_value(666666))
        ("analysisConfig.activeBranches", boost::program_options::value<std::string>()->default_value("*"))
        ("analysisConfig.plotsToDraw", boost::program_options::value<std::string>()->default_value("StandardPlots|UETunePlots|MCGeneratorPlots|JESVariationPlots|SignalModellingPlots|BasicMasses|LightPulls|BJetPulls|MixedPulls|ExtraPlotsFitCombTypeEtc|JetPts|TopPts|FitPts|LeptonJetsExtra|JetDetails|EventObservables"))
        ("analysisConfig.ratioYMin", boost::program_options::value<double>()->default_value(0.49))
        ("analysisConfig.ratioYMax", boost::program_options::value<double>()->default_value(1.51))
        ("analysisConfig.plotPermutations", boost::program_options::value<bool>()->default_value(true))
        ("analysisConfig.renameCombinationTypes", boost::program_options::value<std::string>()->default_value(" correct| wrong| unmatched"))
        ("analysisConfig.redefineCombinationTypeColorShifts", boost::program_options::value<std::string>()->default_value("0|-8|-11"))
        ("analysisConfig.normalizeToData", boost::program_options::value<bool>()->default_value(false))
        ("analysisConfig.signalFactor", boost::program_options::value<double>()->default_value(1.0))
        ("analysisConfig.backgroundSamplesLept", boost::program_options::value<std::string>()->default_value("Summer12_WJets|Summer12_ZJets|Summer12_singleTop|Summer12_QCD|Summer12_Diboson"))
        ("analysisConfig.backgroundFactorsLept", boost::program_options::value<std::string>()->default_value("1|1|1|1|1"))
        ("templates.fSig", boost::program_options::value<double>()->default_value(0.0))
        ("templates.fCP", boost::program_options::value<double>()->default_value(0.0))
        ("templates.fWP", boost::program_options::value<double>()->default_value(0.0))
        ("templates.fUN", boost::program_options::value<double>()->default_value(0.0))
        ("templates.muo_fCP", boost::program_options::value<double>()->default_value(0.0))
        ("templates.muo_fWP", boost::program_options::value<double>()->default_value(0.0))
        ("templates.muo_fUN", boost::program_options::value<double>()->default_value(0.0))
        ("templates.ele_fCP", boost::program_options::value<double>()->default_value(0.0))
        ("templates.ele_fWP", boost::program_options::value<double>()->default_value(0.0))
        ("templates.ele_fUN", boost::program_options::value<double>()->default_value(0.0))
        ("templates.parsCP", boost::program_options::value<std::string>()->default_value(""))
        ("templates.parsWP", boost::program_options::value<std::string>()->default_value(""))
        ("templates.parsUN", boost::program_options::value<std::string>()->default_value(""))
        ("templates.parsCPJES", boost::program_options::value<std::string>()->default_value(""))
        ("templates.parsWPJES", boost::program_options::value<std::string>()->default_value(""))
        ("templates.parsUNJES", boost::program_options::value<std::string>()->default_value(""))
        ("templates.muo_parsCP", boost::program_options::value<std::string>()->default_value(""))
        ("templates.muo_parsWP", boost::program_options::value<std::string>()->default_value(""))
        ("templates.muo_parsUN", boost::program_options::value<std::string>()->default_value(""))
        ("templates.muo_parsCPJES", boost::program_options::value<std::string>()->default_value(""))
        ("templates.muo_parsWPJES", boost::program_options::value<std::string>()->default_value(""))
        ("templates.muo_parsUNJES", boost::program_options::value<std::string>()->default_value(""))
        ("templates.ele_parsCP", boost::program_options::value<std::string>()->default_value(""))
        ("templates.ele_parsWP", boost::program_options::value<std::string>()->default_value(""))
        ("templates.ele_parsUN", boost::program_options::value<std::string>()->default_value(""))
        ("templates.ele_parsCPJES", boost::program_options::value<std::string>()->default_value(""))
        ("templates.ele_parsWPJES", boost::program_options::value<std::string>()->default_value(""))
        ("templates.ele_parsUNJES", boost::program_options::value<std::string>()->default_value(""))
        ("templates.parsBKG", boost::program_options::value<std::string>()->default_value(""))
        ("templates.parsBKGJES", boost::program_options::value<std::string>()->default_value(""))
        ("templates.maxTopMass", boost::program_options::value<double>()->default_value(215.))
        ("calibration.massOffset", boost::program_options::value<std::string>()->default_value("0.0"))
        ("calibration.massSlopeMass", boost::program_options::value<std::string>()->default_value("0.0"))
        ("calibration.massSlopeJES", boost::program_options::value<std::string>()->default_value("0.0"))
        ("calibration.massSlopeMassJES", boost::program_options::value<std::string>()->default_value("0.0"))
        ("calibration.jesOffset", boost::program_options::value<std::string>()->default_value("0.0"))
        ("calibration.jesSlopeMass", boost::program_options::value<std::string>()->default_value("0.0"))
        ("calibration.jesSlopeJES", boost::program_options::value<std::string>()->default_value("0.0"))
        ("calibration.jesSlopeMassJES", boost::program_options::value<std::string>()->default_value("0.0"))
        ("calibration.muo_massOffset", boost::program_options::value<std::string>()->default_value("0.0"))
        ("calibration.muo_massSlopeMass", boost::program_options::value<std::string>()->default_value("0.0"))
        ("calibration.muo_massSlopeJES", boost::program_options::value<std::string>()->default_value("0.0"))
        ("calibration.muo_massSlopeMassJES", boost::program_options::value<std::string>()->default_value("0.0"))
        ("calibration.muo_jesOffset", boost::program_options::value<std::string>()->default_value("0.0"))
        ("calibration.muo_jesSlopeMass", boost::program_options::value<std::string>()->default_value("0.0"))
        ("calibration.muo_jesSlopeJES", boost::program_options::value<std::string>()->default_value("0.0"))
        ("calibration.muo_jesSlopeMassJES", boost::program_options::value<std::string>()->default_value("0.0"))
        ("calibration.ele_massOffset", boost::program_options::value<std::string>()->default_value("0.0"))
        ("calibration.ele_massSlopeMass", boost::program_options::value<std::string>()->default_value("0.0"))
        ("calibration.ele_massSlopeJES", boost::program_options::value<std::string>()->default_value("0.0"))
        ("calibration.ele_massSlopeMassJES", boost::program_options::value<std::string>()->default_value("0.0"))
        ("calibration.ele_jesOffset", boost::program_options::value<std::string>()->default_value("0.0"))
        ("calibration.ele_jesSlopeMass", boost::program_options::value<std::string>()->default_value("0.0"))
        ("calibration.ele_jesSlopeJES", boost::program_options::value<std::string>()->default_value("0.0"))
        ("calibration.ele_jesSlopeMassJES", boost::program_options::value<std::string>()->default_value("0.0"))
        ("skimmer.mcweight", boost::program_options::value<double>()->default_value(-1.0))
        ;

    boost::program_options::store(boost::program_options::parse_command_line(ac, av, desc), *programOptions_);

    std::string channel = programOptions_->operator[]("channel").as<std::string>();
    std::string fileNameSnippet;
    int length = -1;
    if (((length = 8) && !strncmp(channel.c_str(), "electron", length)) ||
        ((length = 4) && !strncmp(channel.c_str(), "muon"    , length)) ||
        ((length = 6) && !strncmp(channel.c_str(), "lepton"  , length))) {
  		fileNameSnippet = "LeptonJets";
   		fileNameSnippet += channel.substr(length);
    }
    else if (((length = 7) && !strncmp(channel.c_str(), "alljets", length))) {
      fileNameSnippet = "AllJets";
      fileNameSnippet += channel.substr(length);
    }
    else if (((length = 7) && !strncmp(channel.c_str(), "hamburg", length))) {
      fileNameSnippet = "Hamburg";
      fileNameSnippet += channel.substr(length);
    }
    else {
      std::cerr << "Stopping analysis! Specified decay channel *" << channel << "* not known!" << std::endl;
      return;
    }
    // Read in additional information to supersede default config, separated by "_"
    // e.g. c lepton_BReg will read in Configuration_LeptonJets_BReg.conf and set values
    // before default config is read in.
    size_t splitPos = fileNameSnippet.find_last_of("_"); //equal to string::npos if not found
    if(fileNameSnippet.find_first_of("_")!=splitPos) std::cerr << "Stopping analysis! Check input decay channel: Only one _ character supported to read in additional information" << std::endl;
    if(splitPos==std::string::npos)splitPos=fileNameSnippet.length();
    
    char configFile[99];
    if(splitPos!=fileNameSnippet.length()){
	    sprintf(configFile, "Configuration_%s.conf", fileNameSnippet.substr(0,fileNameSnippet.length()).c_str());
    	std::cout << "Superseding default setup with information from " << fileNameSnippet.substr(splitPos,fileNameSnippet.length()) << " " << configFile<< std::endl;
	    std::ifstream optionsFile(configFile, std::ifstream::in);
	    boost::program_options::store(boost::program_options::parse_config_file(optionsFile, desc), *programOptions_);
	    optionsFile.close();
    }


    sprintf(configFile, "Configuration_%s.conf", fileNameSnippet.substr(0,splitPos).c_str());
    std::ifstream optionsFile(configFile, std::ifstream::in);
    
    boost::program_options::store(boost::program_options::parse_config_file(optionsFile, desc), *programOptions_);
    boost::program_options::notify(*programOptions_);

    if (programOptions_->count("help")) {
      std::cout << desc << std::endl;
      exit(0);
    }
    /*
  if (programOptions_->count("compression")) {
      std::cout << "Compression level was set to " << programOptions_->operator[]("compression").as<int>() << "." << std::endl;
  }
  else {
    //std::cout << "Compression level was not set.\n";
  }
     */
  }
}

const boost::program_options::variables_map*
ProgramOptionsReader::GetProgramOptions(){
  return 0;
}


std::string ProgramOptionsReader::GetOptionReplaced(std::string whichParameter,std::string replaceHelper){
  std::string tempReturn = GetOption<std::string>(whichParameter);
  if(replaceHelper!=""){
    std::vector<std::string> vsPars;
    boost::split(vsPars, replaceHelper, boost::is_any_of("|"));
    assert(vsPars.size()==2);
    boost::replace_all(tempReturn,  vsPars.at(0), vsPars.at(1));
  }
  return tempReturn;
}
