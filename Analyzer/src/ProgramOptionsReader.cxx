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
        ("method,m", boost::program_options::value<std::string>()->default_value("Ideogram"),
            "Top mass measurement method\n"
            "  GenMatch: \tGaussian fit of correct permutations (MC only)\n"
            "  Ideogram: \tJoint likelihood fit of mt and JES given data sample\n"
            "  RooFit:   \tTemplate fit of mt (JES, fSig) given data sample using RooFit"
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
          "  alljets: \tallhadronic"
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
        ("fsig,f", boost::program_options::value<double>()->default_value(0.460272275), "Input signal fraction for pseudo-experiments")
        ("lumi,L", boost::program_options::value<double>()->default_value(0.0), "Luminosity for each pseudo-experiment")
        ("number,N", boost::program_options::value<int>()->default_value(10000), "Number of pseudo-experiments per job")
        ("walltime,W", boost::program_options::value<double>()->default_value(10), "set walltime limit for pseudo-experiments in minutes")
        ("shape,S", boost::program_options::value<double>()->default_value(1.0), "Background shape scaling factor for gamma")
        ("permu,P", boost::program_options::value<double>()->default_value(0.0), "Change permutation fractions by: fUN-P, fWP+0.5P, fCP+0.5P")
        ("fastsim,F", boost::program_options::value<int>()->default_value(0), "use additional calibration for FastSim")
        ("preliminary,p", boost::program_options::value<int>()->default_value(1), "use \"Preliminary\" label for plots")
        ("cmsenergy,e", boost::program_options::value<int>()->default_value(7), "cms energy to be used (for example for plots)")
        ("pullWidth", boost::program_options::value<double>()->default_value(1.0), "pull width correction factor")
        ("analysisConfig.selection", boost::program_options::value<std::string>())
        ("analysisConfig.selectionCP", boost::program_options::value<std::string>())
        ("analysisConfig.selectionWP", boost::program_options::value<std::string>())
        ("analysisConfig.selectionUN", boost::program_options::value<std::string>())
        ("analysisConfig.samplePath", boost::program_options::value<std::string>())
        ("analysisConfig.var1", boost::program_options::value<std::string>())
        ("analysisConfig.var2", boost::program_options::value<std::string>())
        ("analysisConfig.var3", boost::program_options::value<std::string>())
        ("analysisConfig.var4", boost::program_options::value<std::string>())
        ("analysisConfig.maxPermutations", boost::program_options::value<int>()->default_value(666666))
        ("analysisConfig.activeBranches", boost::program_options::value<std::string>()->default_value("*"))
        ("templates.fSig", boost::program_options::value<double>()->default_value(0.0))
        ("templates.fCP", boost::program_options::value<double>()->default_value(0.0))
        ("templates.fWP", boost::program_options::value<double>()->default_value(0.0))
        ("templates.fUN", boost::program_options::value<double>()->default_value(0.0))
        ("templates.parsCP", boost::program_options::value<std::string>()->default_value(""))
        ("templates.parsWP", boost::program_options::value<std::string>()->default_value(""))
        ("templates.parsUN", boost::program_options::value<std::string>()->default_value(""))
        ("templates.parsCPJES", boost::program_options::value<std::string>()->default_value(""))
        ("templates.parsWPJES", boost::program_options::value<std::string>()->default_value(""))
        ("templates.parsUNJES", boost::program_options::value<std::string>()->default_value(""))
        ("templates.parsBKG", boost::program_options::value<std::string>()->default_value(""))
        ("templates.parsBKGJES", boost::program_options::value<std::string>()->default_value(""))
        ("calibration.massOffset", boost::program_options::value<std::string>()->default_value("0.0"))
        ("calibration.massSlopeMass", boost::program_options::value<std::string>()->default_value("0.0"))
        ("calibration.massSlopeJES", boost::program_options::value<std::string>()->default_value("0.0"))
        ("calibration.massSlopeMassJES", boost::program_options::value<std::string>()->default_value("0.0"))
        ("calibration.jesOffset", boost::program_options::value<std::string>()->default_value("0.0"))
        ("calibration.jesSlopeMass", boost::program_options::value<std::string>()->default_value("0.0"))
        ("calibration.jesSlopeJES", boost::program_options::value<std::string>()->default_value("0.0"))
        ("calibration.jesSlopeMassJES", boost::program_options::value<std::string>()->default_value("0.0"))
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
    else {
      std::cerr << "Stopping analysis! Specified decay channel *" << channel << "* not known!" << std::endl;
      return;
    }
    char configFile[99];
    sprintf(configFile, "Configuration_%s.conf", fileNameSnippet.c_str());
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
