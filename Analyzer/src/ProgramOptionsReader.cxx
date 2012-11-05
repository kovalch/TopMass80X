/*
 * ProgramOptionsReader.cxx
 *
 *  Created on: Oct 29, 2012
 *      Author: eschliec
 */

#include "ProgramOptionsReader.h"

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
    // Declare the supported options.
    boost::program_options::options_description desc("Allowed options");
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
            "  hc: \tPerform any hardcoded tasks in TopMass constructor"
            "  diff: \tDifferential top mass, specify binning option"
        )
        ("input,i", boost::program_options::value<std::string>()->default_value("writeFullHadTree_data_2011"),
            "Identifier of input file to be analyzed: Z2_F11_172_5_sig"
        )
        ("lepton,l", boost::program_options::value<std::string>()->default_value("all"),
          "Channel\n"
          "  electron \t\n"
          "  muon \t\n"
          "  all: \telectron+muon"
          "  none: \tallhadronic"
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
        ("fsig,f", boost::program_options::value<double>()->default_value(0.504), "Input signal fraction for pseudo-experiments")
        ("lumi,L", boost::program_options::value<double>()->default_value(0.0), "Luminosity for each pseudo-experiment")
        ("number,N", boost::program_options::value<int>()->default_value(10000), "Number of pseudo-experiments per job")
        ("walltime,W", boost::program_options::value<double>()->default_value(10), "set walltime limit for pseudo-experiments in minutes")
        ("shape,S", boost::program_options::value<double>()->default_value(1.0), "Background shape scaling factor for gamma")
        ("permu,P", boost::program_options::value<double>()->default_value(0.0), "Change permutation fractions by: fUN-P, fWP+0.5P, fCP+0.5P")
        ("fastsim,F", boost::program_options::value<int>()->default_value(0), "use additional calibration for FastSim")
        ("xmlconfig,X", boost::program_options::value<std::string>()->default_value("/afs/naf.desy.de/group/cms/scratch/eschliec/TopMass_hg_devel/Analyzer/Configuration_alljets.xml"),
            "path to the XML configuration file"
        )
        ;

    boost::program_options::store(boost::program_options::parse_command_line(ac, av, desc), *programOptions_);
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
