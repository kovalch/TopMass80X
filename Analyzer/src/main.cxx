#include "TopMass.h"
#include "XMLConfigReader.h"

#include <boost/program_options.hpp>

#include <iostream>

int main(int ac, char** av)
{
  // Declare the supported options.
  po::options_description desc("Allowed options");
  desc.add_options()
    ("help,h", "produce help message")
    ("method,m", po::value<std::string>()->default_value("Ideogram"),
      "Top mass measurement method\n"
      "  GenMatch: \tGaussian fit of correct permutations (MC only)\n"
      "  Ideogram: \tJoint likelihood fit of mt and JES given data sample"
      "  RooFit:   \tTemplate fit of mt (JES, fSig) given data sample using RooFit"
    )
    ("task,t", po::value<std::string>()->default_value("sm"),
      "Task to be done\n"
      "  sm: \tSingle measurement based on input parameters\n"
      "  pe: \tPerform pseudo-experiments, additional settings necessary\n"
      "  hc: \tPerform any hardcoded tasks in TopMass constructor"
    )
    ("input,i", po::value<std::string>()->default_value("writeFullHadTree_data_2011"),
      "Identifier of input file to be analyzed: Z2_F11_172_5_sig")
    ("bins,b", po::value<int>()->default_value(1), "Number of phasespace bins")
    ("weight,w", po::value<int>()->default_value(0), "Weight type -> think, document")
    ("mass,M", po::value<double>()->default_value(172.5), "Input top mass for pseudo-experiments")
    ("jes,J", po::value<double>()->default_value(1.0), "Input JES for pseudo-experiments")
    ("fsig,f", po::value<double>()->default_value(0.504), "Input signal fraction for pseudo-experiments")
    ("lumi,L", po::value<double>()->default_value(0.0), "Luminosity for each pseudo-experiment")
    ("number,N", po::value<int>()->default_value(10000), "Number of pseudo-experiments per job")
    ("walltime,W", po::value<double>()->default_value(10), "set walltime limit for pseudo-experiments in minutes")
    ("shape,S", po::value<double>()->default_value(1.0), "Background shape scaling factor for gamma")
    ("permu,P", po::value<double>()->default_value(0.0), "Change permutation fractions by: fUN-P, fWP+0.5P, fCP+0.5P")
    ("fastsim,F", po::value<int>()->default_value(0), "use additional calibration for FastSim")
  ;

  po::variables_map vm;
  po::store(po::parse_command_line(ac, av, desc), vm);
  po::notify(vm);    

  if (vm.count("help")) {
      std::cout << desc << std::endl;
      return 1;
  }

  if (vm.count("compression")) {
      std::cout << "Compression level was set to " << vm["compression"].as<int>() << "." << std::endl;
  }
  else {
    //std::cout << "Compression level was not set.\n";
  }
  
  XMLConfigReader xmlConfig = XMLConfigReader();

  TopMass* top = new TopMass(vm);
  delete top;
}
