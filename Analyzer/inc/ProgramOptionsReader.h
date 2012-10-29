/*
 * ProgramOptionsReader.h
 *
 *  Created on: Oct 29, 2012
 *      Author: eschliec
 */

#ifndef PROGRAMOPTIONSREADER_H_
#define PROGRAMOPTIONSREADER_H_

#include <assert.h>
#include <iostream>
#include <string>

#include <boost/program_options.hpp>

class ProgramOptionsReader {
public:
  ProgramOptionsReader(int ac, char** av);
  virtual ~ProgramOptionsReader();

  static const boost::program_options::variables_map* GetProgramOptions();
  template <class T>
  static T GetOption(std::string whichParameter);

private:
  static boost::program_options::variables_map* programOptions_;
  void ReadProgramOptions(int ac, char** av);
};

template <class T>
T
ProgramOptionsReader::GetOption(std::string whichParameter){
  if (programOptions_->count(whichParameter)){
    return programOptions_->operator[](whichParameter).as<T>();
  }
  else{
    std::cerr << "Program option *" << whichParameter << "* not found!" << std::endl;
    assert(0);
  }
}

#endif /* PROGRAMOPTIONSREADER_H_ */
