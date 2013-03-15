/*
 * ProgramOptionsReader.h
 *
 *  Created on: Oct 29, 2012
 *      Author: eschliec
 */

#ifndef PROGRAMOPTIONSREADER_H_
#define PROGRAMOPTIONSREADER_H_

#include "CommandLineOptionsReader.h"
#include "XMLConfigReader.h"

class ProgramOptionsReader {
public:
  ProgramOptionsReader(int ac, char** av);
  virtual ~ProgramOptionsReader();

  template <class T>
  static T GetOption(std::string whichParameter);

private:
  static CommandLineOptionsReader* vm_;
  static XMLConfigReader* xmlConfig_;
};

template <class T>
T
ProgramOptionsReader::GetOption(std::string whichParameter){
  if(vm_->GetOption<T>(whichParameter)) return *vm_->GetOption<T>(whichParameter);
  else if(xmlConfig_->GetParameter(whichParameter)) return *xmlConfig_->GetParameter(whichParameter.c_str());
  else{
    std::cerr << "Parameter *" << whichParameter << "* neither found in command line nor in XML configuration!" << std::endl;
    assert(0);
  }
}

#endif /* PROGRAMOPTIONSREADER_H_ */
