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
#include <boost/algorithm/string/replace.hpp>
#include <boost/type_traits/is_same.hpp>

class ProgramOptionsReader {
public:
  ProgramOptionsReader(int ac, char** av);
  virtual ~ProgramOptionsReader();

  static const boost::program_options::variables_map* GetProgramOptions();
  template <class T>
  static T GetOption(std::string whichParameter);
  static std::string GetOptionReplaced(std::string whichParameter,std::string replaceHelper="");

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

#endif /* PROGRAMOPTIONSREADER_H_ */
