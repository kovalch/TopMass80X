/*
 * XMLConfigReader.cxx
 *
 *  Created on: Oct 24, 2012
 *      Author: eschliec
 */

#include <assert.h>
#include <iostream>

#include "ProgramOptionsReader.h"
#include "XMLConfigReader.h"

typedef ProgramOptionsReader po;

XMLConfigReader::XMLConfigReader() {
  ReadConfigFromXMLFile();
}

XMLConfigReader::~XMLConfigReader() {
  delete config_;
}

txml::XMLDocument* XMLConfigReader::config_(0);

void
XMLConfigReader::ReadConfigFromXMLFile(){
  if(!config_){
    config_ = new txml::XMLDocument();
    TString xmlFilePath = po::GetOption<std::string>("xmlconfig");
    int errorID = config_->LoadFile(xmlFilePath);
    if(errorID) {
      std::cerr << "Parsing of XML file (" << xmlFilePath << ") failed with error " << errorID << "!" << std::endl;
      assert(!errorID);
    }
    txml::XMLElement *analysisConfiguration = 0;

    analysisConfiguration = config_->FirstChildElement("analysisConfig");
    if(!analysisConfiguration){
      std::cerr << "No *analysisConfig* object contained in XMLFile:\n" << xmlFilePath << std::endl;
      assert(0);
    }
    if(analysisConfiguration->NoChildren()){
      std::cerr << "No configuration contained in *" << analysisConfiguration->Value() << "* object in XMLFile:\n" << xmlFilePath << std::endl;
      assert(0);
    }
  }
}

const txml::XMLNode*
XMLConfigReader::GetConfig(){
  return config_->FirstChildElement("analysisConfig");
}

TString
XMLConfigReader::GetParameter(TString whichParameter){
  const txml::XMLElement* element = GetConfig()->FirstChildElement(whichParameter);
  if(element) return TString(element->GetText());
  else{
    std::cerr << "Parameter *" << whichParameter << "* not found in XML configuration!" << std::endl;
    assert(0);
  }
  return TString("");
}
