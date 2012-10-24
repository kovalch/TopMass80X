/*
 * XMLConfigReader.cxx
 *
 *  Created on: Oct 24, 2012
 *      Author: eschliec
 */

#include <assert.h>
#include <iostream>

#include "XMLConfigReader.h"


XMLConfigReader::XMLConfigReader() {
  ReadConfigFromXMLFile();
}

XMLConfigReader::~XMLConfigReader() {
  delete _config;
}

txml::XMLDocument* XMLConfigReader::_config(0);

void
XMLConfigReader::ReadConfigFromXMLFile(){
  if(!_config){
    _config = new txml::XMLDocument();
    TString xmlFilePath = "/afs/naf.desy.de/group/cms/scratch/eschliec/TopMass_hg_devel/Analyzer/Configuration_alljets.xml";
    int errorID = _config->LoadFile(xmlFilePath);
    if(errorID) {
      std::cerr << "Parsing of XML file (" << xmlFilePath << ") failed with error " << errorID << "!" << std::endl;
      assert(!errorID);
    }
    txml::XMLElement *analysisConfiguration = 0;

    analysisConfiguration = _config->FirstChildElement("analysisConfig");
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
  return _config->FirstChildElement("analysisConfig");
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
