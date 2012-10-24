/*
 * XMLConfigReader.h
 *
 *  Created on: Oct 24, 2012
 *      Author: eschliec
 */

#ifndef XMLCONFIGREADER_H_
#define XMLCONFIGREADER_H_

#include "TString.h"

#include "tinyxml2.h"

namespace txml = tinyxml2;

class XMLConfigReader {
public:
  XMLConfigReader();
  ~XMLConfigReader();

  static const txml::XMLNode* GetConfig();
  static TString GetParameter(TString whichParameter);

private:
  static txml::XMLDocument* _config;
  void ReadConfigFromXMLFile();

};

#endif /* XMLCONFIGREADER_H_ */
