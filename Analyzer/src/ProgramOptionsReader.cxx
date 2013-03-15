/*
 * ProgramOptionsReader.cxx
 *
 *  Created on: Oct 29, 2012
 *      Author: eschliec
 */

#include "ProgramOptionsReader.h"

ProgramOptionsReader::ProgramOptionsReader(int ac, char** av) {
  vm_        = new CommandLineOptionsReader(ac,av);
  xmlConfig_ = new XMLConfigReader();
}

CommandLineOptionsReader* ProgramOptionsReader::vm_(0);
XMLConfigReader* ProgramOptionsReader::xmlConfig_(0);

ProgramOptionsReader::~ProgramOptionsReader() {
}
