#include "TopMass.h"
#include "XMLConfigReader.h"
#include "ProgramOptionsReader.h"

#include <iostream>

int main(int ac, char** av)
{
  ProgramOptionsReader vm = ProgramOptionsReader(ac,av);
  XMLConfigReader xmlConfig = XMLConfigReader();

  TopMass* top = new TopMass();
  delete top;
}
