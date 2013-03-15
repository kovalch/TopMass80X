#include "TopMass.h"
#include "CommandLineOptionsReader.h"
#include "XMLConfigReader.h"

#include <iostream>

int main(int ac, char** av)
{
  CommandLineOptionsReader vm = CommandLineOptionsReader(ac,av);
  XMLConfigReader xmlConfig = XMLConfigReader();

  TopMass* top = new TopMass();
  delete top;
}
