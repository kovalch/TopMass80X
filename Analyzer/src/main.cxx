#include "TopMass.h"
#include "ProgramOptionsReader.h"

#include <iostream>

int main(int ac, char** av)
{
  ProgramOptionsReader vm = ProgramOptionsReader(ac,av);

  TopMass* top = new TopMass();
  delete top;
}
