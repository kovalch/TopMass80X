#include "TopMass.h"
#include "ProgramOptionsReader.h"

int main(int ac, char** av)
{
  ProgramOptionsReader por = ProgramOptionsReader(ac,av);

  TopMass* top = new TopMass();
  delete top;
}
