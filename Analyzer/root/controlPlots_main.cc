#include "controlPlots.h"
#include "ProgramOptionsReader.h"

int main(int ac, char** av)
{
  ProgramOptionsReader por = ProgramOptionsReader(ac,av);

  TopMassControlPlots top = TopMassControlPlots();
}
