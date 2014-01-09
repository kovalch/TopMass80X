#include "skimmer.h"
#include "ProgramOptionsReader.h"

int main(int ac, char** av)
{
  ProgramOptionsReader por = ProgramOptionsReader(ac,av);

  Skimmer *skim = new Skimmer();
  delete skim;
}
