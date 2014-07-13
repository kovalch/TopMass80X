#include "infRemoval.h"
#include "ProgramOptionsReader.h"

int main(int ac, char** av)
{
  ProgramOptionsReader por = ProgramOptionsReader(ac,av);

  InfRemover *remover = new InfRemover();
  delete remover;
}
