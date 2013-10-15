#include "duplicateRemoval.h"
#include "ProgramOptionsReader.h"

int main(int ac, char** av)
{
  ProgramOptionsReader por = ProgramOptionsReader(ac,av);

  DuplicateRemover *remover = new DuplicateRemover();
  delete remover;
}
