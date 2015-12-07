#include "AddPDFweights.h"
#include "ProgramOptionsReader.h"

int main(int ac, char** av)
{
  ProgramOptionsReader por = ProgramOptionsReader(ac,av);

  AddPDFweights *adder = new AddPDFweights();
  delete adder;
}
