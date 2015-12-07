#include "ProgramOptionsReader.h"
#include "TemplateDerivation.h"

int main(int ac, char** av)
{
  ProgramOptionsReader por = ProgramOptionsReader(ac,av);

  TemplateDerivation *top = new TemplateDerivation();
  delete top;
}
