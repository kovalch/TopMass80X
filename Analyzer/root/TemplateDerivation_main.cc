#include "ProgramOptionsReader.h"
#include "TemplateDerivation.h"

int main(int ac, char** av) {
  ProgramOptionsReader por = ProgramOptionsReader(ac, av);
  std::cout<<"TemplateDerivation_Main DEBUG before TemplateDerivation Constructor"<<std::endl;
  TemplateDerivation top;
  top.run();
}
