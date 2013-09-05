#include "TopMass.h"
#include "ProgramOptionsReader.h"

#include <iostream>
#include <csignal>

#include "TFile.h"
#include "TROOT.h"

void signalHandler( int signum )
{
  std::cout << "\nInterrupt signal (" << signum << ") received." << std::endl;

  // cleanup and close up stuff here
  // terminate program

  for(int i=0; i<gROOT->GetListOfFiles()->GetSize(); ++i)
    ((TFile*)gROOT->GetListOfFiles()->At(i))->Close();

  exit(signum);

}

int main(int ac, char** av)
{

  signal(SIGINT, signalHandler);

  ProgramOptionsReader por = ProgramOptionsReader(ac,av);

  TopMass* top = new TopMass();
  delete top;
}
