#include "calibration.h"
#include "ProgramOptionsReader.h"

int main(int ac, char** av)
{
  ProgramOptionsReader por = ProgramOptionsReader(ac,av);

  TopMassCalibration *top = new TopMassCalibration();
  delete top;
}
