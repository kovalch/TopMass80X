#ifndef TOPMASS_H
#define TOPMASS_H

#include <string>
#include <vector>

class TopMass {
  private:
    //bool fexists(const char *filename);
    void WriteEnsembleTest(const std::vector<float>& vBinning);
      
    const std::string fBinning_, fTask_;

  public:
    TopMass();
};

#endif
