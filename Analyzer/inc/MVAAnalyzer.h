#include "MassAnalyzer.h"
#include "TF1.h"

class MVAAnalyzer : public MassAnalyzer {
  private:
    TF1* gaus;
  
  public:
    MVAAnalyzer(const std::string& identifier, TTree* tree) : MassAnalyzer(identifier, tree), gaus(0) {};
    void Analyze(const std::string& cuts, int i, int j);

};
