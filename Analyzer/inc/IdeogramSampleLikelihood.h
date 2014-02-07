#include "IdeogramCombLikelihoodAllJets.h"
#include "IdeogramCombLikelihoodLeptonJets.h"
#include "ProgramOptionsReader.h"

#include "Math/IFunction.h"
 
class IdeogramSampleLikelihood : public ROOT::Math::IBaseFunctionMultiDim {
public:
  double DoEval(const double *x) const;
  
  unsigned int NDim() const
   {
     return 2;
   }

  IdeogramSampleLikelihood* Clone() const {
    return new IdeogramSampleLikelihood();
  }

  void AddFunctions(std::vector<std::vector<IdeogramCombLikelihood*>> eventFunctions) {
    eventFunctions_ = eventFunctions;
  };

private:
  std::vector<std::vector<IdeogramCombLikelihood*>> eventFunctions_;
};
