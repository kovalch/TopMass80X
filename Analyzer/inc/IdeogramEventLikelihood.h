#include "IdeogramCombLikelihood.h"

// This class is not used in the likelihood calculation
// but just for drawing illustrative ideogram plots

class IdeogramEventLikelihood {
public:
  double DoEval(const double *x, const double *p) const;

  void SetFunction(std::vector<IdeogramCombLikelihood*> eventFunction) {
    eventFunction_ = eventFunction;
  };

private:
  std::vector<IdeogramCombLikelihood*> eventFunction_;
};
