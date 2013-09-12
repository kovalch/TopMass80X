#ifndef IDEOGRAMCOMBLIKELIHOOD_H
#define IDEOGRAMCOMBLIKELIHOOD_H

class IdeogramCombLikelihood {
public:
  IdeogramCombLikelihood() {};
  virtual ~IdeogramCombLikelihood() {};

  virtual double Evaluate(double *x, double *p) = 0;
};

#endif
