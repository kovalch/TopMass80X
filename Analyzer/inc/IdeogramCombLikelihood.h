#ifndef IDEOGRAMCOMBLIKELIHOOD_H
#define IDEOGRAMCOMBLIKELIHOOD_H

#include <vector>

class IdeogramCombLikelihood {
 public:
  IdeogramCombLikelihood();
  virtual ~IdeogramCombLikelihood() {};

  virtual double Evaluate(double *x, double *p) = 0;

  std::vector<double> readParameters(const char *whichParameter);
  void SetFixedParams(double p0 = 0, double p1 = 0, double p2 = 0, double p3 = 0, double p4 = 0, double p5 = 0, double p6 = 0, double p7 = 0, double p8 = 0);
  double GetFixedParam(int index);
  void SetActive(bool active);
  bool IsActive();

 protected:
  std::vector<double> fp_;
  bool useFixedParams_;
  bool active_;

};

#endif
