#ifndef DATASAMPLE_H
#define DATASAMPLE_H

#include <vector>

class DataSample
{
public:

  void Clear();
  void Fill(double topMass, double wMass, double prob, double weight, unsigned char index);
  template<class T> void Fill(double topMass, double wMass, double prob, double weight, T index) = delete;

  DataSample& operator+=(const DataSample& sample);

  int nEvents;
  double maxWeight;

  std::vector<double> topMasses;
  std::vector<double> wMasses;
  std::vector<double> fitProbs;
  std::vector<double> weights;
  std::vector<unsigned char> indices;
};

#endif
