#ifndef DATASAMPLE_H
#define DATASAMPLE_H

#include <vector>

class DataSample
{
public:

  struct SimpleEvent
  {
    struct Permutation
    {
      double topMass;
      double wMass;
      double prob;
      int bin;
    };

    std::vector<Permutation> permutations;
    int leptonFlavour;
    double weight;
  };
  DataSample() : nEvents(0), maxWeight(0) {}
  void Clear();
  void Fill(double topMass, double wMass, double prob, int leptonFlavour, double weight, int index, int bin);
  //template<class T> void Fill(double topMass, double wMass, double prob, double weight, T index) = delete;
  void AddEvent(const SimpleEvent& event);

  DataSample& operator+=(const DataSample& sample);

  int nEvents;
  double maxWeight;

  std::vector<SimpleEvent> events;
};

#endif
