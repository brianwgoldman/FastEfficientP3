// Brian Goldman

// Implementation of the Parameter-less GA's
// method for removing the population size parameter
// Designed to match the variant in the paper:
// "Parameter-Less Hierarchical BOA"
// by M. Pelikan and T. K. Lin
// Takes a configurable sub_optimizer which
// is kept in a collection of populations, each using a
// different population size.

#ifndef POPLESS_H_
#define POPLESS_H_

#include "Optimizer.h"

// Inherits and implements the optimizer
class Popless : public Optimizer {
 public:
  Popless(Random& _rand, shared_ptr<Evaluator> _evaluator,
          Configuration& _config);
  bool iterate() override;
  create_optimizer(Popless);
 private:
  // collections of populations
  vector<shared_ptr<Optimizer>> pops;
  // records how many iterations each population has performed
  vector<size_t> generations;
  // keeps track of the size of the current population
  size_t largest_size;
  // which population gets to perform the next iteration
  size_t up_next;
};

#endif /* POPLESS_H_ */
