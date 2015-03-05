// Implements the Parameter-less GA method of population
// parameter removal.

#include "Popless.h"
#include "OptimizationCollection.h"

Popless::Popless(Random& _rand, shared_ptr<Evaluator> _evaluator,
                 Configuration& _config)
    : Optimizer(_rand, _evaluator, _config) {
  largest_size = config.get<int>("pop_size");
  // copy the configuration and set up the sub_optimizer
  Configuration lesser(config);
  auto optimizer_method = config.get<optimize::pointer>("sub_optimizer");
  pops.push_back(optimizer_method(rand, evaluator, lesser));
  generations.push_back(0);
  up_next = 0;
}

bool Popless::iterate() {
  // performs an iteration of one population
  pops[up_next]->iterate();
  generations[up_next]++;
  if (generations[up_next] % 2 == 0) {
    up_next++;
    // creates new population as necessary
    if (up_next >= pops.size()) {
      // set up the new population size
      Configuration lesser(config);
      largest_size *= 2;
      lesser.set("pop_size", largest_size);
      auto optimizer_method = config.get<optimize::pointer>("sub_optimizer");
      pops.push_back(optimizer_method(rand, evaluator, lesser));
      generations.push_back(0);
    }
  } else {
    up_next = 0;
  }
  // Effectively terminates populations after they have performed
  // too many iterations
  while (generations[up_next] >= length) {
    up_next++;
  }
  return true;
}
