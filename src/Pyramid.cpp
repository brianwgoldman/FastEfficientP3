// Brian Goldman

// Implemention of the Parameter-less Population Pyramid
// Full description given in our publication:
// "Fast and Efficient Black Box Optimization using the Parameter-less Population Pyramid"
// by B. W. Goldman and W. F. Punch

#include "Pyramid.h"
#include <sstream>
using std::endl;

// Applies crossover between the passed in solution as each level
// of the pyramid
void Pyramid::climb(vector<bool> & solution, float & fitness) {
  // attempts to add this solution to the base level
  add_unique(solution, 0);
  for (size_t level = 0; level < pops.size(); level++) {
    float prev = fitness;
    // Use population clusters and population solutions to make new solution
    pops[level].improve(rand, solution, fitness, cross_counter);
    // add it to the next level if its a strict fitness improvement,
    // or configured to always add solutions
    if (not only_add_improvements or prev < fitness) {
      add_unique(solution, level + 1);
    }
  }
}

// Add the solution to the specified level of the population, only if it is
// unique.  Returns true if added.
bool Pyramid::add_unique(const vector<bool> & solution, size_t level) {
  if (seen.find(solution) == seen.end()) {
    if (pops.size() == level) {
      if (config.get<int>("wait_until_k_modeled"))
        if (pops.size() > 0 and not pops.back().k_modeled()) {
          return false;
        }
      // Create new levels as necessary
      pops.push_back(Population(config, level));
    }
    // Add the solution and rebuild the tree
    pops[level].add(solution);
    pops[level].rebuild_tree(rand);
    // track that this solution is now in the pyramid
    if (config.get<int>("prevent_duplicates")) {
      seen.insert(solution);
    }
    return true;
  }
  return false;
}

// Performs a full iteration of P3
bool Pyramid::iterate() {
  // generate a random solution
  restarts++;
  vector<bool> solution = rand_vector(rand, length);
  float fitness = local_counter->evaluate(solution);
  // perform a local search hill climber
  hill_climber(rand, solution, fitness, local_counter);
  // perform crossover with each level of the pyramid
  climb(solution, fitness);
  // P3 never "converges"
  return true;
}

string Pyramid::finalize() {
  std::ostringstream out;
  out << "# Restarts: "
      << restarts
      << " Hill: "
      // Convert back to Middle_Layer pointers to access the counters
      << std::static_pointer_cast<Middle_Layer>(local_counter)->evaluations
      << " Cross: "
      << std::static_pointer_cast<Middle_Layer>(cross_counter)->evaluations
      << endl;
  // output column headers
  out << "Size\tSuccesses\tTies\tFailures\tFitness\tDonation_Attempts\tDonation_Failures"
      << endl;
  for (const auto& pop : pops) {
    float total = 0;
    for (const auto& solution : pop.solutions) {
      total += evaluator->evaluate(solution);
    }
    out << pop.solutions.size() << "\t" << pop.successes << "\t" << pop.ties
        << "\t" << pop.failures << "\t" << total / pop.solutions.size() << "\t"
        << pop.donation_attempts << "\t" << pop.donation_failures << endl;
  }
  return out.str();
}
