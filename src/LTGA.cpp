// Brian Goldman

// Implements the Linkage Tree Genetic Algorithm
// Allows for features not used in the paper, like
// removing the binary tournament, performing hill climbing, etc

#include "LTGA.h"
#include "MiddleLayer.h"
using std::endl;

// Constructs and evaluates the initial population of solutions
LTGA::LTGA(Random& _rand, shared_ptr<Evaluator> _evaluator,
           Configuration& _config)
    : Optimizer(_rand, _evaluator, _config),
      pop(_config),
      local_counter(new Middle_Layer(config, _evaluator, false)),
      cross_counter(new Middle_Layer(config, _evaluator, false)) {
  pop_size = config.get<int>("pop_size");
  disable_binary_insert = config.get<int>("binary_insert") != 1;

  // optionally applies a hill climber to the initial solutions
  hc = config.get<hill_climb::pointer>("hill_climber");

  float fitness;
  vector<vector<bool>> solutions;
  for (size_t i = 0; i < pop_size; i++) {
    // create and evaluate solutions
    auto solution = rand_vector(rand, length);
    fitness = local_counter->evaluate(solution);
    // Apply hill climber if configured to do so
    hc(rand, solution, fitness, local_counter);
    solutions.push_back(solution);
    fitnesses[solution] = fitness;
  }
  // Inserts the solutions into the population using a binary tournament
  binary_insert(rand, solutions, pop);
}

// Determines which solutions in the population should be used when
// calculating the pairwise entropy.
void LTGA::binary_insert(Random& rand, vector<vector<bool>> & solutions,
                         Population& next_pop) {
  std::shuffle(solutions.begin(), solutions.end(), rand);
  for (size_t i = 0; i + 1 < solutions.size(); i += 2) {
    // compete adjacent solutions in the solutions vector
    if (fitnesses[solutions[i]] < fitnesses[solutions[i + 1]]) {
      next_pop.add(solutions[i + 1]);
      // if disable_binary_insert is true, the worse solution
      // is still used when calculating entropy.
      next_pop.add(solutions[i], disable_binary_insert);
    } else {
      next_pop.add(solutions[i]);
      next_pop.add(solutions[i + 1], disable_binary_insert);
    }
  }
}

// Performs a single generation of evolution
void LTGA::generation() {
  // constructs the required crossovers
  pop.rebuild_tree(rand);
  float fitness;
  Population next_generation(config);
  // temporary storage of the generated solutions
  vector<vector<bool>> solutions;
  for (auto solution : pop.solutions) {
    fitness = fitnesses[solution];
    // uses crossover to improve a copy of the solution
    pop.improve(rand, solution, fitness, cross_counter);
    fitnesses[solution] = fitness;
    solutions.push_back(solution);
  }
  // creates the next_generation from the offspring solutions
  binary_insert(rand, solutions, next_generation);

  // data recording
  float total = 0;
  for (const auto& solution : pop.solutions) {
    total += fitnesses[solution];
  }
  metadata << pop.solutions.size() << "\t" << pop.successes << "\t" << pop.ties
           << "\t" << pop.failures << "\t" << total / pop.solutions.size()
           << "\t" << pop.donation_attempts << "\t" << pop.donation_failures
           << endl;

  pop = next_generation;
}

// Implements the Optimization interface, returns true as long as
// LTGA thinks it can still improve
bool LTGA::iterate() {
  generation();
  // keeps track of if the new population actually has different solutions
  // than the previous generation
  decltype(pop_set) new_set(pop.solutions.begin(), pop.solutions.end());
  if (new_set == pop_set) {
    // all organisms identical after two generation, so stop
    return false;
  }
  pop_set = new_set;
  return true;
}

// Extract end of run logging information to be printed to meta files.
string LTGA::finalize() {
  std::ostringstream out;
  out << "# Restarts: "
      << pop_size
      << " Hill: "
      // Convert back to Middle_Layer pointers to access the counters
      << std::static_pointer_cast<Middle_Layer>(local_counter)->evaluations
      << " Cross: "
      << std::static_pointer_cast<Middle_Layer>(cross_counter)->evaluations
      << endl;
  // output column headers
  out << "Size\tSuccesses\tTies\tFailures\tFitness\tDonation_Attempts\tDonation_Failures"
      << endl;
  out << metadata.str();
  return out.str();
}
