#include "HBOA.h"
#include <sstream>
using namespace std;

Factorizer::Factorizer(size_t n)
    : factors(n + 1, 1) {
  for (size_t i = 2; i <= sqrt(n); i++) {
    // if it is still prime, sieve
    if (factors[i] == 1) {
      for (size_t j = i + i; j <= n; j += i) {
        factors[j] = i;
      }
    }
  }
}

void Factorizer::operator()(vector<size_t>& buckets) const {
  for (size_t i = buckets.size() - 1; i >= 2; i--) {
    if (buckets[i]) {
      // If the number is not a prime
      if (factors[i] != 1) {
        buckets[factors[i]] += buckets[i];
        buckets[i / factors[i]] += buckets[i];
        buckets[i] = 0;
      }
    }
  }
}

void Factorial_Fraction::mul_by_factorial(size_t n) {
  safe_size(n);
  for (size_t i = 2; i <= n; i++) {
    upper[i]++;
  }
}

void Factorial_Fraction::div_by_factorial(size_t n) {
  safe_size(n);
  for (size_t i = 2; i <= n; i++) {
    lower[i]++;
  }
}

void Factorial_Fraction::simplify(const Factorizer& factorizer) {
  factorizer(upper);
  factorizer(lower);
  for (size_t i = 2; i < upper.size(); i++) {
    size_t extra = min(upper[i], lower[i]);
    upper[i] -= extra;
    lower[i] -= extra;
  }
}

Factorial_Fraction::operator double() const {
  double result = 1;
  for (size_t i = 2; i < upper.size(); i++) {
    result *= pow(i, upper[i]);
    result /= pow(i, lower[i]);
  }
  return result;
}

std::ostream& operator<<(std::ostream& out, const Factorial_Fraction& rhs) {
  out << "(1";
  for (size_t i = 2; i < rhs.upper.size(); i++) {
    if (rhs.upper[i]) {
      out << " * " << i << "^" << rhs.upper[i];
    }
  }
  out << ")/(1";
  for (size_t i = 2; i < rhs.lower.size(); i++) {
    if (rhs.lower[i]) {
      out << " * " << i << "^" << rhs.lower[i];
    }
  }
  out << ")";
  return out;
}

void Factorial_Fraction::safe_size(size_t n) {
  if (n >= upper.size()) {
    upper.resize(n + 1, 0);
    lower.resize(n + 1, 0);
  }
}

Factorial_Fraction Factorial_Fraction::operator*(
    const Factorial_Fraction& rhs) const {
  Factorial_Fraction result(*this);
  result.safe_size(rhs.upper.size());
  for (size_t i = 2; i < rhs.upper.size(); i++) {
    result.upper[i] += rhs.upper[i];
    result.lower[i] += rhs.lower[i];
  }
  return result;
}

Factorial_Fraction Factorial_Fraction::operator/(
    const Factorial_Fraction& rhs) const {
  Factorial_Fraction result(*this);
  result.safe_size(rhs.upper.size());
  for (size_t i = 2; i < rhs.upper.size(); i++) {
    result.upper[i] += rhs.lower[i];
    result.lower[i] += rhs.upper[i];
  }
  return result;
}

Bayesian_Tree::Bayesian_Tree(const Bayesian_Tree& rhs)
    : index(rhs.index),
      is_leaf(rhs.is_leaf),
      solutions(rhs.solutions),
      left(nullptr),
      right(nullptr),
      counts(rhs.counts) {
  if (not is_leaf) {
    left = new Bayesian_Tree(*rhs.left);
    right = new Bayesian_Tree(*rhs.right);
  }
}

Bayesian_Tree::~Bayesian_Tree() {
  delete left;
  delete right;
}

void swap(Bayesian_Tree& lhs, Bayesian_Tree& rhs) {
  swap(lhs.index, rhs.index);
  swap(lhs.is_leaf, rhs.is_leaf);
  swap(lhs.solutions, rhs.solutions);
  swap(lhs.left, rhs.left);
  swap(lhs.right, rhs.right);
  swap(lhs.counts, rhs.counts);
}

Bayesian_Tree& Bayesian_Tree::operator=(Bayesian_Tree rhs) {
  swap(*this, rhs);
  return *this;
}

Bayesian_Tree::operator string() const {
  if (is_leaf) {
    return to_string(float(counts[1]) / solutions.size());
  }
  string zero = static_cast<string>(*left);
  string one = static_cast<string>(*right);
  ostringstream out;
  out << "(" << index << ", " << zero << ", " << one << ")";
  return out.str();
}

ostream& operator<<(ostream& out, const Bayesian_Tree& tree) {
  out << static_cast<string>(tree);
  return out;
}

void Bayesian_Tree::add_solution(const vector<bool>& solution) {
  if (is_leaf) {
    solutions.push_back(solution);
    counts[solution[index]]++;
  } else if (solution[index]) {
    right->add_solution(solution);
  } else {
    left->add_solution(solution);
  }
}

void Bayesian_Tree::split(size_t new_index) {
  if (not is_leaf) {
    return;
  }
  is_leaf = false;
  left = new Bayesian_Tree(index);
  right = new Bayesian_Tree(index);
  index = new_index;
  for (const auto& solution : solutions) {
    if (solution[index]) {
      right->add_solution(solution);
    } else {
      left->add_solution(solution);
    }
  }
  solutions.clear();
}

void Bayesian_Tree::join() {
  if (is_leaf) {
    return;
  }
  left->join();
  right->join();
  solutions = left->solutions;
  solutions.insert(solutions.end(), right->solutions.begin(),
                   right->solutions.end());
  index = left->index;
  counts[0] = left->counts[0] + right->counts[0];
  counts[1] = left->counts[1] + right->counts[1];
  delete right;
  delete left;
  right = left = nullptr;
  is_leaf = true;
}

Factorial_Fraction Bayesian_Tree::bde_fraction() {
  if (is_leaf) {
    Factorial_Fraction result;
    result.div_by_factorial(solutions.size() + 1);
    // how many solutions have the index bit set to 0 or 1
    result.mul_by_factorial(counts[0]);
    result.mul_by_factorial(counts[1]);
    return result;
  } else {
    return left->bde_fraction() * right->bde_fraction();
  }
}

Factorial_Fraction Bayesian_Tree::splitless_bde(size_t new_index) {
  array<size_t, 4> splitcount = { 0, 0, 0, 0 };
  for (const auto& solution : solutions) {
    splitcount[(solution[new_index] << 1) | solution[index]]++;
  }
  Factorial_Fraction result;
  // left child
  result.div_by_factorial(splitcount[0] + splitcount[1] + 1);
  // right child
  result.div_by_factorial(splitcount[2] + splitcount[3] + 1);
  // individual pieces
  for (const auto& count : splitcount) {
    result.mul_by_factorial(count);
  }
  return result;
}

void Bayesian_Tree::set_bit(Random& rand, vector<bool>& solution) {
  if (is_leaf) {
    double probability = 0.5;
    if (solutions.size()) {
      probability = static_cast<double>(counts[1]) / solutions.size();
    }
    bernoulli_distribution dist(probability);
    solution[index] = dist(rand);
  } else {
    if (solution[index]) {
      right->set_bit(rand, solution);
    } else {
      left->set_bit(rand, solution);
    }
  }
}

Bayesian_Forest::Bayesian_Forest(size_t n)
    : trees(0),
      ordering(n),
      added(0) {
  for (size_t i = 0; i < n; i++) {
    trees.emplace_back(i);
  }
  // ordering
  iota(ordering.begin(), ordering.end(), 0);
}

void Bayesian_Forest::add_solution(const vector<bool>& solution) {
  added += 1;
  for (auto& tree : trees) {
    tree.add_solution(solution);
  }
}

void Bayesian_Forest::build_forest() {
  // Initially all trees have themselves as both prev and post
  prev.resize(trees.size());
  post.resize(trees.size());
  for (size_t i = 0; i < trees.size(); i++) {
    prev[i] = {i};
    post[i] = {i};
  }
  vector<tuple<double, Bayesian_Tree*, size_t>> options;
  // maximum factor you are going to need is 1 more than all solutions.
  Factorizer factors(added + 1);
  for (auto& tree : trees) {
    // ensure it starts as a leaf.
    tree.join();
    build_options(tree, factors, options);
  }

  // ensures minimum quality for the split
  double threshold = pow(2, -0.5 * log2(added));
  filter(options, threshold);
  while (options.size()) {
    auto best = *min_element(options.begin(), options.end());
    Bayesian_Tree& splitting = *get<1>(best);
    size_t new_index = get<2>(best);
    // update dependency information
    auto* new_req = &(post[splitting.index]);
    for (size_t upstream : prev[new_index]) {
      post[upstream].insert(new_req->begin(), new_req->end());
    }

    new_req = &(prev[new_index]);
    for (size_t down_stream : post[splitting.index]) {
      prev[down_stream].insert(new_req->begin(), new_req->end());
    }

    splitting.split(new_index);

    // build new options for the two new trees
    build_options(*splitting.left, factors, options);
    build_options(*splitting.right, factors, options);

    // clear out any invalid splits
    filter(options, threshold);
  }

  // determine what order bits can be generated to ensure all dependencies
  // are met
  sort(ordering.begin(), ordering.end(),
       [this](size_t a, size_t b) {return prev[a].size() < prev[b].size();});
}

void Bayesian_Forest::generate(Random& rand, vector<bool> & solution) {
  for (const auto& index : ordering) {
    trees[index].set_bit(rand, solution);
  }
}

std::ostream& operator<<(std::ostream& out, const Bayesian_Forest& forest) {
  for (size_t i = 0; i < forest.trees.size(); i++) {
    out << i << ": " << forest.trees[i] << endl;
  }
  return out;
}

void Bayesian_Forest::build_options(
    Bayesian_Tree& tree, const Factorizer& factors,
    vector<std::tuple<double, Bayesian_Tree*, size_t>>& options) const {
  Factorial_Fraction initial = tree.bde_fraction();
  const auto & invalid = post[tree.index];
  for (size_t i = 0; i < trees.size(); i++) {
    // if i is not invalid
    if (invalid.count(i) == 0) {
      Factorial_Fraction with_split = tree.splitless_bde(i);
      Factorial_Fraction result = initial / with_split;
      result.simplify(factors);
      options.emplace_back(result, &tree, i);
    }
  }
}

void Bayesian_Forest::filter(
    vector<std::tuple<double, Bayesian_Tree*, size_t>> & options,
    double threshold) {
  vector<tuple<double, Bayesian_Tree*, size_t>> new_options;
  for (const auto& option : options) {
    // only keep options that use a node which is still a leaf, that can
    // be added without causing a circular dependency, and is better than
    // the limiting threshold.
    if (get<1>(option)->is_leaf
        and post[get<1>(option)->index].count(get<2>(option)) == 0
        and get<0>(option) < threshold) {
      new_options.push_back(option);
    }
  }
  swap(options, new_options);
}

HBOA::HBOA(Random& _rand, Evaluator& _evaluator, Configuration& _config)
    : Optimizer(_rand, _evaluator, _config) {
  size_t pop_size = config.get<int>("pop_size");

  // optionally applies a hill climber to the initial solutions
  hc = config.get<hill_climb::pointer>("hill_climber");

  float fitness;
  for (size_t i = 0; i < pop_size; i++) {
    // create and evaluate solutions
    auto solution = rand_vector(rand, length);
    fitness = evaluator.evaluate(solution);
    // Apply hill climber if configured to do so
    hc(rand, solution, fitness, evaluator);
    solutions.push_back(solution);
    fitnesses.push_back(fitness);
  }
  selection.resize(pop_size);
  iota(selection.begin(), selection.end(), 0);
  rtr_size = min(length, pop_size / 20);
  if (rtr_size == 0) {
    rtr_size = 1;
  }
}

bool HBOA::iterate() {
  Bayesian_Forest model(length);
  // binary tournament
  uniform_int_distribution<int> choose(0, solutions.size() - 1);
  int x, y;
  for (size_t i = 0; i < selection.size(); i++) {
    x = choose(rand);
    y = choose(rand);
    // compete adjacent solutions in the solutions vector
    if (fitnesses[x] < fitnesses[y]) {
      model.add_solution(solutions[y]);
    } else {
      model.add_solution(solutions[x]);
    }
  }
  model.build_forest();

  // storage for the newly generated solution
  float fitness;
  vector<bool> solution(length);
  // track if anything in the current population has been replaced
  bool replaced = false;
  for (size_t i = 0; i < solutions.size(); i++) {
    model.generate(rand, solution);
    fitness = evaluator.evaluate(solution);
    hc(rand, solution, fitness, evaluator);
    int choice = rtr_nearest(solution);
    if (fitnesses[choice] < fitness) {
      solutions[choice] = solution;
      fitnesses[choice] = fitness;
      replaced = true;
    }
  }
  return replaced;
}

int HBOA::rtr_nearest(const vector<bool>& solution) {
  int index;
  int choice = -1;
  int distance = length + 1;
  for (size_t i = 0; i < rtr_size; i++) {
    // only look at unused indices.
    index = std::uniform_int_distribution<int>(i, selection.size() - 1)(rand);
    swap(selection[i], selection[index]);
    int new_dist = hamming_distance(solutions[selection[i]], solution);
    if (new_dist < distance) {
      distance = new_dist;
      choice = selection[i];
    }
  }
  return choice;
}
