// Brian Goldman

// Implementation of the Hierarchical Bayesian Optimization Algorithm
// based on the publication:
// Pelikan, M. and Goldberg, D. (2006). Hierarchical bayesian optimization algorithm. In Scalable
// Optimization via Probabilistic Modeling, volume 33 of Studies in Computational Intelligence, pages
// 63–90. Springer Berlin Heidelberg.

#include "HBOA.h"
#include <sstream>
using namespace std;

// Multiply this fraction by n!
void Factorial_Fraction::mul_by_factorial(size_t n) {
  // ensures the factorial is large enough
  safe_size(n);
  for (size_t i = 2; i <= n; i++) {
    upper[i]++;
  }
}

// Divide this fraction by n!
void Factorial_Fraction::div_by_factorial(size_t n) {
  // ensures the factorial is large enough
  safe_size(n);
  for (size_t i = 2; i <= n; i++) {
    lower[i]++;
  }
}

// Trivial simplification to remove numbers which appear
// in both the top and the bottom of the fraction
void Factorial_Fraction::simplify() {
  for (size_t i = 2; i < upper.size(); i++) {
    // find the common value
    size_t extra = min(upper[i], lower[i]);
    upper[i] -= extra;
    lower[i] -= extra;
  }
}

// Calculate the real value from the fraction
Factorial_Fraction::operator double() const {
  double result = 1;
  // the current value to be added from 't'op and 'b'ottom
  size_t t = 2, b = 2;
  // copy number counts to temporary vectors
  vector<size_t> top(upper);
  vector<size_t> bot(lower);
  // keep integrating numbers from the top to the bottom
  while (t < top.size() and b < bot.size()) {
    // above 1 make smaller
    if (result > 1) {
      // if there are more 'b' numbers to integrate
      if (bot[b]) {
        result /= b;
        bot[b]--;
      } else {
        b++;
      }
    } else {
      // if there are more 't' numbers to integrate
      if (top[t]) {
        result *= t;
        top[t]--;
      } else {
        t++;
      }
    }
  }
  // finish any remaining top
  while (t < top.size()) {
    if (top[t]) {
      result *= t;
      top[t]--;
    } else {
      t++;
    }
  }
  // finish any remaining bottom
  while (b < bot.size()) {
    if (bot[b]) {
      result /= b;
      bot[b]--;
    } else {
      b++;
    }
  }
  return result;
}

// Utility to display the fraction for debugging purposes
std::ostream& operator<<(std::ostream& out, const Factorial_Fraction& rhs) {
  // Output the top of the fraction.
  out << "(1";
  for (size_t i = 2; i < rhs.upper.size(); i++) {
    if (rhs.upper[i]) {
      out << " * " << i << "^" << rhs.upper[i];
    }
  }
  // Output the bottom of the fraction
  out << ")/(1";
  for (size_t i = 2; i < rhs.lower.size(); i++) {
    if (rhs.lower[i]) {
      out << " * " << i << "^" << rhs.lower[i];
    }
  }
  out << ")";
  return out;
}

// Ensure the factorial is large enough to add numbers up to
// (and including) n
void Factorial_Fraction::safe_size(size_t n) {
  if (n >= upper.size()) {
    upper.resize(n + 1, 0);
    lower.resize(n + 1, 0);
  }
}

// Multiply two fractions
Factorial_Fraction Factorial_Fraction::operator*(
    const Factorial_Fraction& rhs) const {
  // copy the left side into the result
  Factorial_Fraction result(*this);
  result.safe_size(rhs.upper.size());
  // add the right side's counts to the result
  for (size_t i = 2; i < rhs.upper.size(); i++) {
    result.upper[i] += rhs.upper[i];
    result.lower[i] += rhs.lower[i];
  }
  return result;
}

// Divide two fraction
Factorial_Fraction Factorial_Fraction::operator/(
    const Factorial_Fraction& rhs) const {
  // copy the numerator into the result
  Factorial_Fraction result(*this);
  result.safe_size(rhs.upper.size());
  // the right side's numerator goes into the results denominator
  // and vice versa
  for (size_t i = 2; i < rhs.upper.size(); i++) {
    result.upper[i] += rhs.lower[i];
    result.lower[i] += rhs.upper[i];
  }
  return result;
}

// Recursive copy for a Bayesian Tree
Bayesian_Tree::Bayesian_Tree(const Bayesian_Tree& rhs)
    : index(rhs.index),
      is_leaf(rhs.is_leaf),
      solutions(rhs.solutions),
      left(nullptr),
      right(nullptr),
      counts(rhs.counts) {
  // Recursively create children trees
  if (not is_leaf) {
    left = new Bayesian_Tree(*rhs.left);
    right = new Bayesian_Tree(*rhs.right);
  }
}

// Cleans up dynamic memory
Bayesian_Tree::~Bayesian_Tree() {
  // recursion automatically stopped by delete nullptr
  delete left;
  delete right;
}

// Swap all member variables to swap two trees
void swap(Bayesian_Tree& lhs, Bayesian_Tree& rhs) {
  swap(lhs.index, rhs.index);
  swap(lhs.is_leaf, rhs.is_leaf);
  swap(lhs.solutions, rhs.solutions);
  swap(lhs.left, rhs.left);
  swap(lhs.right, rhs.right);
  swap(lhs.counts, rhs.counts);
}

// Use the copy and swap idiom to handle assignment
Bayesian_Tree& Bayesian_Tree::operator=(Bayesian_Tree rhs) {
  swap(*this, rhs);
  return *this;
}

// Utility function to convert a tree into a string
Bayesian_Tree::operator string() const {
  if (is_leaf) {
    // display the probability this leaf returns 1
    return to_string(float(counts[1]) / solutions.size());
  }
  // convert the children to strings
  string zero = static_cast<string>(*left);
  string one = static_cast<string>(*right);
  ostringstream out;
  // combine the strings of the children with the splitting index
  out << "(" << index << ", " << zero << ", " << one << ")";
  return out.str();
}

// Utility function to output tree to a stream
ostream& operator<<(ostream& out, const Bayesian_Tree& tree) {
  out << static_cast<string>(tree);
  return out;
}

// Add a solution pointer to a tree to be used when splitting leaves
void Bayesian_Tree::add_solution(const vector<bool>* const solution) {
  // leaves hold on to solution pointers and keep counts
  if (is_leaf) {
    solutions.push_back(solution);
    counts[solution->at(index)]++;
  // if not a leaf, decide which child gets this solution
  } else if (solution->at(index)) {
    right->add_solution(solution);
  } else {
    left->add_solution(solution);
  }
}

// Split a leaf based on an index value
void Bayesian_Tree::split(size_t new_index) {
  // Cannot split leaves.
  if (not is_leaf) {
    return;
  }
  is_leaf = false;
  // create new children nodes
  left = new Bayesian_Tree(index);
  right = new Bayesian_Tree(index);
  index = new_index;
  // move the solutions based on their value for index.
  for (const auto& solution : solutions) {
    if (solution->at(index)) {
      right->add_solution(solution);
    } else {
      left->add_solution(solution);
    }
  }
  // stop storing solutions at this internal node
  solutions.clear();
}

// Join the subtrees of this internal node back to a single node
void Bayesian_Tree::join() {
  // cannot join leaves
  if (is_leaf) {
    return;
  }
  // ensure both children are leaves
  left->join();
  right->join();
  // copy up solutions in children back to this node
  solutions = left->solutions;
  solutions.insert(solutions.end(), right->solutions.begin(),
                   right->solutions.end());
  index = left->index;
  // copy counts from children back to this node
  counts[0] = left->counts[0] + right->counts[0];
  counts[1] = left->counts[1] + right->counts[1];
  // clean up children
  delete right;
  delete left;
  right = left = nullptr;
  is_leaf = true;
}

// Calculate the BDe score of this tree
Factorial_Fraction Bayesian_Tree::bde_fraction() {
  if (is_leaf) {
    Factorial_Fraction result;
    result.div_by_factorial(solutions.size() + 1);
    // how many solutions have the index bit set to 0 or 1
    result.mul_by_factorial(counts[0]);
    result.mul_by_factorial(counts[1]);
    return result;
  } else {
    // Internal nodes are the product of their children's scores
    return left->bde_fraction() * right->bde_fraction();
  }
}

// Calculate the BDe effect of splitting this leaf without actually splitting it
Factorial_Fraction Bayesian_Tree::splitless_bde(size_t new_index) {
  // count how frequently the 2 combinations of new_index and index exist in solutions.
  array<size_t, 4> splitcount = { 0, 0, 0, 0 };
  for (const auto& solution : solutions) {
    splitcount[(solution->at(new_index) << 1) | solution->at(index)]++;
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

// Use the tree and the current solution to determine the probabilistic value of "index"
void Bayesian_Tree::set_bit(Random& rand, vector<bool>& solution) {
  if (is_leaf) {
    // If solutions==0, probability is even
    double probability = 0.5;
    if (solutions.size()) {
      // set based on proportion of solutions which reached this leaf with value 1
      probability = static_cast<double>(counts[1]) / solutions.size();
    }
    bernoulli_distribution dist(probability);
    solution[index] = dist(rand);
  } else {
    // recurse based on internal node's index until you reach a leaf.
    if (solution[index]) {
      right->set_bit(rand, solution);
    } else {
      left->set_bit(rand, solution);
    }
  }
}

// Initialize a collection of n empty trees
Bayesian_Forest::Bayesian_Forest(size_t n)
    : trees(0),
      ordering(n),
      added(0) {
  for (size_t i = 0; i < n; i++) {
    trees.emplace_back(i);
  }
  // acyclic dependency ordering of trees
  iota(ordering.begin(), ordering.end(), 0);
}

// Add the solution pointer to each tree in the forest
void Bayesian_Forest::add_solution(const vector<bool>* solution) {
  added += 1;
  for (auto& tree : trees) {
    tree.add_solution(solution);
  }
}

// Construct a forest from the added solutions
void Bayesian_Forest::build_forest(Random & rand) {
  // Initially all trees have themselves as both prev and post
  prev.resize(trees.size());
  post.resize(trees.size());
  for (size_t i = 0; i < trees.size(); i++) {
    prev[i] = {i};
    post[i] = {i};
  }
  // Keep track of which splits are the most useful
  vector<tuple<double, Bayesian_Tree*, size_t>> options;
  for (auto& tree : trees) {
    // ensure it starts as a leaf.
    tree.join();
    // add to "options" all possible splits of this tree
    build_options(tree, options);
  }

  // Note: The inequality is inverted from the paper
  // ensures minimum quality for the split
  double threshold = pow(2, -0.5 * log2(added));
  // remove all options worse than threshold, and break ties randomly
  filter(options, threshold);
  shuffle(options.begin(), options.end(), rand);
  // while there are still worthwhile splits
  while (options.size()) {
    // get the best split
    auto best = *min_element(options.begin(), options.end());
    Bayesian_Tree& splitting = *get<1>(best);
    size_t new_index = get<2>(best);
    // update dependency information to prevent cycles
    auto* new_req = &(post[splitting.index]);
    for (size_t upstream : prev[new_index]) {
      post[upstream].insert(new_req->begin(), new_req->end());
    }

    new_req = &(prev[new_index]);
    for (size_t down_stream : post[splitting.index]) {
      prev[down_stream].insert(new_req->begin(), new_req->end());
    }
    // perform the chosen split
    splitting.split(new_index);

    // build new options for the two new trees
    build_options(*splitting.left, options);
    build_options(*splitting.right, options);

    // clear out any invalid splits
    filter(options, threshold);
    shuffle(options.begin(), options.end(), rand);
  }

  // determine what order bits can be generated to ensure all dependencies
  // are met
  sort(ordering.begin(), ordering.end(),
       [this](size_t a, size_t b) {return prev[a].size() < prev[b].size();});
}

// Use each tree in the forest to set all bits in the solution
void Bayesian_Forest::generate(Random& rand, vector<bool> & solution) {
  for (const auto& index : ordering) {
    trees[index].set_bit(rand, solution);
  }
}

// Utility function to display all trees in the forest
std::ostream& operator<<(std::ostream& out, const Bayesian_Forest& forest) {
  for (size_t i = 0; i < forest.trees.size(); i++) {
    out << i << ": " << forest.trees[i] << endl;
  }
  return out;
}

// Internal function which determines all possible ways of splitting a tree
// and adds them to options
void Bayesian_Forest::build_options(
    Bayesian_Tree& tree,
    vector<std::tuple<double, Bayesian_Tree*, size_t>>& options) const {
  // BDe score prior to split
  Factorial_Fraction initial = tree.bde_fraction();
  const auto & invalid = post[tree.index];
  for (size_t i = 0; i < trees.size(); i++) {
    // if splitting on i does not create a cycle
    if (invalid.count(i) == 0) {
      // BDe score after split
      Factorial_Fraction with_split = tree.splitless_bde(i);
      // Ratio improvement
      Factorial_Fraction result = initial / with_split;
      result.simplify();
      options.emplace_back(result, &tree, i);
    }
  }
}

// Remove any split which is no longer useful
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

// Set up initial population
HBOA::HBOA(Random& _rand, shared_ptr<Evaluator> _evaluator,
           Configuration& _config)
    : Optimizer(_rand, _evaluator, _config) {
  size_t pop_size = config.get<int>("pop_size");

  // optionally applies a hill climber to the initial solutions
  hc = config.get<hill_climb::pointer>("hill_climber");

  float fitness;
  for (size_t i = 0; i < pop_size; i++) {
    // create and evaluate solutions
    auto solution = rand_vector(rand, length);
    fitness = evaluator->evaluate(solution);
    // Apply hill climber if configured to do so
    hc(rand, solution, fitness, evaluator);
    solutions.push_back(solution);
    fitnesses.push_back(fitness);
  }
  // set up tools for performing selection
  selection.resize(pop_size);
  iota(selection.begin(), selection.end(), 0);
  rtr_size = min(length, pop_size / 20);
  if (rtr_size == 0) {
    rtr_size = 1;
  }
  completed_generations = 0;
}

bool HBOA::iterate() {
  completed_generations++;
  Bayesian_Forest model(length);
  // Perform a binary tournament to select members of the model
  uniform_int_distribution<int> choose(0, solutions.size() - 1);
  int x, y;
  for (size_t i = 0; i < solutions.size(); i++) {
    x = choose(rand);
    y = choose(rand);
    // choose the better of x and y
    if (fitnesses[x] < fitnesses[y]) {
      model.add_solution(&(solutions[y]));
    } else {
      model.add_solution(&(solutions[x]));
    }
  }

  // Build the model
  model.build_forest(rand);

  // storage for the newly generated solution
  float fitness;
  vector<bool> solution(length);
  for (size_t i = 0; i < solutions.size(); i++) {
    // create and evaluate the solution
    model.generate(rand, solution);
    fitness = evaluator->evaluate(solution);

    // perform restricted tournament replacement
    int choice = rtr_nearest(solution);
    if (fitnesses[choice] < fitness) {
      solutions[choice] = solution;
      fitnesses[choice] = fitness;
    }
  }
  // termination is when generations equals length
  return completed_generations < length;
}

// Finds the current population member most similar to solution from
// a random sample
int HBOA::rtr_nearest(const vector<bool>& solution) {
  int index;
  int choice = -1;
  int distance = length + 1;
  for (size_t i = 0; i < rtr_size; i++) {
    // only look at unused indices.
    index = std::uniform_int_distribution<int>(i, selection.size() - 1)(rand);
    swap(selection[i], selection[index]);
    // Find the distance between the selection and the solution
    int new_dist = hamming_distance(solutions[selection[i]], solution);
    // keep it if it is closer than any previous selection
    if (new_dist < distance) {
      distance = new_dist;
      choice = selection[i];
    }
  }
  return choice;
}
