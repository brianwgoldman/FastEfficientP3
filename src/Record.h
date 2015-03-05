// Brian Goldman

// A record stores the result information for a single optimization run,
// recording how many evaluations are required to reach each new fitness level

#ifndef RECORD_H_
#define RECORD_H_

#include <vector>
#include <utility>
#include <array>
#include <stdlib.h>
#include <chrono>
#include "Util.h"
#include "Configuration.h"

using std::vector;
using std::pair;
using std::array;

// Labels used for decoding the "summary" return statistics
enum RecordKey {
  MES,
  MAD,
  FAILURES
};

// Stores optimization history information
class Record {
 public:
  Record(Configuration& config)
      : verbosity(config.get<int>("verbosity")) {
  }
  virtual ~Record() = default;
  // Tracks increases in fitness.  If verbosity is set high enough, outputs
  // information to the screen.
  void add(float fitness, int evals);
  // The best fitness achieved in this optimization run
  const pair<float, int>& best() const;
  // Summarizes the end of run statistics for mutliple records
  static array<int, 3> summarize(const vector<Record>& records,
                                 Configuration& config);
  // access to the stored values:
  const vector<pair<float, int>>& progression() {
    return levels;
  }

  // Used to track total run time.
  void start_clock() {
    start_time = std::chrono::steady_clock::now();
  }

  // Returns the number of seconds elapsed betwee start_clock() and best fitness found.
  float seconds_used() {
    return std::chrono::duration<float>(last_update_time - start_time).count();
  }
  string metadata;
 private:
  // Raw stored data
  vector<pair<float, int>> levels;
  std::chrono::steady_clock::time_point start_time;
  std::chrono::steady_clock::time_point last_update_time;
  int verbosity;
};

#endif /* RECORD_H_ */
