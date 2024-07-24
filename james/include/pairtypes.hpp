#pragma once

#include <cstddef>
#include <functional>
#include <stdexcept>
#include <unordered_map>
#include <vector>

namespace James::Bond {
// Commutative pair class
class Pair {
public:
  int typeA{};
  int typeB{};

  Pair(int typeA, int typeB) : typeA(typeA), typeB(typeB) {}

  // overloaded equality operator
  friend bool operator==(const Pair &p1, const Pair &p2) {
    if ((p1.typeA == p2.typeA && p1.typeB == p2.typeB) ||
        (p1.typeA == p2.typeB && p1.typeB == p2.typeA)) {
      return true;
    } else {
      return false;
    }
  }
};

// Custom hash for the Pair class
// Here we implement the custom hash as a standalone function object.
struct MyHash {
  std::size_t operator()(const Pair &s) const noexcept {
    std::size_t h1 = std::hash<int>{}(s.typeA);
    std::size_t h2 = std::hash<int>{}(s.typeB);
    if (s.typeA > s.typeB) {
      return h1 ^ h2;
    } else {
      return h2 ^ h1;
    }
  }
};

// Create an unordered_map with the Pair as the key and the cutoff as the
// value
inline std::unordered_map<Pair, double, MyHash>
create_pairtype_cutoffs(std::vector<Pair> &pairs,
                        std::vector<double> &cutoffs) {
  auto cutoff_map = std::unordered_map<Pair, double, MyHash>{};

  if (pairs.size() != cutoffs.size()) {
    std::runtime_error("Inconsistent Pair and cutoff values\n");
  }

  for (size_t i = 0; i < pairs.size(); i++) {
    cutoff_map[pairs[i]] = cutoffs[i];
  }

  return cutoff_map;
}
} // namespace James::Bond