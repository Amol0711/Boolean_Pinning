// Linear-work support certificate for planted two-input Boolean networks.
// No random generation occurs here. Rules are g10, g01, g11, with g00 = 0.
#include <algorithm>
#include <cstdint>
#include <vector>

extern "C" int certify_core(std::int64_t n, const std::int64_t* parents,
                            const std::uint8_t* rules, std::uint8_t* alive,
                            std::int64_t* residual_parent, std::int64_t* pins,
                            std::int64_t* stats) noexcept {
  try {
    std::vector<std::int64_t> offset(n + 1, 0), cursor(n), slots(2*n);
    for (std::int64_t e = 0; e < 2*n; ++e) ++offset[parents[e]+1];
    for (std::int64_t i = 0; i < n; ++i) offset[i+1] += offset[i];
    for (std::int64_t i = 0; i < n; ++i) cursor[i] = offset[i];
    for (std::int64_t e = 0; e < 2*n; ++e) slots[cursor[parents[e]]++] = e;
    std::vector<std::uint8_t> mask(n, 3);
    std::vector<std::int64_t> queue;
    queue.reserve(n);
    for (std::int64_t i = 0; i < n; ++i) {
      alive[i] = (rules[3*i] | rules[3*i+1] | rules[3*i+2]) != 0;
      residual_parent[i] = -1;
      if (!alive[i]) queue.push_back(i);
    }
    for (std::size_t head = 0; head < queue.size(); ++head) {
      const auto source = queue[head];
      for (auto j = offset[source]; j < offset[source+1]; ++j) {
        const auto e = slots[j];
        const auto i = e/2;
        if (!alive[i]) continue;
        mask[i] &= static_cast<std::uint8_t>(~(1u << (e%2)));
        const auto m = mask[i];
        const bool supported = (m == 3) || (m == 1 && rules[3*i]) ||
                               (m == 2 && rules[3*i+1]);
        if (!supported) {
          alive[i] = 0;
          queue.push_back(i);
        }
      }
    }
    std::int64_t r = 0, doubles = 0, npins = 0;
    for (std::int64_t i = 0; i < n; ++i) {
      if (!alive[i]) continue;
      ++r;
      const bool first = alive[parents[2*i]];
      const bool second = alive[parents[2*i+1]];
      if (first && second) ++doubles;
      else if (first) residual_parent[i] = parents[2*i];
      else if (second) residual_parent[i] = parents[2*i+1];
      else return 2;  // A nonempty fixed-support row cannot have zero slots.
    }
    if (doubles == 0) {
      std::vector<std::int64_t> seen(n, 0);
      std::vector<std::uint8_t> chosen(n, 0);
      for (std::int64_t i = 0; i < n; ++i) {
        if (!alive[i] || seen[i]) continue;
        const auto stamp = i+1;
        auto v = i;
        while (!seen[v]) {
          seen[v] = stamp;
          v = residual_parent[v];
        }
        if (seen[v] == stamp) {
          const auto entry = v;
          auto smallest = v;
          v = residual_parent[v];
          while (v != entry) {
            smallest = std::min(smallest, v);
            v = residual_parent[v];
          }
          chosen[smallest] = 1;
        }
      }
      for (std::int64_t i = 0; i < n; ++i)
        if (chosen[i]) pins[npins++] = i;
    }
    stats[0] = r;
    stats[1] = doubles;
    stats[2] = doubles ? -1 : npins;
    stats[3] = npins;
    return 0;
  } catch (...) {
    return 1;
  }
}
