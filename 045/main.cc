#include <algorithm>
#include <cassert>
#include <compare>
#include <cstdlib>
#include <functional>
#include <iostream>
#include <iterator>
#include <limits>
#include <numeric>
#include <optional>
#include <ostream>
#include <set>
#include <span>
#include <type_traits>
#include <unordered_set>
#include <vector>

namespace tbrekalo {

// https://prng.di.unimi.it/xoshiro256plusplus.c
class Xoshiro256pp {
  static constexpr std::array<std::uint64_t, 4> kJump128{
      0x180ec6d33cfd0aba, 0xd5a61266f0c9392c, 0xa9582618e03fc9aa,
      0x39abdc4529b1661c};

  static constexpr std::array<std::uint64_t, 4> kJump192{
      0x76e15d3efefdcbbf, 0xc5004e441c522fb3, 0x77710069854ee241,
      0x39109bb02acbe635};

  [[gnu::always_inline]] static constexpr std::uint64_t rotl(
      const std::uint64_t x, int k) {
    return (x << k) | (x >> (64 - k));
  }

  constexpr auto jumpImpl(std::array<std::uint64_t, 4> const& jump)
      -> Xoshiro256pp {
    auto result = *this;
    std::uint64_t s0 = 0;
    std::uint64_t s1 = 0;
    std::uint64_t s2 = 0;
    std::uint64_t s3 = 0;
    for (std::uint64_t j : jump) {
      for (int b = 0; b < 64; ++b) {
        if (j & (static_cast<std::uint64_t>(1) << b)) {
          s0 ^= s_[0];
          s1 ^= s_[1];
          s2 ^= s_[2];
          s3 ^= s_[3];
        }
        next();
      }
    }
    s_[0] = s0;
    s_[1] = s1;
    s_[2] = s2;
    s_[3] = s3;

    return result;
  }

  std::array<std::uint64_t, 4> s_;

 public:
  using result_type = std::uint64_t;
  using value_type = std::array<std::uint64_t, 4>;

  constexpr Xoshiro256pp() = delete;
  constexpr explicit Xoshiro256pp(value_type const& state) : s_(state) {}

  constexpr static auto min() noexcept -> result_type {
    return std::numeric_limits<result_type>::min();
  }

  constexpr static auto max() noexcept -> result_type {
    return std::numeric_limits<result_type>::max();
  }

  constexpr auto next() noexcept -> result_type {
    const std::uint64_t result = rotl(s_[0] + s_[3], 23) + s_[0];
    const std::uint64_t t = s_[1] << 17;

    s_[2] ^= s_[0];
    s_[3] ^= s_[1];
    s_[1] ^= s_[2];
    s_[0] ^= s_[3];

    s_[2] ^= t;

    s_[3] = rotl(s_[3], 45);

    return result;
  }

  constexpr auto operator()() noexcept -> result_type { return next(); }

  /* This is the jump function for the generator. It is equivalent
     to 2^128 calls to next(); it can be used to generate 2^128
     non-overlapping subsequences for parallel computations. */
  constexpr auto jump() noexcept -> Xoshiro256pp { return jumpImpl(kJump128); }

  /* This is the long-jump function for the generator. It is equivalent to
     2^192 calls to next(); it can be used to generate 2^64 starting points,
     from each of which jump() will generate 2^64 non-overlapping
     subsequences for parallel distributed computations. */
  constexpr auto longJump(void) noexcept -> Xoshiro256pp {
    return jumpImpl(kJump192);
  }
};

}  // namespace tbrekalo

namespace tbrekalo {

static constexpr int kMaxXY = 10'000;
static constexpr int kBlockLength = 1'000;

static Xoshiro256pp rng(Xoshiro256pp::value_type{1, 2, 3, 4});

struct Box {
  int lx;
  int rx;
  int ly;
  int ry;
};

struct Problem {
  int n;
  int m;
  int q;
  int l;
  int w;

  std::vector<int> group_sizes;
  std::vector<Box> boxes;
};

struct Road {
  int a;
  int b;
};

struct Group {
  std::vector<int> members;
  std::vector<Road> roads;
};

struct ClusterId {
  int x;
  int y;

  friend constexpr auto operator<=>(ClusterId lhs, ClusterId rhs) noexcept
      -> std::strong_ordering = default;
};

struct Coord {
  int x;
  int y;

  friend constexpr auto operator<=>(Coord, Coord) noexcept -> Coord = default;
  friend constexpr auto operator+(Coord lhs, Coord rhs) noexcept -> Coord {
    return {
        .x = lhs.x + rhs.x,
        .y = lhs.y + rhs.y,
    };
  };
};

static auto LoadProblem(std::istream& istrm) -> Problem {
  Problem dst{};

  istrm >> dst.n >> dst.m >> dst.q >> dst.l >> dst.w;

  dst.group_sizes.resize(dst.m);
  for (int i = 0; i < dst.m; ++i) {
    istrm >> dst.group_sizes[i];
  }

  dst.boxes.resize(dst.n);
  for (int i = 0; i < dst.n; ++i) {
    /* clang-format off */
    istrm >>
      dst.boxes[i].lx >>
      dst.boxes[i].rx >>
      dst.boxes[i].ly >>
      dst.boxes[i].ry;
    /* clang-format on */

    std::clamp(dst.boxes[i].lx, 0, kMaxXY);
    std::clamp(dst.boxes[i].rx, 0, kMaxXY);

    std::clamp(dst.boxes[i].ly, 0, kMaxXY);
    std::clamp(dst.boxes[i].ry, 0, kMaxXY);
  }

  return dst;
};

static auto PrintSolution(std::ostream& ostrm, std::span<Group const> groups) {
  ostrm << "!\n";
  for (auto const& [members, roads] : groups) {
    for (auto it = members.begin(); it != members.end(); ++it) {
      ostrm << *it << (std::next(it) != members.end() ? ' ' : '\n');
    }
    for (auto [a, b] : roads) {
      ostrm << a << ' ' << b << '\n';
    }
  }
  std::flush(ostrm);
}

static auto Query(std::istream& istrm, std::ostream& ostrm, int l,
                  std::span<int const> cities) -> std::vector<Road> {
  assert(cities.size() >= 3);
  ostrm << "? " << l << " ";
  for (auto i = 0; i < cities.size(); ++i) {
    ostrm << cities[i] << (i + 1 != cities.size() ? ' ' : '\n');
  }

  std::flush(ostrm);
  std::vector<Road> dst(cities.size() - 1);
  for (int i = 0; i < dst.size(); ++i) {
    istrm >> dst[i].a >> dst[i].b;
    if (dst[i].a > dst[i].b) {
      std::swap(dst[i].a, dst[i].b);
    }
  }

  return dst;
}

static auto NaiveGroupRoads(std::span<int const> members) -> std::vector<Road> {
  assert(members.size() >= 1);
  if (members.size() == 1) {
    return {};
  }

  std::vector<Road> dst(members.size() - 1);

  auto lhs = members.begin();
  auto rhs = std::next(members.begin());
  for (int i = 0; rhs != members.end(); ++lhs, ++rhs, ++i) {
    assert(i < dst.size());
    dst[i] = Road{.a = *lhs, .b = *rhs};
  }

  return dst;
}

;

static auto BoxAvgCoord(Box box) -> Coord {
  return Coord{.x = box.lx + (box.rx - box.lx) / 2,
               .y = box.ly + (box.ry - box.ly) / 2};
}

static auto BoxClusterId(int block_length, Box box) -> ClusterId {
  assert(block_length != 0);

  return {
      .x = (box.lx + (box.rx - box.lx) / 2) / block_length,
      .y = (box.ly + (box.ry - box.ly) / 2) / block_length,
  };
}

template <class T>
static auto CalcDist(T lhs, T rhs) -> int
  requires(std::is_same_v<T, ClusterId> || std::is_same_v<T, Coord>)
{
  return (lhs.x - rhs.x) * (lhs.x - rhs.x) + (lhs.y - rhs.y) * (lhs.y - rhs.y);
}

static auto NaiveClustering(std::istream&, std::ostream&,
                            Problem const& problem) -> std::vector<Group> {
  struct GroupSpec {
    int id;
    int size;
  };

  struct Edge {
    int from;
    int to;
    int dist;
  };

  struct EdgeCmp {
    constexpr auto operator()(Edge lhs, Edge rhs) const noexcept -> bool {
      return lhs.dist < rhs.dist;
    }
  };

  std::unordered_set<int> ungrouped;
  std::vector<int> grouped(problem.n, -1);
  std::vector<std::vector<int>> dists(
      problem.n, std::vector<int>(problem.n, kMaxXY * kMaxXY));

  for (int i = 1; i < problem.n; ++i) {
    for (int j = 0; j < i; ++j) {
      dists[i][j] = dists[j][i] = CalcDist(BoxAvgCoord(problem.boxes[i]),
                                           BoxAvgCoord(problem.boxes[j]));
    }
  }

  for (int i = 0; i < problem.n; ++i) {
    ungrouped.insert(i);
  }

  std::vector<Group> groups(problem.m);
  std::vector<GroupSpec> group_specs(problem.m);
  for (int i = 0; i < problem.m; ++i) {
    group_specs[i] = {
        .id = i,
        .size = problem.group_sizes[i],
    };
  }

  int q = problem.q / 2;
  std::ranges::sort(group_specs, std::greater<>{}, &GroupSpec::size);

  for (auto [group_id, group_target_size] : group_specs) {
    std::set<Edge, EdgeCmp> edges;
    auto insert_edges = [n = problem.n, &grouped, &dists, &edges](int from) {
      for (int i = 0; i < n; ++i) {
        if (grouped[i] != -1) {
          continue;
        }

        edges.insert(Edge{.from = from, .to = i, .dist = dists[from][i]});
      }
    };

    edges.insert(
        Edge{.from = *ungrouped.begin(), .to = *ungrouped.begin(), .dist = 0});
    while (!edges.empty()) {
      auto [from, city, _] = *edges.begin();
      edges.erase(*edges.begin());

      if (grouped[city] != -1) {
        continue;
      }

      groups[group_id].members.push_back(city);
      grouped[city] = group_id;
      ungrouped.erase(city);

      if (groups[group_id].members.size() == group_target_size) {
        break;
      }

      insert_edges(city);
    }
  }

  return groups;
}

static auto OptimizeRoads(std::istream& istrm, std::ostream& ostrm,
                          std::span<Box const> boxes, std::vector<Group> groups,
                          int q, int l) -> std::vector<Group> {
  assert(l >= 3);
  for (int group_id = 0; group_id < groups.size(); ++group_id) {
    auto const& members = groups[group_id].members;
    auto& roads = groups[group_id].roads;

    for (auto lo = 0; lo + 1 < members.size();) {
      std::vector<Road> insertion;
      auto len = std::min<int>(l, members.size() - lo);
      if (len >= 3 && q > 0) {
        std::vector<int> buffer;
        buffer.reserve(len);
        buffer.insert(buffer.end(), members.begin() + lo,
                      members.begin() + lo + len);

        insertion = Query(istrm, ostrm, len, buffer);
        --q;

      } else {
        insertion =
            NaiveGroupRoads(std::span(members.begin() + lo, members.end()));
        len = members.size() - lo;
      }

      roads.insert(roads.end(), insertion.begin(), insertion.end());
      lo += len - 1;
    }
  }

  return groups;
}

}  // namespace tbrekalo

namespace tb = tbrekalo;

auto main(void) -> int {
  auto problem = tb::LoadProblem(std::cin);
  auto groups = tb::OptimizeRoads(
      std::cin, std::cout, problem.boxes,
      tb::NaiveClustering(std::cin, std::cout, problem), problem.q, problem.l);

  tb::PrintSolution(std::cout, groups);
}
