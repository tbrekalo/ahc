#include <cassert>
#include <compare>
#include <cstdlib>
#include <deque>
#include <iostream>
#include <iterator>
#include <ostream>
#include <unordered_map>
#include <unordered_set>
#include <vector>

namespace tbrekalo {

static constexpr int kMaxXY = 10'000;
static constexpr int kBlockLength = 1'000;

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

struct SpatialId {
  int x;
  int y;

  friend constexpr auto operator<=>(SpatialId lhs, SpatialId rhs) noexcept
      -> std::strong_ordering = default;
};

struct SpatialIdHash {
  constexpr auto operator()(SpatialId id) const noexcept {
    return (id.x << 5) + id.x + id.y;
  }
};

using SpatialMap =
    std::unordered_map<SpatialId, std::unordered_set<int>, SpatialIdHash>;

auto LoadProblem(std::istream& istrm) -> Problem {
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
  }

  return dst;
};

auto PrintSolution(std::ostream& ostrm, std::vector<Group> const& groups) {
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

auto NaiveGroupRoads(std::vector<int> const& members) -> std::vector<Road> {
  assert(members.size() >= 2);
  std::vector<Road> dst(members.size() - 1);

  auto lhs = members.begin();
  auto rhs = std::next(members.begin());
  for (int i = 0; rhs != members.end(); ++lhs, ++rhs, ++i) {
    dst[i] = Road{.a = *lhs, .b = *rhs};
  }

  return dst;
}

auto BoxSpatialId(int block_length, Box box) -> SpatialId {
  return {
      .x = (box.lx + (box.rx - box.lx) / 2) / block_length,
      .y = (box.ly + (box.ry - box.ly) / 2) / block_length,
  };
}

auto ExpandFrontier(std::vector<std::vector<char>> const& visited, SpatialId id)
    -> std::vector<SpatialId> {
  std::vector<SpatialId> dst;
  auto const max_id = visited.size();
  for (int dx = -1; dx <= 1; ++dx) {
    for (int dy = -1; dy <= 1; ++dy) {
      int x = id.x + dx;
      int y = id.y + dy;
      if (dx == 0 && dy == 0 || x < 0 || x >= max_id || y < 0 || y >= max_id ||
          visited[x][y]) {
        continue;
      }
    }
  }
  return dst;
}

auto SpatialGrouping(std::istream&, std::ostream&, Problem const& problem)
    -> std::vector<Group> {
  int const block_length = kBlockLength;
  int const n_blocks = kMaxXY / block_length + 1;

  std::vector<Group> groups(problem.m);
  for (int i = 0; i < groups.size(); ++i) {
    groups[i].members.reserve(problem.group_sizes[i]);
  }

  std::unordered_set<int> ungrouped;
  std::vector<std::vector<std::vector<int>>> areas(
      n_blocks, std::vector<std::vector<int>>(n_blocks));
  std::vector<std::vector<char>> visited(n_blocks,
                                         std::vector<char>(n_blocks, 0));

  for (int i = 0; i < problem.n; ++i) {
    auto [x, y] = BoxSpatialId(kBlockLength, problem.boxes[i]);
    assert(x < n_blocks);
    assert(y < n_blocks);

    areas[x][y].push_back(i);
    ungrouped.insert(i);
  }

  std::vector<SpatialId> stack{SpatialId{.x = n_blocks / 2, .y = n_blocks / 2}};
  for (int g_id = 0; g_id < problem.m;) {
  }

  for (int i = 0; i < problem.m; ++i) {
    groups[i].roads = NaiveGroupRoads(groups[i].members);
  }

  return groups;
}

}  // namespace tbrekalo
namespace tb = tbrekalo;

auto main(void) -> int {
  tb::PrintSolution(std::cout, tb::SpatialGrouping(std::cin, std::cout,
                                                   tb::LoadProblem(std::cin)));
}
