#include <algorithm>
#include <cassert>
#include <compare>
#include <cstdlib>
#include <iostream>
#include <iterator>
#include <ostream>
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
    for (auto it : members) {
      ostrm << it << ' ';
    }
    ostrm << '\n';
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

auto SpatialGrouping(std::istream&, std::ostream&, Problem const& problem)
    -> std::vector<Group> {
  std::vector<Group> groups(problem.m);
  std::vector<int> spatial_nodes(problem.n);
  for (int i = 0; i < spatial_nodes.size(); ++i) {
    spatial_nodes[i] = i;
  }

  std::sort(spatial_nodes.begin(), spatial_nodes.end(),
            [&boxes = problem.boxes](int a, int b) -> bool {
              return BoxSpatialId(kBlockLength, boxes[a]) <
                     BoxSpatialId(kBlockLength, boxes[b]);
            });

  for (int i = 0, j = 0, k = 0; i < problem.n; ++i) {
    groups[j].members.push_back(spatial_nodes[i]);
    if (++k == problem.group_sizes[j]) {
      k = 0;
      ++j;
    }
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
