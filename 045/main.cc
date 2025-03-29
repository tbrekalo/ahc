#include <cassert>
#include <cstdlib>
#include <iostream>
#include <iterator>
#include <ostream>
#include <vector>

namespace tbrekalo {

static constexpr int kMaxXY = 10'000;
static constexpr int kBlockLength = 1'000;
static_assert(kMaxXY % kBlockLength == 0);

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

auto SpatialGrouping(std::istream&, std::ostream&, Problem const& problem)
    -> std::vector<Group> {
  std::vector<Group> dst(problem.m);
  auto const n_blocks = kMaxXY / kBlockLength + 1;
  std::vector<std::vector<std::vector<int>>> grid(
      n_blocks, std::vector<std::vector<int>>(n_blocks));

  for (auto i = 0; i < problem.n; ++i) {
    auto const x =
        problem.boxes[i].lx + (problem.boxes[i].rx - problem.boxes[i].lx) / 2;
    auto const y =
        problem.boxes[i].ly + (problem.boxes[i].ry - problem.boxes[i].ly) / 2;

    grid[x / kBlockLength][y / kBlockLength].push_back(i);
  }

  std::vector<int> node_group(problem.n);
  for (int i = 0, j = 0, k = 0; i < problem.n; ++i) {
    node_group[i] = j;
    dst[j].members.push_back(i);
    if (++k == problem.group_sizes[j]) {
      k = 0;
      ++j;
    }
  }

  for (int i = 0; i < problem.m; ++i) {
    dst[i].roads = NaiveGroupRoads(dst[i].members);
  }

  return dst;
}

}  // namespace tbrekalo
namespace tb = tbrekalo;

auto main(void) -> int {
  tb::PrintSolution(std::cout, tb::SpatialGrouping(std::cin, std::cout,
                                                   tb::LoadProblem(std::cin)));
}
