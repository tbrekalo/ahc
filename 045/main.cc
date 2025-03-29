#include <cstdlib>
#include <iostream>
#include <ostream>
#include <vector>

namespace tbrekalo {

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
    for (int i = 0; i < members.size(); ++i) {
      ostrm << members[i] << (i + 1 != members.size() ? ' ' : '\n');
    }
    for (auto [a, b] : roads) {
      ostrm << a << ' ' << b << '\n';
    }
  }
  std::flush(ostrm);
}

auto Solve(std::istream& istrm, std::ostream& ostrm, Problem const& problem)
    -> std::vector<Group> {
  std::vector<Group> groups(problem.m);
  for (int i = 0, j = 0; j < problem.m; ++j) {
    groups[j].members.resize(problem.group_sizes[j]);
    for (int k = 0; k < problem.group_sizes[j]; ++k) {
      groups[j].members[k] = i;
      if (k > 0) {
        groups[j].roads.push_back(Road{.a = i - 1, .b = i});
      }
      ++i;
    }
  }

  return groups;
}

}  // namespace tbrekalo
namespace tb = tbrekalo;

auto main(void) -> int {
  tb::PrintSolution(std::cout,
                    tb::Solve(std::cin, std::cout, tb::LoadProblem(std::cin)));
}
