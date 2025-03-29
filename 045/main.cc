#include <algorithm>
#include <cassert>
#include <compare>
#include <cstdlib>
#include <functional>
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

    std::clamp(dst.boxes[i].lx, 0, kMaxXY);
    std::clamp(dst.boxes[i].rx, 0, kMaxXY);

    std::clamp(dst.boxes[i].ly, 0, kMaxXY);
    std::clamp(dst.boxes[i].ry, 0, kMaxXY);
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

auto BoxSpatialId(int block_length, Box box) -> SpatialId {
  assert(block_length != 0);

  return {
      .x = (box.lx + (box.rx - box.lx) / 2) / block_length,
      .y = (box.ly + (box.ry - box.ly) / 2) / block_length,
  };
}

auto SpatialIdDist(SpatialId lhs, SpatialId rhs) -> int {
  return (lhs.x - rhs.x) * (lhs.x - rhs.x) + (lhs.y - rhs.y) * (lhs.y - rhs.y);
}

auto SpatialGrouping(std::istream&, std::ostream&, Problem const& problem)
    -> std::vector<Group> {
  struct GroupSpec {
    int id;
    int size;
  };

  int const block_length = kBlockLength;
  std::vector<Group> groups(problem.m);
  std::vector<int> buffer(problem.n);

  for (int city_id = 0; city_id < problem.n; ++city_id) {
    buffer[city_id] = city_id;
  }

  std::vector<GroupSpec> group_specs(problem.m);
  for (int i = 0; i < problem.m; ++i) {
    group_specs[i] = GroupSpec{.id = i, .size = problem.group_sizes[i]};
  }

  std::ranges::sort(group_specs, std::greater<>{}, &GroupSpec::size);
  for (auto [group_id, size] : group_specs) {
    groups[group_id].members.reserve(size);
    std::sort(buffer.begin() + 1, buffer.end(),
              [block_length, &boxes = problem.boxes, ref = *buffer.begin()](
                  int a, int b) {
                return SpatialIdDist(BoxSpatialId(block_length, boxes[a]),
                                     BoxSpatialId(block_length, boxes[ref])) <
                       SpatialIdDist(BoxSpatialId(block_length, boxes[b]),
                                     BoxSpatialId(block_length, boxes[ref]));
              });
    std::copy(buffer.begin(), buffer.begin() + size,
              std::back_inserter(groups[group_id].members));
    buffer.erase(buffer.begin(), buffer.begin() + size);
    groups[group_id].roads = NaiveGroupRoads(groups[group_id].members);
  }

  return groups;
}

}  // namespace tbrekalo

namespace tb = tbrekalo;

auto main(void) -> int {
  tb::PrintSolution(std::cout, tb::SpatialGrouping(std::cin, std::cout,
                                                   tb::LoadProblem(std::cin)));
}
