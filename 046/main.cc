#include <array>
#include <cassert>
#include <compare>
#include <iostream>
#include <limits>
#include <ostream>
#include <span>
#include <utility>
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

static inline constexpr int N = 20;
static inline constexpr int M = 40;

struct Coord {
  int r;
  int c;

  friend constexpr auto operator+(Coord lhs, Coord rhs) noexcept -> Coord;
  friend constexpr auto operator-(Coord lhs, Coord rhs) noexcept -> Coord;
  friend constexpr auto operator<=>(Coord lhs, Coord rhs) noexcept
      -> std::strong_ordering = default;
};

[[gnu::always_inline]] inline constexpr auto InBounds(Coord coord) noexcept
    -> bool {
  /* clang-format off */
  auto [r, c] = coord;
  return r >= 0 && r < N &&
         c >= 0 && c < N;
  /* clang-format on */
}

[[gnu::always_inline]] constexpr auto operator+(Coord lhs, Coord rhs) noexcept
    -> Coord {
  Coord dst{.r = lhs.r + rhs.r, .c = lhs.c + rhs.c};
  assert(InBounds(dst));
  return dst;
}
[[gnu::always_inline]] constexpr auto operator-(Coord lhs, Coord rhs) noexcept
    -> Coord {
  Coord dst{.r = lhs.r - rhs.r, .c = lhs.c - rhs.c};
  assert(InBounds(dst));
  return dst;
}

using Coords = std::vector<Coord>;

using Grid = std::array<std::array<char, N>, N>;

[[gnu::always_inline]] inline constexpr auto CanPlace(Grid const& grid,
                                                      Coord coord) noexcept
    -> bool {
  auto [r, c] = coord;
  return InBounds(coord) && grid[coord.r][coord.c] == 0;
}

enum class Action : char { M, S, A };
static auto operator<<(std::ostream& ostrm, Action action) -> std::ostream& {
  static constexpr std::array<char, 3> ACTION_TO_CHAR{'M', 'S', 'A'};
  return ostrm << ACTION_TO_CHAR[std::to_underlying(action)];
}

enum class Dir : char { U, D, L, R };
static auto operator<<(std::ostream& ostrm, Dir dir) -> std::ostream& {
  static constexpr std::array<char, 4> DIR_TO_CHAR{'U', 'D', 'L', 'R'};
  return ostrm << DIR_TO_CHAR[std::to_underlying(dir)];
}

struct Turn {
  Action action;
  Dir dir;
};

static auto operator<<(std::ostream& ostrm, Turn turn) -> std::ostream& {
  return ostrm << turn.action << ' ' << turn.dir;
}

[[gnu::always_inline]] inline constexpr auto DirToDif(Dir dir) -> Coord {
  /* clang-format off */
  constexpr std::array<Coord, 4> DIR_TO_DIF{Coord{-1, 0}, 
                                            Coord{ 1, 0},
                                            Coord{0,  -1},
                                            Coord{0,   1}};
  /* clang-format on */
  return DIR_TO_DIF[std::to_underlying(dir)];
};

static auto Load(std::istream& istrm) -> Coords {
  int _;
  Coords dst;
  istrm >> _ >> _;
  dst.resize(M);
  for (int i = 0; i < M; ++i) {
    istrm >> dst[i].r >> dst[i].c;
  }

  return dst;
}

static auto Print(std::ostream& ostrm, std::span<Turn> turns) -> void {
  for (int i = 0; i < turns.size(); ++i) {
    ostrm << turns[i] << '\n';
  }
  std::flush(ostrm);
}

static auto Solve(Coords const& coords) -> std::vector<Turn> {
  std::vector<Turn> turns;
  Coord cur = coords[0];
  auto repeat_move = [&turns, &cur](Turn turn, int n) {
    for (int i = 0; i < n; ++i) {
      turns.push_back(turn);
      cur = cur + DirToDif(turn.dir);
    }
  };

  for (int i = 1; i < M; ++i) {
    while (cur != coords[i]) {
      if (cur.r < coords[i].r) {
        repeat_move(Turn{.action = Action::M, .dir = Dir::D},
                    coords[i].r - cur.r);
        continue;
      }
      if (cur.r > coords[i].r) {
        repeat_move(Turn{.action = Action::M, .dir = Dir::U},
                    cur.r - coords[i].r);
        continue;
      }

      if (cur.c > coords[i].c) {
        repeat_move(Turn{.action = Action::M, .dir = Dir::L},
                    cur.c - coords[i].c);
        continue;
      }
      if (cur.c < coords[i].c) {
        repeat_move(Turn{.action = Action::M, .dir = Dir::R},
                    coords[i].c - cur.c);
        continue;
      }
    }
  }

  return turns;
}

}  // namespace tbrekalo

namespace tb = tbrekalo;

auto main(int, char**) -> int {
  auto problem = tb::Load(std::cin);
  auto turns = tb::Solve(problem);
  tb::Print(std::cout, turns);
  return 0;
}
