#include <algorithm>
#include <array>
#include <cassert>
#include <chrono>
#include <compare>
#include <iostream>
#include <limits>
#include <numeric>
#include <ostream>
#include <random>
#include <utility>
#include <vector>

namespace tbrekalo {

static constexpr auto MAX_RUNTIME = std::chrono::milliseconds(1'800);
static const auto INIT_TIME = std::chrono::system_clock::now();

static auto elapsed() -> std::chrono::milliseconds {
  return std::chrono::duration_cast<std::chrono::milliseconds>(
      std::chrono::system_clock::now() - INIT_TIME);
}

}  // namespace tbrekalo

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

static Xoshiro256pp RNG{{42, 9, 7, 1998}};
static std::bernoulli_distribution P(0.25);

}  // namespace tbrekalo

namespace tbrekalo {

static inline constexpr int N = 20;
static inline constexpr int M = 40;
static inline constexpr int MAX_TURNS = 2 * N * M;

struct Coord {
  int row;
  int col;

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
  Coord dst{.row = lhs.row + rhs.row, .col = lhs.col + rhs.col};
  assert(InBounds(dst));
  return dst;
}
[[gnu::always_inline]] constexpr auto operator-(Coord lhs, Coord rhs) noexcept
    -> Coord {
  Coord dst{.row = lhs.row - rhs.row, .col = lhs.col - rhs.col};
  assert(InBounds(dst));
  return dst;
}

using Coords = std::vector<Coord>;

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

struct Solution {
  int score;
  std::vector<Turn> turns;
};

static auto Load(std::istream& istrm) -> Coords {
  int _;
  Coords dst;
  istrm >> _ >> _;
  dst.resize(M);
  for (int i = 0; i < M; ++i) {
    istrm >> dst[i].row >> dst[i].col;
  }

  return dst;
}

static auto Print(std::ostream& ostrm, Solution const& solution) -> void {
  std::cerr << "final-score=" << solution.score << std::endl;
  for (int i = 0; i < solution.turns.size(); ++i) {
    ostrm << solution.turns[i] << '\n';
  }
  std::flush(ostrm);
}

struct Range {
  int lo = 0;
  int hi = N - 1;
};

static auto Solve(Coords const& coords, bool alter) -> Solution {
  std::array<int, N> row_cnt;
  std::array<int, N> col_cnt;
  std::array<Range, N> row_range;
  std::array<Range, N> col_range;

  {
    std::ranges::fill(row_cnt, 0);
    std::ranges::fill(col_cnt, 0);
    for (int i = 1; i < M; ++i) {
      ++row_cnt[coords[i].row];
      ++col_cnt[coords[i].col];
    }
  }

  std::vector<Turn> turns;
  auto can_keep_turning = [&turns] -> bool { return turns.size() < MAX_TURNS; };
  Coord cur = coords[0];
  auto repeat_move = [&turns, &can_keep_turning, &cur](Turn turn, int n) {
    for (int i = 0; i < n && can_keep_turning(); ++i) {
      turns.push_back(turn);
      cur = cur + DirToDif(turn.dir);
    }
  };

  int m = 1;
  for (; m < M && can_keep_turning(); ++m) {
    auto [row, col] = coords[m];
    while (cur != coords[m] && can_keep_turning()) {
      // down
      if (cur.row < row) {
        if (row_range[row].hi - row < row - cur.row) {
          turns.push_back(Turn{.action = Action::S, .dir = Dir::D});
          cur = {.row = row_range[row].hi, .col = cur.col};
          continue;
        }
        repeat_move(Turn{.action = Action::M, .dir = Dir::D}, row - cur.row);
        continue;
      }
      // up
      if (cur.row > row) {
        if (row - row_range[row].lo < cur.row - row) {
          turns.push_back(Turn{.action = Action::S, .dir = Dir::U});
          cur = {.row = row_range[row].lo, .col = cur.col};
          continue;
        }
        repeat_move(Turn{.action = Action::M, .dir = Dir::U}, cur.row - row);
        continue;
      }

      // left
      if (cur.col > col) {
        if (col - col_range[col].lo < cur.col - col) {
          turns.push_back(Turn{.action = Action::S, .dir = Dir::L});
          cur = {.row = cur.row, .col = col_range[col].lo};
          continue;
        }
        repeat_move(Turn{.action = Action::M, .dir = Dir::L}, cur.col - col);
        continue;
      }
      // right
      if (cur.col < col) {
        if (col_range[col].hi - col < col - cur.col) {
          turns.push_back(Turn{.action = Action::S, .dir = Dir::R});
          cur = {.row = cur.row, .col = col_range[col].hi};
          continue;
        }
        repeat_move(Turn{.action = Action::M, .dir = Dir::R}, col - cur.col);
        continue;
      }
    }

    if (alter) {
      --row_cnt[row];
      if (m + 1 < M && row > 0 &&
          std::accumulate(row_cnt.begin(), row_cnt.begin() + row, 0,
                          std::plus<>{}) == 0 &&
          can_keep_turning()) {
        turns.push_back(Turn{.action = Action::A, .dir = Dir::U});
        col_range[col].lo = row - 1;
      }
      if (m + 1 < M && row + 1 < N &&
          std::accumulate(row_cnt.begin() + row + 1, row_cnt.end(), 0,
                          std::plus<>{}) == 0 &&
          can_keep_turning()) {
        turns.push_back(Turn{.action = Action::A, .dir = Dir::D});
        col_range[col].lo = row + 1;
      }

      --col_cnt[col];
      if (m + 1 < M && col > 0 &&
          std::accumulate(col_cnt.begin(), col_cnt.begin() + col, 0,
                          std::plus<>{}) == 0 &&
          can_keep_turning()) {
        turns.push_back(Turn{.action = Action::A, .dir = Dir::L});
        col_range[col].lo = col - 1;
      }
      if (m + 1 < M && col + 1 < N &&
          std::accumulate(col_cnt.begin() + col + 1, col_cnt.end(), 0,
                          std::plus<>{}) == 0 &&
          can_keep_turning()) {
        turns.push_back(Turn{.action = Action::A, .dir = Dir::R});
        col_range[col].lo = col + 1;
      }
    }
  }

  assert(turns.size() <= MAX_TURNS);
  return Solution{
      .score = m < M ? m : M + 2 * N * M - static_cast<int>(turns.size()),
      .turns = std::move(turns),
  };
}

}  // namespace tbrekalo

namespace tb = tbrekalo;

auto main(int, char**) -> int {
  auto coords = tb::Load(std::cin);
  auto solution = tb::Solve(coords, false);
  if (auto altered = tb::Solve(coords, true); altered.score > solution.score) {
    solution = std::move(altered);
  }
  std::cerr << "elapsed=" << tb::elapsed() << std::endl;
  tb::Print(std::cout, solution);
  return 0;
}
