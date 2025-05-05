#include <algorithm>
#include <array>
#include <cassert>
#include <chrono>
#include <cmath>
#include <compare>
#include <iostream>
#include <limits>
#include <optional>
#include <ostream>
#include <random>
#include <span>
#include <type_traits>
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

}  // namespace tbrekalo

namespace tbrekalo {

static inline constexpr int N = 20;
static inline constexpr int M = 40;
static inline constexpr int MAX_TURNS = 2 * N * M;

static constexpr auto MAX_SCORE = 2 * N * M + M;
static constexpr auto BEGIN_TEMP = MAX_SCORE;
static constexpr auto END_TEMP = N;

static auto temperature() -> double {
  return BEGIN_TEMP *
         std::pow(1. * END_TEMP / BEGIN_TEMP, 1. * elapsed() / MAX_RUNTIME);
}

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
  return dst;
}
[[gnu::always_inline]] constexpr auto operator-(Coord lhs, Coord rhs) noexcept
    -> Coord {
  Coord dst{.row = lhs.row - rhs.row, .col = lhs.col - rhs.col};
  return dst;
}

using Coords = std::array<Coord, M>;
using CoordsView = std::span<Coord const, M>;

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

// prefix increment; ++dir
[[gnu::always_inline]] inline auto operator++(Dir& dir) -> Dir& {
  dir = Dir(std::to_underlying(dir) + 1);
  return dir;
}

// postfix increment; dir++
[[gnu::always_inline]] inline auto operator++(Dir& dir, int) -> Dir {
  auto ret = dir;
  ++dir;
  return ret;
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

static auto FromToDir(Coord from, Coord to) noexcept -> std::optional<Dir> {
  if (from.row > to.row) {
    return Dir::U;
  }
  if (from.row < to.row) {
    return Dir::D;
  }
  if (from.col > to.col) {
    return Dir::L;
  }
  if (from.col < to.col) {
    return Dir::R;
  }

  return std::nullopt;
}

[[gnu::always_inline]] inline constexpr auto RowDif(Coord from,
                                                    Coord to) noexcept -> int {
  return from.row < to.row ? to.row - from.row : from.row - to.row;
}

[[gnu::always_inline]] inline constexpr auto ColDif(Coord from,
                                                    Coord to) noexcept -> int {
  return from.col < to.col ? to.col - from.col : from.col - to.col;
}

static auto FromToDif(Coord from, Coord to) noexcept -> int {
  auto opt_dir = FromToDir(from, to);
  if (!opt_dir) {
    return 0;
  }

  int ret = 0;
  switch (*opt_dir) {
    case Dir::U:
      ret = from.row - to.row;
      break;
    case Dir::D:
      ret = to.row - from.row;
      break;
    case Dir::L:
      ret = from.col - to.col;
      break;
    case Dir::R:
      ret = to.col - from.col;
      break;
  }
  assert(0 < ret && ret <= N);
  return ret;
}

enum class Cell : char { EMPTY, TARGET, BLOCK };

template <class T>
using GridT = std::array<std::array<T, N>, N>;

using GridCell = GridT<Cell>;
using GridProb = GridT<double>;

template <class T>
static constexpr auto FillGrid(T value) -> GridT<T> {
  GridT<T> dst;
  std::ranges::for_each(dst, [value](GridT<T>::reference row) -> void {
    std::ranges::fill(row, value);
  });

  return dst;
}

template <class Fn>
static constexpr auto GenerateGrid(Fn&& fn) -> GridT<std::invoke_result_t<Fn>> {
  using Grid = GridT<std::invoke_result_t<Fn>>;

  Grid dst;
  std::ranges::for_each(
      dst, [fn = std::forward<Fn>(fn)](Grid::reference row) -> void {
        std::ranges::generate(row, fn);
      });

  return dst;
}

static auto FuzzProbs(GridProb grid) -> GridProb {
  std::normal_distribution<> normal(0.1, 0.0025);
  for (int r = 0; r < N; ++r) {
    for (int c = 0; c < N; ++c) {
      grid[r][c] = std::clamp(grid[r][c] + normal(RNG), 0.0025, 0.95);
    }
  }
  return grid;
}

static constexpr auto CreateProbGrid(float value) -> GridProb {
  GridProb dst;
  std::ranges::for_each(dst, [value](GridProb::reference row) -> void {
    std::ranges::fill(row, value);
  });
  return dst;
};

struct Solution {
  int score;
  std::vector<Turn> turns;
  GridProb block_probs;
};

static auto Load(std::istream& istrm) -> Coords {
  int _;
  Coords dst;
  istrm >> _ >> _;
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

static auto Solve(CoordsView coords, GridProb p) -> Solution {
  Coord cur = coords[0];
  std::vector<Turn> turns;
  auto can_keep_turning = [&turns] [[gnu::always_inline]] -> bool {
    return turns.size() < MAX_TURNS;
  };
  auto try_turn = [&turns, &can_keep_turning](Turn turn) -> bool {
    if (can_keep_turning()) {
      turns.push_back(turn);
      return true;
    }
    return false;
  };

  GridProb adjusted_p = FillGrid(0.);
  GridCell grid = FillGrid(Cell::EMPTY);
  for (int i = 1; i < coords.size(); ++i) {
    grid[coords[i].row][coords[i].col] = Cell::TARGET;
  }

  auto alter_block = [&grid](Coord target) -> void {
    auto& cell = grid[target.row][target.col];
    assert(cell != Cell::TARGET);
    cell = cell == Cell::EMPTY ? Cell::BLOCK : Cell::EMPTY;
  };

  auto try_block = [p, &grid, alter_block](
                       Coord base, Dir banned_dir) -> std::vector<Turn> {
    std::vector<Turn> dst;
    std::uniform_real_distribution<float> distr(0, 1);
    for (auto dir = Dir::U; dir < Dir::R; ++dir) {
      if (dir == banned_dir) {
        continue;
      }
      if (auto target = base + DirToDif(dir);
          InBounds(target) && grid[target.row][target.col] == Cell::EMPTY &&
          distr(RNG) < p[target.row][target.col]) {
        dst.push_back(Turn{.action = Action::A, .dir = dir});
        alter_block(target);
      }
    }
    return dst;
  };

  auto find_next_to_block = [&grid](Coord pos,
                                    Dir dir) -> std::pair<Coord, Coord> {
    auto nxt = pos + DirToDif(dir);
    while (true) {
      nxt = pos + DirToDif(dir);
      if (!InBounds(nxt) || grid[nxt.row][nxt.col] == Cell::BLOCK) {
        break;
      }
      pos = nxt;
    }
    return {pos, nxt};
  };

  auto calc_slide_cost = [](Coord from, Coord to, Coord block) -> int {
    auto block_to_dist = FromToDif(block, to);
    return block_to_dist;
  };

  auto dir_move = [&adjusted_p, &try_turn, &grid, &alter_block, &try_block,
                   &find_next_to_block,
                   &calc_slide_cost](Coord cur, Coord to) -> Coord {
    auto dir = FromToDir(cur, to);
    if (!dir) {
      return cur;
    }

    auto move_dif = FromToDif(cur, to);
    auto [next_to_block, block] = find_next_to_block(cur, *dir);
    int slide_dif = calc_slide_cost(cur, to, next_to_block);
    if (slide_dif < move_dif) {
      if (InBounds(block)) {
        adjusted_p[block.row][block.col] = 0.95;
      }
      if (try_turn(Turn{.action = Action::S, .dir = *dir})) {
        cur = next_to_block;
        move_dif = FromToDif(cur, to);
        auto opt_dir = FromToDir(cur, to);
        if (!opt_dir) {
          return cur;
        }
        dir = *opt_dir;
      }
    }

    for (int i = 0; i < move_dif; ++i) {
      for (auto block_turn : try_block(cur, *dir)) {
        if (!try_turn(block_turn)) {
          break;
        }
      }
      auto nxt = cur + DirToDif(*dir);
      if (grid[nxt.row][nxt.col] == Cell::BLOCK) {
        if (!try_turn(Turn{.action = Action::A, .dir = *dir})) {
          break;
        }
        alter_block(nxt);
      }
      if (!try_turn(Turn{.action = Action::M, .dir = *dir})) {
        break;
      }
      cur = nxt;
    }

    return cur;
  };

  int m = 1;
  for (; m < M && can_keep_turning(); ++m) {
    while (cur != coords[m] && can_keep_turning()) {
      cur = dir_move(cur, coords[m]);
    }
  }

  assert(turns.size() <= MAX_TURNS);
  return Solution{
      .score = m < M ? m : M + 2 * N * M - static_cast<int>(turns.size()),
      .turns = std::move(turns),
      .block_probs = adjusted_p,
  };
}

}  // namespace tbrekalo

namespace tb = tbrekalo;

auto main(int, char**) -> int {
  using namespace std::literals;
  std::uniform_real_distribution<> p(0, 1);
  auto last_report = tb::elapsed();

  auto coords = tb::Load(std::cin);
  auto best = tb::Solve(coords, tb::FillGrid(0.001));
  for (auto current = best; tb::elapsed() < tb::MAX_RUNTIME;) {
    auto candidate = tb::Solve(coords, tb::FuzzProbs(current.block_probs));
    if (auto elapsed = tb::elapsed(); elapsed - last_report > 100ms) {
      std::cerr << "elapsed=" << elapsed << " best-score=" << best.score
                << " current-score" << current.score
                << " candidate-score=" << candidate.score << std::endl;
      last_report = elapsed;
    }

    auto acc_p = std::exp(static_cast<double>(candidate.score - current.score) /
                          tb::temperature());
    if (candidate.score > current.score || p(tb::RNG) < acc_p) {
      current = std::move(candidate);
    }

    if (current.score > best.score) {
      best = current;
    }
  }
  tb::Print(std::cout, best);
  return 0;
}
