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

}  // namespace tbrekalo

namespace tbrekalo{}

namespace tb = tbrekalo;

auto main(int, char**) -> int { return 0; }
