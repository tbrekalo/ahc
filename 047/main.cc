#include <algorithm>
#include <array>
#include <cassert>
#include <chrono>
#include <cmath>
#include <functional>
#include <iostream>
#include <limits>
#include <memory>
#include <numeric>
#include <random>
#include <ranges>
#include <span>
#include <unordered_set>
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

static constexpr auto N = 36;
static constexpr auto M = 12;
static constexpr auto L = 1'000'000;
static constexpr auto S = 12;
static constexpr auto ALPHA = 6;

struct Word {
  std::string value;
  int preference;
};

using TransProb = std::array<char, M>;

struct Row {
  char label;
  TransProb value;
};

using Model = std::array<Row, M>;

struct K3Mer {
  char a;
  char b;
  char c;

  friend inline constexpr auto operator<=>(K3Mer lhs,
                                           K3Mer rhs) noexcept = default;
};

struct K3MerHash {
  inline constexpr auto operator()(K3Mer kmer) const noexcept -> std::size_t {
    std::size_t hash = 5381;
    auto arr = std::bit_cast<std::array<char, 3>>(kmer);
    for (auto it : arr) {
      hash = ((hash << 5) + hash) + it;
    }

    return hash;
  }
};

class Trie {
  struct Node {
    Node(char value) : value(value) { std::ranges::fill(children, nullptr); };

    char value;
    std::array<std::unique_ptr<Node>, ALPHA + 1> children;
  };

  struct SearchState {
    std::string value;
    Node* node;
  };

  std::shared_ptr<Node> root_;
  std::unordered_set<std::string> words_;
  std::vector<SearchState> search_state_;

  auto embed(std::string_view word) -> void {
    auto node = root_.get();
    for (auto c :
         std::views::transform(word, [](char c) -> int { return c - 'a'; })) {
      if (node->children[c] == nullptr) {
        node->children[c] = std::make_unique<Node>(c);
      }
      node = node->children[c].get();
    }

    if (node->value != '$') {
      node->children[ALPHA] = std::make_unique<Node>('$');
    }
  }

 public:
  Trie(std::span<Word> words) : root_(std::make_shared<Node>('$')) {
    search_state_.push_back({"", root_.get()});
    for (auto const& [word, _] : words) {
      words_.insert(word);
      embed(word);
    }
  }

  auto search(char c) {
    std::vector<SearchState> next_state;
    for (auto [str, node] : search_state_) {
      if (node->children[c] != nullptr) {
      }
    }
  }

  auto clear() -> void { search_state_.clear(); }
};

auto Load(std::istream& istrm) -> std::vector<Word> {
  int _;
  std::vector<Word> dst(N);
  istrm >> _ >> _ >> _;
  for (int i = 0; i < N; ++i) {
    istrm >> dst[i].value >> dst[i].preference;
  }

  std::ranges::sort(dst, std::greater<>{}, &Word::preference);
  return dst;
}

auto Print(std::ostream& ostrm, Model model) -> void {
  for (auto [label, value] : model) {
    assert(std::accumulate(value.begin(), value.end(), 0) == 100);
    ostrm << label << ' ';
    for (auto it : value) {
      ostrm << static_cast<int>(it) << ' ';
    }
    ostrm << '\n';
  }
  std::flush(ostrm);
}

auto GenerateRandomTrans() -> TransProb {
  TransProb dst;
  std::ranges::fill(dst, 0);

  for (int i = 100, j = 0; i > 0; j = (j + 1) % M) {
    auto x = std::uniform_int_distribution<char>(0, i)(RNG);
    dst[j] += x;
    i -= x;
  }

  return dst;
}

auto CreateEmptyModel() -> Model {
  Model dst;
  for (int i = 0; i < M; ++i) {
    dst[i].label = 'a' + (i % ALPHA);
    std::ranges::fill(dst[i].value, 0);
  }

  return dst;
}

auto GenerateInitial() -> Model {
  Model dst = CreateEmptyModel();
  for (auto& row : dst) {
    row.value = GenerateRandomTrans();
  }

  return dst;
}

auto GenerateSequence(Model model) {
  std::string word;
  std::vector<std::discrete_distribution<int>> distr{M};
  for (int i = 0; i < M; ++i) {
    distr[i] = std::discrete_distribution<int>(model[i].value.begin(),
                                               model[i].value.end());
  }

  for (auto s = 0, i = 0; i < L; ++i) {
    word.push_back(model[s].label);
    s = distr[s](RNG);
  }

  return word;
}

auto EvaluateModel(Model model, std::span<Word const> words, int n_samples) {
  for (auto [value, preference] : words) {
  }
}

auto CreateK3MerUniformModel(std::span<Word const> words, int k) -> Model {
  assert(k <= words.size());
  Model model = CreateEmptyModel();
  std::array<std::array<std::uint64_t, M>, M> score;
  for (auto& row : score) {
    std::ranges::fill(row, 0);
  }

  for (int i = 0; i < k; ++i) {
    auto const& [value, preference] = words[i];
    for (int j = 2; j < value.size(); ++j) {
      score[value[j - 2] - 'a'][value[j - 1] - 'a'] += preference;
      score[ALPHA + value[j - 1] - 'a'][ALPHA + value[j] - 'a'] += preference;
    }
  }
  for (int i = 0; i < M; ++i) {
    int prob_sum = 0;
    std::ranges::fill(model[i].value, 0);
    auto freq_sum = std::accumulate(score[i].begin(), score[i].end(), 0);
    if (freq_sum == 0) {
      model[i].value = GenerateRandomTrans();
      continue;
    }

    for (int j = 0; j < M; ++j) {
      int value = 100 * score[i][j] / freq_sum;
      assert(prob_sum + value <= 100);

      model[i].value[j] = value;
      prob_sum += value;
    }

    for (int j = 0; prob_sum < 100; j = (j + 1) % M, ++prob_sum) {
      ++model[i].value[j];
    }
  }

  return model;
}

}  // namespace tbrekalo

namespace tb = tbrekalo;

auto main(int, char**) -> int {
  auto words = tb::Load(std::cin);
  auto trie = tb::Trie(words);

  auto model = tb::CreateK3MerUniformModel(words, tb::N);
  auto sequence = tb::GenerateSequence(model);

  tb::Print(std::cout, model);
  std::cerr << "elapsed=" << tb::elapsed() << std::endl;
  return 0;
}
