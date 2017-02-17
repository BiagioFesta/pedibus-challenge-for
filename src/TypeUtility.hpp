// Copyright 2017 <Biagio Festa>
#ifndef __FOR_CH_TYPE_UTILITY_HPP
#define __FOR_CH_TYPE_UTILITY_HPP

#include <vector>
#include <utility>
#include <map>
#include <cassert>
#include <unordered_map>
#include <algorithm>

using VertexIndex = unsigned;

using EdgeIndex = unsigned;
using EdgeLink = std::pair<VertexIndex, VertexIndex>;

using RealNumber = float;

template<typename T>
using RowMatrix = std::vector<T>;
template<typename T>
using Matrix = std::vector<RowMatrix<T>>;

using RowMatrixReal = RowMatrix<RealNumber>;
using MatrixReal = Matrix<RealNumber>;

struct HashPair {
  template<typename T1, typename T2>
  std::size_t operator()(const std::pair<T1, T2>& p) const noexcept {
    return (std::hash<T1>()(p.first) ^ std::hash<T2>()(p.second));
  }
};

template<typename Key, typename T>
using Map = std::map<Key, T>;

template<typename Key, typename T, typename Hash = std::hash<Key>>
using HashMap = std::unordered_map<Key, T, Hash>;

template<typename Key, typename T>
using HashMapPair = HashMap<Key, T, HashPair>;

static constexpr VertexIndex SCHOOL_INDEX = 0;

namespace for_ch {
class ProblemDatas;
}

struct Solution {
  std::vector<bool> m_active_edges;
  int m_num_leaves = -1;
  RealNumber m_danger = -1;

  bool operator==(const Solution& oth) const noexcept {
    return m_active_edges == oth.m_active_edges;
  }

  bool operator!=(const Solution& oth) const noexcept {
    return m_active_edges != oth.m_active_edges;
  }

  bool operator<(const Solution& oth) const noexcept {
#ifndef NDEBUG
    if (m_active_edges < oth.m_active_edges) {
      assert(*this != oth);
    }
#endif
    return m_active_edges < oth.m_active_edges;
  }

  bool compute_feasibility(const for_ch::ProblemDatas& problem) const noexcept;

  template<typename T>
  static void AddSet(const std::set<T>& source,
                     std::set<T>* destination) noexcept;
};


template<typename T>
void Solution::AddSet(const std::set<T>& source,
                      std::set<T>* destination) noexcept {
  std::for_each(source.cbegin(),
                source.cend(),
                [&destination] (const T& e) {
                  destination->insert(e);
                });
}


#endif  // __FOR_CH_TYPE_UTILITY_HPP
