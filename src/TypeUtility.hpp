// Copyright 2017 <Biagio Festa>
#ifndef __FOR_CH_TYPE_UTILITY_HPP
#define __FOR_CH_TYPE_UTILITY_HPP

#include <vector>
#include <utility>
#include <map>
#include <cassert>
#include <unordered_map>

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

#endif  // __FOR_CH_TYPE_UTILITY_HPP
