// Copyright 2017 <Biagio Festa>
#ifndef __FOR_CH_TYPE_UTILITY_HPP
#define __FOR_CH_TYPE_UTILITY_HPP

#include <vector>
#include <utility>
#include <map>
#include <cassert>

using VertexIndex = unsigned;
using VectVertices = std::vector<VertexIndex>;

using EdgeIndex = unsigned;
using VectEdges = std::vector<EdgeIndex>;
using EdgeLink = std::pair<VertexIndex, VertexIndex>;

using RealNumber = double;

template<typename T>
using RowMatrix = std::vector<T>;
template<typename T>
using Matrix = std::vector<RowMatrix<T>>;

using RowMatrixReal = RowMatrix<RealNumber>;
using MatrixReal = Matrix<RealNumber>;

template<typename Key, typename T>
using Map = std::map<Key, T>;

static constexpr VertexIndex SCHOOL_INDEX = 0;

#endif  // __FOR_CH_TYPE_UTILITY_HPP
