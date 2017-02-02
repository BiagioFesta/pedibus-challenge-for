// Copyright 2017 <Biagio Festa>
#ifndef __FOR_CH_PROBLEM_DATAS__HPP
#define __FOR_CH_PROBLEM_DATAS__HPP
#include <vector>
#include <set>
#include <string>
#include <sstream>
#include "TypeUtility.hpp"

namespace for_ch {

class ProblemDatas {
 public:
  void parse_problem_dat(const std::string& filename);
 public:
  unsigned m_numNodes;
  unsigned m_numEdges;
  RealNumber m_alpha;
  std::vector<int> m_coordX_vertices;
  std::vector<int> m_coordY_vertices;
  MatrixReal m_distance_matrix;
  std::vector<RealNumber> m_max_distances;
  std::vector<RealNumber> m_distances_FromNearest;
  MatrixReal m_dangerousness_matrix;

  HashMap<EdgeIndex, EdgeLink> m_mapEdge_index2link;
  HashMapPair<EdgeLink, EdgeIndex> m_mapEdge_link2index;
  HashMap<VertexIndex, std::set<EdgeIndex>> m_map_vertex2outedges;
  HashMap<VertexIndex, std::set<EdgeIndex>> m_map_vertex2inedges;

  template<typename T>
  static std::vector<T> parse_vector_dat(std::string vector_data);
  template<typename T>
  static Matrix<T> parse_matrix_dat(std::string matrix_data);

  void compute_all_distances();
  void compute_edges_indices();

  static void trim_str(std::string* str);
};

template<typename T>
std::vector<T> ProblemDatas::parse_vector_dat(std::string vector_data) {
  std::vector<T> rtn;
  std::map<unsigned, T> map_vector;
  std::istringstream sstream(vector_data);

  unsigned index, value;
  while (sstream.eof() == false) {
    sstream >> index;
    sstream >> value;
    map_vector[index] = value;
  }

#ifndef NDEBUG
  index = 0;
  for (const auto& p : map_vector) {
    assert(p.first == index++);
  }
#endif

  rtn.resize(map_vector.size());
  for (const auto& p : map_vector) {
    rtn[p.first] = p.second;
  }

  return rtn;
}

template<typename T>
Matrix<T> ProblemDatas::parse_matrix_dat(std::string matrix_data) {
  static constexpr unsigned SIZE_BUFFER = 1024 * 1024;
  Matrix<T> rtn;
  std::map<unsigned, RowMatrix<T>> map_matrix;
  std::istringstream sstream(matrix_data);

  unsigned index_row;
  std::string line;
  T value;
  char buffer[SIZE_BUFFER];

  while (sstream.eof() == false) {
    sstream >> index_row;
    sstream.getline(buffer, SIZE_BUFFER, '\n');
    line = std::string(buffer, sstream.gcount() - 1);
    trim_str(&line);

    std::istringstream linestream(line);
    RowMatrix<T> row;
    while (linestream.eof() == false) {
      linestream >> value;
      row.push_back(value);
    }

    map_matrix[index_row] = std::move(row);
  }

#ifndef NDEBUG
  index_row = 0;
  for (const auto& p : map_matrix) {
    assert(p.first == index_row++);
  }
#endif

  rtn.resize(map_matrix.size());
  for (const auto& p : map_matrix) {
    rtn[p.first] = std::move(p.second);
  }

  return rtn;
}

}  // namespace for_ch


#endif  // __FOR_CH_PROBLEM_DATAS__HPP
