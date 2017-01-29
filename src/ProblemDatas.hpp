// Copyright 2017 <Biagio Festa>
#ifndef __FOR_CH_PROBLEM_DATAS__HPP
#define __FOR_CH_PROBLEM_DATAS__HPP
#include <vector>
#include <set>
#include <string>
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

  HashMap<EdgeIndex, EdgeLink> m_mapEdge_index2link;
  HashMapPair<EdgeLink, EdgeIndex> m_mapEdge_link2index;
  HashMap<VertexIndex, std::set<EdgeIndex>> m_map_vertex2outedges;
  HashMap<VertexIndex, std::set<EdgeIndex>> m_map_vertex2inedges;

  template<typename T>
  static std::vector<T> parse_vector_dat(std::string vector_data);
  void compute_all_distances();
  void compute_edges_indices();

  static void trim_str(std::string* str);
};

template<typename T>
std::vector<T> ProblemDatas::parse_vector_dat(std::string vector_data) {
  std::vector<T> rtn;
  trim_str(&vector_data);
  while (vector_data.size()) {
    std::string line = vector_data.substr(0, vector_data.find('\n'));
    vector_data.erase(0, line.size() + 1);
    trim_str(&line);
    if (line.size()) {
      std::string index = line.substr(0, line.find(' '));
      std::string x = line.substr(line.find(' '));
      trim_str(&index);
      trim_str(&x);
      unsigned nindex = std::stoi(index);
      int nx = std::stoi(x);
      assert(nindex == rtn.size());
      rtn.push_back(nx);
    }
  }
  return rtn;
}

}  // namespace for_ch


#endif  // __FOR_CH_PROBLEM_DATAS__HPP
