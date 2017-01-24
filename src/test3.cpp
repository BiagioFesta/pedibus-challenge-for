// Copyright 2017 <Biagio Festa>
#include <boost/algorithm/string.hpp>
#include <iostream>
#include <string>
#include <vector>
#include <fstream>
#include <utility>
#include <map>
#include <cmath>
#include <algorithm>
#include <set>
#include <iterator>
#include <stack>
#include <queue>

class AlgorithmSolver {
 public:
  using RealNumber = double;
  using VertexType = unsigned;
  using State = VertexType;
  using Action = VertexType;
  template<typename T>
  using RowMatrix = std::vector<T>;
  template<typename T>
  using Matrix = std::vector<RowMatrix<T>>;
  using MatrixR = Matrix<RealNumber>;
  using VectVertices = std::vector<VertexType>;
  static const VertexType SCHOOL_NODE = 0;

  struct Path {
    VectVertices m_pathVertices;
    std::vector<RealNumber> m_routeLenPerLevel;
  };

  struct Heuristic {
    explicit Heuristic(const AlgorithmSolver& problem) :
        mp_problem(&problem) {
    }

    RealNumber operator()(const VertexType& v1,
                          const VertexType& v2) const {
      using PairIndices = std::pair<unsigned, unsigned>;
      // This return an heuristic value which says how much
      // v2 is good respect to v1

      const auto& distanceMatrix = mp_problem->m_distanceMatrix;
      const auto& nearestDistance = mp_problem->m_distanceFromNearest;

      // Check whether v2 is in a existent path
      bool v2_in_a_path_as_leaf = false;
      bool v2_is_a_leaf = false;
      PairIndices indices;
      if (mp_problem->checkVertexIsInAPath(v2, &indices)) {
        if (indices.second != 0) {
          v2_in_a_path_as_leaf = true;
        } else {
          v2_is_a_leaf = true;
        }
      }

      if (v2_is_a_leaf == true) {
        return 0;
      }

      return ((m_A * nearestDistance[v2]) +
              (m_B * distanceMatrix[v1][v2]) +
              (m_C * v2_in_a_path_as_leaf));
    }

    // m_A factor scale is the grade of isolation of a point
    RealNumber m_A = 0.1;

    // m_B factor scale of distance from the evaluator point
    RealNumber m_B = 1.0;

    // m_C factor scale if v2 is already in a path (we don't want to link
    // with a internal node in a existent path BUT is not a leaf
    RealNumber m_C = 100;

    const AlgorithmSolver* mp_problem;
  };

  template<typename HeuristicT>
  struct CompareVertices {
    explicit CompareVertices(const VertexType& v, const HeuristicT& h):
        mp_v(&v), mp_h(&h) {
    }

    bool operator()(const VertexType& v1, const VertexType& v2) const {
      return (*mp_h)(*mp_v, v1) >= (*mp_h)(*mp_v, v2);
    }

    const VertexType* mp_v;
    const HeuristicT* mp_h;
  };

  template<typename HeuristicT>
  bool launchProblemSolver();

  void parse_problem_dat(const std::string& filename);
  void compute_all_distances();

  template<typename HeuristicT>
  bool construct_path();

  template<typename T>
  static std::vector<T> parse_vector_dat(std::string vector_data);

  /// This method takes a path and returns all possibile vertices
  /// you can add to the input path. The vertices you can apply to the path
  /// satisfy all distance constraints
  /// @param [in] current_path   The path you want to know which vertices you
  ///                            can add.
  /// @param [out] nexts         A vector of vertices which contains all
  ///                            vertices you can properly add to that path.
  void find_possible_next_nodes(const Path& current_path,
                                VectVertices* nexts);

  /// This method takes a vertex and adds it to a path
  /// @param [in] v                 The vertex you want to add to the path
  /// @param [in.out] current_path  The path you want to modify
  void apply_next_node_to_path(const VertexType& v, Path* path);

  template<typename HeuristicT>
  void sort_heuristic_accordance(const VertexType& lastInPath,
                                 VectVertices* vertices);

  void write_on_csv_coordinate(std::ostream* out_stream);

  /// Check whether a vertex is already in a path
  ///
  /// @param [out]    ....
  /// @return true whether is already in a path, false otherwise
  ///
  bool checkVertexIsInAPath(const VertexType& v,
                            std::pair<unsigned, unsigned>* outIndices) const;

  /// The number of nodes in the problem
  unsigned m_numNodes;

  /// The alpha problem parameter
  RealNumber m_alpha;

  /// For all nodes in the problem, its X coordinate
  std::vector<int> m_coordX_vertices;

  /// For all nodes in the problem, its Y coordinae
  std::vector<int> m_coordY_vertices;

  /// The problem distance matrix.
  /// E.g. m_distanceMatrix[i][j] returns the distance from i to j
  MatrixR m_distanceMatrix;

  /// For all nodes in the problem, the max distance it can walk
  std::vector<RealNumber> m_maxDistances;

  /// For all nodes in the problem, the distance from the nearest other node
  std::vector<RealNumber> m_distanceFromNearest;

  /// The set of vertices which have been linked in some path
  std::set<VertexType> m_linkedVertices;

  /// All paths founds
  std::vector<Path> m_foundPaths;
};

template<typename HeuristicT>
bool AlgorithmSolver::launchProblemSolver() {
  m_foundPaths.clear();
  m_linkedVertices.clear();

  while (m_linkedVertices.size() != m_numNodes) {
    if (construct_path<HeuristicT>() == false) {
      // if the algorithm cannot construct a path
      // that means there is a problem
      // Not feasible solution found
      return false;
    }
  }

  return true;
}

bool AlgorithmSolver::checkVertexIsInAPath(
    const VertexType& v, std::pair<unsigned, unsigned>* outIndices) const {
  assert(outIndices != nullptr);

  // First of all check if v is linked
  if (m_linkedVertices.find(v) != m_linkedVertices.cend()) {
    for (unsigned i = 0; i < m_foundPaths.size(); ++i) {
      const Path& p = m_foundPaths[i];
      for (unsigned j = 0; j < p.m_pathVertices.size(); ++j) {
        if (p.m_pathVertices[j] == v) {
          *outIndices = std::make_pair(i, j);
          return true;
        }
      }
    }
  }
  return false;
}

template<typename HeuristicT>
bool AlgorithmSolver::construct_path() {
  using OpenList = std::stack<Path>;

  // Get the list of free vertices (not already linked)
  std::vector<VertexType> freeVertices;
  std::vector<RealNumber> distancesFree;
  for (unsigned i = 0; i < m_numNodes; ++i) {
    if (m_linkedVertices.find(i) == m_linkedVertices.cend()) {
      freeVertices.push_back(i);
      distancesFree.push_back(m_distanceMatrix[i][SCHOOL_NODE]);
    }
  }

  // Among the free nodes, find what is most distant from school
  const auto farestDistance = std::max_element(distancesFree.cbegin(),
                                               distancesFree.cend());
  // Get the index of the node farest from school
  const auto farestIndex = std::distance(distancesFree.cbegin(),
                                         farestDistance);
  // Get the vertex most distant from school
  const auto& farestNode = freeVertices[farestIndex];

  // Some local variables
  std::pair<unsigned, unsigned> indices;
  VectVertices nexts;
  OpenList open;

  // Load into the openlist the initial path
  Path initialPath;
  initialPath.m_pathVertices.push_back(farestNode);
  initialPath.m_routeLenPerLevel.push_back(0);
  open.push(std::move(initialPath));

  while (open.empty() == false) {
    const Path current = std::move(open.top());
    open.pop();

    // If the last node in the current path is the school
    // a path has been found.
    if (current.m_pathVertices.back() == SCHOOL_NODE) {
      // All vertices in the path have to be added to
      // the list of linked vertices
      for (const auto& v : current.m_pathVertices) {
        m_linkedVertices.insert(v);
      }

      // Insert the path into the member list
      m_foundPaths.push_back(std::move(current));
      return true;
    }

    // Find all possibile nodes in which the current path can
    // continue into
    find_possible_next_nodes(current, &nexts);

    // All possibile choices must to be sorted in according to the
    // heuristic
    sort_heuristic_accordance<HeuristicT>(current.m_pathVertices.back(),
                                          &nexts);

    // Test all possible choices
    for (const auto& n : nexts) {
      // If n is already in a path the current path must to be linked
      if (checkVertexIsInAPath(n, &indices)) {
        if (indices.second == 0) {
          // In this case the vertex n is a leaf of an existent path
          Path& oldPath = m_foundPaths[indices.first];
          const RealNumber& lenPath = oldPath.m_routeLenPerLevel[0];
          const VertexType& leafPath = oldPath.m_pathVertices[0];
          const VertexType& lastCurrent = current.m_pathVertices.back();
          const RealNumber& lenLink = m_distanceMatrix[lastCurrent][leafPath];

          bool admissible = true;
          for (unsigned i = 0;
               i < current.m_pathVertices.size() && admissible == true;
               ++i) {
            RealNumber newDistance = current.m_routeLenPerLevel[i] + lenLink +
                lenPath;
            if (newDistance > m_maxDistances[current.m_pathVertices[i]]) {
              admissible = false;
            }
          }

          if (admissible == true) {
            std::for_each(   // for each vertices in the current path
                current.m_pathVertices.crbegin(),
                current.m_pathVertices.crend(),
                [this, &oldPath, &lenPath]
                (const VertexType& v) {
                  RealNumber newDistance =
                      m_distanceMatrix[v][oldPath.m_pathVertices.front()] +
                      lenPath;

                  // insert in head position the vertex in the oldPath
                  oldPath.m_pathVertices.insert(
                      oldPath.m_pathVertices.begin(), v);

                  // insert the routeLen in head
                  oldPath.m_routeLenPerLevel.insert(
                      oldPath.m_routeLenPerLevel.begin(), newDistance);

                  // insert the vertex 'v' in the linked nodes
                  m_linkedVertices.insert(v);
                });
            return true;
          }   // if admissible == true
        } else {
          // In this case the vertex n is a inner node of an existent path
          // Note that inner node can be also zero (the SCHOOL NODE)
          const VertexType& lastCurrent = current.m_pathVertices.back();
          const Path & oldPath = m_foundPaths[indices.first];
          const RealNumber& lenPath =
              oldPath.m_routeLenPerLevel[indices.second];
          const VertexType& linkRef = oldPath.m_pathVertices[indices.second];
          const RealNumber& lenLink = m_distanceMatrix[lastCurrent][linkRef];

          // Create new nodes as the list from linkRef to end (of old Path)
          VectVertices newPathNodes;
          std::copy(oldPath.m_pathVertices.cbegin() + indices.second,
                    oldPath.m_pathVertices.cend(),
                    std::back_inserter(newPathNodes));
          std::vector<RealNumber> newRouteDistance;
          std::copy(oldPath.m_routeLenPerLevel.cbegin() + indices.second,
                    oldPath.m_routeLenPerLevel.cend(),
                    std::back_inserter(newRouteDistance));

          // Check if admissible
          bool admissible = true;
          for (unsigned i = 0;
               i < current.m_pathVertices.size() && admissible == true;
               ++i) {
            RealNumber newDistance = current.m_routeLenPerLevel[i] + lenLink +
                lenPath;
            if (newDistance > m_maxDistances[current.m_pathVertices[i]]) {
              admissible = false;
            }
          }

          if (admissible == true) {
            // Create a new path
            Path newPath = current;
            newPath.m_pathVertices.insert(newPath.m_pathVertices.cend(),
                                          newPathNodes.cbegin(),
                                          newPathNodes.cend());
            const auto itNew = newPath.m_routeLenPerLevel.insert(
                newPath.m_routeLenPerLevel.cend(),
                newRouteDistance.cbegin(),
                newRouteDistance.cend());
            std::for_each(newPath.m_routeLenPerLevel.begin(),
                          itNew,
                          [&lenLink, &lenPath] (RealNumber& n) {
                            n = n + lenLink + lenPath;
                          });
            open.push(std::move(newPath));
          }   // if admissible == true
        }
      } else {
        // In this case the node n is free and can be easly
        // added to the current path
        Path newPath = current;
        apply_next_node_to_path(n, &newPath);
        open.push(std::move(newPath));
      }
    }  // for all n in next
  }  // while openList is not empty

  // Is the open list is empty and no path has been found,
  // return false
  return false;
}

template<typename HeuristicT>
void AlgorithmSolver::sort_heuristic_accordance(
    const VertexType& lastInPath, VectVertices* vertices) {
  assert(vertices != nullptr);
  HeuristicT h(*this);
  CompareVertices<HeuristicT> comp(lastInPath, h);
  std::sort(vertices->begin(), vertices->end(), comp);
}

void AlgorithmSolver::apply_next_node_to_path(const VertexType& v,
                                              Path* path) {
  assert(path != nullptr);

  // Assert the node you want to add is not already in the path
  assert(std::find(path->m_pathVertices.cbegin(),
                   path->m_pathVertices.cend(),
                   v) == path->m_pathVertices.cend());

  // Takes the distance from the last in the current path to v
  const RealNumber& distance_last_to_i =
      m_distanceMatrix[path->m_pathVertices.back()][v];

  // Add v as last vertices in the path
  path->m_pathVertices.push_back(v);

  // Now you have to add the walk distance to all nodes in the
  // path
  for (unsigned i = 0;
       i < path->m_routeLenPerLevel.size();
       ++i) {
    path->m_routeLenPerLevel[i] += distance_last_to_i;

    // Assertion the new insertion does not compromize the
    // integrity of the path for no nodes
    assert(path->m_routeLenPerLevel[i] <=
           m_maxDistances[path->m_pathVertices[i]]);
  }

  // For the last added in the path (v) the walk distance is zeros
  path->m_routeLenPerLevel.push_back(0);
}

void AlgorithmSolver::write_on_csv_coordinate(std::ostream* out_stream) {
  assert(out_stream);
  *out_stream << "IndexNode, X, Y";
  for (unsigned i = 0; i < m_numNodes; ++i) {
    *out_stream << "\n" << i << ", " <<
        m_coordX_vertices[i] << ", " <<
        m_coordY_vertices[i];
  }
}

void AlgorithmSolver::find_possible_next_nodes(const Path& current_path,
                                               VectVertices* nexts) {
  /// This algorithm give a node should return all possibile next
  /// nodes.

  /// Clean the output
  nexts->clear();

  const unsigned LEN_PATH = current_path.m_pathVertices.size();

  /// If current_path is empty, then you can select all except the school
  if (LEN_PATH == 0) {
    for (unsigned i = 0; i < m_numNodes; ++i) {
      if (i != SCHOOL_NODE) {
        nexts->push_back(i);
      }
    }
  } else {
    /// In that case there is path

    // Get the last of the path
    const auto& last_node = current_path.m_pathVertices.back();

    for (unsigned i = 0; i < m_numNodes; ++i) {
      // The i node cannot be already in the path
      if (std::find(current_path.m_pathVertices.cbegin(),
                    current_path.m_pathVertices.cend(),
                    i) == current_path.m_pathVertices.cend()) {
        const RealNumber& distance_last_to_i = m_distanceMatrix[last_node][i];
        // Check if adding the i node no vertices is not consistent
        bool consistent = true;
        for (unsigned n = 0; n < LEN_PATH && consistent == true; ++n) {
          const auto& level_node = current_path.m_pathVertices[n];
          if (current_path.m_routeLenPerLevel[n] +
              distance_last_to_i > m_maxDistances[level_node]) {
            consistent = false;
          }
        }
        if (consistent == true) {
          nexts->push_back(i);
        }
      }
    }
  }
}


void AlgorithmSolver::compute_all_distances() {
  assert(m_numNodes > 0);
  assert(m_coordX_vertices.size() == m_numNodes);
  assert(m_coordY_vertices.size() == m_numNodes);

  m_distanceMatrix.resize(m_numNodes);
  m_maxDistances.resize(m_numNodes);
  m_distanceFromNearest.resize(m_numNodes);

  for (unsigned i = 0; i < m_numNodes; ++i) {
    auto& row = m_distanceMatrix[i];
    row.resize(m_numNodes);
    for (unsigned j = 0; j < m_numNodes; ++j) {
      int dx = m_coordX_vertices[i] - m_coordX_vertices[j];
      int dy = m_coordY_vertices[i] - m_coordY_vertices[j];
      RealNumber distance = std::sqrt(std::pow(dx, 2) + std::pow(dy, 2));
      row[j] = distance;
    }

    // compute the max distance
    m_maxDistances[i] = m_alpha * row[SCHOOL_NODE];

    // The distance from the nearest
    // First of all I need to erase the distance with itself
    std::remove_reference<decltype(row)>::type purged_row;
    std::copy_if(row.cbegin(), row.cend(), std::back_inserter(purged_row),
                 [] (const RealNumber& d) {
                   return d != 0;
                 });
    m_distanceFromNearest[i] =
        *(std::min_element(purged_row.cbegin(),
                           purged_row.cend()));
  }
}

template<typename T>
std::vector<T> AlgorithmSolver::parse_vector_dat(std::string vector_data) {
  std::vector<T> rtn;
  boost::trim(vector_data);
  while (vector_data.size()) {
    std::string line = vector_data.substr(0, vector_data.find('\n'));
    vector_data.erase(0, line.size() + 1);
    boost::trim(line);
    if (line.size()) {
      std::string index = line.substr(0, line.find(' '));
      std::string x = line.substr(line.find(' '));
      boost::trim(index);
      boost::trim(x);
      unsigned nindex = std::stoi(index);
      int nx = std::stoi(x);
      assert(nindex == rtn.size());
      rtn.push_back(nx);
    }
  }
  return rtn;
}

void AlgorithmSolver::parse_problem_dat(const std::string& filename) {
  constexpr unsigned BUFFER_SIZE = 1024 * 1024;

  const auto lambda_parse_command = [this](std::string cmd) {
    boost::trim(cmd);
    size_t fin1, fin2;

    fin1 = cmd.find("param ");
    fin2 = cmd.find(":=");
    if (fin1 != std::string::npos) {
      auto param_name = cmd.substr(fin1 + 6, fin2 - fin1 - 6);
      boost::trim(param_name);
      param_name = param_name.substr(0, param_name.find(' '));
      if (param_name == "n") {
        std::string num_vertices = cmd.substr(fin2 + 2);
        boost::trim(num_vertices);
        m_numNodes = std::stoi(num_vertices) + 1;
      } else if (param_name == "alpha") {
        std::string alpha = cmd.substr(fin2 + 2);
        boost::trim(alpha);
        m_alpha = stod(alpha);
      } else if (param_name == "coordX") {
        m_coordX_vertices.clear();
        std::string vector = cmd.substr(fin2 + 2);
        m_coordX_vertices =
            parse_vector_dat<decltype(m_coordX_vertices)::value_type>(
                std::move(vector));
      } else if (param_name == "coordY") {
        m_coordY_vertices.clear();
        std::string vector = cmd.substr(fin2 + 2);
        m_coordY_vertices =
            parse_vector_dat<decltype(m_coordY_vertices)::value_type>(
                std::move(vector));
      }
    }
  };

  std::ifstream file;
  file.open(filename);
  if (file.fail()) {
    throw std::runtime_error(std::string(
        "Cannot open the file '" + filename + "'"));
  }

  char buffer[BUFFER_SIZE];
  std::string command;
  unsigned read_chars;

  while (file.eof() == false) {
    command.clear();
    do {
      file.get(buffer, BUFFER_SIZE, ';');
      read_chars = file.gcount();
      command.append(buffer, read_chars);
      file.get();
    } while (read_chars == BUFFER_SIZE);

    lambda_parse_command(command);
  }

  file.close();
}

int main(int argc, char *argv[]) {
  AlgorithmSolver algorithm;
  algorithm.parse_problem_dat("pedibus_300.dat");
  algorithm.compute_all_distances();
  std::ofstream file;
  file.open("coords20.csv");
  algorithm.write_on_csv_coordinate(&file);
  file.close();

  algorithm.launchProblemSolver<AlgorithmSolver::Heuristic>();

  return 0;
}
