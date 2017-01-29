// Copyright 2017 <Biagio Festa>
#ifndef __FOR_CH__HE_SOLVER__HPP
#define __FOR_CH__HE_SOLVER__HPP
#include <memory>
#include <stack>
#include <algorithm>
#include <iterator>
#include "ProblemDatas.hpp"

namespace for_ch {

class HESolver {
 public:
  explicit HESolver(const std::shared_ptr<ProblemDatas>& problem) noexcept;

  bool run(std::vector<bool>* active_edges);
 private:
  struct Path {
    std::vector<VertexIndex> m_pathVertices;
    std::vector<RealNumber> m_routeLenPerLevel;
  };

  struct Heuristic {
    explicit Heuristic(const HESolver& solver) :
        mp_solver(&solver) {
    }

    RealNumber operator()(const VertexIndex& v1,
                          const VertexIndex& v2) const {
      using PairIndices = EdgeLink;
      // This return an heuristic value which says how much
      // v2 is good respect to v1

      const auto& distanceMatrix =
          mp_solver->mp_problem->m_distance_matrix;
      const auto& nearestDistance =
          mp_solver->mp_problem->m_distances_FromNearest;

      // Check whether v2 is in a existent path
      bool v2_in_a_path_as_inner = false;
      bool v2_is_a_leaf = false;
      PairIndices indices;
      if (mp_solver->checkVertexIsInAPath(v2, &indices)) {
        if (indices.second != 0) {
          v2_in_a_path_as_inner = true;
        } else {
          v2_is_a_leaf = true;
        }
      }

      if (v2_is_a_leaf == true) {
        return 0;
      }

      return ((m_A * nearestDistance[v2]) +
              (m_B * distanceMatrix[v1][v2]) +
              (m_C * v2_in_a_path_as_inner) +
              (m_D * distanceMatrix[v2][SCHOOL_INDEX]));
    }

    // m_A factor scale is the grade of isolation of a point
    RealNumber m_A = 0.1;

    // m_B factor scale of distance from the evaluator point
    RealNumber m_B = 1.0;

    // m_C factor scale if v2 is already in a path (we don't want to link
    // with a internal node in a existent path BUT is not a leaf
    RealNumber m_C = 100;

    // m_D factor scale if v2 is near to school
    RealNumber m_D = 0.4;

    const HESolver* mp_solver;
  };

  template<typename HeuristicT>
  struct CompareVertices {
    explicit CompareVertices(const VertexIndex& v, const HeuristicT& h):
        mp_v(&v), mp_h(&h) {
    }

    bool operator()(const VertexIndex& v1, const VertexIndex& v2) const {
      return (*mp_h)(*mp_v, v1) > (*mp_h)(*mp_v, v2);
    }

    const VertexIndex* mp_v;
    const HeuristicT* mp_h;
  };
  
  std::shared_ptr<ProblemDatas> mp_problem;
  std::vector<Path> m_foundPaths;
  std::set<VertexIndex> m_linkedVertices;

  /// This method takes a path and returns all possibile vertices
  /// you can add to the input path. The vertices you can apply to the path
  /// satisfy all distance constraints
  /// @param [in] current_path   The path you want to know which vertices you
  ///                            can add.
  /// @param [out] nexts         A vector of vertices which contains all
  ///                            vertices you can properly add to that path.
  void find_possible_next_nodes(const Path& current_path,
                                std::vector<VertexIndex>* nexts);
  
  /// This method takes a vertex and adds it to a path
  /// @param [in] v                 The vertex you want to add to the path
  /// @param [in.out] current_path  The path you want to modify
  void apply_next_node_to_path(const VertexIndex& v, Path* path);

  /// Check whether a vertex is already in a path
  ///
  /// @param [out]    ....
  /// @return true whether is already in a path, false otherwise
  ///
  bool checkVertexIsInAPath(const VertexIndex& v,
                            EdgeLink* outIndices) const;

  template<typename HeuristicT>
  bool launchProblemSolver();

  template<typename HeuristicT>
  bool construct_path();

  template<typename HeuristicT>
  void sort_heuristic_accordance(const VertexIndex& lastInPath,
                                 std::vector<VertexIndex>* vertices);
};


template<typename HeuristicT>
bool HESolver::launchProblemSolver() {
  const unsigned num_nodes = mp_problem->m_numNodes;
  m_foundPaths.clear();
  m_linkedVertices.clear();

  while (m_linkedVertices.size() != num_nodes) {
    if (construct_path<HeuristicT>() == false) {
      // if the algorithm cannot construct a path
      // that means there is a problem
      // Not feasible solution found
      return false;
    }
  }

  return true;
}

template<typename HeuristicT>
bool HESolver::construct_path() {
  using OpenList = std::stack<Path>;

  const unsigned num_nodes = mp_problem->m_numNodes;
  
  // Get the list of free vertices (not already linked)
  std::vector<VertexIndex> freeVertices;
  std::vector<RealNumber> distancesFree;
  for (unsigned i = 0; i < num_nodes; ++i) {
    if (m_linkedVertices.find(i) == m_linkedVertices.cend()) {
      freeVertices.push_back(i);
      distancesFree.push_back(
          mp_problem->m_distance_matrix[i][SCHOOL_INDEX]);
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
  EdgeLink indices;
  std::vector<VertexIndex> nexts;
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
    if (current.m_pathVertices.back() == SCHOOL_INDEX) {
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
    sort_heuristic_accordance<HeuristicT>(
        current.m_pathVertices.back(), &nexts);

    // Test all possible choices
    for (const auto& n : nexts) {
      // If n is already in a path the current path must to be linked
      if (checkVertexIsInAPath(n, &indices)) {
        if (indices.second == 0) {
          // In this case the vertex n is a leaf of an existent path
          Path& oldPath = m_foundPaths[indices.first];
          const RealNumber& lenPath = oldPath.m_routeLenPerLevel[0];
          const VertexIndex& leafPath = oldPath.m_pathVertices[0];
          const VertexIndex& lastCurrent = current.m_pathVertices.back();
          const RealNumber& lenLink =
              mp_problem->m_distance_matrix[lastCurrent][leafPath];

          bool admissible = true;
          for (unsigned i = 0;
               i < current.m_pathVertices.size() && admissible == true;
               ++i) {
            RealNumber newDistance = current.m_routeLenPerLevel[i] + lenLink +
                lenPath;
            if (newDistance >
                mp_problem->m_max_distances[current.m_pathVertices[i]]) {
              admissible = false;
            }
          }

          if (admissible == true) {
            std::for_each(   // for each vertices in the current path
                current.m_pathVertices.crbegin(),
                current.m_pathVertices.crend(),
                [this, &oldPath, &lenPath]
                (const VertexIndex& v) {
                  RealNumber newDistance =
                      mp_problem->m_distance_matrix[v][oldPath.m_pathVertices.front()]
                      + lenPath;

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
          const VertexIndex& lastCurrent = current.m_pathVertices.back();
          const Path & oldPath = m_foundPaths[indices.first];
          const RealNumber& lenPath =
              oldPath.m_routeLenPerLevel[indices.second];
          const VertexIndex& linkRef = oldPath.m_pathVertices[indices.second];
          const RealNumber& lenLink =
              mp_problem->m_distance_matrix[lastCurrent][linkRef];

          // Create new nodes as the list from linkRef to end (of old Path)
          std::vector<VertexIndex> newPathNodes;
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
            if (newDistance >
                mp_problem->m_max_distances[current.m_pathVertices[i]]) {
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
void HESolver::sort_heuristic_accordance(
    const VertexIndex& lastInPath,
    std::vector<VertexIndex>* vertices) {
  assert(vertices != nullptr);
  HeuristicT h(*this);
  CompareVertices<HeuristicT> comp(lastInPath, h);
  std::sort(vertices->begin(), vertices->end(), comp);
}



}  // namespace for_ch

#endif  // __FOR_CH__HE_SOLVER__HPP
