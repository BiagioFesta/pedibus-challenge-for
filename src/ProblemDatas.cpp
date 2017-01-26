// Copyright 2017 <Biagio Festa>
#include <utility>
#include <string>
#include <cmath>
#include <fstream>
#include <algorithm>
#include <set>
#include "ProblemDatas.hpp"

void ProblemDatas::parse_problem_dat(const std::string& filename) {
  constexpr unsigned BUFFER_SIZE = 1024 * 1024;

  m_numNodes = 0;
  m_numEdges = 0;
  m_alpha = 0;
  m_coordX_vertices.clear();
  m_coordY_vertices.clear();
  m_distance_matrix.clear();
  m_max_distances.clear();
  m_distances_FromNearest.clear();

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
        m_numEdges = std::pow(m_numNodes - 1, 2);
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

  compute_all_distances();
  compute_edges_indices();
}

void ProblemDatas::compute_all_distances() {
  assert(m_numNodes > 0);
  assert(m_coordX_vertices.size() == m_numNodes);
  assert(m_coordY_vertices.size() == m_numNodes);

  m_distance_matrix.resize(m_numNodes);
  m_max_distances.resize(m_numNodes);
  m_distances_FromNearest.resize(m_numNodes);

  for (unsigned i = 0; i < m_numNodes; ++i) {
    auto& row = m_distance_matrix[i];
    row.resize(m_numNodes);
    for (unsigned j = 0; j < m_numNodes; ++j) {
      int dx = m_coordX_vertices[i] - m_coordX_vertices[j];
      int dy = m_coordY_vertices[i] - m_coordY_vertices[j];
      RealNumber distance = std::sqrt(std::pow(dx, 2) + std::pow(dy, 2));
      row[j] = distance;
    }

    // compute the max distance
    m_max_distances[i] = m_alpha * row[SCHOOL_INDEX];

    // The distance from the nearest
    // First of all I need to erase the distance with itself
    std::remove_reference<decltype(row)>::type purged_row;
    std::copy_if(row.cbegin(), row.cend(), std::back_inserter(purged_row),
                 [] (const RealNumber& d) {
                   return d != 0;
                 });
    m_distances_FromNearest[i] =
        *(std::min_element(purged_row.cbegin(),
                           purged_row.cend()));
  }
}

void ProblemDatas::compute_edges_indices() {
  m_mapEdge_index2link.clear();
  m_mapEdge_link2index.clear();
  m_map_vertex2outedges.clear();
  m_map_vertex2inedges.clear();

  EdgeIndex e = 0;
  for (VertexIndex i = 1; i < m_numNodes; ++i) {
    for (VertexIndex j = 0; j < m_numNodes; ++j) {
      if (j != i) {
        EdgeLink link = {i, j};
        m_mapEdge_index2link[e] = link;
        m_mapEdge_link2index[link] = e;
        m_map_vertex2outedges[i].insert(e);
        m_map_vertex2inedges[j].insert(e);
        ++e;
      }
    }
  }

  m_map_vertex2outedges[SCHOOL_INDEX] = std::set<EdgeIndex>();

  assert(m_mapEdge_index2link.size() == m_numEdges);
  assert(m_mapEdge_link2index.size() == m_numEdges);
  assert(m_map_vertex2outedges.size() == m_numNodes);
  assert(m_map_vertex2inedges.size() == m_numNodes);
}
