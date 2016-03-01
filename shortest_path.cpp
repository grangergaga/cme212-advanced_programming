/**
 * @file shortest_path.cpp
 * Test script for using our templated Graph to determine shortest paths.
 *
 * @brief Reads in two files specified on the command line.
 * First file: 3D Points (one per line) defined by three doubles
 * Second file: Tetrahedra (one per line) defined by 4 indices into the point
 * list
 */

#include <vector>
#include <queue>
#include <fstream>

#include "CME212/SDLViewer.hpp"
#include "CME212/Util.hpp"
#include "CME212/Color.hpp"

#include "Graph.hpp"


/** Comparator that compares the distance from a given point p.
 */
struct MyComparator {
   Point p_;
   MyComparator(const Point& p) : p_(p) {
   };

   template <typename NODE>
   bool operator()(const NODE& node1, const NODE& node2) const {
     Point diff1 = node1.position() - p_;
     Point diff2 = node2.position() - p_;
     if (norm(diff1) < norm(diff2)) return true;
    return false;
  }
};


/** Calculate shortest path lengths in @a g from the nearest node to @a point.
 * @param[in,out] g Input graph
 * @param[in] point Point to find the nearest node to.
 * @post Graph has modified node values indicating the minimum path length
 *           to the nearest node to @a point
 * @post Graph nodes that are unreachable to the nearest node to @a point have
 *           the value() -1.
 * @return The maximum path length found.
 *
 * Finds the nearest node to @a point and treats that as the root node for a
 * breadth first search.
 * This sets node's value() to the length of the shortest path to
 * the root node. The root's value() is 0. Nodes unreachable from
 * the root have value() -1.
 */
int shortest_path_lengths(Graph<int, int>& g, const Point& point) {
  MyComparator mc = MyComparator(point);

  // Find the closest node to the given Point as root.
  Graph<int, int>::node_iterator nifirst = g.node_begin();
  Graph<int, int>::node_iterator nilast = g.node_end();
  Graph<int, int>::node_iterator niroot = std::min_element(nifirst, nilast, mc);
  Graph<int, int>::node_type root = *niroot;

  // Set all the nodes' default values to -1 and the root's value to 0.
  for(; nifirst != nilast; ++nifirst){
    (*nifirst).value() = -1;
  }
  root.value() = 0;

  // Set the current longest distance from root.
  int max = 0;

  /** The queue is to store the nodes needed to be evaluated.
   *  For each node n in the queue, we iterate the adjent nodes and set their values,
   *  push them into the queue and then pop n. In this process, update max.
   */
  std::queue<Graph<int, int>::node_type> waiting;
  waiting.push(root);
  while(!waiting.empty()){
    Graph<int, int>::node_type r = waiting.front();
    int cur = r.value();
    Graph<int, int>::incident_iterator rbegin = r.edge_begin();
    Graph<int, int>::incident_iterator rend = r.edge_end();
    for(; rbegin != rend; ++rbegin){
      if((*rbegin).node2().value() == -1){
        (*rbegin).node2().value() = cur + 1;
        if((*rbegin).node2().value() > max) max = (*rbegin).node2().value();
        waiting.push((*rbegin).node2());
      }
      if((*rbegin).node2().value() > cur + 1){
        (*rbegin).node2().value() = cur + 1;
        if((*rbegin).node2().value() > max) max = (*rbegin).node2().value();
      }
    }
    waiting.pop();
  }
  return max;
}

/** One ColorFunction described in the homework document, it defines the color
 *  based on the node's value comparing to the reference value maxpath_, which
 *  is returned by the shortest_path_lengths function.
 */
struct PathColorFn{
  int maxpath_;
  PathColorFn(const int p) : maxpath_(p) {
  };
  CME212::Color operator()(Graph<int, int>::node_type n){
    float c = float(n.value()) / float(maxpath_);
    return CME212::Color::make_heat(1.0 - c);
  }
};


/** Another ColorFunction based on the node's position. It shows red color
 *  near the middle of the graph and purple far from the middle, in terms of
 *  x-coordinate.
 */
struct PositionColorFn{
  PositionColorFn(){
  }
  CME212::Color operator()(Graph<int, int>::node_type n){
    double nm = norm(n.position());
    double x = n.position().x;
    double c = std::abs(x) / nm;
    return CME212::Color::make_heat(c);
  }
};


int main(int argc, char** argv)
{
  // Check arguments
  if (argc < 3) {
    std::cerr << "Usage: " << argv[0] << " NODES_FILE TETS_FILE\n";
    exit(1);
  }

  // Construct a Graph
  typedef Graph<int, int> GraphType;
  GraphType graph;
  std::vector<GraphType::node_type> nodes;

  // Create a nodes_file from the first input argument
  std::ifstream nodes_file(argv[1]);
  // Interpret each line of the nodes_file as a 3D Point and add to the Graph
  Point p;
  while (CME212::getline_parsed(nodes_file, p))
    nodes.push_back(graph.add_node(p));

  // Create a tets_file from the second input argument
  std::ifstream tets_file(argv[2]);
  // Interpret each line of the tets_file as four ints which refer to nodes
  std::array<int,4> t;
  while (CME212::getline_parsed(tets_file, t))
    for (unsigned i = 1; i < t.size(); ++i)
      for (unsigned j = 0; j < i; ++j)
        graph.add_edge(nodes[t[i]], nodes[t[j]]);

  // Print out the stats
  std::cout << graph.num_nodes() << " " << graph.num_edges() << std::endl;

  // Launch the SDLViewer
  CME212::SDLViewer viewer;
  viewer.launch();
  auto node_map = viewer.empty_node_map(graph);

  Point pref = Point(-1, 0, 1);
  int path = shortest_path_lengths(graph, pref);
  PathColorFn pcf = PathColorFn(path);
  viewer.add_nodes(graph.node_begin(), graph.node_end(), pcf, node_map);

  // Test the PositionColorFn, the color is presented according to the nodes' x coordinats.
  //PositionColorFn pocf = PositionColorFn();
  //viewer.add_nodes(graph.node_begin(), graph.node_end(), pocf, node_map);

  viewer.add_edges(graph.edge_begin(), graph.edge_end(), node_map);
  viewer.center_view();

  return 0;
}
