/**
 * @file subgraph.cpp
 * Test script for viewing a subgraph from our Graph
 *
 * @brief Reads in two files specified on the command line.
 * First file: 3D Points (one per line) defined by three doubles
 * Second file: Tetrahedra (one per line) defined by 4 indices into the point
 * list
 */

#include <fstream>
#include <iterator>

#include "CME212/SDLViewer.hpp"
#include "CME212/Util.hpp"

#include "Graph.hpp"

/** An iterator that skips over elements of another iterator based on whether
 * those elements satisfy a predicate.
 *
 * Given an iterator range [@a first, @a last) and a predicate @a pred,
 * this iterator models a filtered range such that all i with
 * @a first <= i < @a last and @a pred(*i) appear in order of the original range.
 */
template <typename Pred, typename It>
class filter_iterator
    : private equality_comparable<filter_iterator<Pred,It>> {
 public:
  // Get all of the iterator traits and make them our own
  typedef typename std::iterator_traits<It>::value_type        value_type;
  typedef typename std::iterator_traits<It>::pointer           pointer;
  typedef typename std::iterator_traits<It>::reference         reference;
  typedef typename std::iterator_traits<It>::difference_type   difference_type;
  typedef typename std::input_iterator_tag                     iterator_category;

  typedef filter_iterator<Pred,It> self_type;

  // Constructor
  filter_iterator(const Pred& p, const It& first, const It& last)
      : p_(p), it_(first), end_(last) {
    // HW1 #4: YOUR CODE HERE
  }

  value_type operator*() const {
    return *it_;
  }

  self_type& operator++(){
    ++it_;
    if(!p_(it_) and it_ != end_) ++it_;
    return *this;
  }

  bool operator==(const self_type& fi) const {
    return it_ == fi.it_ and end_ == fi.end_;
  }

 private:
  Pred p_;
  It it_;
  It end_;
};

/** Helper function for constructing filter_iterators.
 *
 * Usage:
 * // Construct an iterator that filters odd values out and keeps even values.
 * std::vector<int> a = ...;
 * auto it = make_filtered(a.begin(), a.end(), [](int k) {return k % 2 == 0;});
 */
template <typename Pred, typename Iter>
filter_iterator<Pred,Iter> make_filtered(const Iter& it, const Iter& end,
                                         const Pred& p) {
  return filter_iterator<Pred,Iter>(p, it, end);
}

/** One Predicate Functor.
 *  Filter all the nodes with odd indices.
 */
struct OddPredicate {
  template <typename NODE>
  bool operator()(const NODE& ni) {
    return ni.index() % 2 == 0;
  }
};

/** One Predicate Functor.
 *  Filter all the nodes with z-coordinate < 0.
 */
struct BelowPredicate {
  template <typename NODE>
  bool operator()(const NODE& ni) {
    return ni.position().z >= 0;
  }
};

/** One Predicate Functor
 *  Filter all the nodes with probability 0.9, only 10% of nodes are left.
 */
struct ProbPredicate {
  template <typename NODE>
  bool operator()(const NODE& ni) {
    (void) ni;
    return rand() % 1 > 0.9;
  }
};

/** Test predicate for HW1 #4 */
struct SlicePredicate {
  template <typename NODE>
  bool operator()(const NODE& n) {
    return n.position().x < 0;
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
  typedef Graph<int, double> GraphType;
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

  // Test the OddPredicate.
  //OddPredicate op;
  //filter_iterator<OddPredicate, GraphType::node_iterator> filter_begin = make_filtered(graph.node_begin(), graph.node_end(), op);
  //filter_iterator<OddPredicate, GraphType::node_iterator> filter_end = make_filtered(graph.node_end(), graph.node_end(), op);

  // Test the BelowPredicate.
  // BelowPredicate bp;
  // filter_iterator<BelowPredicate, GraphType::node_iterator> filter_begin = make_filtered(graph.node_begin(), graph.node_end(), bp);
  // filter_iterator<BelowPredicate, GraphType::node_iterator> filter_end = make_filtered(graph.node_end(), graph.node_end(), bp);

  // Test the ProbPredicate.
  ProbPredicate pp;
  filter_iterator<ProbPredicate, GraphType::node_iterator> filter_begin = make_filtered(graph.node_begin(), graph.node_end(), pp);
  filter_iterator<ProbPredicate, GraphType::node_iterator> filter_end = make_filtered(graph.node_end(), graph.node_end(), pp);

  auto node_map = viewer.empty_node_map(graph);
  viewer.add_nodes(filter_begin, filter_end, node_map);
  viewer.add_edges(graph.edge_begin(), graph.edge_end(), node_map);
  viewer.center_view();

  return 0;
}
