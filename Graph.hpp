#ifndef CME212_GRAPH_HPP
#define CME212_GRAPH_HPP

/** @file Graph.hpp
 * @brief An undirected graph type
 */

#include <algorithm>
#include <vector>
#include <map>
#include <cassert>

#include "CME212/Util.hpp"
#include "CME212/Point.hpp"


/** @class Graph
 * @brief A template for 3D undirected graphs.
 *
 * Users can add and retrieve nodes and edges. Edges are unique (there is at
 * most one edge between any pair of distinct nodes).
 */
template <typename V, typename E> class Graph {
 public:

  //
  // PUBLIC TYPE DEFINITIONS
  //

  /** Type of the node value. */
  typedef V node_value_type;

  /** Type of the edge value */
  typedef E edge_value_type;

  /** Type of this graph. */
  typedef Graph graph_type;

  /** Predeclaration of Node type. */
  class Node;
  /** Synonym for Node (following STL conventions). */
  typedef Node node_type;

  /** Predeclaration of Edge type. */
  class Edge;
  /** Synonym for Edge (following STL conventions). */
  typedef Edge edge_type;

  /** Type of node iterators, which iterate over all graph nodes. */
  class NodeIterator;
  /** Synonym for NodeIterator */
  typedef NodeIterator node_iterator;

  /** Type of edge iterators, which iterate over all graph edges. */
  class EdgeIterator;
  /** Synonym for EdgeIterator */
  typedef EdgeIterator edge_iterator;

  /** Type of incident iterators, which iterate incident edges to a node. */
  class IncidentIterator;
  /** Synonym for IncidentIterator */
  typedef IncidentIterator incident_iterator;

  /** Type of indexes and sizes.
      Return type of Graph::Node::index(), Graph::num_nodes(),
      Graph::num_edges(), and argument type of Graph::node(size_type) */
  typedef unsigned size_type;

  //
  // CONSTRUCTORS AND DESTRUCTOR
  //

  /** Construct an empty graph. */
  Graph()
  // Initialize the private attributes of a graph,with no nodes and no edges.
    : nodes_(std::vector<std::pair<Point*, node_value_type>>(0)),
    adjacency_(std::vector<std::vector<std::pair<size_type, edge_value_type>>>(0)){

  }

  /** Default destructor */
  ~Graph() = default;

  //
  // NODES
  //

  /** @class Graph::Node
   * @brief Class representing the graph's nodes.
   *
   * Node objects are used to access information about the Graph's nodes.
   */
  class Node: private totally_ordered <Node> {

   public:
    /** Construct an invalid node.
     *
     * Valid nodes are obtained from the Graph class, but it
     * is occasionally useful to declare an @i invalid node, and assign a
     * valid node to it later.
     */
    Node(){
    }

    /** Return this node's position, whcih can be modified.*/
    Point& position(){
      assert(this->graph_ != NULL);
      return *((*graph_).nodes_[node_id_].first);
    }

    /** Return this node's position. */
    const Point& position() const{
      assert(this->graph_ != NULL);
      return *((*graph_).nodes_[node_id_].first);
    }

    /** Return this node's index, a number in the range [0, graph_size_). */
    size_type index() const {
      assert(this->graph_ != NULL);
      /** Check if the current node is a valid node
       * (otherwise it doesn't make sense to get the index)
       * And also check if the index is in the range of [0, graph_size_). */
      assert(graph_->size() != 0 and graph_->size() > node_id_);
      return node_id_;
    }

    /** Test whether this node and @a n are equal.
     *
     * Equal nodes have the same graph and the same index.
     */
    bool operator==(const Node& n) const {
      if(graph_ == n.graph_ and node_id_ == n.node_id_){
        return true;
      }
      return false;
    }

    /** Test whether this node is less than @a n in a global order.
     *
     * This ordering function is useful for STL containers such as
     * std::map<>. It need not have any geometric meaning.
     *
     * The node ordering relation must obey trichotomy: For any two nodes x
     * and y, exactly one of x == y, x < y, and y < x is true.
     *
     * It may not make sense to compare the pointer of graph, but we need a
     * way to make the ordering relation obey trichotomy.
     */
    bool operator<(const Node& n) const {
      assert(n.graph_ != NULL);
      if(graph_ == n.graph_ and node_id_ < n.node_id_) return true;
      else if(graph_ < n.graph_) return true;
      return false;
    }

    /** Non-const version of this node's value,
     *  get access and set value to it.
     */
    node_value_type& value(){
      assert(this->graph_ != NULL);
      return (*graph_).nodes_[node_id_].second;
    }

    /** Non-const version of this node's value,
     *  only get access to it.
     */
    const node_value_type& value() const {
      assert(this->graph_ != NULL);
      return (*graph_).nodes_[node_id_].second;
    }

    /** Take a look at the adjacency table,
     *  return the number of nodes adjacent to the current node.
     */
    size_type degree() const{
      return graph_->adjacency_[node_id_].size();
    }

    /** Return the incident-iterator related to the current node.
     *  find the begin() map::iterator in the adjacent table.
     *  Note that current node's id is stored as node1 in the related
     *  edge of the incident iterator.
     */
    incident_iterator edge_begin() const{
      incident_iterator ii = IncidentIterator(graph_, node_id_, 0);
      return ii;
    }

    /** Return the incident-iterator related to the current node.
     *  find the end() map::iterator in the adjacent table.
     *  Note that current node's id is stored as node1 in the related
     *  edge of the incident iterator.
     */
    incident_iterator edge_end() const {
      size_type incnum = graph_->adjacency_[node_id_].size();
      incident_iterator ii = IncidentIterator(graph_, node_id_, incnum);
      return ii;
    }

   private:
     /** graph_ is a pointer to the Graph object where the node is located.
      * node_id_ is the index of this node in the graph. */
     Graph* graph_;
     size_type node_id_;

     /** Construct a valid Node using graph and index */
     Node(const Graph* graph, size_type node_id)
         : graph_(const_cast<Graph*>(graph)), node_id_(node_id) {
     }

    // Allow Graph to access Node's private member data and functions.
    friend class Graph;
  };

  /** Return the number of nodes in the graph.
   *
   * Complexity: O(1).
   */
  size_type size() const {
    return nodes_.size();
  }

  /** Synonym for size(). */
  size_type num_nodes() const {
    return size();
  }

  /** Add a node to the graph, returning the added node.
   * @param[in] position The new node's position
   * @post new num_nodes() == old num_nodes() + 1
   * @post result_node.index() == old num_nodes()
   * @post result_node.value() is default value of the node_value_type
   *
   * Complexity: O(1) amortized operations.
   */
  Node add_node(const Point& position) {
    Point *curpoint = new Point;
    *curpoint = position;
    node_value_type innervalue = node_value_type();
    std::pair<Point*, node_value_type> curnode(curpoint, innervalue);
    nodes_.push_back(curnode);
    /** Update the attributes of graph. */
    std::vector<std::pair<size_type, edge_value_type>> curadj;
    adjacency_.push_back(curadj);
    return Node(this, size() - 1);
    return Node();
  }

  /** Add a node to the graph, returning the added node.
   * @param[in] position The new node's position
   * @param[in] innervalue The new node's value
   * @post new num_nodes() == old num_nodes() + 1
   * @post result_node.index() == old num_nodes()
   * @post result_node.value() == innervalue
   *
   * Complexity: O(1) amortized operations.
   */
  Node add_node(const Point& position, const node_value_type& innervalue) {
    Point *curpoint = new Point;
    *curpoint = position;
    std::pair<Point*, node_value_type> curnode(curpoint, innervalue);
    nodes_.push_back(curnode);
    /** Update the attributes of graph. */
    std::vector<std::pair<size_type, edge_value_type>> curadj;
    adjacency_.push_back(curadj);
    return Node(this, size() - 1);
    return Node();
  }

  /** Determine if a Node belongs to this Graph
   * @return True if @a n is currently a Node of this Graph
   *
   * Complexity: O(1).
   */
  bool has_node(const Node& n) const {
    assert(n.graph_ != NULL);
    if(n.graph_ == this and n.node_id_ < size()){
      return true;
    }
    return false;
  }

  /** Return the node with index @a i.
   * @pre 0 <= @a i < num_nodes()
   * @post result_node.index() == i
   *
   * Complexity: O(1).
   */
  Node node(size_type i) const {
    assert(i < size());
    return Node(this, i);
    return Node();
  }

  /** Remove node n from current graph.
    * @post new graph.size() == old graph.size() - 1;
    * @post new num_edges() == old num_edges() - n.degree();
    * complexity: O(degree()^2)
    * Invalidate: nodes with node_id_ == n.node_id_ or node_id_ == old graph.size() - 1;
                  edges whose at least one endpoint incident to n;
                  node_iterators representing the invalid nodes;
                  edge_iterators representing the invalid edges;

  */
  void remove_node(const Node& n){
    assert(has_node(n));
    auto g = n.graph_;
    auto nid = n.node_id_;
    auto lid = g->size() - 1;
    g->nodes_[nid] = g->nodes_[lid];
    g->nodes_.pop_back();
    for(size_type i = 0; i < g->adjacency_[nid].size(); ++i){
      size_type oid = g->adjacency_[nid][i].first;
      for(size_type j = 0; j < g->adjacency_[oid].size(); ++j){
        if(g->adjacency_[oid][j].first == nid){
          g->adjacency_[oid][j] = g->adjacency_[oid].back();
          g->adjacency_[oid].pop_back();
        }
      }
    }
    g->adjacency_[nid] = g->adjacency_[lid];
    g->adjacency_.pop_back();
    for(size_type i = 0; i < g->adjacency_[nid].size(); ++i){
      size_type loid = g->adjacency_[nid][i].first;
      for(size_type j = 0; j < g->adjacency_[loid].size(); ++j){
        if(g->adjacency_[loid][j].first == lid){
          g->adjacency_[loid][j].first = nid;
        }
      }
    }

  }

  node_iterator remove_node(node_iterator n_it){
    assert(n_it != node_end());
    auto n = *n_it;
    remove_node(n);
    return n_it;
  }

  //
  // EDGES
  //

  /** @class Graph::Edge
   * @brief Class representing the graph's edges.
   *
   * Edges are order-insensitive pairs of nodes. Two Edges with the same nodes
   * are considered equal if they connect the same nodes, in either order.
   */
  class Edge: private totally_ordered <Edge> {
   public:
    /** Construct an invalid Edge. */
    Edge(){
    }

    size_type Node2Id() const{
      assert(this->graph_ != NULL);
      return graph_->adjacency_[node1_id_][node2_vecid_].first;
    }

    /** Return a node of this Edge */
    Node node1() const {
      assert(this->graph_ != NULL);
      return Node(graph_, node1_id_);
      return Node();
    }

    /** Return the other node of this Edge */
    Node node2() const {
      assert(this->graph_ != NULL);
      return Node(graph_, Node2Id());
      return Node();
    }

    /** Test whether this edge and @a e are equal.
     *
     * Equal edges represent the same undirected edge between two nodes.
     */
    bool operator==(const Edge& e) const {
      assert(e.graph_ != NULL and this->graph_ != NULL);
      if(graph_ == e.graph_ and
        ((node1_id_ == e.node1_id_ and Node2Id() == e.Node2Id()) or
        (Node2Id() == e.node1_id_ and node1_id_ == e.Node2Id()))) return true;
      return false;
    }

    /** Test whether this edge is less than @a e in a global order.
     * The order is according to the edge index.
     * This ordering function is useful for STL containers such as
     * std::map<>. It need not have any interpretive meaning.
     */
    bool operator<(const Edge& e) const {
      assert(e.graph_ != NULL and this->graph_ != NULL);
      size_type curmin = std::min(node1_id_, Node2Id());
      size_type curmax = std::max(node1_id_, Node2Id());
      size_type emin = std::min(e.node1_id_, e.Node2Id());
      size_type emax = std::max(e.node1_id_, e.Node2Id());
      if ((graph_ < e.graph_) or
  			((graph_ == e.graph_) and ((curmin < emin) or
  		   ((curmin == emin) and (curmax < emax))))) {
  		  return true;
  		}
      return false;
    }

    /** Non-const version of this edge's value,
     *  get access and set value to it.
     */
    edge_value_type& value(){
      assert(this->graph_ != NULL);
      return (*graph_).adjacency_[node1_id_][node2_vecid_].second;
    }

    /** Const version of this edge's value,
     *  only get access to it.
     */
     const edge_value_type& value() const{
       assert(this->graph_ != NULL);
       return (*graph_).adjacency_[node1_id_][node2_vecid_].second;
     }

    double length() const{
      return norm(node1().position() - node2().position());
    }

    Edge dual() const{
      size_type node2id = graph_->adjacency_[node1_id_][node2_vecid_].first;
      size_type p = 0;
      for(; p < graph_->adjacency_[node2id].size(); ++p){
        if(graph_->adjacency_[node2id][p].first == node1_id_){
          return Edge(graph_, node2id, p);
        }
      }
      return Edge();
    }

   private:
     /** graph_ is a pointer to the Graph object where the edge is located.
      * edge_id_ is the index of this edge in the graph. */
     Graph* graph_;
     size_type node1_id_;
     size_type node2_vecid_;

     /** Construct a valid edge */
     Edge(const Graph* graph, size_type node1_id, size_type node2_vecid)
         : graph_(const_cast<Graph*>(graph)), node1_id_(node1_id), node2_vecid_(node2_vecid){
     }

    // Allow Graph to access Edge's private member data and functions.
    friend class Graph;
  };

  /** Return the total number of edges in the graph.
   *
   * Complexity: O(1)
   */
  size_type num_edges() const {
    size_type numedges = 0;
    for(auto &i : adjacency_){
      numedges += i.size();
    }
    return numedges/2;
  }

  /** Return the edge with index @a i.
   * @pre 0 <= @a i < num_edges()
   *
   * Complexity: O(1)
   */
   Edge edge(size_type i) const {
     assert(i < num_edges());
     return *std::next(edge_begin(), i);
   }

  /** Test whether two nodes are connected by an edge.
   * @pre @a a and @a b are valid nodes of this graph
   * @return True if for some @a i, edge(@a i) connects @a a and @a b.
   *
   * Complexity: O(1)
   */
  bool has_edge(const Node& a, const Node& b) const {
    assert(this == a.graph_ and this == b.graph_);
    size_type aid = a.node_id_;
    size_type bid = b.node_id_;
    assert(aid < size() and bid < size());
    size_type deg = a.degree();

    for(size_type i = 0; i < deg; ++i){
      if(adjacency_[aid][i].first == bid) return true;
    }
    return false;
  }

  /** Add an edge to the graph, or return the current edge if it already exists.
   * @pre @a a and @a b are distinct valid nodes of this graph
   * @return an Edge object e with e.node1() == @a a and e.node2() == @a b
   * @post has_edge(@a a, @a b) == true
   * @post If old has_edge(@a a, @a b), new num_edges() == old num_edges().
   *       Else,                        new num_edges() == old num_edges() + 1.
   *
   * Can invalidate edge indexes -- in other words, old edge(@a i) might not
   * equal new edge(@a i). Must not invalidate outstanding Edge objects.
   *
   * Complexity: O(1)
   */
  Edge add_edge(const Node& a, const Node& b) {
    assert(this == a.graph_ and this == b.graph_);
    size_type aid = a.node_id_;
    size_type bid = b.node_id_;
    assert(aid < size() and bid < size() and aid != bid);

    size_type deg = a.degree();

    for(size_type i = 0; i < deg; ++i){
      if(adjacency_[aid][i].first == bid) return Edge(this, aid, i);
    }

    /** If there isn't such an edge, create one, update the attributes of the graph,
     * then return this newly created edge. */
    edge_value_type ev = edge_value_type();
    std::pair<size_type, edge_value_type> bpair(bid, ev);
    std::pair<size_type, edge_value_type> apair(aid, ev);
    adjacency_[aid].push_back(bpair);
    adjacency_[bid].push_back(apair);
    assert(this->has_edge(a, b) == true and this->has_edge(b, a) == true);
    return Edge(this, aid, adjacency_[aid].size() - 1);
    return Edge();
  }

  /** Remove edges between two nodes.
    * @post new has_edge(a, b) == FALSE;
    * @post if old has_edge(a, b) then new num_edges() == old num_edges() - 1;
    * @complexity: O(a.degree() + b.degree())
    * @invalidates: edge(graph, aid, old adjacency_[aid].size() - 1);
                    edge(graph, bid, old adjacency_[bid].size() - 1);
                    edge_iterators representing the invalid edges above;
    */
  bool remove_edge(const Node& a, const Node& b){
    if(!has_edge(a, b)) return false;
    size_type aid = a.node_id_;
    size_type bid = b.node_id_;
    for(size_type i = 0; i < adjacency_[aid].size(); ++i){
      if(adjacency_[aid][i].first == bid){
        adjacency_[aid][i] = adjacency_[aid][adjacency_[aid].size() - 1];
        adjacency_[aid].pop_back();
      }
    }
    for(size_type i = 0; i < adjacency_[bid].size(); ++i){
      if(adjacency_[bid][i].first == aid){
        adjacency_[bid][i] = adjacency_[bid][adjacency_[bid].size() - 1];
        adjacency_[bid].pop_back();
      }
    }
    return !(has_edge(a, b));
  }

  bool remove_edge(Edge& e){
    auto n1 = e.node1();
    auto n2 = e.node2();
    return remove_edge(n1, n2);
  }

  /** Use edge_iterator to remove edge and return an edge_iterator,
    * which should be fixed after removing.
    */
  edge_iterator remove_edge(edge_iterator e_it){
    auto e = *e_it;
    remove_edge(e);
    e_it.fix();
    return e_it;
  }


  /** Remove all nodes and edges from this graph.
   * @post num_nodes() == 0 && num_edges() == 0
   *
   * Invalidates all outstanding Node and Edge objects.
   */
  void clear() {
    nodes_.clear();
    adjacency_.clear();
  }

  //
  // Node Iterator
  //

  /** @class Graph::NodeIterator
   * @brief Iterator class for nodes. A forward iterator. */
  class NodeIterator: private totally_ordered <NodeIterator> {
   public:
    // These type definitions help us use STL's iterator_traits.
    /** Element type. */
    typedef Node value_type;
    /** Type of pointers to elements. */
    typedef Node* pointer;
    /** Type of references to elements. */
    typedef Node& reference;
    /** Iterator category. */
    typedef std::input_iterator_tag iterator_category;
    /** Difference between iterators */
    typedef std::ptrdiff_t difference_type;

    /** Construct an invalid NodeIterator. */
    NodeIterator() {
    }

    Node operator*() const {
      return Node(graph_, nodeindex_);
    }

    node_iterator& operator++() {
      nodeindex_ ++;
      return *this;
    }

    bool operator==(const NodeIterator& ni) const {
      return graph_ == ni.graph_ and nodeindex_ == ni.nodeindex_;
    }

   private:
     /** Store the pointer to graph and the nodeindex to get access to the nodes
      *  and iterate through.
      */
     Graph* graph_;
     size_type nodeindex_;
     NodeIterator(const Graph* graph, size_type nodeindex)
      : graph_(const_cast<Graph*>(graph)), nodeindex_(nodeindex){}

     friend class Graph;

  };

  // Iterator to first node.
  node_iterator node_begin() const {
    node_iterator ni = NodeIterator(this, 0);
    return ni;
  }

  // Iterator to last node.
  node_iterator node_end() const {
    node_iterator ni = NodeIterator(this, size());
    return ni;
  }

  //
  // Edge Iterator
  //

  /** @class Graph::EdgeIterator
   * @brief Iterator class for edges. A forward iterator. */
  class EdgeIterator: private totally_ordered <EdgeIterator> {
   public:
    // These type definitions help us use STL's iterator_traits.
    /** Element type. */
    typedef Edge value_type;
    /** Type of pointers to elements. */
    typedef Edge* pointer;
    /** Type of references to elements. */
    typedef Edge& reference;
    /** Iterator category. */
    typedef std::input_iterator_tag iterator_category;
    /** Difference between iterators */
    typedef std::ptrdiff_t difference_type;

    /** Construct an invalid EdgeIterator. */
    EdgeIterator() {
    }

    Edge operator*() const {
      return Edge(graph_, center_, outid_);
    }

    EdgeIterator& operator++() {
      ++outid_;
      fix();
      return *this;
    }

    bool operator==(const EdgeIterator& ei) const {
      return (graph_ == ei. graph_ and center_ == ei.center_ and outid_ == outid_);
    }

    void fix() {
      while (center_ < graph_->adjacency_.size()) {
        while (outid_ < graph_->adjacency_[center_].size()) {
          if (center_ < graph_->adjacency_[center_][outid_].first) {
            return;
          }
          ++outid_;
        }
        ++center_;
        outid_ = 0;
      }
      return;
    }

   private:
     /** Store the pointer to graph and the edgeindex to get access to the edges
      *  and iterate through.
      */
     Graph* graph_;
     size_type center_;
     size_type outid_;

     /** Construct a valid EdgeIterator. */
     EdgeIterator(const Graph* graph, size_type center, size_type outid)
       : graph_(const_cast<Graph*>(graph)), center_(center), outid_(outid){}
    friend class Graph;
  };

  // Iterate the first edge.
  edge_iterator edge_begin() const {
    edge_iterator ei = EdgeIterator(this, 0, 0);
    return ei;
  }

  // Iterate the last edge.
  edge_iterator edge_end() const {
    edge_iterator ei = EdgeIterator(this, this->adjacency_.size(), 0);
    return ei;
  }

  //
  // Incident Iterator
  //

  /** @class Graph::IncidentIterator
   * @brief Iterator class for edges incident to a node. A forward iterator. */
  class IncidentIterator: private totally_ordered <IncidentIterator> {
   public:
    // These type definitions help us use STL's iterator_traits.
    /** Element type. */
    typedef Edge value_type;
    /** Type of pointers to elements. */
    typedef Edge* pointer;
    /** Type of references to elements. */
    typedef Edge& reference;
    /** Iterator category. */
    typedef std::input_iterator_tag iterator_category;
    /** Difference between iterators */
    typedef std::ptrdiff_t difference_type;

    /** Construct an invalid IncidentIterator. */
    IncidentIterator() {
    }

    // HW1 #5: YOUR CODE HERE
    // Supply definitions AND SPECIFICATIONS for:
    Edge operator*() const {
      return Edge(graph_, center_, outid_);
    }

    IncidentIterator& operator++() {
      ++outid_;
      return *this;
    }

    bool operator==(const IncidentIterator& ii) const {
      return (graph_ == ii.graph_ and center_ == ii.center_ and outid_ == ii.outid_);
    }

   private:
    friend class Graph;

    /** Store the node1_id and node2_id for the edge, where node1 is the
     *  central node and node2 are the adjacent nodes to node1.
     */
    graph_type* graph_;
    size_type center_;
    size_type outid_;

    /** Construct a valid IncidentIterator. */
    IncidentIterator(const Graph* graph, size_type center, size_type outid)
      : graph_(const_cast<Graph*>(graph)), center_(center), outid_(outid){};
  };

 private:


   /** A vector of all nodes' positions, and the element's index is the node_id_. */
   std::vector<std::pair<Point*, node_value_type>> nodes_;

   /** A vector storing every node's connections. */
   std::vector<std::vector<std::pair<size_type, edge_value_type>>> adjacency_;

};

#endif // CME212_GRAPH_HPP
