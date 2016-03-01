/**
 * @file mass_spring.cpp
 * Implementation of mass-spring system using Graph
 *
 * @brief Reads in two files specified on the command line.
 * First file: 3D Points (one per line) defined by three doubles
 * Second file: Tetrahedra (one per line) defined by 4 indices into the point
 * list
 */

#include <fstream>

#include "CME212/SDLViewer.hpp"
#include "CME212/Util.hpp"
#include "CME212/Color.hpp"
#include "CME212/Point.hpp"

#include "Graph.hpp"



// Gravity in meters/sec^2
static constexpr double grav = 9.81;

/** Custom structure of data to store with Nodes */
struct NodeData {
  Point vel;       //< Node velocity
  double mass;     //< Node mass
};

// Custom structure of data to store with Edges. */
struct EdgeData {
  double K;
  double L;
};

// HW2 #1 YOUR CODE HERE
// Define your Graph type
typedef Graph<NodeData, EdgeData> GraphType;
typedef typename GraphType::node_type Node;
typedef typename GraphType::edge_type Edge;


/** Change a graph's nodes according to a step of the symplectic Euler
 *    method with the given node force.
 * @param[in,out] g      Graph
 * @param[in]     t      The current time (useful for time-dependent forces)
 * @param[in]     dt     The time step
 * @param[in]     force  Function object defining the force per node
 * @return the next time step (usually @a t + @a dt)
 *
 * @tparam G::node_value_type supports ???????? YOU CHOOSE
 * @tparam F is a function object called as @a force(n, @a t),
 *           where n is a node of the graph and @a t is the current time.
 *           @a force must return a Point representing the force vector on Node
 *           at time @a t.
 */

/** Version with no constraints. */
template <typename G, typename F>
double symp_euler_step(G& g, double t, double dt, F force) {
  // Compute the t+dt position
  for (auto it = g.node_begin(); it != g.node_end(); ++it) {
    auto n = *it;
    if (n.position() == Point(0, 0, 0) or n.position() == Point(1, 0, 0)){
      n.value().vel = Point(0, 0, 0);
    }
    // Update the position of the node according to its velocity
    // x^{n+1} = x^{n} + v^{n} * dt
    n.position() += n.value().vel * dt;
  }
  // Compute the t+dt velocity
  for (auto it = g.node_begin(); it != g.node_end(); ++it) {
    auto n = *it;
    // v^{n+1} = v^{n} + F(x^{n+1},t) * dt / m
    n.value().vel += force(n, t) * (dt / n.value().mass);
  }

  return t + dt;
}

/** Version with constraints. */
template <typename G, typename F, typename C>
double symp_euler_step(G& g, double t, double dt, F force, C cons) {
  // Compute the t+dt position
  for (auto it = g.node_begin(); it != g.node_end(); ++it) {
    auto n = *it;
    if (n.position() == Point(0, 0, 0) or n.position() == Point(1, 0, 0)){
      n.value().vel = Point(0, 0, 0);
    }
    // Update the position of the node according to its velocity
    // x^{n+1} = x^{n} + v^{n} * dt
    n.position() += n.value().vel * dt;
  }
  // Compute the t+dt velocity
  for (auto it = g.node_begin(); it != g.node_end(); ++it) {
    auto n = *it;
    // v^{n+1} = v^{n} + F(x^{n+1},t) * dt / m
    n.value().vel += force(n, t) * (dt / n.value().mass);
  }
  cons(g);
  return t + dt;
}

/** Force function object for HW2 #1. */
struct Problem1Force {
  /** Return the force applying to @a n at time @a t.
   *
   * For HW2 #1, this is a combination of mass-spring force and gravity,
   * except that points at (0, 0, 0) and (1, 0, 0) never move. We can
   * model that by returning a zero-valued force. */
  /** We need to initialize K_ and L_ first before using this functor. */
  double K_;
  double L_;
  Problem1Force(const double K, const double L) : K_(K), L_(L) {
  };
  template <typename NODE>
  Point operator()(NODE n, double t) {
    (void) t;     // silence compiler warnings
    if(n.position() == Point(0, 0, 0) || n.position() == Point(1, 0, 0)){
      return Point(0, 0, 0);
    }
    Point Fspring = Point(0, 0, 0);
    Point xi = n.position();
    for(auto i = n.edge_begin(); i != n.edge_end(); ++i) {
      Point xj = (*i).node2().position();
      Fspring += (-K_ * (norm(xi - xj) - L_) / norm(xi - xj)) * (xi - xj);
    }
    return Fspring + n.value().mass * Point(0, 0, -grav);
  }
};

/** Force function object for HW2 #2, only the spring force. */
struct Problem2Force {
  /** Return the spring force applying to @a n at time @a t.
   *
   * For HW2 #2, this is a combination of mass-spring force and gravity,
   * except that points at (0, 0, 0) and (1, 0, 0) never move. We can
   * model that by returning a zero-valued force. */
  template <typename NODE>
  Point operator()(NODE n, double t) {
    (void) t;     // silence compiler warnings
    if(n.position() == Point(0, 0, 0) || n.position() == Point(1, 0, 0)){
      return Point(0, 0, 0);
    }
    Point Fspring = Point(0, 0, 0);
    Point xi = n.position();
    for(auto i = n.edge_begin(); i != n.edge_end(); ++i) {
      double K = (*i).value().K;
      double L = (*i).value().L;
      //std::cout << K << "  " << L << std::endl;
      Point xj = (*i).node2().position();
      double len = (*i).length();
      Fspring += (-K * (len - L) / len) * (xi - xj);
    }
    return Fspring + (n.value().mass * Point(0, 0, -grav));
  }
};


/** Force function object for HW2 #3, only the gravity force. */
struct GravityForce {
  /** Return the gravity force. */
  template <typename NODE>
  Point operator()(NODE n, double t) {
    (void) t;     // silence compiler warnings
    return  (n.value().mass * Point(0, 0, -grav));
  }
};

/** Force function object for HW2 #3, only the spring force. */
struct MassSpringForce {
  /** Return the mass spring force. */
  template <typename NODE>
  Point operator()(NODE n, double t) {
    (void) t;
    Point Fspring = Point(0, 0, 0);
    Point xi = n.position();
    for(auto i = n.edge_begin(); i != n.edge_end(); ++i) {
      double K = (*i).value().K;
      double L = (*i).value().L;
      //std::cout << K << "  " << L << std::endl;
      Point xj = (*i).node2().position();
      double len = (*i).length();
      Fspring += (-K * (len - L) / len) * (xi - xj);
    }
    return Fspring;
  }
};

/** Force function object for HW2 #3, only the damping force. */
struct DampingForce {
  double c_;
  DampingForce(const double c) : c_(c) {
  };
  /** Return the damping force. */
  template <typename NODE>
  Point operator()(NODE n, double t) {
    (void) n;
    (void) t;
    return -c_ * n.value().vel;
  }
};

/** Zero Force, used in make_combined when only two parameters are passed. */
struct ZeroForce {
  template <typename NODE>
  Point operator()(NODE n, double t) {
    (void) n;
    (void) t;
    return Point(0, 0, 0);
  }
};

template <typename f1, typename f2, typename f3>
struct make_combined_force {
  f1 f1_;
  f2 f2_;
  f3 f3_;
  /** Combine 3 forces. */
  make_combined_force(f1 force1, f2 force2, f3 force3):
    f1_(force1), f2_(force2), f3_(force3){};
  /** Combine 2 forces. */
  make_combined_force(f1 force1, f2 force2):
    f1_(force1), f2_(force2), f3_(ZeroForce()){};
  /** Return the combined force. */
  template<typename NODE>
  Point operator()(NODE n, double t) {
    return f1_(n, t) + f2_(n, t) + f3_(n, t);
  }
};

/** The constraint of a plane. */
struct PlaneConstraint {
  double l = -0.75;
  template<typename GRAPH>
  void operator()(GRAPH& g){
    for(auto it = g.node_begin(); it != g.node_end(); ++it){
      auto n = *it;
      if(n.position().z < l){
        n.position().z = l;
        n.value().vel.z = 0;
      }
    }
    return;
  }
};

/** The constraint of sphere. */
struct SphereConstraint {
  Point c = Point(0.5, 0.5, -0.5);
  double r = 0.15;
  template<typename GRAPH>
  void operator()(GRAPH& g){
    for(auto it = g.node_begin(); it != g.node_end(); ++it){
      auto n = *it;
      auto x = n.position();
      auto R = (x - c)/norm(x - c);
      if(norm(x - c) < r){
        n.value().vel -= (n.value().vel * R) * R;
        n.position() = c + r * R;
      }
    }
    return;
  }
};

/** The functor to remove nodes which violate the constraint. */
struct SphereRemove {
  Point c = Point(0.5, 0.5, -0.5);
  double r = 0.15;
  template<typename GRAPH>
  void operator()(GRAPH& g){
    auto it = g.node_begin();
    while(it != g.node_end()){
      auto n = *it;
      auto x = n.position();
      if(norm(x - c) < r){
        it = g.remove_node(it);
      }
      else{
        ++it;
      }
    }
    return;
  }
};


int main(int argc, char** argv) {
  // Check arguments
  if (argc < 3) {
    std::cerr << "Usage: " << argv[0] << " NODES_FILE TETS_FILE\n";
    exit(1);
  }

  // Construct an empty graph
  GraphType graph;

  // Create a nodes_file from the first input argument
  std::ifstream nodes_file(argv[1]);
  // Interpret each line of the nodes_file as a 3D Point and add to the Graph
  Point p;
  std::vector<typename GraphType::node_type> nodes;
  while (CME212::getline_parsed(nodes_file, p))
    nodes.push_back(graph.add_node(p));

  // Create a tets_file from the second input argument
  std::ifstream tets_file(argv[2]);
  // Interpret each line of the tets_file as four ints which refer to nodes
  std::array<int,4> t;
  while (CME212::getline_parsed(tets_file, t)) {
    graph.add_edge(nodes[t[0]], nodes[t[1]]);
    graph.add_edge(nodes[t[0]], nodes[t[2]]);
 //#if 0
    // Diagonal edges: include as of HW2 #2
    graph.add_edge(nodes[t[0]], nodes[t[3]]);
    graph.add_edge(nodes[t[1]], nodes[t[2]]);
//#endif
    graph.add_edge(nodes[t[1]], nodes[t[3]]);
    graph.add_edge(nodes[t[2]], nodes[t[3]]);
  }

  // Initialize the node values: mass and velocity.
  for(auto i = graph.node_begin(); i != graph.node_end(); ++i){
    (*i).value().vel = Point(0, 0, 0);
    (*i).value().mass = 1.0 / (double)graph.num_nodes();
  }

  // Initialize the edge values: K and L.
  for(auto i = graph.edge_begin(); i != graph.edge_end(); ++i){
    auto e = *i;
    auto edual = e.dual();
    e.value().K = 100;
    edual.value().K = 100;
    e.value().L = e.length();
    edual.value().L = edual.length();
  }

  // double K = 100;
  // double L = (*(graph.edge_begin())).length();

  // Print out the stats
  std::cout << graph.num_nodes() << " " << graph.num_edges() << std::endl;

  // Launch the SDLViewer
  CME212::SDLViewer viewer;
  auto node_map = viewer.empty_node_map(graph);
  viewer.launch();

  viewer.add_nodes(graph.node_begin(), graph.node_end(), node_map);
  viewer.add_edges(graph.edge_begin(), graph.edge_end(), node_map);

  viewer.center_view();

  // Begin the mass-spring simulation
  double dt = 0.001;
  double t_start = 0;
  double t_end = 5.0;

  for (double t = t_start; t < t_end; t += dt) {
    //std::cout << "t = " << t << std::endl;
    //symp_euler_step(graph, t, dt, Problem1Force(K, L));
    //symp_euler_step(graph, t, dt, Problem2Force());
    //symp_euler_step(graph, t, dt, make_combined_force<GravityForce, MassSpringForce, ZeroForce>(GravityForce(), MassSpringForce()));
    auto force = make_combined_force<GravityForce, MassSpringForce, DampingForce>(GravityForce(), MassSpringForce(), DampingForce(1.0/graph.size()));
    symp_euler_step(graph, t, dt, force, SphereRemove());
    // Update viewer with nodes' new positions
    viewer.clear();
    node_map.clear();
    viewer.add_nodes(graph.node_begin(), graph.node_end(), node_map);
    viewer.add_edges(graph.edge_begin(), graph.edge_end(), node_map);
    viewer.set_label(t);

    // These lines slow down the animation for small graphs, like grid0_*.
    // Feel free to remove them or tweak the constants.
    if (graph.size() < 100)
      CME212::sleep(0.001);
  }

  return 0;
}
