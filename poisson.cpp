/**
 * @file poisson.cpp
 * Test script for treating the Graph as a MTL Matrix
 * and solving a Poisson equation.
 *
 * @brief Reads in two files specified on the command line.
 * First file: 3D Points (one per line) defined by three doubles.
 * Second file: Eges (one per line) defined by 2 indices into the point list
 *              of the first file.
 *
 * Launches an SDLViewer to visualize the solution.
 */

#include <fstream>

#include "CME212/SDLViewer.hpp"
#include "CME212/Util.hpp"
#include "CME212/Point.hpp"
#include "CME212/BoundingBox.hpp"
#include <boost/numeric/mtl/mtl.hpp>
#include <boost/numeric/itl/itl.hpp>

#include "Graph.hpp"

typedef Graph<char,char> GraphType;
typedef GraphType::node_type NodeType;

/** Another implementation of boundary. */
// bool boundary(const NodeType& n, double h) {
//   std::set<size_t> neighbor;
//   neighbor.emplace(n.index());
//   int k = 0;
//   for(auto i = n.edge_begin(); i != n.edge_end(); ++i){
//     auto m = (*i).node2();
//     if(neighbor.find(m.index()) == neighbor.end()){
//       neighbor.emplace(m.index());
//       if(norm(m.position() - n.position()) < 1.5 * h){
//         ++k;
//       }
//     }
//     for(auto j = m.edge_begin(); j != m.edge_end(); ++j){
//       auto l = (*j).node2();
//       if(neighbor.find(l.index()) == neighbor.end()){
//         neighbor.emplace(l.index());
//         if(norm(l.position() - n.position()) < 1.5 * h){
//           ++k;
//         }
//       }
//     }
//   }
//   std::cout << k << std::endl;
//   if(k < 8) return true;
//   return false;
// }

/* Check if the node is on the boundary, */
bool boundary(const NodeType& n) {
  auto p = n.position();
  if(norm_inf(p) == 1) return true;
  if(norm_inf(p - Point(0.6, 0.6, 0)) < 0.2 or norm_inf(p - Point(-0.6, 0.6, 0)) < 0.2
    or norm_inf(p - Point(0.6, -0.6, 0)) < 0.2 or norm_inf(p - Point(-0.6, -0.6, 0)) < 0.2) return true;
  Box3D b(Point(-0.6, -0.2, -1), Point(0.6, 0.2, 1));
  if(b.contains(p)) return true;
  return false;
}

/** A color function. */
struct PositionColorFn{
  PositionColorFn(){
  }
  CME212::Color operator()(GraphType::node_type n){
    double nm = norm(n.position());
    double x = n.position().z;
    double c = std::abs(x) / nm;
    return CME212::Color::make_heat(c);
  }
};

/* The function which returns the Point value according to the Vector x */
template<typename Vector>
struct Positionfunction{
  Vector x_;
  Positionfunction(Vector x) : x_(x){
  }
  Point& operator()(GraphType::node_type n){
    size_t i = n.index();
    Point& p = n.position();
    p.z = x_[i];
    return p;
  }
};

/* GraphSymmetricMatrix which is implemented by graph. */
class GraphSymmetricMatrix {

public:
  GraphSymmetricMatrix(GraphType* graph)
  // Initialize the private attributes of a graph,with no nodes and no edges.
    : graph_(graph){
  }

  double element(size_t i, size_t j){
    // double h = (*(graph_->edge_begin())).length();
    if(i == j and boundary(graph_->node(i))) return 1.0;
    if(i != j and (boundary(graph_->node(i)) or boundary(graph_->node(j)))) return 0;
    else{
      if(i == j) return -double(graph_->node(i).degree());
      auto n1 = graph_->node(i);
      auto n2 = graph_->node(j);
      if(graph_->has_edge(n1, n2)) return 1.0;
      else return 0.0;
    }
  }

  /** Helper function to perform multiplication . Allows for delayed
    * evaluation of results .
    * Assign :: apply (a , b ) resolves to an assignment operation such as
    * a += b , a -= b , or a = b .
    * @pre @a size ( v ) == size ( w ) */
  template <typename VectorIn, typename VectorOut, typename Assign>
  void mult(const VectorIn& v, VectorOut& w, Assign) const{
    size_t s = graph_->size();
    assert(mtl::size(v) == s);
    assert(mtl::size(w) == s);
    for(size_t i = 0; i < s; ++i){
      double current = 0;
      for (size_t j = indp_[i]; j < indp_[i+1]; ++j){
  			current += v[indi_[j]] * elem_[j];
  		}
      Assign::apply(w[i], current);
    }
  }

  /* * Matvec forwards to MTL ’s lazy m at _ c ve c _ mu l t ip l i er operator */
  template <typename Vector>
  mtl::vec::mat_cvec_multiplier<GraphSymmetricMatrix, Vector>
  operator *(const Vector& v) const{
    return {*this, v};
  }
  /* Construct a sparse matrix using the graph's feature. */
  void tosparse(){
    assert(graph_ != NULL);
    indp_.push_back(0);
    for(size_t i = 0; i < num_rows(); ++i){
      for(size_t j = 0; j < num_cols(); ++j){
        if(this->element(i, j) != 0){
          elem_.push_back(this->element(i, j));
          indi_.push_back(j);
        }
      }
      indp_.push_back(indi_.size());
    }
  }

  size_t num_rows() const{
    return graph_->size();
  }

  size_t num_cols() const{
    return graph_->size();
  }

  size_t m_size() const{
    return graph_->size() * graph_->size();
  }

private:
  GraphType* graph_;
  std::vector<double> elem_;
  std::vector<size_t> indp_;
  std::vector<size_t> indi_;
};

/* * The number of rows in the matrix . */
inline std::size_t num_rows(const GraphSymmetricMatrix& A){
  return A.num_rows();
}

/* * The number of columns in the matrix . */
inline std::size_t num_cols(const GraphSymmetricMatrix& A){
  return A.num_cols();
}

/* * The number of elements in the matrix . */
inline std::size_t size(const GraphSymmetricMatrix& A){
  return A.m_size();
}

/* * Traits that MTL uses to determine properties of our IdentityMatrix . */
namespace mtl {
  namespace ashape {
    /* * Define IdentityMatrix to be a non - scalar type . */
    template <>
      struct ashape_aux <GraphSymmetricMatrix> {
      typedef nonscal type ;
      };
  } // end namespace ashape

  /* * IdentityMatrix implements the Collection concept
  * with value_type and size_type */
  template <>
  struct Collection <GraphSymmetricMatrix> {
    typedef double value_type ;
    typedef unsigned size_type ;
  };
} // end namespace mtl

namespace itl{
  template <class Real, typename Viewer, typename Vector>
  class visual_iteration : public cyclic_iteration<Real> {

    typedef cyclic_iteration<Real> super;
    typedef visual_iteration self;

    void visual_iter(){

      // dispaly during iteration
      auto node_map = viewer_->empty_node_map(*graph_);
      viewer_->clear();
      viewer_->add_nodes(graph_->node_begin(), graph_->node_end(), pcf_, Positionfunction<Vector>(*x_), node_map);
      viewer_->add_edges(graph_->edge_begin(), graph_->edge_end(), node_map);
      viewer_->set_label(this->i);
      viewer_->center_view();
    }

  public:

    visual_iteration(const Vector& r0, int max_iter_, Real tol_,
       const Viewer* viewer, const GraphType* graph,
       const PositionColorFn pcf , const Vector* x, Real atol_ = Real(0), int cycle = 10)
     :super(r0, max_iter_, tol_, atol_, cycle), viewer_(const_cast<Viewer*>(viewer)),
      graph_(const_cast<GraphType*>(graph)), pcf_(pcf), x_(const_cast<Vector*>(x)){
        viewer_->launch();
        visual_iter();
      }

    bool finished() {
      visual_iter();
      return super::finished();
    }

    template <typename T>
    bool finished(const T& r)
    {
       bool ret= super::finished(r);
       visual_iter();
       return ret;
    }
  private:
    Viewer* viewer_;
    GraphType* graph_;
    PositionColorFn pcf_;
    Vector* x_;
  };
}


/** Remove all the nodes in graph @a g whose posiiton is within Box3D @a bb.
 * @post For all i, 0 <= i < @a g.num_nodes(),
 *        not bb.contains(g.node(i).position())
 */
void remove_box(GraphType& g, const Box3D& bb) {
  for(auto i = g.node_begin(); i != g.node_end(); ++i){
    auto n = *i;
    if(bb.contains(n.position())){
      g.remove_node(*i);
    }
  }
  return;
}

/* The function g. */
double g(Point& p){
  if(norm_inf(p) == 1) return 0.0;
  if(norm_inf(p - Point(0.6, 0.6, 0)) < 0.2 or norm_inf(p - Point(-0.6, -0.6, 0)) < 0.2
     or norm_inf(p - Point(0.6, -0.6, 0)) < 0.2 or norm_inf(p - Point(-0.6, 0.6, 0)) < 0.2 ) return -0.2;
  Box3D b(Point(-0.6, -0.2, -1), Point(0.6, 0.2, 1));
  if(b.contains(p)) return 1.0;
  return 100.0;
}

/* The function f. */
double f(Point& p){
  return 5 * cos(norm_1(p));
}


int main(int argc, char** argv)
{
  // Check arguments
  if (argc < 2) {
    std::cerr << "Usage: " << argv[0] << " NODES_FILE TETS_FILE\n";
    exit(1);
  }

  // Define an empty Graph
  GraphType graph;

  // Create a nodes_file from the first input argument
  std::ifstream nodes_file(argv[1]);
  // Interpret each line of the nodes_file as a 3D Point and add to the Graph
  std::vector<typename GraphType::node_type> node_vec;
  Point p;
  while (CME212::getline_parsed(nodes_file, p))
    node_vec.push_back(graph.add_node(2*p - Point(1,1,0)));

  // Create a tets_file from the second input argument
  std::ifstream tets_file(argv[2]);
  // Interpret each line of the tets_file as four ints which refer to nodes
  std::array<int,4> t;
  while (CME212::getline_parsed(tets_file, t)) {
    graph.add_edge(node_vec[t[0]], node_vec[t[1]]);
    graph.add_edge(node_vec[t[0]], node_vec[t[2]]);
    graph.add_edge(node_vec[t[1]], node_vec[t[3]]);
    graph.add_edge(node_vec[t[2]], node_vec[t[3]]);
  }

  // Get the edge length, should be the same for each edge
  auto it = graph.edge_begin();
  assert(it != graph.edge_end());
  double h = norm((*it).node1().position() - (*it).node2().position());
  // Make holes in our Graph
  remove_box(graph, Box3D(Point(-0.8+h,-0.8+h,-1), Point(-0.4-h,-0.4-h,1)));
  remove_box(graph, Box3D(Point( 0.4+h,-0.8+h,-1), Point( 0.8-h,-0.4-h,1)));
  remove_box(graph, Box3D(Point(-0.8+h, 0.4+h,-1), Point(-0.4-h, 0.8-h,1)));
  remove_box(graph, Box3D(Point( 0.4+h, 0.4+h,-1), Point( 0.8-h, 0.8-h,1)));
  remove_box(graph, Box3D(Point(-0.6+h,-0.2+h,-1), Point( 0.6-h, 0.2-h,1)));

  // HW3: YOUR CODE HERE
  size_t node_num = graph.size();
  mtl::vec::dense_vector<double> b(node_num, 0.0);
  int kk = 0;
  for(size_t i = 0; i < node_num; ++i){
    auto n = graph.node(i);
    auto x = n.position();
    if(boundary(n)){
      b[i] = g(x);
      ++kk;
    }
    else{
      double gxi = h * h * f(x);
      for(auto ni = n.edge_begin(); ni != n.edge_end(); ++ni){
        auto n2 = (*ni).node2();
        if(boundary(n2)) gxi -= g(n2.position());
      }
      b[i] = gxi;
    }
  }
  std::cout << kk << std::endl;
  std::cout << "norm(b) = " << mtl::two_norm(b) << std::endl;

  GraphSymmetricMatrix A(&graph);
  A.tosparse();
  mtl::vec::dense_vector<double> x(node_num, 0.0);
  //itl::cyclic_iteration<double> iter(b, 1000, 1.e-11, 0.0, 10);
  CME212::SDLViewer viewer;
  itl::visual_iteration<double, CME212::SDLViewer, mtl::vec::dense_vector<double>>
    iter(b, 1000, 1.e-11, &viewer, &graph, PositionColorFn(), &x);
  itl::cg(A, x, b, iter);

  return 0;
}
