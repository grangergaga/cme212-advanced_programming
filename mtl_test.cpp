/**
 * @file mtl_test.cpp
 * Test script for interfacing with MTL4 and it's linear solvers.
 */

// HW3: Need to install/include Boost and MTL in Makefile
#include <boost/numeric/mtl/mtl.hpp>
#include <boost/numeric/itl/itl.hpp>

// HW3: YOUR CODE HERE
// Define a IdentityMatrix that interfaces with MTL

class IdentityMatrix {
  /* * Compute the product of a vector with this identity matrix
  */
public:

  IdentityMatrix(size_t size)
  // Initialize the private attributes of a graph,with no nodes and no edges.
    : size_(size){

  }

  /** Helper function to perform multiplication . Allows for delayed
    * evaluation of results .
    * Assign :: apply (a , b ) resolves to an assignment operation such as
    * a += b , a -= b , or a = b .
    * @pre @a size ( v ) == size ( w ) */
  template <typename VectorIn, typename VectorOut, typename Assign>
  void mult(const VectorIn& v, VectorOut& w, Assign) const{
    assert(mtl::size(v) == size_);
    assert(mtl::size(w) == size_);
    for(size_t i = 0; i < size_; ++i){
      Assign::apply(w[i], v[i]);
    }
  }

  /* * Matvec forwards to MTL ’s lazy m at _ c ve c _ mu l t ip l i er operator */
  template <typename Vector>
  mtl::vec::mat_cvec_multiplier<IdentityMatrix, Vector>
  operator *(const Vector& v) const{
    return {*this, v};
  }

  size_t num_rows() const{
    return size_;
  }

  size_t num_cols() const{
    return size_;
  }

  size_t m_size() const{
    return size_ * size_;
  }

private:
  size_t size_;
};

/* * The number of rows in the matrix . */
inline std::size_t num_rows(const IdentityMatrix& A){
  return A.num_rows();
}

/* * The number of columns in the matrix . */
inline std::size_t num_cols(const IdentityMatrix& A){
  return A.num_cols();
}

/* * The number of elements in the matrix . */
inline std::size_t size(const IdentityMatrix& A){
  return A.m_size();
}

/* * Traits that MTL uses to determine properties of our IdentityMatrix . */
namespace mtl {
  namespace ashape {
    /* * Define IdentityMatrix to be a non - scalar type . */
    template <>
      struct ashape_aux <IdentityMatrix> {
      typedef nonscal type ;
      };
  } // end namespace ashape

  /* * IdentityMatrix implements the Collection concept
  * with value_type and size_type */
  template <>
  struct Collection <IdentityMatrix> {
    typedef double value_type ;
    typedef unsigned size_type ;
  };
} // end namespace mtl

int main()
{
  typedef mtl::vec::dense_vector<double> Vector;
  typedef IdentityMatrix IMatrix;
  const size_t N = 1000;
  IMatrix I(N);

  //Create and Identity Preconditioner (there is no conditioning)
  itl::pc::identity<IMatrix, double> P(I);

  //Set b such that x == 1 is the solution; start with x == 0
  Vector x(N, 1.0), b(N);
  b = I * x;
  x = 0.0;
  //Termination criterion: r < 1e-6 * b or N iterations
  itl::noisy_iteration<double> iter(b, 500, 1.e-6);

  //Solve Ax == b with left preconditioner P
  itl::bicgstab(I, x, b, P, iter);
  return 0;
}
