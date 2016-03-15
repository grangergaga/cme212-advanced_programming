/** @file iterator_adaptor_thrust.cpp
 * @brief A simple example of using Thrust's iterator adaptors.
 * Compare these implementations with those of hand-written iterators and/or
 * hand-written iterator-adpators!
 *
 * Thrust uses another level of indirection (facade and iterator_core_access)
 * to avoid any conflicts with overloading. The user can optionally implement
 *  void increment()
 *  void decrement()
 *  void advance(DifferenceType)
 *  Reference dereference() const
 *  bool equal(OtherIt) const
 *  DifferenceType distance_to(OtherIt) const
 *
 * Which are then used to implement the iterator interface.
 *
 * That said, when possible, it is usually easier to use and compose
 * the predesigned "fancy iterators" like
 *   transform_iterator
 *   counting_iterator
 *   constant_iterator
 *   permutation_iterator
 *   reverse_iterator
 *   zip_iterator
 *
 * For more information and examples, see
 * https://thrust.github.io/doc/classthrust_1_1iterator__adaptor.html
 * https://github.com/thrust/thrust/tree/master/thrust/iterator
 *
 * INSTALL Thrust with
 *   $ git clone git@github.com:thrust/thrust.git
 * and -I./thrust in your Makefile
 */

#include <thrust/iterator/iterator_adaptor.h>
#include <thrust/iterator/counting_iterator.h>
using thrust::iterator_adaptor;

#include <vector>
#include <numeric>


/** Derived class implementing a strided iterator (random access).
 */
struct my_stride_iterator
    : public iterator_adaptor<my_stride_iterator,                // Derived type
                              std::vector<int>::iterator> {      // Base type
  my_stride_iterator(const std::vector<int>::iterator& m)
      : iterator_adaptor(m) {
  }
 private:
  friend class thrust::iterator_core_access;
  void increment() {
    base_reference() += 2;
  }
  void decrement() {
    base_reference() -= 2;
  }
  void advance(typename my_stride_iterator::difference_type n) {
    base_reference() += 2*n;
  }
  typename my_stride_iterator::difference_type
  distance_to(const my_stride_iterator& it) const {
    return (it.base() - base()) / 2;
  }
};


/** Derived class implementing a transformed return value (random access).
 * NOTE: Because the elements are implicit, the dereference returns
 * a value rather a reference. Thus, this is a read-only iterator and
 * we need to specify the Reference type in the adaptor to be a value.
 */
struct my_transform_iterator
    : public iterator_adaptor<my_transform_iterator,             // Derived type
                              std::vector<int>::iterator,        // Base type
                              int,                               // Value type
                              thrust::use_default,               // System
                              std::random_access_iterator_tag,   // Category
                              int> {                             // Reference
  my_transform_iterator(const std::vector<int>::iterator& m)
      : iterator_adaptor(m) {
  }
 private:
  friend class thrust::iterator_core_access;
  int dereference() const {
    return (*base()) * (*base());
  }
};

/** Derived class implementing an iterator with an index into container.
 */
struct my_position_iterator
    : public iterator_adaptor<my_position_iterator,              // Derived type
                              thrust::counting_iterator<int>> {  // Base type
  my_position_iterator(std::vector<int>& c, int pos)
      : iterator_adaptor(thrust::counting_iterator<int>(pos)), c_(c) {
  }
 private:
  friend class thrust::iterator_core_access;
  int& dereference() const {
    return c_[*base()];
  }
  std::vector<int>& c_;
};



#include <algorithm>
#include <iostream>

int main()
{
  std::vector<int> a = {0,2,4,6,1,2,3,4};
  std::copy(a.begin(), a.end()-1, std::ostream_iterator<int>(std::cout, ", "));
  std::cout << *(a.end()-1)  << std::endl;

  auto begin = a.begin();
  auto end   = a.end();

  std::cout << "Normal sum:  " << std::accumulate(begin, end, 0) << " = ";
  std::copy(begin, end-1, std::ostream_iterator<int>(std::cout, " + "));
  std::cout << *(end-1) << std::endl;

  auto pbegin = my_position_iterator(a, 0);
  auto pend   = my_position_iterator(a, a.size());

  std::cout << "By position: " << std::accumulate(pbegin, pend, 0) << " = ";
  std::copy(pbegin, pend-1, std::ostream_iterator<int>(std::cout, " + "));
  std::cout << *(pend-1) << std::endl;

  auto tbegin = my_transform_iterator(a.begin());
  auto tend   = my_transform_iterator(a.end());

  std::cout << "Squared sum: " << std::accumulate(tbegin, tend, 0) << " = ";
  std::copy(tbegin, tend-1, std::ostream_iterator<int>(std::cout, " + "));
  std::cout << *(tend-1) << std::endl;

  auto sbegin = my_stride_iterator(a.begin());
  auto send   = my_stride_iterator(a.end());

  std::cout << "Strided sum:  " << std::accumulate(sbegin, send, 0) << " = ";
  std::copy(sbegin, send-1, std::ostream_iterator<int>(std::cout, " + "));
  std::cout << *(send-1) << std::endl;

  return 0;
}
