#ifndef __COL_VECTOR_HPP__
#define __COL_VECTOR_HPP__

#include "matrix.hpp"

template <typename T> class col_vector : public matrix<T> {//Column vector
public:
  col_vector();
  col_vector(const unsigned int& n, const bool init_0=false);
  col_vector(const unsigned int& size, const T*);
  col_vector(const col_vector&);
  col_vector(const matrix<T>&);

  unsigned int size() const;
  T sum() const;
  col_vector<T>& set_to(const T &value=0);

  T& operator[](const unsigned int&) const;
};

#endif
