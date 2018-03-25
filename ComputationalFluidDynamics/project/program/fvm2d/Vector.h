#ifndef VECTOR_H_
#define VECTOR_H_

#include <algorithm>

template <typename T, int N >
class Vector{
public:
  typedef T value_type;
  typedef int size_type;
  typedef T* iterator;
  typedef const T* const_iterator;
  typedef T& reference;
  typedef const T& const_reference;

private:
  value_type val[N];
  size_type len = N;

public:
  Vector() = default;
  Vector(const T& _v) { 
    for(auto& i : val) i = _v;
  }

  Vector(const Vector<T,N>& _v) {
    (*this) = _v;
  }

  Vector<T,N>& operator=(const Vector<T,N>& _v){
    std::copy(_v.val, _v.val+N, val);
    return *this;
  }

  Vector<T,N>& operator=(const T& _v){
    for(auto& i : val) i = _v;
    return (*this);
  }

  size_type size(){ return len;}
  size_type size() const { return len;}

  inline value_type& operator[](size_type index) { return val[index];} 
  inline const value_type& operator[](size_type index) const { return
    val[index];}

  Vector<T,N>& operator+=(const Vector<T,N> & _vec){
    for(size_type i = 0; i < len; ++i)
      val[i] += _vec.val[i];
    return *this;
  }

  Vector<T,N>& operator-=(const Vector<T,N> & _vec){
    for(size_type i = 0; i < len; ++i)
      val[i] -= _vec.val[i];
    return *this;
  }

  Vector<T,N>& operator*=(const value_type& factor){
    for(size_type i = 0; i < len; ++i)
      val[i] *= factor;
    return *this;
  }

  Vector<T,N>& operator/=(const value_type& factor){
    for(size_type i = 0; i < len; ++i)
      val[i] /= factor;
    return *this;
  }

  iterator begin() { return val; }
  const iterator begin() const { return val;}

  iterator end() { return val + N; }
  const iterator end() const { return val + N;}

};

template <typename ostream, typename T, int N>
ostream& operator<<(ostream& os, const Vector<T,N>& vec){
  for(int i = 0; i < vec.size(); ++i)
    os << vec[i] << " ";
  return os;
}

template <typename T, int N>
Vector<T,N> operator+(const Vector<T,N>& v1, const Vector<T,N>& v2){
  Vector<T,N> v3;
  v3 = v1;
  v3 += v2;
  return v3;
}

template <typename T, int N>
Vector<T,N> operator-(const Vector<T,N>& v1, const Vector<T,N>& v2){
  Vector<T,N> v3;
  v3 = v1;
  v3 -= v2;
  return v3;
}

template <typename T, int N>
Vector<T,N> operator*(const Vector<T,N>& v1, const T& factor){
  Vector<T,N> v3;
  v3 = v1;
  v3 *= factor;
  return v3;
}

template <typename T, int N>
Vector<T,N> operator*(const T& factor, const Vector<T,N>& v1){
  Vector<T,N> v3;
  v3 = v1;
  v3 *= factor;
  return v3;
}

template <typename T, int N>
Vector<T,N> operator/(const Vector<T,N>& v1, const T& factor){
  Vector<T,N> v3;
  v3 = v1;
  v3 /= factor;
  return v3;
}

#endif
