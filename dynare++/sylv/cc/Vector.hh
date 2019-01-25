/* $Header: /var/lib/cvs/dynare_cpp/sylv/cc/Vector.h,v 1.1.1.1 2004/06/04 13:01:13 kamenik Exp $ */

/* Tag $Name:  $ */

#ifndef VECTOR_H
#define VECTOR_H

/* NOTE! Vector and ConstVector have not common super class in order
 * to avoid running virtual method invokation mechanism. Some
 * members, and methods are thus duplicated */

#include <array>
#include <memory>
#include <complex>

#if defined(MATLAB_MEX_FILE) || defined(OCTAVE_MEX_FILE)
# include <dynmex.h>
#endif

class GeneralMatrix;
class ConstVector;

class Vector
{
  friend class ConstVector;
protected:
  int len{0};
  int off{0}; // offset to double* pointer
  int s{1}; // stride (also called "skip" in some places)
  std::shared_ptr<double> data;
public:
  Vector() = default;
  Vector(int l) : len{l}, data{new double[l], [](double *arr) { delete[] arr; }}
  {
  }
  Vector(Vector &v) = default;
  Vector(const Vector &v);
  Vector(Vector &&v) = default;
  // We don't want implict conversion from ConstVector, since it's expensive
  explicit Vector(const ConstVector &v);
  Vector(std::shared_ptr<double> d, int l)
    : len(l), data{std::move(d)}
  {
  }
  Vector(Vector &v, int off_arg, int l);
  Vector(const Vector &v, int off_arg, int l);
  Vector(Vector &v, int off_arg, int skip, int l);
  Vector(const Vector &v, int off_arg, int skip, int l);
#if defined(MATLAB_MEX_FILE) || defined(OCTAVE_MEX_FILE)
  explicit Vector(mxArray *p);
#endif
  Vector &operator=(const Vector &v);
  Vector &operator=(Vector &&v);
  Vector &operator=(const ConstVector &v);
  double &
  operator[](int i)
  {
    return data.get()[off+s*i];
  }
  const double &
  operator[](int i) const
  {
    return data.get()[off+s*i];
  }
  const double *
  base() const
  {
    return data.get() + off;
  }
  double *
  base()
  {
    return data.get() + off;
  }
  int
  length() const
  {
    return len;
  }
  int
  skip() const
  {
    return s;
  }

  /** Exact equality. */
  bool operator==(const Vector &y) const;
  bool operator!=(const Vector &y) const;
  /** Lexicographic ordering. */
  bool operator<(const Vector &y) const;
  bool operator<=(const Vector &y) const;
  bool operator>(const Vector &y) const;
  bool operator>=(const Vector &y) const;

  virtual ~Vector() = default;
  void zeros();
  void nans();
  void infs();
  void rotatePair(double alpha, double beta1, double beta2, int i);
  // Computes this += r*v
  void add(double r, const Vector &v);
  // Computes this += r*v
  void add(double r, const ConstVector &v);
  // Computes this += z*v (where this and v are intepreted as complex vectors)
  void addComplex(const std::complex<double> &z, const Vector &v);
  // Computes this += z*v (where this and v are intepreted as complex vectors)
  void addComplex(const std::complex<double> &z, const ConstVector &v);
  void mult(double r);
  double getNorm() const;
  double getMax() const;
  double getNorm1() const;
  double dot(const Vector &y) const;
  bool isFinite() const;
  void print() const;

  /* multiplies | alpha -beta1|           |b1|   |x1|
     |             |\otimes I .|  | = |  |
     | -beta2 alpha|           |b2|   |x2|
  */
  static void mult2(double alpha, double beta1, double beta2,
                    Vector &x1, Vector &x2,
                    const Vector &b1, const Vector &b2);
  /* same, but adds instead of set */
  static void mult2a(double alpha, double beta1, double beta2,
                     Vector &x1, Vector &x2,
                     const Vector &b1, const Vector &b2);
  /* same, but subtracts instead of add */
  static void
  mult2s(double alpha, double beta1, double beta2,
         Vector &x1, Vector &x2,
         const Vector &b1, const Vector &b2)
  {
    mult2a(-alpha, -beta1, -beta2, x1, x2, b1, b2);
  }
private:
  void copy(const double *d, int inc);
};

class ConstGeneralMatrix;

class ConstVector
{
  friend class Vector;
protected:
  int len;
  int off{0}; // offset to double* pointer
  int s{1}; // stride (also called "skip" in some places)
  std::shared_ptr<const double> data;
public:
  // Implicit conversion from Vector is ok, since it's cheap
  ConstVector(const Vector &v);
  ConstVector(const ConstVector &v) = default;
  ConstVector(ConstVector &&v) = default;
  ConstVector(std::shared_ptr<const double> d, int l) : len{l}, data{std::move(d)}
  {
  }
  ConstVector(const ConstVector &v, int off_arg, int l);
  ConstVector(const ConstVector &v, int off_arg, int skip, int l);
  ConstVector(std::shared_ptr<const double> d, int skip, int l);
#if defined(MATLAB_MEX_FILE) || defined(OCTAVE_MEX_FILE)
  explicit ConstVector(const mxArray *p);
#endif
  virtual ~ConstVector() = default;
  ConstVector &operator=(const ConstVector &v) = delete;
  ConstVector &operator=(ConstVector &&v) = delete;
  const double &
  operator[](int i) const
  {
    return data.get()[off+s*i];
  }
  const double *
  base() const
  {
    return data.get() + off;
  }
  int
  length() const
  {
    return len;
  }
  int
  skip() const
  {
    return s;
  }
  /** Exact equality. */
  bool operator==(const ConstVector &y) const;
  bool
  operator!=(const ConstVector &y) const
  {
    return !operator==(y);
  }
  /** Lexicographic ordering. */
  bool operator<(const ConstVector &y) const;
  bool
  operator<=(const ConstVector &y) const
  {
    return operator<(y) || operator==(y);
  }
  bool
  operator>(const ConstVector &y) const
  {
    return !operator<=(y);
  }
  bool
  operator>=(const ConstVector &y) const
  {
    return !operator<(y);
  }

  double getNorm() const;
  double getMax() const;
  double getNorm1() const;
  double dot(const ConstVector &y) const;
  bool isFinite() const;
  void print() const;
};

class ZeroPad
{
public:
  static const int length = 16;
private:
  std::array<double, length> pad;
public:
  ZeroPad();
  const double *
  getBase() const
  {
    return pad.data();
  }
};

#endif /* VECTOR_H */
