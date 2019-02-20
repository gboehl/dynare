/* $Id: monoms.h 148 2005-04-19 15:12:26Z kamenik $ */
/* Copyright 2004, Ondra Kamenik */

#ifndef MONOMS_H
#define MONOMS_H

#include <vector>
#include <random>
#include <memory>

#include "int_sequence.hh"
#include "gs_tensor.hh"
#include "t_container.hh"
#include "sparse_tensor.hh"
#include "Vector.hh"

class IntGenerator
{
  int maxim{5};
  double probab{0.3};
  std::mt19937 mtgen;
  std::uniform_real_distribution<> dis;
public:
  IntGenerator() = default;
  void init(int nf, int ny, int nv, int nw, int nu, int mx, double prob);
  int get();
};

extern IntGenerator intgen;

class Monom : public IntSequence
{
public:
  Monom(int len); // generate a random monom
  Monom(int len, int item); // generate monom whose items are the given item
  double deriv(const IntSequence &vars) const;
  // this = this*m^ex (in monomial sense)
  void multiplyWith(int ex, const Monom &m);
  void print() const;
};

class Monom2Vector;
class Monom1Vector
{
  friend class Monom2Vector;
  int nx;
  int len;
  std::vector<Monom> x;
public:
  Monom1Vector(int nxx, int l);
  ~Monom1Vector() = default;
  void deriv(const IntSequence &c, Vector &out) const;
  std::unique_ptr<FGSTensor> deriv(int dim) const;
  void print() const;
};

//class Monom3Vector;
class Monom2Vector
{
  int ny, nu;
  int len;
  std::vector<Monom> y, u;
public:
  // generate random vector of monom two vector
  Monom2Vector(int nyy, int nuu, int l);
  // calculate g(x(y,u))
  Monom2Vector(const Monom1Vector &g, const Monom2Vector &xmon);
  ~Monom2Vector() = default;
  void deriv(const Symmetry &s, const IntSequence &c, Vector &out) const;
  std::unique_ptr<FGSTensor> deriv(const Symmetry &s) const;
  FGSContainer *deriv(int maxdim) const;
  void print() const;
};

class Monom4Vector
{
  int len;
  int nx1, nx2, nx3, nx4;
  std::vector<Monom> x1, x2, x3, x4;
public:
  /* random for g(y,u,sigma) */
  Monom4Vector(int l, int ny, int nu);
  /* random for G(y,u,u',sigma) */
  Monom4Vector(int l, int ny, int nu, int nup);
  /* random for f(y+,y,y-,u) */
  Monom4Vector(int l, int nbigg, int ng, int ny, int nu);
  /* substitution f(G(y,u,u',sigma),g(y,u,sigma),y,u) */
  Monom4Vector(const Monom4Vector &f, const Monom4Vector &bigg,
               const Monom4Vector &g);
  ~Monom4Vector() = default;
  FSSparseTensor *deriv(int dim) const;
  std::unique_ptr<FGSTensor> deriv(const Symmetry &s) const;
  void deriv(const Symmetry &s, const IntSequence &coor, Vector &out) const;
  void print() const;
protected:
  void init_random();
};

struct SparseDerivGenerator
{
  int maxdimen;
  FGSContainer bigg;
  FGSContainer g;
  FGSContainer rcont;
  FSSparseTensor **const ts;
  SparseDerivGenerator(int nf, int ny, int nu, int nup, int nbigg, int ng,
                       int mx, double prob, int maxdim);
  ~SparseDerivGenerator();
};

struct DenseDerivGenerator
{
  int maxdimen;
  FGSContainer *xcont;
  FGSContainer *rcont;
  std::vector<std::unique_ptr<FGSTensor>> ts;
  UGSContainer *uxcont;
  std::vector<std::unique_ptr<UGSTensor>> uts;
  DenseDerivGenerator(int ng, int nx, int ny, int nu,
                      int mx, double prob, int maxdim);
  void unfold();
  ~DenseDerivGenerator();
};

#endif
