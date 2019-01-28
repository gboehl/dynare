// Copyright 2005, Ondra Kamenik

// Quadrature.

/* This file defines an interface for one dimensional (non-nested) quadrature
   |OneDQuadrature|, and a parent for all multi-dimensional
   quadratures. This parent class |Quadrature| presents a general concept of
   quadrature, this is
   $$\int f(x){\rm d}x \approx\sum_{i=1}^N w_ix_i$$
   The class |Quadrature| just declares this concept. The concept is
   implemented by class |QuadratureImpl| which paralelizes the
   summation. All implementations therefore wishing to use the parallel
   implementation should inherit from |QuadratureImpl| and integration is
   done.

   The integration concept relies on a point iterator, which goes through
   all $x_i$ and $w_i$ for $i=1,\ldots,N$. All the iterators must be able
   to go through only a portion of the set $i=1,\ldots,N$. This enables
   us to implement paralelism, for two threads for example, one iterator
   goes from the beginning to the (approximately) half, and the other
   goes from the half to the end.

   Besides this concept of the general quadrature, this file defines also
   one dimensional quadrature, which is basically a scheme of points and
   weights for different levels. The class |OneDQuadrature| is a parent
   of all such objects, the classes |GaussHermite| and |GaussLegendre|
   are specific implementations for Gauss--Hermite and Gauss--Legendre
   quadratures resp. */

#ifndef QUADRATURE_H
#define QUADRATURE_H

#include <iostream>
#include <fstream>
#include <iomanip>

#include "vector_function.hh"
#include "int_sequence.hh"
#include "sthread.hh"

/* This pure virtual class represents a concept of one-dimensional
   (non-nested) quadrature. So, one dimensional quadrature must return
   number of levels, number of points in a given level, and then a point
   and a weight in a given level and given order. */

class OneDQuadrature
{
public:
  virtual ~OneDQuadrature() = default;
  virtual int numLevels() const = 0;
  virtual int numPoints(int level) const = 0;
  virtual double point(int level, int i) const = 0;
  virtual double weight(int lelel, int i) const = 0;
};

/* This is a general concept of multidimensional quadrature. at this
   general level, we maintain only a dimension, and declare virtual
   functions for integration. The function take two forms; first takes a
   constant |VectorFunction| as an argument, creates locally
   |VectorFunctionSet| and do calculation, second one takes as an
   argument |VectorFunctionSet|.

   Part of the interface is a method returning a number of evaluations
   for a specific level. Note two things: this assumes that the number of
   evaluations is known apriori and thus it is not applicable for
   adaptive quadratures, second for Monte Carlo type of quadrature, the
   level is a number of evaluations. */

class Quadrature
{
protected:
  int dim;
public:
  Quadrature(int d) : dim(d)
  {
  }
  virtual ~Quadrature() = default;
  int
  dimen() const
  {
    return dim;
  }
  virtual void integrate(const VectorFunction &func, int level,
                         int tn, Vector &out) const = 0;
  virtual void integrate(VectorFunctionSet &fs, int level, Vector &out) const = 0;
  virtual int numEvals(int level) const = 0;
};

/* This is just an integration worker, which works over a given
   |QuadratureImpl|. It also needs the function, level, a specification
   of the subgroup of points, and output vector.

   See |@<|QuadratureImpl| class declaration@>| for details. */

template <typename _Tpit>
class QuadratureImpl;

template <typename _Tpit>
class IntegrationWorker : public sthread::detach_thread
{
  const QuadratureImpl<_Tpit> &quad;
  VectorFunction &func;
  int level;
  int ti;
  int tn;
  Vector &outvec;
public:
  IntegrationWorker(const QuadratureImpl<_Tpit> &q, VectorFunction &f, int l,
                    int tii, int tnn, Vector &out)
    : quad(q), func(f), level(l), ti(tii), tn(tnn), outvec(out)
  {
  }

  /* This integrates the given portion of the integral. We obtain first
     and last iterators for the portion (|beg| and |end|). Then we iterate
     through the portion. and finally we add the intermediate result to the
     result |outvec|.

     This method just everything up as it is coming. This might be imply
     large numerical errors, perhaps in future I will implement something
     smarter. */

  void
  operator()() override
  {
    _Tpit beg = quad.begin(ti, tn, level);
    _Tpit end = quad.begin(ti+1, tn, level);
    Vector tmpall(outvec.length());
    tmpall.zeros();
    Vector tmp(outvec.length());

    // note that since beg came from begin, it has empty signal
    // and first evaluation gets no signal
    for (_Tpit run = beg; run != end; ++run)
      {
        func.eval(run.point(), run.signal(), tmp);
        tmpall.add(run.weight(), tmp);
      }

    {
      sthread::synchro syn(&outvec, "IntegrationWorker");
      outvec.add(1.0, tmpall);
    }
  }
};

/* This is the class which implements the integration. The class is
   templated by the iterator type. We declare a method |begin| returning
   an iterator to the beginnning of the |ti|-th portion out of total |tn|
   portions for a given level.

   In addition, we define a method which saves all the points to a given
   file. Only for debugging purposes. */

template <typename _Tpit>
class QuadratureImpl : public Quadrature
{
  friend class IntegrationWorker<_Tpit>;
public:
  QuadratureImpl(int d) : Quadrature(d)
  {
  }

  /* Just fill a thread group with workes and run it. */
  void
  integrate(VectorFunctionSet &fs, int level, Vector &out) const override
  {
    // todo: out.length()==func.outdim()
    // todo: dim == func.indim()
    out.zeros();
    sthread::detach_thread_group gr;
    for (int ti = 0; ti < fs.getNum(); ti++)
      gr.insert(std::make_unique<IntegrationWorker<_Tpit>>(*this, fs.getFunc(ti),
                                                           level, ti, fs.getNum(), out));
    gr.run();
  }
  void
  integrate(const VectorFunction &func,
            int level, int tn, Vector &out) const override
  {
    VectorFunctionSet fs(func, tn);
    integrate(fs, level, out);
  }

  /* Just for debugging. */
  void
  savePoints(const string &fname, int level) const
  {
    ofstream fd{fname, std::ios::out | std::ios::trunc};
    if (fd.fail())
      {
        // todo: raise
        std::cerr << "Cannot open file " << fname << " for writing." << std::endl;
        exit(EXIT_FAILURE);
      }
    _Tpit beg = begin(0, 1, level);
    _Tpit end = begin(1, 1, level);
    fd << std::setprecision(12);
    for (_Tpit run = beg; run != end; ++run)
      {
        fd << std::setw(16) << run.weight();
        for (int i = 0; i < dimen(); i++)
          fd << '\t' << std::setw(16) << run.point()[i];
        fd << std::endl;
      }
    fd.close();
  }

  _Tpit
  start(int level) const
  {
    return begin(0, 1, level);
  }
  _Tpit
  end(int level) const
  {
    return begin(1, 1, level);
  }
protected:
  virtual _Tpit begin(int ti, int tn, int level) const = 0;
};

/* This is only an interface to a precalculated data in file {\tt
   precalc\_quadrature.dat} which is basically C coded static data. It
   implements |OneDQuadrature|. The data file is supposed to define the
   following data: number of levels, array of number of points at each
   level, an array of weights and array of points. The both latter array
   store data level by level. An offset for a specific level is stored in
   |offsets| integer sequence.

   The implementing subclasses just fill the necessary data from the
   file, the rest is calculated here. */

class OneDPrecalcQuadrature : public OneDQuadrature
{
  int num_levels;
  const int *num_points;
  const double *weights;
  const double *points;
  IntSequence offsets;
public:
  OneDPrecalcQuadrature(int nlevels, const int *npoints,
                        const double *wts, const double *pts)
    : num_levels(nlevels),  num_points(npoints),
      weights(wts), points(pts), offsets(num_levels)
  {
    calcOffsets();
  }
  ~OneDPrecalcQuadrature() override = default;
  int
  numLevels() const override
  {
    return num_levels;
  }
  int
  numPoints(int level) const override
  {
    return num_points[level-1];
  }
  double
  point(int level, int i) const override
  {
    return points[offsets[level-1]+i];
  }
  double
  weight(int level, int i) const override
  {
    return weights[offsets[level-1]+i];
  }
protected:
  void calcOffsets();
};

/* Just precalculated Gauss--Hermite quadrature. This quadrature integrates exactly integrals
   $$\int_{-\infty}^{\infty} x^ke^{-x^2}{\rm d}x$$
   for level $k$.

   Note that if pluging this one-dimensional quadrature to product or
   Smolyak rule in order to integrate a function $f$ through normally
   distributed inputs, one has to wrap $f$ to
   |GaussConverterFunction| and apply the product or Smolyak rule to the
   new function.

   Check {\tt precalc\_quadrature.dat} for available levels. */

class GaussHermite : public OneDPrecalcQuadrature
{
public:
  GaussHermite();
};

/* Just precalculated Gauss--Legendre quadrature. This quadrature integrates exactly integrals
   $$\int_0^1x^k{\rm d}x$$
   for level $k$.

   Check {\tt precalc\_quadrature.dat} for available levels. */

class GaussLegendre : public OneDPrecalcQuadrature
{
public:
  GaussLegendre();
};

#endif
