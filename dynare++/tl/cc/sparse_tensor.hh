// Copyright 2004, Ondra Kamenik

// Sparse tensor.

/* Here we declare a sparse full and general symmetry tensors with the
   multidimensional index along columns. We implement them as a |multimap|
   associating to each sequence of coordinates |IntSequence| a set of
   pairs (row, number). This is very convenient but not optimal in terms
   of memory consumption. So the implementation can be changed.

   The current |multimap| implementation allows insertions.  Another
   advantage of this approach is that we do not need to calculate column
   numbers from the |IntSequence|, since the column is accessed directly
   via the key which is |IntSequence|.

   The only operation we need to do with the full symmetry sparse tensor
   is a left multiplication of a row oriented single column tensor. The
   result of such operation is a column of the same size as the sparse
   tensor. Other important operations are slicing operations. We need to
   do sparse and dense slices of full symmetry sparse tensors. In fact,
   the only constructor of general symmetry sparse tensor is slicing from
   the full symmetry sparse. */

#ifndef SPARSE_TENSOR_H
#define SPARSE_TENSOR_H

#include "symmetry.hh"
#include "tensor.hh"
#include "gs_tensor.hh"
#include "Vector.hh"

#include <map>

using namespace std;

// |ltseq| predicate
struct ltseq
{
  bool
  operator()(const IntSequence &s1, const IntSequence &s2) const
  {
    return s1 < s2;
  }
};

/* This is a super class of both full symmetry and general symmetry
   sparse tensors. It contains a |multimap| and implements insertions. It
   tracks maximum and minimum row, for which there is an item. */

class SparseTensor
{
public:
  typedef pair<int, double> Item;
  typedef multimap<IntSequence, Item, ltseq> Map;
  typedef Map::const_iterator const_iterator;
protected:
  typedef Map::iterator iterator;

  Map m;
  const int dim;
  const int nr;
  const int nc;
  int first_nz_row;
  int last_nz_row;
public:
  SparseTensor(int d, int nnr, int nnc)
    : dim(d), nr(nnr), nc(nnc), first_nz_row(nr), last_nz_row(-1)
  {
  }
  SparseTensor(const SparseTensor &t)
    : m(t.m), dim(t.dim), nr(t.nr), nc(t.nc)
  {
  }
  virtual ~SparseTensor()
  = default;
  void insert(const IntSequence &s, int r, double c);
  const Map &
  getMap() const
  {
    return m;
  }
  int
  dimen() const
  {
    return dim;
  }
  int
  nrows() const
  {
    return nr;
  }
  int
  ncols() const
  {
    return nc;
  }
  double
  getFillFactor() const
  {
    return ((double) m.size())/(nrows()*ncols());
  }
  double getFoldIndexFillFactor() const;
  double getUnfoldIndexFillFactor() const;
  int
  getNumNonZero() const
  {
    return m.size();
  }
  int
  getFirstNonZeroRow() const
  {
    return first_nz_row;
  }
  int
  getLastNonZeroRow() const
  {
    return last_nz_row;
  }
  virtual const Symmetry&getSym() const = 0;
  void print() const;
  bool isFinite() const;
};

/* This is a full symmetry sparse tensor. It implements
   |multColumnAndAdd| and in addition to |sparseTensor|, it has |nv|
   (number of variables), and symmetry (basically it is a dimension). */

class FSSparseTensor : public SparseTensor
{
public:
  typedef SparseTensor::const_iterator const_iterator;
private:
  const int nv;
  const Symmetry sym;
public:
  FSSparseTensor(int d, int nvar, int r);
  FSSparseTensor(const FSSparseTensor &t);
  void insert(const IntSequence &s, int r, double c);
  void multColumnAndAdd(const Tensor &t, Vector &v) const;
  const Symmetry &
  getSym() const
  {
    return sym;
  }
  int
  nvar() const
  {
    return nv;
  }
  void print() const;
};

/* This is a general symmetry sparse tensor. It has |TensorDimens| and
   can be constructed as a slice of the full symmetry sparse tensor. The
   slicing constructor takes the same form as the slicing |FGSTensor|
   constructor from full symmetry sparse tensor. */

class GSSparseTensor : public SparseTensor
{
public:
  typedef SparseTensor::const_iterator const_iterator;
private:
  const TensorDimens tdims;
public:
  GSSparseTensor(const FSSparseTensor &t, const IntSequence &ss,
                 const IntSequence &coor, const TensorDimens &td);
  GSSparseTensor(const GSSparseTensor &t)
     
  = default;
  void insert(const IntSequence &s, int r, double c);
  const Symmetry &
  getSym() const
  {
    return tdims.getSym();
  }
  const TensorDimens &
  getDims() const
  {
    return tdims;
  }
  void print() const;

};

#endif
