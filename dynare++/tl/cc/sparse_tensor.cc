// Copyright 2004, Ondra Kamenik

#include "sparse_tensor.hh"
#include "fs_tensor.hh"
#include "tl_exception.hh"

#include <cmath>

/* This is straightforward. Before we insert anything, we do a few
   checks. Then we reset |first_nz_row| and |last_nz_row| if necessary. */

void
SparseTensor::insert(const IntSequence &key, int r, double c)
{
  TL_RAISE_IF(r < 0 || r >= nr,
              "Row number out of dimension of tensor in SparseTensor::insert");
  TL_RAISE_IF(key.size() != dimen(),
              "Wrong length of key in SparseTensor::insert");
  TL_RAISE_IF(!std::isfinite(c),
              "Insertion of non-finite value in SparseTensor::insert");

  iterator first_pos = m.lower_bound(key);

  // check that pair |key| and |r| is unique
  iterator last_pos = m.upper_bound(key);
  for (iterator it = first_pos; it != last_pos; ++it)
    if ((*it).second.first == r)
      {
        TL_RAISE("Duplicate <key, r> insertion in SparseTensor::insert");
        return;
      }

  m.insert(first_pos, Map::value_type(key, Item(r, c)));
  if (first_nz_row > r)
    first_nz_row = r;
  if (last_nz_row < r)
    last_nz_row = r;
}

/* This returns true if all items are finite (not Nan nor Inf). */

bool
SparseTensor::isFinite() const
{
  bool res = true;
  const_iterator run = m.begin();
  while (res && run != m.end())
    {
      if (!std::isfinite((*run).second.second))
        res = false;
      ++run;
    }
  return res;
}

/* This returns a ratio of a number of non-zero columns in folded
   tensor to the total number of columns. */

double
SparseTensor::getFoldIndexFillFactor() const
{
  int cnt = 0;
  const_iterator start_col = m.begin();
  while (start_col != m.end())
    {
      cnt++;
      const IntSequence &key = (*start_col).first;
      start_col = m.upper_bound(key);
    }

  return ((double) cnt)/ncols();
}

/* This returns a ratio of a number of non-zero columns in unfolded
   tensor to the total number of columns. */

double
SparseTensor::getUnfoldIndexFillFactor() const
{
  int cnt = 0;
  const_iterator start_col = m.begin();
  while (start_col != m.end())
    {
      const IntSequence &key = (*start_col).first;
      Symmetry s(key);
      cnt += Tensor::noverseq(s);
      start_col = m.upper_bound(key);
    }

  return ((double) cnt)/ncols();
}

/* This prints the fill factor and all items. */

void
SparseTensor::print() const
{
  printf("Fill: %3.2f %%\n", 100*getFillFactor());
  const_iterator start_col = m.begin();
  while (start_col != m.end())
    {
      const IntSequence &key = (*start_col).first;
      printf("Column: "); key.print();
      const_iterator end_col = m.upper_bound(key);
      int cnt = 1;
      for (const_iterator run = start_col; run != end_col; ++run, cnt++)
        {
          if ((cnt/7)*7 == cnt)
            printf("\n");
          printf("%d(%6.2g)  ", (*run).second.first, (*run).second.second);
        }
      printf("\n");
      start_col = end_col;
    }
}

FSSparseTensor::FSSparseTensor(int d, int nvar, int r)
  : SparseTensor(d, r, FFSTensor::calcMaxOffset(nvar, d)),
    nv(nvar), sym(d)
{
}

FSSparseTensor::FSSparseTensor(const FSSparseTensor &t)
  : SparseTensor(t),
    nv(t.nvar()), sym(t.sym)
{
}

void
FSSparseTensor::insert(const IntSequence &key, int r, double c)
{
  TL_RAISE_IF(!key.isSorted(),
              "Key is not sorted in FSSparseTensor::insert");
  TL_RAISE_IF(key[key.size()-1] >= nv || key[0] < 0,
              "Wrong value of the key in FSSparseTensor::insert");
  SparseTensor::insert(key, r, c);
}

/* We go through the tensor |t| which is supposed to have single
   column. If the item of |t| is nonzero, we make a key by sorting the
   index, and then we go through all items having the same key (it is its
   column), obtain the row number and the element, and do the
   multiplication.

   The test for non-zero is |a != 0.0|, since there will be items which
   are exact zeros.

   I have also tried to make the loop through the sparse tensor outer, and
   find index of tensor |t| within the loop. Surprisingly, it is little
   slower (for monomial tests with probability of zeros equal 0.3). But
   everything depends how filled is the sparse tensor. */

void
FSSparseTensor::multColumnAndAdd(const Tensor &t, Vector &v) const
{
  // check compatibility of input parameters
  TL_RAISE_IF(v.length() != nrows(),
              "Wrong size of output vector in FSSparseTensor::multColumnAndAdd");
  TL_RAISE_IF(t.dimen() != dimen(),
              "Wrong dimension of tensor in FSSparseTensor::multColumnAndAdd");
  TL_RAISE_IF(t.ncols() != 1,
              "The input tensor is not single-column in FSSparseTensor::multColumnAndAdd");

  for (Tensor::index it = t.begin(); it != t.end(); ++it)
    {
      int ind = *it;
      double a = t.get(ind, 0);
      if (a != 0.0)
        {
          IntSequence key(it.getCoor());
          key.sort();

          // check that |key| is within the range
          TL_RAISE_IF(key[0] < 0 || key[key.size()-1] >= nv,
                      "Wrong coordinates of index in FSSparseTensor::multColumnAndAdd");

          const_iterator first_pos = m.lower_bound(key);
          const_iterator last_pos = m.upper_bound(key);
          for (const_iterator cit = first_pos; cit != last_pos; ++cit)
            {
              int r = (*cit).second.first;
              double c = (*cit).second.second;
              v[r] += c * a;
            }
        }
    }
}

void
FSSparseTensor::print() const
{
  printf("FS Sparse tensor: dim=%d, nv=%d, (%dx%d)\n", dim, nv, nr, nc);
  SparseTensor::print();
}

// |GSSparseTensor| slicing constructor
/* This is the same as |@<|FGSTensor| slicing from |FSSparseTensor|@>|. */
GSSparseTensor::GSSparseTensor(const FSSparseTensor &t, const IntSequence &ss,
                               const IntSequence &coor, const TensorDimens &td)
  : SparseTensor(td.dimen(), t.nrows(), td.calcFoldMaxOffset()),
    tdims(td)
{
  // set |lb| and |ub| to lower and upper bounds of slice indices
  /* This is the same as |@<set |lb| and |ub| to lower and upper bounds
     of indices@>| in {\tt gs\_tensor.cpp}, see that file for details. */
  IntSequence s_offsets(ss.size(), 0);
  for (int i = 1; i < ss.size(); i++)
    s_offsets[i] = s_offsets[i-1] + ss[i-1];

  IntSequence lb(coor.size());
  IntSequence ub(coor.size());
  for (int i = 0; i < coor.size(); i++)
    {
      lb[i] = s_offsets[coor[i]];
      ub[i] = s_offsets[coor[i]] + ss[coor[i]] - 1;
    }

  FSSparseTensor::const_iterator lbi = t.getMap().lower_bound(lb);
  FSSparseTensor::const_iterator ubi = t.getMap().upper_bound(ub);
  for (FSSparseTensor::const_iterator run = lbi; run != ubi; ++run)
    {
      if (lb.lessEq((*run).first) && (*run).first.lessEq(ub))
        {
          IntSequence c((*run).first);
          c.add(-1, lb);
          insert(c, (*run).second.first, (*run).second.second);
        }
    }

}

void
GSSparseTensor::insert(const IntSequence &s, int r, double c)
{
  TL_RAISE_IF(!s.less(tdims.getNVX()),
              "Wrong coordinates of index in GSSparseTensor::insert");
  SparseTensor::insert(s, r, c);
}

void
GSSparseTensor::print() const
{
  printf("GS Sparse tensor: (%dx%d)\nSymmetry: ", nr, nc);
  tdims.getSym().print();
  printf("NVS: ");
  tdims.getNVS().print();
  SparseTensor::print();
}
