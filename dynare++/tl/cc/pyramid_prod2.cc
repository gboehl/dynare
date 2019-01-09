// Copyright 2004, Ondra Kamenik

#include "pyramid_prod2.hh"
#include "rfs_tensor.hh"

/* Here we only call |sp.createPackedColumns(c, cols, unit_flag)| which
   fills |cols| and |unit_flag| for the given column |c|. Then we set
   |end_seq| according to |unit_flag| and columns lengths. */

IrregTensorHeader::IrregTensorHeader(const StackProduct<FGSTensor> &sp,
                                     const IntSequence &c)
  : nv(sp.getAllSize()),
    unit_flag(sp.dimen()),
    cols(new Vector *[sp.dimen()]),
    end_seq(sp.dimen())
{
  sp.createPackedColumns(c, cols, unit_flag);
  for (int i = 0; i < sp.dimen(); i++)
    {
      end_seq[i] = cols[i]->length();
      if (unit_flag[i] != -1)
        end_seq[i] = unit_flag[i]+1;
    }
}

/* Here we have to increment the given integer sequence. We do it by
   the following code, whose pattern is valid for all tensor. The only
   difference is how we increment item of coordinates. */

void
IrregTensorHeader::increment(IntSequence &v) const
{
  TL_RAISE_IF(v.size() != dimen(),
              "Wrong size of coordinates in IrregTensorHeader::increment");

  if (v.size() == 0)
    return;
  int i = v.size()-1;

  // increment |i|-th item in coordinate |v|
  /* Here we increment item of coordinates. Whenever we reached end of
     column coming from matrices, and |unit_flag| is not $-1$, we have to
     jump to that |unit_flag|. */
  v[i]++;
  if (unit_flag[i] != -1 && v[i] == cols[i]->length()-1)
    v[i] = unit_flag[i];

  while (i > 0 && v[i] == end_seq[i])
    {
      v[i] = 0;
      i--;
      // increment |i|-th item in coordinate |v|
      /* Same code as above */
      v[i]++;
      if (unit_flag[i] != -1 && v[i] == cols[i]->length()-1)
        v[i] = unit_flag[i];
    }
}

IrregTensorHeader::~IrregTensorHeader()
{
  for (int i = 0; i < dimen(); i++)
    delete cols[i];
  delete [] cols;
}

/* It is a product of all column lengths. */

int
IrregTensorHeader::calcMaxOffset() const
{
  int res = 1;
  for (int i = 0; i < dimen(); i++)
    res *= cols[i]->length();
  return res;
}

/* Everything is done in |IrregTensorHeader|, only we have to Kronecker
   multiply all columns of the header. */

IrregTensor::IrregTensor(const IrregTensorHeader &h)
  : Tensor(along_row, IntSequence(h.dimen(), 0), h.end_seq,
           h.calcMaxOffset(), 1, h.dimen()),
    header(h)
{
  if (header.dimen() == 1)
    {
      getData() = *(header.cols[0]);
      return;
    }

  auto *last = new Vector(*(header.cols[header.dimen()-1]));
  for (int i = header.dimen()-2; i > 0; i--)
    {
      auto *newlast = new Vector(last->length()*header.cols[i]->length());
      KronProd::kronMult(ConstVector(*(header.cols[i])),
                         ConstVector(*last), *newlast);
      delete last;
      last = newlast;
    }
  KronProd::kronMult(ConstVector(*(header.cols[0])),
                     ConstVector(*last), getData());
  delete last;
}

void
IrregTensor::addTo(FRSingleTensor &out) const
{
  for (index it = begin(); it != end(); ++it)
    {
      IntSequence tmp(it.getCoor());
      tmp.sort();
      Tensor::index ind(&out, tmp);
      out.get(*ind, 0) += get(*it, 0);
    }
}
