/*
 * Copyright © 2004 Ondra Kamenik
 * Copyright © 2019-2024 Dynare Team
 *
 * This file is part of Dynare.
 *
 * Dynare is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * Dynare is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with Dynare.  If not, see <https://www.gnu.org/licenses/>.
 */

// Even more general symmetry tensor.

/* Here we define an abstraction for a tensor, which has a general symmetry,
   but the symmetry is not of what is modelled by Symmetry. This kind of tensor
   comes to existence when we evaluate something like:

   [B_y²u³]_α₁α₂β₁β₂β₃ = … + [g_y³]_γ₁γ₂γ₃ [g_yu]^γ₁_α₁β₃ [g_yu]^γ₂_α₂β₁
                             [g_u]^γ₃_β₂ + …

   If the tensors are unfolded, we obtain a tensor

    g_y³·(g_yu ⊗ g_yu ⊗ g_u)

   Obviously, this tensor can have a symmetry not compatible with ordering
   α₁α₂β₁β₂β₃, (in other words, not compatible with symmetry y²u³). In fact,
   the indices are permuted.

   This kind of tensor must be added to [B_y²u³]. Its dimensions are the same
   as of [B_y²u³], but some coordinates are permuted. The addition is the only
   action we need to do with the tensor.

   Another application where this permuted symmetry tensor appears is a slice
   of a fully symmetric tensor. If the symmetric dimension of the tensor is
   partitioned to continuous parts, and we are interested only in data with a
   given symmetry (permuted) of the partitions, then we have the permuted
   symmetry tensor. For instance, if x is partitioned x=(a,b,c,d), and having
   tensor [f_x³], one can define a slice (subtensor) [f_aca]. The data of this
   tensor are permuted of [f_a²c].

   Here we also define the folded version of permuted symmetry tensor. It
   has permuted symmetry and is partially folded. One can imagine it as a
   product of a few dimensions, each of them is folded and having a few
   variables. The underlying variables are permuted. The product of such
   dimensions is described by PerTensorDimens2. The tensor holding the
   underlying data is FPSTensor. */

#ifndef PS_TENSOR_HH
#define PS_TENSOR_HH

#include "equivalence.hh"
#include "gs_tensor.hh"
#include "kron_prod.hh"
#include "permutation.hh"
#include "sparse_tensor.hh"
#include "tensor.hh"

/* Here we declare a class describing dimensions of permuted symmetry tensor.
   It inherits from TensorDimens and adds a permutation which permutes ‘nvmax’.
   It has two constructors, each corresponds to a context where the tensor
   appears.

   The first constructor calculates the permutation from a given equivalence.

   The second constructor corresponds to dimensions of a slice. Let us take
   [f_aca] as an example. First it calculates TensorDimens of [f_a²c], then it
   calculates a permutation corresponding to ordering of “aca” to “a²c”, and
   applies this permutation on the dimensions as the first constructor. The
   constructor takes only stack sizes (lengths of a, b, c, and d), and
   coordinates of picked partitions.

   Note that inherited methods calcUnfoldColumns() and calcFoldColumns() work,
   since number of columns is independent on the permutation, and
   calcFoldColumns() does not use changed ‘nvmax’, it uses ‘nvs’, so it is
   OK. */

class PerTensorDimens : public TensorDimens
{
private:
  static IntSequence
  sortIntSequence(IntSequence s)
  {
    s.sort();
    return s;
  }

protected:
  Permutation per;

public:
  PerTensorDimens(Symmetry s, IntSequence nvars, const Equivalence& e) :
      TensorDimens(std::move(s), std::move(nvars)), per(e)
  {
    per.apply(nvmax);
  }
  PerTensorDimens(TensorDimens td, const Equivalence& e) : TensorDimens(std::move(td)), per(e)
  {
    per.apply(nvmax);
  }
  PerTensorDimens(TensorDimens td, Permutation p) : TensorDimens(std::move(td)), per(std::move(p))
  {
    per.apply(nvmax);
  }
  PerTensorDimens(const IntSequence& ss, const IntSequence& coor) :
      TensorDimens(ss, sortIntSequence(coor)), per(coor)
  {
    per.apply(nvmax);
  }
  [[nodiscard]] bool
  operator==(const PerTensorDimens& td) const
  {
    return TensorDimens::operator==(td) && per == td.per;
  }
  [[nodiscard]] int
  tailIdentity() const
  {
    return per.tailIdentity();
  }
  [[nodiscard]] const Permutation&
  getPer() const
  {
    return per;
  }
};

/* Here we declare the permuted symmetry unfolded tensor. It has
   PerTensorDimens as a member. It inherits from UTensor which requires to
   implement fold() method. There is no folded counterpart, so in our
   implementation we raise unconditional exception, and return some dummy
   object (just to make it compilable without warnings).

   The class has two sorts of constructors corresponding to a context where it
   appears. The first constructs object from a given matrix, and Kronecker
   product. Within the constructor, all the calculations are performed. Also we
   need to define dimensions, these are the same of the resulting matrix (in
   our example [B_y²u³]) but permuted. The permutation is done in
   PerTensorDimens constructor.

   The second type of constructor is slicing. It makes a slice from
   FSSparseTensor. The slice is given by stack sizes, and coordinates of
   picked stacks.

   There are two algorithms for filling a slice of a sparse tensor. The
   first fillFromSparseOne() works well for more dense tensors, the
   second fillFromSparseTwo() is better for very sparse tensors. We
   provide a static method, which decides what of the two algorithms is
   better. */

class UPSTensor : public UTensor
{
  const PerTensorDimens tdims;

public:
  // UPSTensor constructors from Kronecker product
  /* Here we have four constructors making an UPSTensor from a product
     of matrix and Kronecker product. */

  /* Constructs the tensor from equivalence classes of the given equivalence in
     an order given by the equivalence */
  UPSTensor(TensorDimens td, const Equivalence& e, const ConstTwoDMatrix& a,
            const KronProdAll& kp) :
      UTensor(indor::along_col, PerTensorDimens(td, e).getNVX(), a.nrows(), kp.ncols(), td.dimen()),
      tdims(std::move(td), e)
  {
    kp.mult(a, *this);
  }

  /* Same as the previous one but with optimized KronProdAllOptim, which has a
     different order of matrices than given by the classes in the equivalence.
     This permutation is projected to the permutation of the UPSTensor. */
  UPSTensor(TensorDimens td, const Equivalence& e, const ConstTwoDMatrix& a,
            const KronProdAllOptim& kp) :
      UTensor(indor::along_col, PerTensorDimens(td, Permutation(e, kp.getPer())).getNVX(),
              a.nrows(), kp.ncols(), td.dimen()),
      tdims(std::move(td), Permutation(e, kp.getPer()))
  {
    kp.mult(a, *this);
  }

  /* Same as the first constructor, but the classes of the equivalence are
     permuted by the given permutation. */
  UPSTensor(TensorDimens td, const Equivalence& e, const Permutation& p, const ConstTwoDMatrix& a,
            const KronProdAll& kp) :
      UTensor(indor::along_col, PerTensorDimens(td, Permutation(e, p)).getNVX(), a.nrows(),
              kp.ncols(), td.dimen()),
      tdims(std::move(td), Permutation(e, p))
  {
    kp.mult(a, *this);
  }

  /* Most general constructor. It allows for a permutation of equivalence
     classes, and for optimized KronProdAllOptim, which permutes the permuted
     equivalence classes. */
  UPSTensor(TensorDimens td, const Equivalence& e, const Permutation& p, const ConstTwoDMatrix& a,
            const KronProdAllOptim& kp) :
      UTensor(indor::along_col,
              PerTensorDimens(td, Permutation(e, Permutation(p, kp.getPer()))).getNVX(), a.nrows(),
              kp.ncols(), td.dimen()),
      tdims(std::move(td), Permutation(e, Permutation(p, kp.getPer())))
  {
    kp.mult(a, *this);
  }

  UPSTensor(const FSSparseTensor& t, const IntSequence& ss, const IntSequence& coor,
            PerTensorDimens ptd);
  UPSTensor(const UPSTensor&) = default;
  UPSTensor(UPSTensor&&) = default;

  void increment(IntSequence& v) const override;
  void decrement(IntSequence& v) const override;
  [[nodiscard]] std::unique_ptr<FTensor> fold() const override;

  [[nodiscard]] int getOffset(const IntSequence& v) const override;
  void addTo(FGSTensor& out) const;
  void addTo(UGSTensor& out) const;

  enum class fill_method
  {
    first,
    second
  };
  static fill_method decideFillMethod(const FSSparseTensor& t);

private:
  [[nodiscard]] int tailIdentitySize() const;
  void fillFromSparseOne(const FSSparseTensor& t, const IntSequence& ss, const IntSequence& coor);
  void fillFromSparseTwo(const FSSparseTensor& t, const IntSequence& ss, const IntSequence& coor);
};

/* Here we define an abstraction for the tensor dimension with the symmetry
   like xuv|uv|xu|y|y|x|x|y. These symmetries come as induces symmetries of
   equivalence and some outer symmetry. Thus the underlying variables are
   permuted. One can imagine the dimensions as an unfolded product of
   dimensions which consist of folded products of variables.

   We inherit from PerTensorDimens since we need the permutation implied by the
   equivalence. The new member are the induced symmetries (symmetries of each
   folded dimensions) and ‘ds’ which are sizes of the dimensions. The number of
   folded dimensions is return by ‘numSyms’.

   The object is constructed from outer tensor dimensions and from
   equivalence with optionally permuted classes. */

class PerTensorDimens2 : public PerTensorDimens
{
  InducedSymmetries syms;
  IntSequence ds;

public:
  PerTensorDimens2(const TensorDimens& td, const Equivalence& e, const Permutation& p) :
      PerTensorDimens(td, Permutation(e, p)), syms(e, p, td.getSym()), ds(syms.size())
  {
    setDimensionSizes();
  }
  PerTensorDimens2(const TensorDimens& td, const Equivalence& e) :
      PerTensorDimens(td, e), syms(e, td.getSym()), ds(syms.size())
  {
    setDimensionSizes();
  }
  [[nodiscard]] int
  numSyms() const
  {
    return static_cast<int>(syms.size());
  }
  [[nodiscard]] const Symmetry&
  getSym(int i) const
  {
    return syms[i];
  }
  [[nodiscard]] int
  calcMaxOffset() const
  {
    return ds.mult();
  }
  [[nodiscard]] int calcOffset(const IntSequence& coor) const;
  void print() const;

protected:
  void setDimensionSizes();
};

/* Here we define an abstraction of the permuted symmetry folded tensor. It is
   needed in context of the Faà Di Bruno formula for folded stack container
   multiplied with container of dense folded tensors, or multiplied by one full
   symmetry sparse tensor.

   For example, if we perform the Faà Di Bruno for $F=f(z)$, where
      ⎛g(x,y,u,v)⎞
      ⎢ h(x,y,u) ⎥
    z=⎢     x    ⎥
      ⎝     y    ⎠
   we get for one concrete equivalence:

    [F_x⁴y³u³v²] = … + [f_g²h²x²y]·([g]_xv ⊗ [g]_u²v ⊗ [h]_xu ⊗ [h]_y² ⊗
                                    [I]_x ⊗ [I]_x ⊗ [I]_y)
                     + …

   The class FPSTensor represents the tensor on the right. Its dimension
   corresponds to a product of 7 dimensions with the following symmetries:
   xv|u²v|xu|y²|x|x|y. Such a dimension is described by PerTensorDimens2.

   The tensor is constructed in a context of stack container multiplication,
   so, it is constructed from dimensions ‘td’ (dimensions of the output
   tensor), stack product ‘sp’ (implied symmetries picking tensors from a stack
   container, here it is z), then a sorted integer sequence of the picked
   stacks of the stack product (it is always sorted, here it is
   (0,0,1,1,2,2,3)), then the tensor [f_g²h²x²y] (its symmetry must be the same
   as symmetry given by the ‘istacks’), and finally from the equivalence with
   permuted classes.

   We implement increment() and getOffset() methods, decrement() and unfold()
   raise an exception. Also, we implement addTo() method, which adds the tensor
   data (partially unfolded) to folded general symmetry tensor. */

template<typename _Ttype>
class StackProduct;

class FPSTensor : public FTensor
{
  const PerTensorDimens2 tdims;

public:
  /* As for UPSTensor, we provide four constructors allowing for
     combinations of permuting equivalence classes, and optimization of
     KronProdAllOptim. These constructors multiply with dense general
     symmetry tensor (coming from the dense container, or as a dense slice
     of the full symmetry sparse tensor). In addition to these 4
     constructors, we have one constructor multiplying with general
     symmetry sparse tensor (coming as a sparse slice of the full symmetry
     sparse tensor). */
  FPSTensor(const TensorDimens& td, const Equivalence& e, const ConstTwoDMatrix& a,
            const KronProdAll& kp) :
      FTensor(indor::along_col, PerTensorDimens(td, e).getNVX(), a.nrows(), kp.ncols(), td.dimen()),
      tdims(td, e)
  {
    kp.mult(a, *this);
  }
  FPSTensor(const TensorDimens& td, const Equivalence& e, const ConstTwoDMatrix& a,
            const KronProdAllOptim& kp) :
      FTensor(indor::along_col, PerTensorDimens(td, Permutation(e, kp.getPer())).getNVX(),
              a.nrows(), kp.ncols(), td.dimen()),
      tdims(td, e, kp.getPer())
  {
    kp.mult(a, *this);
  }
  FPSTensor(const TensorDimens& td, const Equivalence& e, const Permutation& p,
            const ConstTwoDMatrix& a, const KronProdAll& kp) :
      FTensor(indor::along_col, PerTensorDimens(td, Permutation(e, p)).getNVX(), a.nrows(),
              kp.ncols(), td.dimen()),
      tdims(td, e, p)
  {
    kp.mult(a, *this);
  }
  FPSTensor(const TensorDimens& td, const Equivalence& e, const Permutation& p,
            const ConstTwoDMatrix& a, const KronProdAllOptim& kp) :
      FTensor(indor::along_col,
              PerTensorDimens(td, Permutation(e, Permutation(p, kp.getPer()))).getNVX(), a.nrows(),
              kp.ncols(), td.dimen()),
      tdims(td, e, Permutation(p, kp.getPer()))
  {
    kp.mult(a, *this);
  }

  FPSTensor(const TensorDimens& td, const Equivalence& e, const Permutation& p,
            const GSSparseTensor& t, const KronProdAll& kp);

  FPSTensor(const FPSTensor&) = default;
  FPSTensor(FPSTensor&&) = default;

  void increment(IntSequence& v) const override;
  void decrement(IntSequence& v) const override;
  [[nodiscard]] std::unique_ptr<UTensor> unfold() const override;

  [[nodiscard]] int getOffset(const IntSequence& v) const override;
  void addTo(FGSTensor& out) const;
};

#endif
