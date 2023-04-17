/*
 * Copyright © 2004 Ondra Kamenik
 * Copyright © 2019-2023 Dynare Team
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

// Stack of containers.

/* Here we develop abstractions for stacked containers of tensors. For
   instance, in perturbation methods for DSGE models, we need the function:

                   ⎛ G(y*,u,u′,σ) ⎞
    z(y*,u,u′,σ) = ⎢  g(y*,u,σ)   ⎥
                   ⎢      y*      ⎥
                   ⎝      u       ⎠

   and we need to calculate one step of Faà Di Bruno formula:

                                         ₗ
    [B_sᵏ]_α₁…αₗ = [f_zˡ]_β₁…βₗ    ∑     ∏  [z_{s^|cₘ|}]_cₘ(α)^βₘ
                                c∈ℳₗ,ₖ ᵐ⁼¹

   where we have containers for derivatives of G and g.

   The main purpose of this file is to define abstractions for stack of
   containers and possibly raw variables, and code multAndAdd() method
   calculating (one step of) the Faà Di Bruno formula for folded and
   unfolded tensors. Note also, that tensors [f_zˡ] are sparse.

   The abstractions are built as follows. At the top, there is an
   interface describing stack of columns. It contains pure virtual
   methods needed for manipulating the container stack. For technical
   reasons it is a template. Both versions (folded, and unfolded) provide
   all interface necessary for implementation of multAndAdd(). The second
   way of inheritance is first general implementation of the interface
   StackContainer, and then specific (ZContainer for our specific z).
   The only method which is virtual also after StackContainer is
   getType(), which is implemented in the specialization and determines
   behaviour of the stack. The complete classes are obtained by
   inheriting from the both branches, as it is drawn below:


                        StackContainerInterface<FGSTensor>
                          ↙                           ↘
            StackContainer<FGSTensor>             FoldedStackContainer
                        ↓
              ZContainer<FGSTensor>                 ↙
                                  ↘
                                   FoldedZContainer


                        StackContainerInterface<UGSTensor>
                          ↙                           ↘
            StackContainer<UGSTensor>           UnfoldedStackContainer
                        ↓
              ZContainer<UGSTensor>                 ↙
                                 ↘
                                  UnfoldedZContainer


   We have also two supporting classes StackProduct and KronProdStack
   and a number of worker classes used as threads. */

#ifndef STACK_CONTAINER_H
#define STACK_CONTAINER_H

#include "int_sequence.hh"
#include "equivalence.hh"
#include "tl_static.hh"
#include "t_container.hh"
#include "kron_prod.hh"
#include "permutation.hh"
#include "sthread.hh"

/* Here is the general interface to stack container. The subclasses
   maintain IntSequence of stack sizes, i.e. size of G, g, y, and
   u. Then a convenience IntSequence of stack offsets. Then vector of
   pointers to containers, in our example G, and g.

   A non-virtual subclass must implement getType() which determines
   dependency of stack items on symmetries. There are three possible types
   for a symmetry. Either the stack item derivative wrt. the symmetry is
   a matrix, or a unit matrix, or zero.

   Method isZero() returns true if the derivative of a given stack item
   wrt. to given symmetry is zero as defined by getType() or the
   derivative is not present in the container. In this way, we can
   implement the formula conditional some of the tensors are zero, which
   is not true (they are only missing).

   Method createPackedColumn() returns a vector of stack derivatives with
   respect to the given symmetry and of the given column, where all zeros
   from zero types, or unit matrices are deleted. See kron_prod.hh for an
   explanation. */

template<class _Ttype>
class StackContainerInterface
{
public:
  using _Ctype = TensorContainer<_Ttype>;
  enum class itype { matrix, unit, zero };
public:
  StackContainerInterface() = default;
  virtual ~StackContainerInterface() = default;
  virtual const IntSequence &getStackSizes() const = 0;
  virtual IntSequence &getStackSizes() = 0;
  virtual const IntSequence &getStackOffsets() const = 0;
  virtual IntSequence &getStackOffsets() = 0;
  virtual int numConts() const = 0;
  virtual const _Ctype &getCont(int i) const = 0;
  virtual itype getType(int i, const Symmetry &s) const = 0;
  virtual int numStacks() const = 0;
  virtual bool isZero(int i, const Symmetry &s) const = 0;
  virtual const _Ttype &getMatrix(int i, const Symmetry &s) const = 0;
  virtual int getLengthOfMatrixStacks(const Symmetry &s) const = 0;
  virtual int getUnitPos(const Symmetry &s) const = 0;
  virtual std::unique_ptr<Vector> createPackedColumn(const Symmetry &s,
                                                     const IntSequence &coor,
                                                     int &iu) const = 0;
  int
  getAllSize() const
  {
    return getStackOffsets()[numStacks()-1]
      + getStackSizes()[numStacks()-1];
  }
};

/* Here is StackContainer, which implements almost all interface
   StackContainerInterface but one method getType() which is left for
   implementation to specializations.

   It does not own its tensors. */

template<class _Ttype>
class StackContainer : virtual public StackContainerInterface<_Ttype>
{
public:
  using _Stype = StackContainerInterface<_Ttype>;
  using _Ctype = typename StackContainerInterface<_Ttype>::_Ctype;
  using itype = typename StackContainerInterface<_Ttype>::itype;
protected:
  int num_conts;
  IntSequence stack_sizes;
  IntSequence stack_offsets;
  std::vector<const _Ctype *> conts;
public:
  StackContainer(int ns, int nc)
    : stack_sizes(ns, 0), stack_offsets(ns, 0),
      conts(nc)
  {
  }
  const IntSequence &
  getStackSizes() const override
  {
    return stack_sizes;
  }
  IntSequence &
  getStackSizes() override
  {
    return stack_sizes;
  }
  const IntSequence &
  getStackOffsets() const override
  {
    return stack_offsets;
  }
  IntSequence &
  getStackOffsets() override
  {
    return stack_offsets;
  }
  int
  numConts() const override
  {
    return conts.size();
  }
  const _Ctype &
  getCont(int i) const override
  {
    return *(conts[i]);
  }
  itype getType(int i, const Symmetry &s) const override = 0;
  int
  numStacks() const override
  {
    return stack_sizes.size();
  }
  bool
  isZero(int i, const Symmetry &s) const override
  {
    TL_RAISE_IF(i < 0 || i >= numStacks(),
                "Wrong index to stack in StackContainer::isZero.");
    return (getType(i, s) == itype::zero
            || (getType(i, s) == itype::matrix && !conts[i]->check(s)));
  }

  const _Ttype &
  getMatrix(int i, const Symmetry &s) const override
  {
    TL_RAISE_IF(isZero(i, s) || getType(i, s) == itype::unit,
                "Matrix is not returned in StackContainer::getMatrix");
    return conts[i]->get(s);
  }

  int
  getLengthOfMatrixStacks(const Symmetry &s) const override
  {
    int res = 0;
    int i = 0;
    while (i < numStacks() && getType(i, s) == itype::matrix)
      res += stack_sizes[i++];
    return res;
  }

  int
  getUnitPos(const Symmetry &s) const override
  {
    if (s.dimen() != 1)
      return -1;
    int i = numStacks()-1;
    while (i >= 0 && getType(i, s) != itype::unit)
      i--;
    return i;
  }

  std::unique_ptr<Vector>
  createPackedColumn(const Symmetry &s,
                     const IntSequence &coor, int &iu) const override
  {
    TL_RAISE_IF(s.dimen() != coor.size(),
                "Incompatible coordinates for symmetry in StackContainer::createPackedColumn");

    int len = getLengthOfMatrixStacks(s);
    iu = -1;
    int i = 0;
    if (-1 != (i = getUnitPos(s)))
      {
        iu = stack_offsets[i] + coor[0];
        len++;
      }

    auto res = std::make_unique<Vector>(len);
    i = 0;
    while (i < numStacks() && getType(i, s) == itype::matrix)
      {
        const _Ttype &t = getMatrix(i, s);
        Tensor::index ind(t, coor);
        Vector subres(*res, stack_offsets[i], stack_sizes[i]);
        subres = ConstGeneralMatrix(t).getCol(*ind);
        i++;
      }
    if (iu != -1)
      (*res)[len-1] = 1;

    return res;
  }

protected:
  void
  calculateOffsets()
  {
    stack_offsets[0] = 0;
    for (int i = 1; i < stack_offsets.size(); i++)
      stack_offsets[i] = stack_offsets[i-1] + stack_sizes[i-1];
  }
};

class WorkerFoldMAADense;
class WorkerFoldMAASparse1;
class WorkerFoldMAASparse2;
class WorkerFoldMAASparse4;
class FoldedStackContainer : virtual public StackContainerInterface<FGSTensor>
{
  friend class WorkerFoldMAADense;
  friend class WorkerFoldMAASparse1;
  friend class WorkerFoldMAASparse2;
  friend class WorkerFoldMAASparse4;
public:
  static constexpr double fill_threshold = 0.00005;
  void
  multAndAdd(int dim, const TensorContainer<FSSparseTensor> &c,
             FGSTensor &out) const
  {
    if (c.check(Symmetry{dim}))
      multAndAdd(c.get(Symmetry{dim}), out);
  }
  void multAndAdd(const FSSparseTensor &t, FGSTensor &out) const;
  void multAndAdd(int dim, const FGSContainer &c, FGSTensor &out) const;
protected:
  void multAndAddSparse1(const FSSparseTensor &t, FGSTensor &out) const;
  void multAndAddSparse2(const FSSparseTensor &t, FGSTensor &out) const;
  void multAndAddSparse3(const FSSparseTensor &t, FGSTensor &out) const;
  void multAndAddSparse4(const FSSparseTensor &t, FGSTensor &out) const;
  void multAndAddStacks(const IntSequence &fi, const FGSTensor &g,
                        FGSTensor &out, std::mutex &mut) const;
  void multAndAddStacks(const IntSequence &fi, const GSSparseTensor &g,
                        FGSTensor &out, std::mutex &mut) const;
};

class WorkerUnfoldMAADense;
class WorkerUnfoldMAASparse1;
class WorkerUnfoldMAASparse2;
class UnfoldedStackContainer : virtual public StackContainerInterface<UGSTensor>
{
  friend class WorkerUnfoldMAADense;
  friend class WorkerUnfoldMAASparse1;
  friend class WorkerUnfoldMAASparse2;
public:
  static constexpr double fill_threshold = 0.00005;
  void
  multAndAdd(int dim, const TensorContainer<FSSparseTensor> &c,
             UGSTensor &out) const
  {
    if (c.check(Symmetry{dim}))
      multAndAdd(c.get(Symmetry{dim}), out);
  }
  void multAndAdd(const FSSparseTensor &t, UGSTensor &out) const;
  void multAndAdd(int dim, const UGSContainer &c, UGSTensor &out) const;
protected:
  void multAndAddSparse1(const FSSparseTensor &t, UGSTensor &out) const;
  void multAndAddSparse2(const FSSparseTensor &t, UGSTensor &out) const;
  void multAndAddStacks(const IntSequence &fi, const UGSTensor &g,
                        UGSTensor &out, std::mutex &mut) const;
};

/* Here is the specialization of the StackContainer. We implement
   here the x needed in DSGE context for welfare assessment. We implement getType() and define a constructor feeding the data and sizes.

   It depends on four variables U(y,u,u',σ), the variable u' being introduced to enable additions with 4-variable tensors*/

template<class _Ttype>
class XContainer : public StackContainer<_Ttype>
{
public:
  using _Tparent = StackContainer<_Ttype>;
  using _Stype = StackContainerInterface<_Ttype>;
  using _Ctype = typename _Tparent::_Ctype;
  using itype = typename _Tparent::itype;
  XContainer(const _Ctype *g, int ng)
    : _Tparent(1, 1)
  {
    _Tparent::stack_sizes = { ng };
    _Tparent::conts[0] = g;
    _Tparent::calculateOffsets();
  }

  /* Here we say, what happens if we derive z. recall the top of the
     file, how z looks, and code is clear. */

  itype
  getType(int i, const Symmetry &s) const override
  {
    if (i==0)
      {
        if (s[2] > 0)
          return itype::zero;
        else
          return itype::matrix;
      }

    TL_RAISE("Wrong stack index in XContainer::getType");
  }

};

class FoldedXContainer : public XContainer<FGSTensor>,
                         public FoldedStackContainer
{
public:
  using _Ctype = TensorContainer<FGSTensor>;
  FoldedXContainer(const _Ctype *g, int ng)
    : XContainer<FGSTensor>(g, ng)
  {
  }
};

class UnfoldedXContainer : public XContainer<UGSTensor>,
                           public UnfoldedStackContainer
{
public:
  using _Ctype = TensorContainer<UGSTensor>;
  UnfoldedXContainer(const _Ctype *g, int ng)
    : XContainer<UGSTensor>(g, ng)
  {
  }
};
/* Here is the specialization of the StackContainer. We implement
   here the z needed in DSGE context. We implement getType() and define
   a constructor feeding the data and sizes.

   Note that it has two containers, the first is dependent on four variables
   G(y*,u,u′,σ), and the second dependent on three variables g(y*,u,σ). So that
   we would be able to stack them, we make the second container g be dependent
   on four variables, the third being u′ a dummy and always returning zero if
   dimension of u′ is positive. */

template<class _Ttype>
class ZContainer : public StackContainer<_Ttype>
{
public:
  using _Tparent = StackContainer<_Ttype>;
  using _Stype = StackContainerInterface<_Ttype>;
  using _Ctype = typename _Tparent::_Ctype;
  using itype = typename _Tparent::itype;
  ZContainer(const _Ctype *gss, int ngss, const _Ctype *g, int ng,
             int ny, int nu)
    : _Tparent(4, 2)
  {
    _Tparent::stack_sizes = { ngss, ng, ny, nu };
    _Tparent::conts[0] = gss;
    _Tparent::conts[1] = g;
    _Tparent::calculateOffsets();
  }

  /* Here we say, what happens if we derive z. recall the top of the
     file, how z looks, and code is clear. */

  itype
  getType(int i, const Symmetry &s) const override
  {
    if (i == 0)
      return itype::matrix;
    if (i == 1)
      {
        if (s[2] > 0)
          return itype::zero;
        else
          return itype::matrix;
      }
    if (i == 2)
      {
        if (s == Symmetry{1, 0, 0, 0})
          return itype::unit;
        else
          return itype::zero;
      }
    if (i == 3)
      {
        if (s == Symmetry{0, 1, 0, 0})
          return itype::unit;
        else
          return itype::zero;
      }

    TL_RAISE("Wrong stack index in ZContainer::getType");
  }

};

class FoldedZContainer : public ZContainer<FGSTensor>,
                         public FoldedStackContainer
{
public:
  using _Ctype = TensorContainer<FGSTensor>;
  FoldedZContainer(const _Ctype *gss, int ngss, const _Ctype *g, int ng,
                   int ny, int nu)
    : ZContainer<FGSTensor>(gss, ngss, g, ng, ny, nu)
  {
  }
};

class UnfoldedZContainer : public ZContainer<UGSTensor>,
                           public UnfoldedStackContainer
{
public:
  using _Ctype = TensorContainer<UGSTensor>;
  UnfoldedZContainer(const _Ctype *gss, int ngss, const _Ctype *g, int ng,
                     int ny, int nu)
    : ZContainer<UGSTensor>(gss, ngss, g, ng, ny, nu)
  {
  }
};

/* Here we have another specialization of container used in context of
   DSGE. We define a container for

    G(y,u,u′,σ)=g**(g*(y,u,σ),u′,σ)

   For some reason, the symmetry of g** has length 4 although it
   is really dependent on three variables (To now the reason, consult
   the ZContainer class declaration). So, it has four stacks, the
   third one is dummy, and always returns zero. The first stack
   corresponds to a container of g*. */

template<class _Ttype>
class GContainer : public StackContainer<_Ttype>
{
public:
  using _Tparent = StackContainer<_Ttype>;
  using _Stype = StackContainerInterface<_Ttype>;
  using _Ctype = typename StackContainer<_Ttype>::_Ctype;
  using itype = typename StackContainer<_Ttype>::itype;
  GContainer(const _Ctype *gs, int ngs, int nu)
    : StackContainer<_Ttype>(4, 1)
  {
    _Tparent::stack_sizes = { ngs, nu, nu, 1 };
    _Tparent::conts[0] = gs;
    _Tparent::calculateOffsets();
  }

  /* Here we define the dependencies in g**(g*(y,u,σ),u′,σ). Also note, that
     first derivative of g* w.r.t. σ is always zero, so we also add this
     information. */

  itype
  getType(int i, const Symmetry &s) const override
  {
    if (i == 0)
      {
        if (s[2] > 0 || s == Symmetry{0, 0, 0, 1})
          return itype::zero;
        else
          return itype::matrix;
      }
    if (i == 1)
      {
        if (s == Symmetry{0, 0, 1, 0})
          return itype::unit;
        else
          return itype::zero;
      }
    if (i == 2)
      return itype::zero;
    if (i == 3)
      {
        if (s == Symmetry{0, 0, 0, 1})
          return itype::unit;
        else
          return itype::zero;
      }

    TL_RAISE("Wrong stack index in GContainer::getType");
  }

};

class FoldedGContainer : public GContainer<FGSTensor>,
                         public FoldedStackContainer
{
public:
  using _Ctype = TensorContainer<FGSTensor>;
  FoldedGContainer(const _Ctype *gs, int ngs, int nu)
    : GContainer<FGSTensor>(gs, ngs, nu)
  {
  }
};

class UnfoldedGContainer : public GContainer<UGSTensor>,
                           public UnfoldedStackContainer
{
public:
  using _Ctype = TensorContainer<UGSTensor>;
  UnfoldedGContainer(const _Ctype *gs, int ngs, int nu)
    : GContainer<UGSTensor>(gs, ngs, nu)
  {
  }
};

/* Here we have a support class for product of StackContainers. It
   only adds a dimension to StackContainer. It selects the symmetries
   according to equivalence classes passed to the constructor. The
   equivalence can have permuted classes by some given
   permutation. Nothing else is interesting. */

template<class _Ttype>
class StackProduct
{
public:
  using _Stype = StackContainerInterface<_Ttype>;
  using _Ctype = typename _Stype::_Ctype;
  using itype = typename _Stype::itype;
protected:
  const _Stype &stack_cont;
  InducedSymmetries syms;
  Permutation per;
public:
  StackProduct(const _Stype &sc, const Equivalence &e,
               const Symmetry &os)
    : stack_cont(sc), syms(e, os), per(e)
  {
  }
  StackProduct(const _Stype &sc, const Equivalence &e,
               const Permutation &p, const Symmetry &os)
    : stack_cont(sc), syms(e, p, os), per(e, p)
  {
  }
  int
  dimen() const
  {
    return syms.size();
  }
  int
  getAllSize() const
  {
    return stack_cont.getAllSize();
  }
  const Symmetry &
  getProdSym(int ip) const
  {
    return syms[ip];
  }
  bool
  isZero(const IntSequence &istacks) const
  {
    TL_RAISE_IF(istacks.size() != dimen(),
                "Wrong istacks coordinates for StackProduct::isZero");

    bool res = false;
    int i = 0;
    while (i < dimen() && !(res = stack_cont.isZero(istacks[i], syms[i])))
      i++;
    return res;
  }

  itype
  getType(int is, int ip) const
  {
    TL_RAISE_IF(is < 0 || is >= stack_cont.numStacks(),
                "Wrong index to stack in StackProduct::getType");
    TL_RAISE_IF(ip < 0 || ip >= dimen(),
                "Wrong index to stack container in StackProduct::getType");
    return stack_cont.getType(is, syms[ip]);
  }

  const _Ttype &
  getMatrix(int is, int ip) const
  {
    return stack_cont.getMatrix(is, syms[ip]);
  }

  std::vector<std::unique_ptr<Vector>>
  createPackedColumns(const IntSequence &coor, IntSequence &iu) const
  {
    TL_RAISE_IF(iu.size() != dimen(),
                "Wrong storage length for unit flags in StackProduct::createPackedColumns");
    TL_RAISE_IF(coor.size() != per.size(),
                "Wrong size of index coor in StackProduct::createPackedColumns");
    IntSequence perindex(coor.size());
    std::vector<std::unique_ptr<Vector>> vs;
    per.apply(coor, perindex);
    int off = 0;
    for (int i = 0; i < dimen(); i++)
      {
        IntSequence percoor(perindex, off, syms[i].dimen() + off);
        vs.push_back(stack_cont.createPackedColumn(syms[i], percoor, iu[i]));
        off += syms[i].dimen();
      }
    return vs;
  }

  int
  getSize(int is) const
  {
    return stack_cont.getStackSizes()[is];
  }

  int
  numMatrices(const IntSequence &istacks) const
  {
    TL_RAISE_IF(istacks.size() != dimen(),
                "Wrong size of stack coordinates in StackContainer::numMatrices");
    int ret = 0;
    int ip = 0;
    while (ip < dimen() && getType(istacks[ip], ip) == _Stype::matrix)
      {
        ret++;
        ip++;
      }
    return ret;
  }
};

/* Here we only inherit from Kronecker product KronProdAllOptim, only to
   allow for a constructor constructing from StackProduct. */

template<class _Ttype>
class KronProdStack : public KronProdAllOptim
{
public:
  using _Ptype = StackProduct<_Ttype>;
  using _Stype = StackContainerInterface<_Ttype>;

  /* Here we construct KronProdAllOptim from StackContainer and given
     selections of stack items from stack containers in the product. We
     only decide whether to insert matrix, or unit matrix.

     At this point, we do not call KronProdAllOptim::optimizeOrder(), so
     the KronProdStack behaves like KronProdAll (i.e. no optimization
     is done). */

  KronProdStack(const _Ptype &sp, const IntSequence &istack)
    : KronProdAllOptim(sp.dimen())
  {
    TL_RAISE_IF(sp.dimen() != istack.size(),
                "Wrong stack product dimension for KronProdStack constructor");

    for (int i = 0; i < sp.dimen(); i++)
      {
        TL_RAISE_IF(sp.getType(istack[i], i) == _Stype::itype::zero,
                    "Attempt to construct KronProdStack from zero matrix");
        if (sp.getType(istack[i], i) == _Stype::itype::unit)
          setUnit(i, sp.getSize(istack[i]));
        if (sp.getType(istack[i], i) == _Stype::itype::matrix)
          {
            const TwoDMatrix &m = sp.getMatrix(istack[i], i);
            TL_RAISE_IF(m.nrows() != sp.getSize(istack[i]),
                        "Wrong size of returned matrix in KronProdStack constructor");
            setMat(i, m);
          }
      }
  }
};

class WorkerFoldMAADense : public sthread::detach_thread
{
  const FoldedStackContainer &cont;
  Symmetry sym;
  const FGSContainer &dense_cont;
  FGSTensor &out;
public:
  WorkerFoldMAADense(const FoldedStackContainer &container,
                     Symmetry s,
                     const FGSContainer &dcontainer,
                     FGSTensor &outten);
  void operator()(std::mutex &mut) override;
};

class WorkerFoldMAASparse1 : public sthread::detach_thread
{
  const FoldedStackContainer &cont;
  const FSSparseTensor &t;
  FGSTensor &out;
  IntSequence coor;
public:
  WorkerFoldMAASparse1(const FoldedStackContainer &container,
                       const FSSparseTensor &ten,
                       FGSTensor &outten, IntSequence c);
  void operator()(std::mutex &mut) override;
};

class WorkerFoldMAASparse2 : public sthread::detach_thread
{
  const FoldedStackContainer &cont;
  const FSSparseTensor &t;
  FGSTensor &out;
  IntSequence coor;
public:
  WorkerFoldMAASparse2(const FoldedStackContainer &container,
                       const FSSparseTensor &ten,
                       FGSTensor &outten, IntSequence c);
  void operator()(std::mutex &mut) override;
};

class WorkerFoldMAASparse4 : public sthread::detach_thread
{
  const FoldedStackContainer &cont;
  const FSSparseTensor &t;
  FGSTensor &out;
  IntSequence coor;
public:
  WorkerFoldMAASparse4(const FoldedStackContainer &container,
                       const FSSparseTensor &ten,
                       FGSTensor &outten, IntSequence c);
  void operator()(std::mutex &mut) override;
};

class WorkerUnfoldMAADense : public sthread::detach_thread
{
  const UnfoldedStackContainer &cont;
  Symmetry sym;
  const UGSContainer &dense_cont;
  UGSTensor &out;
public:
  WorkerUnfoldMAADense(const UnfoldedStackContainer &container,
                       Symmetry s,
                       const UGSContainer &dcontainer,
                       UGSTensor &outten);
  void operator()(std::mutex &mut) override;
};

class WorkerUnfoldMAASparse1 : public sthread::detach_thread
{
  const UnfoldedStackContainer &cont;
  const FSSparseTensor &t;
  UGSTensor &out;
  IntSequence coor;
public:
  WorkerUnfoldMAASparse1(const UnfoldedStackContainer &container,
                         const FSSparseTensor &ten,
                         UGSTensor &outten, IntSequence c);
  void operator()(std::mutex &mut) override;
};

class WorkerUnfoldMAASparse2 : public sthread::detach_thread
{
  const UnfoldedStackContainer &cont;
  const FSSparseTensor &t;
  UGSTensor &out;
  IntSequence coor;
public:
  WorkerUnfoldMAASparse2(const UnfoldedStackContainer &container,
                         const FSSparseTensor &ten,
                         UGSTensor &outten, IntSequence c);
  void operator()(std::mutex &mut) override;
};

#endif
