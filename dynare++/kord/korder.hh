/*
 * Copyright © 2004 Ondra Kamenik
 * Copyright © 2019-2021 Dynare Team
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

// Higher order at deterministic steady

/* The main purpose of this file is to implement a perturbation method
   algorithm for an DSGE model for higher order approximations. The input of
   the algorithm are sparse tensors as derivatives of the dynamic system, then
   dimensions of vector variables, then the first order approximation to the
   decision rule and finally a covariance matrix of exogenous shocks. The
   output are higher order derivatives of decision rule yₜ=g(y*ₜ₋₁,uₜ,σ). The
   class provides also a method for checking a size of residuals of the solved
   equations.

   The algorithm is implemented in KOrder class. The class contains both
   unfolded and folded containers to allow for switching (usually from unfold
   to fold) during the calculations. The algorithm is implemented in a few
   templated methods. To do this, we need some container type traits, which are
   in ‘ctraits’ struct. Also, the KOrder class contains some information
   encapsulated in other classes, which are defined here. These include:
   PartitionY, MatrixA, MatrixS and MatrixB.
*/

#ifndef KORDER_H
#define KORDER_H

#include "int_sequence.hh"
#include "fs_tensor.hh"
#include "gs_tensor.hh"
#include "t_container.hh"
#include "stack_container.hh"
#include "normal_moments.hh"
#include "t_polynomial.hh"
#include "faa_di_bruno.hh"
#include "journal.hh"
#include "pascal_triangle.hh"

#include "kord_exception.hh"
#include "GeneralSylvester.hh"

#include <dynlapack.h>

#include <cmath>
#include <type_traits>

// The enum class passed as template parameter for many data structures
enum class Storage { fold, unfold };

/* In ‘ctraits’ we define a number of types. We have a type for tensor Ttensor,
   and types for each pair of folded/unfolded containers used as a member in
   KOrder. */

class FoldedZXContainer;
class UnfoldedZXContainer;
class FoldedGXContainer;
class UnfoldedGXContainer;

template<Storage type>
class ctraits
{
public:
  using Ttensor = std::conditional_t<type == Storage::fold, FGSTensor, UGSTensor>;
  using Ttensym = std::conditional_t<type == Storage::fold, FFSTensor, UFSTensor>;
  using Tg = std::conditional_t<type == Storage::fold, FGSContainer, UGSContainer>;
  using Tgs = std::conditional_t<type == Storage::fold, FGSContainer, UGSContainer>;
  using Tgss = std::conditional_t<type == Storage::fold, FGSContainer, UGSContainer>;
  using TG = std::conditional_t<type == Storage::fold, FGSContainer, UGSContainer>;
  using TU = std::conditional_t<type == Storage::fold, FGSContainer, UGSContainer>;
  using TW = std::conditional_t<type == Storage::fold, FGSContainer, UGSContainer>;
  using TWrond = std::conditional_t<type == Storage::fold, FGSContainer, UGSContainer>;
  using TZstack = std::conditional_t<type == Storage::fold, FoldedZContainer, UnfoldedZContainer>;
  using TXstack = std::conditional_t<type == Storage::fold, FoldedXContainer, UnfoldedXContainer>;
  using TGstack = std::conditional_t<type == Storage::fold, FoldedGContainer, UnfoldedGContainer>;
  using Tm = std::conditional_t<type == Storage::fold, FNormalMoments, UNormalMoments>;
  using Tpol = std::conditional_t<type == Storage::fold, FTensorPolynomial, UTensorPolynomial>;
  using TZXstack = std::conditional_t<type == Storage::fold, FoldedZXContainer, UnfoldedZXContainer>;
  using TGXstack = std::conditional_t<type == Storage::fold, FoldedGXContainer, UnfoldedGXContainer>;
};

/* The PartitionY class defines the partitioning of state variables y. The
   vector y, and subvectors y* and y** are defined as:

        ⎡ static⎤
        ⎢  pred ⎥
    y = ⎢  both ⎥        ⎡pred⎤         ⎡  both ⎤
        ⎣forward⎦   y* = ⎣both⎦   y** = ⎣forward⎦

   where ‘static’ means variables appearing only at time t, ‘pred’ means
   variables appearing at time t−1, but not at t+1, ‘both’ means variables
   appearing both at t−1 and t+1 (regardless of appearance at t), and ‘forward’
   means variables appearing at t+1 but not at t−1.

   The class maintains the four lengths, and returns the whole length, length
   of y*, and length of y**$.
*/

struct PartitionY
{
  const int nstat;
  const int npred;
  const int nboth;
  const int nforw;
  PartitionY() : PartitionY(0,0,0,0) {}
  PartitionY(int num_stat, int num_pred,
             int num_both, int num_forw)
    : nstat(num_stat), npred(num_pred),
      nboth(num_both), nforw(num_forw)
  {
  }
  int
  ny() const
  {
    return nstat+npred+nboth+nforw;
  }
  int
  nys() const
  {
    return npred+nboth;
  }
  int
  nyss() const
  {
    return nboth+nforw;
  }
};

/* This is an abstraction for a square matrix with attached PLU factorization.
   It can calculate the PLU factorization and apply the inverse with some given
   matrix.

   We use LAPACK PLU decomposition for the inverse. We store the L and U in the
   ‘inv’ array and ‘ipiv’ is the permutation P. */

class PLUMatrix : public TwoDMatrix
{
public:
  PLUMatrix(int n)
    : TwoDMatrix(n, n),
      inv(nrows()*ncols()),
      ipiv(nrows())
  {
  }
  void multInv(TwoDMatrix &m) const;
private:
  Vector inv;
  std::vector<lapack_int> ipiv;
protected:
  void calcPLU();
};

/* The class MatrixA is used for matrix [f_y] + [0 [f_y**₊]·[g**_y*] 0]
   which is central for the perturbation method step. */

class MatrixA : public PLUMatrix
{
public:
  MatrixA(const FSSparseTensor &f, const IntSequence &ss,
          const TwoDMatrix &gy, const PartitionY &ypart);
};

/* The class MatrixS slightly differs from MatrixA. It is used for
   [f_y] + [0 [f_y**₊]·[g**_y*] 0] + [0 0 [f_y**₊]]
   which is needed when recovering g_σᵏ. */

class MatrixS : public PLUMatrix
{
public:
  MatrixS(const FSSparseTensor &f, const IntSequence &ss,
          const TwoDMatrix &gy, const PartitionY &ypart);
};

/* The B matrix is equal to [f_y**₊]. We have just a constructor. */

class MatrixB : public TwoDMatrix
{
public:
  MatrixB(const FSSparseTensor &f, const IntSequence &ss)
    : TwoDMatrix(FGSTensor(f, ss, IntSequence(1, 0),
                           TensorDimens(ss, IntSequence(1, 0))))
  {
  }
};

/* Class for the higher order approximations.

   Most of the code is templated, and template types are calculated in
   ‘ctraits’. So all templated methods get a template argument ‘T’, which
   can be either ‘Storage::fold’, or ‘Storage::unfold’.
*/

class KOrder
{
protected:
  /* Variable sizes: PartitionY struct maintaining partitions of y, see
     PartitionY struct declaration */
  const PartitionY ypart;

  const int ny;
  const int nu;
  const int maxk;

  /* Tensor variable dimension: variable sizes of all tensors in
     containers, sizes of y*, u, u′ and σ */
  IntSequence nvs;

  /* Tensor containers: folded and unfolded containers for g, g_y*, g_y** (the
     latter two collect appropriate subtensors of g, they do not allocate any
     new space), G, G stack, Z stack

     The names are not important because they do not appear anywhere else since
     we access them by template functions. */
  UGSContainer _ug;
  FGSContainer _fg;
  UGSContainer _ugs;
  FGSContainer _fgs;
  UGSContainer _ugss;
  FGSContainer _fgss;
  UGSContainer _uG;
  FGSContainer _fG;
  UnfoldedZContainer _uZstack;
  FoldedZContainer _fZstack;
  UnfoldedGContainer _uGstack;
  FoldedGContainer _fGstack;

  /* Moments: both folded and unfolded normal moment containers, both are
     calculated at initialization */
  UNormalMoments _um;
  FNormalMoments _fm;

  /* Dynamic model derivatives: just a reference to the container of sparse
     tensors of the system derivatives, lives outside the class */
  const TensorContainer<FSSparseTensor> &f;

  /* Matrices: matrix A, matrix S, and matrix B, see MatrixA and MatrixB class
     declarations */
  const MatrixA matA;
  const MatrixS matS;
  const MatrixB matB;

  /* These are the declarations of the template functions accessing the
     containers. We declare template methods for accessing containers depending
     on ‘fold’ and ‘unfold’ flag, we implement their specializations*/
  template<Storage t>
  typename ctraits<t>::Tg &g();
  template<Storage t>
  const typename ctraits<t>::Tg &g() const;

  template<Storage t>
  typename ctraits<t>::Tgs &gs();
  template<Storage t>
  const typename ctraits<t>::Tgs &gs() const;

  template<Storage t>
  typename ctraits<t>::Tgss &gss();
  template<Storage t>
  const typename ctraits<t>::Tgss &gss() const;

  template<Storage t>
  typename ctraits<t>::TG &G();
  template<Storage t>
  const typename ctraits<t>::TG &G() const;

  template<Storage t>
  typename ctraits<t>::TZstack &Zstack();
  template<Storage t>
  const typename ctraits<t>::TZstack &Zstack() const;

  template<Storage t>
  typename ctraits<t>::TGstack &Gstack();
  template<Storage t>
  const typename ctraits<t>::TGstack &Gstack() const;

  template<Storage t>
  typename ctraits<t>::Tm &m();
  template<Storage t>
  const typename ctraits<t>::Tm &m() const;

  Journal &journal;
public:
  KOrder(int num_stat, int num_pred, int num_both, int num_forw,
         const TensorContainer<FSSparseTensor> &fcont,
         const TwoDMatrix &gy, const TwoDMatrix &gu, const TwoDMatrix &v,
         Journal &jr);

  /* Performs k-order step provided that k=2 or the k−1-th step has been
     run, this is the core method */
  template<Storage t>
  void performStep(int order);

  /* Calculates residuals of all solved equations for k-order and reports their
     sizes, it is runnable after k-order performStep() has been run */
  template<Storage t>
  double check(int dim) const;

  template<Storage t>
  Vector calcStochShift(int order, double sigma) const;
  void switchToFolded();
  const PartitionY &
  getPartY() const
  {
    return ypart;
  }
  const FGSContainer &
  getFoldDers() const
  {
    return _fg;
  }
  const UGSContainer &
  getUnfoldDers() const
  {
    return _ug;
  }
  const FGSContainer &
  getFoldDersS() const
  {
    return _fgs;
  }
  static bool
  is_even(int i)
  {
    return i % 2 == 0;
  }
protected:
  /* Inserts a g derivative to the g container and also creates subtensors and
     insert them to g_y* and g_y** containers */
  template<Storage t>
  void insertDerivative(std::unique_ptr<typename ctraits<t>::Ttensor> der);

  /* Solves the sylvester equation (templated fold, and unfold) */
  template<Storage t>
  void sylvesterSolve(typename ctraits<t>::Ttensor &der) const;

  /* Calculates derivatives of F by Faà Di Bruno for the sparse container of
     system derivatives and Z stack container */
  template<Storage t>
  std::unique_ptr<typename ctraits<t>::Ttensor> faaDiBrunoZ(const Symmetry &sym) const;

  /* Calculates derivatives of G by Faà Di Bruno for the dense container g**
     and G stack */
  template<Storage t>
  std::unique_ptr<typename ctraits<t>::Ttensor> faaDiBrunoG(const Symmetry &sym) const;

  // Recovers g_y*ⁱ
  template<Storage t>
  void recover_y(int i);
  // Recovers g_y*ⁱuʲ
  template<Storage t>
  void recover_yu(int i, int j);
  // Recovers g_y*ⁱσʲ
  template<Storage t>
  void recover_ys(int i, int j);
  // Recovers g_y*ⁱuʲσᵏ
  template<Storage t>
  void recover_yus(int i, int j, int k);
  template<Storage t>
  // Recovers g_σⁱ
  void recover_s(int i);
  // Calculates specified derivatives of G and inserts them to the container
  template<Storage t>
  void fillG(int i, int j, int k);

  // Calculates Dᵢⱼₖ
  template<Storage t>
  typename ctraits<t>::Ttensor calcD_ijk(int i, int j, int k) const;
  template<Storage t>
  typename ctraits<t>::Ttensor calcD_ik(int i, int k) const;
  template<Storage t>
  typename ctraits<t>::Ttensor calcD_k(int k) const;

  // Calculates Eᵢⱼₖ
  template<Storage t>
  typename ctraits<t>::Ttensor calcE_ijk(int i, int j, int k) const;
  template<Storage t>
  typename ctraits<t>::Ttensor calcE_ik(int i, int k) const;
  template<Storage t>
  typename ctraits<t>::Ttensor calcE_k(int k) const;
};

/* Here we insert the result to the container. Along the insertion, we
   also create subtensors and insert as well. */

template<Storage t>
void
KOrder::insertDerivative(std::unique_ptr<typename ctraits<t>::Ttensor> der)
{
  auto der_ptr = der.get();
  g<t>().insert(std::move(der));
  gs<t>().insert(std::make_unique<typename ctraits<t>::Ttensor>(ypart.nstat, ypart.nys(), *der_ptr));
  gss<t>().insert(std::make_unique<typename ctraits<t>::Ttensor>(ypart.nstat+ypart.npred,
                                                                 ypart.nyss(), *der_ptr));
}

/* Here we implement Faà Di Bruno formula

    ₖ                        ₗ
    ∑  [f_zˡ]_γ₁…γₗ    ∑     ∏  [z_{s^|cₘ|}]^γₘ
   ˡ⁼¹              c∈ℳₗ,ₖ ᵐ⁼¹

   where s is a given outer symmetry and k is the dimension of the symmetry. */

template<Storage t>
std::unique_ptr<typename ctraits<t>::Ttensor>
KOrder::faaDiBrunoZ(const Symmetry &sym) const
{
  JournalRecordPair pa(journal);
  pa << u8"Faà Di Bruno Z container for " << sym << endrec;
  auto res = std::make_unique<typename ctraits<t>::Ttensor>(ny, TensorDimens(sym, nvs));
  FaaDiBruno bruno(journal);
  bruno.calculate(Zstack<t>(), f, *res);
  return res;
}

/* The same as KOrder::faaDiBrunoZ(), but for g** and G stack. */

template<Storage t>
std::unique_ptr<typename ctraits<t>::Ttensor>
KOrder::faaDiBrunoG(const Symmetry &sym) const
{
  JournalRecordPair pa(journal);
  pa << u8"Faà Di Bruno G container for " << sym << endrec;
  TensorDimens tdims(sym, nvs);
  auto res = std::make_unique<typename ctraits<t>::Ttensor>(ypart.nyss(), tdims);
  FaaDiBruno bruno(journal);
  bruno.calculate(Gstack<t>(), gss<t>(), *res);
  return res;
}

/* Here we solve [F_yⁱ]=0. First we calculate conditional G_yⁱ (it misses l=1
   and l=i since g_yⁱ does not exist yet). Then calculate conditional F_yⁱ and
   we have the right hand side of equation. Since we miss two orders, we solve
   by Sylvester, and insert the solution as the derivative g_yⁱ. Then we need
   to update G_yⁱ running multAndAdd() for both dimensions 1 and i.

   Requires: everything at order ≤ i−1

   Provides: g_yⁱ and G_yⁱ
*/
template<Storage t>
void
KOrder::recover_y(int i)
{
  Symmetry sym{i, 0, 0, 0};
  JournalRecordPair pa(journal);
  pa << "Recovering symmetry " << sym << endrec;

  auto G_yi = faaDiBrunoG<t>(sym);
  auto G_yi_ptr = G_yi.get();
  G<t>().insert(std::move(G_yi));

  auto g_yi = faaDiBrunoZ<t>(sym);
  g_yi->mult(-1.0);

  sylvesterSolve<t>(*g_yi);

  insertDerivative<t>(std::move(g_yi));

  auto &gss_y = gss<t>().get(Symmetry{1, 0, 0, 0});
  gs<t>().multAndAdd(gss_y, *G_yi_ptr);
  auto &gss_yi = gss<t>().get(sym);
  gs<t>().multAndAdd(gss_yi, *G_yi_ptr);
}

/* Here we solve [F_yⁱuʲ]=0 to obtain g_yⁱuʲ for j>0. We calculate conditional
   G_yⁱuʲ (this misses only l=1) and calculate conditional F_yⁱuʲ and we have
   the right hand side. It is solved by multiplication of inversion of A. Then
   we insert the result, and update G_yⁱuʲ by multAndAdd() for l=1.

   Requires: everything at order ≤ i+j−1, G_yⁱ⁺ʲ and g_yⁱ⁺ʲ.

   Provides: g_yⁱuʲ and G_yⁱuʲ
*/
template<Storage t>
void
KOrder::recover_yu(int i, int j)
{
  Symmetry sym{i, j, 0, 0};
  JournalRecordPair pa(journal);
  pa << "Recovering symmetry " << sym << endrec;

  auto G_yiuj = faaDiBrunoG<t>(sym);
  auto G_yiuj_ptr = G_yiuj.get();
  G<t>().insert(std::move(G_yiuj));

  auto g_yiuj = faaDiBrunoZ<t>(sym);
  g_yiuj->mult(-1.0);
  matA.multInv(*g_yiuj);
  insertDerivative<t>(std::move(g_yiuj));

  gs<t>().multAndAdd(gss<t>().get(Symmetry{1, 0, 0, 0}), *G_yiuj_ptr);
}

/* Here we solve [F_yⁱσʲ]+[Dᵢⱼ]+[Eᵢⱼ]=0 to obtain g_yⁱσʲ. We calculate
   conditional G_yⁱσʲ (missing dimensions 1 and i+j), calculate conditional
   F_yⁱσʲ. Before we can calculate Dᵢⱼ and Eᵢⱼ, we have to calculate
   G_yⁱu′ᵐσʲ⁻ᵐ for m=1,…,j. Then we add the Dᵢⱼ and Eᵢⱼ to obtain the right
   hand side. Then we solve the sylvester to obtain g_yⁱσʲ. Then we update
   G_yⁱσʲ for l=1 and l=i+j.

   Requires: everything at order ≤ i+j−1, g_yⁱ⁺ʲ, G_yⁱu′ʲ and g_yⁱuʲ through
     Dᵢⱼ, g_yⁱuᵐσʲ⁻ᵐ for m=1,…,j−1 through Eᵢⱼ.

   Provides: g_yⁱσʲ and G_yⁱσʲ, and finally G_yⁱu′ᵐσʲ⁻ᵐ for m=1,…,j. The latter
     is calculated by fillG() before the actual calculation.
*/
template<Storage t>
void
KOrder::recover_ys(int i, int j)
{
  Symmetry sym{i, 0, 0, j};
  JournalRecordPair pa(journal);
  pa << "Recovering symmetry " << sym << endrec;

  fillG<t>(i, 0, j);

  if (is_even(j))
    {
      auto G_yisj = faaDiBrunoG<t>(sym);
      auto G_yisj_ptr = G_yisj.get();
      G<t>().insert(std::move(G_yisj));

      auto g_yisj = faaDiBrunoZ<t>(sym);

      {
        auto D_ij = calcD_ik<t>(i, j);
        g_yisj->add(1.0, D_ij);
      }

      if (j >= 3)
        {
          auto E_ij = calcE_ik<t>(i, j);
          g_yisj->add(1.0, E_ij);
        }

      g_yisj->mult(-1.0);

      sylvesterSolve<t>(*g_yisj);

      insertDerivative<t>(std::move(g_yisj));

      Gstack<t>().multAndAdd(1, gss<t>(), *G_yisj_ptr);
      Gstack<t>().multAndAdd(i+j, gss<t>(), *G_yisj_ptr);
    }
}

/* Here we solve [F_yⁱuʲσᵏ]+[Dᵢⱼₖ]+[Eᵢⱼₖ]=0 to obtain g_yⁱuʲσᵏ. First we
   calculate conditional G_yⁱuʲσᵏ (missing only for dimension l=1), then we
   evaluate conditional F_yⁱuʲσᵏ. Before we can calculate Dᵢⱼₖ, and Eᵢⱼₖ, we
   need to insert G_yⁱuʲu′ᵐσᵏ⁻ᵐ for m=1,…,k. This is done by fillG(). Then we
   have right hand side and we multiply by A⁻¹ to obtain g_yⁱuʲσᵏ. Finally we
   have to update G_yⁱuʲσᵏ by multAndAdd() for dimension l=1.

   Requires: everything at order ≤ i+j+k, g_yⁱ⁺ʲσᵏ through G_yⁱuʲσᵏ involved in
     right hand side, then g_yⁱuʲ⁺ᵏ through Dᵢⱼₖ, and g_yⁱuʲ⁺ᵐσᵏ⁻ᵐ for m=1,…,k−1
     through Eᵢⱼₖ.

   Provides: g_yⁱuʲσᵏ, G_yⁱuʲσᵏ, and G_yⁱuʲu′ᵐσᵏ⁻ᵐ for m=1,…,k
*/
template<Storage t>
void
KOrder::recover_yus(int i, int j, int k)
{
  Symmetry sym{i, j, 0, k};
  JournalRecordPair pa(journal);
  pa << "Recovering symmetry " << sym << endrec;

  fillG<t>(i, j, k);

  if (is_even(k))
    {
      auto G_yiujsk = faaDiBrunoG<t>(sym);
      auto G_yiujsk_ptr = G_yiujsk.get();
      G<t>().insert(std::move(G_yiujsk));

      auto g_yiujsk = faaDiBrunoZ<t>(sym);

      {
        auto D_ijk = calcD_ijk<t>(i, j, k);
        g_yiujsk->add(1.0, D_ijk);
      }

      if (k >= 3)
        {
          auto E_ijk = calcE_ijk<t>(i, j, k);
          g_yiujsk->add(1.0, E_ijk);
        }

      g_yiujsk->mult(-1.0);

      matA.multInv(*g_yiujsk);
      insertDerivative<t>(std::move(g_yiujsk));

      Gstack<t>().multAndAdd(1, gss<t>(), *G_yiujsk_ptr);
    }
}

/* Here we solve [F_{σⁱ}]+[Dᵢ]+[Eᵢ]=0 to recover g_σⁱ. First we calculate
   conditional G_σⁱ (missing dimension l=1 and l=i), then we calculate
   conditional F_σⁱ. Before we can calculate Dᵢ and Eᵢ, we have to obtain
   G_u′ᵐσⁱ⁻ᵐ for m=1,…,i. Then adding Dᵢ and Eᵢ we have the right hand side. We
   solve by S⁻¹ multiplication and update G_σⁱ by calling multAndAdd() for
   dimension l=1.

   Recall that the solved equation here is:

    [f_y][g_σᵏ]+[f_y**₊][g**_y*][g*_σᵏ]+[f_y**₊][g**_σᵏ] = RHS

   This is a sort of deficient sylvester equation (sylvester equation for
   dimension=0), we solve it by S⁻¹. See the constructor of MatrixS to see how
   S looks like.

   Requires: everything at order ≤ i−1, g_yⁱ and g_yⁱ⁻ʲσʲ, then g_uᵏ through
     F_u′ᵏ, and g_yᵐuʲσᵏ for j=1,…,i−1 and m+j+k=i through F_u′ʲσⁱ⁻ʲ.

   Provides: g_σⁱ, G_σⁱ, and G_u′ᵐσⁱ⁻ᵐ for m=1,…,i
*/
template<Storage t>
void
KOrder::recover_s(int i)
{
  Symmetry sym{0, 0, 0, i};
  JournalRecordPair pa(journal);
  pa << "Recovering symmetry " << sym << endrec;

  fillG<t>(0, 0, i);

  if (is_even(i))
    {
      auto G_si = faaDiBrunoG<t>(sym);
      auto G_si_ptr = G_si.get();
      G<t>().insert(std::move(G_si));

      auto g_si = faaDiBrunoZ<t>(sym);

      {
        auto D_i = calcD_k<t>(i);
        g_si->add(1.0, D_i);
      }

      if (i >= 3)
        {
          auto E_i = calcE_k<t>(i);
          g_si->add(1.0, E_i);
        }

      g_si->mult(-1.0);

      matS.multInv(*g_si);
      insertDerivative<t>(std::move(g_si));

      Gstack<t>().multAndAdd(1, gss<t>(), *G_si_ptr);
      Gstack<t>().multAndAdd(i, gss<t>(), *G_si_ptr);
    }
}

/* Here we calculate and insert G_yⁱuʲu′ᵐσᵏ⁻ᵐ for m=1,…,k. The derivatives are
   inserted only for k−m being even. */

template<Storage t>
void
KOrder::fillG(int i, int j, int k)
{
  for (int m = 1; m <= k; m++)
    if (is_even(k-m))
      {
        auto G_yiujupms = faaDiBrunoG<t>(Symmetry{i, j, m, k-m});
        G<t>().insert(std::move(G_yiujupms));
      }
}

/* Here we calculate:

    [Dᵢⱼₖ]_α₁…αᵢβ₁…βⱼ = [F_yⁱuʲu′ᵏ]_α₁…αᵢβ₁…βⱼγ₁…γₖ [Σ]^γ₁…γₖ

   So it is non zero only for even k. */

template<Storage t>
typename ctraits<t>::Ttensor
KOrder::calcD_ijk(int i, int j, int k) const
{
  typename ctraits<t>::Ttensor res(ny, TensorDimens(Symmetry{i, j, 0, 0}, nvs));
  res.zeros();
  if (is_even(k))
    {
      auto tmp = faaDiBrunoZ<t>(Symmetry{i, j, k, 0});
      tmp->contractAndAdd(2, res, m<t>().get(Symmetry{k}));
    }
  return res;
}

/* Here we calculate
                        ₖ₋₁ ⎛k⎞
    [Eᵢⱼₖ]_α₁…αᵢβ₁…βⱼ =  ∑  ⎝m⎠ [F_yⁱuʲu′ᵐσᵏ⁻ᵐ]_α₁…αᵢβ₁…βⱼγ₁…γₘ [Σ]^γ₁…γₘ
                        ᵐ⁼¹
   The sum can sum only for even m. */

template<Storage t>
typename ctraits<t>::Ttensor
KOrder::calcE_ijk(int i, int j, int k) const
{
  typename ctraits<t>::Ttensor res(ny, TensorDimens(Symmetry{i, j, 0, 0}, nvs));
  res.zeros();
  for (int n = 2; n <= k-1; n += 2)
    {
      auto tmp = faaDiBrunoZ<t>(Symmetry{i, j, n, k-n});
      tmp->mult(static_cast<double>(PascalTriangle::noverk(k, n)));
      tmp->contractAndAdd(2, res, m<t>().get(Symmetry{n}));
    }
  return res;
}

template<Storage t>
typename ctraits<t>::Ttensor
KOrder::calcD_ik(int i, int k) const
{
  return calcD_ijk<t>(i, 0, k);
}

template<Storage t>
typename ctraits<t>::Ttensor
KOrder::calcD_k(int k) const
{
  return calcD_ijk<t>(0, 0, k);
}

template<Storage t>
typename ctraits<t>::Ttensor
KOrder::calcE_ik(int i, int k) const
{
  return calcE_ijk<t>(i, 0, k);
}

template<Storage t>
typename ctraits<t>::Ttensor
KOrder::calcE_k(int k) const
{
  return calcE_ijk<t>(0, 0, k);
}

/* Here is the core routine. It calls methods recovering derivatives in the
   right order. Recall, that the code, namely Faà Di Bruno’s formula, is
   implemented as to be run conditionally on the current contents of
   containers. So, if some call of Faà Di Bruno evaluates derivatives, and some
   derivatives are not present in the container, then it is considered to be
   zero. So, we have to be very careful to put everything in the right order.
   The order here can be derived from dependencies, or it is in the paper.

   The method recovers all the derivatives of the given ‘order’.

   The precondition of the method is that all tensors of order ‘order-1’, which
   are not zero, exist (including G). The postcondition of of the method is
   derivatives of g and G of order ‘order’ are calculated and stored in the
   containers. Responsibility of precondition lays upon the constructor (for
   ‘order==2’), or upon the previous call of performStep().

   From the code, it is clear, that all g are calculated. If one goes through
   all the recovering methods, he should find out that also all G are
   provided. */

template<Storage t>
void
KOrder::performStep(int order)
{
  KORD_RAISE_IF(order-1 != g<t>().getMaxDim(),
                "Wrong order for KOrder::performStep");
  JournalRecordPair pa(journal);
  pa << "Performing step for order = " << order << endrec;

  recover_y<t>(order);

  for (int i = 0; i < order; i++)
    recover_yu<t>(i, order-i);

  for (int j = 1; j < order; j++)
    {
      for (int i = j-1; i >= 1; i--)
        recover_yus<t>(order-j, i, j-i);
      recover_ys<t>(order-j, j);
    }

  for (int i = order-1; i >= 1; i--)
    recover_yus<t>(0, i, order-i);

  recover_s<t>(order);
}

/* Here we check for residuals of all the solved equations at the given order.
   The method returns the largest residual size. Each check simply evaluates
   the equation. */

template<Storage t>
double
KOrder::check(int dim) const
{
  KORD_RAISE_IF(dim > g<t>().getMaxDim(),
                "Wrong dimension for KOrder::check");
  JournalRecordPair pa(journal);
  pa << "Checking residuals for order = " << dim << endrec;

  double maxerror = 0.0;

  // Check for F_yⁱuʲ=0
  for (int i = 0; i <= dim; i++)
    {
      Symmetry sym{dim-i, i, 0, 0};
      auto r = faaDiBrunoZ<t>(sym);
      double err = r->getData().getMax();
      JournalRecord(journal) << "\terror for symmetry " << sym << "\tis " << err << endrec;
      maxerror = std::max(err, maxerror);
    }

  // Check for F_yⁱuʲu′ᵏ+Dᵢⱼₖ+Eᵢⱼₖ=0
  for (auto &si : SymmetrySet(dim, 3))
    {
      int i = si[0];
      int j = si[1];
      int k = si[2];
      if (i+j > 0 && k > 0)
        {
          Symmetry sym{i, j, 0, k};
          auto r = faaDiBrunoZ<t>(sym);
          auto D_ijk = calcD_ijk<t>(i, j, k);
          r->add(1.0, D_ijk);
          auto E_ijk = calcE_ijk<t>(i, j, k);
          r->add(1.0, E_ijk);
          double err = r->getData().getMax();
          JournalRecord(journal) << "\terror for symmetry " << sym << "\tis " << err << endrec;
          maxerror = std::max(err, maxerror);
        }
    }

  // Check for F_σⁱ+Dᵢ+Eᵢ=0
  auto r = faaDiBrunoZ<t>(Symmetry{0, 0, 0, dim});
  auto D_k = calcD_k<t>(dim);
  r->add(1.0, D_k);
  auto E_k = calcE_k<t>(dim);
  r->add(1.0, E_k);
  double err = r->getData().getMax();
  Symmetry sym{0, 0, 0, dim};
  JournalRecord(journal) << "\terror for symmetry " << sym << "\tis " << err << endrec;
  maxerror = std::max(err, maxerror);

  return maxerror;
}

template<Storage t>
Vector
KOrder::calcStochShift(int order, double sigma) const
{
  Vector res(ny);
  res.zeros();
  int jfac = 1;
  for (int j = 1; j <= order; j++, jfac *= j)
    if (is_even(j))
      {
        auto ten = calcD_k<t>(j);
        res.add(std::pow(sigma, j)/jfac, ten.getData());
      }
  return res;
}

#endif
