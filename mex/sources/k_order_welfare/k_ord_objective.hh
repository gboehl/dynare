/*
 * Copyright © 2021-2024 Dynare Team
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

#ifndef K_ORD_OBJECTIVE_HH
#define K_ORD_OBJECTIVE_HH

#include "k_ord_dynare.hh"
#include "objective_m.hh"

class KordwDynare;

class KordwDynare
{
public:
  KordpDynare& model;
  const ConstVector& NNZD;

private:
  Journal& journal;
  Vector& params;
  Vector resid;
  TensorContainer<FSSparseTensor> ud; // planner's objective derivatives, in Dynare++ form
  std::vector<int> dynppToDyn;        // Maps Dynare++ jacobian variable indices to Dynare ones
  std::vector<int> dynToDynpp;        // Maps Dynare jacobian variable indices to Dynare++ ones
  std::unique_ptr<ObjectiveMFile> objectiveFile;

public:
  KordwDynare(KordpDynare& m, ConstVector& NNZD_arg, Journal& jr, Vector& inParams,
              std::unique_ptr<ObjectiveMFile> objectiveFile_arg, const std::vector<int>& varOrder);
  void calcDerivativesAtSteady();
  void populateDerivativesContainer(const std::vector<TwoDMatrix>& dyn_ud, int ord);
  [[nodiscard]] const TensorContainer<FSSparseTensor>&
  getPlannerObjDerivatives() const
  {
    return ud;
  }
  [[nodiscard]] const KordpDynare&
  getModel() const
  {
    return model;
  }
  [[nodiscard]] const Vector&
  getResid() const
  {
    return resid;
  }

private:
  /* Computes the permutations mapping back and forth between Dynare and
     Dynare++ orderings of variables */
  void computeJacobianPermutation(const std::vector<int>& var_order);
};

class KOrderWelfare
{
public:
  const PartitionY ypart;
  const int ny;
  const int nu;
  const int maxk;
  const int order;
  const double discount_factor;
  /* Tensor variable dimension: variable sizes of all tensors in
     containers, sizes of y, u and σ */
  IntSequence nvs;
  UGSContainer _uU;
  FGSContainer _fU;
  UGSContainer _uW;
  FGSContainer _fW;
  UGSContainer _uWrond;
  FGSContainer _fWrond;
  UGSContainer _ug;
  FGSContainer _fg;
  UGSContainer _ugs;
  FGSContainer _fgs;
  UnfoldedXContainer _uXstack;
  FoldedXContainer _fXstack;
  UnfoldedGContainer _uGstack;
  FoldedGContainer _fGstack;

  /* Moments: both folded and unfolded normal moment containers, both are
     calculated at initialization */
  UNormalMoments _um;
  FNormalMoments _fm;

  /* Planner objective derivatives: just a reference to the container of sparse
     tensors of the derivatives, lives outside the class */
  const TensorContainer<FSSparseTensor>& u;

  /* These are the declarations of the template functions accessing the
     containers. We declare template methods for accessing containers depending
     on ‘fold’ and ‘unfold’ flag, we implement their specializations*/
  template<Storage t>
  typename ctraits<t>::Tg& g();
  template<Storage t>
  const typename ctraits<t>::Tg& g() const;
  template<Storage t>
  typename ctraits<t>::Tgs& gs();
  template<Storage t>
  const typename ctraits<t>::Tgs& gs() const;
  template<Storage t>
  typename ctraits<t>::TU& U();
  template<Storage t>
  const typename ctraits<t>::TU& U() const;
  template<Storage t>
  typename ctraits<t>::TW& W();
  template<Storage t>
  const typename ctraits<t>::TW& W() const;
  template<Storage t>
  typename ctraits<t>::TWrond& Wrond();
  template<Storage t>
  const typename ctraits<t>::TWrond& Wrond() const;
  template<Storage t>
  typename ctraits<t>::TXstack& Xstack();
  template<Storage t>
  const typename ctraits<t>::TXstack& Xstack() const;
  template<Storage t>
  typename ctraits<t>::TGstack& Gstack();
  template<Storage t>
  const typename ctraits<t>::TGstack& Gstack() const;
  template<Storage t>
  typename ctraits<t>::Tm& m();
  template<Storage t>
  const typename ctraits<t>::Tm& m() const;

  Journal& journal;

public:
  KOrderWelfare(int num_stat, int num_pred, int num_both, int num_forw, int nu, int ord,
                double discount_factor, const TensorContainer<FSSparseTensor>& ucont,
                FGSContainer g, FGSContainer gs, const TwoDMatrix& v, Journal& jr);

  /* Performs k-order step provided that k=2 or the k−1-th step has been
     run, this is the core method */
  template<Storage t>
  void performStep(int order);
  /* Calculates residuals of all solved equations for k-order and reports their
     sizes, it is runnable after k-order performStep() has been run */
  template<Storage t>
  [[nodiscard]] double check(int dim) const;

  [[nodiscard]] const FGSContainer&
  getFoldU() const
  {
    return _fU;
  }
  [[nodiscard]] const UGSContainer&
  getUnfoldU() const
  {
    return _uU;
  }
  [[nodiscard]] const FGSContainer&
  getFoldW() const
  {
    return _fW;
  }
  [[nodiscard]] const UGSContainer&
  getUnfoldW() const
  {
    return _uW;
  }
  static bool
  is_even(int i)
  {
    return i % 2 == 0;
  }

protected:
  /* Calculates derivatives of U by Faà Di Bruno for the sparse container of
     planner's objective derivatives and the container for the decision rule derivatives*/
  template<Storage t>
  std::unique_ptr<typename ctraits<t>::Ttensor> faaDiBrunoU(const Symmetry& sym) const;
  /* Calculates derivatives of the compounded functions W and Gstack using the Faà Di Bruno
   * formula*/
  template<Storage t>
  std::unique_ptr<typename ctraits<t>::Ttensor> faaDiBrunoW(const Symmetry& sym) const;

  /* Solves the sylvester equation (templated fold, and unfold) */
  template<Storage t>
  void sylvesterSolve(typename ctraits<t>::Ttensor& der) const;

  // Recovers W_y*ⁱ
  template<Storage t>
  void recover_y(int i);
  // Recovers W_y*ⁱuʲ
  template<Storage t>
  void recover_yu(int i, int j);
  // Recovers W_y*ⁱσʲ
  template<Storage t>
  void recover_ys(int i, int j);
  // Recovers W_y*ⁱuʲσᵏ
  template<Storage t>
  void recover_yus(int i, int j, int k);
  // Recovers W_σⁱ
  template<Storage t>
  void recover_s(int i);

  // Calculates specified derivatives of Wrond and inserts them to the container
  template<Storage t>
  void fillWrond(int i, int j, int k);

  // Calculates Dᵢⱼₖ
  template<Storage t>
  typename ctraits<t>::Ttensor calcH_ijk(int i, int j, int k) const;
  template<Storage t>
  typename ctraits<t>::Ttensor calcH_ik(int i, int k) const;
  template<Storage t>
  typename ctraits<t>::Ttensor calcH_k(int k) const;

  // Calculates Eᵢⱼₖ
  template<Storage t>
  typename ctraits<t>::Ttensor calcJ_ijk(int i, int j, int k) const;
  template<Storage t>
  typename ctraits<t>::Ttensor calcJ_ik(int i, int k) const;
  template<Storage t>
  typename ctraits<t>::Ttensor calcJ_k(int k) const;
};

/* Here we implement Faà Di Bruno formula

    ₖ                        ₗ
    ∑  [u_zˡ]_γ₁…γₗ    ∑     ∏  [h_{s^|cₘ|}]^γₘ
   ˡ⁼¹              c∈ℳₗ,ₖ ᵐ⁼¹

   where s is a given outer symmetry and k is the dimension of the symmetry. */
template<Storage t>
std::unique_ptr<typename ctraits<t>::Ttensor>
KOrderWelfare::faaDiBrunoU(const Symmetry& sym) const
{
  JournalRecordPair pa(journal);
  pa << "Faà Di Bruno U container for " << sym << endrec;
  auto res = std::make_unique<typename ctraits<t>::Ttensor>(1, TensorDimens(sym, nvs));
  FaaDiBruno bruno(journal);
  bruno.calculate(Xstack<t>(), u, *res);
  return res;
}

/* The same as KOrder::faaDiBrunoW(), but for W and G stack. */
template<Storage t>
std::unique_ptr<typename ctraits<t>::Ttensor>
KOrderWelfare::faaDiBrunoW(const Symmetry& sym) const
{
  JournalRecordPair pa(journal);
  pa << "Faà Di Bruno G container for " << sym << endrec;
  TensorDimens tdims(sym, nvs);
  auto res = std::make_unique<typename ctraits<t>::Ttensor>(1, tdims);
  FaaDiBruno bruno(journal);
  bruno.calculate(Gstack<t>(), W<t>(), *res);
  return res;
}

template<Storage t>
void
KOrderWelfare::performStep(int order)
{
  KORD_RAISE_IF(order - 1 != W<t>().getMaxDim() and order > 1,
                "Wrong order for KOrder::performStep");
  JournalRecordPair pa(journal);
  pa << "Performing step for order = " << order << endrec;

  recover_y<t>(order);

  for (int i = 0; i < order; i++)
    recover_yu<t>(i, order - i);

  recover_ys<t>(order - 1, 1);

  for (int j = 2; j < order; j++)
    {
      for (int i = 1; i <= j - 1; i++)
        recover_yus<t>(order - j, i, j - i);
      recover_ys<t>(order - j, j);
      recover_yus<t>(0, order - j, j);
    }

  recover_s<t>(order);
}

/* Here we solve [F_yⁱ]=0. First we calculate conditional W_yⁱ (it misses l=i since W_yⁱ does not
   exist yet). Then calculate conditional F_yⁱ and we have the right hand side of equation. Since we
   miss two orders, we solve by Sylvester, and  the solution as the derivative g_yⁱ. Then we need to
   update G_yⁱ running multAndAdd() for both dimensions 1 and i.

   Requires: everything at order ≤ i−1

   Provides: W_yⁱ and Wrond_yⁱ
*/
template<Storage t>
void
KOrderWelfare::recover_y(int i)
{
  Symmetry sym {i, 0, 0, 0};
  JournalRecordPair pa(journal);
  pa << "Conditional welfare: recovering symmetry " << sym << "\n" << endrec;

  auto Wrond_yi = faaDiBrunoW<t>(sym);
  auto Wrond_yi_ptr = Wrond_yi.get();
  Wrond<t>().insert(std::move(Wrond_yi));

  auto W_yi = U<t>().get(sym);
  W_yi.add(discount_factor, *Wrond_yi_ptr);

  sylvesterSolve<t>(W_yi);

  W<t>().insert(std::make_unique<typename ctraits<t>::Ttensor>(W_yi));

  Gstack<t>().multAndAdd(i, W<t>(), *Wrond_yi_ptr);
}

/* Here we solve [F_yⁱuʲ]=0 to obtain W_yⁱuʲ for j>0. We calculate conditional
   Wrond_yⁱuʲ and calculate conditional F_yⁱuʲ and we have
   the right hand side.

   Requires: everything at order ≤ i+j−1, Wrond_yⁱ⁺ʲ and W_yⁱ⁺ʲ.

   Provides: W_yⁱuʲ and Wrond_yⁱuʲ
*/
template<Storage t>
void
KOrderWelfare::recover_yu(int i, int j)
{
  Symmetry sym {i, j, 0, 0};
  JournalRecordPair pa(journal);
  pa << "Conditional welfare: recovering symmetry " << sym << endrec;

  auto Wrond_yiuj = faaDiBrunoW<t>(sym);
  auto Wrond_yiuj_ptr = Wrond_yiuj.get();
  Wrond<t>().insert(std::move(Wrond_yiuj));

  auto W_yiuj = U<t>().get(sym);
  W_yiuj.add(discount_factor, *Wrond_yiuj_ptr);

  W<t>().insert(std::make_unique<typename ctraits<t>::Ttensor>(W_yiuj));
}

/* Here we solve [F_yⁱσʲ]+[Dᵢⱼ]+[Eᵢⱼ]=0 to obtain W_yⁱσʲ. We calculate
   conditional Wrond_yⁱσʲ (missing dimensions i+j), calculate conditional
   F_yⁱσʲ. Before we can calculate Dᵢⱼ and Eᵢⱼ, we have to calculate
   Wrond_yⁱu′ᵐσʲ⁻ᵐ for m=1,…,j. Then we add the Dᵢⱼ and Eᵢⱼ to obtain the right
   hand side. Then we solve the sylvester to obtain g_yⁱσʲ. Then we update
   G_yⁱσʲ for l=1 and l=i+j.

   Requires: everything at order ≤ i+j−1, g_yⁱ⁺ʲ, G_yⁱu′ʲ and g_yⁱuʲ through
     Dᵢⱼ, g_yⁱuᵐσʲ⁻ᵐ for m=1,…,j−1 through Eᵢⱼ.

   Provides: g_yⁱσʲ and G_yⁱσʲ, and finally G_yⁱu′ᵐσʲ⁻ᵐ for m=1,…,j. The latter
     is calculated by fillG() before the actual calculation.
*/
template<Storage t>
void
KOrderWelfare::recover_ys(int i, int j)
{
  Symmetry sym {i, 0, 0, j};
  JournalRecordPair pa(journal);
  pa << "Conditional welfare: recovering symmetry " << sym << endrec;

  fillWrond<t>(i, 0, j);

  if (is_even(j))
    {
      auto Wrond_yisj = faaDiBrunoW<t>(sym);
      auto Wrond_yisj_ptr = Wrond_yisj.get();
      Wrond<t>().insert(std::move(Wrond_yisj));

      auto W_yisj = U<t>().get(sym);
      W_yisj.add(discount_factor, *Wrond_yisj_ptr);
      {
        auto H_ij = calcH_ik<t>(i, j);
        W_yisj.add(-1.0, H_ij);
      }

      if (j >= 3)
        {
          auto J_ij = calcJ_ik<t>(i, j);
          W_yisj.add(-1.0, J_ij);
        }

      sylvesterSolve<t>(W_yisj);

      W<t>().insert(std::make_unique<typename ctraits<t>::Ttensor>(W_yisj));

      Gstack<t>().multAndAdd(i + j, W<t>(), *Wrond_yisj_ptr);
    }
}

/* Here we solve [F_yⁱuʲσᵏ]+[Dᵢⱼₖ]+[Eᵢⱼₖ]=0 to obtain g_yⁱuʲσᵏ. First we
   calculate conditional G_yⁱuʲσᵏ (missing only for dimension l=1), then we
   evaluate conditional F_yⁱuʲσᵏ. Before we can calculate Dᵢⱼₖ, and Eᵢⱼₖ, we
   need to  G_yⁱuʲu′ᵐσᵏ⁻ᵐ for m=1,…,k. This is done by fillG(). Then we
   have right hand side and we multiply by A⁻¹ to obtain g_yⁱuʲσᵏ. Finally we
   have to update G_yⁱuʲσᵏ by multAndAdd() for dimension l=1.

   Requires: everything at order ≤ i+j+k, g_yⁱ⁺ʲσᵏ through G_yⁱuʲσᵏ involved in
     right hand side, then g_yⁱuʲ⁺ᵏ through Dᵢⱼₖ, and g_yⁱuʲ⁺ᵐσᵏ⁻ᵐ for m=1,…,k−1
     through Eᵢⱼₖ.

   Provides: g_yⁱuʲσᵏ, G_yⁱuʲσᵏ, and G_yⁱuʲu′ᵐσᵏ⁻ᵐ for m=1,…,k
*/
template<Storage t>
void
KOrderWelfare::recover_yus(int i, int j, int k)
{
  Symmetry sym {i, j, 0, k};
  JournalRecordPair pa(journal);
  pa << "Conditional welfare: recovering symmetry " << sym << endrec;

  fillWrond<t>(i, j, k);

  if (is_even(k))
    {
      auto Wrond_yiujsk = faaDiBrunoW<t>(sym);
      auto Wrond_yiujsk_ptr = Wrond_yiujsk.get();
      Wrond<t>().insert(std::move(Wrond_yiujsk));

      auto W_yiujsk = U<t>().get(sym);
      W_yiujsk.add(discount_factor, *Wrond_yiujsk_ptr);

      auto H_ijk = calcH_ijk<t>(i, j, k);
      W_yiujsk.add(-1.0, H_ijk);

      if (k >= 3)
        {
          auto J_ijk = calcJ_ijk<t>(i, j, k);
          W_yiujsk.add(-1.0, J_ijk);
        }

      W<t>().insert(std::make_unique<typename ctraits<t>::Ttensor>(W_yiujsk));

      Gstack<t>().multAndAdd(i + j + k, W<t>(), *Wrond_yiujsk_ptr);
    }
}

/* Here we solve [F_{σⁱ}]+[Dᵢ]+[Eᵢ]=0 to recover g_σⁱ. First we calculate
   conditional Wrond_σⁱ (missing dimension and l=i), then we calculate
   conditional F_σⁱ. Before we can calculate Dᵢ and Eᵢ, we have to obtain
   G_u′ᵐσⁱ⁻ᵐ for m=1,…,i. Then adding Dᵢ and Eᵢ we have the right hand side. We
   solve by S⁻¹ multiplication and update G_σⁱ by calling multAndAdd() for
   dimension l=1.

   Recall that the solved equation here is:
 ny(ny), nu(nu), F_u′ʲσⁱ⁻ʲ.

   Provides: g_σⁱ, G_σⁱ, and G_u′ᵐσⁱ⁻ᵐ for m=1,…,i
*/
template<Storage t>
void
KOrderWelfare::recover_s(int i)
{
  Symmetry sym {0, 0, 0, i};
  JournalRecordPair pa(journal);
  pa << "Conditional welfare: recovering symmetry " << sym << endrec;

  fillWrond<t>(0, 0, i);

  if (is_even(i))
    {
      auto Wrond_si = faaDiBrunoW<t>(sym);
      auto Wrond_si_ptr = Wrond_si.get();
      Wrond<t>().insert(std::move(Wrond_si));

      auto W_si = U<t>().get(sym);
      W_si.add(discount_factor, *Wrond_si_ptr);

      {
        auto H_i = calcH_k<t>(i);
        W_si.add(-1.0, H_i);
      }

      if (i >= 3)
        {
          auto J_i = calcJ_k<t>(i);
          W_si.add(-1.0, J_i);
        }

      W_si.mult(1 / (1 - discount_factor));

      W<t>().insert(std::make_unique<typename ctraits<t>::Ttensor>(W_si));

      Gstack<t>().multAndAdd(i, W<t>(), *Wrond_si_ptr);
    }
}

/* Here we calculate and insert Wrond_yⁱuʲu′ᵐσᵏ⁻ᵐ for m=1,…,k. The derivatives are inserted only for
 * k−m being even. */
template<Storage t>
void
KOrderWelfare::fillWrond(int i, int j, int k)
{
  for (int m = 1; m <= k; m++)
    if (is_even(k - m))
      {
        auto Wrond_yiujupms = faaDiBrunoW<t>(Symmetry {i, j, m, k - m});
        Wrond<t>().insert(std::move(Wrond_yiujupms));
      }
}

/* Here we calculate:

    [Hᵢⱼₖ]_α₁…αᵢβ₁…βⱼ = [F_yⁱuʲu′ᵏ]_α₁…αᵢβ₁…βⱼγ₁…γₖ [Σ]^γ₁…γₖ
                      = -β [Wrond_yⁱuʲu′ᵏ]_α₁…αᵢβ₁…βⱼγ₁…γₖ [Σ]^γ₁…γₖ
   So it is non zero only for even k. */
template<Storage t>
typename ctraits<t>::Ttensor
KOrderWelfare::calcH_ijk(int i, int j, int k) const
{
  typename ctraits<t>::Ttensor res(1, TensorDimens(Symmetry {i, j, 0, 0}, nvs));
  res.zeros();
  if (is_even(k))
    {
      auto tmp = faaDiBrunoW<t>(Symmetry {i, j, k, 0});
      tmp->contractAndAdd(2, res, m<t>().get(Symmetry {k}));
    }
  res.mult(-discount_factor);
  return res;
}

/* Here we calculate
                        ₖ₋₁ ⎛k⎞
    [Jᵢⱼₖ]_α₁…αᵢβ₁…βⱼ =  ∑  ⎝m⎠ [F_yⁱuʲu′ᵐσᵏ⁻ᵐ]_α₁…αᵢβ₁…βⱼγ₁…γₘ [Σ]^γ₁…γₘ
                        ᵐ⁼¹
                           ₖ₋₁ ⎛k⎞
                      = -β  ∑  ⎝m⎠ [Wrond_yⁱuʲu′ᵐσᵏ⁻ᵐ]_α₁…αᵢβ₁…βⱼγ₁…γₘ [Σ]^γ₁…γₘ
                           ᵐ⁼¹
The sum can sum only for even m. */
template<Storage t>
typename ctraits<t>::Ttensor
KOrderWelfare::calcJ_ijk(int i, int j, int k) const
{
  typename ctraits<t>::Ttensor res(1, TensorDimens(Symmetry {i, j, 0, 0}, nvs));
  res.zeros();
  for (int n = 2; n <= k - 1; n += 2)
    {
      auto tmp = faaDiBrunoW<t>(Symmetry {i, j, n, k - n});
      tmp->mult(static_cast<double>(PascalTriangle::noverk(k, n)));
      tmp->contractAndAdd(2, res, m<t>().get(Symmetry {n}));
    }
  res.mult(-discount_factor);
  return res;
}

template<Storage t>
typename ctraits<t>::Ttensor
KOrderWelfare::calcH_ik(int i, int k) const
{
  return calcH_ijk<t>(i, 0, k);
}

template<Storage t>
typename ctraits<t>::Ttensor
KOrderWelfare::calcH_k(int k) const
{
  return calcH_ijk<t>(0, 0, k);
}

template<Storage t>
typename ctraits<t>::Ttensor
KOrderWelfare::calcJ_ik(int i, int k) const
{
  return calcJ_ijk<t>(i, 0, k);
}

template<Storage t>
typename ctraits<t>::Ttensor
KOrderWelfare::calcJ_k(int k) const
{
  return calcJ_ijk<t>(0, 0, k);
}

/* Here we check for residuals of all the solved equations at the given order.
   The method returns the largest residual size. Each check simply evaluates
   the equation. */

template<Storage t>
double
KOrderWelfare::check(int dim) const
{
  KORD_RAISE_IF(dim > W<t>().getMaxDim(), "Wrong dimension for KOrderWelfare::check");
  JournalRecordPair pa(journal);
  pa << "Checking residuals for order = " << dim << endrec;

  double maxerror = 0.0;

  // Check for F_yⁱuʲ=0
  for (int i = 0; i <= dim; i++)
    {
      Symmetry sym {dim - i, i, 0, 0};
      auto r = W<t>().get(sym);
      r.add(-1.0, U<t>().get(sym));
      r.add(-discount_factor, Wrond<t>().get(sym));
      double err = r.getData().getMax();
      JournalRecord(journal) << "\terror for symmetry " << sym << "\tis " << err << endrec;
      maxerror = std::max(err, maxerror);
    }

  // Check for F_yⁱuʲu′ᵏ+Hᵢⱼₖ+Jᵢⱼₖ=0
  for (auto& si : SymmetrySet(dim, 3))
    {
      int i = si[0];
      int j = si[1];
      int k = si[2];
      if (i + j > 0 && k > 0 && is_even(k))
        {
          Symmetry sym {i, j, 0, k};
          auto r = W<t>().get(sym);
          r.add(-1.0, U<t>().get(sym));
          r.add(-discount_factor, Wrond<t>().get(sym));
          auto H_ijk = calcH_ijk<t>(i, j, k);
          r.add(1.0, H_ijk);
          auto J_ijk = calcJ_ijk<t>(i, j, k);
          r.add(1.0, J_ijk);
          double err = r.getData().getMax();
          JournalRecord(journal) << "\terror for symmetry " << sym << "\tis " << err << endrec;
          maxerror = std::max(err, maxerror);
        }
    }

  // Check for F_σⁱ+Dᵢ+Eᵢ=0
  if (is_even(dim))
    {
      Symmetry sym {0, 0, 0, dim};
      auto r = W<t>().get(sym);
      r.add(-1.0, U<t>().get(sym));
      r.add(-discount_factor, Wrond<t>().get(sym));
      auto H_k = calcH_k<t>(dim);
      r.add(1.0, H_k);
      auto J_k = calcJ_k<t>(dim);
      r.add(1.0, J_k);
      double err = r.getData().getMax();
      JournalRecord(journal) << "\terror for symmetry " << sym << "\tis " << err << endrec;
      maxerror = std::max(err, maxerror);
    }

  return maxerror;
}

#endif
