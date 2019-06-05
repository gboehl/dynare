// Copyright 2004, Ondra Kamenik

// First order at deterministic steady state

#ifndef FIRST_ORDER_H
#define FIRST_ORDER_H

#include "korder.hh"

#include <mutex>

template<Storage>
class FirstOrderDerivs;
class FirstOrder
{
  template <Storage>
  friend class FirstOrderDerivs;
  PartitionY ypart;
  int nu;
  TwoDMatrix gy;
  TwoDMatrix gu;
  bool bk_cond;
  double b_error;
  int sdim;
  Vector alphar;
  Vector alphai;
  Vector beta;
  double qz_criterium;
  Journal &journal;

  // Passed to LAPACK's DGGES
  static lapack_int order_eigs(const double *alphar, const double *alphai, const double *beta);

  // The value of qz_criterium_global used by the order_eigs function
  /* NB: we have no choice but to use a global variable, since LAPACK won't
     take a closure */
  static double qz_criterium_global;

  // Protects the static qz_criterium_global
  static std::mutex mut;
public:
  FirstOrder(int num_stat, int num_pred, int num_both, int num_forw,
             int num_u, const FSSparseTensor &f, Journal &jr, double qz_crit)
    : ypart(num_stat, num_pred, num_both, num_forw),
      nu(num_u),
      gy(ypart.ny(), ypart.nys()),
      gu(ypart.ny(), nu),
      alphar(ypart.ny()+ypart.nboth),
      alphai(ypart.ny()+ypart.nboth),
      beta(ypart.ny()+ypart.nboth),
      qz_criterium(qz_crit),
      journal(jr)
  {
    solve(FFSTensor(f));
  }
  bool
  isStable() const
  {
    return bk_cond;
  }
  const TwoDMatrix &
  getGy() const
  {
    return gy;
  }
  const TwoDMatrix &
  getGu() const
  {
    return gu;
  }
protected:
  void solve(const TwoDMatrix &f);
  void journalEigs();
};

/* This class only converts the derivatives g_y* and gᵤ to a folded or unfolded
   container. */

template <Storage t>
class FirstOrderDerivs : public ctraits<t>::Tg
{
public:
  FirstOrderDerivs(const FirstOrder &fo)
    : ctraits<t>::Tg(4)
  {
    IntSequence nvs{fo.ypart.nys(), fo.nu, fo.nu, 1};
    auto ten = std::make_unique<typename ctraits<t>::Ttensor>(fo.ypart.ny(), TensorDimens(Symmetry{1, 0, 0, 0}, nvs));
    ten->zeros();
    ten->add(1.0, fo.gy);
    this->insert(std::move(ten));
    ten = std::make_unique<typename ctraits<t>::Ttensor>(fo.ypart.ny(), TensorDimens(Symmetry{0, 1, 0, 0}, nvs));
    ten->zeros();
    ten->add(1.0, fo.gu);
    this->insert(std::move(ten));
  }
};

#endif
