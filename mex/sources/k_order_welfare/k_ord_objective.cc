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

#include "k_ord_objective.hh"

#include <cassert>
#include <utility>

KordwDynare::KordwDynare(KordpDynare& m, Journal& jr, Vector& inParams,
                         std::unique_ptr<ObjectiveMFile> objectiveFile_arg,
                         const std::vector<int>& dr_order) :
    model {m},
    journal {jr},
    params {inParams},
    resid(1),
    ud {1},
    objectiveFile {std::move(objectiveFile_arg)}
{
  dynToDynpp.resize(model.ny());
  for (int i = 0; i < model.ny(); i++)
    dynToDynpp[dr_order[i]] = i;
}

void
KordwDynare::calcDerivativesAtSteady()
{

  assert(ud.begin() == ud.end());

  Vector xx(model.nexog());
  xx.zeros();
  resid.zeros();
  objectiveFile->eval(model.getSteady(), xx, params, resid, dynToDynpp, ud);
}

template<>
ctraits<Storage::unfold>::Tg&
KOrderWelfare::g<Storage::unfold>()
{
  return _ug;
}
template<>
const ctraits<Storage::unfold>::Tg&
KOrderWelfare::g<Storage::unfold>() const
{
  return _ug;
}
template<>
ctraits<Storage::fold>::Tg&
KOrderWelfare::g<Storage::fold>()
{
  return _fg;
}
template<>
const ctraits<Storage::fold>::Tg&
KOrderWelfare::g<Storage::fold>() const
{
  return _fg;
}
template<>
ctraits<Storage::unfold>::Tgs&
KOrderWelfare::gs<Storage::unfold>()
{
  return _ugs;
}
template<>
const ctraits<Storage::unfold>::Tgs&
KOrderWelfare::gs<Storage::unfold>() const
{
  return _ugs;
}
template<>
ctraits<Storage::fold>::Tgs&
KOrderWelfare::gs<Storage::fold>()
{
  return _fgs;
}
template<>
const ctraits<Storage::fold>::Tgs&
KOrderWelfare::gs<Storage::fold>() const
{
  return _fgs;
}

template<>
ctraits<Storage::unfold>::TU&
KOrderWelfare::U<Storage::unfold>()
{
  return _uU;
}
template<>
const ctraits<Storage::unfold>::TU&
KOrderWelfare::U<Storage::unfold>() const
{
  return _uU;
}
template<>
ctraits<Storage::fold>::TU&
KOrderWelfare::U<Storage::fold>()
{
  return _fU;
}
template<>
const ctraits<Storage::fold>::TU&
KOrderWelfare::U<Storage::fold>() const
{
  return _fU;
}
template<>
ctraits<Storage::unfold>::TW&
KOrderWelfare::W<Storage::unfold>()
{
  return _uW;
}
template<>
const ctraits<Storage::unfold>::TW&
KOrderWelfare::W<Storage::unfold>() const
{
  return _uW;
}
template<>
ctraits<Storage::fold>::TW&
KOrderWelfare::W<Storage::fold>()
{
  return _fW;
}
template<>
const ctraits<Storage::fold>::TW&
KOrderWelfare::W<Storage::fold>() const
{
  return _fW;
}
template<>
ctraits<Storage::unfold>::TWrond&
KOrderWelfare::Wrond<Storage::unfold>()
{
  return _uWrond;
}
template<>
const ctraits<Storage::unfold>::TWrond&
KOrderWelfare::Wrond<Storage::unfold>() const
{
  return _uWrond;
}
template<>
ctraits<Storage::fold>::TWrond&
KOrderWelfare::Wrond<Storage::fold>()
{
  return _fWrond;
}
template<>
const ctraits<Storage::fold>::TWrond&
KOrderWelfare::Wrond<Storage::fold>() const
{
  return _fWrond;
}
template<>
ctraits<Storage::unfold>::TGstack&
KOrderWelfare::Gstack<Storage::unfold>()
{
  return _uGstack;
}
template<>
const ctraits<Storage::unfold>::TGstack&
KOrderWelfare::Gstack<Storage::unfold>() const
{
  return _uGstack;
}
template<>
ctraits<Storage::fold>::TGstack&
KOrderWelfare::Gstack<Storage::fold>()
{
  return _fGstack;
}
template<>
const ctraits<Storage::fold>::TGstack&
KOrderWelfare::Gstack<Storage::fold>() const
{
  return _fGstack;
}
template<>
ctraits<Storage::unfold>::TXstack&
KOrderWelfare::Xstack<Storage::unfold>()
{
  return _uXstack;
}
template<>
const ctraits<Storage::unfold>::TXstack&
KOrderWelfare::Xstack<Storage::unfold>() const
{
  return _uXstack;
}
template<>
ctraits<Storage::fold>::TXstack&
KOrderWelfare::Xstack<Storage::fold>()
{
  return _fXstack;
}
template<>
const ctraits<Storage::fold>::TXstack&
KOrderWelfare::Xstack<Storage::fold>() const
{
  return _fXstack;
}
template<>
ctraits<Storage::unfold>::Tm&
KOrderWelfare::m<Storage::unfold>()
{
  return _um;
}
template<>
const ctraits<Storage::unfold>::Tm&
KOrderWelfare::m<Storage::unfold>() const
{
  return _um;
}
template<>
ctraits<Storage::fold>::Tm&
KOrderWelfare::m<Storage::fold>()
{
  return _fm;
}
template<>
const ctraits<Storage::fold>::Tm&
KOrderWelfare::m<Storage::fold>() const
{
  return _fm;
}

KOrderWelfare::KOrderWelfare(int num_stat, int num_pred, int num_both, int num_forw, int nu,
                             int ord, double discount_factor,
                             const TensorContainer<FSSparseTensor>& ucont, FGSContainer g_arg,
                             FGSContainer gs_arg, const TwoDMatrix& v, Journal& jr) :
    ypart(num_stat, num_pred, num_both, num_forw),
    ny(ypart.ny()),
    nu(nu),
    maxk(ucont.getMaxDim()),
    order(ord),
    discount_factor {discount_factor},
    nvs {ypart.nys(), nu, nu, 1},
    _uU(4),
    _fU(4),
    _uW(4),
    _fW(4),
    _uWrond(4),
    _fWrond(4),
    _ug(4),
    _fg(std::move(g_arg)),
    _ugs(4),
    _fgs(std::move(gs_arg)),
    _uXstack(&_ug, ny),
    _fXstack(&_fg, ny),
    _uGstack(&_ugs, ypart.nys(), nu),
    _fGstack(&_fgs, ypart.nys(), nu),
    _um(maxk, v),
    _fm(_um),
    u(ucont),
    journal(jr)
{
  KORD_RAISE_IF(v.ncols() != nu, "Wrong number of columns of Vcov in KOrderWelfare constructor");
  KORD_RAISE_IF(nu != v.nrows(), "Wrong number of rows of Vcov in KOrderWelfare constructor");
  for (int ord = 1; ord <= order; ord++)
    {
      JournalRecordPair pa(journal);
      pa << "Unconditional welfare : performing step for order = " << ord << "\n" << endrec;
      for (int j = 0; j <= ord; j++)
        for (int i = 0; i <= j; i++)
          {
            Symmetry sym {ord - j, i, 0, j - i};
            pa << "Recovering symmetry " << sym << "\n" << endrec;
            auto U_sym = faaDiBrunoU<Storage::fold>(sym);
            U<Storage::fold>().insert(std::move(U_sym));
          }
    }
  U<Storage::unfold>() = UGSContainer(U<Storage::fold>());
  g<Storage::unfold>() = UGSContainer(g<Storage::fold>());
  gs<Storage::unfold>() = UGSContainer(gs<Storage::fold>());
}

// KOrderWelfare::sylvesterSolve() unfolded specialization
/* Here we have an unfolded specialization of sylvesterSolve(). We simply
   create the sylvester object and solve it. Note that the W_y is not
   continuous in memory as assumed by the sylvester code, so we make a
   temporary copy and pass it as matrix C.

   If the B matrix is empty, in other words there are now forward looking
   variables, then the system becomes AX=D which is solved by simple
   matA.multInv().

   If one wants to display the diagnostic messages from the Sylvester module,
   then after the sylv.solve() one needs to call sylv.getParams().print(""). */

template<>
void
KOrderWelfare::sylvesterSolve<Storage::unfold>(ctraits<Storage::unfold>::Ttensor& der) const
{
  JournalRecordPair pa(journal);
  pa << "Sylvester equation for dimension = " << der.getSym()[0] << endrec;
  KORD_RAISE_IF(!der.isFinite(), "RHS of Sylverster is not finite");
  TwoDMatrix gs_y(gs<Storage::unfold>().get(Symmetry {1, 0, 0, 0}));
  TwoDMatrix A(1, 1);
  A.unit();
  TwoDMatrix B(1, 1);
  B.unit();
  B.mult(-discount_factor);
  GeneralSylvester sylv(der.getSym()[0], 1, ypart.nys(), 0, A.getData(), B.getData(),
                        gs_y.getData(), der.getData());
  sylv.solve();
}

// KOrder::sylvesterSolve() folded specialization
/* Here is the folded specialization of sylvester. We unfold the right hand
   side. Then we solve it by the unfolded version of sylvesterSolve(), and fold
   it back and copy to output vector. */

template<>
void
KOrderWelfare::sylvesterSolve<Storage::fold>(ctraits<Storage::fold>::Ttensor& der) const
{
  ctraits<Storage::unfold>::Ttensor tmp(der);
  sylvesterSolve<Storage::unfold>(tmp);
  ctraits<Storage::fold>::Ttensor ftmp(tmp);
  der.getData() = const_cast<const Vector&>(ftmp.getData());
}
