/*
 * Copyright © 2005 Ondra Kamenik
 * Copyright © 2019 Dynare Team
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

#include "korder_stoch.hh"

/* Same as MatrixA constructor, but the submatrix ‘gss_ys’ is passed
   directly. */

MatrixAA::MatrixAA(const FSSparseTensor& f, const IntSequence& ss, const TwoDMatrix& gss_ys,
                   const PartitionY& ypart) :
    PLUMatrix(ypart.ny())
{
  zeros();

  IntSequence c {1};
  FGSTensor f_y(f, ss, c, TensorDimens(ss, c));
  add(1.0, f_y);

  c[0] = 0;
  FGSTensor f_yss(f, ss, c, TensorDimens(ss, c));
  TwoDMatrix sub(*this, ypart.nstat, ypart.nys());
  sub.multAndAdd(f_yss, gss_ys);

  calcPLU();
}

// KOrderStoch folded constructor code
KOrderStoch::KOrderStoch(const PartitionY& yp, int nu, const TensorContainer<FSSparseTensor>& fcont,
                         const FGSContainer& hh, Journal& jr) :
    nvs {yp.nys(), nu, nu, 1},
    ypart(yp),
    journal(jr),
    _ug(4),
    _fg(4),
    _ugs(4),
    _fgs(4),
    _uG(4),
    _fG(4),
    _uh(nullptr),
    _fh(&hh),
    _uZstack(&_uG, ypart.nyss(), &_ug, ypart.ny(), ypart.nys(), nu),
    _fZstack(&_fG, ypart.nyss(), &_fg, ypart.ny(), ypart.nys(), nu),
    _uGstack(&_ugs, ypart.nys(), nu),
    _fGstack(&_fgs, ypart.nys(), nu),
    f(fcont),
    matA(fcont.get(Symmetry {1}), _uZstack.getStackSizes(), hh.get(Symmetry {1, 0, 0, 0}), ypart)
{
}

// KOrderStoch unfolded constructor code
KOrderStoch::KOrderStoch(const PartitionY& yp, int nu, const TensorContainer<FSSparseTensor>& fcont,
                         const UGSContainer& hh, Journal& jr) :
    nvs {yp.nys(), nu, nu, 1},
    ypart(yp),
    journal(jr),
    _ug(4),
    _fg(4),
    _ugs(4),
    _fgs(4),
    _uG(4),
    _fG(4),
    _uh(&hh),
    _fh(nullptr),
    _uZstack(&_uG, ypart.nyss(), &_ug, ypart.ny(), ypart.nys(), nu),
    _fZstack(&_fG, ypart.nyss(), &_fg, ypart.ny(), ypart.nys(), nu),
    _uGstack(&_ugs, ypart.nys(), nu),
    _fGstack(&_fgs, ypart.nys(), nu),
    f(fcont),
    matA(fcont.get(Symmetry {1}), _uZstack.getStackSizes(), hh.get(Symmetry {1, 0, 0, 0}), ypart)
{
}

// KOrderStoch convenience method specializations
template<>
ctraits<Storage::unfold>::Tg&
KOrderStoch::g<Storage::unfold>()
{
  return _ug;
}
template<>
const ctraits<Storage::unfold>::Tg&
KOrderStoch::g<Storage::unfold>() const
{
  return _ug;
}
template<>
ctraits<Storage::fold>::Tg&
KOrderStoch::g<Storage::fold>()
{
  return _fg;
}
template<>
const ctraits<Storage::fold>::Tg&
KOrderStoch::g<Storage::fold>() const
{
  return _fg;
}
template<>
ctraits<Storage::unfold>::Tgs&
KOrderStoch::gs<Storage::unfold>()
{
  return _ugs;
}
template<>
const ctraits<Storage::unfold>::Tgs&
KOrderStoch::gs<Storage::unfold>() const
{
  return _ugs;
}
template<>
ctraits<Storage::fold>::Tgs&
KOrderStoch::gs<Storage::fold>()
{
  return _fgs;
}
template<>
const ctraits<Storage::fold>::Tgs&
KOrderStoch::gs<Storage::fold>() const
{
  return _fgs;
}
template<>
const ctraits<Storage::unfold>::Tgss&
KOrderStoch::h<Storage::unfold>() const
{
  return *_uh;
}
template<>
const ctraits<Storage::fold>::Tgss&
KOrderStoch::h<Storage::fold>() const
{
  return *_fh;
}
template<>
ctraits<Storage::unfold>::TG&
KOrderStoch::G<Storage::unfold>()
{
  return _uG;
}
template<>
const ctraits<Storage::unfold>::TG&
KOrderStoch::G<Storage::unfold>() const
{
  return _uG;
}
template<>
ctraits<Storage::fold>::TG&
KOrderStoch::G<Storage::fold>()
{
  return _fG;
}
template<>
const ctraits<Storage::fold>::TG&
KOrderStoch::G<Storage::fold>() const
{
  return _fG;
}
template<>
ctraits<Storage::unfold>::TZXstack&
KOrderStoch::Zstack<Storage::unfold>()
{
  return _uZstack;
}
template<>
const ctraits<Storage::unfold>::TZXstack&
KOrderStoch::Zstack<Storage::unfold>() const
{
  return _uZstack;
}
template<>
ctraits<Storage::fold>::TZXstack&
KOrderStoch::Zstack<Storage::fold>()
{
  return _fZstack;
}
template<>
const ctraits<Storage::fold>::TZXstack&
KOrderStoch::Zstack<Storage::fold>() const
{
  return _fZstack;
}
template<>
ctraits<Storage::unfold>::TGXstack&
KOrderStoch::Gstack<Storage::unfold>()
{
  return _uGstack;
}
template<>
const ctraits<Storage::unfold>::TGXstack&
KOrderStoch::Gstack<Storage::unfold>() const
{
  return _uGstack;
}
template<>
ctraits<Storage::fold>::TGXstack&
KOrderStoch::Gstack<Storage::fold>()
{
  return _fGstack;
}
template<>
const ctraits<Storage::fold>::TGXstack&
KOrderStoch::Gstack<Storage::fold>() const
{
  return _fGstack;
}
