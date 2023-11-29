/*
 * Copyright © 2004 Ondra Kamenik
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

#include "korder.hh"
#include "kord_exception.hh"

/* Here we set ‘ipiv’ and ‘inv’ members of the PLUMatrix depending on its
   content. It is assumed that subclasses will call this method at the end of
   their constructors. */

void
PLUMatrix::calcPLU()
{
  lapack_int info;
  lapack_int rows = nrows(), lda = ld;
  inv = const_cast<const Vector&>(getData());
  dgetrf(&rows, &rows, inv.base(), &lda, ipiv.data(), &info);
}

/* Here we just call the LAPACK machinery to multiply by the inverse. */

void
PLUMatrix::multInv(TwoDMatrix& m) const
{
  KORD_RAISE_IF(m.nrows() != ncols(), "The matrix is not square in PLUMatrix::multInv");
  lapack_int info;
  lapack_int lda = ld;
  lapack_int mcols = m.ncols();
  lapack_int mrows = m.nrows();
  lapack_int ldb = m.getLD();
  dgetrs("N", &mrows, &mcols, inv.base(), &lda, ipiv.data(), m.getData().base(), &ldb, &info);
  KORD_RAISE_IF(info != 0, "Info!=0 in PLUMatrix::multInv");
}

/* Here we construct the matrix A. Its dimension is ‘ny’, and it is

    A=[f_y] + [0 [f_y**₊]·[g**_y*] 0],

   where the first zero spans ‘nstat’ columns, and last zero spans ‘nforw’
   columns. */

MatrixA::MatrixA(const FSSparseTensor& f, const IntSequence& ss, const TwoDMatrix& gy,
                 const PartitionY& ypart) :
    PLUMatrix(ypart.ny())
{
  zeros();

  IntSequence c {1};
  FGSTensor f_y(f, ss, c, TensorDimens(ss, c));
  add(1.0, f_y);

  ConstTwoDMatrix gss_ys(ypart.nstat + ypart.npred, ypart.nyss(), gy);
  c[0] = 0;
  FGSTensor f_yss(f, ss, c, TensorDimens(ss, c));
  TwoDMatrix sub(*this, ypart.nstat, ypart.nys());
  sub.multAndAdd(ConstTwoDMatrix(f_yss), gss_ys);

  calcPLU();
}

/* Here we construct the matrix S. Its dimension is ‘ny’, and it is

    S = [f_y] + [0 [f_y**₊]·[g**_y*] 0] + [0 0 [f_y**₊]]

   It is, in fact, the matrix A plus the third summand. The first zero in the
   summand spans ‘nstat’ columns, the second zero spans ‘npred’ columns. */

MatrixS::MatrixS(const FSSparseTensor& f, const IntSequence& ss, const TwoDMatrix& gy,
                 const PartitionY& ypart) :
    PLUMatrix(ypart.ny())
{
  zeros();

  IntSequence c {1};
  FGSTensor f_y(f, ss, c, TensorDimens(ss, c));
  add(1.0, f_y);

  ConstTwoDMatrix gss_ys(ypart.nstat + ypart.npred, ypart.nyss(), gy);
  c[0] = 0;
  FGSTensor f_yss(f, ss, c, TensorDimens(ss, c));
  TwoDMatrix sub(*this, ypart.nstat, ypart.nys());
  sub.multAndAdd(ConstTwoDMatrix(f_yss), gss_ys);

  TwoDMatrix sub2(*this, ypart.nstat + ypart.npred, ypart.nyss());
  sub2.add(1.0, f_yss);

  calcPLU();
}

// KOrder member access method specializations
/* These are the specializations of container access methods. Nothing
   interesting here. */

template<>
ctraits<Storage::unfold>::Tg&
KOrder::g<Storage::unfold>()
{
  return _ug;
}
template<>
const ctraits<Storage::unfold>::Tg&
KOrder::g<Storage::unfold>() const
{
  return _ug;
}
template<>
ctraits<Storage::fold>::Tg&
KOrder::g<Storage::fold>()
{
  return _fg;
}
template<>
const ctraits<Storage::fold>::Tg&
KOrder::g<Storage::fold>() const
{
  return _fg;
}
template<>
ctraits<Storage::unfold>::Tgs&
KOrder::gs<Storage::unfold>()
{
  return _ugs;
}
template<>
const ctraits<Storage::unfold>::Tgs&
KOrder::gs<Storage::unfold>() const
{
  return _ugs;
}
template<>
ctraits<Storage::fold>::Tgs&
KOrder::gs<Storage::fold>()
{
  return _fgs;
}
template<>
const ctraits<Storage::fold>::Tgs&
KOrder::gs<Storage::fold>() const
{
  return _fgs;
}
template<>
ctraits<Storage::unfold>::Tgss&
KOrder::gss<Storage::unfold>()
{
  return _ugss;
}
template<>
const ctraits<Storage::unfold>::Tgss&
KOrder::gss<Storage::unfold>() const
{
  return _ugss;
}
template<>
ctraits<Storage::fold>::Tgss&
KOrder::gss<Storage::fold>()
{
  return _fgss;
}
template<>
const ctraits<Storage::fold>::Tgss&
KOrder::gss<Storage::fold>() const
{
  return _fgss;
}
template<>
ctraits<Storage::unfold>::TG&
KOrder::G<Storage::unfold>()
{
  return _uG;
}
template<>
const ctraits<Storage::unfold>::TG&
KOrder::G<Storage::unfold>() const
{
  return _uG;
}
template<>
ctraits<Storage::fold>::TG&
KOrder::G<Storage::fold>()
{
  return _fG;
}
template<>
const ctraits<Storage::fold>::TG&
KOrder::G<Storage::fold>() const
{
  return _fG;
}
template<>
ctraits<Storage::unfold>::TZstack&
KOrder::Zstack<Storage::unfold>()
{
  return _uZstack;
}
template<>
const ctraits<Storage::unfold>::TZstack&
KOrder::Zstack<Storage::unfold>() const
{
  return _uZstack;
}
template<>
ctraits<Storage::fold>::TZstack&
KOrder::Zstack<Storage::fold>()
{
  return _fZstack;
}
template<>
const ctraits<Storage::fold>::TZstack&
KOrder::Zstack<Storage::fold>() const
{
  return _fZstack;
}
template<>
ctraits<Storage::unfold>::TGstack&
KOrder::Gstack<Storage::unfold>()
{
  return _uGstack;
}
template<>
const ctraits<Storage::unfold>::TGstack&
KOrder::Gstack<Storage::unfold>() const
{
  return _uGstack;
}
template<>
ctraits<Storage::fold>::TGstack&
KOrder::Gstack<Storage::fold>()
{
  return _fGstack;
}
template<>
const ctraits<Storage::fold>::TGstack&
KOrder::Gstack<Storage::fold>() const
{
  return _fGstack;
}
template<>
ctraits<Storage::unfold>::Tm&
KOrder::m<Storage::unfold>()
{
  return _um;
}
template<>
const ctraits<Storage::unfold>::Tm&
KOrder::m<Storage::unfold>() const
{
  return _um;
}
template<>
ctraits<Storage::fold>::Tm&
KOrder::m<Storage::fold>()
{
  return _fm;
}
template<>
const ctraits<Storage::fold>::Tm&
KOrder::m<Storage::fold>() const
{
  return _fm;
}

/* Here is the constructor of the KOrder class. We pass what we have to. The
   partitioning of the y vector, a sparse container with model derivatives,
   then the first order approximation, these are g_y and gᵤ matrices, and
   covariance matrix of exogenous shocks ‘v’.

   We build the members, it is nothing difficult. Note that we do not make a
   physical copy of sparse tensors, so during running the class, the outer
   world must not change them.

   In the body, we have to set ‘nvs’ array, and initialize g and G containers
   to comply to preconditions of performStep(). */

KOrder::KOrder(int num_stat, int num_pred, int num_both, int num_forw,
               const TensorContainer<FSSparseTensor>& fcont, const TwoDMatrix& gy,
               const TwoDMatrix& gu, const TwoDMatrix& v, Journal& jr) :
    ypart(num_stat, num_pred, num_both, num_forw),
    ny(ypart.ny()), nu(gu.ncols()), maxk(fcont.getMaxDim()), nvs {ypart.nys(), nu, nu, 1}, _ug(4),
    _fg(4), _ugs(4), _fgs(4), _ugss(4), _fgss(4), _uG(4), _fG(4),
    _uZstack(&_uG, ypart.nyss(), &_ug, ny, ypart.nys(), nu),
    _fZstack(&_fG, ypart.nyss(), &_fg, ny, ypart.nys(), nu), _uGstack(&_ugs, ypart.nys(), nu),
    _fGstack(&_fgs, ypart.nys(), nu), _um(maxk, v), _fm(_um), f(fcont),
    matA(f.get(Symmetry {1}), _uZstack.getStackSizes(), gy, ypart),
    matS(f.get(Symmetry {1}), _uZstack.getStackSizes(), gy, ypart),
    matB(f.get(Symmetry {1}), _uZstack.getStackSizes()), journal(jr)
{
  KORD_RAISE_IF(gy.ncols() != ypart.nys(), "Wrong number of columns in gy in KOrder constructor");
  KORD_RAISE_IF(v.ncols() != nu, "Wrong number of columns of Vcov in KOrder constructor");
  KORD_RAISE_IF(nu != v.nrows(), "Wrong number of rows of Vcov in KOrder constructor");
  KORD_RAISE_IF(maxk < 2, "Order of approximation must be at least 2 in KOrder constructor");
  KORD_RAISE_IF(gy.nrows() != ypart.ny(), "Wrong number of rows in gy in KOrder constructor");
  KORD_RAISE_IF(gu.nrows() != ypart.ny(), "Wrong number of rows in gu in KOrder constructor");
  KORD_RAISE_IF(gu.ncols() != nu, "Wrong number of columns in gu in KOrder constructor");

  // Put g_y and gᵤ in the container
  /* Note that g_σ is zero by construction and we do not insert it to the
     container. We insert a new physical copies. */
  auto tgy = std::make_unique<UGSTensor>(ny, TensorDimens(Symmetry {1, 0, 0, 0}, nvs));
  tgy->getData() = gy.getData();
  insertDerivative<Storage::unfold>(std::move(tgy));
  auto tgu = std::make_unique<UGSTensor>(ny, TensorDimens(Symmetry {0, 1, 0, 0}, nvs));
  tgu->getData() = gu.getData();
  insertDerivative<Storage::unfold>(std::move(tgu));

  // Put G_y, G_u and G_u′ in the container
  /* Also note that since g_σ is zero, so is G_σ. */
  auto tGy = faaDiBrunoG<Storage::unfold>(Symmetry {1, 0, 0, 0});
  G<Storage::unfold>().insert(std::move(tGy));
  auto tGu = faaDiBrunoG<Storage::unfold>(Symmetry {0, 1, 0, 0});
  G<Storage::unfold>().insert(std::move(tGu));
  auto tGup = faaDiBrunoG<Storage::unfold>(Symmetry {0, 0, 1, 0});
  G<Storage::unfold>().insert(std::move(tGup));
}

// KOrder::sylvesterSolve() unfolded specialization
/* Here we have an unfolded specialization of sylvesterSolve(). We simply
   create the sylvester object and solve it. Note that the g*_y is not
   continuous in memory as assumed by the sylvester code, so we make a
   temporary copy and pass it as matrix C.

   If the B matrix is empty, in other words there are now forward looking
   variables, then the system becomes AX=D which is solved by simple
   matA.multInv().

   If one wants to display the diagnostic messages from the Sylvester module,
   then after the sylv.solve() one needs to call sylv.getParams().print(""). */

template<>
void
KOrder::sylvesterSolve<Storage::unfold>(ctraits<Storage::unfold>::Ttensor& der) const
{
  JournalRecordPair pa(journal);
  pa << "Sylvester equation for dimension = " << der.getSym()[0] << endrec;
  if (ypart.nys() > 0 && ypart.nyss() > 0)
    {
      KORD_RAISE_IF(!der.isFinite(), "RHS of Sylverster is not finite");
      TwoDMatrix gs_y(gs<Storage::unfold>().get(Symmetry {1, 0, 0, 0}));
      GeneralSylvester sylv(der.getSym()[0], ny, ypart.nys(), ypart.nstat + ypart.npred,
                            matA.getData(), matB.getData(), gs_y.getData(), der.getData());
      sylv.solve();
    }
  else if (ypart.nys() > 0 && ypart.nyss() == 0)
    matA.multInv(der);
}

// KOrder::sylvesterSolve() folded specialization
/* Here is the folded specialization of sylvester. We unfold the right hand
   side. Then we solve it by the unfolded version of sylvesterSolve(), and fold
   it back and copy to output vector. */

template<>
void
KOrder::sylvesterSolve<Storage::fold>(ctraits<Storage::fold>::Ttensor& der) const
{
  ctraits<Storage::unfold>::Ttensor tmp(der);
  sylvesterSolve<Storage::unfold>(tmp);
  ctraits<Storage::fold>::Ttensor ftmp(tmp);
  der.getData() = const_cast<const Vector&>(ftmp.getData());
}

void
KOrder::switchToFolded()
{
  JournalRecordPair pa(journal);
  pa << "Switching from unfolded to folded" << endrec;

  int maxdim = g<Storage::unfold>().getMaxDim();
  for (int dim = 1; dim <= maxdim; dim++)
    for (auto& si : SymmetrySet(dim, 4))
      {
        if (si[2] == 0 && g<Storage::unfold>().check(si))
          {
            auto ft = std::make_unique<FGSTensor>(g<Storage::unfold>().get(si));
            insertDerivative<Storage::fold>(std::move(ft));
            if (dim > 1)
              {
                gss<Storage::unfold>().remove(si);
                gs<Storage::unfold>().remove(si);
                g<Storage::unfold>().remove(si);
              }
          }
        if (G<Storage::unfold>().check(si))
          {
            auto ft = std::make_unique<FGSTensor>(G<Storage::unfold>().get(si));
            G<Storage::fold>().insert(std::move(ft));
            if (dim > 1)
              G<Storage::fold>().remove(si);
          }
      }
}
