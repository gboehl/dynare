// Copyright 2004, Ondra Kamenik

#include "kord_exception.hh"
#include "korder.hh"

PLUMatrix::PLUMatrix(const PLUMatrix &plu)
  : TwoDMatrix(plu), inv(plu.inv), ipiv(new lapack_int[nrows()])
{
  memcpy(ipiv, plu.ipiv, nrows()*sizeof(lapack_int));
}

/* Here we set |ipiv| and |inv| members of the |PLUMatrix| depending on
   its content. It is assumed that subclasses will call this method at
   the end of their constructors. */

void
PLUMatrix::calcPLU()
{
  lapack_int info;
  lapack_int rows = nrows(), lda = ld;
  inv = (const Vector &) getData();
  dgetrf(&rows, &rows, inv.base(), &lda, ipiv, &info);
}

/* Here we just call the LAPACK machinery to multiply by the inverse. */

void
PLUMatrix::multInv(TwoDMatrix &m) const
{
  KORD_RAISE_IF(m.nrows() != ncols(),
                "The matrix is not square in PLUMatrix::multInv");
  lapack_int info;
  lapack_int lda = ld;
  lapack_int mcols = m.ncols();
  lapack_int mrows = m.nrows();
  lapack_int ldb = m.getLD();
  double *mbase = m.getData().base();
  dgetrs("N", &mrows, &mcols, inv.base(), &lda, ipiv,
         mbase, &ldb, &info);
  KORD_RAISE_IF(info != 0,
                "Info!=0 in PLUMatrix::multInv");
}

/* Here we construct the matrix $A$. Its dimension is |ny|, and it is
   $$A=\left[f_{y}\right]+
   \left[0 \left[f_{y^{**}_+}\right]\cdot\left[g^{**}_{y^*}\right] 0\right]$$,
   where the first zero spans |nstat| columns, and last zero spans
   |nforw| columns. */

MatrixA::MatrixA(const FSSparseTensor &f, const IntSequence &ss,
                 const TwoDMatrix &gy, const PartitionY &ypart)
  : PLUMatrix(ypart.ny())
{
  zeros();

  IntSequence c{1};
  FGSTensor f_y(f, ss, c, TensorDimens(ss, c));
  add(1.0, f_y);

  ConstTwoDMatrix gss_ys(ypart.nstat+ypart.npred, ypart.nyss(), gy);
  c[0] = 0;
  FGSTensor f_yss(f, ss, c, TensorDimens(ss, c));
  TwoDMatrix sub(*this, ypart.nstat, ypart.nys());
  sub.multAndAdd(ConstTwoDMatrix(f_yss), gss_ys);

  calcPLU();
}

/* Here we construct the matrix $S$. Its dimension is |ny|, and it is
   $$S=\left[f_{y}\right]+
   \left[0\quad\left[f_{y^{**}_+}\right]\cdot\left[g^{**}_{y^*}\right]\quad
   0\right]+ \left[0\quad 0\quad\left[f_{y^{**}_+}\right]\right]$$
   It is, in fact, the matrix $A$ plus the third summand. The first zero
   in the summand spans |nstat| columns, the second zero spans |npred|
   columns. */

MatrixS::MatrixS(const FSSparseTensor &f, const IntSequence &ss,
                 const TwoDMatrix &gy, const PartitionY &ypart)
  : PLUMatrix(ypart.ny())
{
  zeros();

  IntSequence c{1};
  FGSTensor f_y(f, ss, c, TensorDimens(ss, c));
  add(1.0, f_y);

  ConstTwoDMatrix gss_ys(ypart.nstat+ypart.npred, ypart.nyss(), gy);
  c[0] = 0;
  FGSTensor f_yss(f, ss, c, TensorDimens(ss, c));
  TwoDMatrix sub(*this, ypart.nstat, ypart.nys());
  sub.multAndAdd(ConstTwoDMatrix(f_yss), gss_ys);

  TwoDMatrix sub2(*this, ypart.nstat+ypart.npred, ypart.nyss());
  sub2.add(1.0, f_yss);

  calcPLU();
}

// |KOrder| member access method specializations
/* These are the specializations of container access methods. Nothing
   interesting here. */

template<>
ctraits<KOrder::unfold>::Tg& KOrder::g<KOrder::unfold>()
{
  return _ug;
}
template<>
const ctraits<KOrder::unfold>::Tg& KOrder::g<KOrder::unfold>() const
{ return _ug;}
template<>
ctraits<KOrder::fold>::Tg& KOrder::g<KOrder::fold>()
{
  return _fg;
}
template<>
const ctraits<KOrder::fold>::Tg& KOrder::g<KOrder::fold>() const
{ return _fg;}
template<>
ctraits<KOrder::unfold>::Tgs& KOrder::gs<KOrder::unfold>()
{
  return _ugs;
}
template<>
const ctraits<KOrder::unfold>::Tgs& KOrder::gs<KOrder::unfold>() const
{ return _ugs;}
template<>
ctraits<KOrder::fold>::Tgs& KOrder::gs<KOrder::fold>()
{
  return _fgs;
}
template<>
const ctraits<KOrder::fold>::Tgs& KOrder::gs<KOrder::fold>() const
{ return _fgs;}
template<>
ctraits<KOrder::unfold>::Tgss& KOrder::gss<KOrder::unfold>()
{
  return _ugss;
}
template<>
const ctraits<KOrder::unfold>::Tgss& KOrder::gss<KOrder::unfold>() const
{ return _ugss;}
template<>
ctraits<KOrder::fold>::Tgss& KOrder::gss<KOrder::fold>()
{
  return _fgss;
}
template<>
const ctraits<KOrder::fold>::Tgss& KOrder::gss<KOrder::fold>() const
{ return _fgss;}
template<>
ctraits<KOrder::unfold>::TG& KOrder::G<KOrder::unfold>()
{
  return _uG;
}
template<>
const ctraits<KOrder::unfold>::TG& KOrder::G<KOrder::unfold>() const
{ return _uG;}
template<>
ctraits<KOrder::fold>::TG& KOrder::G<KOrder::fold>()
{
  return _fG;
}
template<>
const ctraits<KOrder::fold>::TG& KOrder::G<KOrder::fold>() const
{ return _fG;}
template<>
ctraits<KOrder::unfold>::TZstack& KOrder::Zstack<KOrder::unfold>()
{
  return _uZstack;
}
template<>
const ctraits<KOrder::unfold>::TZstack& KOrder::Zstack<KOrder::unfold>() const
{ return _uZstack;}
template<>
ctraits<KOrder::fold>::TZstack& KOrder::Zstack<KOrder::fold>()
{
  return _fZstack;
}
template<>
const ctraits<KOrder::fold>::TZstack& KOrder::Zstack<KOrder::fold>() const
{ return _fZstack;}
template<>
ctraits<KOrder::unfold>::TGstack& KOrder::Gstack<KOrder::unfold>()
{
  return _uGstack;
}
template<>
const ctraits<KOrder::unfold>::TGstack& KOrder::Gstack<KOrder::unfold>() const
{ return _uGstack;}
template<>
ctraits<KOrder::fold>::TGstack& KOrder::Gstack<KOrder::fold>()
{
  return _fGstack;
}
template<>
const ctraits<KOrder::fold>::TGstack& KOrder::Gstack<KOrder::fold>() const
{ return _fGstack;}
template<>
ctraits<KOrder::unfold>::Tm& KOrder::m<KOrder::unfold>()
{
  return _um;
}
template<>
const ctraits<KOrder::unfold>::Tm& KOrder::m<KOrder::unfold>() const
{ return _um;}
template<>
ctraits<KOrder::fold>::Tm& KOrder::m<KOrder::fold>()
{
  return _fm;
}
template<>
const ctraits<KOrder::fold>::Tm& KOrder::m<KOrder::fold>() const
{ return _fm;}

/* Here is the constructor of the |KOrder| class. We pass what we have
   to. The partitioning of the $y$ vector, a sparse container with model
   derivatives, then the first order approximation, these are $g_y$ and
   $g_u$ matrices, and covariance matrix of exogenous shocks |v|.

   We build the members, it is nothing difficult. Note that we do not make
   a physical copy of sparse tensors, so during running the class, the
   outer world must not change them.

   In the body, we have to set |nvs| array, and initialize $g$ and $G$
   containers to comply to preconditions of |performStep|. */

KOrder::KOrder(int num_stat, int num_pred, int num_both, int num_forw,
               const TensorContainer<FSSparseTensor> &fcont,
               const TwoDMatrix &gy, const TwoDMatrix &gu, const TwoDMatrix &v,
               Journal &jr)
  : ypart(num_stat, num_pred, num_both, num_forw),
    ny(ypart.ny()), nu(gu.ncols()), maxk(fcont.getMaxDim()),
    nvs{ypart.nys(), nu, nu, 1},
    _ug(4), _fg(4), _ugs(4), _fgs(4), _ugss(4), _fgss(4),
    _uG(4), _fG(4),
    _uZstack(&_uG, ypart.nyss(), &_ug, ny, ypart.nys(), nu),
    _fZstack(&_fG, ypart.nyss(), &_fg, ny, ypart.nys(), nu),
    _uGstack(&_ugs, ypart.nys(), nu),
    _fGstack(&_fgs, ypart.nys(), nu),
    _um(maxk, v), _fm(_um), f(fcont),
    matA(f.get(Symmetry{1}), _uZstack.getStackSizes(), gy, ypart),
    matS(f.get(Symmetry{1}), _uZstack.getStackSizes(), gy, ypart),
    matB(f.get(Symmetry{1}), _uZstack.getStackSizes()),
    journal(jr)
{
  KORD_RAISE_IF(gy.ncols() != ypart.nys(),
                "Wrong number of columns in gy in KOrder constructor");
  KORD_RAISE_IF(v.ncols() != nu,
                "Wrong number of columns of Vcov in KOrder constructor");
  KORD_RAISE_IF(nu != v.nrows(),
                "Wrong number of rows of Vcov in KOrder constructor");
  KORD_RAISE_IF(maxk < 2,
                "Order of approximation must be at least 2 in KOrder constructor");
  KORD_RAISE_IF(gy.nrows() != ypart.ny(),
                "Wrong number of rows in gy in KOrder constructor");
  KORD_RAISE_IF(gu.nrows() != ypart.ny(),
                "Wrong number of rows in gu in KOrder constructor");
  KORD_RAISE_IF(gu.ncols() != nu,
                "Wrong number of columns in gu in KOrder constructor");

  // put $g_y$ and $g_u$ to the container
  /* Note that $g_\sigma$ is zero by the nature and we do not insert it to
     the container. We insert a new physical copies. */
  auto tgy = std::make_unique<UGSTensor>(ny, TensorDimens(Symmetry{1, 0, 0, 0}, nvs));
  tgy->getData() = gy.getData();
  insertDerivative<unfold>(std::move(tgy));
  auto tgu = std::make_unique<UGSTensor>(ny, TensorDimens(Symmetry{0, 1, 0, 0}, nvs));
  tgu->getData() = gu.getData();
  insertDerivative<unfold>(std::move(tgu));

  // put $G_y$, $G_u$ and $G_{u'}$ to the container
  /* Also note that since $g_\sigma$ is zero, so $G_\sigma$. */
  auto tGy = faaDiBrunoG<unfold>(Symmetry{1, 0, 0, 0});
  G<unfold>().insert(std::move(tGy));
  auto tGu = faaDiBrunoG<unfold>(Symmetry{0, 1, 0, 0});
  G<unfold>().insert(std::move(tGu));
  auto tGup = faaDiBrunoG<unfold>(Symmetry{0, 0, 1, 0});
  G<unfold>().insert(std::move(tGup));
}

// |KOrder::sylvesterSolve| unfolded specialization
/* Here we have an unfolded specialization of |sylvesterSolve|. We
   simply create the sylvester object and solve it. Note that the $g^*_y$
   is not continuous in memory as assumed by the sylvester code, so we
   make a temporary copy and pass it as matrix $C$.

   If the $B$ matrix is empty, in other words there are now forward
   looking variables, then the system becomes $AX=D$ which is solved by
   simple |matA.multInv()|.

   If one wants to display the diagnostic messages from the Sylvester
   module, then after the |sylv.solve()| one needs to call
   |sylv.getParams().print("")|. */

template<>
void
KOrder::sylvesterSolve<KOrder::unfold>(ctraits<unfold>::Ttensor &der) const
{
  JournalRecordPair pa(journal);
  pa << "Sylvester equation for dimension = " << der.getSym()[0] << endrec;
  if (ypart.nys() > 0 && ypart.nyss() > 0)
    {
      KORD_RAISE_IF(!der.isFinite(),
                    "RHS of Sylverster is not finite");
      TwoDMatrix gs_y(gs<unfold>().get(Symmetry{1, 0, 0, 0}));
      GeneralSylvester sylv(der.getSym()[0], ny, ypart.nys(),
                            ypart.nstat+ypart.npred,
                            matA.getData(), matB.getData(),
                            gs_y.getData(), der.getData());
      sylv.solve();
    }
  else if (ypart.nys() > 0 && ypart.nyss() == 0)
    {
      matA.multInv(der);
    }
}

// |KOrder::sylvesterSolve| folded specialization
/* Here is the folded specialization of sylvester. We unfold the right
   hand side. Then we solve it by the unfolded version of
   |sylvesterSolve|, and fold it back and copy to output vector. */

template<>
void
KOrder::sylvesterSolve<KOrder::fold>(ctraits<fold>::Ttensor &der) const
{
  ctraits<unfold>::Ttensor tmp(der);
  sylvesterSolve<unfold>(tmp);
  ctraits<fold>::Ttensor ftmp(tmp);
  der.getData() = (const Vector &) (ftmp.getData());
}

void
KOrder::switchToFolded()
{
  JournalRecordPair pa(journal);
  pa << "Switching from unfolded to folded" << endrec;

  int maxdim = g<unfold>().getMaxDim();
  for (int dim = 1; dim <= maxdim; dim++)
    for (auto &si : SymmetrySet(dim, 4))
      {
        if (si[2] == 0 && g<unfold>().check(si))
          {
            auto ft = std::make_unique<FGSTensor>(g<unfold>().get(si));
            insertDerivative<fold>(std::move(ft));
            if (dim > 1)
              {
                gss<unfold>().remove(si);
                gs<unfold>().remove(si);
                g<unfold>().remove(si);
              }
          }
        if (G<unfold>().check(si))
          {
            auto ft = std::make_unique<FGSTensor>(G<unfold>().get(si));
            G<fold>().insert(std::move(ft));
            if (dim > 1)
              {
                G<fold>().remove(si);
              }
          }
      }
}
