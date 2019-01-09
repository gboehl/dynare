// Copyright 2005, Ondra Kamenik

// Higher order at stochastic steady

/* This file defines a number of classes of which |KOrderStoch| is the
   main purpose. Basically, |KOrderStoch| calculates first and higher
   order Taylor expansion of a policy rule at $\sigma>0$ with explicit
   forward $g^{**}$. More formally, we have to solve a policy rule $g$
   from
   $$E_t\left[f(g^{**}(g^*(y^*_t,u_t,\sigma),u_{t+1},\sigma),g(y^*,u_t,\sigma),y^*,u_t)\right]$$
   As an introduction in {\tt approximation.hweb} argues, $g^{**}$ at
   tine $t+1$ must be given from outside. Let the explicit
   $E_t(g^{**}(y^*,u_{t+1},\sigma)$ be equal to $h(y^*,\sigma)$. Then we
   have to solve
   $$f(h(g^*(y^*,u,\sigma),\sigma),g(y,u,\sigma),y,u),$$
   which is much easier than fully implicit system for $\sigma=0$.

   Besides the class |KOrderStoch|, we declare here also classes for the
   new containers corresponding to
   $f(h(g^*(y^*,u,\sigma),\sigma),g(y,u,\sigma),y,u)$. Further, we
   declare |IntegDerivs| and |StochForwardDerivs| classes which basically
   calculate $h$ as an extrapolation based on an approximation to $g$ at
   lower $\sigma$. */

#include "korder.hh"
#include "faa_di_bruno.hh"
#include "journal.hh"

/* This class is a container, which has a specialized constructor
   integrating the policy rule at given $\sigma$. */

template <int t>
class IntegDerivs : public ctraits<t>::Tgss
{
public:
  IntegDerivs(int r, const IntSequence &nvs, const _Tgss &g, const __Tm &mom,
              double at_sigma);
};

/* This constructor integrates a rule (namely its $g^{**}$ part) with
   respect to $u=\tilde\sigma\eta$, and stores to the object the
   derivatives of this integral $h$ at $(y^*,u,\sigma)=(\tilde
   y^*,0,\tilde\sigma)$. The original container of $g^{**}$, the moments of
   the stochastic shocks |mom| and the $\tilde\sigma$ are input.

   The code follows the following derivation
   \def\lims{\vbox{\baselineskip=0pt\lineskip=1pt
   \setbox0=\hbox{$\scriptstyle n+k=p$}\hbox to\wd0{\hss$\scriptstyle m=0$\hss}\box0}}
   $$
   \eqalign{h(y,\sigma)&=E_t\left[g(y,u',\sigma)\right]=\cr
   &=\tilde y+\sum_{d=1}{1\over d!}\sum_{i+j+k=d}\pmatrix{d\cr i,j,k}\left[g_{y^iu^j\sigma^k}\right]
   (y^*-\tilde y^*)^i\sigma^j\Sigma^j(\sigma-\tilde\sigma)^k\cr
   &=\tilde y+\sum_{d=1}{1\over d!}\sum_{i+m+n+k=d}\pmatrix{d\cr i,m+n,k}
   \left[g_{y^iu^{m+n}\sigma^k}\right]
   \hat y^{*i}\Sigma^{m+n}\pmatrix{m+n\cr m,n}{\tilde\sigma}^m\hat\sigma^{k+n}\cr
   &=\tilde y+\sum_{d=1}{1\over d!}\sum_{i+m+n+k=d}\pmatrix{d\cr i,m,n,k}
   \left[g_{y^iu^{m+n}\sigma^k}\right]
   \Sigma^{m+n}{\tilde\sigma}^m\hat y^{*i}\hat\sigma^{k+n}\cr
   &=\tilde y+\sum_{d=1}{1\over d!}\sum_{i+p=d}\sum_{\lims}\pmatrix{d\cr i,m,n,k}
   \left[g_{y^iu^{m+n}\sigma^k}\right]
   \Sigma^{m+n}{\tilde\sigma}^m\hat y^{*i}\hat\sigma^{k+n}\cr
   &=\tilde y+\sum_{d=1}{1\over d!}\sum_{i+p=d}\pmatrix{d\cr i,p}
   \left[\sum_{\lims}\pmatrix{p\cr n,k}{1\over m!}
   \left[g_{y^iu^{m+n}\sigma^k}\right]
   \Sigma^{m+n}{\tilde\sigma}^m\right]\hat y^{*i}\hat\sigma^{k+n},
   }
   $$
   where $\pmatrix{a\cr b_1,\ldots, b_n}$ is a generalized combination
   number, $p=k+n$, $\hat\sigma=\sigma-\tilde\sigma$, $\hat
   y^*=y^*-\tilde y$, and we dropped writing the multidimensional indexes
   in Einstein summation.

   This implies that
   $$h_{y^i\sigma^p}=\sum_{\lims}\pmatrix{p\cr n,k}{1\over m!}
   \left[g_{y^iu^{m+n}\sigma^k}\right]
   \Sigma^{m+n}{\tilde\sigma}^m$$
   and this is exactly what the code does. */

template <int t>
IntegDerivs<t>::IntegDerivs(int r, const IntSequence &nvs, const _Tgss &g, const __Tm &mom,
                            double at_sigma)
  : ctraits<t>::Tgss(4)
{
  int maxd = g.getMaxDim();
  for (int d = 1; d <= maxd; d++)
    {
      for (int i = 0; i <= d; i++)
        {
          int p = d-i;
          Symmetry sym(i, 0, 0, p);
          _Ttensor *ten = new _Ttensor(r, TensorDimens(sym, nvs));

          // calculate derivative $h_{y^i\sigma^p}$
          /* This code calculates
             $$h_{y^i\sigma^p}=\sum_{\lims}\pmatrix{p\cr n,k}{1\over m!}
             \left[g_{y^iu^{m+n}\sigma^k}\right]
             \Sigma^{m+n}{\tilde\sigma}^m$$
             and stores it in |ten|. */
          ten->zeros();
          for (int n = 0; n <= p; n++)
            {
              int k = p-n;
              int povern = Tensor::noverk(p, n);
              int mfac = 1;
              for (int m = 0; i+m+n+k <= maxd; m++, mfac *= m)
                {
                  double mult = (pow(at_sigma, m)*povern)/mfac;
                  Symmetry sym_mn(i, m+n, 0, k);
                  if (m+n == 0 && g.check(sym_mn))
                    ten->add(mult, *(g.get(sym_mn)));
                  if (m+n > 0 && KOrder::is_even(m+n) && g.check(sym_mn))
                    {
                      _Ttensor gtmp(*(g.get(sym_mn)));
                      gtmp.mult(mult);
                      gtmp.contractAndAdd(1, *ten, *(mom.get(Symmetry(m+n))));
                    }
                }
            }

          this->insert(ten);
        }
    }
}

/* This class calculates an extrapolation of expectation of forward
   derivatives. It is a container, all calculations are done in a
   constructor.

   The class calculates derivatives of $E[g(y*,u,\sigma)]$ at $(\bar
   y^*,\bar\sigma)$. The derivatives are extrapolated based on
   derivatives at $(\tilde y^*,\tilde\sigma)$. */

template <int t>
class StochForwardDerivs : public ctraits<t>::Tgss
{
public:
  StochForwardDerivs(const PartitionY &ypart, int nu,
                     const _Tgss &g, const __Tm &m,
                     const Vector &ydelta, double sdelta,
                     double at_sigma);
};

/* This is the constructor which performs the integration and the
   extrapolation. Its parameters are: |g| is the container of derivatives
   at $(\tilde y,\tilde\sigma)$; |m| are the moments of stochastic
   shocks; |ydelta| is a difference of the steady states $\bar y-\tilde
   y$; |sdelta| is a difference between new sigma and old sigma
   $\bar\sigma-\tilde\sigma$, and |at_sigma| is $\tilde\sigma$. There is
   no need of inputing the $\tilde y$.

   We do the operation in four steps:
   \orderedlist
   \li Integrate $g^{**}$, the derivatives are at $(\tilde y,\tilde\sigma)$
   \li Form the (full symmetric) polynomial from the derivatives stacking
   $\left[\matrix{y^*\cr\sigma}\right]$
   \li Centralize this polynomial about $(\bar y,\bar\sigma)$
   \li Recover general symmetry tensors from the (full symmetric) polynomial
   \endorderedlist */

template <int t>
StochForwardDerivs<t>::StochForwardDerivs(const PartitionY &ypart, int nu,
                                          const _Tgss &g, const __Tm &m,
                                          const Vector &ydelta, double sdelta,
                                          double at_sigma)
  : ctraits<t>::Tgss(4)
{
  int maxd = g.getMaxDim();
  int r = ypart.nyss();

  // make |g_int| be integral of $g^{**}$ at $(\tilde y,\tilde\sigma)$
  /* This simply constructs |IntegDerivs| class. Note that the |nvs| of
     the tensors has zero dimensions for shocks, this is because we need to
     make easily stacks of the form $\left[\matrix{y^*\cr\sigma}\right]$ in
     the next step. */
  IntSequence nvs(4);
  nvs[0] = ypart.nys(); nvs[1] = 0; nvs[2] = 0; nvs[3] = 1;
  IntegDerivs<t> g_int(r, nvs, g, m, at_sigma);

  // make |g_int_sym| be full symmetric polynomial from |g_int|
  /* Here we just form a polynomial whose unique variable corresponds to
     $\left[\matrix{y^*\cr\sigma}\right]$ stack. */
  _Tpol g_int_sym(r, ypart.nys()+1);
  for (int d = 1; d <= maxd; d++)
    {
      auto *ten = new _Ttensym(r, ypart.nys()+1, d);
      ten->zeros();
      for (int i = 0; i <= d; i++)
        {
          int k = d-i;
          if (g_int.check(Symmetry(i, 0, 0, k)))
            ten->addSubTensor(*(g_int.get(Symmetry(i, 0, 0, k))));
        }
      g_int_sym.insert(ten);
    }

  // make |g_int_cent| the centralized polynomial about $(\bar y,\bar\sigma)$
  /* Here we centralize the polynomial to $(\bar y,\bar\sigma)$ knowing
     that the polynomial was centralized about $(\tilde
     y,\tilde\sigma)$. This is done by derivating and evaluating the
     derivated polynomial at $(\bar y-\tilde
     y,\bar\sigma-\tilde\sigma)$. The stack of this vector is |delta| in
     the code. */
  Vector delta(ypart.nys()+1);
  Vector dy(delta, 0, ypart.nys());
  ConstVector dy_in(ydelta, ypart.nstat, ypart.nys());
  dy = dy_in;
  delta[ypart.nys()] = sdelta;
  _Tpol g_int_cent(r, ypart.nys()+1);
  for (int d = 1; d <= maxd; d++)
    {
      g_int_sym.derivative(d-1);
      _Ttensym *der = g_int_sym.evalPartially(d, delta);
      g_int_cent.insert(der);
    }

  // pull out general symmetry tensors from |g_int_cent|
  /* Here we only recover the general symmetry derivatives from the full
     symmetric polynomial. Note that the derivative get the true |nvs|. */
  IntSequence ss(4);
  ss[0] = ypart.nys(); ss[1] = 0; ss[2] = 0; ss[3] = 1;
  IntSequence pp(4);
  pp[0] = 0; pp[1] = 1; pp[2] = 2; pp[3] = 3;
  IntSequence true_nvs(nvs);
  true_nvs[1] = nu; true_nvs[2] = nu;
  for (int d = 1; d <= maxd; d++)
    {
      if (g_int_cent.check(Symmetry(d)))
        {
          for (int i = 0; i <= d; i++)
            {
              Symmetry sym(i, 0, 0, d-i);
              IntSequence coor(sym, pp);
              _Ttensor *ten = new _Ttensor(*(g_int_cent.get(Symmetry(d))), ss, coor,
                                           TensorDimens(sym, true_nvs));
              this->insert(ten);
            }
        }
    }
}

/* This container corresponds to $h(g^*(y,u,\sigma),\sigma)$. Note that
   in our application, the $\sigma$ as a second argument to $h$ will be
   its fourth variable in symmetry, so we have to do four member stack
   having the second and third stack dummy. */

template <class _Ttype>
class GXContainer : public GContainer<_Ttype>
{
public:
  typedef StackContainerInterface<_Ttype> _Stype;
  typedef typename StackContainer<_Ttype>::_Ctype _Ctype;
  typedef typename StackContainer<_Ttype>::itype itype;
  GXContainer(const _Ctype *gs, int ngs, int nu)
    : GContainer<_Ttype>(gs, ngs, nu)
  {
  }
  itype getType(int i, const Symmetry &s) const override;
};

/* This routine corresponds to this stack:
   $$\left[\matrix{g^*(y,u,\sigma)\cr dummy\cr dummy\cr\sigma}\right]$$ */

template <class _Ttype>
typename GXContainer<_Ttype>::itype
GXContainer<_Ttype>::getType(int i, const Symmetry &s) const
{
  if (i == 0)
    if (s[2] > 0)
      return _Stype::zero;
    else
      return _Stype::matrix;
  if (i == 1)
    return _Stype::zero;
  if (i == 2)
    return _Stype::zero;
  if (i == 3)
    if (s == Symmetry(0, 0, 0, 1))
      return _Stype::unit;
    else
      return _Stype::zero;

  KORD_RAISE("Wrong stack index in GXContainer::getType");
}

/* This container corresponds to $f(H(y,u,\sigma),g(y,u,sigma),y,u)$,
   where the $H$ has the size (number of rows) as $g^{**}$. Since it is
   very simmilar to |ZContainer|, we inherit form it and override only
   |getType| method. */

template <class _Ttype>
class ZXContainer : public ZContainer<_Ttype>
{
public:
  typedef StackContainerInterface<_Ttype> _Stype;
  typedef typename StackContainer<_Ttype>::_Ctype _Ctype;
  typedef typename StackContainer<_Ttype>::itype itype;
  ZXContainer(const _Ctype *gss, int ngss, const _Ctype *g, int ng, int ny, int nu)
    : ZContainer<_Ttype>(gss, ngss, g, ng, ny, nu)
  {
  }
  itype getType(int i, const Symmetry &s) const override;
};

/* This |getType| method corresponds to this stack:
   $$\left[\matrix{H(y,u,\sigma)\cr g(y,u,\sigma)\cr y\cr u}\right]$$ */

template <class _Ttype>
typename ZXContainer<_Ttype>::itype
ZXContainer<_Ttype>::getType(int i, const Symmetry &s) const
{
  if (i == 0)
    if (s[2] > 0)
      return _Stype::zero;
    else
      return _Stype::matrix;
  if (i == 1)
    if (s[2] > 0)
      return _Stype::zero;
    else
      return _Stype::matrix;
  if (i == 2)
    if (s == Symmetry(1, 0, 0, 0))
      return _Stype::unit;
    else
      return _Stype::zero;
  if (i == 3)
    if (s == Symmetry(0, 1, 0, 0))
      return _Stype::unit;
    else
      return _Stype::zero;

  KORD_RAISE("Wrong stack index in ZXContainer::getType");
}

class UnfoldedGXContainer : public GXContainer<UGSTensor>, public UnfoldedStackContainer
{
public:
  typedef TensorContainer<UGSTensor> _Ctype;
  UnfoldedGXContainer(const _Ctype *gs, int ngs, int nu)
    : GXContainer<UGSTensor>(gs, ngs, nu)
  {
  }
};

class FoldedGXContainer : public GXContainer<FGSTensor>, public FoldedStackContainer
{
public:
  typedef TensorContainer<FGSTensor> _Ctype;
  FoldedGXContainer(const _Ctype *gs, int ngs, int nu)
    : GXContainer<FGSTensor>(gs, ngs, nu)
  {
  }
};

class UnfoldedZXContainer : public ZXContainer<UGSTensor>, public UnfoldedStackContainer
{
public:
  typedef TensorContainer<UGSTensor> _Ctype;
  UnfoldedZXContainer(const _Ctype *gss, int ngss, const _Ctype *g, int ng, int ny, int nu)
    : ZXContainer<UGSTensor>(gss, ngss, g, ng, ny, nu)
  {
  }
};

class FoldedZXContainer : public ZXContainer<FGSTensor>, public FoldedStackContainer
{
public:
  typedef TensorContainer<FGSTensor> _Ctype;
  FoldedZXContainer(const _Ctype *gss, int ngss, const _Ctype *g, int ng, int ny, int nu)
    : ZXContainer<FGSTensor>(gss, ngss, g, ng, ny, nu)
  {
  }
};

/* This matrix corresponds to
   $$\left[f_{y}\right]+ \left[0
   \left[f_{y^{**}_+}\right]\cdot\left[h^{**}_{y^*}\right] 0\right]$$
   This is very the same as |MatrixA|, the only difference that the
   |MatrixA| is constructed from whole $h_{y^*}$, not only from
   $h^{**}_{y^*}$, hence the new abstraction. */

class MatrixAA : public PLUMatrix
{
public:
  MatrixAA(const FSSparseTensor &f, const IntSequence &ss,
           const TwoDMatrix &gyss, const PartitionY &ypart);
};

/* This class calculates derivatives of $g$ given implicitly by
   $f(h(g^*(y,u,\sigma),\sigma),g(y,u,\sigma),y,u)$, where $h(y,\sigma)$
   is given from outside.

   Structurally, the class is very similar to |KOrder|, but calculations
   are much easier. The two constructors construct an object from sparse
   derivatives of $f$, and derivatives of $h$. The caller must ensure
   that the both derivatives are done at the same point.

   The calculation for order $k$ (including $k=1$) is done by a call
   |performStep(k)|. The derivatives can be retrived by |getFoldDers()|
   or |getUnfoldDers()|. */

class KOrderStoch
{
protected:
  IntSequence nvs;
  PartitionY ypart;
  Journal &journal;
  UGSContainer _ug;
  FGSContainer _fg;
  UGSContainer _ugs;
  FGSContainer _fgs;
  UGSContainer _uG;
  FGSContainer _fG;
  const UGSContainer *_uh;
  const FGSContainer *_fh;
  UnfoldedZXContainer _uZstack;
  FoldedZXContainer _fZstack;
  UnfoldedGXContainer _uGstack;
  FoldedGXContainer _fGstack;
  const TensorContainer<FSSparseTensor> &f;
  MatrixAA matA;
public:
  KOrderStoch(const PartitionY &ypart, int nu, const TensorContainer<FSSparseTensor> &fcont,
              const FGSContainer &hh, Journal &jr);
  KOrderStoch(const PartitionY &ypart, int nu, const TensorContainer<FSSparseTensor> &fcont,
              const UGSContainer &hh, Journal &jr);
  template <int t>
  void performStep(int order);
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
protected:
  template <int t>
  _Ttensor *faaDiBrunoZ(const Symmetry &sym) const;
  template <int t>
  _Ttensor *faaDiBrunoG(const Symmetry &sym) const;

  // convenience access methods
  template<int t>
  _Tg&g();
  template<int t>
  const _Tg&g() const;
  template<int t>
  _Tgs&gs();
  template<int t>
  const _Tgs&gs() const;
  template<int t>
  const _Tgss&h() const;
  template<int t>
  _TG&G();
  template<int t>
  const _TG&G() const;
  template<int t>
  _TZXstack&Zstack();
  template<int t>
  const _TZXstack&Zstack() const;
  template<int t>
  _TGXstack&Gstack();
  template<int t>
  const _TGXstack&Gstack() const;
};

/* This calculates a derivative of $f(G(y,u,\sigma),g(y,u,\sigma),y,u)$
   of a given symmetry. */

template <int t>
_Ttensor *
KOrderStoch::faaDiBrunoZ(const Symmetry &sym) const
{
  JournalRecordPair pa(journal);
  pa << "Faa Di Bruno ZX container for " << sym << endrec;
  _Ttensor *res = new _Ttensor(ypart.ny(), TensorDimens(sym, nvs));
  FaaDiBruno bruno(journal);
  bruno.calculate(Zstack<t>(), f, *res);
  return res;
}

/* This calculates a derivative of
   $G(y,u,\sigma)=h(g^*(y,u,\sigma),\sigma)$ of a given symmetry. */

template <int t>
_Ttensor *
KOrderStoch::faaDiBrunoG(const Symmetry &sym) const
{
  JournalRecordPair pa(journal);
  pa << "Faa Di Bruno GX container for " << sym << endrec;
  TensorDimens tdims(sym, nvs);
  auto *res = new _Ttensor(ypart.nyss(), tdims);
  FaaDiBruno bruno(journal);
  bruno.calculate(Gstack<t>(), h<t>(), *res);
  return res;
}

/* This retrieves all $g$ derivatives of a given dimension from implicit
   $f(h(g^*(y,u,\sigma),\sigma),g(y,u,\sigma),y,u)$. It suppose that all
   derivatives of smaller dimensions have been retrieved.

   So, we go through all symmetries $s$, calculate $G_s$ conditional on
   $g_s=0$, insert the derivative to the $G$ container, then calculate
   $F_s$ conditional on $g_s=0$. This is a righthand side. The left hand
   side is $matA\cdot g_s$. The $g_s$ is retrieved as
   $$g_s=-matA^{-1}\cdot RHS.$$ Finally we have to update $G_s$ by
   calling |Gstack<t>().multAndAdd(1, h<t>(), *G_sym)|. */

template <int t>
void
KOrderStoch::performStep(int order)
{
  int maxd = g<t>().getMaxDim();
  KORD_RAISE_IF(order-1 != maxd && (order != 1 || maxd != -1),
                "Wrong order for KOrderStoch::performStep");
  SymmetrySet ss(order, 4);
  for (symiterator si(ss); !si.isEnd(); ++si)
    {
      if ((*si)[2] == 0)
        {
          JournalRecordPair pa(journal);
          pa << "Recovering symmetry " << *si << endrec;

          _Ttensor *G_sym = faaDiBrunoG<t>(*si);
          G<t>().insert(G_sym);

          _Ttensor *g_sym = faaDiBrunoZ<t>(*si);
          g_sym->mult(-1.0);
          matA.multInv(*g_sym);
          g<t>().insert(g_sym);
          gs<t>().insert(new _Ttensor(ypart.nstat, ypart.nys(), *g_sym));

          Gstack<t>().multAndAdd(1, h<t>(), *G_sym);
        }
    }
}
