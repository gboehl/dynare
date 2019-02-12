/* $Id: tests.cpp 148 2005-04-19 15:12:26Z kamenik $ */
/* Copyright 2004, Ondra Kamenik */

#include "SylvException.hh"
#include "tl_exception.hh"
#include "gs_tensor.hh"
#include "factory.hh"
#include "monoms.hh"
#include "t_container.hh"
#include "stack_container.hh"
#include "t_polynomial.hh"
#include "rfs_tensor.hh"
#include "ps_tensor.hh"
#include "tl_static.hh"

#include <string>
#include <utility>
#include <memory>
#include <vector>
#include <ctime>
#include <cstdlib>
#include <iostream>
#include <iomanip>

class TestRunnable
{
public:
  const std::string name;
  int dim; // dimension of the solved problem
  int nvar; // number of variable of the solved problem
  TestRunnable(std::string name_arg, int d, int nv)
    : name{std::move(name_arg)}, dim(d), nvar(nv)
  {
  }
  virtual ~TestRunnable() = default;
  bool test() const;
  virtual bool run() const = 0;
protected:
  template<class _Ttype>
  static bool index_forward(const Symmetry &s, const IntSequence &nvs);

  template <class _Ttype>
  static bool index_backward(const Symmetry &s, const IntSequence &nvs);

  template <class _Ttype>
  static bool index_offset(const Symmetry &s, const IntSequence &nvs);

  static bool fold_unfold(const FTensor *folded);
  static bool
  fs_fold_unfold(int r, int nv, int dim)
  {
    Factory f;
    FTensor *folded = f.make<FFSTensor>(r, nv, dim);
    return fold_unfold(folded); // folded deallocated in fold_unfold
  }
  static bool
  r_fold_unfold(int r, int nv, int dim)
  {
    Factory f;
    FTensor *folded = f.make<FRTensor>(r, nv, dim);
    return fold_unfold(folded); // folded deallocated in fold_unfold
  }
  static bool
  gs_fold_unfold(int r, const Symmetry &s, const IntSequence &nvs)
  {
    Factory f;
    FTensor *folded = f.make<FGSTensor>(r, s, nvs);
    return fold_unfold(folded); // folded deallocated in fold_unfold
  }

  static bool dense_prod(const Symmetry &bsym, const IntSequence &bnvs,
                         int hdim, int hnv, int rows);

  static bool folded_monomial(int ng, int nx, int ny, int nu, int dim);

  static bool unfolded_monomial(int ng, int nx, int ny, int nu, int dim);

  static bool fold_zcont(int nf, int ny, int nu, int nup, int nbigg,
                         int ng, int dim);

  static bool unfold_zcont(int nf, int ny, int nu, int nup, int nbigg,
                           int ng, int dim);

  static bool folded_contraction(int r, int nv, int dim);

  static bool unfolded_contraction(int r, int nv, int dim);

  static bool poly_eval(int r, int nv, int maxdim);

};

bool
TestRunnable::test() const
{
  std::cout << "Running test <" << name << ">" << std::endl;
  clock_t start = clock();
  bool passed = run();
  clock_t end = clock();
  std::cout << "CPU time " << static_cast<double>(end-start)/CLOCKS_PER_SEC
            << " (CPU seconds)..................";
  if (passed)
    std::cout << "passed\n\n";
  else
    std::cout << "FAILED\n\n";
  return passed;
}

/****************************************************/
/*     definition of TestRunnable static methods    */
/****************************************************/
template <class _Ttype>
bool
TestRunnable::index_forward(const Symmetry &s, const IntSequence &nvs)
{
  int fails = 0;
  int ndecr = 0;
  int nincr = 0;
  _Ttype dummy(0, TensorDimens(s, nvs));
  typename _Ttype::index run = dummy.end();
  do
    {
      --run;
      ndecr++;
      typename _Ttype::index run2 = dummy.begin();
      for (int i = 0; i < *run; i++)
        {
          ++run2;
          nincr++;
        }
      if (!(run == run2))
        fails++;
    }
  while (run != dummy.begin());

  std::cout << "\tnumber of columns    = " << dummy.ncols() << '\n'
            << "\tnumber of increments = " << nincr << '\n'
            << "\tnumber of decrements = " << ndecr << '\n'
            << "\tnumber of failures   = " << fails << '\n';

  return fails == 0;
}

template <class _Ttype>
bool
TestRunnable::index_backward(const Symmetry &s, const IntSequence &nvs)
{
  int fails = 0;
  int ndecr = 0;
  int nincr = 0;
  _Ttype dummy(0, TensorDimens(s, nvs));
  typename _Ttype::index run = dummy.begin();
  while (run != dummy.end())
    {
      typename _Ttype::index run2 = dummy.end();
      for (int i = 0; i < dummy.ncols() - *run; i++)
        {
          --run2;
          ndecr++;
        }
      if (!(run == run2))
        fails++;
      ++run;
      nincr++;
    }

  std::cout << "\tnumber of columns    = " << dummy.ncols() << '\n'
            << "\tnumber of increments = " << nincr << '\n'
            << "\tnumber of decrements = " << ndecr << '\n'
            << "\tnumber of failures   = " << fails << '\n';

  return fails == 0;
}

template <class _Ttype>
bool
TestRunnable::index_offset(const Symmetry &s, const IntSequence &nvs)
{
  int fails = 0;
  int nincr = 0;
  _Ttype dummy(0, TensorDimens(s, nvs));
  for (typename _Ttype::index run = dummy.begin();
       run != dummy.end(); ++run, nincr++)
    {
      typename _Ttype::index run2(&dummy, run.getCoor());
      if (!(run == run2))
        fails++;
    }

  std::cout << "\tnumber of columns    = " << dummy.ncols() << '\n'
            << "\tnumber of increments = " << nincr << '\n'
            << "\tnumber of failures   = " << fails << '\n';

  return fails == 0;
}

bool
TestRunnable::fold_unfold(const FTensor *folded)
{
  UTensor *unfolded = &(folded->unfold());
  FTensor *folded2 = &(unfolded->fold());
  folded2->add(-1.0, *folded);
  double normInf = folded2->getNormInf();
  double norm1 = folded2->getNorm1();
  std::cout << "\tfolded size:       (" << folded->nrows() << ", " << folded->ncols() << ")\n"
            << "\tunfolded size:     (" << unfolded->nrows() << ", " << unfolded->ncols() << ")\n"
            << "\tdifference normInf: " << normInf << '\n'
            << "\tdifference norm1:   " << norm1 << '\n';

  delete folded;
  delete unfolded;
  delete folded2;

  return normInf < 1.0e-15;
}

bool
TestRunnable::dense_prod(const Symmetry &bsym, const IntSequence &bnvs,
                         int hdim, int hnv, int rows)
{
  Factory f;
  FGSContainer *cont
    = f.makeCont<FGSTensor, FGSContainer>(hnv, bnvs, bsym.dimen()-hdim+1);
  auto *fh
    = f.make<FGSTensor>(rows, Symmetry{hdim}, IntSequence(1, hnv));
  UGSTensor uh(*fh);
  FGSTensor fb(rows, TensorDimens(bsym, bnvs));
  fb.getData().zeros();
  clock_t s1 = clock();
  cont->multAndAdd(uh, fb);
  clock_t s2 = clock();
  UGSContainer ucont(*cont);
  clock_t s3 = clock();
  UGSTensor ub(rows, fb.getDims());
  ub.getData().zeros();
  clock_t s4 = clock();
  ucont.multAndAdd(uh, ub);
  clock_t s5 = clock();

  UGSTensor btmp(fb);
  btmp.add(-1, ub);
  double norm = btmp.getData().getMax();
  double norm1 = btmp.getNorm1();
  double normInf = btmp.getNormInf();

  std::cout << "\ttime for folded product:     " << static_cast<double>(s2-s1)/CLOCKS_PER_SEC << '\n'
            << "\ttime for unfolded product:   " << static_cast<double>(s5-s4)/CLOCKS_PER_SEC << '\n'
            << "\ttime for container convert:  " << static_cast<double>(s3-s2)/CLOCKS_PER_SEC << '\n'
            << "\tunfolded difference normMax: " << norm << '\n'
            << "\tunfolded difference norm1:   " << norm1 << '\n'
            << "\tunfolded difference normInf: " << normInf << '\n';

  delete cont;
  delete fh;

  return norm < 1.e-13;
}

bool
TestRunnable::folded_monomial(int ng, int nx, int ny, int nu, int dim)
{
  clock_t gen_time = clock();
  DenseDerivGenerator gen(ng, nx, ny, nu, 5, 0.3, dim);
  gen_time = clock()-gen_time;
  std::cout << "\ttime for monom generation: "
            << static_cast<double>(gen_time)/CLOCKS_PER_SEC << '\n';
  IntSequence nvs{ny, nu};
  double maxnorm = 0;
  for (int ydim = 0; ydim <= dim; ydim++)
    {
      Symmetry s{ydim, dim-ydim};
      std::cout << "\tSymmetry: ";
      s.print();
      FGSTensor res(ng, TensorDimens(s, nvs));
      res.getData().zeros();
      clock_t stime = clock();
      for (int d = 1; d <= dim; d++)
        gen.xcont->multAndAdd(*(gen.ts[d-1]), res);
      stime = clock() - stime;
      std::cout << "\t\ttime for symmetry: "
                << static_cast<double>(stime)/CLOCKS_PER_SEC << '\n';
      const FGSTensor *mres = gen.rcont->get(s);
      res.add(-1.0, *mres);
      double normtmp = res.getData().getMax();
      std::cout << "\t\terror normMax:     " << normtmp << '\n';
      if (normtmp > maxnorm)
        maxnorm = normtmp;
    }
  return maxnorm < 1.0e-10;
}

bool
TestRunnable::unfolded_monomial(int ng, int nx, int ny, int nu, int dim)
{
  clock_t gen_time = clock();
  DenseDerivGenerator gen(ng, nx, ny, nu, 5, 0.3, dim);
  gen_time = clock()-gen_time;
  std::cout << "\ttime for monom generation: "
            << static_cast<double>(gen_time)/CLOCKS_PER_SEC << '\n';
  clock_t u_time = clock();
  gen.unfold();
  u_time = clock() - u_time;
  std::cout << "\ttime for monom unfolding:  "
            << static_cast<double>(u_time)/CLOCKS_PER_SEC << '\n';
  IntSequence nvs{ny, nu};
  double maxnorm = 0;
  for (int ydim = 0; ydim <= dim; ydim++)
    {
      Symmetry s{ydim, dim-ydim};
      std::cout << "\tSymmetry: ";
      s.print();
      UGSTensor res(ng, TensorDimens(s, nvs));
      res.getData().zeros();
      clock_t stime = clock();
      for (int d = 1; d <= dim; d++)
        gen.uxcont->multAndAdd(*(gen.uts[d-1]), res);
      stime = clock() - stime;
      std::cout << "\t\ttime for symmetry: "
                << static_cast<double>(stime)/CLOCKS_PER_SEC << '\n';
      const FGSTensor *mres = gen.rcont->get(s);
      FGSTensor foldres(res);
      foldres.add(-1.0, *mres);
      double normtmp = foldres.getData().getMax();
      std::cout << "\t\terror normMax:     " << normtmp << '\n';
      if (normtmp > maxnorm)
        maxnorm = normtmp;
    }
  return maxnorm < 1.0e-10;
}

bool
TestRunnable::fold_zcont(int nf, int ny, int nu, int nup, int nbigg,
                         int ng, int dim)
{
  clock_t gen_time = clock();
  SparseDerivGenerator dg(nf, ny, nu, nup, nbigg, ng,
                          5, 0.55, dim);
  gen_time = clock()-gen_time;
  for (int d = 1; d <= dim; d++)
    std::cout << "\tfill of dim=" << d << " tensor:     "
              << std::setprecision(2) << std::fixed << std::setw(6)
              << 100*dg.ts[d-1]->getFillFactor()
              << std::setprecision(6) << std::defaultfloat << " %\n";
  std::cout << "\ttime for monom generation: "
            << static_cast<double>(gen_time)/CLOCKS_PER_SEC << '\n';

  IntSequence nvs{ny, nu, nup, 1};
  double maxnorm = 0.0;

  // form ZContainer
  FoldedZContainer zc(&dg.bigg, nbigg, &dg.g, ng, ny, nu);

  for (int d = 2; d <= dim; d++)
    for (auto &si : SymmetrySet(d, 4))
      {
        std::cout << "\tSymmetry: ";
        si.print();
        FGSTensor res(nf, TensorDimens(si, nvs));
        res.getData().zeros();
        clock_t stime = clock();
        for (int l = 1; l <= si.dimen(); l++)
          zc.multAndAdd(*(dg.ts[l-1]), res);
        stime = clock() - stime;
        std::cout << "\t\ttime for symmetry: "
                  << static_cast<double>(stime)/CLOCKS_PER_SEC << '\n';
        const FGSTensor *mres = dg.rcont.get(si);
        res.add(-1.0, *mres);
        double normtmp = res.getData().getMax();
        std::cout << "\t\terror normMax:     " << normtmp << '\n';
        if (normtmp > maxnorm)
          maxnorm = normtmp;
      }

  return maxnorm < 1.0e-10;
}

bool
TestRunnable::unfold_zcont(int nf, int ny, int nu, int nup, int nbigg,
                           int ng, int dim)
{
  clock_t gen_time = clock();
  SparseDerivGenerator dg(nf, ny, nu, nup, nbigg, ng,
                          5, 0.55, dim);
  gen_time = clock()-gen_time;
  for (int d = 1; d <= dim; d++)
    std::cout << "\tfill of dim=" << d << " tensor:     "
              << std::setprecision(2) << std::fixed << std::setw(6)
              << 100*dg.ts[d-1]->getFillFactor()
              << std::setprecision(6) << std::defaultfloat << " %\n";
  std::cout << "\ttime for monom generation: "
            << static_cast<double>(gen_time)/CLOCKS_PER_SEC << '\n';

  clock_t con_time = clock();
  UGSContainer uG_cont(dg.bigg);
  UGSContainer ug_cont(dg.g);
  con_time = clock()-con_time;
  std::cout << "\ttime for container unfold: "
            << static_cast<double>(con_time)/CLOCKS_PER_SEC << '\n';

  IntSequence nvs{ny, nu, nup, 1};
  double maxnorm = 0.0;

  // form ZContainer
  UnfoldedZContainer zc(&uG_cont, nbigg, &ug_cont, ng, ny, nu);

  for (int d = 2; d <= dim; d++)
    for (auto &si : SymmetrySet(d, 4))
      {
        std::cout << "\tSymmetry: ";
        si.print();
        UGSTensor res(nf, TensorDimens(si, nvs));
        res.getData().zeros();
        clock_t stime = clock();
        for (int l = 1; l <= si.dimen(); l++)
          zc.multAndAdd(*(dg.ts[l-1]), res);
        stime = clock() - stime;
        std::cout << "\t\ttime for symmetry: "
                  << static_cast<double>(stime)/CLOCKS_PER_SEC << '\n';
        FGSTensor fold_res(res);
        const FGSTensor *mres = dg.rcont.get(si);
        fold_res.add(-1.0, *mres);
        double normtmp = fold_res.getData().getMax();
        std::cout << "\t\terror normMax:     " << normtmp << '\n';
        if (normtmp > maxnorm)
          maxnorm = normtmp;
      }

  return maxnorm < 1.0e-10;
}

bool
TestRunnable::folded_contraction(int r, int nv, int dim)
{
  Factory fact;
  Vector *x = fact.makeVector(nv);

  auto *forig = fact.make<FFSTensor>(r, nv, dim);
  auto *f = new FFSTensor(*forig);
  clock_t ctime = clock();
  for (int d = dim-1; d > 0; d--)
    {
      FFSTensor *fnew = new FFSTensor(*f, ConstVector(*x));
      delete f;
      f = fnew;
    }
  ctime = clock() - ctime;
  Vector res(forig->nrows());
  res.zeros();
  f->multaVec(res, *x);

  UFSTensor u(*forig);
  clock_t utime = clock();
  URSingleTensor ux(*x, dim);
  Vector v(u.nrows());
  v.zeros();
  u.multaVec(v, ux.getData());
  utime = clock() - utime;

  v.add(-1.0, res);
  std::cout << "\ttime for folded contraction: "
            << static_cast<double>(ctime)/CLOCKS_PER_SEC << '\n'
            << "\ttime for unfolded power:     "
            << static_cast<double>(utime)/CLOCKS_PER_SEC << '\n'
            << "\terror normMax:     " << v.getMax() << '\n'
            << "\terror norm1:       " << v.getNorm1() << '\n';

  delete f;
  delete x;

  return (v.getMax() < 1.e-10);
}

bool
TestRunnable::unfolded_contraction(int r, int nv, int dim)
{
  Factory fact;
  Vector *x = fact.makeVector(nv);

  auto *forig = fact.make<FFSTensor>(r, nv, dim);
  UFSTensor uorig(*forig);
  delete forig;
  auto *u = new UFSTensor(uorig);
  clock_t ctime = clock();
  for (int d = dim-1; d > 0; d--)
    {
      UFSTensor *unew = new UFSTensor(*u, ConstVector(*x));
      delete u;
      u = unew;
    }
  ctime = clock() - ctime;
  Vector res(uorig.nrows());
  res.zeros();
  u->multaVec(res, *x);

  clock_t utime = clock();
  URSingleTensor ux(*x, dim);
  Vector v(uorig.nrows());
  v.zeros();
  uorig.multaVec(v, ux.getData());
  utime = clock() - utime;

  v.add(-1.0, res);
  std::cout << "\ttime for unfolded contraction: "
            << static_cast<double>(ctime)/CLOCKS_PER_SEC << '\n'
            << "\ttime for unfolded power:       "
            << static_cast<double>(utime)/CLOCKS_PER_SEC << '\n'
            << "\terror normMax:     " << v.getMax() << '\n'
            << "\terror norm1:       " << v.getNorm1() << '\n';

  delete u;
  delete x;

  return (v.getMax() < 1.e-10);
}

bool
TestRunnable::poly_eval(int r, int nv, int maxdim)
{
  Factory fact;
  Vector *x = fact.makeVector(nv);

  Vector out_ft(r);
  out_ft.zeros();
  Vector out_fh(r);
  out_fh.zeros();
  Vector out_ut(r);
  out_ut.zeros();
  Vector out_uh(r);
  out_uh.zeros();

  UTensorPolynomial *up;
  {
    FTensorPolynomial *fp = fact.makePoly<FFSTensor, FTensorPolynomial>(r, nv, maxdim);

    clock_t ft_cl = clock();
    fp->evalTrad(out_ft, *x);
    ft_cl = clock() - ft_cl;
    std::cout << "\ttime for folded power eval:    "
              << static_cast<double>(ft_cl)/CLOCKS_PER_SEC << '\n';

    clock_t fh_cl = clock();
    fp->evalHorner(out_fh, *x);
    fh_cl = clock() - fh_cl;
    std::cout << "\ttime for folded horner eval:   "
              << static_cast<double>(fh_cl)/CLOCKS_PER_SEC << '\n';

    up = new UTensorPolynomial(*fp);
    delete fp;
  }

  clock_t ut_cl = clock();
  up->evalTrad(out_ut, *x);
  ut_cl = clock() - ut_cl;
  std::cout << "\ttime for unfolded power eval:  "
            << static_cast<double>(ut_cl)/CLOCKS_PER_SEC << '\n';

  clock_t uh_cl = clock();
  up->evalHorner(out_uh, *x);
  uh_cl = clock() - uh_cl;
  std::cout << "\ttime for unfolded horner eval: "
            << static_cast<double>(uh_cl)/CLOCKS_PER_SEC << '\n';

  out_ft.add(-1.0, out_ut);
  double max_ft = out_ft.getMax();
  out_fh.add(-1.0, out_ut);
  double max_fh = out_fh.getMax();
  out_uh.add(-1.0, out_ut);
  double max_uh = out_uh.getMax();

  std::cout << "\tfolded power error norm max:     " << max_ft << '\n'
            << "\tfolded horner error norm max:    " << max_fh << '\n'
            << "\tunfolded horner error norm max:  " << max_uh << '\n';

  delete up;
  delete x;
  return (max_ft+max_fh+max_uh < 1.0e-10);
}

/****************************************************/
/*     definition of TestRunnable subclasses        */
/****************************************************/
class SmallIndexForwardFold : public TestRunnable
{
public:
  SmallIndexForwardFold()
    : TestRunnable("small index forward for fold (44)(222)", 5, 4)
  {
  }
  bool
  run() const override
  {
    Symmetry s{2, 3};
    IntSequence nvs{4, 2};
    return index_forward<FGSTensor>(s, nvs);
  }
};

class SmallIndexForwardUnfold : public TestRunnable
{
public:
  SmallIndexForwardUnfold()
    : TestRunnable("small index forward for unfold (44)(222)", 5, 4)
  {
  }
  bool
  run() const override
  {
    Symmetry s{2, 3};
    IntSequence nvs{4, 2};
    return index_forward<UGSTensor>(s, nvs);
  }
};

class IndexForwardFold : public TestRunnable
{
public:
  IndexForwardFold()
    : TestRunnable("index forward for fold (55)(222)(22)", 7, 5)
  {
  }
  bool
  run() const override
  {
    Symmetry s{2, 3, 2};
    IntSequence nvs{5, 2, 2};
    return index_forward<FGSTensor>(s, nvs);
  }
};

class IndexForwardUnfold : public TestRunnable
{
public:
  IndexForwardUnfold()
    : TestRunnable("index forward for unfold (55)(222)(22)", 7, 5)
  {
  }
  bool
  run() const override
  {
    Symmetry s{2, 3, 2};
    IntSequence nvs{5, 2, 2};
    return index_forward<UGSTensor>(s, nvs);
  }
};

class SmallIndexBackwardFold : public TestRunnable
{
public:
  SmallIndexBackwardFold()
    : TestRunnable("small index backward for fold (3)(3)(222)", 5, 3)
  {
  }
  bool
  run() const override
  {
    Symmetry s{1, 1, 3};
    IntSequence nvs{3, 3, 2};
    return index_backward<FGSTensor>(s, nvs);
  }
};

class IndexBackwardFold : public TestRunnable
{
public:
  IndexBackwardFold()
    : TestRunnable("index backward for fold (44)(222)(44)", 7, 4)
  {
  }
  bool
  run() const override
  {
    Symmetry s{2, 3, 2};
    IntSequence nvs{4, 2, 4};
    return index_backward<FGSTensor>(s, nvs);
  }
};

class SmallIndexBackwardUnfold : public TestRunnable
{
public:
  SmallIndexBackwardUnfold()
    : TestRunnable("small index backward for unfold (3)(3)(222)", 5, 3)
  {
  }
  bool
  run() const override
  {
    Symmetry s{1, 1, 3};
    IntSequence nvs{3, 3, 2};
    return index_backward<UGSTensor>(s, nvs);
  }
};

class IndexBackwardUnfold : public TestRunnable
{
public:
  IndexBackwardUnfold()
    : TestRunnable("index backward for unfold (44)(222)(44)", 7, 4)
  {
  }
  bool
  run() const override
  {
    Symmetry s{2, 3, 2};
    IntSequence nvs{4, 2, 4};
    return index_backward<UGSTensor>(s, nvs);
  }
};

class SmallIndexOffsetFold : public TestRunnable
{
public:
  SmallIndexOffsetFold()
    : TestRunnable("small index offset for fold (44)(222)", 5, 4)
  {
  }
  bool
  run() const override
  {
    Symmetry s{2, 3};
    IntSequence nvs{4, 2};
    return index_offset<FGSTensor>(s, nvs);
  }
};

class SmallIndexOffsetUnfold : public TestRunnable
{
public:
  SmallIndexOffsetUnfold()
    : TestRunnable("small index offset for unfold (44)(222)", 5, 4)
  {
  }
  bool
  run() const override
  {
    Symmetry s{2, 3};
    IntSequence nvs{4, 2};
    return index_offset<UGSTensor>(s, nvs);
  }
};

class IndexOffsetFold : public TestRunnable
{
public:
  IndexOffsetFold()
    : TestRunnable("index offset for fold (55)(222)(22)", 5, 5)
  {
  }
  bool
  run() const override
  {
    Symmetry s{2, 3, 2};
    IntSequence nvs{5, 2, 2};
    return index_offset<FGSTensor>(s, nvs);
  }
};

class IndexOffsetUnfold : public TestRunnable
{
public:
  IndexOffsetUnfold()
    : TestRunnable("index offset for unfold (55)(222)(22)", 7, 5)
  {
  }
  bool
  run() const override
  {
    Symmetry s{2, 3, 2};
    IntSequence nvs{5, 2, 2};
    return index_offset<UGSTensor>(s, nvs);
  }
};

class SmallFoldUnfoldFS : public TestRunnable
{
public:
  SmallFoldUnfoldFS()
    : TestRunnable("small fold-unfold for full symmetry (444)", 3, 4)
  {
  }
  bool
  run() const override
  {
    return fs_fold_unfold(5, 4, 3);
  }
};

class SmallFoldUnfoldGS : public TestRunnable
{
public:
  SmallFoldUnfoldGS()
    : TestRunnable("small fold-unfold for gen symmetry (3)(33)(22)", 5, 3)
  {
  }
  bool
  run() const override
  {
    Symmetry s{1, 2, 2};
    IntSequence nvs{3, 3, 2};
    return gs_fold_unfold(5, s, nvs);
  }
};

class FoldUnfoldFS : public TestRunnable
{
public:
  FoldUnfoldFS()
    : TestRunnable("fold-unfold for full symmetry (9999)", 4, 9)
  {
  }
  bool
  run() const override
  {
    return fs_fold_unfold(5, 9, 4);
  }
};

class FoldUnfoldGS : public TestRunnable
{
public:
  FoldUnfoldGS()
    : TestRunnable("fold-unfold for gen symmetry (66)(2)(66)", 5, 6)
  {
  }
  bool
  run() const override
  {
    Symmetry s{2, 1, 2};
    IntSequence nvs{6, 2, 6};
    return gs_fold_unfold(5, s, nvs);
  }
};

class SmallFoldUnfoldR : public TestRunnable
{
public:
  SmallFoldUnfoldR()
    : TestRunnable("small fold-unfold for row full symmetry (333)", 3, 3)
  {
  }
  bool
  run() const override
  {
    return r_fold_unfold(5, 3, 3);
  }
};

class FoldUnfoldR : public TestRunnable
{
public:
  FoldUnfoldR()
    : TestRunnable("fold-unfold for row full symmetry (66666)", 5, 6)
  {
  }
  bool
  run() const override
  {
    return r_fold_unfold(5, 6, 5);
  }
};

class SmallDenseProd : public TestRunnable
{
public:
  SmallDenseProd()
    : TestRunnable("small dense prod bsym=1-2,nvs=3-2,h=2-3,r=2", 3, 3)
  {
  }
  bool
  run() const override
  {
    IntSequence bnvs{3, 2};
    return dense_prod(Symmetry{1, 2}, bnvs, 2, 3, 2);
  }
};

class DenseProd : public TestRunnable
{
public:
  DenseProd()
    : TestRunnable("dense prod bsym=2-3,nvs=10-7,h=3-15,r=10", 5, 15)
  {
  }
  bool
  run() const override
  {
    IntSequence bnvs{10, 7};
    return dense_prod(Symmetry{2, 3}, bnvs, 3, 15, 10);
  }
};

class BigDenseProd : public TestRunnable
{
public:
  BigDenseProd()
    : TestRunnable("dense prod bsym=3-2,nvs=13-11,h=3-20,r=20", 6, 20)
  {
  }
  bool
  run() const override
  {
    IntSequence bnvs{13, 11};
    return dense_prod(Symmetry{3, 2}, bnvs, 3, 20, 20);
  }
};

class SmallFoldedMonomial : public TestRunnable
{
public:
  SmallFoldedMonomial()
    : TestRunnable("folded vrs. monoms (g,x,y,u)=(10,4,5,3), dim=4", 4, 8)
  {
  }
  bool
  run() const override
  {
    return folded_monomial(10, 4, 5, 3, 4);
  }
};

class FoldedMonomial : public TestRunnable
{
public:
  FoldedMonomial()
    : TestRunnable("folded vrs. monoms (g,x,y,u)=(20,12,10,5), dim=4", 4, 15)
  {
  }
  bool
  run() const override
  {
    return folded_monomial(20, 12, 10, 5, 4);
  }
};

class SmallUnfoldedMonomial : public TestRunnable
{
public:
  SmallUnfoldedMonomial()
    : TestRunnable("unfolded vrs. monoms (g,x,y,u)=(10,4,5,3), dim=4", 4, 8)
  {
  }
  bool
  run() const override
  {
    return unfolded_monomial(10, 4, 5, 3, 4);
  }
};

class UnfoldedMonomial : public TestRunnable
{
public:
  UnfoldedMonomial()
    : TestRunnable("unfolded vrs. monoms (g,x,y,u)=(20,12,10,5), dim=4", 4, 15)
  {
  }
  bool
  run() const override
  {
    return unfolded_monomial(20, 12, 10, 5, 4);
  }
};

class FoldedContractionSmall : public TestRunnable
{
public:
  FoldedContractionSmall()
    : TestRunnable("folded contraction small (r=5, nv=4, dim=3)", 3, 4)
  {
  }
  bool
  run() const override
  {
    return folded_contraction(5, 4, 3);
  }
};

class FoldedContractionBig : public TestRunnable
{
public:
  FoldedContractionBig()
    : TestRunnable("folded contraction big (r=20, nv=12, dim=5)", 5, 12)
  {
  }
  bool
  run() const override
  {
    return folded_contraction(20, 12, 5);
  }
};

class UnfoldedContractionSmall : public TestRunnable
{
public:
  UnfoldedContractionSmall()
    : TestRunnable("unfolded contraction small (r=5, nv=4, dim=3)", 3, 4)
  {
  }
  bool
  run() const override
  {
    return unfolded_contraction(5, 4, 3);
  }
};

class UnfoldedContractionBig : public TestRunnable
{
public:
  UnfoldedContractionBig()
    : TestRunnable("unfolded contraction big (r=20, nv=12, dim=5)", 5, 12)
  {
  }
  bool
  run() const override
  {
    return unfolded_contraction(20, 12, 5);
  }
};

class PolyEvalSmall : public TestRunnable
{
public:
  PolyEvalSmall()
    : TestRunnable("polynomial evaluation small (r=4, nv=5, maxdim=4)", 4, 5)
  {
  }
  bool
  run() const override
  {
    return poly_eval(4, 5, 4);
  }
};

class PolyEvalBig : public TestRunnable
{
public:
  PolyEvalBig()
    : TestRunnable("polynomial evaluation big (r=244, nv=97, maxdim=2)", 2, 97)
  {
  }
  bool
  run() const override
  {
    return poly_eval(244, 97, 2);
  }
};

class FoldZContSmall : public TestRunnable
{
public:
  FoldZContSmall()
    : TestRunnable("folded Z container (r=3,ny=2,nu=2,nup=1,G=2,g=2,dim=3)",
                   3, 8)
  {
  }
  bool
  run() const override
  {
    return fold_zcont(3, 2, 2, 1, 2, 2, 3);
  }
};

class FoldZCont : public TestRunnable
{
public:
  FoldZCont()
    : TestRunnable("folded Z container (r=13,ny=5,nu=7,nup=4,G=6,g=7,dim=4)",
                   4, 25)
  {
  }
  bool
  run() const override
  {
    return fold_zcont(13, 5, 7, 4, 6, 7, 4);
  }
};

class UnfoldZContSmall : public TestRunnable
{
public:
  UnfoldZContSmall()
    : TestRunnable("unfolded Z container (r=3,ny=2,nu=2,nup=1,G=2,g=2,dim=3)",
                   3, 8)
  {
  }
  bool
  run() const override
  {
    return unfold_zcont(3, 2, 2, 1, 2, 2, 3);
  }
};

class UnfoldZCont : public TestRunnable
{
public:
  UnfoldZCont()
    : TestRunnable("unfolded Z container (r=13,ny=5,nu=7,nup=4,G=6,g=7,dim=4",
                   4, 25)
  {
  }
  bool
  run() const override
  {
    return unfold_zcont(13, 5, 7, 4, 6, 7, 4);
  }
};

int
main()
{
  std::vector<std::unique_ptr<TestRunnable>> all_tests;
  // fill in vector of all tests
  all_tests.push_back(std::make_unique<SmallIndexForwardFold>());
  all_tests.push_back(std::make_unique<SmallIndexForwardUnfold>());
  all_tests.push_back(std::make_unique<IndexForwardFold>());
  all_tests.push_back(std::make_unique<IndexForwardUnfold>());
  all_tests.push_back(std::make_unique<SmallIndexBackwardFold>());
  all_tests.push_back(std::make_unique<IndexBackwardFold>());
  all_tests.push_back(std::make_unique<SmallIndexBackwardUnfold>());
  all_tests.push_back(std::make_unique<IndexBackwardUnfold>());
  all_tests.push_back(std::make_unique<SmallIndexOffsetFold>());
  all_tests.push_back(std::make_unique<SmallIndexOffsetUnfold>());
  all_tests.push_back(std::make_unique<IndexOffsetFold>());
  all_tests.push_back(std::make_unique<IndexOffsetUnfold>());
  all_tests.push_back(std::make_unique<SmallFoldUnfoldFS>());
  all_tests.push_back(std::make_unique<SmallFoldUnfoldGS>());
  all_tests.push_back(std::make_unique<FoldUnfoldFS>());
  all_tests.push_back(std::make_unique<FoldUnfoldGS>());
  all_tests.push_back(std::make_unique<SmallFoldUnfoldR>());
  all_tests.push_back(std::make_unique<FoldUnfoldR>());
  all_tests.push_back(std::make_unique<SmallDenseProd>());
  all_tests.push_back(std::make_unique<DenseProd>());
  all_tests.push_back(std::make_unique<BigDenseProd>());
  all_tests.push_back(std::make_unique<SmallFoldedMonomial>());
  all_tests.push_back(std::make_unique<FoldedMonomial>());
  all_tests.push_back(std::make_unique<SmallUnfoldedMonomial>());
  all_tests.push_back(std::make_unique<UnfoldedMonomial>());
  all_tests.push_back(std::make_unique<FoldedContractionSmall>());
  all_tests.push_back(std::make_unique<FoldedContractionBig>());
  all_tests.push_back(std::make_unique<UnfoldedContractionSmall>());
  all_tests.push_back(std::make_unique<UnfoldedContractionBig>());
  all_tests.push_back(std::make_unique<PolyEvalSmall>());
  all_tests.push_back(std::make_unique<PolyEvalBig>());
  all_tests.push_back(std::make_unique<FoldZContSmall>());
  all_tests.push_back(std::make_unique<FoldZCont>());
  all_tests.push_back(std::make_unique<UnfoldZContSmall>());
  all_tests.push_back(std::make_unique<UnfoldZCont>());

  // find maximum dimension and maximum nvar
  int dmax = 0;
  int nvmax = 0;
  for (const auto &test : all_tests)
    {
      if (dmax < test->dim)
        dmax = test->dim;
      if (nvmax < test->nvar)
        nvmax = test->nvar;
    }
  tls.init(dmax, nvmax); // initialize library

  // launch the tests
  int success = 0;
  for (const auto &test : all_tests)
    {
      try
        {
          if (test->test())
            success++;
        }
      catch (const TLException &e)
        {
          std::cout << "Caught TL exception in <" << test->name << ">:\n";
          e.print();
        }
      catch (SylvException &e)
        {
          std::cout << "Caught Sylv exception in <" << test->name << ">:\n";
          e.printMessage();
        }
    }

  int nfailed = all_tests.size() - success;
  std::cout << "There were " << nfailed << " tests that failed out of "
            << all_tests.size() << " tests run." << std::endl;

  if (nfailed)
    return EXIT_FAILURE;
  else
    return EXIT_SUCCESS;
}
