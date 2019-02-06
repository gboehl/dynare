/* $Id: tests.cpp 431 2005-08-16 15:41:01Z kamenik $ */
/* Copyright 2005, Ondra Kamenik */

#include "GeneralMatrix.hh"
#include <dynlapack.h>
#include "SylvException.hh"

#include "rfs_tensor.hh"
#include "normal_moments.hh"

#include "vector_function.hh"
#include "quadrature.hh"
#include "smolyak.hh"
#include "product.hh"
#include "quasi_mcarlo.hh"

#include <iomanip>
#include <chrono>
#include <cmath>
#include <iostream>
#include <utility>
#include <array>
#include <memory>
#include <cstdlib>

// evaluates unfolded (Dx)^k power, where x is a vector, D is a
// Cholesky factor (lower triangular)
class MomentFunction : public VectorFunction
{
  GeneralMatrix D;
  int k;
public:
  MomentFunction(const GeneralMatrix &inD, int kk)
    : VectorFunction(inD.numRows(), UFSTensor::calcMaxOffset(inD.numRows(), kk)),
      D(inD), k(kk)
  {
  }
  MomentFunction(const MomentFunction &func) = default;

  std::unique_ptr<VectorFunction>
  clone() const override
  {
    return std::make_unique<MomentFunction>(*this);
  }
  void eval(const Vector &point, const ParameterSignal &sig, Vector &out) override;
};

void
MomentFunction::eval(const Vector &point, const ParameterSignal &sig, Vector &out)
{
  if (point.length() != indim() || out.length() != outdim())
    {
      std::cerr << "Wrong length of vectors in MomentFunction::eval" << std::endl;
      std::exit(EXIT_FAILURE);
    }
  Vector y(point);
  y.zeros();
  D.multaVec(y, point);
  URSingleTensor ypow(y, k);
  out.zeros();
  out.add(1.0, ypow.getData());
}

class TensorPower : public VectorFunction
{
  int k;
public:
  TensorPower(int nvar, int kk)
    : VectorFunction(nvar, UFSTensor::calcMaxOffset(nvar, kk)), k(kk)
  {
  }
  TensorPower(const TensorPower &func) = default;

  std::unique_ptr<VectorFunction>
  clone() const override
  {
    return std::make_unique<TensorPower>(*this);
  }
  void eval(const Vector &point, const ParameterSignal &sig, Vector &out) override;
};

void
TensorPower::eval(const Vector &point, const ParameterSignal &sig, Vector &out)
{
  if (point.length() != indim() || out.length() != outdim())
    {
      std::cerr << "Wrong length of vectors in TensorPower::eval" << std::endl;
      std::exit(EXIT_FAILURE);
    }
  URSingleTensor ypow(point, k);
  out.zeros();
  out.add(1.0, ypow.getData());
}

// evaluates (1+1/d)^d*(x_1*...*x_d)^(1/d), its integral over <0,1>^d
// is 1.0, and its variation grows exponetially
// d = dim
class Function1 : public VectorFunction
{
  int dim;
public:
  Function1(int d)
    : VectorFunction(d, 1), dim(d)
  {
  }
  Function1(const Function1 &f)
    : VectorFunction(f.indim(), f.outdim()), dim(f.dim)
  {
  }
  std::unique_ptr<VectorFunction>
  clone() const override
  {
    return std::make_unique<Function1>(*this);
  }
  void eval(const Vector &point, const ParameterSignal &sig, Vector &out) override;
};

void
Function1::eval(const Vector &point, const ParameterSignal &sig, Vector &out)
{
  if (point.length() != dim || out.length() != 1)
    {
      std::cerr << "Wrong length of vectors in Function1::eval" << std::endl;
      std::exit(EXIT_FAILURE);
    }
  double r = 1;
  for (int i = 0; i < dim; i++)
    r *= point[i];
  r = pow(r, 1.0/dim);
  r *= pow(1.0 + 1.0/dim, (double) dim);
  out[0] = r;
}

// evaluates Function1 but with transformation x_i=0.5(y_i+1)
// this makes the new function integrate over <-1,1>^d to 1.0
class Function1Trans : public Function1
{
public:
  Function1Trans(int d)
    : Function1(d)
  {
  }
  Function1Trans(const Function1Trans &func) = default;

  std::unique_ptr<VectorFunction>
  clone() const override
  {
    return std::make_unique<Function1Trans>(*this);
  }
  void eval(const Vector &point, const ParameterSignal &sig, Vector &out) override;
};

void
Function1Trans::eval(const Vector &point, const ParameterSignal &sig, Vector &out)
{
  Vector p(point.length());
  for (int i = 0; i < p.length(); i++)
    p[i] = 0.5*(point[i]+1);
  Function1::eval(p, sig, out);
  out.mult(pow(0.5, indim()));
}

// WallTimer class. Constructor saves the wall time, destructor
// cancels the current time from the saved, and prints the message
// with time information
class WallTimer
{
  std::string mes;
  std::chrono::time_point<std::chrono::high_resolution_clock> start;
  bool new_line;
public:
  WallTimer(std::string m, bool nl = true)
    : mes{m}, start{std::chrono::high_resolution_clock::now()}, new_line{nl}
  {
  }
  ~WallTimer()
  {
    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> duration = end - start;
    std::cout << mes << std::setw(8) << std::setprecision(4) << duration.count();
    if (new_line)
      std::cout << std::endl;
  }
};

/****************************************************/
/*     declaration of TestRunnable class            */
/****************************************************/
class TestRunnable
{
public:
  const std::string name;
  int dim; // dimension of the solved problem
  int nvar; // number of variable of the solved problem
  TestRunnable(std::string name_arg, int d, int nv)
    : name{move(name_arg)}, dim(d), nvar(nv)
  {
  }
  virtual ~TestRunnable() = default;
  bool test() const;
  virtual bool run() const = 0;
protected:
  static bool smolyak_normal_moments(const GeneralMatrix &m, int imom, int level);
  static bool product_normal_moments(const GeneralMatrix &m, int imom, int level);
  static bool qmc_normal_moments(const GeneralMatrix &m, int imom, int level);
  static bool smolyak_product_cube(const VectorFunction &func, const Vector &res,
                                   double tol, int level);
  static bool qmc_cube(const VectorFunction &func, double res, double tol, int level);
};

bool
TestRunnable::test() const
{
  std::cout << "Running test <" << name << ">" << std::endl;
  bool passed;
  {
    WallTimer tim("Wall clock time ", false);
    passed = run();
  }
  if (passed)
    {
      std::cout << "............................ passed" << std::endl << std::endl;
      return passed;
    }
  else
    {
      std::cout << "............................ FAILED" << std::endl << std::endl;
      return passed;
    }
}

/****************************************************/
/*     definition of TestRunnable static methods    */
/****************************************************/
bool
TestRunnable::smolyak_normal_moments(const GeneralMatrix &m, int imom, int level)
{
  // first make m*m' and then Cholesky factor
  GeneralMatrix mtr(m, "transpose");
  GeneralMatrix msq(m, mtr);

  // make vector function
  int dim = m.numRows();
  TensorPower tp(dim, imom);
  GaussConverterFunction func(tp, msq);

  // smolyak quadrature
  Vector smol_out(UFSTensor::calcMaxOffset(dim, imom));
  {
    WallTimer tim("\tSmolyak quadrature time:         ");
    GaussHermite gs;
    SmolyakQuadrature quad(dim, level, gs);
    quad.integrate(func, level, sthread::detach_thread_group::max_parallel_threads, smol_out);
    std::cout << "\tNumber of Smolyak evaluations:    " << quad.numEvals(level) << std::endl;
  }

  // check against theoretical moments
  UNormalMoments moments(imom, msq);
  smol_out.add(-1.0, (moments.get(Symmetry(imom)))->getData());
  std::cout << "\tError:                         " << std::setw(16) << std::setprecision(12) << smol_out.getMax() << std::endl;
  return smol_out.getMax() < 1.e-7;
}

bool
TestRunnable::product_normal_moments(const GeneralMatrix &m, int imom, int level)
{
  // first make m*m' and then Cholesky factor
  GeneralMatrix mtr(m, "transpose");
  GeneralMatrix msq(m, mtr);

  // make vector function
  int dim = m.numRows();
  TensorPower tp(dim, imom);
  GaussConverterFunction func(tp, msq);

  // product quadrature
  Vector prod_out(UFSTensor::calcMaxOffset(dim, imom));
  {
    WallTimer tim("\tProduct quadrature time:         ");
    GaussHermite gs;
    ProductQuadrature quad(dim, gs);
    quad.integrate(func, level, sthread::detach_thread_group::max_parallel_threads, prod_out);
    std::cout << "\tNumber of product evaluations:    " << quad.numEvals(level) << std::endl;
  }

  // check against theoretical moments
  UNormalMoments moments(imom, msq);
  prod_out.add(-1.0, (moments.get(Symmetry(imom)))->getData());
  std::cout << "\tError:                         " << std::setw(16) << std::setprecision(12) << prod_out.getMax() << std::endl;
  return prod_out.getMax() < 1.e-7;
}

bool
TestRunnable::smolyak_product_cube(const VectorFunction &func, const Vector &res,
                                   double tol, int level)
{
  if (res.length() != func.outdim())
    {
      std::cerr << "Incompatible dimensions of check value and function." << std::endl;
      std::exit(EXIT_FAILURE);
    }

  GaussLegendre glq;
  Vector out(func.outdim());
  double smol_error;
  double prod_error;
  {
    WallTimer tim("\tSmolyak quadrature time:         ");
    SmolyakQuadrature quad(func.indim(), level, glq);
    quad.integrate(func, level, sthread::detach_thread_group::max_parallel_threads, out);
    out.add(-1.0, res);
    smol_error = out.getMax();
    std::cout << "\tNumber of Smolyak evaluations:    " << quad.numEvals(level) << std::endl;
    std::cout << "\tError:                            " << std::setw(16) << std::setprecision(12) << smol_error << std::endl;
  }
  {
    WallTimer tim("\tProduct quadrature time:         ");
    ProductQuadrature quad(func.indim(), glq);
    quad.integrate(func, level, sthread::detach_thread_group::max_parallel_threads, out);
    out.add(-1.0, res);
    prod_error = out.getMax();
    std::cout << "\tNumber of product evaluations:    " << quad.numEvals(level) << std::endl;
    std::cout << "\tError:                            " << std::setw(16) << std::setprecision(12) << prod_error << std::endl;
  }

  return smol_error < tol && prod_error < tol;
}

bool
TestRunnable::qmc_cube(const VectorFunction &func, double res, double tol, int level)
{
  Vector r(1);
  double error1;
  {
    WallTimer tim("\tQuasi-Monte Carlo (Warnock scrambling) time:  ");
    WarnockPerScheme wps;
    QMCarloCubeQuadrature qmc(func.indim(), level, wps);
    //		qmc.savePoints("warnock.txt", level);
    qmc.integrate(func, level, sthread::detach_thread_group::max_parallel_threads, r);
    error1 = std::max(res - r[0], r[0] - res);
    std::cout << "\tQuasi-Monte Carlo (Warnock scrambling) error: " << std::setw(16) << std::setprecision(12) << error1 << std::endl;
  }
  double error2;
  {
    WallTimer tim("\tQuasi-Monte Carlo (reverse scrambling) time:  ");
    ReversePerScheme rps;
    QMCarloCubeQuadrature qmc(func.indim(), level, rps);
    //		qmc.savePoints("reverse.txt", level);
    qmc.integrate(func, level, sthread::detach_thread_group::max_parallel_threads, r);
    error2 = std::max(res - r[0], r[0] - res);
    std::cout << "\tQuasi-Monte Carlo (reverse scrambling) error: " << std::setw(16) << std::setprecision(12) << error2 << std::endl;
  }
  double error3;
  {
    WallTimer tim("\tQuasi-Monte Carlo (no scrambling) time:       ");
    IdentityPerScheme ips;
    QMCarloCubeQuadrature qmc(func.indim(), level, ips);
    //		qmc.savePoints("identity.txt", level);
    qmc.integrate(func, level, sthread::detach_thread_group::max_parallel_threads, r);
    error3 = std::max(res - r[0], r[0] - res);
    std::cout << "\tQuasi-Monte Carlo (no scrambling) error:      " << std::setw(16) << std::setprecision(12) << error3 << std::endl;
  }

  return error1 < tol && error2 < tol && error3 < tol;
}

/****************************************************/
/*     definition of TestRunnable subclasses        */
/****************************************************/
class SmolyakNormalMom1 : public TestRunnable
{
public:
  SmolyakNormalMom1()
    : TestRunnable("Smolyak normal moments (dim=2, level=4, order=4)", 4, 2)
  {
  }

  bool
  run() const override
  {
    GeneralMatrix m(2, 2);
    m.zeros(); m.get(0, 0) = 1; m.get(1, 1) = 1;
    return smolyak_normal_moments(m, 4, 4);
  }
};

class SmolyakNormalMom2 : public TestRunnable
{
public:
  SmolyakNormalMom2()
    : TestRunnable("Smolyak normal moments (dim=3, level=8, order=8)", 8, 3)
  {
  }

  bool
  run() const override
  {
    GeneralMatrix m(3, 3);
    m.zeros();
    m.get(0, 0) = 1; m.get(0, 2) = 0.5; m.get(1, 1) = 1;
    m.get(1, 0) = 0.5; m.get(2, 2) = 2; m.get(2, 1) = 4;
    return smolyak_normal_moments(m, 8, 8);
  }
};

class ProductNormalMom1 : public TestRunnable
{
public:
  ProductNormalMom1()
    : TestRunnable("Product normal moments (dim=2, level=4, order=4)", 4, 2)
  {
  }

  bool
  run() const override
  {
    GeneralMatrix m(2, 2);
    m.zeros(); m.get(0, 0) = 1; m.get(1, 1) = 1;
    return product_normal_moments(m, 4, 4);
  }
};

class ProductNormalMom2 : public TestRunnable
{
public:
  ProductNormalMom2()
    : TestRunnable("Product normal moments (dim=3, level=8, order=8)", 8, 3)
  {
  }

  bool
  run() const override
  {
    GeneralMatrix m(3, 3);
    m.zeros();
    m.get(0, 0) = 1; m.get(0, 2) = 0.5; m.get(1, 1) = 1;
    m.get(1, 0) = 0.5; m.get(2, 2) = 2; m.get(2, 1) = 4;
    return product_normal_moments(m, 8, 8);
  }
};

// note that here we pass 1,1 to tls since smolyak has its own PascalTriangle
class F1GaussLegendre : public TestRunnable
{
public:
  F1GaussLegendre()
    : TestRunnable("Function1 Gauss-Legendre (dim=6, level=13", 1, 1)
  {
  }

  bool
  run() const override
  {
    Function1Trans f1(6);
    Vector res(1); res[0] = 1.0;
    return smolyak_product_cube(f1, res, 1e-2, 13);
  }
};

class F1QuasiMCarlo : public TestRunnable
{
public:
  F1QuasiMCarlo()
    : TestRunnable("Function1 Quasi-Monte Carlo (dim=6, level=1000000)", 1, 1)
  {
  }

  bool
  run() const override
  {
    Function1 f1(6);
    return qmc_cube(f1, 1.0, 1.e-4, 1000000);
  }
};

int
main()
{
  std::vector<std::unique_ptr<TestRunnable>> all_tests;
  // fill in vector of all tests
  all_tests.push_back(std::make_unique<SmolyakNormalMom1>());
  all_tests.push_back(std::make_unique<SmolyakNormalMom2>());
  all_tests.push_back(std::make_unique<ProductNormalMom1>());
  all_tests.push_back(std::make_unique<ProductNormalMom2>());
  all_tests.push_back(std::make_unique<F1GaussLegendre>());
  all_tests.push_back(std::make_unique<F1QuasiMCarlo>());

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
          std::cout << "Caught TL exception in <" << test->name << ">:" << std::endl;
          e.print();
        }
      catch (SylvException &e)
        {
          std::cout << "Caught Sylv exception in <" << test->name << ">:" << std::endl;
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
