/* $Header: /var/lib/cvs/dynare_cpp/sylv/cc/GeneralSylvester.cpp,v 1.1.1.1 2004/06/04 13:00:20 kamenik Exp $ */

/* Tag $Name:  $ */

#include "GeneralSylvester.hh"
#include "SchurDecomp.hh"
#include "SylvException.hh"
#include "TriangularSylvester.hh"
#include "IterativeSylvester.hh"
#include "int_power.hh"

#include <ctime>

GeneralSylvester::GeneralSylvester(int ord, int n, int m, int zero_cols,
                                   const ConstVector &da, const ConstVector &db,
                                   const ConstVector &dc, const ConstVector &dd,
                                   const SylvParams &ps)
  : pars(ps),
    order(ord), a(Vector{da}, n),
    b(Vector{db}, n, n-zero_cols), c(Vector{dc}, m), d(Vector{dd}, n, power(m, order)),
    solved(false)
{
  init();
}

GeneralSylvester::GeneralSylvester(int ord, int n, int m, int zero_cols,
                                   const ConstVector &da, const ConstVector &db,
                                   const ConstVector &dc, Vector &dd,
                                   const SylvParams &ps)
  : pars(ps),
    order(ord), a(Vector{da}, n),
    b(Vector{db}, n, n-zero_cols), c(Vector{dc}, m), d(dd, n, power(m, order)),
    solved(false)
{
  init();
}

GeneralSylvester::GeneralSylvester(int ord, int n, int m, int zero_cols,
                                   const ConstVector &da, const ConstVector &db,
                                   const ConstVector &dc, const ConstVector &dd,
                                   bool alloc_for_check)
  : pars(alloc_for_check),
    order(ord), a(Vector{da}, n),
    b(Vector{db}, n, n-zero_cols), c(Vector{dc}, m), d(Vector{dd}, n, power(m, order)),
    solved(false)
{
  init();
}

GeneralSylvester::GeneralSylvester(int ord, int n, int m, int zero_cols,
                                   const ConstVector &da, const ConstVector &db,
                                   const ConstVector &dc, Vector &dd,
                                   bool alloc_for_check)
  : pars(alloc_for_check),
    order(ord), a(Vector{da}, n),
    b(Vector{db}, n, n-zero_cols), c(Vector{dc}, m), d(dd, n, power(m, order)),
    solved(false)
{
  init();
}

void
GeneralSylvester::init()
{
  GeneralMatrix ainvb(b);
  double rcond1;
  double rcondinf;
  a.multInvLeft2(ainvb, d, rcond1, rcondinf);
  pars.rcondA1 = rcond1;
  pars.rcondAI = rcondinf;
  bdecomp = std::make_unique<SchurDecompZero>(ainvb);
  cdecomp = std::make_unique<SimilarityDecomp>(c.getData(), c.numRows(), *(pars.bs_norm));
  cdecomp->check(pars, c);
  cdecomp->infoToPars(pars);
  if (*(pars.method) == SylvParams::solve_method::recurse)
    sylv = std::make_unique<TriangularSylvester>(*bdecomp, *cdecomp);
  else
    sylv = std::make_unique<IterativeSylvester>(*bdecomp, *cdecomp);
}

void
GeneralSylvester::solve()
{
  if (solved)
    throw SYLV_MES_EXCEPTION("Attempt to run solve() more than once.");

  clock_t start = clock();
  // multiply d
  d.multLeftITrans(bdecomp->getQ());
  d.multRightKron(cdecomp->getQ(), order);
  // convert to KronVector
  KronVector dkron(d.getData(), getM(), getN(), order);
  // solve
  sylv->solve(pars, dkron);
  // multiply d back
  d.multLeftI(bdecomp->getQ());
  d.multRightKron(cdecomp->getInvQ(), order);
  clock_t end = clock();
  pars.cpu_time = ((double) (end-start))/CLOCKS_PER_SEC;

  solved = true;
}

void
GeneralSylvester::check(const ConstVector &ds)
{
  if (!solved)
    throw SYLV_MES_EXCEPTION("Cannot run check on system, which is not solved yet.");

  // calculate xcheck = AX+BXC^i-D
  SylvMatrix dcheck(d.numRows(), d.numCols());
  dcheck.multLeft(b.numRows()-b.numCols(), b, d);
  dcheck.multRightKron(c, order);
  dcheck.multAndAdd(a, d);
  dcheck.getData().add(-1.0, ds);
  // calculate relative norms
  pars.mat_err1 = dcheck.getNorm1()/d.getNorm1();
  pars.mat_errI = dcheck.getNormInf()/d.getNormInf();
  pars.mat_errF = dcheck.getData().getNorm()/d.getData().getNorm();
  pars.vec_err1 = dcheck.getData().getNorm1()/d.getData().getNorm1();
  pars.vec_errI = dcheck.getData().getMax()/d.getData().getMax();
}
