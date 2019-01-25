/* $Header: /var/lib/cvs/dynare_cpp/sylv/cc/GeneralSylvester.h,v 1.1.1.1 2004/06/04 13:00:20 kamenik Exp $ */

/* Tag $Name:  $ */

#ifndef GENERAL_SYLVESTER_H
#define GENERAL_SYLVESTER_H

#include "SylvMatrix.hh"
#include "SimilarityDecomp.hh"
#include "SylvesterSolver.hh"

#include <memory>

class GeneralSylvester
{
  SylvParams pars;
  int order;
  const SqSylvMatrix a;
  const SylvMatrix b;
  const SqSylvMatrix c;
  SylvMatrix d;
  bool solved;
  std::unique_ptr<SchurDecompZero> bdecomp;
  std::unique_ptr<SimilarityDecomp> cdecomp;
  std::unique_ptr<SylvesterSolver> sylv;
public:
  /* construct with my copy of d*/
  GeneralSylvester(int ord, int n, int m, int zero_cols,
                   const ConstVector &da, const ConstVector &db,
                   const ConstVector &dc, const ConstVector &dd,
                   const SylvParams &ps);
  GeneralSylvester(int ord, int n, int m, int zero_cols,
                   const ConstVector &da, const ConstVector &db,
                   const ConstVector &dc, const ConstVector &dd,
                   bool alloc_for_check = false);
  /* construct with provided storage for d */
  GeneralSylvester(int ord, int n, int m, int zero_cols,
                   const ConstVector &da, const ConstVector &db,
                   const ConstVector &dc, Vector &dd,
                   bool alloc_for_check = false);
  GeneralSylvester(int ord, int n, int m, int zero_cols,
                   const ConstVector &da, const ConstVector &db,
                   const ConstVector &dc, Vector &dd,
                   const SylvParams &ps);
  virtual ~GeneralSylvester() = default;
  int
  getM() const
  {
    return c.numRows();
  }
  int
  getN() const
  {
    return a.numRows();
  }
  const double *
  getResult() const
  {
    return d.base();
  }
  const SylvParams &
  getParams() const
  {
    return pars;
  }
  SylvParams &
  getParams()
  {
    return pars;
  }
  void solve();
  void check(const ConstVector &ds);
private:
  void init();
};

#endif /* GENERAL_SYLVESTER_H */
