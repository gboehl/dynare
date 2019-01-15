/* $Header: /var/lib/cvs/dynare_cpp/sylv/cc/SchurDecomp.h,v 1.1.1.1 2004/06/04 13:00:44 kamenik Exp $ */

/* Tag $Name:  $ */

#ifndef SCHUR_DECOMP_H
#define SCHUR_DECOMP_H

#include "SylvMatrix.hh"
#include "QuasiTriangular.hh"

#include <memory>

class QuasiTriangular;
class SchurDecomp
{
  SqSylvMatrix q;
  // Stores t if is owned
  std::unique_ptr<QuasiTriangular> t_storage;
  QuasiTriangular *t;
public:
  SchurDecomp(const SqSylvMatrix &m);
  SchurDecomp(const QuasiTriangular &tr);
  SchurDecomp(QuasiTriangular &tr);
  const SqSylvMatrix &
  getQ() const
  {
    return q;
  }
  const QuasiTriangular &
  getT() const
  {
    return *t;
  }
  SqSylvMatrix &
  getQ()
  {
    return q;
  }
  QuasiTriangular &
  getT()
  {
    return *t;
  }
  virtual int getDim() const;
  virtual ~SchurDecomp() = default;
};

class SchurDecompZero : public SchurDecomp
{
  GeneralMatrix ru; /* right upper matrix */
public:
  SchurDecompZero(const GeneralMatrix &m);
  const GeneralMatrix &
  getRU() const
  {
    return ru;
  }
  int getDim() const override;
  int
  getZeroCols() const
  {
    return ru.numRows();
  }
};

#endif /* SCHUR_DECOMP_H */
