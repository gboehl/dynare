// $Id: dynare3.h 1764 2008-03-31 14:30:55Z kamenik $
// Copyright 2005, Ondra Kamenik

#ifndef DYNARE3_H
#define DYNARE3_H

#include "../tl/cc/t_container.hh"
#include "../tl/cc/sparse_tensor.hh"
#include "../kord/decision_rule.hh"
#include "../kord/dynamic_model.hh"

#include "dynare_model.hh"
#include "nlsolve.hh"

#include <vector>

#include <matio.h>

class Dynare;

class DynareNameList : public NameList
{
  std::vector<const char *> names;
public:
  DynareNameList(const Dynare &dynare);
  int
  getNum() const override
  {
    return (int) names.size();
  }
  const char *
  getName(int i) const override
  {
    return names[i];
  }
  /** This for each string of the input vector calculates its index
   * in the names. And returns the resulting vector of indices. If
   * the name cannot be found, then an exception is raised. */
  std::vector<int> selectIndices(const std::vector<const char *> &ns) const;
};

class DynareExogNameList : public NameList
{
  std::vector<const char *> names;
public:
  DynareExogNameList(const Dynare &dynare);
  int
  getNum() const override
  {
    return (int) names.size();
  }
  const char *
  getName(int i) const override
  {
    return names[i];
  }
};

class DynareStateNameList : public NameList
{
  std::vector<const char *> names;
public:
  DynareStateNameList(const Dynare &dynare, const DynareNameList &dnl,
                      const DynareExogNameList &denl);
  int
  getNum() const override
  {
    return (int) names.size();
  }
  const char *
  getName(int i) const override
  {
    return names[i];
  }
};

// The following only implements DynamicModel with help of ogdyn::DynareModel

class DynareJacobian;
class Dynare : public DynamicModel
{
  friend class DynareNameList;
  friend class DynareExogNameList;
  friend class DynareStateNameList;
  friend class DynareJacobian;
  Journal &journal;
  ogdyn::DynareModel *model;
  Vector *ysteady;
  TensorContainer<FSSparseTensor> md;
  DynareNameList *dnl;
  DynareExogNameList *denl;
  DynareStateNameList *dsnl;
  ogp::FormulaEvaluator *fe;
  ogp::FormulaDerEvaluator *fde;
  const double ss_tol;
public:
  /** Parses the given model file and uses the given order to
   * override order from the model file (if it is != -1). */
  Dynare(const char *modname, int ord, double sstol, Journal &jr);
  /** Parses the given equations with explicitly given names. */
  Dynare(const char **endo, int num_endo,
         const char **exo, int num_exo,
         const char **par, int num_par,
         const char *equations, int len, int ord,
         double sstol, Journal &jr);
  /** Makes a deep copy of the object. */
  Dynare(const Dynare &dyn);
  DynamicModel *
  clone() const override
  {
    return new Dynare(*this);
  }
  
  ~Dynare() override;
  int
  nstat() const override
  {
    return model->getAtoms().nstat();
  }
  int
  nboth() const override
  {
    return model->getAtoms().nboth();
  }
  int
  npred() const override
  {
    return model->getAtoms().npred();
  }
  int
  nforw() const override
  {
    return model->getAtoms().nforw();
  }
  int
  nexog() const override
  {
    return model->getAtoms().nexo();
  }
  int
  nys() const
  {
    return model->getAtoms().nys();
  }
  int
  nyss() const
  {
    return model->getAtoms().nyss();
  }
  int
  ny() const
  {
    return model->getAtoms().ny();
  }
  int
  order() const override
  {
    return model->getOrder();
  }

  const NameList &
  getAllEndoNames() const override
  {
    return *dnl;
  }
  const NameList &
  getStateNames() const override
  {
    return *dsnl;
  }
  const NameList &
  getExogNames() const override
  {
    return *denl;
  }

  TwoDMatrix &
  getVcov()
  {
    return model->getVcov();
  }
  const TwoDMatrix &
  getVcov() const override
  {
    return model->getVcov();
  }
  Vector &
  getParams()
  {
    return model->getParams();
  }
  const Vector &
  getParams() const
  {
    return model->getParams();
  }
  void
  setInitOuter(const Vector &x)
  {
    model->setInitOuter(x);
  }

  const TensorContainer<FSSparseTensor> &
  getModelDerivatives() const override
  {
    return md;
  }
  const Vector &
  getSteady() const override
  {
    return *ysteady;
  }
  Vector &
  getSteady() override
  {
    return *ysteady;
  }
  const ogdyn::DynareModel &
  getModel() const
  {
    return *model;
  }

  // here is true public interface
  void solveDeterministicSteady(Vector &steady);
  void
  solveDeterministicSteady() override
  {
    solveDeterministicSteady(*ysteady);
  }
  void evaluateSystem(Vector &out, const ConstVector &yy, const Vector &xx) override;
  void evaluateSystem(Vector &out, const ConstVector &yym, const ConstVector &yy,
                      const ConstVector &yyp, const Vector &xx) override;
  void calcDerivatives(const Vector &yy, const Vector &xx);
  void calcDerivativesAtSteady() override;

  void writeMat(mat_t *fd, const char *prefix) const;
  void writeDump(const std::string &basename) const;
private:
  void writeModelInfo(Journal &jr) const;
};

class DynareEvalLoader : public ogp::FormulaEvalLoader, public Vector
{
public:
  DynareEvalLoader(const ogp::FineAtoms &a, Vector &out);
  void
  load(int i, double res) override
  {
    operator[](i) = res;
  }
};

class DynareDerEvalLoader : public ogp::FormulaDerEvalLoader
{
protected:
  const ogp::FineAtoms &atoms;
  TensorContainer<FSSparseTensor> &md;
public:
  DynareDerEvalLoader(const ogp::FineAtoms &a, TensorContainer<FSSparseTensor> &mod_ders,
                      int order);
  void load(int i, int iord, const int *vars, double res) override;
};

class DynareJacobian : public ogu::Jacobian, public ogp::FormulaDerEvalLoader
{
protected:
  Dynare &d;
public:
  DynareJacobian(Dynare &dyn);
  ~DynareJacobian()
  override = default;
  void load(int i, int iord, const int *vars, double res) override;
  void eval(const Vector &in) override;
};

class DynareVectorFunction : public ogu::VectorFunction
{
protected:
  Dynare &d;
public:
  DynareVectorFunction(Dynare &dyn)
    : d(dyn)
  {
  }
  ~DynareVectorFunction()
  override = default;
  int
  inDim() const override
  {
    return d.ny();
  }
  int
  outDim() const override
  {
    return d.ny();
  }
  void eval(const ConstVector &in, Vector &out) override;
};

#endif
