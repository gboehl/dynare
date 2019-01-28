// Copyright 2004, Ondra Kamenik

#include "kord_exception.hh"
#include "decision_rule.hh"
#include "dynamic_model.hh"

#include "SymSchurDecomp.hh"

#include <dynlapack.h>

#include <limits>
#include <utility>
#include <memory>

template <>
int DRFixPoint<KOrder::fold>::max_iter = 10000;
template <>
int DRFixPoint<KOrder::unfold>::max_iter = 10000;
template <>
double DRFixPoint<KOrder::fold>::tol = 1.e-10;
template <>
double DRFixPoint<KOrder::unfold>::tol = 1.e-10;
template <>
int DRFixPoint<KOrder::fold>::max_newton_iter = 50;
template <>
int DRFixPoint<KOrder::unfold>::max_newton_iter = 50;
template <>
int DRFixPoint<KOrder::fold>::newton_pause = 100;
template <>
int DRFixPoint<KOrder::unfold>::newton_pause = 100;

// |FoldDecisionRule| conversion from |UnfoldDecisionRule|
FoldDecisionRule::FoldDecisionRule(const UnfoldDecisionRule &udr)
  : DecisionRuleImpl<KOrder::fold>(ctraits<KOrder::fold>::Tpol(udr.nrows(), udr.nvars()),
                                   udr.ypart, udr.nu, udr.ysteady)
{
  for (const auto & it : udr)
    {
      insert(new ctraits<KOrder::fold>::Ttensym(*(it.second)));
    }
}

// |UnfoldDecisionRule| conversion from |FoldDecisionRule|
UnfoldDecisionRule::UnfoldDecisionRule(const FoldDecisionRule &fdr)
  : DecisionRuleImpl<KOrder::unfold>(ctraits<KOrder::unfold>::Tpol(fdr.nrows(), fdr.nvars()),
                                     fdr.ypart, fdr.nu, fdr.ysteady)
{
  for (const auto & it : fdr)
    {
      insert(new ctraits<KOrder::unfold>::Ttensym(*(it.second)));
    }
}

SimResults::~SimResults()
{
  for (int i = 0; i < getNumSets(); i++)
    {
      delete data[i];
      delete shocks[i];
    }
}

/* This runs simulations with an output to journal file. Note that we
   report how many simulations had to be thrown out due to Nan or Inf. */

void
SimResults::simulate(int num_sim, const DecisionRule &dr, const Vector &start,
                     const TwoDMatrix &vcov, Journal &journal)
{
  JournalRecordPair paa(journal);
  paa << "Performing " << num_sim << " stochastic simulations for "
      << num_per << " periods burning " << num_burn << " initial periods"  << endrec;
  simulate(num_sim, dr, start, vcov);
  int thrown = num_sim - data.size();
  if (thrown > 0)
    {
      JournalRecord rec(journal);
      rec << "I had to throw " << thrown << " simulations away due to Nan or Inf" << endrec;
    }
}

/* This runs a given number of simulations by creating
   |SimulationWorker| for each simulation and inserting them to the
   thread group. */

void
SimResults::simulate(int num_sim, const DecisionRule &dr, const Vector &start,
                     const TwoDMatrix &vcov)
{
  std::vector<RandomShockRealization> rsrs;
  rsrs.reserve(num_sim);

  sthread::detach_thread_group gr;
  for (int i = 0; i < num_sim; i++)
    {
      RandomShockRealization sr(vcov, system_random_generator.int_uniform());
      rsrs.push_back(sr);
      gr.insert(std::make_unique<SimulationWorker>(*this, dr, DecisionRule::horner,
                                                   num_per+num_burn, start, rsrs.back()));
    }
  gr.run();
}

/* This adds the data with the realized shocks. It takes only periods
   which are not to be burnt. If the data is not finite, the both data
   and shocks are thrown away. */

bool
SimResults::addDataSet(TwoDMatrix *d, ExplicitShockRealization *sr)
{
  KORD_RAISE_IF(d->nrows() != num_y,
                "Incompatible number of rows for SimResults::addDataSets");
  KORD_RAISE_IF(d->ncols() != num_per+num_burn,
                "Incompatible number of cols for SimResults::addDataSets");
  bool ret = false;
  if (d->isFinite())
    {
      data.push_back(new TwoDMatrix((const TwoDMatrix &) (*d), num_burn, num_per));
      shocks.push_back(new ExplicitShockRealization(
                                                    ConstTwoDMatrix(sr->getShocks(), num_burn, num_per)));
      ret = true;
    }

  delete d;
  delete sr;
  return ret;
}

void
SimResults::writeMat(const char *base, const char *lname) const
{
  char matfile_name[100];
  sprintf(matfile_name, "%s.mat", base);
  mat_t *matfd = Mat_Create(matfile_name, nullptr);
  if (matfd != nullptr)
    {
      writeMat(matfd, lname);
      Mat_Close(matfd);
    }
}

/* This save the results as matrices with given prefix and with index
   appended. If there is only one matrix, the index is not appended. */

void
SimResults::writeMat(mat_t *fd, const char *lname) const
{
  char tmp[100];
  for (int i = 0; i < getNumSets(); i++)
    {
      if (getNumSets() > 1)
        sprintf(tmp, "%s_data%d", lname, i+1);
      else
        sprintf(tmp, "%s_data", lname);
      ConstTwoDMatrix m(*(data[i]));
      m.writeMat(fd, tmp);
    }
}

void
SimResultsStats::simulate(int num_sim, const DecisionRule &dr,
                          const Vector &start,
                          const TwoDMatrix &vcov, Journal &journal)
{
  SimResults::simulate(num_sim, dr, start, vcov, journal);
  {
    JournalRecordPair paa(journal);
    paa << "Calculating means from the simulations." << endrec;
    calcMean();
  }
  {
    JournalRecordPair paa(journal);
    paa << "Calculating covariances from the simulations." << endrec;
    calcVcov();
  }
}

/* Here we do not save the data itself, we save only mean and vcov. */
void
SimResultsStats::writeMat(mat_t *fd, const char *lname) const
{
  char tmp[100];
  sprintf(tmp, "%s_mean", lname);
  ConstTwoDMatrix m(num_y, 1, mean);
  m.writeMat(fd, tmp);
  sprintf(tmp, "%s_vcov", lname);
  ConstTwoDMatrix(vcov).writeMat(fd, tmp);
}

void
SimResultsStats::calcMean()
{
  mean.zeros();
  if (data.size()*num_per > 0)
    {
      double mult = 1.0/data.size()/num_per;
      for (auto & i : data)
        {
          for (int j = 0; j < num_per; j++)
            {
              ConstVector col{i->getCol(j)};
              mean.add(mult, col);
            }
        }
    }
}

void
SimResultsStats::calcVcov()
{
  if (data.size()*num_per > 1)
    {
      vcov.zeros();
      double mult = 1.0/(data.size()*num_per - 1);
      for (auto & i : data)
        {
          const TwoDMatrix &d = *i;
          for (int j = 0; j < num_per; j++)
            {
              for (int m = 0; m < num_y; m++)
                {
                  for (int n = m; n < num_y; n++)
                    {
                      double s = (d.get(m, j)-mean[m])*(d.get(n, j)-mean[n]);
                      vcov.get(m, n) += mult*s;
                      if (m != n)
                        vcov.get(n, m) += mult*s;
                    }
                }
            }
        }
    }
  else
    {
      vcov.infs();
    }
}

void
SimResultsDynamicStats::simulate(int num_sim, const DecisionRule &dr,
                                 const Vector &start,
                                 const TwoDMatrix &vcov, Journal &journal)
{
  SimResults::simulate(num_sim, dr, start, vcov, journal);
  {
    JournalRecordPair paa(journal);
    paa << "Calculating means of the conditional simulations." << endrec;
    calcMean();
  }
  {
    JournalRecordPair paa(journal);
    paa << "Calculating variances of the conditional simulations." << endrec;
    calcVariance();
  }
}

void
SimResultsDynamicStats::writeMat(mat_t *fd, const char *lname) const
{
  char tmp[100];
  sprintf(tmp, "%s_cond_mean", lname);
  ConstTwoDMatrix(mean).writeMat(fd, tmp);
  sprintf(tmp, "%s_cond_variance", lname);
  ConstTwoDMatrix(variance).writeMat(fd, tmp);
}

void
SimResultsDynamicStats::calcMean()
{
  mean.zeros();
  if (data.size() > 0)
    {
      double mult = 1.0/data.size();
      for (int j = 0; j < num_per; j++)
        {
          Vector meanj{mean.getCol(j)};
          for (auto & i : data)
            {
              ConstVector col{i->getCol(j)};
              meanj.add(mult, col);
            }
        }
    }
}

void
SimResultsDynamicStats::calcVariance()
{
  if (data.size() > 1)
    {
      variance.zeros();
      double mult = 1.0/(data.size()-1);
      for (int j = 0; j < num_per; j++)
        {
          ConstVector meanj{mean.getCol(j)};
          Vector varj{variance.getCol(j)};
          for (auto & i : data)
            {
              Vector col{i->getCol(j)};
              col.add(-1.0, meanj);
              for (int k = 0; k < col.length(); k++)
                col[k] = col[k]*col[k];
              varj.add(mult, col);
            }
        }
    }
  else
    {
      variance.infs();
    }
}

void
SimResultsIRF::simulate(const DecisionRule &dr, Journal &journal)
{
  JournalRecordPair paa(journal);
  paa << "Performing " << control.getNumSets() << " IRF simulations for "
      << num_per << " periods; shock=" << ishock << ", impulse=" << imp << endrec;
  simulate(dr);
  int thrown = control.getNumSets() - data.size();
  if (thrown > 0)
    {
      JournalRecord rec(journal);
      rec << "I had to throw " << thrown
          << " simulations away due to Nan or Inf" << endrec;
    }
  calcMeans();
  calcVariances();
}

void
SimResultsIRF::simulate(const DecisionRule &dr)
{
  sthread::detach_thread_group gr;
  for (int idata = 0; idata < control.getNumSets(); idata++)
    gr.insert(std::make_unique<SimulationIRFWorker>(*this, dr, DecisionRule::horner,
                                                    num_per, idata, ishock, imp));
  gr.run();
}

void
SimResultsIRF::calcMeans()
{
  means.zeros();
  if (data.size() > 0)
    {
      for (auto & i : data)
        means.add(1.0, *i);
      means.mult(1.0/data.size());
    }
}

void
SimResultsIRF::calcVariances()
{
  if (data.size() > 1)
    {
      variances.zeros();
      for (auto & i : data)
        {
          TwoDMatrix d((const TwoDMatrix &)(*i));
          d.add(-1.0, means);
          for (int j = 0; j < d.nrows(); j++)
            for (int k = 0; k < d.ncols(); k++)
              variances.get(j, k) += d.get(j, k)*d.get(j, k);
          d.mult(1.0/(data.size()-1));
        }
    }
  else
    {
      variances.infs();
    }
}

void
SimResultsIRF::writeMat(mat_t *fd, const char *lname) const
{
  char tmp[100];
  sprintf(tmp, "%s_mean", lname);
  means.writeMat(fd, tmp);
  sprintf(tmp, "%s_var", lname);
  variances.writeMat(fd, tmp);
}

void
RTSimResultsStats::simulate(int num_sim, const DecisionRule &dr, const Vector &start,
                            const TwoDMatrix &v, Journal &journal)
{
  JournalRecordPair paa(journal);
  paa << "Performing " << num_sim << " real-time stochastic simulations for "
      << num_per << " periods" << endrec;
  simulate(num_sim, dr, start, v);
  mean = nc.getMean();
  mean.add(1.0, dr.getSteady());
  nc.getVariance(vcov);
  if (thrown_periods > 0)
    {
      JournalRecord rec(journal);
      rec << "I had to throw " << thrown_periods << " periods away due to Nan or Inf" << endrec;
      JournalRecord rec1(journal);
      rec1 << "This affected " << incomplete_simulations << " out of "
           << num_sim << " simulations" << endrec;
    }
}

void
RTSimResultsStats::simulate(int num_sim, const DecisionRule &dr, const Vector &start,
                            const TwoDMatrix &vcov)
{
  std::vector<RandomShockRealization> rsrs;
  rsrs.reserve(num_sim);

  sthread::detach_thread_group gr;
  for (int i = 0; i < num_sim; i++)
    {
      RandomShockRealization sr(vcov, system_random_generator.int_uniform());
      rsrs.push_back(sr);
      gr.insert(std::make_unique<RTSimulationWorker>(*this, dr, DecisionRule::horner,
                                                     num_per, start, rsrs.back()));
    }
  gr.run();
}

void
RTSimResultsStats::writeMat(mat_t *fd, const char *lname)
{
  char tmp[100];
  sprintf(tmp, "%s_rt_mean", lname);
  ConstTwoDMatrix m(nc.getDim(), 1, mean);
  m.writeMat(fd, tmp);
  sprintf(tmp, "%s_rt_vcov", lname);
  ConstTwoDMatrix(vcov).writeMat(fd, tmp);
}

IRFResults::IRFResults(const DynamicModel &mod, const DecisionRule &dr,
                       const SimResults &control, vector<int> ili,
                       Journal &journal)
  : model(mod), irf_list_ind(std::move(ili))
{
  int num_per = control.getNumPer();
  JournalRecordPair pa(journal);
  pa << "Calculating IRFs against control for " << (int) irf_list_ind.size() << " shocks and for "
     << num_per << " periods" << endrec;
  const TwoDMatrix &vcov = mod.getVcov();
  for (int ishock : irf_list_ind)
    {
      double stderror = sqrt(vcov.get(ishock, ishock));
      irf_res.push_back(new SimResultsIRF(control, model.numeq(), num_per,
                                          ishock, stderror));
      irf_res.push_back(new SimResultsIRF(control, model.numeq(), num_per,
                                          ishock, -stderror));
    }

  for (unsigned int ii = 0; ii < irf_list_ind.size(); ii++)
    {
      irf_res[2*ii]->simulate(dr, journal);
      irf_res[2*ii+1]->simulate(dr, journal);
    }
}

IRFResults::~IRFResults()
{
  for (auto & irf_re : irf_res)
    delete irf_re;
}

void
IRFResults::writeMat(mat_t *fd, const char *prefix) const
{
  for (unsigned int i = 0; i < irf_list_ind.size(); i++)
    {
      char tmp[100];
      int ishock = irf_list_ind[i];
      const char *shockname = model.getExogNames().getName(ishock);
      sprintf(tmp, "%s_irfp_%s", prefix, shockname);
      irf_res[2*i]->writeMat(fd, tmp);
      sprintf(tmp, "%s_irfm_%s", prefix, shockname);
      irf_res[2*i+1]->writeMat(fd, tmp);
    }
}

void
SimulationWorker::operator()()
{
  auto *esr = new ExplicitShockRealization(sr, np);
  TwoDMatrix *m = dr.simulate(em, np, st, *esr);
  {
    sthread::synchro syn(&res, "simulation");
    res.addDataSet(m, esr);
  }
}

/* Here we create a new instance of |ExplicitShockRealization| of the
   corresponding control, add the impulse, and simulate. */

void
SimulationIRFWorker::operator()()
{
  auto *esr
    = new ExplicitShockRealization(res.control.getShocks(idata));
  esr->addToShock(ishock, 0, imp);
  const TwoDMatrix &data = res.control.getData(idata);
  Vector st{data.getCol(res.control.getNumBurn())};
  TwoDMatrix *m = dr.simulate(em, np, st, *esr);
  m->add(-1.0, res.control.getData(idata));
  {
    sthread::synchro syn(&res, "simulation");
    res.addDataSet(m, esr);
  }
}

void
RTSimulationWorker::operator()()
{
  NormalConj nc(res.nc.getDim());
  const PartitionY &ypart = dr.getYPart();
  int nu = dr.nexog();
  const Vector &ysteady = dr.getSteady();

  // initialize vectors and subvectors for simulation
  Vector dyu(ypart.nys()+nu);
  ConstVector ystart_pred(ystart, ypart.nstat, ypart.nys());
  ConstVector ysteady_pred(ysteady, ypart.nstat, ypart.nys());
  Vector dy(dyu, 0, ypart.nys());
  Vector u(dyu, ypart.nys(), nu);
  Vector y(nc.getDim());
  ConstVector ypred(y, ypart.nstat, ypart.nys());

  // simulate the first real-time period
  int ip = 0;
  dy = ystart_pred;
  dy.add(-1.0, ysteady_pred);
  sr.get(ip, u);
  dr.eval(em, y, dyu);
  if (ip >= res.num_burn)
    nc.update(y);

  // simulate other real-time periods@>=
  while (y.isFinite() && ip < res.num_burn + res.num_per)
    {
      ip++;
      dy = ypred;
      sr.get(ip, u);
      dr.eval(em, y, dyu);
      if (ip >= res.num_burn)
        nc.update(y);
    }
  {
    sthread::synchro syn(&res, "rtsimulation");
    res.nc.update(nc);
    if (res.num_per-ip > 0)
      {
        res.incomplete_simulations++;
        res.thrown_periods += res.num_per-ip;
      }
  }
}

/* This calculates factorization $FF^T=V$ in the Cholesky way. It does
   not work for semidefinite matrices. */

void
RandomShockRealization::choleskyFactor(const ConstTwoDMatrix &v)
{
  factor = v;
  lapack_int rows = factor.nrows(), lda = factor.getLD();
  for (int i = 0; i < rows; i++)
    for (int j = i+1; j < rows; j++)
      factor.get(i, j) = 0.0;
  lapack_int info;

  dpotrf("L", &rows, factor.base(), &lda, &info);
  KORD_RAISE_IF(info != 0,
                "Info!=0 in RandomShockRealization::choleskyFactor");
}

/* This calculates $FF^T=V$ factorization by symmetric Schur
   decomposition. It works for semidifinite matrices. */

void
RandomShockRealization::schurFactor(const ConstTwoDMatrix &v)
{
  SymSchurDecomp ssd(v);
  ssd.getFactor(factor);
}

void
RandomShockRealization::get(int n, Vector &out)
{
  KORD_RAISE_IF(out.length() != numShocks(),
                "Wrong length of out vector in RandomShockRealization::get");
  Vector d(out.length());
  for (int i = 0; i < d.length(); i++)
    {
      d[i] = mtwister.normal();
    }
  out.zeros();
  factor.multaVec(out, ConstVector(d));
}

ExplicitShockRealization::ExplicitShockRealization(ShockRealization &sr,
                                                   int num_per)
  : shocks(sr.numShocks(), num_per)
{
  for (int j = 0; j < num_per; j++)
    {
      Vector jcol{shocks.getCol(j)};
      sr.get(j, jcol);
    }
}

void
ExplicitShockRealization::get(int n, Vector &out)
{
  KORD_RAISE_IF(out.length() != numShocks(),
                "Wrong length of out vector in ExplicitShockRealization::get");
  int i = n % shocks.ncols();
  ConstVector icol{shocks.getCol(i)};
  out = icol;
}

void
ExplicitShockRealization::addToShock(int ishock, int iper, double val)
{
  KORD_RAISE_IF(ishock < 0 || ishock > numShocks(),
                "Wrong index of shock in ExplicitShockRealization::addToShock");
  int j = iper % shocks.ncols();
  shocks.get(ishock, j) += val;
}

void
GenShockRealization::get(int n, Vector &out)
{
  KORD_RAISE_IF(out.length() != numShocks(),
                "Wrong length of out vector in GenShockRealization::get");
  ExplicitShockRealization::get(n, out);
  Vector r(numShocks());
  RandomShockRealization::get(n, r);
  for (int j = 0; j < numShocks(); j++)
    if (!isfinite(out[j]))
      out[j] = r[j];
}
