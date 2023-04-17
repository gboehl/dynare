/*
 * Copyright © 2005 Ondra Kamenik
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

#include "smolyak.hh"
#include "symmetry.hh"

#include <iostream>
#include <iomanip>

/* This constructs a beginning of ‘isum’ summand in ‘smolq’. We must be careful
   here, since ‘isum’ can be past-the-end, so no reference to vectors in
   ‘smolq’ by ‘isum’ must be done in this case. */

smolpit::smolpit(const SmolyakQuadrature &q, unsigned int isum)
  : smolq(q), isummand(isum),
    jseq(q.dimen(), 0),
    sig{q.dimen()},
    p{q.dimen()}
{
  if (isummand < q.numSummands())
    setPointAndWeight();
}

bool
smolpit::operator==(const smolpit &spit) const
{
  return &smolq == &spit.smolq && isummand == spit.isummand && jseq == spit.jseq;
}

/* We first try to increase index within the current summand. If we are at
   maximum, we go to a subsequent summand. Note that in this case all indices
   in ‘jseq’ will be zero, so no change is needed. */

smolpit &
smolpit::operator++()
{
  const IntSequence &levpts = smolq.levpoints[isummand];
  int i = smolq.dimen()-1;
  jseq[i]++;
  while (i >= 0 && jseq[i] == levpts[i])
    {
      jseq[i] = 0;
      i--;
      if (i >= 0)
        jseq[i]++;
    }
  sig.signalAfter(std::max(i, 0));

  if (i < 0)
    isummand++;

  if (isummand < smolq.numSummands())
    setPointAndWeight();

  return *this;
}

/* Here we set the point coordinates according to ‘jseq’ and
   ‘isummand’. Also the weight is set here. */

void
smolpit::setPointAndWeight()
{
  // todo: raise if isummand ≥ smolq.numSummands()
  int l = smolq.level;
  int d = smolq.dimen();
  int sumk = (smolq.levels[isummand]).sum();
  int m1exp = l + d - sumk - 1;
  w = (2*(m1exp/2) == m1exp) ? 1.0 : -1.0;
  w *= PascalTriangle::noverk(d-1, sumk-l);
  for (int i = 0; i < d; i++)
    {
      int ki = (smolq.levels[isummand])[i];
      p[i] = (smolq.uquad).point(ki, jseq[i]);
      w *= (smolq.uquad).weight(ki, jseq[i]);
    }
}

/* Debug print. */
void
smolpit::print() const
{
  auto ff = std::cout.flags();
  std::cout << "isum=" << std::left << std::setw(3) << isummand << std::right << ": [";
  for (int i = 0; i < smolq.dimen(); i++)
    std::cout << std::setw(2) << (smolq.levels[isummand])[i] << ' ';
  std::cout << "] j=[";
  for (int i = 0; i < smolq.dimen(); i++)
    std::cout << std::setw(2) << jseq[i] << ' ';
  std::cout << std::showpos << std::fixed << std::setprecision(3)
            << "] " << std::setw(4) << w << "*(";
  for (int i = 0; i < smolq.dimen()-1; i++)
    std::cout << std::setw(4) << p[i] << ' ';
  std::cout << std::setw(4) << p[smolq.dimen()-1] << ')' << std::endl;
  std::cout.flags(ff);
}

/* Here is the constructor of SmolyakQuadrature. We have to setup ‘levels’,
   ‘levpoints’ and ‘cumevals’. We have to go through all d-dimensional
   sequences k, such that l≤|k|≤l+d−1 and all kᵢ are positive integers. This is
   equivalent to going through all k such that l−d≤|k|≤l−1 and all kᵢ are
   non-negative integers. This is equivalent to going through d+1 dimensional
   sequences (k,x) such that |(k,x)|=l−1 and x=0,…,d−1. The resulting sequence
   of positive integers is obtained by adding 1 to all kᵢ. */

SmolyakQuadrature::SmolyakQuadrature(int d, int l, const OneDQuadrature &uq)
  : QuadratureImpl<smolpit>(d), level(l), uquad(uq)
{
  // TODO: check l>1, l≥d
  // TODO: check l≥uquad.miLevel(), l≤uquad.maxLevel()
  int cum = 0;
  for (const auto &si : SymmetrySet(l-1, d+1))
    {
      if (si[d] <= d-1)
        {
          IntSequence lev(si, 0, d);
          lev.add(1);
          levels.push_back(lev);
          IntSequence levpts(d);
          for (int i = 0; i < d; i++)
            levpts[i] = uquad.numPoints(lev[i]);
          levpoints.push_back(levpts);
          cum += levpts.mult();
          cumevals.push_back(cum);
        }
    }
}

/* Here we return a number of evalutions of the quadrature for the given level.
   If the given level is the current one, we simply return the maximum
   cumulative number of evaluations. Otherwise we call costly
   calcNumEvaluations() method. */

int
SmolyakQuadrature::numEvals(int l) const
{
  if (l != level)
    return calcNumEvaluations(l);
  else
    return cumevals[numSummands()-1];
}

/* This divides all the evaluations to ‘tn’ approximately equal groups, and
   returns the beginning of the specified group ‘ti’. The granularity of
   divisions are summands as listed by ‘levels’. */

smolpit
SmolyakQuadrature::begin(int ti, int tn, int l) const
{
  // TODO: raise is level≠l
  if (ti == tn)
    return smolpit(*this, numSummands());

  int totevals = cumevals[numSummands()-1];
  int evals = (totevals*ti)/tn;
  unsigned int isum = 0;
  while (isum+1 < numSummands() && cumevals[isum+1] < evals)
    isum++;
  return smolpit(*this, isum);
}

/* This is the same in a structure as SmolyakQuadrature constructor. We have to
   go through all summands and calculate a number of evaluations in each
   summand. */

int
SmolyakQuadrature::calcNumEvaluations(int lev) const
{
  int cum = 0;
  for (const auto &si : SymmetrySet(lev-1, dim+1))
    {
      if (si[dim] <= dim-1)
        {
          IntSequence lev(si, 0, dim);
          lev.add(1);
          IntSequence levpts(dim);
          for (int i = 0; i < dim; i++)
            levpts[i] = uquad.numPoints(lev[i]);
          cum += levpts.mult();
        }
    }
  return cum;
}

/* This returns a maximum level such that the number of evaluations is less
   than the given number. */

void
SmolyakQuadrature::designLevelForEvals(int max_evals, int &lev, int &evals) const
{
  int last_evals;
  evals = 1;
  lev = 1;
  do
    {
      lev++;
      last_evals = evals;
      evals = calcNumEvaluations(lev);
    }
  while (lev < uquad.numLevels() && evals <= max_evals);
  lev--;
  evals = last_evals;
}
