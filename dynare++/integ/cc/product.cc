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
 * along with Dynare.  If not, see <http://www.gnu.org/licenses/>.
 */

#include "product.hh"
#include "symmetry.hh"

#include <iostream>
#include <iomanip>

/* This constructs a product iterator corresponding to index (j0,0,…,0). */

prodpit::prodpit(const ProductQuadrature &q, int j0, int l)
  : prodq(q), level(l), npoints(q.uquad.numPoints(l)),
    jseq(q.dimen(), 0),
    end_flag(false),
    sig{q.dimen()},
    p{q.dimen()}
{
  if (j0 < npoints)
    {
      jseq[0] = j0;
      setPointAndWeight();
    }
  else
    end_flag = true;
}

bool
prodpit::operator==(const prodpit &ppit) const
{
  return &prodq == &ppit.prodq && end_flag == ppit.end_flag && jseq == ppit.jseq;
}

prodpit &
prodpit::operator++()
{
  int i = prodq.dimen()-1;
  jseq[i]++;
  while (i >= 0 && jseq[i] == npoints)
    {
      jseq[i] = 0;
      i--;
      if (i >= 0)
        jseq[i]++;
    }
  sig.signalAfter(std::max(i, 0));

  if (i == -1)
    end_flag = true;

  if (!end_flag)
    setPointAndWeight();

  return *this;
}

/* This calculates the weight and sets point coordinates from the indices. */

void
prodpit::setPointAndWeight()
{
  w = 1.0;
  for (int i = 0; i < prodq.dimen(); i++)
    {
      p[i] = (prodq.uquad).point(level, jseq[i]);
      w *= (prodq.uquad).weight(level, jseq[i]);
    }
}

/* Debug print. */

void
prodpit::print() const
{
  auto ff = std::cout.flags();
  std::cout << "j=[";
  for (int i = 0; i < prodq.dimen(); i++)
    std::cout << std::setw(2) << jseq[i];
  std::cout << std::showpos << std::fixed << std::setprecision(3)
            << "] " << std::setw(4) << w << "*(";
  for (int i = 0; i < prodq.dimen()-1; i++)
    std::cout << std::setw(4) << p[i] << ' ';
  std::cout << std::setw(4) << p[prodq.dimen()-1] << ')' << std::endl;
  std::cout.flags(ff);
}

ProductQuadrature::ProductQuadrature(int d, const OneDQuadrature &uq)
  : QuadratureImpl<prodpit>(d), uquad(uq)
{
  // TODO: check d≥1
}

/* This calls prodpit constructor to return an iterator which points
   approximatelly at ‘ti’-th portion out of ‘tn’ portions. First we find
   out how many points are in the level, and then construct an interator
   (j0,0,…,0) where j0=ti·npoints/tn. */

prodpit
ProductQuadrature::begin(int ti, int tn, int l) const
{
  // TODO: raise if l<dimen()
  // TODO: check l ≤ uquad.numLevels()
  int npoints = uquad.numPoints(l);
  return prodpit(*this, ti*npoints/tn, l);
}

/* This just starts at the first level and goes to a higher level as long as a
   number of evaluations (which is nₖᵈ for k being the level) is less than the
   given number of evaluations. */

void
ProductQuadrature::designLevelForEvals(int max_evals, int &lev, int &evals) const
{
  int last_evals;
  evals = 1;
  lev = 1;
  do
    {
      lev++;
      last_evals = evals;
      evals = numEvals(lev);
    }
  while (lev < uquad.numLevels()-2 && evals < max_evals);
  lev--;
  evals = last_evals;
}
