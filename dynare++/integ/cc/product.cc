// Copyright 2005, Ondra Kamenik

#include "product.hh"
#include "symmetry.hh"

#include <iostream>
#include <iomanip>

/* This constructs a product iterator corresponding to index $(j0,0\ldots,0)$. */

prodpit::prodpit(const ProductQuadrature &q, int j0, int l)
  : prodq(q), level(l), npoints(q.uquad.numPoints(l)),
    jseq{q.dimen(), 0},
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
  // todo: throw if |prodq==NULL| or |jseq==NULL| or |sig==NULL| or |end_flag==true|
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
  // todo: raise if |prodq==NULL| or |jseq==NULL| or |sig==NULL| or
  // |p==NULL| or |end_flag==true|
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
  // todo: check |d>=1|
}

/* This calls |prodpit| constructor to return an iterator which points
   approximatelly at |ti|-th portion out of |tn| portions. First we find
   out how many points are in the level, and then construct an interator
   $(j0,0,\ldots,0)$ where $j0=$|ti*npoints/tn|. */

prodpit
ProductQuadrature::begin(int ti, int tn, int l) const
{
  // todo: raise is |l<dimen()|
  // todo: check |l<=uquad.numLevels()|
  int npoints = uquad.numPoints(l);
  return prodpit(*this, ti*npoints/tn, l);
}

/* This just starts at the first level and goes to a higher level as
   long as a number of evaluations (which is $n_k^d$ for $k$ being the
   level) is less than the given number of evaluations. */

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
