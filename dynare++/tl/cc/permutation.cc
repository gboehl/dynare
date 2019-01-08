// Copyright 2004, Ondra Kamenik

#include "permutation.hh"
#include "tl_exception.hh"

/* This is easy, we simply apply the map in the fashion $s\circ m$.. */

void
Permutation::apply(const IntSequence &src, IntSequence &tar) const
{
  TL_RAISE_IF(src.size() != permap.size() || tar.size() != permap.size(),
              "Wrong sizes of input or output in Permutation::apply");
  for (int i = 0; i < permap.size(); i++)
    tar[i] = src[permap[i]];
}

void
Permutation::apply(IntSequence &tar) const
{
  IntSequence tmp(tar);
  apply(tmp, tar);
}

void
Permutation::inverse()
{
  IntSequence former(permap);
  for (int i = 0; i < size(); i++)
    permap[former[i]] = i;
}

/* Here we find a number of trailing indices which are identical with
   the permutation. */

int
Permutation::tailIdentity() const
{
  int i = permap.size();
  while (i > 0 && permap[i-1] == i-1)
    i--;
  return permap.size() - i;
}

/* This calculates a map which corresponds to sorting in the following
   sense: $(\hbox{sorted }s)\circ m = s$, where $s$ is a given sequence.

   We go through |s| and find an the same item in sorted |s|. We
   construct the |permap| from the found pair of indices. We have to be
   careful, to not assign to two positions in |s| the same position in
   sorted |s|, so we maintain a bitmap |flag|, in which we remember
   indices from the sorted |s| already assigned. */

void
Permutation::computeSortingMap(const IntSequence &s)
{
  IntSequence srt(s);
  srt.sort();
  IntSequence flags(s.size(), 0);

  for (int i = 0; i < s.size(); i++)
    {
      int j = 0;
      while (j < s.size() && (flags[j] || srt[j] != s[i]))
        j++;
      TL_RAISE_IF(j == s.size(),
                  "Internal algorithm error in Permutation::computeSortingMap");
      flags[j] = 1;
      permap[i] = j;
    }
}

PermutationSet::PermutationSet()
  : order(1), size(1), pers(new const Permutation *[size])
{
  pers[0] = new Permutation(1);
}

PermutationSet::PermutationSet(const PermutationSet &sp, int n)
  : order(n), size(n*sp.size),
    pers(new const Permutation *[size])
{
  for (int i = 0; i < size; i++)
    pers[i] = NULL;

  TL_RAISE_IF(n != sp.order+1,
              "Wrong new order in PermutationSet constructor");

  int k = 0;
  for (int i = 0; i < sp.size; i++)
    {
      for (int j = 0; j < order; j++, k++)
        {
          pers[k] = new Permutation(*(sp.pers[i]), j);
        }
    }
}

PermutationSet::~PermutationSet()
{
  for (int i = 0; i < size; i++)
    if (pers[i])
      delete pers[i];
  delete [] pers;
}

vector<const Permutation *>
PermutationSet::getPreserving(const IntSequence &s) const
{
  TL_RAISE_IF(s.size() != order,
              "Wrong sequence length in PermutationSet::getPreserving");

  vector<const Permutation *> res;
  IntSequence tmp(s.size());
  for (int i = 0; i < size; i++)
    {
      pers[i]->apply(s, tmp);
      if (s == tmp)
        {
          res.push_back(pers[i]);
        }
    }

  return res;
}

PermutationBundle::PermutationBundle(int nmax)
{
  nmax = max(nmax, 1);
  generateUpTo(nmax);
}

PermutationBundle::~PermutationBundle()
{
  for (unsigned int i = 0; i < bundle.size(); i++)
    delete bundle[i];
}

const PermutationSet &
PermutationBundle::get(int n) const
{
  if (n > (int) (bundle.size()) || n < 1)
    {
      TL_RAISE("Permutation set not found in PermutationSet::get");
      return *(bundle[0]);
    }
  else
    {
      return *(bundle[n-1]);
    }
}

void
PermutationBundle::generateUpTo(int nmax)
{
  if (bundle.size() == 0)
    bundle.push_back(new PermutationSet());

  int curmax = bundle.size();
  for (int n = curmax+1; n <= nmax; n++)
    {
      bundle.push_back(new PermutationSet(*(bundle.back()), n));
    }
}
