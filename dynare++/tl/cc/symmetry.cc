// Copyright (C) 2004-2011, Ondra Kamenik

#include "symmetry.hh"
#include "permutation.hh"

#include <iostream>

/* Construct symmetry as numbers of successively equal items in the sequence. */

Symmetry::Symmetry(const IntSequence &s)
  : IntSequence(s.getNumDistinct(), 0)
{
  int p = 0;
  if (s.size() > 0)
    operator[](p) = 1;
  for (int i = 1; i < s.size(); i++)
    {
      if (s[i] != s[i-1])
        p++;
      operator[](p)++;
    }
}

/* Find a class of the symmetry containing a given index. */

int
Symmetry::findClass(int i) const
{
  int j = 0;
  int sum = 0;
  do
    {
      sum += operator[](j);
      j++;
    }
  while (j < size() && sum <= i);

  return j-1;
}

/* The symmetry is full if it allows for any permutation of indices. It
   means, that there is at most one non-zero index. */

bool
Symmetry::isFull() const
{
  int count = 0;
  for (int i = 0; i < num(); i++)
    if (operator[](i) != 0)
      count++;
  return count <= 1;
}

/* Construct a symiterator of given dimension, starting from the given
   symmetry. */

symiterator::symiterator(int dim_arg, Symmetry run_arg)
  : dim{dim_arg}, run(std::move(run_arg))
{
}

/* Here we move to the next symmetry. We do so only, if we are not at
   the end. If length is 2, we increase lower index and decrease upper
   index, otherwise we increase the subordinal symmetry. If we got to the
   end, we recreate the subordinal symmetry set and set the subordinal
   iterator to the beginning. */

symiterator &
symiterator::operator++()
{
  if (run[0] == dim)
    run[0]++; // Jump to the past-the-end iterator
  else if (run.size() == 2)
    {
      run[0]++;
      run[1]--;
    }
  else
    {
      symiterator subit{dim-run[0], Symmetry(run, run.size()-1)};
      ++subit;
      if (run[1] == dim-run[0]+1)
        {
          run[0]++;
          run[1] = 0;
          /* subit is equal to the past-the-end iterator, so the range
             2..(size()-1) is already set to 0 */
          run[run.size()-1] = dim-run[0];
        }
    }
  return *this;
}

InducedSymmetries::InducedSymmetries(const Equivalence &e, const Symmetry &s)
{
  for (const auto & i : e)
    emplace_back(s, i);
}

// |InducedSymmetries| permuted constructor code
InducedSymmetries::InducedSymmetries(const Equivalence &e, const Permutation &p,
                                     const Symmetry &s)
{
  for (int i = 0; i < e.numClasses(); i++)
    {
      auto it = e.find(p.getMap()[i]);
      emplace_back(s, *it);
    }
}

/* Debug print. */

void
InducedSymmetries::print() const
{
  std::cout << "Induced symmetries: " << size() << std::endl;
  for (unsigned int i = 0; i < size(); i++)
    operator[](i).print();
}
