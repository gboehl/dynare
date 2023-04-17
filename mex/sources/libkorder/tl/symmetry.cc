/*
 * Copyright © 2004-2011 Ondra Kamenik
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

#include "symmetry.hh"
#include "permutation.hh"
#include "tl_exception.hh"

#include <iostream>

/* This constructs an implied symmetry from a more general symmetry and
   equivalence class. For example, let the general symmetry be y³u² and the
   equivalence class is {0,4} picking up first and fifth variable, we
   calculate symmetry corresponding to the picked variables. These are “yu”.
   Thus the constructed sequence must be (1,1), meaning that we picked one
   y and one u. */

Symmetry::Symmetry(const Symmetry &sy, const OrdSequence &cl)
  : IntSequence(sy.num())
{
  const std::vector<int> &se = cl.getData();
  TL_RAISE_IF(sy.dimen() <= se[se.size()-1],
              "Sequence is not reachable by symmetry in IntSequence()");
  for (int i = 0; i < size(); i++)
    operator[](i) = 0;

  for (int i : se)
    operator[](sy.findClass(i))++;
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
             2…(size()−1) is already set to 0 */
          run[run.size()-1] = dim-run[0];
        }
    }
  return *this;
}

InducedSymmetries::InducedSymmetries(const Equivalence &e, const Symmetry &s)
{
  for (const auto &i : e)
    emplace_back(s, i);
}

// InducedSymmetries permuted constructor code
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
