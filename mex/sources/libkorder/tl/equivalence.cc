/*
 * Copyright © 2004 Ondra Kamenik
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

#include "equivalence.hh"
#include "permutation.hh"
#include "tl_exception.hh"

#include <iostream>

int
OrdSequence::operator[](int i) const
{
  TL_RAISE_IF((i < 0 || i >= length()),
              "Index out of range in OrdSequence::operator[]");
  return data[i];
}

/* Here we implement the ordering. It can be changed, or various
   orderings can be used for different problem sizes. We order them
   according to the average, and then according to the first item. */

bool
OrdSequence::operator<(const OrdSequence &s) const
{
  double ta = average();
  double sa = s.average();
  return (ta < sa || ((ta == sa) && (operator[](0) > s[0])));
}

bool
OrdSequence::operator==(const OrdSequence &s) const
{
  if (length() != s.length())
    return false;

  int i = 0;
  while (i < length() && operator[](i) == s[i])
    i++;

  return (i == length());
}

/* The first add() adds a given integer to the class, the second
   iterates through a given sequence and adds everything found in the
   given class. */

void
OrdSequence::add(int i)
{
  auto vit = data.begin();
  while (vit != data.end() && *vit < i)
    ++vit;
  if (vit != data.end() && *vit == i)
    return;
  data.insert(vit, i);
}

void
OrdSequence::add(const OrdSequence &s)
{
  auto vit = s.data.begin();
  while (vit != s.data.end())
    {
      add(*vit);
      ++vit;
    }
}

/* Answers true if a given number is in the class. */
bool
OrdSequence::has(int i) const
{
  auto vit = data.begin();
  while (vit != data.end())
    {
      if (*vit == i)
        return true;
      ++vit;
    }
  return false;
}

/* Return an average of the class. */
double
OrdSequence::average() const
{
  double res = 0;
  for (int i : data)
    res += i;
  TL_RAISE_IF(data.size() == 0,
              "Attempt to take average of empty class in OrdSequence::average");
  return res/data.size();
}

/* Debug print. */
void
OrdSequence::print(const std::string &prefix) const
{
  std::cout << prefix;
  for (int i : data)
    std::cout << i << ' ';
  std::cout << '\n';
}

Equivalence::Equivalence(int num)
  : n(num)
{
  for (int i = 0; i < num; i++)
    {
      OrdSequence s;
      s.add(i);
      classes.push_back(s);
    }
}

Equivalence::Equivalence(int num, const std::string &dummy)
  : n(num)
{
  OrdSequence s;
  for (int i = 0; i < num; i++)
    s.add(i);
  classes.push_back(s);
}

/* Copy constructor that also glues a given couple. */

Equivalence::Equivalence(const Equivalence &e, int i1, int i2)
  : n(e.n),
    classes(e.classes)
{
  auto s1 = find(i1);
  auto s2 = find(i2);
  if (s1 != s2)
    {
      OrdSequence ns(*s1);
      ns.add(*s2);
      classes.erase(s1);
      classes.erase(s2);
      insert(ns);
    }
}

bool
Equivalence::operator==(const Equivalence &e) const
{
  if (!std::operator==(classes, e.classes))
    return false;

  if (n != e.n)
    return false;

  return true;
}

/* Return an iterator pointing to a class having a given integer. */

Equivalence::const_seqit
Equivalence::findHaving(int i) const
{
  auto si = classes.begin();
  while (si != classes.end())
    {
      if (si->has(i))
        return si;
      ++si;
    }
  TL_RAISE_IF(si == classes.end(),
              "Couldn't find equivalence class in Equivalence::findHaving");
  return si;
}

Equivalence::seqit
Equivalence::findHaving(int i)
{
  auto si = classes.begin();
  while (si != classes.end())
    {
      if (si->has(i))
        return si;
      ++si;
    }
  TL_RAISE_IF(si == classes.end(),
              "Couldn't find equivalence class in Equivalence::findHaving");
  return si;
}

/* Find j-th class for a given j. */

Equivalence::const_seqit
Equivalence::find(int j) const
{
  auto si = classes.begin();
  int i = 0;
  while (si != classes.end() && i < j)
    {
      ++si;
      i++;
    }
  TL_RAISE_IF(si == classes.end(),
              "Couldn't find equivalence class in Equivalence::find");
  return si;
}

Equivalence::seqit
Equivalence::find(int j)
{
  auto si = classes.begin();
  int i = 0;
  while (si != classes.end() && i < j)
    {
      ++si;
      i++;
    }
  TL_RAISE_IF(si == classes.end(),
              "Couldn't find equivalence class in Equivalence::find");
  return si;
}

/* Insert a new class yielding the ordering. */
void
Equivalence::insert(const OrdSequence &s)
{
  auto si = classes.begin();
  while (si != classes.end() && *si < s)
    ++si;
  classes.insert(si, s);
}

/* Trace the equivalence into the integer sequence. The classes are in
   some order (described earlier), and items within classes are ordered,
   so this implies, that the data can be linearized. This method
   “prints” them to the sequence. We allow for tracing only a given
   number of classes from the beginning. */

void
Equivalence::trace(IntSequence &out, int num) const
{
  int i = 0;
  int nc = 0;
  for (auto it = begin(); it != end() && nc < num; ++it, ++nc)
    for (int j = 0; j < it->length(); j++, i++)
      {
        TL_RAISE_IF(i >= out.size(),
                    "Wrong size of output sequence in Equivalence::trace");
        out[i] = (*it)[j];
      }
}

void
Equivalence::trace(IntSequence &out, const Permutation &per) const
{
  TL_RAISE_IF(out.size() != n,
              "Wrong size of output sequence in Equivalence::trace");
  TL_RAISE_IF(per.size() != numClasses(),
              "Wrong permutation for permuted Equivalence::trace");
  int i = 0;
  for (int iclass = 0; iclass < numClasses(); iclass++)
    {
      auto itper = find(per.getMap()[iclass]);
      for (int j = 0; j < itper->length(); j++, i++)
        out[i] = (*itper)[j];
    }
}

/* Debug print. */
void
Equivalence::print(const std::string &prefix) const
{
  int i = 0;
  for (auto it = classes.begin();
       it != classes.end();
       ++it, i++)
    {
      std::cout << prefix << "class " << i << ": ";
      it->print("");
    }
}

/* Here we construct a set of all equivalences over n-element set. The
   construction proceeds as follows. We maintain a list of added equivalences.
   At each iteration we pop front of the list, try to add all parents of the
   popped equivalence. This action adds new equivalences to the object and also
   to the added list. We finish the iterations when the added list is empty.

   In the beginning we start with { {0}, {1}, …, {n-1} }. Adding of parents is
   an action which for a given equivalence tries to glue all possible couples
   and checks whether a new equivalence is already in the equivalence set. This
   is not effective, but we will do the construction only ones.

   In this way we breath-first search a lattice of all equivalences. Note
   that the lattice is modular, that is why the result of a construction
   is a list with a property that between two equivalences with the same
   number of classes there are only equivalences with that number of
   classes. Obviously, the list is decreasing in a number of classes
   (since it is constructed by gluing attempts). */

EquivalenceSet::EquivalenceSet(int num)
  : n(num)
{
  std::list<Equivalence> added;
  Equivalence first(n);
  equis.push_back(first);
  addParents(first, added);
  while (!added.empty())
    {
      addParents(added.front(), added);
      added.pop_front();
    }
  if (n > 1)
    equis.emplace_back(n, "");
}

/* This method is used in addParents() and returns true if the object
   already has that equivalence. We trace list of equivalences in reverse
   order since equivalences are ordered in the list from the most
   primitive (nothing equivalent) to maximal (all is equivalent). Since
   we will have much more results of has() method as true, and
   operator==() between equivalences is quick if number of classes
   differ, and in time we will compare with equivalences with less
   classes, then it is more efficient to trace the equivalences from less
   classes to more classes. hence the reverse order. */

bool
EquivalenceSet::has(const Equivalence &e) const
{
  auto rit = equis.rbegin();
  while (rit != equis.rend() && *rit != e)
    ++rit;
  if (rit != equis.rend())
    return true;
  return false;
}

/* Responsibility of this methods is to try to glue all possible
   couples within a given equivalence and add those which are not in the
   list yet. These are added also to the ‘added’ list.

   If number of classes is 2 or 1, we exit, because there is nothing to
   be added. */

void
EquivalenceSet::addParents(const Equivalence &e,
                           std::list<Equivalence> &added)
{
  if (e.numClasses() == 2 || e.numClasses() == 1)
    return;

  for (int i1 = 0; i1 < e.numClasses(); i1++)
    for (int i2 = i1+1; i2 < e.numClasses(); i2++)
      {
        Equivalence ns(e, i1, i2);
        if (!has(ns))
          {
            added.push_back(ns);
            equis.push_back(std::move(ns));
          }
      }
}

/* Debug print. */
void
EquivalenceSet::print(const std::string &prefix) const
{
  int i = 0;
  for (auto it = equis.begin();
       it != equis.end();
       ++it, i++)
    {
      std::cout << prefix << "equivalence " << i << ":(classes "
                << it->numClasses() << ")\n";
      it->print(prefix + "    ");
    }
}

/* Construct the bundle. nmax is a maximum size of underlying set. */
EquivalenceBundle::EquivalenceBundle(int nmax)
{
  nmax = std::max(nmax, 1);
  generateUpTo(nmax);
}

/* Remember, that the first item is EquivalenceSet(1). */
const EquivalenceSet &
EquivalenceBundle::get(int n) const
{
  TL_RAISE_IF(n > static_cast<int>(bundle.size()) || n < 1,
              "Equivalence set not found in EquivalenceBundle::get");
  return bundle[n-1];
}

/* Get ‘curmax’ which is a maximum size in the bundle, and generate for
   all sizes from curmax+1 up to nmax. */

void
EquivalenceBundle::generateUpTo(int nmax)
{
  int curmax = bundle.size();
  for (int i = curmax+1; i <= nmax; i++)
    bundle.emplace_back(i);
}
