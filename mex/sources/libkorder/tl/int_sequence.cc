/*
 * Copyright © 2004 Ondra Kamenik
 * Copyright © 2019-2024 Dynare Team
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

#include "int_sequence.hh"
#include "pascal_triangle.hh"
#include "symmetry.hh"
#include "tl_exception.hh"

#include <iostream>
#include <limits>
#include <numeric>

IntSequence
IntSequence::unfold(const Symmetry& sy) const
{
  IntSequence r(sy.dimen());
  int k = 0;
  for (int i = 0; i < sy.num(); i++)
    for (int j = 0; j < sy[i]; j++, k++)
      r[k] = operator[](i);
  return r;
}

Symmetry
IntSequence::getSymmetry() const
{
  Symmetry r(getNumDistinct());
  int p = 0;
  if (size() > 0)
    r[p] = 1;
  for (int i = 1; i < size(); i++)
    {
      if (operator[](i) != operator[](i - 1))
        p++;
      r[p]++;
    }
  return r;
}

IntSequence
IntSequence::insert(int i) const
{
  IntSequence r(size() + 1);
  int j;
  for (j = 0; j < size() && operator[](j) < i; j++)
    r[j] = operator[](j);
  r[j] = i;
  for (; j < size(); j++)
    r[j + 1] = operator[](j);
  return r;
}

IntSequence
IntSequence::insert(int i, int pos) const
{
  TL_RAISE_IF(pos < 0 || pos > size(), "Wrong position for IntSequence::insert()");
  IntSequence r(size() + 1);
  int j;
  for (j = 0; j < pos; j++)
    r[j] = operator[](j);
  r[j] = i;
  for (; j < size(); j++)
    r[j + 1] = operator[](j);
  return r;
}

IntSequence&
IntSequence::operator=(const IntSequence& s)
{
  TL_RAISE_IF(length != s.length, "Wrong length for in-place IntSequence::operator=");
  std::copy_n(s.data, length, data);
  return *this;
}

IntSequence&
IntSequence::operator=(IntSequence&& s)
{
  TL_RAISE_IF(length != s.length, "Wrong length for in-place IntSequence::operator=");
  std::copy_n(s.data, length, data);
  return *this;
}

bool
IntSequence::operator==(const IntSequence& s) const
{
  return std::equal(data, data + length, s.data, s.data + s.length);
}

std::strong_ordering
IntSequence::operator<=>(const IntSequence& s) const
{
  return std::lexicographical_compare_three_way(data, data + length, s.data, s.data + s.length);
}

bool
IntSequence::lessEq(const IntSequence& s) const
{
  TL_RAISE_IF(size() != s.size(), "Sequence with different lengths in IntSequence::lessEq");

  int i = 0;
  while (i < size() && operator[](i) <= s[i])
    i++;
  return (i == size());
}

bool
IntSequence::less(const IntSequence& s) const
{
  TL_RAISE_IF(size() != s.size(), "Sequence with different lengths in IntSequence::less");

  int i = 0;
  while (i < size() && operator[](i) < s[i])
    i++;
  return (i == size());
}

void
IntSequence::sort()
{
  std::sort(data, data + length);
}

void
IntSequence::monotone()
{
  for (int i = 1; i < length; i++)
    if (operator[](i - 1) > operator[](i))
      operator[](i) = operator[](i - 1);
}

void
IntSequence::pmonotone(const Symmetry& s)
{
  int cum = 0;
  for (int i = 0; i < s.num(); i++)
    {
      for (int j = cum + 1; j < cum + s[i]; j++)
        if (operator[](j - 1) > operator[](j))
          operator[](j) = operator[](j - 1);
      cum += s[i];
    }
}

int
IntSequence::sum() const
{
  return std::accumulate(data, data + length, 0);
}

int
IntSequence::mult(int i1, int i2) const
{
  return std::accumulate(data + i1, data + i2, 1, std::multiplies<>());
}

int
IntSequence::getPrefixLength() const
{
  int i = 0;
  while (i + 1 < size() && operator[](i + 1) == operator[](0))
    i++;
  return i + 1;
}

int
IntSequence::getNumDistinct() const
{
  int res = 0;
  if (length > 0)
    res++;
  for (int i = 1; i < length; i++)
    if (operator[](i) != operator[](i - 1))
      res++;
  return res;
}

int
IntSequence::getMax() const
{
  if (length == 0)
    return std::numeric_limits<int>::min();
  return *std::max_element(data, data + length);
}

void
IntSequence::add(int i)
{
  for (int j = 0; j < size(); j++)
    operator[](j) += i;
}

void
IntSequence::add(int f, const IntSequence& s)
{
  TL_RAISE_IF(size() != s.size(), "Wrong sequence length in IntSequence::add");
  for (int j = 0; j < size(); j++)
    operator[](j) += f * s[j];
}

bool
IntSequence::isPositive() const
{
  return std::all_of(data, data + length, [](int x) { return x >= 0; });
}

bool
IntSequence::isConstant() const
{
  if (length < 2)
    return true;
  return std::all_of(data + 1, data + length, [this](int x) { return x == operator[](0); });
}

bool
IntSequence::isSorted() const
{
  return std::is_sorted(data, data + length);
}

/* Debug print. */

void
IntSequence::print() const
{
  std::cout << '[';
  for (int i = 0; i < size(); i++)
    std::cout << operator[](i) << ' ';
  std::cout << ']' << std::endl;
}

/* Here we calculate the multinomial coefficients
   ⎛   a   ⎞
   ⎝b₁,…,bₙ⎠ where a=b₁+…+bₙ

   (the notation gives the name to the function: “n over seq(ence)”)

   See:
    https://en.wikipedia.org/wiki/Binomial_coefficient#Generalization_to_multinomials
    https://en.wikipedia.org/wiki/Multinomial_theorem

   For n=1, the coefficient is equal to 1.
   For n=2, the multinomial coeffs correspond to the binomial coeffs, i.e. the binomial
    ⎛a⎞                             ⎛  a  ⎞
    ⎝b⎠ is equal to the multinomial ⎝b,a−b⎠

   For n≥3, we have the identity
   ⎛   a   ⎞ ⎛b₁+b₂⎞ ⎛      a      ⎞
   ⎝b₁,…,bₙ⎠=⎝  b₁ ⎠·⎝b₁+b₂,b₃,…,bₙ⎠
   (where the first factor on the right hand side is to be interpreted as a
    binomial coefficient)

   This number is exactly a number of unfolded indices corresponding to one
   folded index, where the sequence b₁,…,bₙ is the symmetry of the index. This
   can be easily seen if the multinomial coefficient is interpreted as the
   number of unique permutations of a word, where ‘a’ is the length of the
   word, ‘n’ is the number of distinct letters, and the bᵢ are the number of
   repetitions of each letter. For example, for a symmetry of the form y⁴u²v³,
   we want to compute the number of permutations of the word ‘yyyyuuvvv’.
                                                ⎛  9  ⎞
   This is equal to the multinomial coefficient ⎝4,2,3⎠.
*/
int
IntSequence::noverseq()
{
  if (size() == 0 || size() == 1)
    return 1;
  data[1] += data[0];
  return PascalTriangle::noverk(data[1], data[0]) * IntSequence(*this, 1, size()).noverseq();
}
