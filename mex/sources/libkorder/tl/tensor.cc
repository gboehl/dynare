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

#include "tensor.hh"
#include "pascal_triangle.hh"
#include "tl_exception.hh"
#include "tl_static.hh"

/* Here we increment a given sequence within full symmetry given by ‘nv’, which
   is number of variables in each dimension. The underlying tensor is unfolded,
   so we increase the rightmost by one, and if it is ‘nv’ we zero it and
   increase the next one to the left. */

void
UTensor::increment(IntSequence& v, int nv)
{
  if (v.size() == 0)
    return;
  int i = v.size() - 1;
  v[i]++;
  while (i > 0 && v[i] == nv)
    {
      v[i] = 0;
      v[--i]++;
    }
}

/* This is dual to UTensor::increment(IntSequence& v, int nv). */

void
UTensor::decrement(IntSequence& v, int nv)
{
  if (v.size() == 0)
    return;
  int i = v.size() - 1;
  v[i]--;
  while (i > 0 && v[i] == -1)
    {
      v[i] = nv - 1;
      v[--i]--;
    }
}

/* Here we increment index for general symmetry for unfolded storage. The
   sequence ‘nvmx’ assigns for each coordinate a number of variables. Since the
   storage is unfolded, we do not need information about what variables are
   symmetric, everything necessary is given by ‘nvmx’. */

void
UTensor::increment(IntSequence& v, const IntSequence& nvmx)
{
  if (v.size() == 0)
    return;
  int i = v.size() - 1;
  v[i]++;
  while (i > 0 && v[i] == nvmx[i])
    {
      v[i] = 0;
      v[--i]++;
    }
}

/* This is a dual code to UTensor::increment(IntSequence& v, const IntSequence&
   nvmx). */

void
UTensor::decrement(IntSequence& v, const IntSequence& nvmx)
{
  if (v.size() == 0)
    return;
  int i = v.size() - 1;
  v[i]--;
  while (i > 0 && v[i] == -1)
    {
      v[i] = nvmx[i] - 1;
      v[--i]--;
    }
}

/* Here we return an offset for a given coordinates of unfolded full symmetry
   tensor. This is easy. */

int
UTensor::getOffset(const IntSequence& v, int nv)
{
  int res = 0;
  for (int i = 0; i < v.size(); i++)
    {
      res *= nv;
      res += v[i];
    }
  return res;
}

/* Also easy. */

int
UTensor::getOffset(const IntSequence& v, const IntSequence& nvmx)
{
  int res = 0;
  for (int i = 0; i < v.size(); i++)
    {
      res *= nvmx[i];
      res += v[i];
    }
  return res;
}

/* Decrementing of coordinates of folded index is not that easy. Note that if a
   trailing part of coordinates is (b,a,a,a) (for instance) with b<a, then a
   preceding coordinates are (b,a−1,n−1,n−1), where n is a number of variables
   ‘nv’. So we find the left most element which is equal to the last element,
   decrease it by one, and then set all elements to the right to n−1. */

void
FTensor::decrement(IntSequence& v, int nv)
{
  int i = v.size() - 1;
  while (i > 0 && v[i - 1] == v[i])
    i--;
  v[i]--;
  for (int j = i + 1; j < v.size(); j++)
    v[j] = nv - 1;
}

/* This calculates order of the given index of our ordering of indices. In
   order to understand how it works, let us take number of variables n and
   dimension k, and write down all the possible combinations of indices in our
   ordering. For example for n=4 and k=3, the sequence looks as the following
   (in column-major order):

    000  111  222  333
    001  112  223
    002  113  233
    003  122
    011  123
    012  133
    013
    022
    023
    033

   Now observe, that a number of sequences starting with zero is the same as
   the total number of sequences with the same number of variables but with
   dimension minus one. More generally, if Sₙ,ₖ denotes the number of indices
   of n variables and dimension k, then the number of indices beginning with m
   is exactly Sₙ₋ₘ,ₖ₋₁. This is because m can be subtracted from all items, and
   we obtain the sequence of indices of n−m variables. So we have the formula:

    Sₙ,ₖ = Sₙ,ₖ₋₁ + Sₙ₋₁,ₖ₋₁ + … + S₁,ₖ₋₁

   Now it is easy to calculate the offset of an index of the form (m,…,m). It
   is the sum of all above it, that is Sₙ,ₖ₋₁ + … + Sₙ₋ₘ,ₖ₋₁. Also, since Sₙ,ₖ
   is the number of ways to choose k elements from a set of n elements with
   repetitions allowed, it is well known that:
           ⎛n+k−1⎞
    Sₙ,ₖ = ⎝  k  ⎠.
   Using the above formula, one can therefore calculate the offset of (m,…,m):
    ⎛n+k−1⎞ ⎛n−m+k−1⎞
    ⎝  k  ⎠−⎝   k   ⎠

   The offset of general index (m₁,m₂,…,mₖ) is calculated recursively, since it
   is the offset of (m₁,…,m₁) for n variables plus the offset of
   (m₂−m₁,m₃−m₁,…,mₖ−m₁) for n−m₁ variables. */

int
FTensor::getOffsetRecurse(IntSequence& v, int nv)
{
  if (v.size() == 0)
    return 0;
  int prefix = v.getPrefixLength();
  int m = v[0];
  int k = v.size();
  int s1 = PascalTriangle::noverk(nv + k - 1, k) - PascalTriangle::noverk(nv - m + k - 1, k);
  IntSequence subv(v, prefix, k);
  subv.add(-m);
  int s2 = getOffsetRecurse(subv, nv - m);
  return s1 + s2;
}
