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

#include "quasi_mcarlo.hh"

#include <cmath>
#include <iostream>
#include <iomanip>
#include <array>

/* Here in the constructor, we have to calculate a maximum length of ‘coeff’
   array for a given ‘base’ and given maximum ‘maxn’. After allocation, we
   calculate the coefficients. */

RadicalInverse::RadicalInverse(int n, int b, int mxn)
  : num(n), base(b), maxn(mxn),
    coeff(static_cast<int>(floor(log(static_cast<double>(maxn))/log(static_cast<double>(b)))+2), 0)
{
  int nr = num;
  j = -1;
  do
    {
      j++;
      coeff[j] = nr % base;
      nr = nr / base;
    }
  while (nr > 0);
}

/* This evaluates the radical inverse. If there was no permutation, we have to
   calculate:

    c₀   c₁        cⱼ
    ── + ── + … + ────
     b   b²       bʲ⁺¹

   which is evaluated as:

    ⎛ ⎛⎛cⱼ 1   cⱼ₋₁⎞ 1   cⱼ₋₂⎞ ⎞ 1   c₀
    ⎢…⎢⎢──·─ + ────⎥·─ + ────⎥…⎥·─ + ──
    ⎝ ⎝⎝ b b     b ⎠ b     b ⎠ ⎠ b    b


   Now with permutation π, we have:

    ⎛ ⎛⎛π(cⱼ) 1   π(cⱼ₋₁)⎞ 1   π(cⱼ₋₂)⎞ ⎞ 1   π(c₀)
    ⎢…⎢⎢─────·─ + ───────⎥·─ + ───────⎥…⎥·─ + ─────
    ⎝ ⎝⎝  b   b       b  ⎠ b       b  ⎠ ⎠ b     b
*/

double
RadicalInverse::eval(const PermutationScheme &p) const
{
  double res = 0;
  for (int i = j; i >= 0; i--)
    {
      int cper = p.permute(i, base, coeff[i]);
      res = (cper + res)/base;
    }
  return res;
}

/* We just add 1 to the lowest coefficient and check for overflow with respect
   to the base. */
void
RadicalInverse::increase()
{
  // TODO: raise if num+1 > maxn
  num++;
  int i = 0;
  coeff[i]++;
  while (coeff[i] == base)
    {
      coeff[i] = 0;
      coeff[++i]++;
    }
  if (i > j)
    j = i;
}

/* Debug print. */
void
RadicalInverse::print() const
{
  std::cout << "n=" << num << " b=" << base << " c=";
  coeff.print();
}

/* Here we have the first 170 primes. This means that we are not able to
   integrate dimensions greater than 170. */

std::array<int, 170> HaltonSequence::primes =
  {
   2, 3, 5, 7, 11, 13, 17, 19, 23, 29,
   31, 37, 41, 43, 47, 53, 59, 61, 67, 71,
   73, 79, 83, 89, 97, 101, 103, 107, 109, 113,
   127, 131, 137, 139, 149, 151, 157, 163, 167, 173,
   179, 181, 191, 193, 197, 199, 211, 223, 227, 229,
   233, 239, 241, 251, 257, 263, 269, 271, 277, 281,
   283, 293, 307, 311, 313, 317, 331, 337, 347, 349,
   353, 359, 367, 373, 379, 383, 389, 397, 401, 409,
   419, 421, 431, 433, 439, 443, 449, 457, 461, 463,
   467, 479, 487, 491, 499, 503, 509, 521, 523, 541,
   547, 557, 563, 569, 571, 577, 587, 593, 599, 601,
   607, 613, 617, 619, 631, 641, 643, 647, 653, 659,
   661, 673, 677, 683, 691, 701, 709, 719, 727, 733,
   739, 743, 751, 757, 761, 769, 773, 787, 797, 809,
   811, 821, 823, 827, 829, 839, 853, 857, 859, 863,
   877, 881, 883, 887, 907, 911, 919, 929, 937, 941,
   947, 953, 967, 971, 977, 983, 991, 997, 1009, 1013
  };

/* This takes first ‘dim’ primes and constructs ‘dim’ radical inverses and
   calls eval(). */

HaltonSequence::HaltonSequence(int n, int mxn, int dim, const PermutationScheme &p)
  : num(n), maxn(mxn), per(p), pt(dim)
{
  // TODO: raise if dim > num_primes
  // TODO: raise if n > mxn
  for (int i = 0; i < dim; i++)
    ri.emplace_back(num, primes[i], maxn);
  eval();
}

/* This calls RadicalInverse::increase() for all radical inverses and calls
   eval(). */

void
HaltonSequence::increase()
{
  for (auto &i : ri)
    i.increase();
  num++;
  if (num <= maxn)
    eval();
}

/* This sets point ‘pt’ to radical inverse evaluations in each dimension. */
void
HaltonSequence::eval()
{
  for (unsigned int i = 0; i < ri.size(); i++)
    pt[i] = ri[i].eval(per);
}

/* Debug print. */
void
HaltonSequence::print() const
{
  auto ff = std::cout.flags();
  for (const auto &i : ri)
    i.print();
  std::cout << "point=[ "
            << std::fixed << std::setprecision(6);
  for (unsigned int i = 0; i < ri.size(); i++)
    std::cout << std::setw(7) << pt[i] << ' ';
  std::cout << ']' << std::endl;
  std::cout.flags(ff);
}

qmcpit::qmcpit(const QMCSpecification &s, int n)
  : spec(s), halton{n, s.level(), s.dimen(), s.getPerScheme()},
    sig{s.dimen()}
{
}

bool
qmcpit::operator==(const qmcpit &qpit) const
{
  return &spec == &qpit.spec && halton.getNum() == qpit.halton.getNum();
}

qmcpit &
qmcpit::operator++()
{
  halton.increase();
  return *this;
}

double
qmcpit::weight() const
{
  return 1.0/spec.level();
}

int
WarnockPerScheme::permute(int i, int base, int c) const
{
  return (c+i) % base;
}

int
ReversePerScheme::permute(int i, int base, int c) const
{
  return (base-c) % base;
}
