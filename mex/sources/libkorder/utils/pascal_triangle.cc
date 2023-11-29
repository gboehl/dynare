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

#include "pascal_triangle.hh"

#include <iostream>
#include <mutex>

void
PascalRow::setFromPrevious(const PascalRow& prev)
{
  k = prev.k + 1;
  clear();
  prolong(prev);
}

/* This prolongs the PascalRow. If it is empty, we set the first item
                    ⎛k+1⎞
   to k+1, which is ⎝ k ⎠, which is the second item in the real
                                 ⎛k⎞
   pascal row, which starts from ⎝k⎠=1. Then we calculate
   other items from the provided row which must be the one with k-1.
 */
void
PascalRow::prolong(const PascalRow& prev)
{
  if (size() == 0)
    push_back(k + 1);
  int last = back();
  for (unsigned int i = size(); i < prev.size(); i++)
    {
      last += prev[i];
      push_back(last);
    }
}

void
PascalRow::prolongFirst(int n)
{
  // TODO: check n = 1
  for (int i = static_cast<int>(size()) + 2; i <= n; i++)
    push_back(i);
}

void
PascalRow::print() const
{
  std::cout << "k=" << k << std::endl;
  for (unsigned int i = 0; i < size(); i++)
    std::cout << operator[](i) << ' ';
  std::cout << std::endl;
}

namespace PascalTriangle
{
namespace // Anonymous namespace that is a functional equivalent of “private”
{
std::vector<PascalRow> tr(1);
std::mutex mut; // For protecting the triangle from concurrent modifications

int
max_n()
{
  return static_cast<int>(tr[0].size() + 1);
}

int
max_k()
{
  return static_cast<int>(tr.size());
}
}

void
ensure(int n, int k)
{
  // Add along n
  if (n > max_n())
    {
      std::lock_guard<std::mutex> lk {mut};
      tr[0].prolongFirst(n);
      for (int i = 2; i <= max_k(); i++)
        tr[i - 1].prolong(tr[i - 2]);
    }

  if (k > max_k())
    {
      std::lock_guard<std::mutex> lk {mut};
      for (int i = max_k() + 1; i <= k; i++)
        {
          tr.emplace_back();
          tr.back().setFromPrevious(tr[i - 2]);
        }
    }
}

int
noverk(int n, int k)
{
  // TODO: raise exception if out of bounds
  if (n - k < k)
    k = n - k;
  if (k == 0)
    return 1;
  ensure(n, k);
  return (tr[k - 1])[n - 1 - k];
}

void
print()
{
  for (const auto& i : tr)
    i.print();
}
}
