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

#include "fine_container.hh"

#include <cmath>

/* Here we construct the vector of new sizes of containers (before
   nc) and copy all remaining sizes behind nc. */

SizeRefinement::SizeRefinement(const IntSequence &s, int nc, int max)
{
  new_nc = 0;
  for (int i = 0; i < nc; i++)
    {
      int nr = s[i]/max;
      if (s[i] % max != 0)
        nr++;
      int ss = (nr > 0) ? static_cast<int>(round(static_cast<double>(s[i])/nr)) : 0;
      for (int j = 0; j < nr - 1; j++)
        {
          rsizes.push_back(ss);
          ind_map.push_back(i);
          new_nc++;
        }
      rsizes.push_back(s[i]-(nr-1)*ss);
      ind_map.push_back(i);
      new_nc++;
    }

  for (int i = nc; i < s.size(); i++)
    {
      rsizes.push_back(s[i]);
      ind_map.push_back(i);
    }
}
