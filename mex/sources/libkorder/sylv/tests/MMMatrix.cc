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

#include "MMMatrix.hh"

#include <fstream>
#include <iomanip>

MMMatrixIn::MMMatrixIn(const std::string& fname)
{
  std::ifstream fd {fname};
  if (fd.fail())
    throw MMException("Cannot open file " + fname + " for reading\n");

  // jump over initial comments
  while (fd.peek() == '%')
    fd.ignore(std::numeric_limits<std::streamsize>::max(), '\n');

  // read in number of rows and cols
  fd >> rows >> cols;
  if (fd.fail())
    throw MMException("Couldn't parse rows and cols\n");
  // read in data
  int len = rows * cols;
  data.resize(len);
  int i = 0;
  while (!fd.eof() && i < len)
    {
      fd >> data[i];
      if (fd.fail())
        throw MMException("Couldn't parse float number\n");
      i++;
    }
  if (i < len)
    throw MMException("Couldn't read all " + std::to_string(len) + " elements, read "
                      + std::to_string(i) + " so far\n");
  fd.close();
}

void
MMMatrixOut::write(const std::string& fname, const GeneralMatrix& m)
{
  std::ofstream fd {fname, std::ios::out | std::ios::trunc};
  if (fd.fail())
    throw MMException("Cannot open file " + fname + " for writing\n");

  fd << "%%%%MatrixMarket matrix array real general" << std::endl
     << m.nrows() << ' ' << m.ncols() << std::endl
     << std::setprecision(35);
  for (int i = 0; i < m.ncols(); i++)
    for (int j = 0; j < m.nrows(); j++)
      fd << std::setw(40) << m.get(i, j);
  fd.close();
}
