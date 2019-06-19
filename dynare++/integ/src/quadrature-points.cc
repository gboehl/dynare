/*
 * Copyright © 2008-2011 Ondra Kamenik
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

#include "parser/cc/matrix_parser.hh"
#include "utils/cc/exception.hh"
#include "sylv/cc/GeneralMatrix.hh"
#include "sylv/cc/Vector.hh"
#include "sylv/cc/SymSchurDecomp.hh"
#include "sylv/cc/SylvException.hh"
#include "integ/cc/quadrature.hh"
#include "integ/cc/smolyak.hh"
#include "integ/cc/product.hh"

#include <getopt.h>

#include <cmath>
#include <cstdlib>
#include <iostream>
#include <fstream>
#include <sstream>
#include <memory>
#include <string>

struct QuadParams
{
  std::string outname;
  std::string vcovname;
  int max_level{3};
  double discard_weight{0.0};
  QuadParams(int argc, char **argv);
  void check_consistency() const;
private:
  enum class opt {max_level, discard_weight, vcov};
};

QuadParams::QuadParams(int argc, char **argv)
{
  if (argc == 1)
    {
      // Print the help and exit
      std::cerr << "Usage: " << argv[0] << " [--max-level INTEGER] [--discard-weight FLOAT] [--vcov FILENAME] OUTPUT_FILENAME" << std::endl;
      std::exit(EXIT_FAILURE);
    }

  outname = argv[argc-1];
  argc--;

  struct option const opts [] = {
    {"max-level", required_argument, nullptr, static_cast<int>(opt::max_level)},
    {"discard-weight", required_argument, nullptr, static_cast<int>(opt::discard_weight)},
    {"vcov", required_argument, nullptr, static_cast<int>(opt::vcov)},
    {nullptr, 0, nullptr, 0}
  };

  int ret;
  int index;
  while (-1 != (ret = getopt_long(argc, argv, "", opts, &index)))
    {
      if (ret == '?')
        {
          std::cerr << "Unknown option, ignored\n";
          continue;
        }

      switch (static_cast<opt>(ret))
        {
        case opt::max_level:
          try
            {
              max_level = std::stoi(optarg);
            }
          catch (const std::invalid_argument &e)
            {
              std::cerr << "Couldn't parse integer " << optarg << ", ignored" << std::endl;
            }
          break;
        case opt::discard_weight:
          try
            {
              discard_weight = std::stod(optarg);
            }
          catch (const std::invalid_argument &e)
            {
              std::cerr << "Couldn't parse float " << optarg << ", ignored" << std::endl;
            }
          break;
        case opt::vcov:
          vcovname = optarg;
          break;
        }
    }

  check_consistency();
}

void
QuadParams::check_consistency() const
{
  if (outname.empty())
    {
      std::cerr << "Error: output name not set" << std::endl;
      std::exit(EXIT_FAILURE);
    }

  if (vcovname.empty())
    {
      std::cerr << "Error: vcov file name not set" << std::endl;
      std::exit(EXIT_FAILURE);
    }
}

int
main(int argc, char **argv)
{
  QuadParams params(argc, argv);

  // Open output file for writing
  std::ofstream fout{params.outname, std::ios::out | std::ios::trunc};
  if (fout.fail())
    {
      std::cerr << "Could not open " << params.outname << " for writing" << std::endl;
      std::exit(EXIT_FAILURE);
    }

  try
    {
      std::ifstream f{params.vcovname};
      std::ostringstream buffer;
      buffer << f.rdbuf();
      std::string contents{buffer.str()};

      // Parse the vcov matrix
      ogp::MatrixParser mp;
      mp.parse(contents);
      if (mp.nrows() != mp.ncols())
        throw ogu::Exception(__FILE__, __LINE__,
                             "VCOV matrix not square");
      // And put to the GeneralMatrix
      GeneralMatrix vcov(mp.nrows(), mp.ncols());
      vcov.zeros();
      for (ogp::MPIterator it = mp.begin(); it != mp.end(); ++it)
        vcov.get(it.row(), it.col()) = *it;

      // Calculate the factor A of vcov, so that A·Aᵀ=VCOV
      GeneralMatrix A(vcov.nrows(), vcov.nrows());
      SymSchurDecomp ssd(vcov);
      ssd.getFactor(A);

      // Construct Gauss-Hermite quadrature
      GaussHermite ghq;
      // Construct Smolyak quadrature
      int level = params.max_level;
      SmolyakQuadrature sq(vcov.nrows(), level, ghq);

      std::cout << "Dimension:                " << vcov.nrows() << std::endl
                << "Maximum level:            " << level << std::endl
                << "Total number of nodes:    " << sq.numEvals(level) << std::endl;

      // Put the points to the vector
      std::vector<std::unique_ptr<Vector>> points;
      for (smolpit qit = sq.start(level); qit != sq.end(level); ++qit)
        points.push_back(std::make_unique<Vector>(const_cast<const Vector &>(qit.point())));
      // Sort and uniq
      std::sort(points.begin(), points.end(), [](auto &a, auto &b) { return a.get() < b.get(); });
      auto new_end = std::unique(points.begin(), points.end());
      points.erase(new_end, points.end());

      std::cout << "Duplicit nodes removed:   " << static_cast<unsigned long>(sq.numEvals(level)-points.size())
                << std::endl;

      // Calculate weights and mass
      double mass = 0.0;
      std::vector<double> weights;
      for (auto & point : points)
        {
          weights.push_back(std::exp(-point->dot(*point)));
          mass += weights.back();
        }

      // Calculate discarded mass
      double discard_mass = 0.0;
      for (double weight : weights)
        if (weight/mass < params.discard_weight)
          discard_mass += weight;

      std::cout << "Total mass discarded:     " << std::fixed << discard_mass/mass << std::endl;

      // Dump the results
      int npoints = 0;
      double upscale_weight = 1/(mass-discard_mass);
      Vector x(vcov.nrows());
      fout << std::setprecision(16);
      for (int i = 0; i < static_cast<int>(weights.size()); i++)
        if (weights[i]/mass >= params.discard_weight)
          {
            // Print the upscaled weight
            fout << std::setw(20) << upscale_weight*weights[i];
            // Multiply point with the factor A and √2
            A.multVec(0.0, x, std::sqrt(2.), *(points[i]));
            // Print the coordinates
            for (int j = 0; j < x.length(); j++)
              fout << ' ' << std::setw(20) << x[j];
            fout << std::endl;
            npoints++;
          }

      std::cout << "Final number of points:   " << npoints << std::endl;

      fout.close();
    }
  catch (const SylvException &e)
    {
      e.printMessage();
      return EXIT_FAILURE;
    }
  catch (const ogu::Exception &e)
    {
      e.print();
      return EXIT_FAILURE;
    }

  return EXIT_SUCCESS;
}
