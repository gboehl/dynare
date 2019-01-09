// Copyright (C) 2008-2011, Ondra Kamenik

#include "parser/cc/matrix_parser.hh"
#include "utils/cc/memory_file.hh"
#include "utils/cc/exception.hh"
#include "sylv/cc/GeneralMatrix.hh"
#include "sylv/cc/Vector.hh"
#include "sylv/cc/SymSchurDecomp.hh"
#include "sylv/cc/SylvException.hh"
#include "integ/cc/quadrature.hh"
#include "integ/cc/smolyak.hh"
#include "integ/cc/product.hh"

#include <getopt.h>
#include <cstdio>

#include <cmath>

struct QuadParams
{
  const char *outname;
  const char *vcovname;
  int max_level;
  double discard_weight;
  QuadParams(int argc, char **argv);
  void check_consistency() const;
private:
  enum {opt_max_level, opt_discard_weight, opt_vcov};
};

QuadParams::QuadParams(int argc, char **argv)
  : outname(NULL), vcovname(NULL), max_level(3), discard_weight(0.0)
{
  if (argc == 1)
    {
      // print the help and exit
      exit(1);
    }

  outname = argv[argc-1];
  argc--;

  struct option const opts [] = {
    {"max-level", required_argument, NULL, opt_max_level},
    {"discard-weight", required_argument, NULL, opt_discard_weight},
    {"vcov", required_argument, NULL, opt_vcov},
    {NULL, 0, NULL, 0}
  };

  int ret;
  int index;
  while (-1 != (ret = getopt_long(argc, argv, "", opts, &index)))
    {
      switch (ret)
        {
        case opt_max_level:
          if (1 != sscanf(optarg, "%d", &max_level))
            fprintf(stderr, "Couldn't parse integer %s, ignored\n", optarg);
          break;
        case opt_discard_weight:
          if (1 != sscanf(optarg, "%lf", &discard_weight))
            fprintf(stderr, "Couldn't parse float %s, ignored\n", optarg);
          break;
        case opt_vcov:
          vcovname = optarg;
          break;
        }
    }

  check_consistency();
}

void
QuadParams::check_consistency() const
{
  if (outname == NULL)
    {
      fprintf(stderr, "Error: output name not set\n");
      exit(1);
    }

  if (vcovname == NULL)
    {
      fprintf(stderr, "Error: vcov file name not set\n");
      exit(1);
    }
}

/** Utility class for ordering pointers to vectors according their
 * ordering. */
struct OrderVec
{
  bool
  operator()(const Vector *a, const Vector *b) const
  {
    return *a < *b;
  }
};

int
main(int argc, char **argv)
{
  QuadParams params(argc, argv);

  // open output file for writing
  FILE *fout;
  if (NULL == (fout = fopen(params.outname, "w")))
    {
      fprintf(stderr, "Could not open %s for writing\n", params.outname);
      exit(1);
    }

  try
    {

      // open memory file for vcov
      ogu::MemoryFile vcov_mf(params.vcovname);

      // parse the vcov matrix
      ogp::MatrixParser mp;
      mp.parse(vcov_mf.length(), vcov_mf.base());
      if (mp.nrows() != mp.ncols())
        throw ogu::Exception(__FILE__, __LINE__,
                             "VCOV matrix not square");
      // and put to the GeneralMatrix
      GeneralMatrix vcov(mp.nrows(), mp.ncols());
      vcov.zeros();
      for (ogp::MPIterator it = mp.begin(); it != mp.end(); ++it)
        vcov.get(it.row(), it.col()) = *it;

      // calculate the factor A of vcov, so that A*A^T=VCOV
      GeneralMatrix A(vcov.numRows(), vcov.numRows());
      SymSchurDecomp ssd(vcov);
      ssd.getFactor(A);

      // construct Gauss-Hermite quadrature
      GaussHermite ghq;
      // construct Smolyak quadrature
      int level = params.max_level;
      SmolyakQuadrature sq(vcov.numRows(), level, ghq);

      printf("Dimension:                %d\n", vcov.numRows());
      printf("Maximum level:            %d\n", level);
      printf("Total number of nodes:    %d\n", sq.numEvals(level));

      // put the points to the vector
      std::vector<Vector *> points;
      for (smolpit qit = sq.start(level); qit != sq.end(level); ++qit)
        points.push_back(new Vector((const Vector &) qit.point()));
      // sort and uniq
      OrderVec ordvec;
      std::sort(points.begin(), points.end(), ordvec);
      auto new_end = std::unique(points.begin(), points.end());
      for (auto it = new_end; it != points.end(); ++it)
        delete *it;
      points.erase(new_end, points.end());

      printf("Duplicit nodes removed:   %lu\n", (unsigned long) (sq.numEvals(level)-points.size()));

      // calculate weights and mass
      double mass = 0.0;
      std::vector<double> weights;
      for (auto & point : points)
        {
          weights.push_back(std::exp(-point->dot(*point)));
          mass += weights.back();
        }

      // calculate discarded mass
      double discard_mass = 0.0;
      for (double weight : weights)
        if (weight/mass < params.discard_weight)
          discard_mass += weight;

      printf("Total mass discarded:     %f\n", discard_mass/mass);

      // dump the results
      int npoints = 0;
      double upscale_weight = 1/(mass-discard_mass);
      Vector x(vcov.numRows());
      for (int i = 0; i < (int) weights.size(); i++)
        if (weights[i]/mass >= params.discard_weight)
          {
            // print the upscaled weight
            fprintf(fout, "%20.16g", upscale_weight*weights[i]);
            // multiply point with the factor A and sqrt(2)
            A.multVec(0.0, x, std::sqrt(2.), *(points[i]));
            // print the coordinates
            for (int j = 0; j < x.length(); j++)
              fprintf(fout, " %20.16g", x[j]);
            fprintf(fout, "\n");
            npoints++;
          }

      printf("Final number of points:   %d\n", npoints);

      fclose(fout);

    }
  catch (const SylvException &e)
    {
      e.printMessage();
      return 1;
    }
  catch (const ogu::Exception &e)
    {
      e.print();
      return 1;
    }

  return 0;
}
