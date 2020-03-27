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
 * along with Dynare.  If not, see <http://www.gnu.org/licenses/>.
 */

/*
  along shocks: m    mult    max_evals
  ellipse:      m    mult    max_evals  (10·m) (0.5·mult)
  simul:        m            max_evals  (10·m)

  --check-scale 2.0 --check-evals 1000 --check-num 10 --check PES
*/

#include <vector>
#include <string>

struct DynareParams
{
  std::string modname;
  std::string basename;
  int num_per;
  int num_burn;
  int num_sim;
  int num_rtper;
  int num_rtsim;
  int num_condper;
  int num_condsim;
  int num_threads;
  int num_steps;
  std::string prefix;
  int seed;
  int order;
  /* Tolerance used for steady state calcs. */
  double ss_tol;
  bool check_along_path;
  bool check_along_shocks;
  bool check_on_ellipse;
  int check_evals;
  int check_num;
  double check_scale;
  /* Flag for doing IRFs even if the irf_list is empty. */
  bool do_irfs_all;
  /* List of shocks for which IRF will be calculated. */
  std::vector<std::string> irf_list;
  bool do_centralize;
  double qz_criterium;
  bool help;
  bool version;
  DynareParams(int argc, char **argv);
  void printHelp() const;
  int
  getCheckShockPoints() const
  {
    return check_num;
  }
  double
  getCheckShockScale() const
  {
    return check_scale;
  }
  int
  getCheckEllipsePoints() const
  {
    return 10*check_num;
  }
  double
  getCheckEllipseScale() const
  {
    return 0.5*check_scale;
  }
  int
  getCheckPathPoints() const
  {
    return 10*check_num;
  }
private:
  enum class opt { per, burn, sim, rtper, rtsim, condper, condsim,
                   prefix, threads,
                   steps, seed, order, ss_tol, check,
                   check_evals, check_scale, check_num, noirfs, irfs,
                   help, version, centralize, no_centralize, qz_criterium };
  void processCheckFlags(const std::string &flags);
  /* This gathers strings from argv[optind] and on not starting with '-' to the
     irf_list. It stops one item before the end, since this is the model
     file. */
  void processIRFList(int argc, char **argv);
};
