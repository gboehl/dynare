/*
 * Copyright © 2004-2011 Ondra Kamenik
 * Copyright © 2019-2022 Dynare Team
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

#include "dynare_params.hh"

#include "sthread.hh"

#include <getopt.h>
#include <string>
#include <iostream>

DynareParams::DynareParams(int argc, char **argv)
  : num_per(100), num_burn(0), num_sim(80),
    num_rtper(0), num_rtsim(0),
    num_condper(0), num_condsim(0),
    num_threads(sthread::default_threads_number()), num_steps(0),
    prefix("dyn"), seed(934098), order(-1), ss_tol(1.e-13),
    check_along_path(false), check_along_shocks(false),
    check_on_ellipse(false), check_evals(1000), check_num(10), check_scale(2.0),
    do_irfs_all(true), do_centralize(true), qz_criterium(1.0+1e-6),
    help(false), version(false)
{
  using namespace std::string_literals;
  if (argc == 1 || argv[1] == "--help"s)
    {
      help = true;
      return;
    }
  if (argc == 1 || argv[1] == "--version"s)
    {
      version = true;
      return;
    }

  modname = argv[argc-1];
  argc--;

  struct option const opts[] =
    {
     {"periods", required_argument, nullptr, static_cast<int>(opt::per)},
     {"per", required_argument, nullptr, static_cast<int>(opt::per)},
     {"burn", required_argument, nullptr, static_cast<int>(opt::burn)},
     {"simulations", required_argument, nullptr, static_cast<int>(opt::sim)},
     {"sim", required_argument, nullptr, static_cast<int>(opt::sim)},
     {"rtperiods", required_argument, nullptr, static_cast<int>(opt::rtper)},
     {"rtper", required_argument, nullptr, static_cast<int>(opt::rtper)},
     {"rtsimulations", required_argument, nullptr, static_cast<int>(opt::rtsim)},
     {"rtsim", required_argument, nullptr, static_cast<int>(opt::rtsim)},
     {"condperiods", required_argument, nullptr, static_cast<int>(opt::condper)},
     {"condper", required_argument, nullptr, static_cast<int>(opt::condper)},
     {"condsimulations", required_argument, nullptr, static_cast<int>(opt::condsim)},
     {"condsim", required_argument, nullptr, static_cast<int>(opt::condsim)},
     {"prefix", required_argument, nullptr, static_cast<int>(opt::prefix)},
     {"threads", required_argument, nullptr, static_cast<int>(opt::threads)},
     {"steps", required_argument, nullptr, static_cast<int>(opt::steps)},
     {"seed", required_argument, nullptr, static_cast<int>(opt::seed)},
     {"order", required_argument, nullptr, static_cast<int>(opt::order)},
     {"ss-tol", required_argument, nullptr, static_cast<int>(opt::ss_tol)},
     {"check", required_argument, nullptr, static_cast<int>(opt::check)},
     {"check-scale", required_argument, nullptr, static_cast<int>(opt::check_scale)},
     {"check-evals", required_argument, nullptr, static_cast<int>(opt::check_evals)},
     {"check-num", required_argument, nullptr, static_cast<int>(opt::check_num)},
     {"qz-criterium", required_argument, nullptr, static_cast<int>(opt::qz_criterium)},
     {"no-irfs", no_argument, nullptr, static_cast<int>(opt::noirfs)},
     {"irfs", no_argument, nullptr, static_cast<int>(opt::irfs)},
     {"centralize", no_argument, nullptr, static_cast<int>(opt::centralize)},
     {"no-centralize", no_argument, nullptr, static_cast<int>(opt::no_centralize)},
     {"help", no_argument, nullptr, static_cast<int>(opt::help)},
     {"version", no_argument, nullptr, static_cast<int>(opt::version)},
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

      try
        {
          switch (static_cast<opt>(ret))
            {
            case opt::per:
              num_per = std::stoi(optarg);
              break;
            case opt::burn:
              num_burn = std::stoi(optarg);
              break;
            case opt::sim:
              num_sim = std::stoi(optarg);
              break;
            case opt::rtper:
              num_rtper = std::stoi(optarg);
              break;
            case opt::rtsim:
              num_rtsim = std::stoi(optarg);
              break;
            case opt::condper:
              num_condper = std::stoi(optarg);
              break;
            case opt::condsim:
              num_condsim = std::stoi(optarg);
              break;
            case opt::prefix:
              prefix = optarg;
              break;
            case opt::threads:
              num_threads = std::stoi(optarg);
              break;
            case opt::steps:
              num_steps = std::stoi(optarg);
              break;
            case opt::seed:
              seed = std::stoi(optarg);
              break;
            case opt::order:
              order = std::stoi(optarg);
              break;
            case opt::ss_tol:
              ss_tol = std::stod(optarg);
              break;
            case opt::check:
              processCheckFlags(optarg);
              break;
            case opt::check_scale:
              check_scale = std::stod(optarg);
              break;
            case opt::check_evals:
              check_evals = std::stoi(optarg);
              break;
            case opt::check_num:
              check_num = std::stoi(optarg);
              break;
            case opt::noirfs:
              irf_list.clear();
              do_irfs_all = false;
              break;
            case opt::irfs:
              processIRFList(argc, argv);
              if (irf_list.empty())
                do_irfs_all = true;
              else
                do_irfs_all = false;
              break;
            case opt::centralize:
              do_centralize = true;
              break;
            case opt::no_centralize:
              do_centralize = false;
              break;
            case opt::qz_criterium:
              qz_criterium = std::stod(optarg);
              break;
            case opt::help:
              help = true;
              break;
            case opt::version:
              version = true;
              break;
            }
        }
      catch (std::invalid_argument &)
        {
          std::cerr << "Couldn't parse option " << optarg << ", ignored\n";
        }
      catch (std::out_of_range &)
        {
          std::cerr << "Out-of-range value " << optarg << ", ignored\n";
        }
    }

  // make basename (get rid of the directory and the extension)
  basename = modname;
  auto pos = basename.find_last_of(R"(/\)");
  if (pos != std::string::npos)
    basename = basename.substr(pos+1);
  pos = basename.find_last_of('.');
  if (pos != std::string::npos)
    basename.erase(pos);
}

void
DynareParams::printHelp() const
{
  std::cout << "usage: dynare++ [--help] [--version] [options] <model file>\n"
    "\n"
    "    --help               print this message and return\n"
    "    --version            print version and return\n"
    "\n"
    "options:\n"
    "    --per <num>          number of periods simulated after burnt [100]\n"
    "    --burn <num>         number of periods burnt [0]\n"
    "    --sim <num>          number of simulations [80]\n"
    "    --rtper <num>        number of RT periods simulated after burnt [0]\n"
    "    --rtsim <num>        number of RT simulations [0]\n"
    "    --condper <num>      number of periods in cond. simulations [0]\n"
    "    --condsim <num>      number of conditional simulations [0]\n"
    "    --steps <num>        steps towards stoch. SS [0=deter.]\n"
    "    --centralize         centralize the rule [do centralize]\n"
    "    --no-centralize      do not centralize the rule [do centralize]\n"
    "    --prefix <string>    prefix of variables in Mat-4 file [\"dyn\"]\n"
    "    --seed <num>         random number generator seed [934098]\n"
    "    --order <num>        order of approximation [no default]\n"
    "    --threads <num>      number of max parallel threads [1/2 * nb. of logical CPUs]\n"
    "    --ss-tol <num>       steady state calcs tolerance [1.e-13]\n"
    "    --check pesPES       check model residuals [no checks]\n"
    "                         lower/upper case switches off/on\n"
    "                           pP  checking along simulation path\n"
    "                           eE  checking on ellipse\n"
    "                           sS  checking along shocks\n"
    "    --check-evals <num>  max number of evals per residual [1000]\n"
    "    --check-num <num>    number of checked points [10]\n"
    "    --check-scale <num>  scaling of checked points [2.0]\n"
    "    --no-irfs            shuts down IRF simulations [do IRFs]\n"
    "    --irfs               performs IRF simulations [do IRFs]\n"
    "    --qz-criterium <num> threshold for stable eigenvalues [1.000001]\n"
    "\n\n";
}

void
DynareParams::processCheckFlags(const std::string &flags)
{
  for (char c : flags)
    switch (c)
      {
      case 'p':
        check_along_path = false;
        break;
      case 'P':
        check_along_path = true;
        break;
      case 'e':
        check_on_ellipse = false;
        break;
      case 'E':
        check_on_ellipse = true;
        break;
      case 's':
        check_along_shocks = false;
        break;
      case 'S':
        check_along_shocks = true;
        break;
      default:
        std::cerr << "Unknown check type selection character <" << c << ">, ignored.\n";
      }
}

void
DynareParams::processIRFList(int argc, char **argv)
{
  irf_list.clear();
  while (optind < argc && *(argv[optind]) != '-')
    {
      irf_list.push_back(argv[optind]);
      optind++;
    }
}
