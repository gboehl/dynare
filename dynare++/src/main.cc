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
 * along with Dynare.  If not, see <http://www.gnu.org/licenses/>.
 */

#include "dynare3.hh"
#include "dynare_exception.hh"
#include "dynare_params.hh"

#include "utils/cc/exception.hh"
#include "parser/cc/parser_exception.hh"
#include "../sylv/cc/SylvException.hh"
#include "../kord/seed_generator.hh"
#include "../kord/global_check.hh"
#include "../kord/approximation.hh"

#include <fstream>
#include <iostream>
#include <cstdlib>

int
main(int argc, char **argv)
{
  DynareParams params(argc, argv);
  if (params.help)
    {
      params.printHelp();
      return EXIT_SUCCESS;
    }
  if (params.version)
    {
      std::cout << "Dynare++ v. " << VERSION << '\n'
                << '\n'
                << u8"Copyright © 2004-2011 Ondra Kamenik\n"
                << u8"Copyright © 2019 Dynare Team\n"
                << "Dynare++ comes with ABSOLUTELY NO WARRANTY and is distributed under the GNU GPL,\n"
                << "version 3 or later (see https://www.gnu.org/licenses/gpl.html)\n";
      return EXIT_SUCCESS;
    }
  sthread::detach_thread_group::max_parallel_threads = params.num_threads;

  try
    {
      // make journal
      Journal journal(params.basename + ".jnl");

      // make dynare object
      Dynare dynare(params.modname, params.order, params.ss_tol, journal);
      // make list of shocks for which we will do IRFs
      std::vector<int> irf_list_ind;
      if (params.do_irfs_all)
        for (int i = 0; i < dynare.nexog(); i++)
          irf_list_ind.push_back(i);
      else
        irf_list_ind = static_cast<const DynareNameList &>(dynare.getExogNames()).selectIndices(params.irf_list);

      // write matlab files
      std::string mfile1(params.basename + "_f.m");
      std::ofstream mfd{mfile1, std::ios::out | std::ios::trunc};
      if (mfd.fail())
        {
          std::cerr << "Couldn't open " << mfile1 << " for writing.\n";
          std::exit(EXIT_FAILURE);
        }
      ogdyn::MatlabSSWriter writer0(dynare.getModel(), params.basename);
      writer0.write_der0(mfd);
      mfd.close();

      std::string mfile2(params.basename + "_ff.m");
      mfd.open(mfile2, std::ios::out | std::ios::trunc);
      if (mfd.fail())
        {
          std::cerr << "Couldn't open " << mfile2 << " for writing.\n";
          std::exit(EXIT_FAILURE);
        }
      ogdyn::MatlabSSWriter writer1(dynare.getModel(), params.basename);
      writer1.write_der1(mfd);
      mfd.close();

      // open mat file
      std::string matfile(params.basename + ".mat");
      mat_t *matfd = Mat_Create(matfile.c_str(), nullptr);
      if (!matfd)
        {
          std::cerr << "Couldn't open " << matfile << " for writing.\n";
          std::exit(EXIT_FAILURE);
        }

      // write info about the model (dimensions and variables)
      dynare.writeMat(matfd, params.prefix);
      // write the dump file corresponding to the input
      dynare.writeDump(params.basename);

      seed_generator::set_meta_seed(static_cast<std::mt19937::result_type>(params.seed));

      TLStatic::init(dynare.order(),
                     dynare.nstat()+2*dynare.npred()+3*dynare.nboth()
                     +2*dynare.nforw()+dynare.nexog());

      Approximation app(dynare, journal, params.num_steps, params.do_centralize, params.qz_criterium);
      try
        {
          app.walkStochSteady();
        }
      catch (const KordException &e)
        {
          // tell about the exception and continue
          std::cout << "Caught (not yet fatal) Kord exception: ";
          e.print();
          JournalRecord rec(journal);
          rec << "Solution routine not finished (" << e.get_message()
              << "), see what happens" << endrec;
        }

      std::string ss_matrix_name(params.prefix + "_steady_states");
      ConstTwoDMatrix(app.getSS()).writeMat(matfd, ss_matrix_name);

      // check the approximation
      if (params.check_along_path || params.check_along_shocks
          || params.check_on_ellipse)
        {
          GlobalChecker gcheck(app, sthread::detach_thread_group::max_parallel_threads, journal);
          if (params.check_along_shocks)
            gcheck.checkAlongShocksAndSave(matfd, params.prefix,
                                           params.getCheckShockPoints(),
                                           params.getCheckShockScale(),
                                           params.check_evals);
          if (params.check_on_ellipse)
            gcheck.checkOnEllipseAndSave(matfd, params.prefix,
                                         params.getCheckEllipsePoints(),
                                         params.getCheckEllipseScale(),
                                         params.check_evals);
          if (params.check_along_path)
            gcheck.checkAlongSimulationAndSave(matfd, params.prefix,
                                               params.getCheckPathPoints(),
                                               params.check_evals);
        }

      // write the folded decision rule to the Mat-4 file
      app.getFoldDecisionRule().writeMat(matfd, params.prefix);

      // simulate conditional
      if (params.num_condper > 0 && params.num_condsim > 0)
        {
          SimResultsDynamicStats rescond(dynare.numeq(), params.num_condper, 0);
          Vector det_ss{app.getSS().getCol(0)};
          rescond.simulate(params.num_condsim, app.getFoldDecisionRule(), det_ss, dynare.getVcov(), journal);
          rescond.writeMat(matfd, params.prefix);
        }

      // simulate unconditional
      //const DecisionRule& dr = app.getUnfoldDecisionRule();
      const DecisionRule &dr = app.getFoldDecisionRule();
      if (params.num_per > 0 && params.num_sim > 0)
        {
          SimResultsStats res(dynare.numeq(), params.num_per, params.num_burn);
          res.simulate(params.num_sim, dr, dynare.getSteady(), dynare.getVcov(), journal);
          res.writeMat(matfd, params.prefix);

          // impulse response functions
          if (!irf_list_ind.empty())
            {
              IRFResults irf(dynare, dr, res, irf_list_ind, journal);
              irf.writeMat(matfd, params.prefix);
            }
        }

      // simulate with real-time statistics
      if (params.num_rtper > 0 && params.num_rtsim > 0)
        {
          RTSimResultsStats rtres(dynare.numeq(), params.num_rtper, params.num_burn);
          rtres.simulate(params.num_rtsim, dr, dynare.getSteady(), dynare.getVcov(), journal);
          rtres.writeMat(matfd, params.prefix);
        }

      Mat_Close(matfd);
    }
  catch (const KordException &e)
    {
      std::cout << "Caught Kord exception: ";
      e.print();
      return e.code();
    }
  catch (const TLException &e)
    {
      std::cout << "Caught TL exception: ";
      e.print();
      return 255;
    }
  catch (SylvException &e)
    {
      std::cout << "Caught Sylv exception: ";
      e.printMessage();
      return 255;
    }
  catch (const DynareException &e)
    {
      std::cout << "Caught Dynare exception: " << e.message() << '\n';
      return 255;
    }
  catch (const ogu::Exception &e)
    {
      std::cout << "Caught ogu::Exception: ";
      e.print();
      return 255;
    }
  catch (const ogp::ParserException &e)
    {
      std::cout << "Caught parser exception: " << e.message() << '\n';
      return 255;
    }

  return EXIT_SUCCESS;
}
