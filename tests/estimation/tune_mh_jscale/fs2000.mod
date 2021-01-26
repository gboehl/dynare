/*
 * Copyright © 2004-2020 Dynare Team
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

@#include "fs2000.inc"

estimation(order=1, datafile='../fsdat_simul', nobs=192, loglinear, mh_replic=10000, mh_nblocks=1, mh_tune_jscale=0.33,mh_tune_guess=0.7,plot_priors=0);

mhdata = load('fs2000/metropolis/fs2000_mh_history_0.mat');

if any(abs(mhdata.record.AcceptanceRatio-options_.mh_tune_jscale.target)>options_.mh_tune_jscale.c2)
    error('Automagic tuning of the MCMC proposal scale parameter did not work as expected!')
end
