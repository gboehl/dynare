function [nvar,vartan,CorrFileNumber] = dsge_simulated_theoretical_correlation(SampleSize,nar,M_,options_,oo_,type)
% function [nvar,vartan,CorrFileNumber] = dsge_simulated_theoretical_correlation(SampleSize,nar,M_,options_,oo_,type)
% This function computes the posterior or prior distribution of the endogenous
% variables' second order moments. Actual computations are done in
% dsge_simulated_theoretical_covariance, see https://git.dynare.org/Dynare/dynare/-/issues/1769
%
% INPUTS
%   SampleSize   [integer]          scalar, number of simulations.
%   nar          [integer]          maximum number of autocorrelations to
%                                   consider
%   M_           [structure]        Dynare structure describing the model.
%   options_     [structure]        Dynare structure defining global options
%   oo_          [structure]        Dynare structure where the results are saved.
%   type         [string]           'prior' or 'posterior'
%
% OUTPUTS
%   nvar           [integer]        nvar is the number of stationary variables.
%   vartan         [char]           array of characters (with nvar rows).
%   CorrFileNumber [integer]        scalar, number of prior or posterior data files (for correlation).

% Copyright (C) 2007-2021 Dynare Team
%
% This file is part of Dynare.
%
% Dynare is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
%
% Dynare is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with Dynare.  If not, see <https://www.gnu.org/licenses/>.

nvar = length(ivar);
[ivar,vartan, options_] = get_variables_list(options_, M_);

% Get informations about the _posterior_draws files.
if strcmpi(type,'posterior')
    CorrFileNumber = length(dir([M_.dname '/metropolis/' M_.fname '_PosteriorCorrelations*']));
elseif strcmpi(type,'prior')
    CorrFileNumber = length(dir([M_.dname '/prior/moments/' M_.fname '_PriorCorrelations*']));
else
    disp('dsge_simulated_theoretical_correlation:: Unknown type!');
    error()
end
