function [fval,grad,hess,exit_flag,info,PHI,SIGMAu,iXX,prior] = dsge_var_likelihood(xparam1,DynareDataset,DynareOptions,Model,EstimatedParameters,BayesInfo,DynareResults)
% Evaluates the posterior kernel of the bvar-dsge model. Deprecated interface.
%
% INPUTS
%   o xparam1       [double]     Vector of model's parameters.
%   o gend          [integer]    Number of observations (without conditionning observations for the lags).
%
% OUTPUTS
%   o fval          [double]     Value of the posterior kernel at xparam1.
%   o cost_flag     [integer]    Zero if the function returns a penalty, one otherwise.
%   o info          [integer]    Vector of informations about the penalty.
%   o PHI           [double]     Stacked BVAR-DSGE autoregressive matrices (at the mode associated to xparam1).
%   o SIGMAu        [double]     Covariance matrix of the BVAR-DSGE (at the mode associated to xparam1).
%   o iXX           [double]     inv(X'X).
%   o prior         [double]     a matlab structure describing the dsge-var prior.
%
% SPECIAL REQUIREMENTS
%   None.

% Copyright (C) 2006-2012 Dynare Team
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
% along with Dynare.  If not, see <http://www.gnu.org/licenses/>.

[fval,info,exit_flag,grad,hess,SteadyState,trend_coeff,PHI,SIGMAu,iXX,prior] = ...
    dsge_var_likelihood_1(xparam1,DynareDataset,DynareInfo,DynareOptions,Model,...
                          EstimatedParameters,BayesInfo,BoundsInfo,DynareResults);