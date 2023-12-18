function pfm = setup_stochastic_perfect_foresight_model_solver(M_,options_,oo_)
% pfm = setup_stochastic_perfect_foresight_model_solver(M_,options_,oo_)
% INPUTS
%  o M_                     [struct]    Dynare's model structure
%  o options_               [struct]    Dynare's options structure
%  o oo_                    [struct]    Dynare's results structure

% Copyright © 2013-2023 Dynare Team
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

pfm.lead_lag_incidence = M_.lead_lag_incidence;
pfm.ny = M_.endo_nbr;
pfm.Sigma = M_.Sigma_e;
if det(pfm.Sigma) > 0
    pfm.Omega = chol(pfm.Sigma,'upper'); % Sigma = Omega'*Omega
end
pfm.number_of_shocks = length(pfm.Sigma);
pfm.stochastic_order = options_.ep.stochastic.order;
pfm.max_lag = M_.maximum_endo_lag;
if pfm.max_lag > 0
    pfm.nyp = nnz(pfm.lead_lag_incidence(1,:));
    pfm.iyp = find(pfm.lead_lag_incidence(1,:)>0);
else
    pfm.nyp = 0;
    pfm.iyp = [];
end
pfm.ny0 = nnz(pfm.lead_lag_incidence(pfm.max_lag+1,:));
pfm.iy0 = find(pfm.lead_lag_incidence(pfm.max_lag+1,:)>0);
if M_.maximum_endo_lead
    pfm.nyf = nnz(pfm.lead_lag_incidence(pfm.max_lag+2,:));
    pfm.iyf = find(pfm.lead_lag_incidence(pfm.max_lag+2,:)>0);
else
    pfm.nyf = 0;
    pfm.iyf = [];
end
pfm.nd = pfm.nyp+pfm.ny0+pfm.nyf;
pfm.nrc = pfm.nyf+1;
pfm.isp = 1:pfm.nyp;
pfm.is = pfm.nyp+1:pfm.ny+pfm.nyp;
pfm.isf = pfm.iyf+pfm.nyp;
pfm.isf1 = pfm.nyp+pfm.ny+1:pfm.nyf+pfm.nyp+pfm.ny+1;
pfm.iz = 1:pfm.ny+pfm.nyp+pfm.nyf;
pfm.periods = options_.ep.periods;
pfm.steady_state = oo_.steady_state;
pfm.params = M_.params;
if M_.maximum_endo_lead
    pfm.i_cols_1 = nonzeros(pfm.lead_lag_incidence(pfm.max_lag+(1:2),:)');
    pfm.i_cols_A1 = find(pfm.lead_lag_incidence(pfm.max_lag+(1:2),:)');
else
    pfm.i_cols_1 = nonzeros(pfm.lead_lag_incidence(pfm.max_lag+1,:)');
    pfm.i_cols_A1 = find(pfm.lead_lag_incidence(pfm.max_lag+1,:)');
end
if pfm.max_lag > 0
    pfm.i_cols_T = nonzeros(pfm.lead_lag_incidence(1:2,:)');
else
    pfm.i_cols_T = nonzeros(pfm.lead_lag_incidence(1,:)');
end
pfm.i_cols_j = 1:pfm.nd;
pfm.i_upd = pfm.ny+(1:pfm.periods*pfm.ny);
if ~options_.bytecode
    pfm.dynamic_model = str2func([M_.fname,'.dynamic']);
end
pfm.verbose = options_.ep.verbosity;
pfm.maxit_ = options_.simul.maxit;
pfm.tolerance = options_.dynatol.f;
