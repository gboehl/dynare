function [initialconditions, samplesize, innovations, DynareOptions, DynareModel, DynareOutput, nx, ny1, iy1, jdx, model_dynamic, y] = simul_backward_model_init(varargin)

% Initialization of the routines simulating backward models.    

% Copyright (C) 2012-2017 Dynare Team
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
    
initialconditions = varargin{1};
samplesize = varargin{2};

DynareOptions = varargin{3};
DynareModel = varargin{4};
DynareOutput = varargin{5};

if DynareModel.maximum_lead
    error('simul_backward_nonlinear_model:: The specified model is not backward looking!')
end

if nargin<6
    % Set the covariance matrix of the structural innovations.
    variances = diag(DynareModel.Sigma_e);
    number_of_shocks = length(DynareModel.Sigma_e);
    positive_var_indx = find(variances>0);
    effective_number_of_shocks = length(positive_var_indx);
    covariance_matrix = DynareModel.Sigma_e(positive_var_indx,positive_var_indx);
    covariance_matrix_upper_cholesky = chol(covariance_matrix);
    % Set seed to its default state.
    if DynareOptions.bnlms.set_dynare_seed_to_default
        set_dynare_seed('default');
    end
    % Simulate structural innovations.
    switch DynareOptions.bnlms.innovation_distribution
      case 'gaussian'
        DynareOutput.bnlms.shocks = randn(samplesize,effective_number_of_shocks)*covariance_matrix_upper_cholesky;
      otherwise
        error(['simul_backward_nonlinear_model:: ' DynareOption.bnlms.innovation_distribution ' distribution for the structural innovations is not (yet) implemented!'])
    end
    % Put the simulated innovations in DynareOutput.exo_simul.
    DynareOutput.exo_simul = zeros(samplesize,number_of_shocks);
    DynareOutput.exo_simul(:,positive_var_indx) = DynareOutput.bnlms.shocks;
    if isfield(DynareModel,'exo_histval') && ~ isempty(DynareModel.exo_histval)
        DynareOutput.exo_simul = [transpose(DynareModel.exo_histval); DynareOutput.exo_simul];
    else
        DynareOutput.exo_simul = [zeros(1,number_of_shocks); DynareOutput.exo_simul];
    end
    innovations = DynareOutput.exo_simul;
else
    innovations = varargin{6};
    DynareOutput.exo_simul = innovations; % innovations
end

if nargout>6
   nx = size(DynareOutput.exo_simul,2);
   ny0 = nnz(DynareModel.lead_lag_incidence(2,:));
   ny1 = nnz(DynareModel.lead_lag_incidence(1,:));
   iy1 = find(DynareModel.lead_lag_incidence(1,:)>0);
   idx = 1:DynareModel.endo_nbr;
   jdx = idx+ny1;
   % Get the name of the dynamic model routine.
   model_dynamic = str2func([DynareModel.fname,'_dynamic']);
   % initialization of vector y.
   y = NaN(length(idx)+ny1,1);
end