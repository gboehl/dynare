function [h_minus_1, h, h_plus_1, h_exo, resid] = get_deriv(M_, ys_)
% function [h_minus_1, h, h_plus_1, h_exo, resid] = get_deriv(M_, ys_)
% Computes dynamic Jacobian and writes it to conformable matrices
%
% INPUTS
% - M_         [struct]     Definition of the model.
% - ys         vector       steady state
%
% OUTPUTS
% - h_minus_1  [N by N]     derivative matrix with respect to lagged endogenous variables
% - h          [N by N]     derivative matrix with respect to contemporanous endogenous variables
% - h_plus_1   [N by N]     derivative matrix with respect to leaded endogenous variables
% - h_exo      [N by N_exo] derivative matrix with respect to exogenous variables
% - resid      [N by 1]     vector of residuals

% Copyright (C) 2021 Dynare Team
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

x = zeros(M_.maximum_lag + M_.maximum_lead + 1,M_.exo_nbr);

iyv = M_.lead_lag_incidence';
iyr0 = find(iyv(:)) ;
z=repmat(ys_,1,M_.maximum_lag + M_.maximum_lead + 1);

[resid,g1]=feval([M_.fname,'.dynamic'],z(iyr0),x, M_.params, ys_, M_.maximum_exo_lag+1);

% Initialize matrices
h_minus_1=zeros(M_.endo_nbr);
h = h_minus_1;
h_plus_1 = h_minus_1;

% build h_minus_1
if M_.maximum_lag
    lag_columns=find(iyv(:,1));
    n_lag_columns=length(lag_columns);
    h_minus_1(:,lag_columns) = g1(:,1:n_lag_columns);
else
    n_lag_columns=0;
end

% build h
contemporaneous_columns=find(iyv(:,M_.maximum_lag+1));
n_contemporaneous_columns = length(contemporaneous_columns);
h(:,contemporaneous_columns) = g1(:,1+n_lag_columns:n_lag_columns+n_contemporaneous_columns);

%build h_plus_1
if M_.maximum_lead
    lead_columns=find(iyv(:,end));
    n_lead_columns = length(lead_columns);
    h_plus_1(:,lead_columns) = g1(:,n_lag_columns+n_contemporaneous_columns+1:n_lag_columns+n_contemporaneous_columns+n_lead_columns);
else
    n_lead_columns=0;
end

h_exo =g1(:,n_lag_columns+n_contemporaneous_columns+n_lead_columns+1:end);
