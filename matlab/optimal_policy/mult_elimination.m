function dr=mult_elimination(varlist,M_, options_, oo_)
% function mult_elimination()
% replaces Lagrange multipliers in Ramsey policy by lagged value of state
% and shock variables
%
% INPUT
%   none
%
% OUTPUT
%   dr: a structure with the new decision rule
%
% SPECIAL REQUIREMENTS
%   none

% Copyright © 2003-2024 Dynare Team
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

dr = oo_.dr;

nstatic = M_.nstatic;
nspred = M_.nspred;
order_var = dr.order_var;
nstates = M_.endo_names(order_var(nstatic+(1:nspred)));

il = strmatch('MULT_',nstates);
nil = setdiff(1:nspred,il);
m_nbr = length(il);
nm_nbr = length(nil);

AA1 = dr.ghx(:,nil);
AA2 = dr.ghx(:,il);
A1 = dr.ghx(nstatic+(1:nspred),nil);
A2 = dr.ghx(nstatic+(1:nspred),il);
B = dr.ghu(nstatic+(1:nspred),:);
A11 = A1(nil,:);
A21 = A1(il,:);
A12 = A2(nil,:);
A22 = A2(il,:);
B1 = B(nil,:);
B2 = B(il,:);

[Q1,R1] = qr([A12; A22]);
n1 = sum(abs(diag(R1)) > 1e-8);

Q1_12 = Q1(1:nm_nbr,n1+1:end);
Q1_22 = Q1(nm_nbr+(1:m_nbr),n1+1:end);
[Q2,R2,E2] = qr(Q1_22');
n2 = sum(abs(diag(R2)) > 1e-8);

R2_1 = inv(R2(1:n2,1:n2));

M1 = AA1 - AA2*E2*[R2_1*Q2(:,1:n2)'*Q1_12'; zeros(m_nbr-n2,nm_nbr)];
M2 = AA2*E2*[R2_1*Q2(:,1:n2)'*[Q1_12' Q1_22']*[A11;A21]; zeros(m_nbr-n2,length(nil))];
M3 = dr.ghu;
M4 = AA2*E2*[R2_1*Q2(:,1:n2)'*[Q1_12' Q1_22']*[B1;B2]; zeros(m_nbr-n2,size(B,2))];

k1 = nstatic+(1:nspred);
k1 = k1(nil);

endo_nbr = M_.orig_endo_nbr;
exo_nbr = M_.exo_nbr;

lead_lag_incidence = M_.lead_lag_incidence(:,1:endo_nbr+exo_nbr);
lead_lag_incidence1 = lead_lag_incidence(1,:) > 0;
maximum_lag = M_.maximum_lag;
for i=1:maximum_lag-1
    lead_lag_incidence1 = [lead_lag_incidence1; lead_lag_incidence(i,:)| ...
                        lead_lag_incidence(i+1,:)];
end
lead_lag_incidence1 = [lead_lag_incidence1; ...
                    lead_lag_incidence(M_.maximum_lag,:) > 0];
k = find(lead_lag_incidence1');
lead_lag_incidence1 = zeros(size(lead_lag_incidence1'));
lead_lag_incidence1(k) = 1:length(k);
lead_lag_incidence1 = lead_lag_incidence1';

dr.M1 = M1;
dr.M2 = M2;
dr.M3 = M3;
dr.M4 = M4;

nvar = length(varlist);

if nvar > 0 && ~options_.noprint
    res_table = zeros(2*(nm_nbr+M_.exo_nbr),nvar);
    headers = {'Variables'};
    for i=1:length(varlist)
        k = strmatch(varlist{i}, M_.endo_names(dr.order_var), 'exact');
        headers = vertcat(headers, varlist{i});

        res_table(1:nm_nbr,i) = M1(k,:)';
        res_table(nm_nbr+(1:nm_nbr),i) = M2(k,:)';
        res_table(2*nm_nbr+(1:M_.exo_nbr),i) = M3(k,:)';
        res_table(2*nm_nbr+M_.exo_nbr+(1:M_.exo_nbr),i) = M4(k,:)';
    end

    my_title='ELIMINATION OF THE MULTIPLIERS';
    lab = nstates(nil);
    labels = cellfun(@(x) horzcat(x, '(-1)'), nstates(nil), 'UniformOutput', false);
    labels = vertcat(labels, cellfun(@(x) horzcat(x, '(-2)'), nstates(nil), 'UniformOutput', false));
    labels = vertcat(labels, M_.exo_names);
    labels = vertcat(labels, cellfun(@(x) horzcat(x, '(-1)'), M_.exo_names, 'UniformOutput', false));
    lh = size(labels,2)+2;
    dyntable(options_, my_title, headers, labels, res_table, lh, 10, 6);
    skipline()
end
