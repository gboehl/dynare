function [dr, info] = dyn_first_order_solver(jacobia, M_, dr, options_, task)
% [dr, info] = dyn_first_order_solver(jacobia, M_, dr, options_, task)
% Computes the first order reduced form of a DSGE model.
%
% INPUTS
% - jacobia       [double]    matrix, the jacobian of the dynamic model.
% - M_            [struct]    Matlab's structre describing the model
% - dr            [struct]    Matlab's structure describing the reduced form model.
% - options_      [struct]    Matlab's structure containing the current state of the options
% - task          [integer]   scalar, if task = 0 then decision rules are computed and if task = 1 then only eigenvales are computed.
%
% OUTPUTS
% - dr            [struct]    Matlab's structure describing the reduced form model.
% - info          [integer]   scalar, error code. Possible values are:
%
%                                     info=0 -> no error,
%                                     info=1 -> the model doesn't determine the current variables uniquely,
%                                     info=2 -> mjdgges dll returned an error,
%                                     info=3 -> Blanchard and Kahn conditions are not satisfied: no stable equilibrium,
%                                     info=4 -> Blanchard and Kahn conditions are not satisfied: indeterminacy,
%                                     info=5 -> Blanchard and Kahn conditions are not satisfied: indeterminacy due to rank failure,
%                                     info=7 -> One of the eigenvalues is close to 0/0 (infinity of complex solutions)

% Copyright © 2001-2024 Dynare Team
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

persistent reorder_jacobian_columns innovations_idx index_s index_m index_c
persistent index_p row_indx index_0m index_0p k1 k2 state_var
persistent ndynamic nstatic nfwrd npred nboth nd nsfwrd n_current index_d
persistent index_e index_d1 index_d2 index_e1 index_e2 row_indx_de_1
persistent row_indx_de_2 cols_b


if ~nargin
    if nargout
        error('dyn_first_order_solver:: Initialization mode returns zero argument!')
    end
    reorder_jacobian_columns = [];
    return
end

exo_nbr = M_.exo_nbr;

if isempty(reorder_jacobian_columns)

    maximum_lag = M_.maximum_endo_lag;
    nfwrd    = M_.nfwrd;
    nboth    = M_.nboth;
    npred    = M_.npred;
    nstatic  = M_.nstatic;
    ndynamic = M_.ndynamic;
    nsfwrd   = M_.nsfwrd;

    k1 = 1:(npred+nboth);
    k2 = 1:(nfwrd+nboth);

    order_var = dr.order_var;
    nd = npred+nfwrd+2*nboth;
    lead_lag_incidence = M_.lead_lag_incidence;
    nz = nnz(lead_lag_incidence);

    lead_id = find(lead_lag_incidence(maximum_lag+2,:));
    if maximum_lag
        lag_id = find(lead_lag_incidence(1,:));
    else
        lag_id = [];
    end

    both_id = intersect(lead_id,lag_id);
    if maximum_lag
        no_both_lag_id = setdiff(lag_id,both_id);
    else
        no_both_lag_id = lag_id;
    end
    innovations_idx = nz+(1:exo_nbr);
    state_var  = [no_both_lag_id, both_id];

    index_c  = nonzeros(lead_lag_incidence(maximum_lag+1,:));             % Index of all endogenous variables present at time=t
    n_current = length(index_c);

    index_s  = npred+nboth+(1:nstatic);     % Index of all static
                                            % endogenous variables
                                            % present at time=t
    indexi_0 = npred+nboth;

    npred0 = nnz(lead_lag_incidence(maximum_lag+1,no_both_lag_id));
    index_0m = indexi_0+nstatic+(1:npred0);
    nfwrd0 = nnz(lead_lag_incidence(maximum_lag+1,lead_id));
    index_0p = indexi_0+nstatic+npred0+(1:nfwrd0);
    index_m  = 1:(npred+nboth);
    index_p  = npred+nboth+n_current+(1:nfwrd+nboth);
    row_indx_de_1 = 1:ndynamic;
    row_indx_de_2 = ndynamic+(1:nboth);
    row_indx = nstatic+row_indx_de_1;
    index_d = [index_0m index_p];
    llx = lead_lag_incidence(:,order_var);
    index_d1 = [find(llx(maximum_lag+1,nstatic+(1:npred))), npred+nboth+(1:nsfwrd) ];
    index_d2 = npred+(1:nboth);

    index_e = [index_m index_0p];
    index_e1 = [1:npred+nboth, npred+nboth+find(llx(maximum_lag+1,nstatic+npred+(1: ...
                                                      nsfwrd)))];
    index_e2 = npred+nboth+(1:nboth);

    [~,cols_b] = find(lead_lag_incidence(maximum_lag+1, order_var));

    reorder_jacobian_columns = [nonzeros(lead_lag_incidence(:,order_var)'); ...
                        nz+(1:exo_nbr)'];
end

dr.ghx = [];
dr.ghu = [];

dr.state_var = state_var;

jacobia = jacobia(:,reorder_jacobian_columns);

if nstatic > 0
    [Q, ~] = qr(jacobia(:,index_s));
    aa = Q'*jacobia;
else
    aa = jacobia;
end

A = aa(:,index_m);  % Jacobian matrix for lagged endogenous variables
B(:,cols_b) = aa(:,index_c);  % Jacobian matrix for contemporaneous endogenous variables
C = aa(:,index_p);  % Jacobian matrix for led endogenous variables

info = 0;
if task ~= 1 && (options_.dr_cycle_reduction || options_.dr_logarithmic_reduction)
    if n_current < M_.endo_nbr
        if options_.dr_cycle_reduction
            error(['The cycle reduction algorithm can''t be used when the ' ...
                   'coefficient matrix for current variables isn''t invertible'])
        elseif options_.dr_logarithmic_reduction
            error(['The logarithmic reduction algorithm can''t be used when the ' ...
                   'coefficient matrix for current variables isn''t invertible'])
        end
    end
    A1 = [aa(row_indx,index_m ) zeros(ndynamic,nfwrd)];
    B1 = [aa(row_indx,index_0m) aa(row_indx,index_0p) ];
    C1 = [zeros(ndynamic,npred) aa(row_indx,index_p)];
    if options_.dr_cycle_reduction
        [ghx, info] = cycle_reduction(A1, B1, C1, options_.dr_cycle_reduction_tol);
    else
        [ghx, info] = logarithmic_reduction(C1, B1, A1, options_.dr_logarithmic_reduction_tol, options_.dr_logarithmic_reduction_maxiter);
    end
    if info
        % cycle_reduction or logarithmic reduction failed and set info
        return
    end
    ghx = ghx(:,index_m);
    hx = ghx(1:npred+nboth,:);
    gx = ghx(1+npred:end,:);
else
    D = zeros(ndynamic+nboth);
    E = zeros(ndynamic+nboth);
    D(row_indx_de_1,index_d1) = aa(row_indx,index_d);
    D(row_indx_de_2,index_d2) = eye(nboth);
    E(row_indx_de_1,index_e1) = -aa(row_indx,index_e);
    E(row_indx_de_2,index_e2) = eye(nboth);

    [ss, tt, w, sdim, dr.eigval, info1] = mjdgges(E, D, options_.qz_criterium, options_.qz_zero_threshold);

    if info1
        if info1 == -30
            % one eigenvalue is close to 0/0
            info(1) = 7;
        else
            info(1) = 2;
            info(2) = info1;
            info(3) = size(E,2);
        end
        return
    end

    dr.sdim = sdim;                      % Number of stable eigenvalues.
    dr.edim = length(dr.eigval)-sdim;    % Number of exposive eigenvalues.

    nba = nd-sdim;

    if task==1
        if rcond(w(npred+nboth+1:end,npred+nboth+1:end)) < 1e-9
            dr.full_rank = 0;
        else
            dr.full_rank = 1;
        end
    end

    if nba ~= nsfwrd
        temp = sort(abs(dr.eigval));
        if nba > nsfwrd
            temp = temp(nd-nba+1:nd-nsfwrd)-1-options_.qz_criterium;
            info(1) = 3;
        elseif nba < nsfwrd
            temp = temp(nd-nsfwrd+1:nd-nba)-1-options_.qz_criterium;
            info(1) = 4;
        end
        info(2) = temp'*temp;
        return
    end

    if task==1, return, end

    %First order approximation
    indx_stable_root = 1: (nd - nsfwrd);         %=> index of stable roots
    indx_explosive_root = npred + nboth + 1:nd;  %=> index of explosive roots
                                                 % derivatives with respect to dynamic state variables
                                                 % forward variables
    Z = w';
    Z11 = Z(indx_stable_root,    indx_stable_root);
    Z21  = Z(indx_explosive_root, indx_stable_root);
    Z22  = Z(indx_explosive_root, indx_explosive_root);
    opts.TRANSA = false; % needed by Octave 4.0.0
    [minus_gx,rc] = linsolve(Z22,Z21,opts);
    if rc < 1e-9
        % Z22 is near singular
        info(1) = 5;
        info(2) = -log(rc);
        return
    end
    gx  = -minus_gx;
    % predetermined variables
    opts.UT = true;
    opts.TRANSA = true;
    hx1 = linsolve(tt(indx_stable_root, indx_stable_root),Z11,opts)';
    opts.UT = false;      % needed by Octave 4.0.0
    opts.TRANSA = false;  % needed by Octave 4.0.0
    hx2 = linsolve(Z11,ss(indx_stable_root, indx_stable_root)',opts)';
    hx =  hx1*hx2;
    ghx = [hx(k1,:); gx(k2(nboth+1:end),:)];
end

if nstatic > 0
    B_static = B(:,1:nstatic);  % submatrix containing the derivatives w.r. to static variables
else
    B_static = [];
end
%static variables, backward variable, mixed variables and forward variables
B_pred = B(:,nstatic+1:nstatic+npred+nboth);
B_fyd = B(:,nstatic+npred+nboth+1:end);

% static variables
if nstatic > 0
    temp = - C(1:nstatic,:)*gx*hx;
    b(:,cols_b) = aa(:,index_c);
    b10 = b(1:nstatic, 1:nstatic);
    b11 = b(1:nstatic, nstatic+1:end);
    temp(:,index_m) = temp(:,index_m)-A(1:nstatic,:);
    temp = b10\(temp-b11*ghx);
    if options_.debug
        if any(any(~isfinite(temp)))
            fprintf('\ndyn_first_order_solver: infinite/NaN elements encountered when solving for the static variables\n')
            fprintf('dyn_first_order_solver: This often arises if there is a singularity.\n\n')
        end
    end
    ghx = [temp; ghx];
end

A_ = real([B_static C*gx+B_pred B_fyd]); % The state_variable of the block are located at [B_pred B_both]

if exo_nbr
    if nstatic > 0
        fu = Q' * jacobia(:,innovations_idx);
    else
        fu = jacobia(:,innovations_idx);
    end

    ghu = - A_ \ fu;
else
    ghu = [];
end

dr.ghx = ghx;
dr.ghu = ghu;

if options_.aim_solver ~= 1
    % Necessary when using Sims' routines for QZ
    dr.ghx = real(ghx);
    dr.ghu = real(ghu);
    hx = real(hx);
end

% non-predetermined variables
dr.gx = gx;
%predetermined (endogenous state) variables, square transition matrix
dr.Gy = hx;
