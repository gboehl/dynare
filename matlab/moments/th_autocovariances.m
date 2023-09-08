function [Gamma_y,stationary_vars] = th_autocovariances(dr,ivar,M_,options_,nodecomposition)
% Computes the theoretical auto-covariances, Gamma_y, for an AR(p) process
% with coefficients dr.ghx and dr.ghu and shock variances Sigma_e
% for a subset of variables ivar.
% Theoretical HP-filtering and band-pass filtering is available as an option
%
% INPUTS
%   dr:               [structure]    Reduced form solution of the DSGE model  (decisions rules)
%   ivar:             [integer]      Vector of indices for a subset of variables.
%   M_                [structure]    Global dynare's structure, description of the DSGE model.
%   options_          [structure]    Global dynare's structure.
%   nodecomposition   [integer]      Scalar, if different from zero the variance decomposition is not triggered.
%
% OUTPUTS
%   Gamma_y           [cell]         Matlab cell of nar+3 (second order approximation) or nar+2 (first order approximation) arrays,
%                                    where nar is the order of the autocorrelation function.
%                                      Gamma_y{1}       [double]  Covariance matrix.
%                                      Gamma_y{i+1}     [double]  Autocorrelation function (for i=1,...,options_.nar).
%                                      Gamma_y{nar+2}   [double]  Variance decomposition.
%                                      Gamma_y{nar+3}   [double]  Expectation of the endogenous variables associated with a second
%                                                                 order approximation.
%   stationary_vars   [integer]      Vector of indices of stationary variables (as a subset of 1:length(ivar))
%
% SPECIAL REQUIREMENTS
%
% Algorithms
%   The means at order=2 are based on the pruned state space as
%   in Kim, Kim, Schaumburg, Sims (2008): Calculating and using second-order accurate
%   solutions of discrete time dynamic equilibrium models.
%   The solution at second order can be written as:
%   \[
%   \hat x_t = g_x \hat x_{t - 1} + g_u u_t + \frac{1}{2}\left( g_{\sigma\sigma} \sigma^2 + g_{xx}\hat x_t^2 + g_{uu} u_t^2 \right)
%   \]
%   Taking expectations on both sides requires to compute E(x^2)=Var(x), which
%   can be obtained up to second order from the first order solution
%   \[
%       \hat x_t = g_x \hat x_{t - 1} + g_u u_t
%   \]
%   by solving the corresponding Lyapunov equation.
%   Given Var(x), the above equation can be solved for E(x_t) as
%   \[
%   E(x_t) = (I - {g_x}\right)^{- 1} 0.5\left( g_{\sigma\sigma} \sigma^2 + g_{xx} Var(\hat x_t) + g_{uu} Var(u_t) \right)
%   \]
%
% Copyright © 2001-2023 Dynare Team
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

if nargin<5
    nodecomposition = 0;
end

if options_.order >= 3
    error('Theoretical moments not implemented above 2nd order')
end

local_order = options_.order;
if local_order~=1 && M_.hessian_eq_zero
    local_order = 1;
end

if local_order>1 && (options_.hp_filter || options_.bandpass.indicator)
    error('Theoretical filtered moments not implemented above 1st order')
end

endo_nbr = M_.endo_nbr;
if isoctave
    warning('off', 'Octave:divide-by-zero')
end
nar = options_.ar;
Gamma_y = cell(nar+2,1);
if isempty(ivar)
    ivar = (1:endo_nbr)';
end
nvar = size(ivar,1);

ghx = dr.ghx;
ghu = dr.ghu;

inv_order_var = dr.inv_order_var;
index_states = M_.nstatic+(1:M_.nspred)';
ghu_states_only = zeros(M_.nspred,M_.exo_nbr);
ghu_states_only(1:M_.nspred,:) = ghu(index_states,:); %get shock impact on states only

% state space representation for state variables only
[A,B] = kalman_transition_matrix(dr,index_states,1:M_.nspred);
% Compute stationary variables for unfiltered moments (filtering will remove unit roots)
if options_.hp_filter ~= 0 || options_.bandpass.indicator
    % By construction, all variables are stationary when filtered
    index_stationary_vars = inv_order_var(ivar);
    stationary_vars = (1:length(ivar))';
else
    [variance_states, unit_root_Schur_vector] =  lyapunov_symm(A,B*M_.Sigma_e*B',options_.lyapunov_fixed_point_tol,options_.qz_criterium,options_.lyapunov_complex_threshold,[],options_.debug);
    index_stationary_vars = inv_order_var(ivar);
    stationary_vars = (1:length(ivar))';
    if ~isempty(unit_root_Schur_vector)
        x = abs(ghx*unit_root_Schur_vector);
        stationary_vars = find(all(x(inv_order_var(ivar),:) < options_.schur_vec_tol,2));
        index_stationary_vars = inv_order_var(ivar(stationary_vars));
    end
end

%compute 2nd order mean correction on stationary variables 
if local_order == 2         % mean correction for 2nd order with no filters; other cases are error out above
    if ~isempty(index_states)
        Ex = (dr.ghs2(index_states)+dr.ghxx(index_states,:)*variance_states(:)+dr.ghuu(index_states,:)*M_.Sigma_e(:))/2;
        Ex = (eye(M_.nspred)-ghx(index_states,:))\Ex;
        Gamma_y{nar+3} = NaN*ones(nvar, 1);
        Gamma_y{nar+3}(stationary_vars) = ghx(index_stationary_vars,:)*Ex+(dr.ghs2(index_stationary_vars)+dr.ghxx(index_stationary_vars,:)*variance_states(:)+...
            dr.ghuu(index_stationary_vars,:)*M_.Sigma_e(:))/2;
    else %no static and no predetermined
        Gamma_y{nar+3} = NaN*ones(nvar, 1);
        Gamma_y{nar+3}(stationary_vars) = (dr.ghs2(index_stationary_vars)+ dr.ghuu(index_stationary_vars,:)*M_.Sigma_e(:))/2;
    end
end

if options_.hp_filter == 0 && ~options_.bandpass.indicator
    aa = ghx(index_stationary_vars,:);
    bb = ghu(index_stationary_vars,:);
    % unconditional variance
    v = NaN*ones(nvar,nvar);
    v(stationary_vars,stationary_vars) = aa*variance_states*aa'+ bb*M_.Sigma_e*bb';
    v(abs(v) < 1e-12) = 0;
    Gamma_y{1} = v;
    % autocorrelations
    if nar > 0
        vxy = (A*variance_states*aa'+ghu_states_only*M_.Sigma_e*bb');
        sy = sqrt(diag(Gamma_y{1}));
        sy = sy(stationary_vars);
        sy = sy *sy';
        Gamma_y{2} = NaN*ones(nvar,nvar);
        Gamma_y{2}(stationary_vars,stationary_vars) = aa*vxy./sy;
        for i=2:nar
            vxy = A*vxy;
            Gamma_y{i+1} = NaN*ones(nvar,nvar);
            Gamma_y{i+1}(stationary_vars,stationary_vars) = aa*vxy./sy;
        end
    end
    % variance decomposition
    if ~nodecomposition && M_.exo_nbr > 0 && size(stationary_vars, 1) > 0
        Gamma_y{nar+2}=compute_variance_decomposition(M_,options_,diag(v(stationary_vars,stationary_vars)),A,aa,ghu,ghu_states_only,stationary_vars,index_stationary_vars,nvar);
    end
else% ==> Theoretical filters.
    aa = ghx(index_stationary_vars,:); %R in Uhlig (2001)
    bb = ghu(index_stationary_vars,:); %S in Uhlig (2001)

    ngrid = options_.filtered_theoretical_moments_grid;
    freqs = 0 : ((2*pi)/ngrid) : (2*pi*(1 - .5/ngrid)); %[0,2*pi)
    tpos  = exp( sqrt(-1)*freqs); %positive frequencies
    tneg  =  exp(-sqrt(-1)*freqs); %negative frequencies
    if options_.bandpass.indicator
        filter_gain = zeros(1,ngrid);
        lowest_periodicity=options_.bandpass.passband(2);
        highest_periodicity=options_.bandpass.passband(1);
        highest_periodicity=max(2,highest_periodicity); % restrict to upper bound of pi
        filter_gain(freqs>=2*pi/lowest_periodicity & freqs<=2*pi/highest_periodicity)=1;
        filter_gain(freqs<=-2*pi/lowest_periodicity+2*pi & freqs>=-2*pi/highest_periodicity+2*pi)=1;
    else
        lambda = options_.hp_filter;
        filter_gain = 4*lambda*(1 - cos(freqs)).^2 ./ (1 + 4*lambda*(1 - cos(freqs)).^2);   %HP transfer function
    end
    mathp_col = NaN(ngrid,length(ivar)^2);
    IA = eye(size(A,1));
    IE = eye(M_.exo_nbr);
    for ig = 1:ngrid
        if filter_gain(ig)==0
            f_hp = zeros(length(ivar),length(ivar));
        else
            f_omega  =(1/(2*pi))*([(IA-A*tneg(ig))\ghu_states_only;IE]...
                                  *M_.Sigma_e*[ghu_states_only'/(IA-A'*tpos(ig)) IE]); % spectral density of state variables; top formula Uhlig (2001), p. 20 with N=0
            g_omega = [aa*tneg(ig) bb]*f_omega*[aa'*tpos(ig); bb']; % spectral density of selected variables; middle formula Uhlig (2001), p. 20; only middle block, i.e. y_t'
            f_hp = filter_gain(ig)^2*g_omega; % spectral density of selected filtered series; top formula Uhlig (2001), p. 21;
        end
        mathp_col(ig,:) = (f_hp(:))';    % store as matrix row for ifft
    end
    % Covariance of filtered series
    imathp_col = real(ifft(mathp_col))*(2*pi); % Inverse Fast Fourier Transformation; middle formula Uhlig (2001), p. 21;
    Gamma_y{1} = reshape(imathp_col(1,:),nvar,nvar);
    % Autocorrelations
    if nar > 0
        sy = sqrt(diag(Gamma_y{1}));
        sy = sy *sy';
        for i=1:nar
            Gamma_y{i+1} = reshape(imathp_col(i+1,:),nvar,nvar)./sy;
        end
    end
    % Variance decomposition
    if ~nodecomposition && M_.exo_nbr > 0
        if M_.exo_nbr == 1
            Gamma_y{nar+2} = ones(nvar,1);
        else
            Gamma_y{nar+2} = zeros(nvar,M_.exo_nbr);
            cs = get_lower_cholesky_covariance(M_.Sigma_e); %make sure Covariance matrix is positive definite
            SS = cs*cs';
            b1 = ghu_states_only;
            b2 = ghu(index_stationary_vars,:);
            mathp_col = NaN(ngrid,length(ivar)^2);
            IA = eye(size(A,1));
            IE = eye(M_.exo_nbr);
            for ig = 1:ngrid
                if filter_gain(ig)==0
                    f_hp = zeros(length(ivar),length(ivar));
                else
                    f_omega  =(1/(2*pi))*( [(IA-A*tneg(ig))\b1;IE]...
                                           *SS*[b1'/(IA-A'*tpos(ig)) IE]); % spectral density of state variables; top formula Uhlig (2001), p. 20 with N=0
                    g_omega = [aa*tneg(ig) b2]*f_omega*[aa'*tpos(ig); b2']; % spectral density of selected variables; middle formula Uhlig (2001), p. 20; only middle block, i.e. y_t'
                    f_hp = filter_gain(ig)^2*g_omega;  % spectral density of selected filtered series; top formula Uhlig (2001), p. 21;
                end
                mathp_col(ig,:) = (f_hp(:))';    % store as matrix row for ifft
            end
            imathp_col = real(ifft(mathp_col))*(2*pi);
            vv = diag(reshape(imathp_col(1,:),nvar,nvar));
            for i=1:M_.exo_nbr
                mathp_col = NaN(ngrid,length(ivar)^2);
                SSi = cs(:,i)*cs(:,i)';
                for ig = 1:ngrid
                    if filter_gain(ig)==0
                        f_hp = zeros(length(ivar),length(ivar));
                    else
                        f_omega  =(1/(2*pi))*( [(IA-A*tneg(ig))\b1;IE]...
                                               *SSi*[b1'/(IA-A'*tpos(ig)) IE]); % spectral density of state variables; top formula Uhlig (2001), p. 20 with N=0
                        g_omega = [aa*tneg(ig) b2]*f_omega*[aa'*tpos(ig); b2']; % spectral density of selected variables; middle formula Uhlig (2001), p. 20; only middle block, i.e. y_t'
                        f_hp = filter_gain(ig)^2*g_omega; % spectral density of selected filtered series; top formula Uhlig (2001), p. 21;
                    end
                    mathp_col(ig,:) = (f_hp(:))';    % store as matrix row for ifft
                end
                imathp_col = real(ifft(mathp_col))*(2*pi);
                Gamma_y{nar+2}(:,i) = abs(diag(reshape(imathp_col(1,:),nvar,nvar)))./vv;
            end
        end
    end
end
if isoctave
    warning('on', 'Octave:divide-by-zero')
end