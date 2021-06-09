function pruned_state_space = pruned_state_space_system(M, options, dr, indy, nlags, useautocorr, compute_derivs)
% Set up the pruned state space ABCD representation:
%   z =      c + A*z(-1) + B*inov
%   y = ys + d + C*z(-1) + D*inov
% References: 
% - Andreasen, Martin M., Jesús Fernández-Villaverde and Juan F. Rubio-Ramírez (2018):
%   "The Pruned State-Space System for Non-Linear DSGE Models: Theory and Empirical Applications",
%   Review of Economic Studies, Volume 85, Issue 1, Pages 1–49.
% - Mutschler, Willi (2018): "Higher-order statistics for DSGE models",
%   Econometrics and Statistics, Volume 6, Pages 44-56.
% =========================================================================
% INPUTS
%   M:              [structure] storing the model information
%   options:        [structure] storing the options
%   dr:             [structure] storing the results from perturbation approximation
%   indy:           [vector]    index of control variables in DR order
%   nlags:          [integer]   number of lags in autocovariances and autocorrelations
%   useautocorr:    [boolean]   true: compute autocorrelations
%   compute_derivs: [boolean]   true: compute derivatives
% -------------------------------------------------------------------------
% OUTPUTS
% pruned_state_space: [structure] with the following fields:
%   indx:       [x_nbr by 1]
%                 index of state variables
%   indy:       [y_nbr by 1]
%                 index of control variables
%   A:          [z_nbr by z_nbr]
%                 state space transition matrix A mapping previous states to current states
%   B:          [z_nbr by inov_nbr]
%                 state space transition matrix B mapping current inovations to current states
%   c:          [z_nbr by 1]
%                 state space transition matrix c mapping constants to current states
%   C:          [y_nbr by z_nbr]
%                 state space measurement matrix C mapping previous states to current controls
%   D:          [y_nbr by inov_nbr]
%                 state space measurement matrix D mapping current inovations to current controls
%   d:          [y_nbr by 1]
%                 state space measurement matrix d mapping constants to current controls
%   Var_inov    [inov_nbr by inov_nbr]
%                 contemporenous covariance matrix of innovations, i.e. E[inov*inov']
%   Var_z       [z_nbr by z_nbr]
%                 contemporenous covariance matrix of states z
%   Var_y       [y_nbr by y_nbr]
%                 contemporenous covariance matrix of controls y
%   Var_yi      [y_nbr by y_nbr by nlags]
%                 autocovariance matrix of controls y
%   Corr_y      [y_nbr by y_nbr]
%                 contemporenous correlation matrix of controls y
%   Corr_yi     [y_nbr by y_nbr by nlags]
%                 autocorrelation matrix of controls y
%   E_y         [y_nbr by 1] 
%                 unconditional theoretical mean of control variables y
%
% if compute_derivs == 1, then the following additional fields are outputed:
%   dA:         [z_nbr by z_nbr by totparam_nbr]
%                 parameter Jacobian of A
%   dB:         [z_nbr by inov_nbr by totparam_nbr]
%                 parameter Jacobian of B
%   dc:         [z_nbr by totparam_nbr]
%                 parameter Jacobian of c
%   dC:         [y_nbr by z_nbr by totparam_nbr]
%                 parameter Jacobian of C
%   dD:         [y_nbr by inov_nbr by totparam_nbr]
%                 parameter Jacobian of D
%   dd:         [y_nbr by totparam_nbr]
%                 parameter Jacobian of d
%   dVar_inov   [inov_nbr by inov_nbr by totparam_nbr]
%                 parameter Jacobian of Var_inov
%   dVar_z      [z_nbr by z_nbr by totparam_nbr]
%                 parameter Jacobian of Var_z
%   dVar_y      [y_nbr by y_nbr by totparam_nbr]
%                 parameter Jacobian of Var_y
%   dVar_yi     [y_nbr by y_nbr by nlags by totparam_nbr]
%                 parameter Jacobian of Var_yi
%   dCorr_y     [y_nbr by y_nbr by totparam_nbr]
%                 parameter Jacobian of Corr_y
%   dCorr_yi    [y_nbr by y_nbr by nlags by totparam_nbr]
%                 parameter Jacobian of Corr_yi
%   dE_y        [y_nbr by totparam_nbr]
%                 parameter Jacobian of E_y
% -------------------------------------------------------------------------
% This function is called by
%   * get_identification_jacobians.m
%   * identification_numerical_objective.m
% -------------------------------------------------------------------------
% This function calls
%   * allVL1.m
%   * commutation.m
%   * disclyap_fast (MEX)
%   * duplication.m
%   * lyapunov_symm.m
%   * prodmom
%   * prodmom_deriv
%   * Q6_plication
%   * quadruplication.m
%   * vec.m
% =========================================================================
% Copyright (C) 2019-2020 Dynare Team
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
% =========================================================================

%% MAIN IDEA:
%   Decompose the state vector x into first-order effects xf, second-order 
%   effects xs, and third-order effects xrd, i.e. x=xf+xs+xrd. Then, Dynare's
%   perturbation approximation for the state vector up to third order
%   (with Gaussian innovations u, i.e. no odd moments, hxxs=huus=hxus=hsss=0) is:
%   x = hx*( xf(-1)+xs(-1)+xrd(-1) )
%     + hu*u
%     + 1/2*hxx*kron( xf(-1)+xs(-1)+xrd(-1) , xf(-1)+xs(-1)+xrd(-1) )
%     + hxu*kron( xf(-1)+xs(-1)+xrd(-1) , u )
%     + 1/2*huu*kron( u , u )
%     + 1/2*hss*sig^2
%     + 1/6*hxxx*kron( kron( xf(-1)+xs(-1)+xrd(-1) , xf(-1)+xs(-1)+xrd(-1) ) , xf(-1)+xs(-1)+xrd(-1) )
%     + 1/6*huuu*kron( kron( u , u ) , u )
%     + 3/6*hxxu*kron( kron( xf(-1)+xs(-1)+xrd(-1) , xf(-1)+xs(-1)+xrd(-1) ) , u )
%     + 3/6*hxuu*kron( kron( xf(-1)+xs(-1)+xrd(-1) , u ) , u)
%     + 3/6*hxss*( xf(-1)+xs(-1)+xrd(-1) )*sig^2
%     + 3/6*huss*u*sig^2
%   where:
%     hx  = dr.ghx(indx,:);    hu  = dr.ghu(indx,:);
%     hxx = dr.ghxx(indx,:);   hxu = dr.ghxu(indx,:);   huu = dr.ghuu(indx,:);   hss = dr.ghs2(indx,:);
%     hxxx = dr.ghxxx(indx,:); huuu = dr.ghuuu(indx,:); hxxu = dr.ghxxu(indx,:); hxuu = dr.ghxxu(indx,:); hxss = dr.ghxss(indx,:); huss = dr.ghuss(indx,:);
%     and similarly for control variables:
%   y = gx*( xf(-1)+xs(-1)+xrd(-1) )
%     + gu*u
%     + 1/2*gxx*kron( xf(-1)+xs(-1)+xrd(-1) , xf(-1)+xs(-1)+xrd(-1) )
%     + gxu*kron( xf(-1)+xs(-1)+xrd(-1) , u )
%     + 1/2*guu*kron( u , u )
%     + 1/2*gss*sig^2
%     + 1/6*gxxx*kron( kron( xf(-1)+xs(-1)+xrd(-1) , xf(-1)+xs(-1)+xrd(-1) ) , xf(-1)+xs(-1)+xrd(-1) )
%     + 1/6*guuu*kron( kron( u , u ) , u )
%     + 3/6*gxxu*kron( kron( xf(-1)+xs(-1)+xrd(-1) , xf(-1)+xs(-1)+xrd(-1) ) , u )
%     + 3/6*gxuu*kron( kron( xf(-1)+xs(-1)+xrd(-1) , u ) , u)
%     + 3/6*gxss*( xf(-1)+xs(-1)+xrd(-1) )*sig^2
%     + 3/6*guss*u*sig^2
%   where:
%     gx  = dr.ghx(indy,:);    gu  = dr.ghu(indy,:);
%     gxx = dr.ghxx(indy,:);   gxu = dr.ghxu(indy,:);   guu = dr.ghuu(indy,:);   gss = dr.ghs2(indy,:);
%     gxxx = dr.ghxxx(indy,:); guuu = dr.ghuuu(indy,:); gxxu = dr.ghxxu(indy,:); gxuu = dr.ghxxu(indy,:); gxss = dr.ghxss(indy,:); guss = dr.ghuss(indy,:);
%
%   PRUNING means getting rid of terms higher than the approximation order, i.e.
%         - involving fourth-order effects:  kron(xf,xrd), kron(xs,xs), kron(xrd,xf), kron(xrd,u), 
%                                            kron(kron(xf,xf),xs), kron(kron(xf,xs),xf), kron(kron(xs,xf),xf)
%                                            kron(kron(xf,xs),u), kron(kron(xs,xf),u)
%                                            kron(kron(xs,u),u)
%                                            xs*sig^2
%         - involving fifth-order effects:   kron(xs,xrd), kron(xrd,xs),
%                                            kron(kron(xf,xf),xrd), kron(kron(xf,xs),xs), kron(kron(xf,xrd),xf), kron(kron(xs,xf),xs), kron(kron(xs,xs),xf), kron(kron(xrd,xf),xf)
%                                            kron(kron(xf,xrd),u), kron(kron(xs,xs),u), kron(kron(xrd,xf),u)
%                                            kron(kron(xrd,u),u)
%                                            xrd*sig^2
%         - involving sixth-order effects:   kron(xrd,xrd),
%                                            kron(kron(xf,xs),xrd), kron(kron(xf,xrd),xs), kron(kron(xs,xrd),xrd), kron(kron(xs,xs),xs), kron(kron(xs,xrd),xf), kron(kron(xf,xf),xs), kron(kron(xrd,xs),xf)
%                                            kron(kron(xs,xrd),u), kron(kron(xrd,xs),u)
%         - involving seventh-order effects: kron(kron(xf,xrd),xrd), kron(kron(xs,xs),xrd), kron(kron(xs,xrd),xs), kron(kron(xrd,xf),xrd), kron(kron(xrd,xs),xs), kron(kron(xrd,xrd),xf)
%                                            kron(kron(xrd,xrd),u)
%         - involving eighth-order effects:  kron(kron(xs,xrd),xrd), kron(kron(xrd,xs),xrd), kron(kron(xrd,xrd),xs)
%         - involving ninth-order effects:   kron(kron(xrd,xrd),xrd)
%   Note that u is treated as a first-order effect and the perturbation parameter sig as a variable.
%
% SUMMARY OF LAW OF MOTIONS: Set up the law of motions for the individual effects, but keep only effects of same order
%   Notation: I_n=eye(n) and K_m_n=commutation(m,n)
%
%   First-order effects: keep xf and u
%       xf = hx*xf(-1) + hu*u
%       Note that we 
%
%   Second-order effects: keep xs, kron(xf,xf), kron(u,u), kron(xf,u), and sig^2
%       xs = hx*xs(-1) + 1/2*hxx*kron(xf(-1),xf(-1)) + 1/2*huu*(kron(u,u)-Sigma_e(:)+Sigma_e(:)) + hxu*kron(xf(-1),u) + 1/2*hss*sig^2%     
%
%   Third-order effects: keep xrd, kron(xf,xs), kron(xs,xf), kron(xs,u), kron(kron(xf,xf),xf), kron(kron(u,u),u), kron(kron(xf,xf),u), kron(kron(xf,u),u), xf*sig^2, u*sig^2
%       xrd = hx*xrd(-1) + 1/2*hxx*(kron(xf(-1),xs(-1))+kron(xs(-1),xf(-1))) + hxu*kron(xs(-1),u) + 1/6*hxxx*kron(xf(-1),kron(xf(-1),xf(-1))) + 1/6*huuu*kron(u,kron(u,u)) + 3/6*hxxu*kron(xf(-1),kron(xf(-1),u)) + 3/6*hxuu*kron(xf(-1),kron(u,u)) + 3/6*hxss*xf(-1)*sig^2 + 3/6*huss*u*sig^2
%     Simplified (due to symmetry in hxx):
%       xrd = hx*xrd(-1) + hxx*(kron(xf(-1),xs(-1)) + hxu*kron(xs(-1),u) + 1/6*hxxx*kron(xf(-1),kron(xf(-1),xf(-1))) + 1/6*huuu*kron(u,kron(u,u)) + 3/6*hxxu*kron(xf(-1),kron(xf(-1),u)) + 3/6*hxuu*kron(xf(-1),kron(u,u)) + 3/6*hxss*xf(-1)*sig^2 + 3/6*huss*u*sig^2%     
%
%   Auxiliary equation kron(xf,xf) to set up the VAR(1) pruned state space system
%       kron(xf,xf) = kron(hx,hx)*kron(xf(-1),xf(-1)) + kron(hu,hu)*(kron(u,u)-Sigma_e(:)+Sigma_e(:)) + kron(hx,u)*kron(xf(-1),u) + kron(u,hx)*kron(u,xf(-1))
%     Simplified using commutation matrix:
%       kron(xf,xf) = kron(hx,hx)*kron(xf(-1),xf(-1)) + (I_xx+K_x_x)*kron(hx,hu)*kron(xf(-1),u) + kron(hu,hu)*kron(u,u)
%
%   Auxiliary equation kron(xf,xs) to set up the VAR(1) pruned state space system
%       kron(xf,xs) = kron(hx,hx)*kron(xf(-1),xs(-1)) + kron(hu,hx)*kron(u,xs(-1))
%                   + kron(hx,1/2*hxx)*kron(kron(xf(-1),xf(-1)),xf(-1)) + kron(hu,1/2*hxx)*kron(kron(u,xf(-1)),xf(-1))
%                   + kron(hx,1/2*huu)*kron(kron(xf(-1),u),u) + kron(hu,1/2*huu)*kron(kron(u,u),u)
%                   + kron(hx,hxu)*kron(kron(xf(-1),xf(-1)),u) + kron(hu,hxu)*kron(kron(u,xf(-1)),u)
%                   + kron(hx,1/2*hss)*xf(-1)*sig^2 + kron(hu,1/2*hss)*u*sig^2
%     Simplified using commutation matrix:
%       kron(xf,xs) = kron(hx,hx)*kron(xf(-1),xs(-1))
%                   + K_x_x*kron(hx,hu)*kron(xs(-1),u)
%                   + kron(hx,1/2*hxx)*kron(kron(xf(-1),xf(-1)),xf(-1))
%                   + ( kron(hx,hxu) + K_x_x*kron(1/2*hxx,hu) )*kron(kron(xf(-1),xf(-1)),u)
%                   + ( kron(hx,1/2*huu) + K_x_x*kron(hxu,hu) )*kron(kron(xf(-1),u),u)
%                   + kron(hu,1/2*huu)*kron(kron(u,u),u)
%                   + kron(hx,1/2*hss)*xf(-1)*sig^2
%                   + kron(hu,1/2*hss)*u*sig^2
%
%   Auxiliary equation kron(kron(xf,xf),xf) to set up the VAR(1) pruned state space system
%       kron(kron(xf,xf),xf) = kron(kron(hx,hx),hx)*kron(kron(xf(-1),xf(-1)),xf(-1))
%                            + kron(kron(hx,hu),hx)*kron(kron(xf(-1),u),xf(-1))
%                            + kron(kron(hx,hx),hu)*kron(kron(xf(-1),xf(-1)),u)
%                            + kron(kron(hx,hu),hu)*kron(kron(xf(-1),u),u)
%                            + kron(kron(hu,hx),hx)*kron(kron(u,xf(-1)),xf(-1))
%                            + kron(kron(hu,hu),hx)*kron(kron(u,u),xf(-1))
%                            + kron(kron(hu,hx),hu)*kron(kron(u,xf(-1)),u)
%                            + kron(kron(hu,hu),hu)*kron(kron(u,u),u)
%     Simplified using commutation matrix:
%       kron(kron(xf,xf),xf) = kron(kron(hx,hx),hx)*kron(kron(xf(-1),xf(-1)),xf(-1))
%                            + ( kron(kron(hx,hx),hu) + K_xx_x*kron(hx,(I_xx+K_x_x)*kron(hx,hu)) )*kron(kron(xf(-1),xf(-1)),u)
%                            + ( kron((I_xx+K_x_x)*kron(hx,hu),hu) + K_xx_x*kron(kron(hx,hu),hu) )*kron(kron(xf(-1),u),u)
%                            + kron(kron(hu,hu),hu)*kron(kron(u,u),u)
%
%   Law of motion for control variables y (either VAROBS variables or if no VAROBS statement is given then for all endogenous variables)
%       y = steady_state(y)
%         + gx*( xf(-1)+xs(-1)+xrd(-1) )
%         + gu*u
%         + 1/2*gxx*kron(xf(-1),xf(-1)) + gxx*kron(xf(-1),xs(-1))
%         + gxu*kron(xf(-1),u) + gxu*kron(xs(-1),u)
%         + 1/2*guu*(kron(u,u)-Sigma_e+Sigma_e)
%         + 1/2*gss*sig^2
%         + 1/6*gxxx*kron(kron(xf(-1),xf(-1)),xf(-1))
%         + 1/6*guuu*kron(kron(u,u),u)
%         + 3/6*gxxu*kron(kron(xf(-1),xf(-1),u)
%         + 3/6*gxuu*kron(kron(xf(-1),u),u)
%         + 3/6*gxss*xf(-1)*sig^2
%         + 3/6*guss*u*sig^2
%
% See code below how z and inov are defined at first, second, and third order,
% and how to set up A, B, C, D and compute unconditional first and second moments of inov, z and y

persistent QPu COMBOS4 Q6Pu COMBOS6 K_u_xx K_u_ux K_xx_x

%% Auxiliary indices and objects
order = options.order;
if isempty(options.qz_criterium)
    % set default value for qz_criterium: if there are no unit roots one can use 1.0
    % If they are possible, you may have have multiple unit roots and the accuracy 
    % decreases when computing the eigenvalues in lyapunov_symm. Hence, we normally use 1+1e-6
    % Note that unit roots are only possible at first-order, at higher order we set it to 1
    options.qz_criterium = 1+1e-6;
end
indx = [M.nstatic+(1:M.nspred) M.endo_nbr+(1:size(dr.ghx,2)-M.nspred)]';
if isempty(indy)
    indy = (1:M.endo_nbr)'; %by default select all variables in DR order
end
u_nbr    = M.exo_nbr;
x_nbr    = length(indx);
y_nbr    = length(indy);
Yss      = dr.ys(dr.order_var);
hx       = dr.ghx(indx,:);
gx       = dr.ghx(indy,:);
hu       = dr.ghu(indx,:);
gu       = dr.ghu(indy,:);
E_uu     = M.Sigma_e; %this is E[u*u']

if compute_derivs
    stderrparam_nbr = length(dr.derivs.indpstderr);
    corrparam_nbr   = size(dr.derivs.indpcorr,1);
    modparam_nbr    = length(dr.derivs.indpmodel);
    totparam_nbr    = stderrparam_nbr+corrparam_nbr+modparam_nbr;
    dYss   = dr.derivs.dYss;
	dhx    = dr.derivs.dghx(indx,:,:);
    dgx    = dr.derivs.dghx(indy,:,:);
	dhu    = dr.derivs.dghu(indx,:,:);
    dgu    = dr.derivs.dghu(indy,:,:);
	dE_uu  = dr.derivs.dSigma_e;
end

% first-order approximation indices for extended state vector z and extended innovations vector inov 
id_z1_xf    = (1:x_nbr);
id_inov1_u  = (1:u_nbr);
if order > 1
    % second-order approximation indices for extended state vector z and extended innovations vector inov 
    id_z2_xs      = id_z1_xf(end)     + (1:x_nbr);
    id_z3_xf_xf   = id_z2_xs(end)     + (1:x_nbr*x_nbr);
    id_inov2_u_u  = id_inov1_u(end)   + (1:u_nbr*u_nbr);
    id_inov3_xf_u = id_inov2_u_u(end) + (1:x_nbr*u_nbr);

    hxx = dr.ghxx(indx,:);
    gxx = dr.ghxx(indy,:);
    hxu = dr.ghxu(indx,:);
    gxu = dr.ghxu(indy,:);
    huu = dr.ghuu(indx,:);
    guu = dr.ghuu(indy,:);
    hss = dr.ghs2(indx,:);
    gss = dr.ghs2(indy,:);
    if compute_derivs        
        dhxx = dr.derivs.dghxx(indx,:,:);
        dgxx = dr.derivs.dghxx(indy,:,:);
        dhxu = dr.derivs.dghxu(indx,:,:);
        dgxu = dr.derivs.dghxu(indy,:,:);
        dhuu = dr.derivs.dghuu(indx,:,:);
        dguu = dr.derivs.dghuu(indy,:,:);
        dhss = dr.derivs.dghs2(indx,:);
        dgss = dr.derivs.dghs2(indy,:);
    end
end
if order > 2
    % third-order approximation indices for extended state vector z and extended innovations vector inov 
    id_z4_xrd        = id_z3_xf_xf(end)      + (1:x_nbr);
    id_z5_xf_xs      = id_z4_xrd(end)        + (1:x_nbr*x_nbr);
    id_z6_xf_xf_xf   = id_z5_xf_xs(end)      + (1:x_nbr*x_nbr*x_nbr);
    id_inov4_xs_u    = id_inov3_xf_u(end)    + (1:x_nbr*u_nbr);
    id_inov5_xf_xf_u = id_inov4_xs_u(end)    + (1:x_nbr*x_nbr*u_nbr);
    id_inov6_xf_u_u  = id_inov5_xf_xf_u(end) + (1:x_nbr*u_nbr*u_nbr);
    id_inov7_u_u_u   = id_inov6_xf_u_u(end)  + (1:u_nbr*u_nbr*u_nbr);

    hxxx = dr.ghxxx(indx,:);
    gxxx = dr.ghxxx(indy,:);
    hxxu = dr.ghxxu(indx,:);
    gxxu = dr.ghxxu(indy,:);
    hxuu = dr.ghxuu(indx,:);
    gxuu = dr.ghxuu(indy,:);
    huuu = dr.ghuuu(indx,:);
    guuu = dr.ghuuu(indy,:);
    hxss = dr.ghxss(indx,:);
    gxss = dr.ghxss(indy,:);
    huss = dr.ghuss(indx,:);
    guss = dr.ghuss(indy,:);
    if compute_derivs
        dhxxx = dr.derivs.dghxxx(indx,:,:);
        dgxxx = dr.derivs.dghxxx(indy,:,:);
        dhxxu = dr.derivs.dghxxu(indx,:,:);
        dgxxu = dr.derivs.dghxxu(indy,:,:);
       	dhxuu = dr.derivs.dghxuu(indx,:,:);
        dgxuu = dr.derivs.dghxuu(indy,:,:);
        dhuuu = dr.derivs.dghuuu(indx,:,:);
        dguuu = dr.derivs.dghuuu(indy,:,:);
        dhxss = dr.derivs.dghxss(indx,:,:);
        dgxss = dr.derivs.dghxss(indy,:,:);
        dhuss = dr.derivs.dghuss(indx,:,:);
        dguss = dr.derivs.dghuss(indy,:,:);
    end
end

%% First-order state space system
% Auxiliary state vector z is defined by:          z    = [xf]
% Auxiliary innovations vector inov is defined by: inov = [u]
z_nbr       = x_nbr;
inov_nbr    = M.exo_nbr;
A           = hx;
B           = hu;
c           = zeros(x_nbr,1);
C           = gx;
D           = gu;
d           = zeros(y_nbr,1);
Varinov     = E_uu;
E_inovzlag1  = zeros(inov_nbr,z_nbr); %at first-order E[inov*z(-1)'] = 0
Om_z        = B*Varinov*B';
E_xf        = zeros(x_nbr,1);

lyapunov_symm_method = 1; %method=1 to initialize persistent variables
[Var_z,Schur_u] = lyapunov_symm(A, Om_z,... %at first-order this algorithm is well established and also used in th_autocovariances.m
                                options.lyapunov_fixed_point_tol, options.qz_criterium, options.lyapunov_complex_threshold,...
                                lyapunov_symm_method,...       
                                options.debug); %we use Schur_u to take care of (possible) nonstationary VAROBS variables in moment computations
%find stationary vars
stationary_vars = (1:y_nbr)';
if ~isempty(Schur_u)
    %base this only on first order, because if first-order is stable so are the higher-order pruned systems
    x = abs(gx*Schur_u);
    stationary_vars = find(all(x < options.schur_vec_tol,2));
end

if compute_derivs == 1
    dA          = dhx;
    dB          = dhu;
    dc          = zeros(x_nbr,totparam_nbr);
    dC          = dgx;
    dD          = dgu;
    dd          = zeros(y_nbr,totparam_nbr);
    dVarinov    = dE_uu;
    dE_xf       = zeros(x_nbr,totparam_nbr);
    dE_inovzlag1 = zeros(z_nbr,inov_nbr,totparam_nbr);
    dVar_z   = zeros(z_nbr,z_nbr,totparam_nbr);
    lyapunov_symm_method = 2;%to spare a lot of computing time while not repeating Schur every time
    for jp1 = 1:totparam_nbr
        if jp1 <= stderrparam_nbr+corrparam_nbr
            dOm_z_jp1 = B*dVarinov(:,:,jp1)*B';
            dVar_z(:,:,jp1) = lyapunov_symm(A, dOm_z_jp1,...
                                           options.lyapunov_fixed_point_tol,options.qz_criterium,options.lyapunov_complex_threshold,...
                                           lyapunov_symm_method,...
                                           options.debug);
        else
            dOm_z_jp1  = dB(:,:,jp1)*Varinov*B' + B*Varinov*dB(:,:,jp1)';
            dVar_z(:,:,jp1) = lyapunov_symm(A, dA(:,:,jp1)*Var_z*A' + A*Var_z*dA(:,:,jp1)' + dOm_z_jp1,...
                                           options.lyapunov_fixed_point_tol,options.qz_criterium,options.lyapunov_complex_threshold,...
                                           lyapunov_symm_method,...
                                           options.debug);
        end
    end
end

if order > 1
    options.qz_criterium = 1; %pruned state space has no unit roots
    % Some common and useful objects for order > 1
    E_xfxf   = Var_z;
    if compute_derivs
        dE_xfxf = dVar_z;
    end
    hx_hx    = kron(hx,hx);
    hx_hu    = kron(hx,hu);
    hu_hu    = kron(hu,hu);
    I_xx     = eye(x_nbr^2);
    K_x_x    = commutation(x_nbr,x_nbr,1);
    invIx_hx = (eye(x_nbr)-hx)\eye(x_nbr);

    %Compute unique fourth order product moments of u, i.e. unique(E[kron(kron(kron(u,u),u),u)],'stable')
    u_nbr4    = u_nbr*(u_nbr+1)/2*(u_nbr+2)/3*(u_nbr+3)/4;
    if isempty(QPu)
        QPu       = quadruplication(u_nbr);
        COMBOS4   = flipud(allVL1(u_nbr, 4)); %all possible (unique) combinations of powers that sum up to four
    end
    E_u_u_u_u = zeros(u_nbr4,1); %only unique entries
    if compute_derivs && (stderrparam_nbr+corrparam_nbr>0)
        dE_u_u_u_u = zeros(u_nbr4,stderrparam_nbr+corrparam_nbr);
    end
    for j4 = 1:size(COMBOS4,1)
        if compute_derivs && (stderrparam_nbr+corrparam_nbr>0)
            [E_u_u_u_u(j4), dE_u_u_u_u(j4,:)] = prodmom_deriv(E_uu, 1:u_nbr, COMBOS4(j4,:), dE_uu(:,:,1:(stderrparam_nbr+corrparam_nbr)), dr.derivs.dCorrelation_matrix(:,:,1:(stderrparam_nbr+corrparam_nbr)));
        else
            E_u_u_u_u(j4) = prodmom(E_uu, 1:u_nbr, COMBOS4(j4,:));
        end
    end
    E_xfxf_uu = kron(E_xfxf,E_uu');

%% Second-order pruned state space system
% Auxiliary state vector z is defined by:          z    = [xf;xs;kron(xf,xf)]
% Auxiliary innovations vector inov is defined by: inov = [u;kron(u,u)-E_uu(:);kron(xf,u)]
    z_nbr    = x_nbr + x_nbr + x_nbr^2;
    inov_nbr = u_nbr + u_nbr^2 + x_nbr*u_nbr;

    A = zeros(z_nbr, z_nbr);
    A(id_z1_xf    , id_z1_xf   ) = hx;
    A(id_z2_xs    , id_z2_xs   ) = hx;
    A(id_z2_xs    , id_z3_xf_xf) = 1/2*hxx;
    A(id_z3_xf_xf , id_z3_xf_xf) = hx_hx;

    B = zeros(z_nbr, inov_nbr);
    B(id_z1_xf    , id_inov1_u   ) = hu;
    B(id_z2_xs    , id_inov2_u_u ) = 1/2*huu;
    B(id_z2_xs    , id_inov3_xf_u) = hxu;
    B(id_z3_xf_xf , id_inov2_u_u ) = hu_hu;
    B(id_z3_xf_xf , id_inov3_xf_u) = (I_xx+K_x_x)*hx_hu;

    c = zeros(z_nbr, 1);
    c(id_z2_xs    , 1) = 1/2*hss + 1/2*huu*E_uu(:);
    c(id_z3_xf_xf , 1) = hu_hu*E_uu(:);

    C = zeros(y_nbr, z_nbr);
    C(: , id_z1_xf   ) = gx;
    C(: , id_z2_xs   ) = gx;
    C(: , id_z3_xf_xf) = 1/2*gxx;

    D = zeros(y_nbr, inov_nbr);
    D(: , id_inov1_u   ) = gu;
    D(: , id_inov2_u_u ) = 1/2*guu;
    D(: , id_inov3_xf_u) = gxu;

    d = 1/2*gss + 1/2*guu*E_uu(:);

    Varinov = zeros(inov_nbr,inov_nbr);
    Varinov(id_inov1_u    , id_inov1_u)    = E_uu;
   %Varinov(id_inov1_u    , id_inov2_u_u ) = zeros(u_nbr,u_nbr^2);
   %Varinov(id_inov1_u    , id_inov3_xf_u) = zeros(u_nbr,x_nbr*u_nbr);
   %Varinov(id_inov2_u_u  , id_inov1_u   ) = zeros(u_nbr^2,u_nbr);
    Varinov(id_inov2_u_u  , id_inov2_u_u ) = reshape(QPu*E_u_u_u_u,u_nbr^2,u_nbr^2)-E_uu(:)*E_uu(:)';
   %Varinov(id_inov2_u_u  , id_inov3_xf_u) = zeros(u_nbr^2,x_nbr*u_nbr);
   %Varinov(id_inov3_xf_u , id_inov1_u   ) = zeros(x_nbr*u_nbr,u_nbr);
   %Varinov(id_inov3_xf_u , id_inov2_u_u ) = zeros(x_nbr*u_nbr,u_nbr^2);
    Varinov(id_inov3_xf_u , id_inov3_xf_u) = E_xfxf_uu;

    E_xs        = invIx_hx*(1/2*hxx*E_xfxf(:) + c(id_z2_xs,1));
    E_inovzlag1 = zeros(inov_nbr,z_nbr); %at second-order E[z(-1)*inov'] = 0
    Om_z        = B*Varinov*transpose(B);

    lyapunov_symm_method = 1; %method=1 to initialize persistent variables (if errorflag)
    [Var_z, errorflag] = disclyap_fast(A,Om_z,options.lyapunov_doubling_tol);
    if errorflag %use Schur-based method
        fprintf('PRUNED_STATE_SPACE_SYSTEM: error flag in disclyap_fast at order=2, use lyapunov_symm\n');
        Var_z = lyapunov_symm(A,Om_z,...
                              options.lyapunov_fixed_point_tol,options.qz_criterium,options.lyapunov_complex_threshold,...
                              lyapunov_symm_method,...
                              options.debug);
        lyapunov_symm_method = 2; %in the following we can reuse persistent variables
    end
    % Make sure some stuff is zero due to Gaussianity
    Var_z(id_z1_xf    , id_z2_xs   ) = zeros(x_nbr,x_nbr);
    Var_z(id_z1_xf    , id_z3_xf_xf) = zeros(x_nbr,x_nbr^2);
    Var_z(id_z2_xs    , id_z1_xf   ) = zeros(x_nbr,x_nbr);
    Var_z(id_z3_xf_xf , id_z1_xf   ) = zeros(x_nbr^2,x_nbr);

    if compute_derivs
        dA           = zeros(z_nbr,z_nbr,totparam_nbr);
        dB           = zeros(z_nbr,inov_nbr,totparam_nbr);
        dc           = zeros(z_nbr,totparam_nbr);
        dC           = zeros(y_nbr,z_nbr,totparam_nbr);
        dD           = zeros(y_nbr,inov_nbr,totparam_nbr);
        dd           = zeros(y_nbr,totparam_nbr);
        dVarinov     = zeros(inov_nbr,inov_nbr,totparam_nbr);
        dE_xs        = zeros(x_nbr,totparam_nbr);
        dE_inovzlag1 = zeros(inov_nbr,z_nbr,totparam_nbr);
        dVar_z       = zeros(z_nbr,z_nbr,totparam_nbr);
        
        for jp2 = 1:totparam_nbr
            if jp2 <= (stderrparam_nbr+corrparam_nbr)
                dE_uu_jp2      = dE_uu(:,:,jp2);
                dE_u_u_u_u_jp2 = QPu*dE_u_u_u_u(:,jp2);                
            else
                dE_uu_jp2      = zeros(u_nbr,u_nbr);
                dE_u_u_u_u_jp2 = zeros(u_nbr^4,1);            
            end
            dhx_jp2        = dhx(:,:,jp2);
            dhu_jp2        = dhu(:,:,jp2);
            dhxx_jp2       = dhxx(:,:,jp2);
            dhxu_jp2       = dhxu(:,:,jp2);
            dhuu_jp2       = dhuu(:,:,jp2);
            dhss_jp2       = dhss(:,jp2);            
            dgx_jp2        = dgx(:,:,jp2);
            dgu_jp2        = dgu(:,:,jp2);
            dgxx_jp2       = dgxx(:,:,jp2);
            dgxu_jp2       = dgxu(:,:,jp2);
            dguu_jp2       = dguu(:,:,jp2);
            dgss_jp2       = dgss(:,jp2);
            dhx_hx_jp2     = kron(dhx_jp2,hx) + kron(hx,dhx_jp2);
            dhu_hu_jp2     = kron(dhu_jp2,hu) + kron(hu,dhu_jp2);
            dhx_hu_jp2     = kron(dhx_jp2,hu) + kron(hx,dhu_jp2);
            dE_xfxf_jp2    = dE_xfxf(:,:,jp2);
            dE_xfxf_uu_jp2 = kron(dE_xfxf_jp2,E_uu) + kron(E_xfxf,dE_uu_jp2);

            dA(id_z1_xf    , id_z1_xf    , jp2) = dhx_jp2;
            dA(id_z2_xs    , id_z2_xs    , jp2) = dhx_jp2;
            dA(id_z2_xs    , id_z3_xf_xf , jp2) = 1/2*dhxx_jp2;
            dA(id_z3_xf_xf , id_z3_xf_xf , jp2) = dhx_hx_jp2;

            dB(id_z1_xf    , id_inov1_u    , jp2) = dhu_jp2;
            dB(id_z2_xs    , id_inov2_u_u  , jp2) = 1/2*dhuu_jp2;
            dB(id_z2_xs    , id_inov3_xf_u , jp2) = dhxu_jp2;
            dB(id_z3_xf_xf , id_inov2_u_u  , jp2) = dhu_hu_jp2;
            dB(id_z3_xf_xf , id_inov3_xf_u , jp2) = (I_xx+K_x_x)*dhx_hu_jp2;

            dc(id_z2_xs    , jp2) = 1/2*dhss_jp2 + 1/2*dhuu_jp2*E_uu(:) + 1/2*huu*dE_uu_jp2(:);
            dc(id_z3_xf_xf , jp2) = dhu_hu_jp2*E_uu(:) + hu_hu*dE_uu_jp2(:);

            dC(: , id_z1_xf    , jp2) = dgx_jp2;
            dC(: , id_z2_xs    , jp2) = dgx_jp2;
            dC(: , id_z3_xf_xf , jp2) = 1/2*dgxx_jp2;

            dD(: , id_inov1_u    , jp2) = dgu_jp2;
            dD(: , id_inov2_u_u  , jp2) = 1/2*dguu_jp2;
            dD(: , id_inov3_xf_u , jp2) = dgxu_jp2;

            dd(:,jp2) = 1/2*dgss_jp2 + 1/2*guu*dE_uu_jp2(:) + 1/2*dguu_jp2*E_uu(:);

            dVarinov(id_inov1_u    , id_inov1_u    , jp2) = dE_uu_jp2;
            dVarinov(id_inov2_u_u  , id_inov2_u_u  , jp2) = reshape(dE_u_u_u_u_jp2,u_nbr^2,u_nbr^2) - dE_uu_jp2(:)*E_uu(:)' - E_uu(:)*dE_uu_jp2(:)';
            dVarinov(id_inov3_xf_u , id_inov3_xf_u , jp2) = dE_xfxf_uu_jp2;
            
            dE_xs(:,jp2) = invIx_hx*( dhx_jp2*E_xs + 1/2*dhxx_jp2*E_xfxf(:) + 1/2*hxx*dE_xfxf_jp2(:) + dc(id_z2_xs,jp2) );
            dOm_z_jp2    = dB(:,:,jp2)*Varinov*B' + B*dVarinov(:,:,jp2)*B' + B*Varinov*dB(:,:,jp2)';
            
            [dVar_z(:,:,jp2), errorflag] = disclyap_fast(A, dA(:,:,jp2)*Var_z*A' + A*Var_z*dA(:,:,jp2)' + dOm_z_jp2, options.lyapunov_doubling_tol);
            if errorflag
                dVar_z(:,:,jp2) = lyapunov_symm(A, dA(:,:,jp2)*Var_z*A' + A*Var_z*dA(:,:,jp2)' + dOm_z_jp2,...
                                               options.lyapunov_fixed_point_tol,options.qz_criterium,options.lyapunov_complex_threshold,...
                                               lyapunov_symm_method,...
                                               options.debug);
                if lyapunov_symm_method == 1
                    lyapunov_symm_method = 2; %now we can reuse persistent schur
                end
            end
            % Make sure some stuff is zero due to Gaussianity
            dVar_z(id_z1_xf    , id_z2_xs    , jp2) = zeros(x_nbr,x_nbr);
            dVar_z(id_z1_xf    , id_z3_xf_xf , jp2) = zeros(x_nbr,x_nbr^2);    
            dVar_z(id_z2_xs    , id_z1_xf    , jp2) = zeros(x_nbr,x_nbr);    
            dVar_z(id_z3_xf_xf , id_z1_xf    , jp2) = zeros(x_nbr^2,x_nbr);
        end
    end

    if order > 2
        % Some common and useful objects for order > 2
        if isempty(K_u_xx)
            K_u_xx   = commutation(u_nbr,x_nbr^2,1);
            K_u_ux   = commutation(u_nbr,u_nbr*x_nbr,1);
            K_xx_x   = commutation(x_nbr^2,x_nbr);
        end
        hx_hss2  = kron(hx,1/2*hss);
        hu_hss2  = kron(hu,1/2*hss);
        hx_hxx2  = kron(hx,1/2*hxx);
        hxx2_hu  = kron(1/2*hxx,hu);
        hx_hxu   = kron(hx,hxu);
        hxu_hu   = kron(hxu,hu);
        hx_huu2  = kron(hx,1/2*huu);
        hu_huu2  = kron(hu,1/2*huu);
        hx_hx_hx = kron(hx,hx_hx);
        hx_hx_hu = kron(hx_hx,hu);
        hu_hx_hx = kron(hu,hx_hx);
        hu_hu_hu = kron(hu_hu,hu);
        hx_hu_hu = kron(hx,hu_hu);
        hu_hx_hu = kron(hu,hx_hu);
        invIxx_hx_hx = (eye(x_nbr^2)-hx_hx)\eye(x_nbr^2);

        % Reuse second-order results
       %E_xfxf       = Var_z(id_z1_xf, id_z1_xf      );                        %this is E[xf*xf'], we already have that
       %E_xfxs       = Var_z(id_z1_xf, id_z2_xs      );                        %this is E[xf*xs']=0 due to gaussianity
       %E_xfxf_xf    = Var_z(id_z1_xf, id_z3_xf_xf   );                        %this is E[xf*kron(xf_xf)']=0 due to gaussianity
       %E_xsxf       = Var_z(id_z2_xs, id_z1_xf      );                        %this is E[xs*xf']=0 due to gaussianity
        E_xsxs       = Var_z(id_z2_xs, id_z2_xs      ) + E_xs*transpose(E_xs); %this is E[xs*xs']
        E_xsxf_xf    = Var_z(id_z2_xs, id_z3_xf_xf   ) + E_xs*E_xfxf(:)';      %this is E[xs*kron(xf,xf)']
       %E_xf_xfxf    = Var_z(id_z3_xf_xf, id_z1_xf   );                        %this is E[kron(xf,xf)*xf']=0 due to gaussianity
        E_xf_xfxs    = Var_z(id_z3_xf_xf, id_z2_xs   ) + E_xfxf(:)*E_xs';      %this is E[kron(xf,xf)*xs']
        E_xf_xfxf_xf = Var_z(id_z3_xf_xf, id_z3_xf_xf) + E_xfxf(:)*E_xfxf(:)'; %this is E[kron(xf,xf)*kron(xf,xf)']
        E_xrdxf = reshape(invIxx_hx_hx*vec(...
                                             hxx*reshape( commutation(x_nbr^2,x_nbr,1)*E_xf_xfxs(:), x_nbr^2,x_nbr)*hx'...
                                           + hxu*kron(E_xs,E_uu)*hu'...
                                           + 1/6*hxxx*reshape(E_xf_xfxf_xf,x_nbr^3,x_nbr)*hx'...
                                           + 1/6*huuu*reshape(QPu*E_u_u_u_u,u_nbr^3,u_nbr)*hu'...
                                           + 3/6*hxxu*kron(E_xfxf(:),E_uu)*hu'...
                                           + 3/6*hxuu*kron(E_xfxf,E_uu(:))*hx'...
                                           + 3/6*hxss*E_xfxf*hx'...
                                           + 3/6*huss*E_uu*hu'...
                                           ),...
                         x_nbr,x_nbr); %this is E[xrd*xf']
        if compute_derivs
            dE_xsxs       = zeros(x_nbr,x_nbr,totparam_nbr);
            dE_xsxf_xf    = zeros(x_nbr,x_nbr^2,totparam_nbr);
            dE_xf_xfxs    = zeros(x_nbr^2,x_nbr,totparam_nbr);
            dE_xf_xfxf_xf = zeros(x_nbr^2,x_nbr^2,totparam_nbr);
            dE_xrdxf      = zeros(x_nbr,x_nbr,totparam_nbr);
            for jp2 = 1:totparam_nbr
                if jp2 < (stderrparam_nbr+corrparam_nbr)
                    dE_u_u_u_u_jp2 = QPu*dE_u_u_u_u(:,jp2);
                else
                    dE_u_u_u_u_jp2 = zeros(u_nbr^4,1);
                end
                dE_xsxs(:,:,jp2)       = dVar_z(id_z2_xs    , id_z2_xs    , jp2) + dE_xs(:,jp2)*transpose(E_xs) + E_xs*transpose(dE_xs(:,jp2));
                dE_xsxf_xf(:,:,jp2)    = dVar_z(id_z2_xs    , id_z3_xf_xf , jp2) + dE_xs(:,jp2)*E_xfxf(:)' + E_xs*vec(dE_xfxf(:,:,jp2))';
                dE_xf_xfxs(:,:,jp2)    = dVar_z(id_z3_xf_xf , id_z2_xs    , jp2) + vec(dE_xfxf(:,:,jp2))*E_xs' + E_xfxf(:)*dE_xs(:,jp2)';
                dE_xf_xfxf_xf(:,:,jp2) = dVar_z(id_z3_xf_xf , id_z3_xf_xf , jp2) + vec(dE_xfxf(:,:,jp2))*E_xfxf(:)' + E_xfxf(:)*vec(dE_xfxf(:,:,jp2))';
                dE_xrdxf(:,:,jp2) = reshape(invIxx_hx_hx*vec(...
                    dhx(:,:,jp2)*E_xrdxf*hx' + hx*E_xrdxf*dhx(:,:,jp2)'...
                  + dhxx(:,:,jp2)*reshape( commutation(x_nbr^2,x_nbr,1)*E_xf_xfxs(:), x_nbr^2,x_nbr)*hx' + hxx*reshape( commutation(x_nbr^2,x_nbr,1)*vec(dE_xf_xfxs(:,:,jp2)), x_nbr^2,x_nbr)*hx' + hxx*reshape( commutation(x_nbr^2,x_nbr,1)*E_xf_xfxs(:), x_nbr^2,x_nbr)*dhx(:,:,jp2)'...
                  + dhxu(:,:,jp2)*kron(E_xs,E_uu)*hu' + hxu*kron(dE_xs(:,jp2),E_uu)*hu' + hxu*kron(E_xs,dE_uu(:,:,jp2))*hu' + hxu*kron(E_xs,E_uu)*dhu(:,:,jp2)'...
                  + 1/6*dhxxx(:,:,jp2)*reshape(E_xf_xfxf_xf,x_nbr^3,x_nbr)*hx' + 1/6*hxxx*reshape(dE_xf_xfxf_xf(:,:,jp2),x_nbr^3,x_nbr)*hx' + 1/6*hxxx*reshape(E_xf_xfxf_xf,x_nbr^3,x_nbr)*dhx(:,:,jp2)'...
                  + 1/6*dhuuu(:,:,jp2)*reshape(QPu*E_u_u_u_u,u_nbr^3,u_nbr)*hu' + 1/6*huuu*reshape(dE_u_u_u_u_jp2,u_nbr^3,u_nbr)*hu' + 1/6*huuu*reshape(QPu*E_u_u_u_u,u_nbr^3,u_nbr)*dhu(:,:,jp2)'...
                  + 3/6*dhxxu(:,:,jp2)*kron(E_xfxf(:),E_uu)*hu' + 3/6*hxxu*kron(vec(dE_xfxf(:,:,jp2)),E_uu)*hu' + 3/6*hxxu*kron(E_xfxf(:),dE_uu(:,:,jp2))*hu' + 3/6*hxxu*kron(E_xfxf(:),E_uu)*dhu(:,:,jp2)'...
                  + 3/6*dhxuu(:,:,jp2)*kron(E_xfxf,E_uu(:))*hx' + 3/6*hxuu*kron(dE_xfxf(:,:,jp2),E_uu(:))*hx' + 3/6*hxuu*kron(E_xfxf,vec(dE_uu(:,:,jp2)))*hx' + 3/6*hxuu*kron(E_xfxf,E_uu(:))*dhx(:,:,jp2)'...
                  + 3/6*dhxss(:,:,jp2)*E_xfxf*hx' + 3/6*hxss*dE_xfxf(:,:,jp2)*hx' + 3/6*hxss*E_xfxf*dhx(:,:,jp2)'...
                  + 3/6*dhuss(:,:,jp2)*E_uu*hu' + 3/6*huss*dE_uu(:,:,jp2)*hu' + 3/6*huss*E_uu*dhu(:,:,jp2)'...
                  ), x_nbr, x_nbr);
            end
        end

        % Compute unique sixth-order product moments of u, i.e. unique(E[kron(kron(kron(kron(kron(u,u),u),u),u),u)],'stable')
        u_nbr6        = u_nbr*(u_nbr+1)/2*(u_nbr+2)/3*(u_nbr+3)/4*(u_nbr+4)/5*(u_nbr+5)/6;       
        if isempty(Q6Pu)
            Q6Pu          = Q6_plication(u_nbr);
            COMBOS6       = flipud(allVL1(u_nbr, 6)); %all possible (unique) combinations of powers that sum up to six
        end
        E_u_u_u_u_u_u = zeros(u_nbr6,1); %only unique entries
        if compute_derivs && (stderrparam_nbr+corrparam_nbr>0)
            dE_u_u_u_u_u_u = zeros(u_nbr6,stderrparam_nbr+corrparam_nbr);
        end
        for j6 = 1:size(COMBOS6,1)
            if compute_derivs && (stderrparam_nbr+corrparam_nbr>0)
                [E_u_u_u_u_u_u(j6), dE_u_u_u_u_u_u(j6,:)] = prodmom_deriv(E_uu, 1:u_nbr, COMBOS6(j6,:), dE_uu(:,:,1:(stderrparam_nbr+corrparam_nbr)), dr.derivs.dCorrelation_matrix(:,:,1:(stderrparam_nbr+corrparam_nbr)));
            else
                E_u_u_u_u_u_u(j6) = prodmom(E_uu, 1:u_nbr, COMBOS6(j6,:));
            end
        end

%% Third-order pruned state space system
% Auxiliary state vector z is defined by:          z    = [xf; xs; kron(xf,xf); xrd; kron(xf,xs); kron(kron(xf,xf),xf)]
% Auxiliary innovations vector inov is defined by: inov = [u; kron(u,u)-E_uu(:); kron(xf,u); kron(xs,u); kron(kron(xf,xf),u); kron(kron(xf,u),u); kron(kron(u,u),u))]
        z_nbr    = x_nbr + x_nbr + x_nbr^2 + x_nbr + x_nbr^2 + x_nbr^3;
        inov_nbr = u_nbr + u_nbr^2 + x_nbr*u_nbr + x_nbr*u_nbr + x_nbr^2*u_nbr + x_nbr*u_nbr^2 + u_nbr^3;

        A = zeros(z_nbr,z_nbr);
        A(id_z1_xf       , id_z1_xf      ) = hx;
        A(id_z2_xs       , id_z2_xs      ) = hx;
        A(id_z2_xs       , id_z3_xf_xf   ) = 1/2*hxx;
        A(id_z3_xf_xf    , id_z3_xf_xf   ) = hx_hx;
        A(id_z4_xrd      , id_z1_xf      ) = 3/6*hxss;
        A(id_z4_xrd      , id_z4_xrd     ) = hx;
        A(id_z4_xrd      , id_z5_xf_xs   ) = hxx;
        A(id_z4_xrd      , id_z6_xf_xf_xf) = 1/6*hxxx;
        A(id_z5_xf_xs    , id_z1_xf      ) = hx_hss2;
        A(id_z5_xf_xs    , id_z5_xf_xs   ) = hx_hx;
        A(id_z5_xf_xs    , id_z6_xf_xf_xf) = hx_hxx2;
        A(id_z6_xf_xf_xf , id_z6_xf_xf_xf) = hx_hx_hx;

        B = zeros(z_nbr,inov_nbr);
        B(id_z1_xf       , id_inov1_u      ) = hu;
        B(id_z2_xs       , id_inov2_u_u    ) = 1/2*huu;
        B(id_z2_xs       , id_inov3_xf_u   ) = hxu;
        B(id_z3_xf_xf    , id_inov2_u_u    ) = hu_hu;
        B(id_z3_xf_xf    , id_inov3_xf_u   ) = (I_xx+K_x_x)*hx_hu;
        B(id_z4_xrd      , id_inov1_u      ) = 3/6*huss;
        B(id_z4_xrd      , id_inov4_xs_u   ) = hxu;
        B(id_z4_xrd      , id_inov5_xf_xf_u) = 3/6*hxxu;
        B(id_z4_xrd      , id_inov6_xf_u_u ) = 3/6*hxuu;
        B(id_z4_xrd      , id_inov7_u_u_u  ) = 1/6*huuu;
        B(id_z5_xf_xs    , id_inov1_u      ) = hu_hss2;
        B(id_z5_xf_xs    , id_inov4_xs_u   ) = K_x_x*hx_hu;
        B(id_z5_xf_xs    , id_inov5_xf_xf_u) = hx_hxu + K_x_x*hxx2_hu;
        B(id_z5_xf_xs    , id_inov6_xf_u_u ) = hx_huu2 + K_x_x*hxu_hu;
        B(id_z5_xf_xs    , id_inov7_u_u_u  ) = hu_huu2;
        B(id_z6_xf_xf_xf , id_inov5_xf_xf_u) = hx_hx_hu + kron(hx,K_x_x*hx_hu) + hu_hx_hx*K_u_xx;
        B(id_z6_xf_xf_xf , id_inov6_xf_u_u ) = hx_hu_hu + hu_hx_hu*K_u_ux + kron(hu,K_x_x*hx_hu)*K_u_ux;
        B(id_z6_xf_xf_xf , id_inov7_u_u_u  ) = hu_hu_hu;

        c = zeros(z_nbr, 1);
        c(id_z2_xs    , 1) = 1/2*hss + 1/2*huu*E_uu(:);
        c(id_z3_xf_xf , 1) = hu_hu*E_uu(:);

        C = zeros(y_nbr,z_nbr);
        C(: , id_z1_xf      ) = gx + 3/6*gxss;
        C(: , id_z2_xs      ) = gx;
        C(: , id_z3_xf_xf   ) = 1/2*gxx;
        C(: , id_z4_xrd     ) = gx;
        C(: , id_z5_xf_xs   ) = gxx;
        C(: , id_z6_xf_xf_xf) = 1/6*gxxx;

        D = zeros(y_nbr,inov_nbr);
        D(: , id_inov1_u      ) = gu + 3/6*guss;
        D(: , id_inov2_u_u    ) = 1/2*guu;
        D(: , id_inov3_xf_u   ) = gxu;
        D(: , id_inov4_xs_u   ) = gxu;
        D(: , id_inov5_xf_xf_u) = 3/6*gxxu;
        D(: , id_inov6_xf_u_u)  = 3/6*gxuu;
        D(: , id_inov7_u_u_u )  = 1/6*guuu;

        d = 1/2*gss + 1/2*guu*E_uu(:);

        Varinov = zeros(inov_nbr,inov_nbr);
        Varinov(id_inov1_u       , id_inov1_u      ) = E_uu;
       %Varinov(id_inov1_u       , id_inov2_u_u    ) = zeros(u_nbr,u_nbr^2);
       %Varinov(id_inov1_u       , id_inov3_xf_u   ) = zeros(u_nbr,x_nbr*u_nbr);
        Varinov(id_inov1_u       , id_inov4_xs_u   ) = kron(E_xs',E_uu);
        Varinov(id_inov1_u       , id_inov5_xf_xf_u) = kron(E_xfxf(:)',E_uu);
       %Varinov(id_inov1_u       , id_inov6_xf_u_u ) = zeros(u_nbr,x_nbr*u_nbr^2);
        Varinov(id_inov1_u       , id_inov7_u_u_u  ) = reshape(QPu*E_u_u_u_u,u_nbr,u_nbr^3);

       %Varinov(id_inov2_u_u     , id_inov1_u      ) = zeros(u_nbr^2,u_nbr);
        Varinov(id_inov2_u_u     , id_inov2_u_u    ) = reshape(QPu*E_u_u_u_u,u_nbr^2,u_nbr^2)-E_uu(:)*E_uu(:)';
       %Varinov(id_inov2_u_u     , id_inov3_xf_u   ) = zeros(u_nbr^2,x_nbr*u_nbr);
       %Varinov(id_inov2_u_u     , id_inov4_xs_u   ) = zeros(u_nbr^2,x_nbr*u_nbr);
       %Varinov(id_inov2_u_u     , id_inov5_xf_xf_u) = zeros(u_nbr^2,x_nbr^2,u_nbr);
       %Varinov(id_inov2_u_u     , id_inov6_xf_u_u ) = zeros(u_nbr^2,x_nbr*u_nbr^2);
       %Varinov(id_inov2_u_u     , id_inov7_u_u_u  ) = zeros(u_nbr^2,u_nbr^3);

       %Varinov(id_inov3_xf_u    , id_inov1_u      ) = zeros(x_nbr*u_nbr,u_nbr);
       %Varinov(id_inov3_xf_u    , id_inov2_u_u    ) = zeros(x_nbr*u_nbr,u_nbr^2);
        Varinov(id_inov3_xf_u    , id_inov3_xf_u   ) = E_xfxf_uu;
       %Varinov(id_inov3_xf_u    , id_inov4_xs_u   ) = zeros(x_nbr*u_nbr,x_nbr*u_nbr);
       %Varinov(id_inov3_xf_u    , id_inov5_xf_xf_u) = zeros(x_nbr*u_nbr,x_nbr^2*u_nbr);
       %Varinov(id_inov3_xf_u    , id_inov6_xf_u_u ) = zeros(x_nbr*u_nbr,x_nbr*u_nbr^2);
       %Varinov(id_inov3_xf_u    , id_inov7_u_u_u   ) = zeros(x_nbr*u_nbr,u_nbr^3);

        Varinov(id_inov4_xs_u    , id_inov1_u      ) = kron(E_xs,E_uu);
       %Varinov(id_inov4_xs_u    , id_inov2_u_u    ) = zeros(x_nbr*u_nbr,u_nbr^2);
       %Varinov(id_inov4_xs_u    , id_inov3_xf_u   ) = zeros(x_nbr*u_nbr,x_nbr*u_nbr);
        Varinov(id_inov4_xs_u    , id_inov4_xs_u   ) = kron(E_xsxs,E_uu);
        Varinov(id_inov4_xs_u    , id_inov5_xf_xf_u) = kron(E_xsxf_xf, E_uu);
       %Varinov(id_inov4_xs_u    , id_inov6_xf_u_u ) = zeros(x_nbr*u_nbr,x_nbr*u_nbr^2);
        Varinov(id_inov4_xs_u    , id_inov7_u_u_u  ) = kron(E_xs,reshape(QPu*E_u_u_u_u,u_nbr,u_nbr^3));

        Varinov(id_inov5_xf_xf_u , id_inov1_u      ) = kron(E_xfxf(:),E_uu);
       %Varinov(id_inov5_xf_xf_u , id_inov2_u_u    ) = zeros(x_nbr^2*u_nbr,u_nbr^2);
       %Varinov(id_inov5_xf_xf_u , id_inov3_xf_u   ) = zeros(x_nbr^2*u_nbr,x_nbr*u_nbr);
        Varinov(id_inov5_xf_xf_u , id_inov4_xs_u   ) = kron(E_xf_xfxs,E_uu);
        Varinov(id_inov5_xf_xf_u , id_inov5_xf_xf_u) = kron(E_xf_xfxf_xf,E_uu);
       %Varinov(id_inov5_xf_xf_u , id_inov6_xf_u_u ) = zeros(x_nbr^2*u_nbr,x_nbr*u_nbr^2);
        Varinov(id_inov5_xf_xf_u , id_inov7_u_u_u  ) = kron(E_xfxf(:),reshape(QPu*E_u_u_u_u,u_nbr,u_nbr^3));

       %Varinov(id_inov6_xf_u_u  , id_inov1_u      ) = zeros(x_nbr*u_nbr^2,u_nbr);
       %Varinov(id_inov6_xf_u_u  , id_inov2_u_u    ) = zeros(x_nbr*u_nbr^2,u_nbr^2);
       %Varinov(id_inov6_xf_u_u  , id_inov3_xf_u   ) = zeros(x_nbr*u_nbr^2,x_nbr*u_nbr);
       %Varinov(id_inov6_xf_u_u  , id_inov4_xs_u   ) = zeros(x_nbr*u_nbr^2,x_nbr*u_nbr);
       %Varinov(id_inov6_xf_u_u  , id_inov5_xf_xf_u) = zeros(x_nbr*u_nbr^2,x_nbr^2*u_nbr);
        Varinov(id_inov6_xf_u_u  , id_inov6_xf_u_u ) = kron(E_xfxf,reshape(QPu*E_u_u_u_u,u_nbr^2,u_nbr^2));
       %Varinov(id_inov6_xf_u_u  , id_inov7_u_u_u  ) = zeros(x_nbr*u_nbr^2,u_nbr^3);

        Varinov(id_inov7_u_u_u   , id_inov1_u      ) = reshape(QPu*E_u_u_u_u,u_nbr^3,u_nbr);
       %Varinov(id_inov7_u_u_u   , id_inov2_u_u    ) = zeros(u_nbr^3,u_nbr^2);
       %Varinov(id_inov7_u_u_u   , id_inov3_xf_u   ) = zeros(u_nbr^3,x_nbr*u_nbr);
        Varinov(id_inov7_u_u_u   , id_inov4_xs_u   ) = kron(E_xs',reshape(QPu*E_u_u_u_u,u_nbr^3,u_nbr));
        Varinov(id_inov7_u_u_u   , id_inov5_xf_xf_u) = kron(transpose(E_xfxf(:)),reshape(QPu*E_u_u_u_u,u_nbr^3,u_nbr));
       %Varinov(id_inov7_u_u_u   , id_inov6_xf_u_u ) = zeros(u_nbr^3,x_nbr*u_nbr^2);
        Varinov(id_inov7_u_u_u   , id_inov7_u_u_u  ) = reshape(Q6Pu*E_u_u_u_u_u_u,u_nbr^3,u_nbr^3);

        E_xrd = zeros(x_nbr,1);%due to gaussianity

        E_inovzlag1 = zeros(inov_nbr,z_nbr); % Attention: E[inov*z(-1)'] is not equal to zero for a third-order approximation due to kron(kron(xf(-1),u),u)
        E_inovzlag1(id_inov6_xf_u_u , id_z1_xf       ) = kron(E_xfxf,E_uu(:));
        E_inovzlag1(id_inov6_xf_u_u , id_z4_xrd      ) = kron(E_xrdxf',E_uu(:));
        E_inovzlag1(id_inov6_xf_u_u , id_z5_xf_xs    ) = kron(reshape(K_xx_x*vec(E_xsxf_xf),x_nbr,x_nbr^2),vec(E_uu)) ;
        E_inovzlag1(id_inov6_xf_u_u , id_z6_xf_xf_xf ) = kron(reshape(E_xf_xfxf_xf,x_nbr,x_nbr^3),E_uu(:));

        Binovzlag1A= B*E_inovzlag1*transpose(A);
        Om_z = B*Varinov*transpose(B) + Binovzlag1A + transpose(Binovzlag1A);

        lyapunov_symm_method = 1; %method=1 to initialize persistent variables
        [Var_z, errorflag] = disclyap_fast(A,Om_z,options.lyapunov_doubling_tol);        
        if errorflag %use Schur-based method
            fprintf('PRUNED_STATE_SPACE_SYSTEM: error flag in disclyap_fast at order=3, use lyapunov_symm\n');
            Var_z = lyapunov_symm(A,Om_z,...
                                    options.lyapunov_fixed_point_tol,options.qz_criterium,options.lyapunov_complex_threshold,...
                                    lyapunov_symm_method,...
                                    options.debug);
            lyapunov_symm_method = 2; %we can now make use of persistent variables from shur
        end
        %make sure some stuff is zero due to Gaussianity
        Var_z(id_z1_xf       , id_z2_xs)       = zeros(x_nbr,x_nbr);
        Var_z(id_z1_xf       , id_z3_xf_xf)    = zeros(x_nbr,x_nbr^2);
        Var_z(id_z2_xs       , id_z1_xf)       = zeros(x_nbr,x_nbr);
        Var_z(id_z2_xs       , id_z4_xrd)      = zeros(x_nbr,x_nbr);
        Var_z(id_z2_xs       , id_z5_xf_xs)    = zeros(x_nbr,x_nbr^2);
        Var_z(id_z2_xs       , id_z6_xf_xf_xf) = zeros(x_nbr,x_nbr^3);
        Var_z(id_z3_xf_xf    , id_z1_xf)       = zeros(x_nbr^2,x_nbr);
        Var_z(id_z3_xf_xf    , id_z4_xrd)      = zeros(x_nbr^2,x_nbr);
        Var_z(id_z3_xf_xf    , id_z5_xf_xs)    = zeros(x_nbr^2,x_nbr^2);
        Var_z(id_z3_xf_xf    , id_z6_xf_xf_xf) = zeros(x_nbr^2,x_nbr^3);
        Var_z(id_z4_xrd      , id_z2_xs)       = zeros(x_nbr,x_nbr);
        Var_z(id_z4_xrd      , id_z3_xf_xf)    = zeros(x_nbr,x_nbr^2);
        Var_z(id_z5_xf_xs    , id_z2_xs)       = zeros(x_nbr^2,x_nbr);
        Var_z(id_z5_xf_xs    , id_z3_xf_xf)    = zeros(x_nbr^2,x_nbr^2);
        Var_z(id_z6_xf_xf_xf , id_z2_xs)       = zeros(x_nbr^3,x_nbr);
        Var_z(id_z6_xf_xf_xf , id_z3_xf_xf)    = zeros(x_nbr^3,x_nbr^2);
        
        if compute_derivs
            dA           = zeros(z_nbr,z_nbr,totparam_nbr);
            dB           = zeros(z_nbr,inov_nbr,totparam_nbr);
            dc           = zeros(z_nbr,totparam_nbr);
            dC           = zeros(y_nbr,z_nbr,totparam_nbr);
            dD           = zeros(y_nbr,inov_nbr,totparam_nbr);
            dd           = zeros(y_nbr,totparam_nbr);
            dVarinov     = zeros(inov_nbr,inov_nbr,totparam_nbr);
            dE_xrd       = zeros(x_nbr,totparam_nbr);
            dE_inovzlag1 = zeros(inov_nbr,z_nbr,totparam_nbr);
            dVar_z       = zeros(z_nbr,z_nbr,totparam_nbr);
            
            for jp3 = 1:totparam_nbr
                if jp3 <= (stderrparam_nbr+corrparam_nbr)
                    dE_uu_jp3          = dE_uu(:,:,jp3);
                    dE_u_u_u_u_jp3     = QPu*dE_u_u_u_u(:,jp3);
                    dE_u_u_u_u_u_u_jp3 = Q6Pu*dE_u_u_u_u_u_u(:,jp3);
                else
                    dE_uu_jp3          = zeros(u_nbr,u_nbr);
                    dE_u_u_u_u_jp3     = zeros(u_nbr^4,1);
                    dE_u_u_u_u_u_u_jp3 = zeros(u_nbr^6,1);
                end
                dhx_jp3       = dhx(:,:,jp3);
                dhu_jp3       = dhu(:,:,jp3);
                dhxx_jp3      = dhxx(:,:,jp3);
                dhxu_jp3      = dhxu(:,:,jp3);
                dhuu_jp3      = dhuu(:,:,jp3);
                dhss_jp3      = dhss(:,jp3);
                dhxxx_jp3     = dhxxx(:,:,jp3);
                dhxxu_jp3     = dhxxu(:,:,jp3);
                dhxuu_jp3     = dhxuu(:,:,jp3);
                dhuuu_jp3     = dhuuu(:,:,jp3);
                dhxss_jp3     = dhxss(:,:,jp3);
                dhuss_jp3     = dhuss(:,:,jp3);
                dgx_jp3       = dgx(:,:,jp3);
                dgu_jp3       = dgu(:,:,jp3);
                dgxx_jp3      = dgxx(:,:,jp3);
                dgxu_jp3      = dgxu(:,:,jp3);
                dguu_jp3      = dguu(:,:,jp3);
                dgss_jp3      = dgss(:,jp3);
                dgxxx_jp3     = dgxxx(:,:,jp3);
                dgxxu_jp3     = dgxxu(:,:,jp3);
                dgxuu_jp3     = dgxuu(:,:,jp3);
                dguuu_jp3     = dguuu(:,:,jp3);
                dgxss_jp3     = dgxss(:,:,jp3);
                dguss_jp3     = dguss(:,:,jp3);
                
                dhx_hx_jp3    = kron(dhx_jp3,hx) + kron(hx,dhx_jp3);
                dhx_hu_jp3    = kron(dhx_jp3,hu) + kron(hx,dhu_jp3);                
                dhu_hu_jp3    = kron(dhu_jp3,hu) + kron(hu,dhu_jp3);
                dhx_hss2_jp3  = kron(dhx_jp3,1/2*hss) + kron(hx,1/2*dhss_jp3);
                dhu_hss2_jp3  = kron(dhu_jp3,1/2*hss) + kron(hu,1/2*dhss_jp3);
                dhx_hxx2_jp3  = kron(dhx_jp3,1/2*hxx) + kron(hx,1/2*dhxx_jp3);
                dhxx2_hu_jp3  = kron(1/2*dhxx_jp3,hu) + kron(1/2*hxx,dhu_jp3);
                dhx_hxu_jp3   = kron(dhx_jp3,hxu) + kron(hx,dhxu_jp3);
                dhxu_hu_jp3   = kron(dhxu_jp3,hu) + kron(hxu,dhu_jp3);
                dhx_huu2_jp3  = kron(dhx_jp3,1/2*huu) + kron(hx,1/2*dhuu_jp3);
                dhu_huu2_jp3  = kron(dhu_jp3,1/2*huu) + kron(hu,1/2*dhuu_jp3);
                dhx_hx_hx_jp3 = kron(dhx_jp3,hx_hx) + kron(hx,dhx_hx_jp3);
                dhx_hx_hu_jp3 = kron(dhx_hx_jp3,hu) + kron(hx_hx,dhu_jp3);
                dhu_hx_hx_jp3 = kron(dhu_jp3,hx_hx) + kron(hu,dhx_hx_jp3);
                dhu_hu_hu_jp3 = kron(dhu_hu_jp3,hu) + kron(hu_hu,dhu_jp3);
                dhx_hu_hu_jp3 = kron(dhx_jp3,hu_hu) + kron(hx,dhu_hu_jp3);
                dhu_hx_hu_jp3 = kron(dhu_jp3,hx_hu) + kron(hu,dhx_hu_jp3);

                dE_xs_jp3         = dE_xs(:,jp3);
                dE_xfxf_jp3       = dE_xfxf(:,:,jp3);
                dE_xsxs_jp3       = dE_xsxs(:,:,jp3);
                dE_xsxf_xf_jp3    = dE_xsxf_xf(:,:,jp3);
                dE_xfxf_uu_jp3    = kron(dE_xfxf_jp3,E_uu) + kron(E_xfxf,dE_uu_jp3);
                dE_xf_xfxs_jp3    = dE_xf_xfxs(:,:,jp3);
                dE_xf_xfxf_xf_jp3 = dE_xf_xfxf_xf(:,:,jp3);
                dE_xrdxf_jp3      = dE_xrdxf(:,:,jp3);
                
                dA(id_z1_xf       , id_z1_xf       , jp3) = dhx_jp3;
                dA(id_z2_xs       , id_z2_xs       , jp3) = dhx_jp3;
                dA(id_z2_xs       , id_z3_xf_xf    , jp3) = 1/2*dhxx_jp3;
                dA(id_z3_xf_xf    , id_z3_xf_xf    , jp3) = dhx_hx_jp3;
                dA(id_z4_xrd      , id_z1_xf       , jp3) = 3/6*dhxss_jp3;
                dA(id_z4_xrd      , id_z4_xrd      , jp3) = dhx_jp3;
                dA(id_z4_xrd      , id_z5_xf_xs    , jp3) = dhxx_jp3;
                dA(id_z4_xrd      , id_z6_xf_xf_xf , jp3) = 1/6*dhxxx_jp3;
                dA(id_z5_xf_xs    , id_z1_xf       , jp3) = dhx_hss2_jp3;
                dA(id_z5_xf_xs    , id_z5_xf_xs    , jp3) = dhx_hx_jp3;
                dA(id_z5_xf_xs    , id_z6_xf_xf_xf , jp3) = dhx_hxx2_jp3;
                dA(id_z6_xf_xf_xf , id_z6_xf_xf_xf , jp3) = dhx_hx_hx_jp3;

                dB(id_z1_xf       , id_inov1_u       , jp3) = dhu_jp3;
                dB(id_z2_xs       , id_inov2_u_u     , jp3) = 1/2*dhuu_jp3;
                dB(id_z2_xs       , id_inov3_xf_u    , jp3) = dhxu_jp3;
                dB(id_z3_xf_xf    , id_inov2_u_u     , jp3) = dhu_hu_jp3;
                dB(id_z3_xf_xf    , id_inov3_xf_u    , jp3) = (I_xx+K_x_x)*dhx_hu_jp3;
                dB(id_z4_xrd      , id_inov1_u       , jp3) = 3/6*dhuss_jp3;
                dB(id_z4_xrd      , id_inov4_xs_u    , jp3) = dhxu_jp3;
                dB(id_z4_xrd      , id_inov5_xf_xf_u , jp3) = 3/6*dhxxu_jp3;
                dB(id_z4_xrd      , id_inov6_xf_u_u  , jp3) = 3/6*dhxuu_jp3;
                dB(id_z4_xrd      , id_inov7_u_u_u   , jp3) = 1/6*dhuuu_jp3;
                dB(id_z5_xf_xs    , id_inov1_u       , jp3) = dhu_hss2_jp3;
                dB(id_z5_xf_xs    , id_inov4_xs_u    , jp3) = K_x_x*dhx_hu_jp3;
                dB(id_z5_xf_xs    , id_inov5_xf_xf_u , jp3) = dhx_hxu_jp3 + K_x_x*dhxx2_hu_jp3;
                dB(id_z5_xf_xs    , id_inov6_xf_u_u  , jp3) = dhx_huu2_jp3 + K_x_x*dhxu_hu_jp3;
                dB(id_z5_xf_xs    , id_inov7_u_u_u   , jp3) = dhu_huu2_jp3;
                dB(id_z6_xf_xf_xf , id_inov5_xf_xf_u , jp3) = dhx_hx_hu_jp3 + kron(dhx_jp3,K_x_x*hx_hu) + kron(hx,K_x_x*dhx_hu_jp3) + dhu_hx_hx_jp3*K_u_xx;
                dB(id_z6_xf_xf_xf , id_inov6_xf_u_u  , jp3) = dhx_hu_hu_jp3 + dhu_hx_hu_jp3*K_u_ux + kron(dhu_jp3,K_x_x*hx_hu)*K_u_ux + kron(hu,K_x_x*dhx_hu_jp3)*K_u_ux;
                dB(id_z6_xf_xf_xf , id_inov7_u_u_u   , jp3) = dhu_hu_hu_jp3;

                dc(id_z2_xs    , jp3) = 1/2*dhss_jp3 + 1/2*dhuu_jp3*E_uu(:) + 1/2*huu*dE_uu_jp3(:);
                dc(id_z3_xf_xf , jp3) = dhu_hu_jp3*E_uu(:) + hu_hu*dE_uu_jp3(:);

                dC(: , id_z1_xf       , jp3) = dgx_jp3 + 3/6*dgxss_jp3;
                dC(: , id_z2_xs       , jp3) = dgx_jp3;
                dC(: , id_z3_xf_xf    , jp3) = 1/2*dgxx_jp3;
                dC(: , id_z4_xrd      , jp3) = dgx_jp3;
                dC(: , id_z5_xf_xs    , jp3) = dgxx_jp3;
                dC(: , id_z6_xf_xf_xf , jp3) = 1/6*dgxxx_jp3;

                dD(: , id_inov1_u       , jp3) = dgu_jp3 + 3/6*dguss_jp3;
                dD(: , id_inov2_u_u     , jp3) = 1/2*dguu_jp3;
                dD(: , id_inov3_xf_u    , jp3) = dgxu_jp3;
                dD(: , id_inov4_xs_u    , jp3) = dgxu_jp3;
                dD(: , id_inov5_xf_xf_u , jp3) = 3/6*dgxxu_jp3;
                dD(: , id_inov6_xf_u_u  , jp3) = 3/6*dgxuu_jp3;
                dD(: , id_inov7_u_u_u   , jp3) = 1/6*dguuu_jp3;

                dd(:,jp3) = 1/2*dgss_jp3 + 1/2*dguu_jp3*E_uu(:) + 1/2*guu*dE_uu_jp3(:);

                dVarinov(id_inov1_u       , id_inov1_u       , jp3) = dE_uu_jp3;
                dVarinov(id_inov1_u       , id_inov4_xs_u    , jp3) = kron(dE_xs_jp3',E_uu) + kron(E_xs',dE_uu_jp3);
                dVarinov(id_inov1_u       , id_inov5_xf_xf_u , jp3) = kron(dE_xfxf_jp3(:)',E_uu) + kron(E_xfxf(:)',dE_uu_jp3);
                dVarinov(id_inov1_u       , id_inov7_u_u_u   , jp3) = reshape(dE_u_u_u_u_jp3,u_nbr,u_nbr^3);
                dVarinov(id_inov2_u_u     , id_inov2_u_u     , jp3) = reshape(dE_u_u_u_u_jp3,u_nbr^2,u_nbr^2) - dE_uu_jp3(:)*E_uu(:)' - E_uu(:)*dE_uu_jp3(:)';
                dVarinov(id_inov3_xf_u    , id_inov3_xf_u    , jp3) = dE_xfxf_uu_jp3;
                dVarinov(id_inov4_xs_u    , id_inov1_u       , jp3) = kron(dE_xs_jp3,E_uu) + kron(E_xs,dE_uu_jp3);
                dVarinov(id_inov4_xs_u    , id_inov4_xs_u    , jp3) = kron(dE_xsxs_jp3,E_uu) + kron(E_xsxs,dE_uu_jp3);
                dVarinov(id_inov4_xs_u    , id_inov5_xf_xf_u , jp3) = kron(dE_xsxf_xf_jp3, E_uu) + kron(E_xsxf_xf, dE_uu_jp3);
                dVarinov(id_inov4_xs_u    , id_inov7_u_u_u   , jp3) = kron(dE_xs_jp3,reshape(QPu*E_u_u_u_u,u_nbr,u_nbr^3)) + kron(E_xs,reshape(dE_u_u_u_u_jp3,u_nbr,u_nbr^3));
                dVarinov(id_inov5_xf_xf_u , id_inov1_u       , jp3) = kron(dE_xfxf_jp3(:),E_uu) + kron(E_xfxf(:),dE_uu_jp3);
                dVarinov(id_inov5_xf_xf_u , id_inov4_xs_u    , jp3) = kron(dE_xf_xfxs_jp3,E_uu) + kron(E_xf_xfxs,dE_uu_jp3);
                dVarinov(id_inov5_xf_xf_u , id_inov5_xf_xf_u , jp3) = kron(dE_xf_xfxf_xf_jp3,E_uu) + kron(E_xf_xfxf_xf,dE_uu_jp3);
                dVarinov(id_inov5_xf_xf_u , id_inov7_u_u_u   , jp3) = kron(dE_xfxf_jp3(:),reshape(QPu*E_u_u_u_u,u_nbr,u_nbr^3)) + kron(E_xfxf(:),reshape(dE_u_u_u_u_jp3,u_nbr,u_nbr^3));
                dVarinov(id_inov6_xf_u_u  , id_inov6_xf_u_u  , jp3) = kron(dE_xfxf_jp3,reshape(QPu*E_u_u_u_u,u_nbr^2,u_nbr^2)) + kron(E_xfxf,reshape(dE_u_u_u_u_jp3,u_nbr^2,u_nbr^2));
                dVarinov(id_inov7_u_u_u   , id_inov1_u       , jp3) = reshape(dE_u_u_u_u_jp3,u_nbr^3,u_nbr);
                dVarinov(id_inov7_u_u_u   , id_inov4_xs_u    , jp3) = kron(dE_xs_jp3',reshape(QPu*E_u_u_u_u,u_nbr^3,u_nbr)) + kron(E_xs',reshape(dE_u_u_u_u_jp3,u_nbr^3,u_nbr));
                dVarinov(id_inov7_u_u_u   , id_inov5_xf_xf_u , jp3) = kron(transpose(dE_xfxf_jp3(:)),reshape(QPu*E_u_u_u_u,u_nbr^3,u_nbr)) + kron(transpose(E_xfxf(:)),reshape(dE_u_u_u_u_jp3,u_nbr^3,u_nbr));
                dVarinov(id_inov7_u_u_u   , id_inov7_u_u_u   , jp3) = reshape(dE_u_u_u_u_u_u_jp3,u_nbr^3,u_nbr^3);

                dE_inovzlag1(id_inov6_xf_u_u , id_z1_xf       , jp3) = kron(dE_xfxf_jp3,E_uu(:)) + kron(E_xfxf,dE_uu_jp3(:));
                dE_inovzlag1(id_inov6_xf_u_u , id_z4_xrd      , jp3) = kron(dE_xrdxf_jp3',E_uu(:)) + kron(E_xrdxf',dE_uu_jp3(:));
                dE_inovzlag1(id_inov6_xf_u_u , id_z5_xf_xs    , jp3) = kron(reshape(K_xx_x*vec(dE_xsxf_xf_jp3),x_nbr,x_nbr^2),vec(E_uu)) + kron(reshape(K_xx_x*vec(E_xsxf_xf),x_nbr,x_nbr^2),vec(dE_uu_jp3)) ;
                dE_inovzlag1(id_inov6_xf_u_u , id_z6_xf_xf_xf , jp3) = kron(reshape(dE_xf_xfxf_xf_jp3,x_nbr,x_nbr^3),E_uu(:)) + kron(reshape(E_xf_xfxf_xf,x_nbr,x_nbr^3),dE_uu_jp3(:));

                dBinovzlag1A_jp3 = dB(:,:,jp3)*E_inovzlag1*transpose(A) + B*dE_inovzlag1(:,:,jp3)*transpose(A) + B*E_inovzlag1*transpose(dA(:,:,jp3));
                dOm_z_jp3 = dB(:,:,jp3)*Varinov*transpose(B) + B*dVarinov(:,:,jp3)*transpose(B) + B*Varinov*transpose(dB(:,:,jp3)) + dBinovzlag1A_jp3 + transpose(dBinovzlag1A_jp3);
                
                [dVar_z(:,:,jp3), errorflag] = disclyap_fast(A, dA(:,:,jp3)*Var_z*A' + A*Var_z*dA(:,:,jp3)' + dOm_z_jp3, options.lyapunov_doubling_tol);
                if errorflag
                    dVar_z(:,:,jp3) = lyapunov_symm(A, dA(:,:,jp3)*Var_z*A' + A*Var_z*dA(:,:,jp3)' + dOm_z_jp3,...
                                                    options.lyapunov_fixed_point_tol,options.qz_criterium,options.lyapunov_complex_threshold,...
                                                    lyapunov_symm_method,...
                                                    options.debug);
                    if lyapunov_symm_method == 1
                        lyapunov_symm_method = 2; %now we can reuse persistent schur
                    end
                end
                %make sure some stuff is zero due to Gaussianity
                dVar_z(id_z1_xf       , id_z2_xs       , jp3) = zeros(x_nbr,x_nbr);
                dVar_z(id_z1_xf       , id_z3_xf_xf    , jp3) = zeros(x_nbr,x_nbr^2);
                dVar_z(id_z2_xs       , id_z1_xf       , jp3) = zeros(x_nbr,x_nbr);
                dVar_z(id_z2_xs       , id_z4_xrd      , jp3) = zeros(x_nbr,x_nbr);
                dVar_z(id_z2_xs       , id_z5_xf_xs    , jp3) = zeros(x_nbr,x_nbr^2);
                dVar_z(id_z2_xs       , id_z6_xf_xf_xf , jp3) = zeros(x_nbr,x_nbr^3);
                dVar_z(id_z3_xf_xf    , id_z1_xf       , jp3) = zeros(x_nbr^2,x_nbr);
                dVar_z(id_z3_xf_xf    , id_z4_xrd      , jp3) = zeros(x_nbr^2,x_nbr);
                dVar_z(id_z3_xf_xf    , id_z5_xf_xs    , jp3) = zeros(x_nbr^2,x_nbr^2);
                dVar_z(id_z3_xf_xf    , id_z6_xf_xf_xf , jp3) = zeros(x_nbr^2,x_nbr^3);
                dVar_z(id_z4_xrd      , id_z2_xs       , jp3) = zeros(x_nbr,x_nbr);
                dVar_z(id_z4_xrd      , id_z3_xf_xf    , jp3) = zeros(x_nbr,x_nbr^2);
                dVar_z(id_z5_xf_xs    , id_z2_xs       , jp3) = zeros(x_nbr^2,x_nbr);
                dVar_z(id_z5_xf_xs    , id_z3_xf_xf    , jp3) = zeros(x_nbr^2,x_nbr^2);
                dVar_z(id_z6_xf_xf_xf , id_z2_xs       , jp3) = zeros(x_nbr^3,x_nbr);
                dVar_z(id_z6_xf_xf_xf , id_z3_xf_xf    , jp3) = zeros(x_nbr^3,x_nbr^2);
            end
        end
    end
end

%% Covariance/Correlation of control variables
Var_y = NaN*ones(y_nbr,y_nbr);
if order < 3
    Var_y(stationary_vars,stationary_vars) = C(stationary_vars,:)*Var_z*C(stationary_vars,:)'...
                                           + D(stationary_vars,:)*Varinov*D(stationary_vars,:)';
else
    Var_y(stationary_vars,stationary_vars) = C(stationary_vars,:)*Var_z*C(stationary_vars,:)'...
                                           + D(stationary_vars,:)*E_inovzlag1*C(stationary_vars,:)'...
                                           + C(stationary_vars,:)*transpose(E_inovzlag1)*D(stationary_vars,:)'...
                                           + D(stationary_vars,:)*Varinov*D(stationary_vars,:)';
end

Var_y(abs(Var_y) < 1e-12) = 0; %find values that are numerical zero
if useautocorr
    sdy = sqrt(diag(Var_y)); %theoretical standard deviation
    sdy = sdy(stationary_vars);
    sy = sdy*sdy';           %cross products of standard deviations
    Corr_y = NaN*ones(y_nbr,y_nbr);
    Corr_y(stationary_vars,stationary_vars) = Var_y(stationary_vars,stationary_vars)./sy;
    Corr_yi = NaN*ones(y_nbr,y_nbr,nlags);
end

if compute_derivs
    dVar_y = NaN*ones(y_nbr,y_nbr,totparam_nbr);
    if useautocorr
        dCorr_y  = NaN*ones(y_nbr,y_nbr,totparam_nbr);
        dCorr_yi = NaN*ones(y_nbr,y_nbr,nlags,totparam_nbr);
    end
    for jpV=1:totparam_nbr
        if order < 3
            dVar_y_tmp = dC(stationary_vars,:,jpV)*Var_z*C(stationary_vars,:)' + C(stationary_vars,:)*dVar_z(:,:,jpV)*C(stationary_vars,:)' + C(stationary_vars,:)*Var_z*dC(stationary_vars,:,jpV)'...
                                                        + dD(stationary_vars,:,jpV)*Varinov*D(stationary_vars,:)' + D(stationary_vars,:)*dVarinov(:,:,jpV)*D(stationary_vars,:)' + D(stationary_vars,:)*Varinov*dD(stationary_vars,:,jpV)';
        else
            dVar_y_tmp = dC(stationary_vars,:,jpV)*Var_z*C(stationary_vars,:)' + C(stationary_vars,:)*dVar_z(:,:,jpV)*C(stationary_vars,:)' + C(stationary_vars,:)*Var_z*dC(stationary_vars,:,jpV)'...
                                                        + dD(stationary_vars,:,jpV)*E_inovzlag1*C(stationary_vars,:)' + D(stationary_vars,:)*dE_inovzlag1(:,:,jpV)*C(stationary_vars,:)' + D(stationary_vars,:)*E_inovzlag1*dC(stationary_vars,:,jpV)'...
                                                        + dC(stationary_vars,:,jpV)*transpose(E_inovzlag1)*D(stationary_vars,:)' + C(stationary_vars,:)*transpose(dE_inovzlag1(:,:,jpV))*D(stationary_vars,:)' + C(stationary_vars,:)*transpose(E_inovzlag1)*dD(stationary_vars,:,jpV)'...
                                                        + dD(stationary_vars,:,jpV)*Varinov*D(stationary_vars,:)' + D(stationary_vars,:)*dVarinov(:,:,jpV)*D(stationary_vars,:)' + D(stationary_vars,:)*Varinov*dD(stationary_vars,:,jpV)';
        end
        dVar_y_tmp(abs(dVar_y_tmp) < 1e-12) = 0; %find values that are numerical zero       
        dVar_y(stationary_vars,stationary_vars,jpV) = dVar_y_tmp;
        if useautocorr
            dsy = 1/2./sdy.*diag(dVar_y(:,:,jpV));
            dsy = dsy(stationary_vars);
            dsy = dsy*sdy'+sdy*dsy';
            dCorr_y(stationary_vars,stationary_vars,jpV) = (dVar_y(stationary_vars,stationary_vars,jpV).*sy-dsy.*Var_y(stationary_vars,stationary_vars))./(sy.*sy);
            dCorr_y(stationary_vars,stationary_vars,jpV) = dCorr_y(stationary_vars,stationary_vars,jpV)-diag(diag(dCorr_y(stationary_vars,stationary_vars,jpV)))+diag(diag(dVar_y(stationary_vars,stationary_vars,jpV)));
        end
    end
end

%% Autocovariances/autocorrelations of lagged control variables
Var_yi = NaN*ones(y_nbr,y_nbr,nlags);
Ai = eye(z_nbr); %this is A^0
hxi = eye(x_nbr);
E_inovzlagi = E_inovzlag1;
Var_zi = Var_z;
if order <= 2
    tmp = A*Var_z*C(stationary_vars,:)' + B*Varinov*D(stationary_vars,:)';
else
    tmp = A*E_inovzlag1'*D(stationary_vars,:)' + B*Varinov*D(stationary_vars,:)';
end
for i = 1:nlags
    if order <= 2
        Var_yi(stationary_vars,stationary_vars,i) = C(stationary_vars,:)*Ai*tmp;
    else
        Var_zi = A*Var_zi + B*E_inovzlagi;
        hxi = hx*hxi;
        E_inovzlagi = zeros(inov_nbr,z_nbr);
        E_inovzlagi(id_inov6_xf_u_u , id_z1_xf       ) = kron(hxi*E_xfxf,E_uu(:));
        E_inovzlagi(id_inov6_xf_u_u , id_z4_xrd      ) = kron(hxi*E_xrdxf',E_uu(:));
        E_inovzlagi(id_inov6_xf_u_u , id_z5_xf_xs    ) = kron(hxi*reshape(K_xx_x*vec(E_xsxf_xf),x_nbr,x_nbr^2),vec(E_uu));
        E_inovzlagi(id_inov6_xf_u_u , id_z6_xf_xf_xf ) = kron(hxi*reshape(E_xf_xfxf_xf,x_nbr,x_nbr^3),E_uu(:));
        Var_yi(stationary_vars,stationary_vars,i)      = C(stationary_vars,:)*Var_zi*C(stationary_vars,:)' + C(stationary_vars,:)*Ai*tmp + D(stationary_vars,:)*E_inovzlagi*C(stationary_vars,:)';
    end    
    if useautocorr
        Corr_yi(stationary_vars,stationary_vars,i) = Var_yi(stationary_vars,stationary_vars,i)./sy;
    end
    Ai = Ai*A; %note that this is A^(i-1)
end

if compute_derivs
    dVar_yi = NaN*ones(y_nbr,y_nbr,nlags,totparam_nbr);
    for jpVi=1:totparam_nbr        
        Ai          = eye(z_nbr);   dAi_jpVi          = zeros(z_nbr,z_nbr);
        hxi         = eye(x_nbr);   dhxi_jpVi         = zeros(x_nbr,x_nbr);
        E_inovzlagi = E_inovzlag1;  dE_inovzlagi_jpVi = dE_inovzlag1(:,:,jpVi);
        Var_zi      = Var_z;        dVar_zi_jpVi      = dVar_z(:,:,jpVi);
        if order <= 2            
            dtmp_jpVi = dA(:,:,jpVi)*Var_z*C(stationary_vars,:)' + A*dVar_z(:,:,jpVi)*C(stationary_vars,:)' + A*Var_z*dC(stationary_vars,:,jpVi)'...
                      + dB(:,:,jpVi)*Varinov*D(stationary_vars,:)' + B*dVarinov(:,:,jpVi)*D(stationary_vars,:)' + B*Varinov*dD(stationary_vars,:,jpVi)';
        else
            dtmp_jpVi = dA(:,:,jpVi)*E_inovzlag1'*D(stationary_vars,:)' + A*dE_inovzlag1(:,:,jpVi)'*D(stationary_vars,:)' + A*E_inovzlag1'*dD(stationary_vars,:,jpVi)'...
                      + dB(:,:,jpVi)*Varinov*D(stationary_vars,:)' + B*dVarinov(:,:,jpVi)*D(stationary_vars,:)' + B*Varinov*dD(stationary_vars,:,jpVi)';
        end

        for i = 1:nlags
            if order <= 2
                dVar_yi(stationary_vars,stationary_vars,i,jpVi) = dC(stationary_vars,:,jpVi)*Ai*tmp + C(stationary_vars,:)*dAi_jpVi*tmp + C(stationary_vars,:)*Ai*dtmp_jpVi;
            else
                Var_zi       = A*Var_zi + B*E_inovzlagi;
                dVar_zi_jpVi = dA(:,:,jpVi)*Var_zi + A*dVar_zi_jpVi + dB(:,:,jpVi)*E_inovzlagi + + B*dE_inovzlagi_jpVi;
                dhxi_jpVi = dhx(:,:,jpVi)*hxi + hx*dhxi_jpVi;
                hxi = hx*hxi;                
                E_inovzlagi = zeros(inov_nbr,z_nbr);
                E_inovzlagi(id_inov6_xf_u_u , id_z1_xf       ) = kron(hxi*E_xfxf,E_uu(:));
                E_inovzlagi(id_inov6_xf_u_u , id_z4_xrd      ) = kron(hxi*E_xrdxf',E_uu(:));
                E_inovzlagi(id_inov6_xf_u_u , id_z5_xf_xs    ) = kron(hxi*reshape(K_xx_x*vec(E_xsxf_xf),x_nbr,x_nbr^2),vec(E_uu));
                E_inovzlagi(id_inov6_xf_u_u , id_z6_xf_xf_xf ) = kron(hxi*reshape(E_xf_xfxf_xf,x_nbr,x_nbr^3),E_uu(:));
                dE_inovzlagi_jpVi = zeros(inov_nbr,z_nbr);
                dE_inovzlagi_jpVi(id_inov6_xf_u_u , id_z1_xf       ) = kron(dhxi_jpVi*E_xfxf,E_uu(:)) + kron(hxi*dE_xfxf(:,:,jpVi),E_uu(:)) + kron(hxi*E_xfxf,vec(dE_uu(:,:,jpVi)));
                dE_inovzlagi_jpVi(id_inov6_xf_u_u , id_z4_xrd      ) = kron(dhxi_jpVi*E_xrdxf',E_uu(:)) + kron(hxi*dE_xrdxf(:,:,jpVi)',E_uu(:)) + kron(hxi*E_xrdxf',vec(dE_uu(:,:,jpVi)));
                dE_inovzlagi_jpVi(id_inov6_xf_u_u , id_z5_xf_xs    ) = kron(dhxi_jpVi*reshape(K_xx_x*vec(E_xsxf_xf),x_nbr,x_nbr^2),vec(E_uu)) + kron(hxi*reshape(K_xx_x*vec(dE_xsxf_xf(:,:,jpVi)),x_nbr,x_nbr^2),vec(E_uu)) + kron(hxi*reshape(K_xx_x*vec(E_xsxf_xf),x_nbr,x_nbr^2),vec(dE_uu(:,:,jpVi)));
                dE_inovzlagi_jpVi(id_inov6_xf_u_u , id_z6_xf_xf_xf ) = kron(dhxi_jpVi*reshape(E_xf_xfxf_xf,x_nbr,x_nbr^3),E_uu(:)) + kron(hxi*reshape(dE_xf_xfxf_xf(:,:,jpVi),x_nbr,x_nbr^3),E_uu(:)) + kron(hxi*reshape(E_xf_xfxf_xf,x_nbr,x_nbr^3),vec(dE_uu(:,:,jpVi)));
                dVar_yi(stationary_vars,stationary_vars,i,jpVi) = dC(stationary_vars,:,jpVi)*Var_zi*C(stationary_vars,:)' + C(stationary_vars,:)*dVar_zi_jpVi*C(stationary_vars,:)' + C(stationary_vars,:)*Var_zi*dC(stationary_vars,:,jpVi)'...
                                                                + dC(stationary_vars,:,jpVi)*Ai*tmp + C(stationary_vars,:)*dAi_jpVi*tmp + C(stationary_vars,:)*Ai*dtmp_jpVi...
                                                                + dD(stationary_vars,:,jpVi)*E_inovzlagi*C(stationary_vars,:)' + D(stationary_vars,:)*dE_inovzlagi_jpVi*C(stationary_vars,:)' + D(stationary_vars,:)*E_inovzlagi*dC(stationary_vars,:,jpVi)';
            end
            if useautocorr
                dsy = 1/2./sdy.*diag(dVar_y(:,:,jpVi));
                dsy = dsy(stationary_vars);
                dsy = dsy*sdy'+sdy*dsy';
                dCorr_yi(stationary_vars,stationary_vars,i,jpVi) = (dVar_yi(stationary_vars,stationary_vars,i,jpVi).*sy-dsy.*Var_yi(stationary_vars,stationary_vars,i))./(sy.*sy);                
            end
            dAi_jpVi = dAi_jpVi*A + Ai*dA(:,:,jpVi);
            Ai = Ai*A;
        end
    end
end    


%% Mean of control variables
E_z = E_xf;
if order > 1
    E_z = [E_xf;E_xs;E_xfxf(:)];
end
if order > 2
    E_xf_xs = zeros(x_nbr^2,1);
    E_xf_xf_xf = zeros(x_nbr^3,1);
    E_z = [E_xf;E_xs;E_xfxf(:);E_xrd;E_xf_xs;E_xf_xf_xf];
end
E_y  = Yss(indy,:) + C*E_z + d;

if compute_derivs
    dE_y = zeros(y_nbr,totparam_nbr);
    for jpE = 1:totparam_nbr
        if order == 1
            dE_z_jpE = dE_xf(:,jpE);
        elseif order == 2
            dE_z_jpE = [dE_xf(:,jpE);dE_xs(:,jpE);vec(dE_xfxf(:,:,jpE))];
        elseif order == 3
            dE_xf_xs_jpE    = zeros(x_nbr^2,1);
            dE_xf_xf_xf_jpE = zeros(x_nbr^3,1);
            dE_z_jpE        = [dE_xf(:,jpE);dE_xs(:,jpE);vec(dE_xfxf(:,:,jpE)); dE_xrd(:,jpE); dE_xf_xs_jpE; dE_xf_xf_xf_jpE];
        end
        dE_y(:,jpE) = dC(:,:,jpE)*E_z + C*dE_z_jpE + dd(:,jpE);
        if jpE > (stderrparam_nbr+corrparam_nbr)
            dE_y(:,jpE) = dE_y(:,jpE) + dYss(indy,jpE-stderrparam_nbr-corrparam_nbr); %add steady state
        end
    end
end
non_stationary_vars = ~ismember((1:y_nbr)',stationary_vars);
E_y(non_stationary_vars) = NaN;
if compute_derivs
    dE_y(non_stationary_vars,:) = NaN;
end

%% Store into output structure
pruned_state_space.indx = indx;
pruned_state_space.indy = indy;
pruned_state_space.A = A;
pruned_state_space.B = B;
pruned_state_space.C = C;
pruned_state_space.D = D;
pruned_state_space.c = c;
pruned_state_space.d = d;
pruned_state_space.Varinov = Varinov;
% pruned_state_space.Var_z  = Var_z; %
pruned_state_space.Var_y  = Var_y;
pruned_state_space.Var_yi = Var_yi;
if useautocorr
    pruned_state_space.Corr_y  = Corr_y;
    pruned_state_space.Corr_yi = Corr_yi;
end
pruned_state_space.E_y = E_y;

if compute_derivs == 1
    pruned_state_space.dA = dA;
    pruned_state_space.dB = dB;
    pruned_state_space.dC = dC;
    pruned_state_space.dD = dD;
    pruned_state_space.dc = dc;
    pruned_state_space.dd = dd;
    pruned_state_space.dVarinov = dVarinov;
    pruned_state_space.dVar_y   = dVar_y;
    pruned_state_space.dVar_yi  = dVar_yi;
    if useautocorr
        pruned_state_space.dCorr_y  = dCorr_y;
        pruned_state_space.dCorr_yi = dCorr_yi;
    end
    pruned_state_space.dE_y = dE_y;
end
