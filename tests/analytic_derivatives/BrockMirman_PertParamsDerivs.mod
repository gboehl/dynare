% "Augmented" Brock Mirman model
% True policy functions and their exact derivatives (ghx,ghu,ghxx,ghxu,ghuu,ghs2,ghxxx,ghxxu,ghxuu,ghuuu,ghxss,ghuss) are computed using Matlab's symbolic toolbox and saved to a mat file
% Created by @wmutschl (Willi Mutschler, willi@mutschler.eu)

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


@#define CREATE_SYMBOLIC = 0
@#define ORDER = 3

%define parameter values which are used for calibration and estimated_params block
@#define value_SE_A     = 0.3
@#define value_SE_X     = 0.1
@#define value_SE_Y     = 0.9
@#define value_rho_AX   = 0.4
@#define value_alph     = 0.35
@#define value_betta    = 0.99
@#define value_rhoA     = 0.9
@#define value_sigA     = 0.6
@#define value_sigX     = 0.25
@#define value_Xss      = 1.2
@#define value_dumA     = 2
@#define value_dumK     = 1
@#define value_dumepsA  = 4
@#define value_dumepsX  = 5


%define parameter values which are used for estimated_params block, note that all are different
@#define prior_value_SE_A     = 0.8
@#define prior_value_SE_X     = 0.6
@#define prior_value_SE_Y     = 0.1
@#define prior_value_rho_AX   = 0.1
@#define prior_value_alph     = 0.15
@#define prior_value_betta    = 0.92
@#define prior_value_rhoA     = 0.3
@#define prior_value_sigA     = 0.3
@#define prior_value_sigX     = 0.15
@#define prior_value_Xss      = 1.1
@#define prior_value_dumA     = 2
@#define prior_value_dumK     = 1
@#define prior_value_dumepsA  = 5
@#define prior_value_dumepsX  = 6

var Y X K C A W Z; %note that declaration order is different from DR order on purpose
varobs C K A;
varexo epsA epsX epsY;
parameters alph betta rhoA sigA sigX Xss dumA dumK dumepsA dumepsX;

%Calibration
alph    = @{value_alph};
betta   = @{value_betta};
rhoA    = @{value_rhoA};
sigA    = @{value_sigA};
sigX    = @{value_sigX};
Xss     = @{value_Xss};
dumA    = @{value_dumA};
dumK    = @{value_dumK};
dumepsA = @{value_dumepsA};
dumepsX = @{value_dumepsX};

model;
%this is original Brock-Mirman model equations with AR(1) shock
C^(-1) = alph*betta*C(+1)^(-1)*A(+1)*K^(alph-1);
K = A*K(-1)^alph - C;
log(A) = rhoA*log(A(-1)) + sigA*epsA;
Y = A*K(-1)^alph + epsY;
%augmented auxiliary equations to get nonzero ghs2, ghxss and ghuss, they have no economic meaning
log(X) = log(Xss) + sigX*epsX;
W = X(+1)*exp(dumA*A(-1))*exp(dumK*K(-1)); %this will get a nonzero ghs2 and ghxss
Z = X(+1)*exp(dumepsA*epsA)*exp(dumepsX*epsX); %this will get a nonzero ghs2 and ghuss
end;

shocks;
var epsA = (@{value_SE_A})^2;
var epsX = (@{value_SE_X})^2;
var epsA, epsX = @{value_rho_AX}*@{value_SE_A}*@{value_SE_X};
var epsY = (@{value_SE_Y})^2;
end;

steady_state_model;
X = Xss;
A = 1;
K = (alph*betta*A)^(1/(1-alph));
C = A*K^alph-K;
Y = A*K^alph;
W = Xss*exp(dumA*A)*exp(dumK*K);
Z = Xss*exp(dumepsA*0)*exp(dumepsX*0);
end;

estimated_params;%note that the parameter ordering is different than declaration order
rhoA,            normal_pdf, @{prior_value_rhoA},    0.1;
betta,           normal_pdf, @{prior_value_betta},   0.1;
alph,            normal_pdf, @{prior_value_alph},    0.1;
corr epsA, epsX, normal_pdf, @{prior_value_rho_AX},  0.1;
sigA,            normal_pdf, @{prior_value_sigA},    0.1;
stderr epsA,     normal_pdf, @{prior_value_SE_A},    0.1;
stderr epsX,     normal_pdf, @{prior_value_SE_X},    0.1;
stderr epsY,     normal_pdf, @{prior_value_SE_Y},    0.1;
sigX,            normal_pdf, @{prior_value_sigX},    0.1;
Xss,             normal_pdf, @{prior_value_Xss},     0.1;
dumepsX,         normal_pdf, @{prior_value_dumepsX}, 0.1;
dumK,            normal_pdf, @{prior_value_dumK},    0.1;
dumepsA,         normal_pdf, @{prior_value_dumepsA}, 0.1;
dumA,            normal_pdf, @{prior_value_dumA},    0.1;
end;

%save calibration parameters as these get overwritten through stoch_simul and identification command
calib_params = M_.params;
calib_Sigma_e = M_.Sigma_e;

stoch_simul(order=@{ORDER},nograph,irf=0,periods=0);
identification(order=@{ORDER},nograph,no_identification_strength);

indpmodel  = estim_params_.param_vals(:,1);
indpstderr = estim_params_.var_exo(:,1);
indpcorr   = estim_params_.corrx(:,1:2);
[I,~]      = find(M_.lead_lag_incidence');

%% Parameter derivatives of perturbation 
@#if CREATE_SYMBOLIC == 1
    syms Y_ Y0 Yp C_ C0 Cp K_ K0 Kp A_ A0 Ap X_ X0 Xp W_ W0 Wp Z_ Z0 Zp;
    syms epsA0 epsX0 epsY0;
    syms RHO_AX SE_A SE_X SE_Y ALPH BETTA RHOA SIGA SIGX XSS DUMA DUMK DUMEPSA DUMEPSX;
    syms sig;
    
    SYM.corr_params        = [RHO_AX];
    SYM.stderr_params_decl = [SE_A SE_X SE_Y];
    SYM.modparams_decl     = [ALPH BETTA RHOA SIGA SIGX XSS DUMA DUMK DUMEPSA DUMEPSX];
    SYM.endo_decl_  = [Y_ X_ K_ C_ A_ W_ Z_];
    SYM.endo_decl0  = [Y0 X0 K0 C0 A0 W0 Z0];
    SYM.endo_declp  = [Yp Xp Kp Cp Ap Wp Zp];
    SYM.ex0         = [epsA0 epsX0 epsY0];
    SYM.model_eqs   = [C0^(-1)-ALPH*BETTA*Cp^(-1)*Ap*K0^(ALPH-1);
                       K0-A0*K_^ALPH+C0;
                       log(A0)-RHOA*log(A_)-SIGA*epsA0;
                       Y0 - A0*K_^ALPH - epsY0;
                       log(X0)-log(XSS)-SIGX*epsX0;
                       W0 - Xp*exp(DUMA*A_)*exp(DUMK*K_);
                       Z0 - Xp*exp(DUMEPSA*epsA0)*exp(DUMEPSX*epsX0);
                       ];
    SYM.Sigma_e  = [SE_A^2 RHO_AX*SE_A*SE_X 0; RHO_AX*SE_A*SE_X  SE_X^2 0; 0 0 SE_Y^2];
    SYM.Correlation_matrix  = [1 RHO_AX 0; RHO_AX 1 0; 0 0 1];
    %steady states
    SYM.epsAbar = sym(0);
    SYM.epsXbar = sym(0);
    SYM.epsYbar = sym(0);
    SYM.Abar    = sym(1);
    SYM.Xbar    = XSS;
    SYM.Kbar    = (ALPH*BETTA*SYM.Abar)^(1/(1-ALPH));
    SYM.Cbar    = SYM.Abar*SYM.Kbar^ALPH-SYM.Kbar;
    SYM.Ybar    = SYM.Abar*SYM.Kbar^ALPH;
    SYM.Wbar    = XSS*exp(DUMA*SYM.Abar)*exp(DUMK*SYM.Kbar);
    SYM.Zbar    = XSS*exp(DUMEPSA*0)*exp(DUMEPSX*0);
    SYM.endo_bar_decl = [SYM.Ybar SYM.Xbar SYM.Kbar SYM.Cbar SYM.Abar SYM.Wbar SYM.Zbar];
    SYM.ex0bar  = [SYM.epsAbar SYM.epsXbar SYM.epsYbar];
    %True policy functions (in declaration order): gs is used for derivatives wrt to perturbation parameter
    SYM.g = [A_^RHOA*K_^ALPH*exp(SIGA*epsA0) + epsY0;               %Y
             XSS*exp(SIGX*epsX0);                                   %X
             ALPH*BETTA*A_^RHOA*K_^ALPH*exp(SIGA*epsA0);            %K
             (1-ALPH*BETTA)*A_^RHOA*K_^ALPH*exp(SIGA*epsA0);        %C
             A_^RHOA*exp(SIGA*epsA0);                               %A
             XSS*exp(SIGX*0)*exp(DUMA*A_)*exp(DUMK*K_);             %W
             XSS*exp(SIGX*0)*exp(DUMEPSA*epsA0)*exp(DUMEPSX*epsX0); %Z
             ];
   SYM.gs = [A_^RHOA*K_^ALPH*exp(SIGA*sig*epsA0) + sig*epsY0;                       %Y
             XSS*exp(SIGX*sig*epsX0);                                               %X
             ALPH*BETTA*A_^RHOA*K_^ALPH*exp(SIGA*sig*epsA0);                        %K
             (1-ALPH*BETTA)*A_^RHOA*K_^ALPH*exp(SIGA*sig*epsA0);                    %C
             A_^RHOA*exp(SIGA*sig*epsA0);                                           %A
             XSS*(1+1/2*SIGX^2*sig^2*SE_X^2)*exp(DUMA*A_)*exp(DUMK*K_);             %W
             XSS*(1+1/2*SIGX^2*sig^2*SE_X^2)*exp(DUMEPSA*epsA0)*exp(DUMEPSX*epsX0); %Z
             ];
    %put into in DR-order
    SYM.g = SYM.g(oo_.dr.order_var);
    SYM.gs = SYM.gs(oo_.dr.order_var);
    SYM.endo_DR_ = SYM.endo_decl_(oo_.dr.order_var);
    SYM.endo_DR0 = SYM.endo_decl0(oo_.dr.order_var);
    SYM.endo_DRp = SYM.endo_declp(oo_.dr.order_var);
    SYM.Yss = SYM.endo_bar_decl(oo_.dr.order_var);
    SYM.yy0 = [SYM.endo_decl_(M_.lead_lag_incidence(1,:)~=0) SYM.endo_decl0(M_.lead_lag_incidence(2,:)~=0) SYM.endo_declp(M_.lead_lag_incidence(3,:)~=0)];
    SYM.yy0ex0 = [SYM.yy0 SYM.ex0];
    SYM.yy0ex0bar = [SYM.endo_bar_decl(I) SYM.ex0bar];
    SYM.x = SYM.endo_DR_(M_.nstatic + (1:M_.nspred));
    SYM.xbar = SYM.Yss(M_.nstatic + (1:M_.nspred));
    SYM.stderr_params = SYM.stderr_params_decl(indpstderr);
    SYM.modparams = SYM.modparams_decl(indpmodel);
    SYM.params = [SYM.stderr_params SYM.corr_params SYM.modparams];
    SYM.yy0ex0_nbr = length(SYM.yy0ex0);

    %% Parameter derivatives
    SYM.dYss = jacobian(SYM.Yss,SYM.modparams);
    SYM.d2Yss = jacobian(SYM.dYss(:),SYM.modparams);
    SYM.g1 = jacobian(SYM.model_eqs,SYM.yy0ex0);
    SYM.g2 = reshape(jacobian(vec(SYM.g1),SYM.yy0ex0),M_.endo_nbr,SYM.yy0ex0_nbr^2);
    for j = 1:M_.endo_nbr
        SYM.g3(j,:) = vec(transpose(jacobian(vec(SYM.g2(j,:)),SYM.yy0ex0))); %dynare ordering
    end
    SYM.dg1 = jacobian(vec(subs(SYM.g1,SYM.yy0ex0,SYM.yy0ex0bar)),SYM.modparams);
    SYM.dg2 = jacobian(vec(subs(SYM.g2,SYM.yy0ex0,SYM.yy0ex0bar)),SYM.modparams);
    SYM.dg3 = jacobian(vec(subs(SYM.g3,SYM.yy0ex0,SYM.yy0ex0bar)),SYM.modparams);
    SYM.dSigma_e = jacobian(SYM.Sigma_e(:),SYM.params);
    SYM.dCorrelation_matrix = jacobian(SYM.Correlation_matrix(:),SYM.params);
    SYM.ghx = jacobian(SYM.g,SYM.x);
    SYM.KalmanA = sym(zeros(M_.endo_nbr,M_.endo_nbr));
    SYM.KalmanA(:,M_.nstatic + (1:M_.nspred)) = SYM.ghx;
    SYM.dghx = jacobian(vec(subs(SYM.ghx, [SYM.x SYM.ex0], [SYM.xbar SYM.ex0bar])),SYM.params);
    SYM.dKalmanA = jacobian(vec(subs(SYM.KalmanA, [SYM.x SYM.ex0], [SYM.xbar SYM.ex0bar])),SYM.params);
    SYM.d2KalmanA = jacobian(SYM.dKalmanA(:),SYM.params);
    SYM.ghu = jacobian(SYM.g,SYM.ex0);
    SYM.dghu = jacobian(vec(subs(SYM.ghu, [SYM.x SYM.ex0], [SYM.xbar SYM.ex0bar])),SYM.params);
    SYM.Om = SYM.ghu*SYM.Sigma_e*transpose(SYM.ghu);
    SYM.dOm = jacobian(vec(subs(SYM.Om,[SYM.x SYM.ex0], [SYM.xbar SYM.ex0bar])),SYM.params);
    SYM.d2Om = jacobian(SYM.dOm(:),SYM.params);
    SYM.ghxx = jacobian(SYM.ghx(:),SYM.x);
    SYM.ghxx = reshape(SYM.ghxx,M_.endo_nbr,M_.nspred^2);
    SYM.dghxx = jacobian(vec(subs(SYM.ghxx, [SYM.x SYM.ex0], [SYM.xbar SYM.ex0bar])),SYM.params);
    SYM.ghxu = jacobian(SYM.ghu(:),SYM.x);
    SYM.ghxu = reshape(SYM.ghxu,M_.endo_nbr,M_.nspred*M_.exo_nbr);
    SYM.dghxu = jacobian(vec(subs(SYM.ghxu, [SYM.x SYM.ex0], [SYM.xbar SYM.ex0bar])),SYM.params);
    SYM.ghuu = jacobian(SYM.ghu(:),SYM.ex0);
    SYM.ghuu = reshape(SYM.ghuu,M_.endo_nbr,M_.exo_nbr*M_.exo_nbr);
    SYM.dghuu = jacobian(vec(subs(SYM.ghuu, [SYM.x SYM.ex0], [SYM.xbar SYM.ex0bar])),SYM.params);
    SYM.ghs2 = jacobian(jacobian(SYM.gs,sig),sig);
    SYM.dghs2 = jacobian(vec(subs(SYM.ghs2, [SYM.x SYM.ex0], [SYM.xbar SYM.ex0bar])),SYM.params);
    for j = 1:M_.endo_nbr
        SYM.ghxxx(j,:) = vec(transpose(jacobian(vec(SYM.ghxx(j,:)),SYM.x)));   %dynare ordering
        SYM.ghxxu(j,:) = vec(transpose(jacobian(vec(SYM.ghxx(j,:)),SYM.ex0))); %dynare ordering
        SYM.ghxuu(j,:) = vec(transpose(jacobian(vec(SYM.ghxu(j,:)),SYM.ex0))); %dynare ordering
        SYM.ghuuu(j,:) = vec(transpose(jacobian(vec(SYM.ghuu(j,:)),SYM.ex0))); %dynare ordering
        SYM.ghxss(j,:) = vec(transpose(jacobian(vec(SYM.ghs2(j,:)),SYM.x)));   %dynare ordering
        SYM.ghuss(j,:) = vec(transpose(jacobian(vec(SYM.ghs2(j,:)),SYM.ex0))); %dynare ordering
    end
    SYM.dghxxx = jacobian(vec(subs(SYM.ghxxx, [SYM.x SYM.ex0], [SYM.xbar SYM.ex0bar])),SYM.params);
    SYM.dghxxu = jacobian(vec(subs(SYM.ghxxu, [SYM.x SYM.ex0], [SYM.xbar SYM.ex0bar])),SYM.params);
    SYM.dghxuu = jacobian(vec(subs(SYM.ghxuu, [SYM.x SYM.ex0], [SYM.xbar SYM.ex0bar])),SYM.params);
    SYM.dghuuu = jacobian(vec(subs(SYM.ghuuu, [SYM.x SYM.ex0], [SYM.xbar SYM.ex0bar])),SYM.params);
    SYM.dghxss = jacobian(vec(subs(SYM.ghxss, [SYM.x SYM.ex0], [SYM.xbar SYM.ex0bar])),SYM.params);
    SYM.dghuss = jacobian(vec(subs(SYM.ghuss, [SYM.x SYM.ex0], [SYM.xbar SYM.ex0bar])),SYM.params);

    % Evaluate numerically
    for jj = 1:2
        if jj == 1
            RHO_AX = @{prior_value_rho_AX};
            SE_A  = @{prior_value_SE_A};
            SE_X  = @{prior_value_SE_X};
            SE_Y  = @{prior_value_SE_Y};
            ALPH = @{prior_value_alph};
            BETTA = @{prior_value_betta};
            RHOA = @{prior_value_rhoA};
            SIGA = @{prior_value_sigA};
            SIGX = @{prior_value_sigX};
            XSS = @{prior_value_Xss};
            DUMA = @{prior_value_dumA};
            DUMK = @{prior_value_dumK};
            DUMEPSA = @{prior_value_dumepsA};
            DUMEPSX = @{prior_value_dumepsX};
        elseif jj == 2
            RHO_AX = @{value_rho_AX};
            SE_A  = @{value_SE_A};
            SE_X  = @{value_SE_X};
            SE_Y  = @{value_SE_Y};
            ALPH = @{value_alph};
            BETTA = @{value_betta};
            RHOA = @{value_rhoA};
            SIGA = @{value_sigA};
            SIGX = @{value_sigX};
            XSS = @{value_Xss};
            DUMA = @{value_dumA};
            DUMK = @{value_dumK};
            DUMEPSA = @{value_dumepsA};
            DUMEPSX = @{value_dumepsX};
        end
        sig   = 1;

        nSYM.Yss   = transpose(double(subs(SYM.Yss)));
        nSYM.dYss  = reshape(double(subs(SYM.dYss)), M_.endo_nbr, length(SYM.modparams));
        nSYM.d2Yss = double(reshape(subs(SYM.d2Yss), [M_.endo_nbr length(SYM.modparams) length(SYM.modparams)]));
        nSYM.g1    = double(subs(subs(SYM.g1, SYM.yy0ex0, SYM.yy0ex0bar)));
        nSYM.dg1   = reshape(double(subs(SYM.dg1)),M_.endo_nbr, SYM.yy0ex0_nbr, length(SYM.modparams));
        nSYM.g2    = sparse(double(subs(subs(SYM.g2, SYM.yy0ex0, SYM.yy0ex0bar))));
        nSYM.dg2   = reshape(sparse(double(subs(SYM.dg2))), M_.endo_nbr, SYM.yy0ex0_nbr^2*length(SYM.modparams));
        nSYM.g3    = sparse(double(subs(subs(SYM.g3, SYM.yy0ex0, SYM.yy0ex0bar))));
        nSYM.dg3   = reshape(sparse(double(subs(SYM.dg3))), M_.endo_nbr, SYM.yy0ex0_nbr^3*length(SYM.modparams));
        nSYM.Sigma_e             = double(subs(SYM.Sigma_e));
        nSYM.dSigma_e            = double(reshape(subs(SYM.dSigma_e), M_.exo_nbr, M_.exo_nbr,length(SYM.params)));
        nSYM.Correlation_matrix  = double(subs(SYM.Correlation_matrix));
        nSYM.dCorrelation_matrix = double(reshape(subs(SYM.dCorrelation_matrix), M_.exo_nbr, M_.exo_nbr,length(SYM.params)));
        nSYM.ghx       = double(subs(subs(SYM.ghx, SYM.yy0ex0, SYM.yy0ex0bar)));
        nSYM.dghx      = reshape(double(subs(SYM.dghx)), M_.endo_nbr, M_.nspred, length(SYM.params));
        nSYM.ghu       = double(subs(subs(SYM.ghu, SYM.yy0ex0, SYM.yy0ex0bar)));
        nSYM.dghu      = reshape(double(subs(SYM.dghu)), M_.endo_nbr, M_.exo_nbr, length(SYM.params));
        nSYM.Om        = double(subs(subs(SYM.Om, SYM.yy0ex0, SYM.yy0ex0bar)));
        nSYM.dOm       = reshape(double(subs(SYM.dOm)), M_.endo_nbr, M_.endo_nbr, length(SYM.params));
        nSYM.d2Om      = reshape(sparse(double(subs(SYM.d2Om))), M_.endo_nbr^2, length(SYM.params)^2);
        SYM.indp2      = reshape(1:length(SYM.params)^2,length(SYM.params),length(SYM.params));
        SYM.indx2      = reshape(1:M_.endo_nbr^2, M_.endo_nbr, M_.endo_nbr);
        nSYM.d2Om      = nSYM.d2Om(dyn_vech(SYM.indx2),dyn_vech(SYM.indp2));
        nSYM.KalmanA   = double(subs(subs(SYM.KalmanA, SYM.yy0ex0, SYM.yy0ex0bar)));
        nSYM.dKalmanA  = reshape(double(subs(SYM.dKalmanA)), M_.endo_nbr, M_.endo_nbr, length(SYM.params));
        nSYM.d2KalmanA = reshape(double(subs(SYM.d2KalmanA)), M_.endo_nbr^2, length(SYM.params)^2);
        nSYM.d2KalmanA = nSYM.d2KalmanA(:,dyn_vech(SYM.indp2));
        nSYM.ghxx      = double(subs(subs(SYM.ghxx, SYM.yy0ex0, SYM.yy0ex0bar)));
        nSYM.dghxx     = reshape(double(subs(SYM.dghxx)), M_.endo_nbr, M_.nspred^2, length(SYM.params));
        nSYM.ghxu      = double(subs(subs(SYM.ghxu, SYM.yy0ex0, SYM.yy0ex0bar)));
        nSYM.dghxu     = reshape(double(subs(SYM.dghxu)), M_.endo_nbr, M_.nspred*M_.exo_nbr, length(SYM.params));
        nSYM.ghuu      = double(subs(subs(SYM.ghuu, SYM.yy0ex0, SYM.yy0ex0bar)));
        nSYM.dghuu     = reshape(double(subs(SYM.dghuu)), M_.endo_nbr, M_.exo_nbr^2, length(SYM.params));
        nSYM.ghs2      = double(subs(subs(SYM.ghs2, SYM.yy0ex0, SYM.yy0ex0bar)));
        nSYM.dghs2     = reshape(double(subs(SYM.dghs2)), M_.endo_nbr, length(SYM.params));
        nSYM.ghxxx     = double(subs(subs(SYM.ghxxx, SYM.yy0ex0, SYM.yy0ex0bar)));
        nSYM.dghxxx    = reshape(double(subs(SYM.dghxxx)), M_.endo_nbr, M_.nspred^3, length(SYM.params));
        nSYM.ghxxu     = double(subs(subs(SYM.ghxxu, SYM.yy0ex0, SYM.yy0ex0bar)));
        nSYM.dghxxu    = reshape(double(subs(SYM.dghxxu)), M_.endo_nbr, M_.nspred^2*M_.exo_nbr, length(SYM.params));
        nSYM.ghxuu     = double(subs(subs(SYM.ghxuu, SYM.yy0ex0, SYM.yy0ex0bar)));
        nSYM.dghxuu    = reshape(double(subs(SYM.dghxuu)), M_.endo_nbr, M_.nspred*M_.exo_nbr^2, length(SYM.params));
        nSYM.ghuuu     = double(subs(subs(SYM.ghuuu, SYM.yy0ex0, SYM.yy0ex0bar)));
        nSYM.dghuuu    = reshape(double(subs(SYM.dghuuu)), M_.endo_nbr, M_.exo_nbr^3, length(SYM.params));
        nSYM.ghxss     = double(subs(subs(SYM.ghxss, SYM.yy0ex0, SYM.yy0ex0bar)));
        nSYM.dghxss    = reshape(double(subs(SYM.dghxss)), M_.endo_nbr, M_.nspred, length(SYM.params));
        nSYM.ghuss     = double(subs(subs(SYM.ghuss, SYM.yy0ex0, SYM.yy0ex0bar)));
        nSYM.dghuss    = reshape(double(subs(SYM.dghuss)), M_.endo_nbr, M_.exo_nbr, length(SYM.params));
        if jj == 1
            nSYMprior = nSYM;
            save('nBrockMirmanSYM.mat','nSYMprior')
        elseif jj==2
            nSYMcalib = nSYM;
            save('nBrockMirmanSYM.mat','nSYMcalib','-append')        
        end
    end
@#endif


clc;
tol_vars.Yss                = 1e-14;
tol_vars.Sigma_e            = 1e-15;
tol_vars.Correlation_matrix = 1e-16;
tol_vars.g1                 = 1e-13;
tol_vars.ghx                = 1e-13;
tol_vars.ghu                = 1e-13;
@#if ORDER > 1
tol_vars.g2                 = 1e-12;
tol_vars.ghxx               = 1e-12;
tol_vars.ghxu               = 1e-12;
tol_vars.ghuu               = 1e-12;
tol_vars.ghs2               = 1e-12;
@#endif
@#if ORDER > 2
tol_vars.g3                 = 1e-11;
tol_vars.ghxxx              = 1e-11;
tol_vars.ghxxu              = 1e-11;
tol_vars.ghxuu              = 1e-11;
tol_vars.ghuuu              = 1e-11;
tol_vars.ghxss              = 1e-11;
tol_vars.ghuss              = 1e-11;
@#endif

tol_dvars.dYss =                [1e-9  1e-9 1e-14 1e-14];
tol_dvars.dSigma_e =            [1e-9  1e-15 1e-14 1e-14];
tol_dvars.dCorrelation_matrix = [1e-9  1e-15 1e-14 1e-14];
tol_dvars.dg1 =                 [1e-6  1e-6  1e-13 1e-13];
tol_dvars.dghx =                [1e-8  1e-8  1e-13 1e-13];
tol_dvars.dghu =                [1e-8  1e-8  1e-13 1e-13];
@#if ORDER > 1
tol_dvars.dg2 =                 [1e-5  1e-5  1e-11 1e-11];
tol_dvars.dghxx =               [1e-6  1e-6  1e-12 1e-12];
tol_dvars.dghxu =               [1e-6  1e-6  1e-12 1e-12];
tol_dvars.dghuu =               [1e-6  1e-6  1e-12 1e-12];
tol_dvars.dghs2 =               [1e-6  1e-6  1e-12 1e-12];
@#endif
@#if ORDER > 2
tol_dvars.dg3 =                 [1e-3  1e-3  1e-9 1e-9];
tol_dvars.dghxxx =              [1e-5  1e-5  1e-11 1e-11];
tol_dvars.dghxxu =              [1e-5  1e-5  1e-11 1e-11];
tol_dvars.dghxuu =              [1e-5  1e-5  1e-11 1e-11];
tol_dvars.dghuuu =              [1e-5  1e-5  1e-11 1e-11];
tol_dvars.dghxss =              [1e-5  1e-5  1e-11 1e-11];
tol_dvars.dghuss =              [1e-5  1e-5  1e-11 1e-11];
@#endif
tol_dvars.d2KalmanA =           [1e-3 1e-3 1e-13 1e-13];
tol_dvars.d2Om =                [1e-3 1e-3 1e-13 1e-13];
tol_dvars.d2Yss =               [1e-3 1e-3 1e-13 1e-13];

options_.dynatol.x = eps.^(1/3); %set numerical differentiation step in fjaco.m

for jj = 1:2
    lst_vars = {'Yss', 'Sigma_e', 'Correlation_matrix','g1','ghx','ghu'};
                 @#if ORDER > 1
    lst_vars =   [lst_vars, 'g2','ghxx','ghxu','ghuu','ghs2'];
                 @#endif
                 @#if ORDER > 2
    lst_vars =   [lst_vars, 'g3','ghxxx','ghxxu','ghxuu','ghuuu','ghxss','ghuss'];
                 @#endif
    lst_dvars = {'dYss','dSigma_e','dg1','dghx','dghu'};
                 @#if ORDER > 1 
    lst_dvars = [lst_dvars, 'dg2','dghxx','dghxu','dghuu','dghs2'];
                 @#endif
                 @#if ORDER > 2
    lst_dvars = [lst_dvars, 'dg3','dghxxx','dghxxu','dghxuu','dghuuu','dghxss','dghuss'];
                 @#endif
    load('nBrockMirmanSYM.mat');
     if jj==1
        strparamset = 'PRIOR';
        nSYM = nSYMprior;
        xparam_prior = set_prior(estim_params_,M_,options_);
        M_ = set_all_parameters(xparam_prior,estim_params_,M_);
    elseif jj==2
        strparamset = 'CALIBRATION';
        nSYM = nSYMcalib;
        xparam1_calib = [];
        for j = 1:length(indpstderr)
            xparam1_calib = [xparam1_calib; sqrt(calib_Sigma_e(j,j))];            
        end
        for j = 1:size(indpcorr,1)
            xparam1_calib = [xparam1_calib; calib_Sigma_e(indpcorr(j,1),indpcorr(j,2))/( sqrt(calib_Sigma_e(indpcorr(j,1),indpcorr(j,1))) * sqrt(calib_Sigma_e(indpcorr(j,2),indpcorr(j,2))) )];
        end
        xparam1_calib = [xparam1_calib; calib_params(indpmodel)];
        M_ = set_all_parameters(xparam1_calib,estim_params_,M_);
    end
    [~,info,M_,options_,oo_] = resol(0,M_, options_, oo_);
    %For convenience we save the objects to compare into oo_.dr.
    oo_.dr.Yss  = oo_.dr.ys(oo_.dr.order_var);
    oo_.dr.Sigma_e = M_.Sigma_e;
    oo_.dr.Correlation_matrix = M_.Correlation_matrix;
    ex0 = oo_.exo_steady_state';
    [~, oo_.dr.g1, oo_.dr.g2, oo_.dr.g3] = feval([M_.fname,'.dynamic'], oo_.dr.ys(I), oo_.exo_steady_state', M_.params, oo_.dr.ys, 1);
    oo_.dr.g3 = unfold_g3(oo_.dr.g3, length(oo_.dr.ys(I))+length(oo_.exo_steady_state')); %add symmetric elements to g3

    fprintf('***** %s: SOME COMMON OBJECTS *****\n', strparamset)
    for id_var = 1:size(lst_vars,2)
        dx = norm( nSYM.(sprintf('%s',lst_vars{id_var})) - oo_.dr.(sprintf('%s',lst_vars{id_var})), Inf);
        fprintf('Max absolute deviation for %s: %e\n', lst_vars{id_var}, dx);
        if dx > tol_vars.(sprintf('%s',lst_vars{id_var}))
            error('Something wrong in steady state computation, solution algorithm or preprocessor')
        end        
    end


    for d2flag = [0 1];
        if d2flag
            lst_dvars = [lst_dvars {'d2KalmanA', 'd2Om', 'd2Yss'}];
        end    
        KRONFLAG = [-1 -2 0 1];    
        for id_kronflag = 1:length(KRONFLAG)
            fprintf('***** %s: d2flag=%d and kronflag=%d *****\n',strparamset, d2flag,KRONFLAG(id_kronflag))
            options_.analytic_derivation_mode = KRONFLAG(id_kronflag);        
            DERIVS = get_perturbation_params_derivs(M_, options_, estim_params_, oo_, indpmodel, indpstderr, indpcorr, d2flag);
            for id_var = 1:size(lst_dvars,2)
                dx = norm( vec(nSYM.(sprintf('%s',lst_dvars{id_var}))) - vec(DERIVS.(sprintf('%s',lst_dvars{id_var}))), Inf);
                fprintf('Max absolute deviation for %s: %e\n', lst_dvars{id_var}, dx);
                if dx > tol_dvars.(sprintf('%s',lst_dvars{id_var}))(id_kronflag)
                    error('Something wrong in get_perturbation_params_derivs.m')
                end
            end
        end
    end
end
