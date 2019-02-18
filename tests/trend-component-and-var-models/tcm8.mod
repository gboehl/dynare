// --+ options: stochastic,json=compute +--

var U2_Q_YED
    U2_G_YER
    U2_STN
    U2_EHIC
    U2_ESTN
    U2_G_EYER
    U2_HH_OCOR
    U2_HH_COR
    U2_H_Q_YER400 ;

varexo res_U2_Q_YED
       res_U2_G_YER
       res_U2_STN
       res_U2_EHIC
       res_U2_ESTN
       res_U2_HH_OCOR
       res_U2_H_Q_YER
       res_ez ;

parameters u2_q_yed_ecm_u2_q_yed_L1
           u2_q_yed_ecm_u2_stn_L1
           u2_q_yed_u2_q_yed_L1
           u2_q_yed_u2_g_yer_L1
           u2_q_yed_u2_stn_L1
           u2_g_yer_ecm_u2_q_yed_L1
           u2_g_yer_ecm_u2_stn_L1
           u2_g_yer_u2_q_yed_L1
           u2_g_yer_u2_g_yer_L1
           u2_g_yer_u2_stn_L1
           u2_stn_ecm_u2_q_yed_L1
           u2_stn_ecm_u2_stn_L1
           u2_stn_u2_q_yed_L1
           u2_stn_u2_g_yer_L1
           u2_q_yed_ecm_u2_g_yer_L1
           u2_g_yer_ecm_u2_g_yer_L1
           u2_stn_ecm_u2_g_yer_L1
           u2_hh_ocor_ecm_u2_q_yed_L1
           u2_hh_ocor_ecm_u2_stn_L1
           u2_hh_ocor_ecm_u2_g_yer_L1
           u2_hh_ocor_u2_q_yed_L1
           u2_hh_ocor_u2_g_yer_L1
           u2_hh_ocor_u2_stn_L1
           u2_hh_ocor_ecm_u2_hh_ocor_L1
           u2_hh_ocor_u2_hh_ocor_L1
           beta
           ecm_pac
           u2_hh_cor_pac_u2_hh_cor_L1 ;


beta = 0.98 ;
ecm_pac = 0.2;
u2_hh_cor_pac_u2_hh_cor_L1 = 0.4;

model;

[name='U2_G_YER']
diff(U2_G_YER) =   u2_g_yer_ecm_u2_q_yed_L1 * (U2_Q_YED(-1) - U2_EHIC(-1))
                 + u2_g_yer_ecm_u2_stn_L1   * (U2_STN(-1)   - U2_ESTN(-1))
                 + u2_g_yer_ecm_u2_g_yer_L1 * (U2_G_YER(-1) - U2_G_EYER(-1))
                 + u2_g_yer_u2_q_yed_L1     * diff(U2_Q_YED(-1))
                 + u2_g_yer_u2_g_yer_L1     * diff(U2_G_YER(-1))
                 + u2_g_yer_u2_stn_L1       * diff(U2_STN(-1))
                 + res_U2_G_YER                                           ;


[name='U2_Q_YED']
diff(U2_Q_YED) =   u2_q_yed_ecm_u2_q_yed_L1 * (U2_Q_YED(-1) - U2_EHIC(-1))
                 + u2_q_yed_ecm_u2_stn_L1   * (U2_STN(-1)   - U2_ESTN(-1))
                 + u2_q_yed_ecm_u2_g_yer_L1 * (U2_G_YER(-1)- U2_G_EYER(-1))
                 + u2_q_yed_u2_q_yed_L1     * diff(U2_Q_YED(-1))
                 + u2_q_yed_u2_g_yer_L1     * diff(U2_G_YER(-1))
                 + u2_q_yed_u2_stn_L1       * diff(U2_STN(-1))
                 + res_U2_Q_YED                                           ;


[name='U2_ESTN']
U2_ESTN        =  U2_ESTN(-1) + res_U2_ESTN                               ;

[name='U2_EHIC']
U2_EHIC        =  U2_EHIC(-1) + res_U2_EHIC                               ;

[name='U2_G_EYER']
U2_G_EYER        =  0;

[name='U2_HH_OCOR']
diff(diff(log(U2_HH_OCOR))) =  u2_hh_ocor_ecm_u2_q_yed_L1   * (U2_Q_YED(-1) - U2_EHIC(-1))
                 + u2_hh_ocor_ecm_u2_stn_L1     * (U2_STN(-1)   - U2_ESTN(-1))
                 + u2_hh_ocor_ecm_u2_g_yer_L1   * (U2_G_YER(-1) - U2_G_EYER(-1))
                 + u2_hh_ocor_ecm_u2_hh_ocor_L1    * (diff(log(U2_HH_OCOR(-1))) - U2_H_Q_YER400(-1))
                 + u2_hh_ocor_u2_q_yed_L1       * diff(U2_Q_YED(-1))
                 + u2_hh_ocor_u2_g_yer_L1       * diff(U2_G_YER(-1))
                 + u2_hh_ocor_u2_stn_L1         * diff(U2_STN(-1))
                 + u2_hh_ocor_u2_hh_ocor_L1        * diff(diff(log(U2_HH_OCOR(-1))))
                 + res_U2_HH_OCOR                                        ;

[name='U2_H_Q_YER400']
U2_H_Q_YER400 = U2_H_Q_YER400(-1) + res_U2_H_Q_YER;

[name='U2_STN']
diff(U2_STN)   =   u2_stn_ecm_u2_q_yed_L1   * (U2_Q_YED(-1) - U2_EHIC(-1))
                 + u2_stn_ecm_u2_stn_L1     * (U2_STN(-1)   - U2_ESTN(-1))
                 + u2_stn_ecm_u2_g_yer_L1   * (U2_G_YER(-1) - U2_G_EYER(-1))
                 + u2_stn_u2_q_yed_L1       * diff(U2_Q_YED(-1))
                 + u2_stn_u2_g_yer_L1       * diff(U2_G_YER(-1))
                 + res_U2_STN                                             ;

[name='zpac']
diff(log(U2_HH_COR)) = ecm_pac*(log(U2_HH_COR(-1))-log(U2_HH_OCOR(-1))) +
                  u2_hh_cor_pac_u2_hh_cor_L1*diff(log(U2_HH_COR(-1))) +
                  pac_expectation(pacman)   +
                  res_ez;
end;

// Use the same calibration as in tcm6.mod.
tcm6 = load('tcm6_data.mat');
M_.params = tcm6.params;

trend_component_model(model_name=toto, eqtags=['U2_Q_YED', 'U2_G_YER', 'U2_STN', 'U2_EHIC', 'U2_ESTN', 'U2_HH_OCOR', 'U2_H_Q_YER400'], targets=['U2_EHIC', 'U2_ESTN', 'U2_H_Q_YER400']);
trend_component_model(model_name=titi, eqtags=['U2_Q_YED', 'U2_G_YER', 'U2_STN', 'U2_EHIC', 'U2_ESTN', 'U2_HH_OCOR', 'U2_H_Q_YER400'], targets=['U2_H_Q_YER400', 'U2_ESTN', 'U2_EHIC']);

% toto
AR  = zeros(4,4,1);
A0  = zeros(4,4,1);
A0s = zeros(4,3,1);

AR(1,:) = [M_.params(strcmp(M_.param_names, 'u2_q_yed_u2_q_yed_L1')), M_.params(strcmp(M_.param_names, 'u2_q_yed_u2_g_yer_L1')), M_.params(strcmp(M_.param_names, 'u2_q_yed_u2_stn_L1')), 0];
AR(2,:) = [M_.params(strcmp(M_.param_names, 'u2_g_yer_u2_q_yed_L1')), M_.params(strcmp(M_.param_names, 'u2_g_yer_u2_g_yer_L1')), M_.params(strcmp(M_.param_names, 'u2_g_yer_u2_stn_L1')), 0];
AR(3,:) = [M_.params(strcmp(M_.param_names, 'u2_stn_u2_q_yed_L1')), M_.params(strcmp(M_.param_names, 'u2_stn_u2_g_yer_L1')), 0, 0];
AR(4,:) = [M_.params(strcmp(M_.param_names, 'u2_hh_ocor_u2_q_yed_L1')), M_.params(strcmp(M_.param_names, 'u2_hh_ocor_u2_g_yer_L1')), M_.params(strcmp(M_.param_names, 'u2_hh_ocor_u2_stn_L1')), M_.params(strcmp(M_.param_names, 'u2_hh_ocor_u2_hh_ocor_L1'))];

A0(1,:) = [M_.params(strcmp(M_.param_names, 'u2_q_yed_ecm_u2_q_yed_L1')), M_.params(strcmp(M_.param_names, 'u2_q_yed_ecm_u2_g_yer_L1')), M_.params(strcmp(M_.param_names, 'u2_q_yed_ecm_u2_stn_L1')), 0];
A0(2,:) = [M_.params(strcmp(M_.param_names, 'u2_g_yer_ecm_u2_q_yed_L1')), M_.params(strcmp(M_.param_names, 'u2_g_yer_ecm_u2_g_yer_L1')), M_.params(strcmp(M_.param_names, 'u2_g_yer_ecm_u2_stn_L1')), 0];
A0(3,:) = [M_.params(strcmp(M_.param_names, 'u2_stn_ecm_u2_q_yed_L1')), M_.params(strcmp(M_.param_names, 'u2_stn_ecm_u2_g_yer_L1')), M_.params(strcmp(M_.param_names, 'u2_stn_ecm_u2_stn_L1')), 0];
A0(4,:) = [M_.params(strcmp(M_.param_names, 'u2_hh_ocor_ecm_u2_q_yed_L1')), M_.params(strcmp(M_.param_names, 'u2_hh_ocor_ecm_u2_g_yer_L1')), M_.params(strcmp(M_.param_names, 'u2_hh_ocor_ecm_u2_stn_L1')), M_.params(strcmp(M_.param_names, 'u2_hh_ocor_ecm_u2_hh_ocor_L1'))];

A0s(1,:) = [M_.params(strcmp(M_.param_names, 'u2_q_yed_ecm_u2_q_yed_L1')), M_.params(strcmp(M_.param_names, 'u2_q_yed_ecm_u2_stn_L1')), 0];
A0s(2,:) = [M_.params(strcmp(M_.param_names, 'u2_g_yer_ecm_u2_q_yed_L1')), M_.params(strcmp(M_.param_names, 'u2_g_yer_ecm_u2_stn_L1')), 0];
A0s(3,:) = [M_.params(strcmp(M_.param_names, 'u2_stn_ecm_u2_q_yed_L1')), M_.params(strcmp(M_.param_names, 'u2_stn_ecm_u2_stn_L1')), 0];
A0s(4,:) = [M_.params(strcmp(M_.param_names, 'u2_hh_ocor_ecm_u2_q_yed_L1')), M_.params(strcmp(M_.param_names, 'u2_hh_ocor_ecm_u2_stn_L1')), M_.params(strcmp(M_.param_names, 'u2_hh_ocor_ecm_u2_hh_ocor_L1'))];

[ARt, A0t, A0st] = tcm8.trend_component_ar_a0('toto', M_.params);

if sum(sum(abs(AR-ARt))) ~= 0
    error('Problem with AR')
end

if sum(sum(abs(A0-A0t))) ~= 0
    error('Problem with A0')
end

if sum(sum(abs(A0s-A0st))) ~= 0
    error('Problem with A0star')
end

[ARtt, A0tt, A0stt] = tcm8.trend_component_ar_a0('titi', M_.params);

if sum(sum(abs(AR-ARtt))) ~= 0
    error('Problem with AR')
end

if sum(sum(abs(A0-A0tt))) ~= 0
    error('Problem with A0')
end

if sum(sum(abs(A0s-A0stt(:,[3,2,1])))) ~= 0
    error('Problem with A0star')
end

pac_model(auxiliary_model_name=toto, discount=beta, model_name=pacman, growth = U2_H_Q_YER400);
pac.initialize('pacman');
C0 = oo_.trend_component.toto.CompanionMatrix;

pac_model(auxiliary_model_name=titi, discount=beta, model_name=pacman1, growth = U2_H_Q_YER400);
pac.initialize('pacman1');
C1 = oo_.trend_component.titi.CompanionMatrix;

if any(abs(C0(:)-tcm6.C0(:)))>1e-12
   error('Companion matrix is not independent of the ordering of the equations.')
end

if any(abs(C1(:)-tcm6.C1(:)))>1e-12
   error('Companion matrix is not independent of the ordering of the equations.')
end