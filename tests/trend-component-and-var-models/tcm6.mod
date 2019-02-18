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
       res_U2_G_EYER
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

[name='U2_Q_YED']
diff(U2_Q_YED) =   u2_q_yed_ecm_u2_q_yed_L1 * (U2_Q_YED(-1) - U2_EHIC(-1))
                 + u2_q_yed_ecm_u2_stn_L1   * (U2_STN(-1)   - U2_ESTN(-1))
                 + u2_q_yed_ecm_u2_g_yer_L1 * (U2_G_YER(-1)- U2_G_EYER(-1))
                 + u2_q_yed_u2_q_yed_L1     * diff(U2_Q_YED(-1))
                 + u2_q_yed_u2_g_yer_L1     * diff(U2_G_YER(-1))
                 + u2_q_yed_u2_stn_L1       * diff(U2_STN(-1))
                 + res_U2_Q_YED                                           ;

[name='U2_G_YER']
diff(U2_G_YER) =   u2_g_yer_ecm_u2_q_yed_L1 * (U2_Q_YED(-1) - U2_EHIC(-1))
                 + u2_g_yer_ecm_u2_stn_L1   * (U2_STN(-1)   - U2_ESTN(-1))
                 + u2_g_yer_ecm_u2_g_yer_L1 * (U2_G_YER(-1) - U2_G_EYER(-1))
                 + u2_g_yer_u2_q_yed_L1     * diff(U2_Q_YED(-1))
                 + u2_g_yer_u2_g_yer_L1     * diff(U2_G_YER(-1))
                 + u2_g_yer_u2_stn_L1       * diff(U2_STN(-1))
                 + res_U2_G_YER                                           ;

[name='U2_STN']
diff(U2_STN)   =   u2_stn_ecm_u2_q_yed_L1   * (U2_Q_YED(-1) - U2_EHIC(-1))
                 + u2_stn_ecm_u2_stn_L1     * (U2_STN(-1)   - U2_ESTN(-1))
                 + u2_stn_ecm_u2_g_yer_L1   * (U2_G_YER(-1) - U2_G_EYER(-1))
                 + u2_stn_u2_q_yed_L1       * diff(U2_Q_YED(-1))
                 + u2_stn_u2_g_yer_L1       * diff(U2_G_YER(-1))
                 + res_U2_STN                                             ;

[name='U2_ESTN']
U2_ESTN        =  U2_ESTN(-1) + res_U2_ESTN                               ;

[name='U2_EHIC']
U2_EHIC        =  U2_EHIC(-1) + res_U2_EHIC                               ;

[name='U2_G_EYER']
U2_G_EYER        =  U2_G_EYER(-1) + res_U2_G_EYER                               ;

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

[name='zpac']
diff(log(U2_HH_COR)) = ecm_pac*(log(U2_HH_COR(-1))-log(U2_HH_OCOR(-1))) +
                  u2_hh_cor_pac_u2_hh_cor_L1*diff(log(U2_HH_COR(-1))) +
                  pac_expectation(pacman)   +
                  res_ez;
end;

// Random calibration for the TREND_COMPONENT_MODEL's parameters.
M_.params(1:25) = randn(25, 1);
params = M_.params;

trend_component_model(model_name=toto, eqtags=['U2_Q_YED', 'U2_G_YER', 'U2_STN', 'U2_EHIC', 'U2_G_EYER', 'U2_ESTN', 'U2_HH_OCOR', 'U2_H_Q_YER400'], targets=['U2_EHIC', 'U2_G_EYER', 'U2_ESTN', 'U2_H_Q_YER400']);
pac_model(auxiliary_model_name=toto, discount=beta, model_name=pacman, growth = U2_H_Q_YER400);
pac.initialize('pacman');
C0 = oo_.trend_component.toto.CompanionMatrix;

trend_component_model(model_name=titi, eqtags=['U2_Q_YED', 'U2_G_YER', 'U2_STN', 'U2_EHIC', 'U2_G_EYER', 'U2_ESTN', 'U2_HH_OCOR', 'U2_H_Q_YER400'], targets=['U2_G_EYER', 'U2_H_Q_YER400', 'U2_ESTN', 'U2_EHIC']);
pac_model(auxiliary_model_name=titi, discount=beta, model_name=pacman1, growth = U2_H_Q_YER400);
pac.initialize('pacman1');
C1 = oo_.trend_component.titi.CompanionMatrix;

save('tcm6_data.mat', 'C0', 'C1', 'params');

[AR0, A00, A00star] = feval([M_.fname '.trend_component_ar_a0'], 'toto', M_.params);
[AR1, A01, A01star] = feval([M_.fname '.trend_component_ar_a0'], 'titi', M_.params);

if sum(sum(abs(AR0-AR1))) ~= 0
    error('Problem with AR matrices')
end

if sum(sum(abs(A00-A01))) ~= 0
    error('Problem with A0 matrices')
end

if sum(sum(abs(A00star-A01star(:,[4,1,3,2])))) ~= 0
    error('Problem with A0star matrices')
end

ar = zeros(4, 4, 1);
ar(1,:) = [M_.params(strcmp(M_.param_names, 'u2_q_yed_u2_q_yed_L1')),   M_.params(strcmp(M_.param_names, 'u2_q_yed_u2_g_yer_L1')),   M_.params(strcmp(M_.param_names, 'u2_q_yed_u2_stn_L1')),   0];
ar(2,:) = [M_.params(strcmp(M_.param_names, 'u2_g_yer_u2_q_yed_L1')),   M_.params(strcmp(M_.param_names, 'u2_g_yer_u2_g_yer_L1')),   M_.params(strcmp(M_.param_names, 'u2_g_yer_u2_stn_L1')),   0];
ar(3,:) = [M_.params(strcmp(M_.param_names, 'u2_stn_u2_q_yed_L1')),     M_.params(strcmp(M_.param_names, 'u2_stn_u2_g_yer_L1')),     0,                                                         0];
ar(4,:) = [M_.params(strcmp(M_.param_names, 'u2_hh_ocor_u2_q_yed_L1')), M_.params(strcmp(M_.param_names, 'u2_hh_ocor_u2_g_yer_L1')), M_.params(strcmp(M_.param_names, 'u2_hh_ocor_u2_stn_L1')), M_.params(strcmp(M_.param_names, 'u2_hh_ocor_u2_hh_ocor_L1'))];

if sum(sum(abs(AR0-ar))) ~= 0
    error('Problem with AR matrix')
end

a0 = zeros(4, 4, 1);
a0(1,:) = [M_.params(strcmp(M_.param_names, 'u2_q_yed_ecm_u2_q_yed_L1')),   M_.params(strcmp(M_.param_names, 'u2_q_yed_ecm_u2_g_yer_L1')),   M_.params(strcmp(M_.param_names, 'u2_q_yed_ecm_u2_stn_L1')),   0];
a0(2,:) = [M_.params(strcmp(M_.param_names, 'u2_g_yer_ecm_u2_q_yed_L1')),   M_.params(strcmp(M_.param_names, 'u2_g_yer_ecm_u2_g_yer_L1')),   M_.params(strcmp(M_.param_names, 'u2_g_yer_ecm_u2_stn_L1')),   0];
a0(3,:) = [M_.params(strcmp(M_.param_names, 'u2_stn_ecm_u2_q_yed_L1')),     M_.params(strcmp(M_.param_names, 'u2_stn_ecm_u2_g_yer_L1')),     M_.params(strcmp(M_.param_names, 'u2_stn_ecm_u2_stn_L1')),     0];
a0(4,:) = [M_.params(strcmp(M_.param_names, 'u2_hh_ocor_ecm_u2_q_yed_L1')), M_.params(strcmp(M_.param_names, 'u2_hh_ocor_ecm_u2_g_yer_L1')), M_.params(strcmp(M_.param_names, 'u2_hh_ocor_ecm_u2_stn_L1')), M_.params(strcmp(M_.param_names, 'u2_hh_ocor_ecm_u2_hh_ocor_L1'))];

if sum(sum(abs(A00star-a0))) ~= 0
    error('Problem with A0star matrix')
end
