// --+ options: json=compute +--

/* REMARK
** ------
**
** You need to have the first line on top of the mod file. The options defined on this line are passed
** to the dynare command (you can add other options, separated by spaces or commas). The option defined
** here is mandatory for the decomposition. It forces Dynare to output another representation of the
** model in JSON file (additionaly to the matlab files) which is used here to manipulate the equations.
*/

var
U2_Q_YED
U2_G_YER
U2_STN
U2_ESTN
U2_EHIC
DE_Q_YED
DE_G_YER
DE_EHIC

;

varexo
res_U2_Q_YED
res_U2_G_YER
res_U2_STN
res_U2_ESTN
res_U2_EHIC
res_DE_Q_YED
res_DE_G_YER
res_DE_EHIC
;

parameters
u2_q_yed_ecm_u2_q_yed_L1
u2_q_yed_ecm_u2_stn_L1
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
u2_estn_u2_estn_L1
u2_ehic_u2_ehic_L1

de_q_yed_ecm_de_q_yed_L1
de_q_yed_ecm_u2_stn_L1
de_q_yed_de_g_yer_L1
de_q_yed_u2_stn_L1
de_g_yer_ecm_de_q_yed_L1
de_g_yer_ecm_u2_stn_L1
de_g_yer_de_q_yed_L1
de_g_yer_de_g_yer_L1
de_g_yer_u2_stn_L1
de_ehic_de_ehic_L1


;

u2_q_yed_ecm_u2_q_yed_L1  = -0.82237516589315   ;
u2_q_yed_ecm_u2_stn_L1    = -0.323715338568976  ;
u2_q_yed_u2_g_yer_L1      =  0.0401361895021084 ;
u2_q_yed_u2_stn_L1        =  0.058397703958446  ;
u2_g_yer_ecm_u2_q_yed_L1  =  0.0189896046977421 ;
u2_g_yer_ecm_u2_stn_L1    = -0.109597659887432  ;
u2_g_yer_u2_q_yed_L1      =  0.0037667967632025 ;
u2_g_yer_u2_g_yer_L1      =  0.480506381923644  ;
u2_g_yer_u2_stn_L1        = -0.0722359286123494 ;
u2_stn_ecm_u2_q_yed_L1    = -0.0438500662608356 ;
u2_stn_ecm_u2_stn_L1      = -0.153283917138772  ;
u2_stn_u2_q_yed_L1        =  0.0328744983772825 ;
u2_stn_u2_g_yer_L1        =  0.292121949736756  ;
u2_estn_u2_estn_L1        =  1                  ;
u2_ehic_u2_ehic_L1        =  1                  ;

de_q_yed_ecm_de_q_yed_L1  = -0.822375165893149  ;
de_q_yed_ecm_u2_stn_L1    = -0.323715338568977  ;
de_q_yed_de_g_yer_L1      =  0.0401361895021082 ;
de_q_yed_u2_stn_L1        =  0.0583977039584461 ;
de_g_yer_ecm_de_q_yed_L1  =  0.0189896046977422 ;
de_g_yer_ecm_u2_stn_L1    = -0.109597659887433  ;
de_g_yer_de_q_yed_L1      =  0.00376679676320256;
de_g_yer_de_g_yer_L1      =  0.480506381923643  ;
de_g_yer_u2_stn_L1        = -0.0722359286123494 ;
de_ehic_de_ehic_L1        =  1                  ;


model(linear);
[name = 'eq1']
diff(U2_Q_YED) =   u2_q_yed_ecm_u2_q_yed_L1 * (U2_Q_YED(-1) - U2_EHIC(-1))
                 + u2_q_yed_ecm_u2_stn_L1   * (U2_STN(-1)   - U2_ESTN(-1))
                 + u2_q_yed_u2_g_yer_L1     * diff(U2_G_YER(-1))
                 + u2_q_yed_u2_stn_L1       * diff(U2_STN(-1))
                 + res_U2_Q_YED                                           ;
[name = 'eq2']
diff(U2_G_YER) =   u2_g_yer_ecm_u2_q_yed_L1 * (U2_Q_YED(-1) - U2_EHIC(-1))
                 + u2_g_yer_ecm_u2_stn_L1   * (U2_STN(-1)   - U2_ESTN(-1))
                 + u2_g_yer_u2_q_yed_L1     * diff(U2_Q_YED(-1))
                 + u2_g_yer_u2_g_yer_L1     * diff(U2_G_YER(-1))
                 + u2_g_yer_u2_stn_L1       * diff(U2_STN(-1))
                 + res_U2_G_YER                                           ;
[name = 'eq3']
diff(U2_STN)   =   u2_stn_ecm_u2_q_yed_L1   * (U2_Q_YED(-1) - U2_EHIC(-1))
                 + u2_stn_ecm_u2_stn_L1     * (U2_STN(-1)   - U2_ESTN(-1))
                 + u2_stn_u2_q_yed_L1       * diff(U2_Q_YED(-1))
                 + u2_stn_u2_g_yer_L1       * diff(U2_G_YER(-1))
                 + res_U2_STN                                             ;
[name = 'eq4']
U2_ESTN        =   u2_estn_u2_estn_L1       * U2_ESTN(-1)
                 + res_U2_ESTN                                            ;
[name = 'eq5']
U2_EHIC        =   u2_ehic_u2_ehic_L1       * U2_EHIC(-1)
                 + res_U2_EHIC                                            ;
[name = 'eq6']
diff(DE_Q_YED) =   de_q_yed_ecm_de_q_yed_L1 * (DE_Q_YED(-1) - DE_EHIC(-1))
                 + de_q_yed_ecm_u2_stn_L1   * (U2_STN(-1)   - U2_ESTN(-1))
                 + de_q_yed_de_g_yer_L1     * diff(DE_G_YER(-1))
                 + de_q_yed_u2_stn_L1       * diff(U2_STN(-1))
                 + res_DE_Q_YED                                           ;
[name = 'eq7']
diff(DE_G_YER) =   de_g_yer_ecm_de_q_yed_L1 * (DE_Q_YED(-1) - DE_EHIC(-1))
                 + de_g_yer_ecm_u2_stn_L1   * (U2_STN(-1)   - U2_ESTN(-1))
                 + de_g_yer_de_q_yed_L1     * diff(DE_Q_YED(-1))
                 + de_g_yer_de_g_yer_L1     * diff(DE_G_YER(-1))
                 + de_g_yer_u2_stn_L1       * diff(U2_STN(-1))
                 + res_DE_G_YER                                           ;
[name = 'eq8']
DE_EHIC        =   de_ehic_de_ehic_L1       * DE_EHIC(-1)
                 + res_DE_EHIC                                            ;



end;

shocks;
var res_U2_Q_YED = 0.005;
var res_U2_G_YER = 0.005;
var res_U2_STN = 0.005;
var res_U2_ESTN = 0.005;
var res_U2_EHIC = 0.005;
var res_DE_Q_YED = 0.005;
var res_DE_G_YER = 0.005;
var res_DE_EHIC = 0.005;
end;

NSIMS = 1000;

calibrated_values = M_.params;
verbatim;
Sigma_e = M_.Sigma_e;
end;

options_.bnlms.set_dynare_seed_to_default = false;

nparampool = length(M_.params);
BETA = zeros(NSIMS, nparampool);
for i=1:NSIMS
    firstobs = rand(3, length(M_.endo_names));
    M_.params = calibrated_values;
    M_.Sigma_e = Sigma_e; 
    simdata = simul_backward_model(dseries(firstobs, dates('1995Q1'), M_.endo_names), 10000);
    simdata = simdata(simdata.dates(5001:6000));
    names=regexp(simdata.name, 'res\w*');
    idxs = [];
    for j=1:length(names)
        if isempty(names{j})
            idxs = [idxs j];
        end
    end
    pooled_fgls(simdata{idxs}, ...
        {'de','u2'}, ...
        {'*_q_yed_ecm_*_q_yed_L1', ...
        '*_q_yed_ecm_u2_stn_L1', ...
        '*_q_yed_*_g_yer_L1', ...
        '*_q_yed_u2_stn_L1', ...
        '*_g_yer_ecm_*_q_yed_L1', ...
        '*_g_yer_ecm_u2_stn_L1', ...
        '*_g_yer_*_q_yed_L1', ...
        '*_g_yer_*_g_yer_L1', ...
        '*_g_yer_u2_stn_L1', ...
        '*_ehic_*_ehic_L1'});
    BETA(i, :) = M_.params';
end

tmp = mean(BETA)' - calibrated_values;
good = [-0.117075725840299
  -0.415772666038150
  -0.164200084092597
  -0.590183263021335
   0.145687952199179
   0.095087188729817
  -0.109667429232965
  -0.281691072046120
  -0.386631489016812
   0.147073439164064
  -0.582630982519539
  -0.768471596296734
  -0.280530746352647
  -0.622354568301842
  -0.213391752393322
  -0.117075725841298
  -0.415772666037151
  -0.164200084092396
  -0.590183263021439
   0.145687952199079
   0.095087188730816
  -0.109667429233025
  -0.281691072045120
  -0.386631489016812
  -0.213391752393322]*1.0e-03;
if sum(abs(tmp-good)) > 1e-14
    error(['sum of tmp - good was: ' num2str(sum(abs(tmp-good)))]);
end

for i=1:nparampool
    figure
    hold on
    title(strrep(M_.param_names(i,:), '_', '\_'));
    histogram(BETA(:,i),50);
    line([calibrated_values(i) calibrated_values(i)], [0 NSIMS/10], 'LineWidth', 2, 'Color', 'r');
    hold off
end
