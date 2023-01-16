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
DE_EHIC        = DE_EHIC(-1) + res_DE_EHIC                                ;



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

NSIMS = 1;

options_.noprint = 1;
calibrated_values = M_.params;
Sigma_e = M_.Sigma_e;
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
    simdata = sur(simdata{idxs});
    BETA(i, :) = M_.params';
    oo_ = rmfield(oo_, 'sur');
end

if NSIMS > 1
    if max(abs(mean(BETA)' - calibrated_values)) > 1e-2
        error(['sum(abs(mean(BETA)'' - calibrated_values)) ' num2str(sum(abs(mean(BETA)' - calibrated_values)))]);
    end
else
    if isoctave
        good = [-8.352553221721005e-01
                -3.197321434625048e-01
                1.833371158471533e-02
                5.535617069191032e-02
                1.561983042003726e-02
                -8.993909626264719e-02
                1.711785670581707e-02
                4.757720040796121e-01
                -8.571091170339767e-02
                -2.950787496937523e-02
                -1.364615168219547e-01
                3.188305976033711e-03
                2.978809614745770e-01
                1.000152497234268e+00
                1.000581562867493e+00
                -8.285078166013811e-01
                -3.282738244438672e-01
                1.707491564143498e-02
                5.924757545427200e-02
                4.811130814635710e-02
                -1.066481566003144e-01
                4.462314484614635e-03
                4.862610106651586e-01
                -7.133910570135164e-02
                1.000000000000000e+00];
    else
        good = [-0.826686196809409
                -0.346753563700393
                0.063013991583949
                0.074802596658698
                -0.017440119721953
                -0.127090614348862
                0.025293280404460
                0.524290302468866
                -0.117611206771440
                -0.027776224547132
                -0.156590828735908
                0.054039707976331
                0.276257666502046
                1.000417289621684
                0.999336865129450
                -0.803258152338916
                -0.309594948488168
                0.051602756230521
                0.039275481081030
                0.024897596371662
                -0.096310133845385
                -0.022630284059365
                0.461683465196454
                -0.110278113383114
                1.000000000000000];
    end
    if max(abs(BETA' - good)) > 1e-14
        error(['sum of BETA'' - good was: ' num2str(sum(abs(BETA' - good)))]);
    end
    return
end

for i=1:nparampool
    figure
    hold on
    title(strrep(M_.param_names(i,:), '_', '\_'));
    histogram(BETA(:,i),50);
    line([calibrated_values(i) calibrated_values(i)], [0 NSIMS/10], 'LineWidth', 2, 'Color', 'r');
    hold off
end
