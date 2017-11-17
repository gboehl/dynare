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

diff(U2_Q_YED) =   u2_q_yed_ecm_u2_q_yed_L1 * (U2_Q_YED(-1) - U2_EHIC(-1))
                 + u2_q_yed_ecm_u2_stn_L1   * (U2_STN(-1)   - U2_ESTN(-1))
                 + u2_q_yed_u2_g_yer_L1     * diff(U2_G_YER(-1))
                 + u2_q_yed_u2_stn_L1       * diff(U2_STN(-1))
                 + res_U2_Q_YED                                           ;

diff(U2_G_YER) =   u2_g_yer_ecm_u2_q_yed_L1 * (U2_Q_YED(-1) - U2_EHIC(-1))
                 + u2_g_yer_ecm_u2_stn_L1   * (U2_STN(-1)   - U2_ESTN(-1))
                 + u2_g_yer_u2_q_yed_L1     * diff(U2_Q_YED(-1))
                 + u2_g_yer_u2_g_yer_L1     * diff(U2_G_YER(-1))
                 + u2_g_yer_u2_stn_L1       * diff(U2_STN(-1))
                 + res_U2_G_YER                                           ;

diff(U2_STN)   =   u2_stn_ecm_u2_q_yed_L1   * (U2_Q_YED(-1) - U2_EHIC(-1))
                 + u2_stn_ecm_u2_stn_L1     * (U2_STN(-1)   - U2_ESTN(-1))
                 + u2_stn_u2_q_yed_L1       * diff(U2_Q_YED(-1))
                 + u2_stn_u2_g_yer_L1       * diff(U2_G_YER(-1))
                 + res_U2_STN                                             ;

U2_ESTN        =   u2_estn_u2_estn_L1       * U2_ESTN(-1)
                 + res_U2_ESTN                                            ;

U2_EHIC        =   u2_ehic_u2_ehic_L1       * U2_EHIC(-1)
                 + res_U2_EHIC                                            ;

diff(DE_Q_YED) =   de_q_yed_ecm_de_q_yed_L1 * (DE_Q_YED(-1) - DE_EHIC(-1))
                 + de_q_yed_ecm_u2_stn_L1   * (U2_STN(-1)   - U2_ESTN(-1))
                 + de_q_yed_de_g_yer_L1     * diff(DE_G_YER(-1))
                 + de_q_yed_u2_stn_L1       * diff(U2_STN(-1))
                 + res_DE_Q_YED                                           ;

diff(DE_G_YER) =   de_g_yer_ecm_de_q_yed_L1 * (DE_Q_YED(-1) - DE_EHIC(-1))
                 + de_g_yer_ecm_u2_stn_L1   * (U2_STN(-1)   - U2_ESTN(-1))
                 + de_g_yer_de_q_yed_L1     * diff(DE_Q_YED(-1))
                 + de_g_yer_de_g_yer_L1     * diff(DE_G_YER(-1))
                 + de_g_yer_u2_stn_L1       * diff(U2_STN(-1))
                 + res_DE_G_YER                                           ;

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

