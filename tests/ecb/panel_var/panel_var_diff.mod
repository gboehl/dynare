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
FR_Q_YED
FR_G_YER
FR_EHIC
IT_Q_YED
IT_G_YER
IT_EHIC
ES_Q_YED
ES_G_YER
ES_EHIC
NL_Q_YED
NL_G_YER
NL_EHIC
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
res_FR_Q_YED
res_FR_G_YER
res_FR_EHIC
res_IT_Q_YED
res_IT_G_YER
res_IT_EHIC
res_ES_Q_YED
res_ES_G_YER
res_ES_EHIC
res_NL_Q_YED
res_NL_G_YER
res_NL_EHIC
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

fr_q_yed_ecm_fr_q_yed_L1
fr_q_yed_ecm_u2_stn_L1
fr_q_yed_fr_g_yer_L1
fr_q_yed_u2_stn_L1
fr_g_yer_ecm_fr_q_yed_L1
fr_g_yer_ecm_u2_stn_L1
fr_g_yer_fr_q_yed_L1
fr_g_yer_fr_g_yer_L1
fr_g_yer_u2_stn_L1
fr_ehic_fr_ehic_L1

it_q_yed_ecm_it_q_yed_L1
it_q_yed_ecm_u2_stn_L1
it_q_yed_it_g_yer_L1
it_q_yed_u2_stn_L1
it_g_yer_ecm_it_q_yed_L1
it_g_yer_ecm_u2_stn_L1
it_g_yer_it_q_yed_L1
it_g_yer_it_g_yer_L1
it_g_yer_u2_stn_L1
it_ehic_it_ehic_L1

es_q_yed_ecm_es_q_yed_L1
es_q_yed_ecm_u2_stn_L1
es_q_yed_es_g_yer_L1
es_q_yed_u2_stn_L1
es_g_yer_ecm_es_q_yed_L1
es_g_yer_ecm_u2_stn_L1
es_g_yer_es_q_yed_L1
es_g_yer_es_g_yer_L1
es_g_yer_u2_stn_L1
es_ehic_es_ehic_L1

nl_q_yed_ecm_nl_q_yed_L1
nl_q_yed_ecm_u2_stn_L1
nl_q_yed_nl_g_yer_L1
nl_q_yed_u2_stn_L1
nl_g_yer_ecm_nl_q_yed_L1
nl_g_yer_ecm_u2_stn_L1
nl_g_yer_nl_q_yed_L1
nl_g_yer_nl_g_yer_L1
nl_g_yer_u2_stn_L1
nl_ehic_nl_ehic_L1

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

fr_q_yed_ecm_fr_q_yed_L1  = -0.822375165893149  ;
fr_q_yed_ecm_u2_stn_L1    = -0.323715338568977  ;
fr_q_yed_fr_g_yer_L1      =  0.0401361895021082 ;
fr_q_yed_u2_stn_L1        =  0.0583977039584461 ;
fr_g_yer_ecm_fr_q_yed_L1  =  0.0189896046977422 ;
fr_g_yer_ecm_u2_stn_L1    = -0.109597659887433  ;
fr_g_yer_fr_q_yed_L1      =  0.00376679676320256;
fr_g_yer_fr_g_yer_L1      =  0.480506381923643  ;
fr_g_yer_u2_stn_L1        = -0.0722359286123494 ;
fr_ehic_fr_ehic_L1        =  1                  ;

it_q_yed_ecm_it_q_yed_L1  = -0.822375165893149  ;
it_q_yed_ecm_u2_stn_L1    = -0.323715338568977  ;
it_q_yed_it_g_yer_L1      =  0.0401361895021082 ;
it_q_yed_u2_stn_L1        =  0.0583977039584461 ;
it_g_yer_ecm_it_q_yed_L1  =  0.0189896046977422 ;
it_g_yer_ecm_u2_stn_L1    = -0.109597659887433  ;
it_g_yer_it_q_yed_L1      =  0.00376679676320256;
it_g_yer_it_g_yer_L1      =  0.480506381923643  ;
it_g_yer_u2_stn_L1        = -0.0722359286123494 ;
it_ehic_it_ehic_L1        =  1                  ;

es_q_yed_ecm_es_q_yed_L1  = -0.822375165893149  ;
es_q_yed_ecm_u2_stn_L1    = -0.323715338568977  ;
es_q_yed_es_g_yer_L1      =  0.0401361895021082 ;
es_q_yed_u2_stn_L1        =  0.0583977039584461 ;
es_g_yer_ecm_es_q_yed_L1  =  0.0189896046977422 ;
es_g_yer_ecm_u2_stn_L1    = -0.109597659887433  ;
es_g_yer_es_q_yed_L1      =  0.00376679676320256;
es_g_yer_es_g_yer_L1      =  0.480506381923643  ;
es_g_yer_u2_stn_L1        = -0.0722359286123494 ;
es_ehic_es_ehic_L1        =  1                  ;

nl_q_yed_ecm_nl_q_yed_L1  = -0.822375165893149  ;
nl_q_yed_ecm_u2_stn_L1    = -0.323715338568977  ;
nl_q_yed_nl_g_yer_L1      =  0.0401361895021082 ;
nl_q_yed_u2_stn_L1        =  0.0583977039584461 ;
nl_g_yer_ecm_nl_q_yed_L1  =  0.0189896046977422 ;
nl_g_yer_ecm_u2_stn_L1    = -0.109597659887433  ;
nl_g_yer_nl_q_yed_L1      =  0.00376679676320256;
nl_g_yer_nl_g_yer_L1      =  0.480506381923643  ;
nl_g_yer_u2_stn_L1        = -0.0722359286123494 ;
nl_ehic_nl_ehic_L1        =  1                  ;


model;

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

U2_ESTN        =   u2_estn_u2_estn_L1       * U2_ESTN
                 + res_U2_ESTN                                            ;

U2_EHIC        =   u2_ehic_u2_ehic_L1       * U2_EHIC
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

DE_EHIC        =   de_ehic_de_ehic_L1       * DE_EHIC
                 + res_DE_EHIC                                            ;

diff(FR_Q_YED) =   fr_q_yed_ecm_fr_q_yed_L1 * (FR_Q_YED(-1) - FR_EHIC(-1))
                 + fr_q_yed_ecm_u2_stn_L1   * (U2_STN(-1)   - U2_ESTN(-1))
                 + fr_q_yed_fr_g_yer_L1     * diff(FR_G_YER(-1))
                 + fr_q_yed_u2_stn_L1       * diff(U2_STN(-1))
                 + res_FR_Q_YED                                           ;

diff(FR_G_YER) =   fr_g_yer_ecm_fr_q_yed_L1 * (FR_Q_YED(-1) - FR_EHIC(-1))
                 + fr_g_yer_ecm_u2_stn_L1   * (U2_STN(-1)   - U2_ESTN(-1))
                 + fr_g_yer_fr_q_yed_L1     * diff(FR_Q_YED(-1))
                 + fr_g_yer_fr_g_yer_L1     * diff(FR_G_YER(-1))
                 + fr_g_yer_u2_stn_L1       * diff(U2_STN(-1))
                 + res_FR_G_YER                                           ;

FR_EHIC        =   fr_ehic_fr_ehic_L1       * FR_EHIC
                 + res_FR_EHIC                                            ;

diff(IT_Q_YED) =   it_q_yed_ecm_it_q_yed_L1 * (IT_Q_YED(-1) - IT_EHIC(-1))
                 + it_q_yed_ecm_u2_stn_L1   * (U2_STN(-1)   - U2_ESTN(-1))
                 + it_q_yed_it_g_yer_L1     * diff(IT_G_YER(-1))
                 + it_q_yed_u2_stn_L1       * diff(U2_STN(-1))
                 + res_IT_Q_YED                                           ;

diff(IT_G_YER) =   it_g_yer_ecm_it_q_yed_L1 * (IT_Q_YED(-1) - IT_EHIC(-1))
                 + it_g_yer_ecm_u2_stn_L1   * (U2_STN(-1)   - U2_ESTN(-1))
                 + it_g_yer_it_q_yed_L1     * diff(IT_Q_YED(-1))
                 + it_g_yer_it_g_yer_L1     * diff(IT_G_YER(-1))
                 + it_g_yer_u2_stn_L1       * diff(U2_STN(-1))
                 + res_IT_G_YER                                           ;

IT_EHIC        =   it_ehic_it_ehic_L1       * IT_EHIC
                 + res_IT_EHIC                                            ;

diff(ES_Q_YED) =   es_q_yed_ecm_es_q_yed_L1 * (ES_Q_YED(-1) - ES_EHIC(-1))
                 + es_q_yed_ecm_u2_stn_L1   * (U2_STN(-1)   - U2_ESTN(-1))
                 + es_q_yed_es_g_yer_L1     * diff(ES_G_YER(-1))
                 + es_q_yed_u2_stn_L1       * diff(U2_STN(-1))
                 + res_ES_Q_YED                                           ;

diff(ES_G_YER) =   es_g_yer_ecm_es_q_yed_L1 * (ES_Q_YED(-1) - ES_EHIC(-1))
                 + es_g_yer_ecm_u2_stn_L1   * (U2_STN(-1)   - U2_ESTN(-1))
                 + es_g_yer_es_q_yed_L1     * diff(ES_Q_YED(-1))
                 + es_g_yer_es_g_yer_L1     * diff(ES_G_YER(-1))
                 + es_g_yer_u2_stn_L1       * diff(U2_STN(-1))
                 + res_ES_G_YER                                           ;

ES_EHIC        =   es_ehic_es_ehic_L1       * ES_EHIC
                 + res_ES_EHIC                                            ;

diff(NL_Q_YED) =   nl_q_yed_ecm_nl_q_yed_L1 * (NL_Q_YED(-1) - NL_EHIC(-1))
                 + nl_q_yed_ecm_u2_stn_L1   * (U2_STN(-1)   - U2_ESTN(-1))
                 + nl_q_yed_nl_g_yer_L1     * diff(NL_G_YER(-1))
                 + nl_q_yed_u2_stn_L1       * diff(U2_STN(-1))
                 + res_NL_Q_YED                                           ;

diff(NL_G_YER) =   nl_g_yer_ecm_nl_q_yed_L1 * (NL_Q_YED(-1) - NL_EHIC(-1))
                 + nl_g_yer_ecm_u2_stn_L1   * (U2_STN(-1)   - U2_ESTN(-1))
                 + nl_g_yer_nl_q_yed_L1     * diff(NL_Q_YED(-1))
                 + nl_g_yer_nl_g_yer_L1     * diff(NL_G_YER(-1))
                 + nl_g_yer_u2_stn_L1       * diff(U2_STN(-1))
                 + res_NL_G_YER                                           ;

NL_EHIC        =   nl_ehic_nl_ehic_L1       * NL_EHIC
                 + res_NL_EHIC                                            ;


end;

mydseries = dseries(randn(40,20), 1, {'U2_Q_YED', ...
                                      'U2_G_YER', ...
                                      'U2_STN', ...
                                      'U2_ESTN', ...
                                      'U2_EHIC', ...
                                      'DE_Q_YED', ...
                                      'DE_G_YER', ...
                                      'DE_EHIC', ...
                                      'FR_Q_YED', ...
                                      'FR_G_YER', ...
                                      'FR_EHIC', ...
                                      'IT_Q_YED', ...
                                      'IT_G_YER', ...
                                      'IT_EHIC', ...
                                      'ES_Q_YED', ...
                                      'ES_G_YER', ...
                                      'ES_EHIC', ...
                                      'NL_Q_YED', ...
                                      'NL_G_YER', ...
                                      'NL_EHIC'});

pooled_ols(mydseries, ...
           {'de','u2','fr', 'it', 'es', 'nl'}, ...
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
