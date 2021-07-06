// --+ options: stochastic,json=compute +--

var foo z x y;
varexo e_x e_y e_z;
parameters a b c d e f beta ;

a =  .9;
b = -.2;
c =  .3;
f =  .8;
d =  .5;
e =  .4;

beta = 1/(1+.02);

// Define a VAR model from a subset of equations in the model block.
var_model(model_name = toto, eqtags = [ 'X' 'Y' 'Z' ]);

/* Define a VAR_EXPECTATION_MODEL
** ------------------------------
**
** model_name:       the name of the VAR_EXPECTATION_MODEL (mandatory).
** var_model_name:   the name of the VAR model used for the expectations (mandatory).
** variable:         the name of the variable to be forecasted (mandatory).
** horizon:          the horizon forecast (mandatory).
** discount:         the discount factor, which can be a value or a declared parameter (default is 1.0, no discounting).
** time_shift:       shifts the information set to the past, must be a non positive scalar. By default, expectations 
**                   about `variable` in period `t+horizon` are formed in period t (time_shift=0) 
**
**
** The `horizon` parameter can be an integer in which case the (discounted) `horizon` step ahead forecast
** is computed using the VAR model `var_model_name`. Alternatively, `horizon` can be a range. In this case
** VAR_EXPECTATION_MODEL returns a discounted sum of expected values. If `horizon` is set equal to the range
** 0:Inf, then VAR_EXPECTATION_MODEL computes:
**
**                                       ∑ βʰ Eₜ[yₜ₊ₕ]
**
** where the sum is over h=0,…,∞ and the conditional expectations are computed with VAR model `var_model_name`.
*/

var_expectation_model(model_name = varexp, expression = x, auxiliary_model_name = toto, horizon = 1, time_shift=1, discount = beta)  ;


model;
[ name = 'X' ]
x = a*x(-1) + b*x(-2) + c*z(-2) + e_x;
[ name = 'Z' ]
z = f*z(-1) + e_z;
[ name = 'Y' ]
y = d*y(-2) + e*z(-1) + e_y;

foo = .5*foo(-1) + var_expectation(varexp);
end;


// Initialize the VAR expectation model, will build the companion matrix of the VAR.
var_expectation.initialize('varexp')

// Update VAR_EXPECTATION reduced form parameters
try
    var_expectation.update('varexp');
    if isfield(M_.var_expectation.varexp, 'time_shift')
        error('Preprocessor should not allow positive values for time_shift option.')
    end
catch
end
