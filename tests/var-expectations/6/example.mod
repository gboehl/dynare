// --+ options: stochastic,json=compute +--

var foo x1 x2 x1bar x2bar;

varexo ex1 ex2 ex1bar ex2bar;

parameters a_x1_0 a_x1_0_ a_x1_1 a_x1_2 a_x1_x2_1 a_x1_x2_2
	   a_x2_0 a_x2_1 a_x2_2 a_x2_x1_1 a_x2_x1_2
       beta ;

a_x1_0 =  -.9;
a_x1_0_ =  -.8;
a_x1_1 =  .4;
a_x1_2 =  .3;
a_x1_x2_1 = .1;
a_x1_x2_2 = .2;


a_x2_0 =  -.9;
a_x2_1 =   .2;
a_x2_2 =  -.1;
a_x2_x1_1 = -.1;
a_x2_x1_2 = .2;

beta = 1/(1+.02);

// Define a TREND_COMPONENT model from a subset of equations in the model block.
trend_component_model(model_name=toto, eqtags=['eq:x1', 'eq:x2', 'eq:x1bar', 'eq:x2bar'], targets=['eq:x1bar', 'eq:x2bar']);

/* Define a VAR_EXPECTATION_MODEL
** ------------------------------
**
** model_name:             the name of the VAR_EXPECTATION_MODEL (mandatory).
** auxiliary_model_name:   the name of the VAR model used for the expectations (mandatory).
** variable:               the name of the variable to be forecasted (mandatory).
** horizon:                the horizon forecast (mandatory).
** discount:               the discount factor, which can be a value or a declared parameter (default is 1.0, no discounting).
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

var_expectation_model(model_name = varexp, variable = x1, auxiliary_model_name = toto, horizon = 15:50, discount = beta)  ;


model;

[name='eq:x1', data_type='nonstationary']
diff(x1) = a_x1_0*(x1(-1)-x1bar(-1))+a_x1_0_*(x2(-1)-x2bar(-1)) + a_x1_1*diff(x1(-1)) + a_x1_2*diff(x1(-2)) + + a_x1_x2_1*diff(x2(-1)) + a_x1_x2_2*diff(x2(-2)) + ex1;     

[name='eq:x2', data_type='nonstationary']
diff(x2) = a_x2_0*(x2(-1)-x2bar(-1)) + a_x2_1*diff(x1(-1)) + a_x2_2*diff(x1(-2)) + + a_x2_x1_1*diff(x2(-1)) + a_x2_x1_2*diff(x2(-2)) + ex2;     

[name='eq:x1bar', data_type='nonstationary']
x1bar = x1bar(-1) + ex1bar;

[name='eq:x2bar', data_type='nonstationary']
x2bar = x2bar(-1) + ex2bar;

foo = .5*foo(-1) + var_expectation(varexp);
end;

// Initialize the VAR expectation model, will build the companion matrix of the VAR.
var_expectation.initialize('varexp')

// Update VAR_EXPECTATION reduced form parameters
var_expectation.update('varexp');

weights = M_.params(M_.var_expectation.varexp.param_indices);

if ~all(weights(1:6)) || ~all(weights(9:10)) || weights(7) || weights(8) || weights(11) || weights(12)
   error('Wrong reduced form parameter for VAR_EXPECTATION_MODEL')
end