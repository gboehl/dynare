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

// Define a VAR_EXPECTATION_MODEL
var_expectation_model(model_name = varexp, variable = x, auxiliary_model_name = toto, horizon = 1, discount = beta)  ;

model;
[ name = 'X' ]
diff(log(x)) = a*diff(log(x(-1))) + b*diff(log(x(-2))) + c*diff(z(-2)) + e_x;
[ name = 'Z' ]
diff(z) = f*diff(z(-1)) + e_z;
[ name = 'Y' ]
log(y) = d*log(y(-2)) + e*diff(z(-1)) + e_y;

foo = .5*foo(-1) + var_expectation(varexp);
end;


// Initialize the VAR expectation model, will build the companion matrix of the VAR.
var_expectation.initialize('varexp')

// Update VAR_EXPECTATION reduced form parameters
var_expectation.update('varexp');

var_expectation.print('varexp');

verbatim;
  set_dynare_seed('default');
  y = zeros(M_.endo_nbr,1);
  y(1:M_.orig_endo_nbr) = rand(M_.orig_endo_nbr, 1);
  x = randn(M_.exo_nbr,1);
  y = example.set_auxiliary_variables(y, x, M_.params);
  y = [y(find(M_.lead_lag_incidence(1,:))); y];
  [residual, g1] = example.dynamic(y, x', M_.params, oo_.steady_state, 1);
  save('example.mat', 'residual', 'g1');
end;