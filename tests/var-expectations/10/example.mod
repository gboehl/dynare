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
var_expectation_model(model_name = varexp, expression = diff(log(x)), auxiliary_model_name = toto, horizon = 1, discount = beta)  ;

model;
[ name = 'X' ]
diff(log(x)) = a*diff(log(x(-1))) + b*diff(log(x(-2))) + c*diff(z(-2)) + e_x;
[ name = 'Z' ]
diff(z) = f*diff(z(-1)) + e_z;
[ name = 'Y' ]
log(y) = d*log(y(-2)) + e*diff(z(-1)) + e_y;

foo = var_expectation(varexp);
end;

// Initialize the VAR expectation model, will build the companion matrix of the VAR.
var_expectation.initialize('varexp')

// Update VAR_EXPECTATION reduced form parameters
var_expectation.update('varexp');

// Print expanded VAR_EXPECTATION expression in a file (to be included in substitution.mod).
var_expectation.print('varexp');

shocks;
  var e_x = .01;
  var e_y = .01;
  var e_z = .01;
end;

verbatim;
  initialconditions =zeros(3,4);
  initialconditions(3,1) = .1; % foo(-1)
  initialconditions(:,2) = .2; % y(-1)
  initialconditions(3,3) = .3; % z(-1)
  initialconditions(2,3) = .4; % z(-2)
  initialconditions(3,4) = .5; % x(-1)
  initialconditions(2,4) = .6; % x(-2)
  initialconditions(1,4) = .7; % x(-3)
  initialconditions = ...
  dseries(initialconditions, dates('2000Q1'), {'foo', 'y','z', 'x'});
    set_dynare_seed('default');
  ts = simul_backward_model(initialconditions, 15);
  foo = ts.foo.data;
  % Evaluate the (VAR) expectation term
  ts = example.var_expectations.evaluate_varexp(ts);
  % Check tthat the evaluation is correct.
  range = dates('2000Q4'):dates('2004Q2');
  if max(abs(ts(range).foo.data-ts(range).varexp.data))>1e-5
     error('Expectation term evaluations do not match!')
  end
end;