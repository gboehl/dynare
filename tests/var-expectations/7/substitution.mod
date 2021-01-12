// --+ options: stochastic,transform_unary_ops,json=compute +--

var foo z x y;
varexo e_x e_y e_z;
parameters a b c d e f beta ;

a =  .9;
b = -.2;
c =  .3;
f =  .8;
d =  .5;
e =  .4;

@#include "example/model/var-expectations/varexp-parameters.inc"

beta = 1/(1+.02);

// Define a VAR_EXPECTATION_MODEL
var_model(model_name = toto, eqtags = [ 'X' 'Y' 'Z' ]);

model;
[ name = 'X' ]
diff(x) = a*diff(x(-1)) + b*diff(x(-2)) + c*z(-2) + e_x;
[ name = 'Z' ]
z = f*z(-1) + e_z;
[ name = 'Y' ]
log(y) = d*log(y(-2)) + e*z(-1) + e_y;

foo = .5*foo(-1) +
@#include "example/model/var-expectations/varexp-expression.inc"
;
end;

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
  ts = simul_backward_model(initialconditions, 100);
  ex = load('example.mat');
end;

delete('example.mat')

if max(abs(ex.foo-ts.foo.data))>1e-12
   error('Simulations do not match!')
end

