// --+ options: stochastic,transform_unary_ops,json=compute +--

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

@#include "example/model/var-expectations/varexp-parameters.inc"

beta = 1/(1+.02);

model;

[name='eq:x1', data_type='nonstationary']
diff(x1) = a_x1_0*(x1(-1)-x1bar(-1))+a_x1_0_*(x2(-1)-x2bar(-1)) + a_x1_1*diff(x1(-1)) + a_x1_2*diff(x1(-2)) + + a_x1_x2_1*diff(x2(-1)) + a_x1_x2_2*diff(x2(-2)) + ex1;

[name='eq:x2', data_type='nonstationary']
diff(x2) = a_x2_0*(x2(-1)-x2bar(-1)) + a_x2_1*diff(x1(-1)) + a_x2_2*diff(x1(-2)) + + a_x2_x1_1*diff(x2(-1)) + a_x2_x1_2*diff(x2(-2)) + ex2;

[name='eq:x1bar', data_type='nonstationary']
x1bar = x1bar(-1) + ex1bar;

[name='eq:x2bar', data_type='nonstationary']
x2bar = x2bar(-1) + ex2bar;

foo = .5*foo(-1) +
@#include "example/model/var-expectations/varexp-expression.inc"
;
end;

shocks;
  var ex1 = .01;
  var ex2 = .01;
  var ex1bar = .02;
  var ex2bar = .02;
end;


verbatim;
  initialconditions =zeros(3,5);
  initialconditions(3,1) = .10; % foo(-1)
  initialconditions(3,2) = .20; % x1(-1)
  initialconditions(2,2) = .22; % x1(-2)
  initialconditions(1,2) = .24; % x1(-3)
  initialconditions(3,3) = .30; % x2(-1)
  initialconditions(2,3) = .32; % x2(-2)
  initialconditions(1,3) = .34; % x2(-3)
  initialconditions(3,4) = .25; % x1bar(-1)
  initialconditions(3,5) = .25; % x2bar(-1)
  initialconditions = ...
  dseries(initialconditions, dates('2000Q1'), {'foo', 'x1','x2', 'x1bar', 'x2bar'});
    set_dynare_seed('default');
  ts = simul_backward_model(initialconditions, 100);
  ex = load('example.mat');
  delete('example.mat')
  if max(abs(ex.foo-ts.foo.data))>1e-12
     error('Simulations do not match!')
  end
end;