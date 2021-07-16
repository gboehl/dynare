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
var_model(structural, model_name = toto, eqtags = [ 'X' 'Y' 'Z' ]);

// Define a VAR_EXPECTATION_MODEL
var_expectation_model(model_name = varexp, expression = diff(log(x)), auxiliary_model_name = toto, horizon = 1, discount = beta)  ;

model;
[ name = 'X' ]
diff(log(x)) = b*diff(z) + a*diff(log(x(-1))) + (1-a)*diff(log(x(-2))) + c*diff(z(-2)) + e_x;
[ name = 'Z' ]
diff(z) = f*(diff(z(-1))-diff(log(x)))+c*diff(z(-2)) + e_z;
[ name = 'Y' ]
log(y) = diff(log(x)) + d*log(y(-2)) + e*diff(z(-1)) + e_y;

foo = var_expectation(varexp);
end;

[ar, a0] = example1.varmatrices('toto', M_.params);

assert(isequal(diag(a0), ones(3,1)), 'Diagonal of a0 is wrong.')

assert(a0(1,3)==-b, 'Element (1,3) in A0 is wrong.')
assert(a0(1,2)==0, 'Element (1,2) in A0 is wrong.')
assert(a0(2,1)==-1, 'Element (2,1) in A0 is wrong.')
assert(a0(2,3)==0, 'Element (2,3) in A0 is wrong.')
assert(a0(3,1)==f, 'Element (3,1) in A0 is wrong.')
assert(a0(3,2)==0, 'Element (3,1) in A0 is wrong.')

assert(isequal(ar(:,:,1), [a 0 0; 0 0 e; 0 0 f]), 'First autoregressive matrix is wrong');
assert(isequal(ar(:,:,2), [1-a 0 c; 0 d 0; 0 0 c]), 'Second autoregressive matrix is wrong');

