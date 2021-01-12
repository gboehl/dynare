var y z;

varexo e x u;

parameters rho;

rho = .9;


model(linear);
y = z +  e ;
z = rho*z(-1) + u;
end;

shocks;
var e; stderr .01;
var u; stderr .01;
end;

varobs y;

varexobs x;

assert(length(options_.varobs)==1);
assert(length(options_.varexobs)==1);
assert(length(options_.varobs_id)==1);
assert(length(options_.varexobs_id)==1);
assert(options_.varobs_id==1);
assert(options_.varexobs_id==2);
assert(options_.varobs{1}=='y');
assert(options_.varexobs{1}=='x');

