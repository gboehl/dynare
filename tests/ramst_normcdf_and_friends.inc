// Tests the normcdf(), normpdf(), erf(), erfc(), trigonometric and hyperbolic functions

var c k t u v w sin_x cos_x tan_x asin_y acos_y atan_y sinh_x cosh_x tanh_x asinh_zz acosh_zz atanh_z;
varexo x y z zz;

parameters alph gam delt bet aa;
alph=0.5;
gam=0.5;
delt=0.02;
bet=0.05;
aa=0.5;

@#if model_options == ""
model;
@#else
model(@{model_options});
@#endif
c + k - aa*x*k(-1)^alph - (1-delt)*k(-1);
c^(-gam) - (1+bet)^(-1)*(aa*alph*x(+1)*k^(alph-1) + 1 - delt)*c(+1)^(-gam);
t = normcdf(x, 2, 3);
u = normpdf(x, 1, 0.5);
v = erf(x);
w = erfc(x);
sin_x = sin(x);
cos_x = cos(x);
tan_x = tan(x);
asin_y = asin(y);
acos_y = acos(y);
atan_y = atan(y);
sinh_x = sinh(x);
cosh_x = cosh(x);
tanh_x = tanh(x);
asinh_zz = asinh(zz);
acosh_zz = acosh(zz);
atanh_z = atanh(z);
end;

initval;
x = 1;
y = 0;
z = 0.1;
zz = 1.1;
k = ((delt+bet)/(1.0*aa*alph))^(1/(alph-1));
c = aa*k^alph-delt*k;
t = 0;
u = 0;
v = 0;
w = 0;
end;

steady;

check;

shocks;
var x;
periods 1;
values 1.2;
var y;
periods 1,2;
values 0.5, 1;
var z;
periods 1;
values 0.5;
var zz;
periods 1;
values 1.2;
end;

simul(periods=20);

if(abs(oo_.steady_state(5) - erf(1)) > 1e-10)
   error('Test failed in static @{mode_name} for erf')
end

if (abs(oo_.endo_simul(5, 2) - erf(1.2)) > 1e-10)
   error('Test failed in dynamic @{mode_name} for erf')
end

if(abs(oo_.steady_state(6) - erfc(1)) > 1e-10)
   error('Test failed in static @{mode_name} for erfc')
end

if (abs(oo_.endo_simul(6, 2) - erfc(1.2)) > 1e-10)
   error('Test failed in dynamic @{mode_name} for erfc')
end

if(abs(oo_.steady_state(4) - normpdf(1, 1, 0.5)) > 1e-10)
   error('Test failed in static @{mode_name} for normpdf')
end

if (abs(oo_.endo_simul(4, 2) - normpdf(1.2, 1, 0.5)) > 1e-10)
   error('Test failed in dynamic @{mode_name} for normpdf')
end

if (abs(oo_.steady_state(3) - normcdf(1, 2, 3)) > 1e-10)
   error('Test failed in static @{mode_name} for normcdf')
end

if (abs(oo_.endo_simul(3, 2) - normcdf(1.2, 2, 3)) > 1e-10)
   error('Test failed in dynamic @{mode_name} for normcdf')
end


sin_x_pos=strmatch('sin_x',M_.endo_names,'exact');
if(abs(oo_.steady_state(sin_x_pos) - sin(1)) > 1e-10)
   error('Test failed in static @{mode_name} for sin')
end

if (abs(oo_.endo_simul(sin_x_pos, 2) - sin(1.2)) > 1e-10)
   error('Test failed in dynamic @{mode_name} for sin')
end

cos_x_pos=strmatch('cos_x',M_.endo_names,'exact');
if(abs(oo_.steady_state(cos_x_pos) - cos(1)) > 1e-10)
   error('Test failed in static @{mode_name} for cos')
end

if (abs(oo_.endo_simul(cos_x_pos, 2) - cos(1.2)) > 1e-10)
   error('Test failed in dynamic @{mode_name} for cos')
end

tan_x_pos=strmatch('tan_x',M_.endo_names,'exact');
if(abs(oo_.steady_state(tan_x_pos) - tan(1)) > 1e-10)
   error('Test failed in static @{mode_name} for tan')
end

if (abs(oo_.endo_simul(tan_x_pos, 2) - tan(1.2)) > 1e-10)
   error('Test failed in dynamic @{mode_name} for tan')
end

asin_y_pos=strmatch('asin_y',M_.endo_names,'exact');
if(abs(oo_.steady_state(asin_y_pos) - asin(0)) > 1e-10)
   error('Test failed in static @{mode_name} for asin')
end

if (abs(oo_.endo_simul(asin_y_pos,1) - asin(0)) > 1e-10) || ...
        (abs(oo_.endo_simul(asin_y_pos, 2) - asin(0.5)) > 1e-10) || ...
        (abs(oo_.endo_simul(asin_y_pos, 3) - asin(1)) > 1e-10)
   error('Test failed in dynamic @{mode_name} for asin')
end

acos_y_pos=strmatch('acos_y',M_.endo_names,'exact');
if(abs(oo_.steady_state(acos_y_pos) - acos(0)) > 1e-10)
   error('Test failed in static @{mode_name} for acos')
end

if (abs(oo_.endo_simul(acos_y_pos,1) - acos(0)) > 1e-10) || ...
        (abs(oo_.endo_simul(acos_y_pos, 2) - acos(0.5)) > 1e-10) || ...
        (abs(oo_.endo_simul(acos_y_pos, 3) - acos(1)) > 1e-10)
   error('Test failed in dynamic @{mode_name} for acos')
end

atan_y_pos=strmatch('atan_y',M_.endo_names,'exact');
if(abs(oo_.steady_state(atan_y_pos) - atan(0)) > 1e-10)
   error('Test failed in static @{mode_name} for atan')
end

if (abs(oo_.endo_simul(atan_y_pos,1) - atan(0)) > 1e-10) || ...
        (abs(oo_.endo_simul(atan_y_pos, 2) - atan(0.5)) > 1e-10) || ...
        (abs(oo_.endo_simul(atan_y_pos, 3) - atan(1)) > 1e-10)
   error('Test failed in dynamic @{mode_name} for atan')
end

// Hyperbolic

sinh_x_pos=strmatch('sinh_x',M_.endo_names,'exact');
if(abs(oo_.steady_state(sinh_x_pos) - sinh(1)) > 1e-10)
   error('Test failed in static @{mode_name} for sinh')
end

if (abs(oo_.endo_simul(sinh_x_pos, 2) - sinh(1.2)) > 1e-10)
   error('Test failed in dynamic @{mode_name} for sinh')
end

cosh_x_pos=strmatch('cosh_x',M_.endo_names,'exact');
if(abs(oo_.steady_state(cosh_x_pos) - cosh(1)) > 1e-10)
   error('Test failed in static @{mode_name} for cosh')
end

if (abs(oo_.endo_simul(cosh_x_pos, 2) - cosh(1.2)) > 1e-10)
   error('Test failed in dynamic @{mode_name} for cosh')
end

tanh_x_pos=strmatch('tanh_x',M_.endo_names,'exact');
if(abs(oo_.steady_state(tanh_x_pos) - tanh(1)) > 1e-10)
   error('Test failed in static @{mode_name} for tanh')
end

if (abs(oo_.endo_simul(tanh_x_pos, 2) - tanh(1.2)) > 1e-10)
   error('Test failed in dynamic @{mode_name} for tanh')
end

asinh_zz_pos=strmatch('asinh_zz',M_.endo_names,'exact');
if(abs(oo_.steady_state(asinh_zz_pos) - asinh(1.1)) > 1e-10)
   error('Test failed in static @{mode_name} for asinh')
end

if (abs(oo_.endo_simul(asinh_zz_pos,1) - asinh(1.1)) > 1e-10) || ...
        (abs(oo_.endo_simul(asinh_zz_pos, 2) - asinh(1.2)) > 1e-10)
   error('Test failed in dynamic @{mode_name} for asinh')
end

acosh_zz_pos=strmatch('acosh_zz',M_.endo_names,'exact');
if(abs(oo_.steady_state(acosh_zz_pos) - acosh(1.1)) > 1e-10)
   error('Test failed in static @{mode_name} for acosh')
end

if (abs(oo_.endo_simul(acosh_zz_pos,1) - acosh(1.1)) > 1e-10) || ...
        (abs(oo_.endo_simul(acosh_zz_pos, 2) - acosh(1.2)) > 1e-10)
   error('Test failed in dynamic @{mode_name} for acosh')
end

atanh_z_pos=strmatch('atanh_z',M_.endo_names,'exact');
if(abs(oo_.steady_state(atanh_z_pos) - atanh(0.1)) > 1e-10)
   error('Test failed in static @{mode_name} for atanh')
end

if (abs(oo_.endo_simul(atanh_z_pos,1) - atanh(0.1)) > 1e-10) || ...
        (abs(oo_.endo_simul(atanh_z_pos, 2) - atanh(0.5)) > 1e-10)
   error('Test failed in dynamic @{mode_name} for atanh')
end

/*
 Local Variables:
   mode: dynare
 End:
 */
