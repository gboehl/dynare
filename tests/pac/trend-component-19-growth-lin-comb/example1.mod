// --+ options: json=compute, stochastic +--

var x1 x2 x1bar x2bar z y gg;

varexo ex1 ex2 ex1bar ex2bar ez ey g;

parameters
       rho_1 rho_2
       a_x1_0 a_x1_1 a_x1_2 a_x1_x2_1 a_x1_x2_2
	   a_x2_0 a_x2_1 a_x2_2 a_x2_x1_1 a_x2_x1_2
	   e_c_m c_z_1 c_z_2 beta ;

rho_1 = .9;
rho_2 = -.2;

a_x1_0 =  -.9;
a_x1_1 =  .4;
a_x1_2 =  .3;
a_x1_x2_1 = .1;
a_x1_x2_2 = .2;

a_x2_0 =  -.9;
a_x2_1 =   .2;
a_x2_2 =  -.1;
a_x2_x1_1 = -.1;
a_x2_x1_2 = .2;

beta  =  .2;
e_c_m =  .5;
c_z_1 =  .2;
c_z_2 = -.1;

trend_component_model(model_name=toto, eqtags=['eq:x1', 'eq:x2', 'eq:x1bar', 'eq:x2bar'], targets=['eq:x1bar', 'eq:x2bar']);

pac_model(auxiliary_model_name=toto, discount=beta, growth=0.5*gg(-1)+beta+ex1, model_name=pacman);

model;

[name='eq:gg']
gg = g;

[name='eq:y']
y = rho_1*y(-1) + rho_2*y(-2) + ey;

[name='eq:x1']
diff(x1) = a_x1_0*(x1(-1)-x1bar(-1)) + a_x1_1*diff(x1(-1)) + a_x1_2*diff(x1(-2)) + a_x1_x2_1*diff(x2(-1)) + a_x1_x2_2*diff(x2(-2)) + ex1;     

[name='eq:x2']
diff(x2) = a_x2_0*(x2(-1)-x2bar(-1)) + a_x2_1*diff(x1(-1)) + a_x2_2*diff(x1(-2)) + a_x2_x1_1*diff(x2(-1)) + a_x2_x1_2*diff(x2(-2)) + ex2;     

[name='eq:x1bar']
x1bar = x1bar(-1) + ex1bar;

[name='eq:x2bar']
x2bar = x2bar(-1) + ex2bar;

[name='zpac']
diff(z) = e_c_m*(x1(-1)-z(-1)) + c_z_1*diff(z(-1))  + c_z_2*diff(z(-2)) + pac_expectation(pacman) + ez;

end;

// Initialize the PAC model (build the Companion VAR representation for the auxiliary model).
pac.initialize('pacman');

// Update the parameters of the PAC expectation model (h0 and h1 vectors).
pac.update.expectation('pacman');

// Set initial conditions to zero. Please use more sensible values if any...
initialconditions = dseries(zeros(10, M_.endo_nbr+M_.exo_nbr), 2000Q1, vertcat(M_.endo_names,M_.exo_names));

// Simulate the model for 5000 periods
TrueData = simul_backward_model(initialconditions, 5000);

// NLS estimation
// Define a structure describing the parameters to be estimated (with initial conditions). 
clear eparams
eparams.e_c_m = .5;
eparams.c_z_1 = .2;
eparams.c_z_2 =-.1;

edata = TrueData;
edata.ez = dseries(NaN(TrueData.nobs, 1), 2000Q1, 'ez');
pac.estimate.nls('zpac', eparams, edata, 2005Q1:2000Q1+200, 'lsqnonlin');

e_c_m_nls = M_.params(strmatch('e_c_m', M_.param_names, 'exact'));
c_z_1_nls = M_.params(strmatch('c_z_1', M_.param_names, 'exact'));
c_z_2_nls = M_.params(strmatch('c_z_2', M_.param_names, 'exact'));
resid_nls = oo_.pac.pacman.equations.eq0.residual;

fprintf('Estimate of e_c_m: %f \n', e_c_m_nls)
fprintf('Estimate of c_z_1: %f \n', c_z_1_nls)
fprintf('Estimate of c_z_2: %f \n', c_z_2_nls)

skipline(2)

// Iterative OLS estimation using estimates from NLS
// Define a structure describing the parameters to be estimated (with initial conditions). 
clear eparams
eparams.e_c_m = e_c_m_nls;
eparams.c_z_1 = c_z_1_nls;
eparams.c_z_2 = c_z_2_nls;

// Define the dataset used for estimation
edata = TrueData;
edata.ez = dseries(NaN(TrueData.nobs, 1), 2000Q1, 'ez');
pac.estimate.iterative_ols('zpac', eparams, edata, 2005Q1:2000Q1+200);

// Test printing of PAC expectations
pac.print('pacman','zpac');

// Print equations where the variable appears in
fprintf('x1bar is in: \n')
print_equations('x1bar')
fprintf('\n')

fprintf('x2bar is in: \n')
print_equations('x2bar', true);
fprintf('\n')

e_c_m_iterative_ols = M_.params(strmatch('e_c_m', M_.param_names, 'exact'));
c_z_1_iterative_ols = M_.params(strmatch('c_z_1', M_.param_names, 'exact'));
c_z_2_iterative_ols = M_.params(strmatch('c_z_2', M_.param_names, 'exact'));
resid_iterative_ols = oo_.pac.pacman.equations.eq0.residual;

fprintf('Estimate of e_c_m: %f \n', e_c_m_iterative_ols)
fprintf('Estimate of c_z_1: %f \n', c_z_1_iterative_ols)
fprintf('Estimate of c_z_2: %f \n', c_z_2_iterative_ols)

if any(abs(resid_nls-resid_iterative_ols).data > 1e-4)
   error('Iterative OLS and NLS do not provide consistent results.')
end
