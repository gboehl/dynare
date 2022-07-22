% Test solve_algo=0 (fsolve), and its companion “fsolve_options” option

var c k;
varexo x;

parameters alph gam delt bet aa;
alph=0.5;
gam=0.5;
delt=0.02;
bet=0.05;
aa=0.5;


model;
c + k - aa*x*k(-1)^alph - (1-delt)*k(-1);
c^(-gam) - (1+bet)^(-1)*(aa*alph*x(+1)*k^(alph-1) + 1 - delt)*c(+1)^(-gam);
end;

initval;
x = 1;
k = 1;
c = 0.1;
end;

if isoctave
steady(solve_algo=0, fsolve_options = ('AutoScaling', 'on'));
elseif user_has_matlab_license('optimization_toolbox')
%% NB: The Display option is accepted but not honoured under Octave (as of version 7)
steady(solve_algo=0, fsolve_options = ('Display', 'iter'));
end
