/* Tests perfect_foresight_with_expectation_errors_{setup,solver}
   using a CSV datafile */

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
k = ((delt+bet)/(1.0*aa*alph))^(1/(alph-1));
c = aa*k^alph-delt*k;
end;

steady;

check;

perfect_foresight_with_expectation_errors_setup(periods = 7, datafile = 'pfwee.csv');

// First simulation with default options
perfect_foresight_with_expectation_errors_solver;
pfwee1 = oo_.endo_simul;

// Second simulation with alternative guess values
perfect_foresight_with_expectation_errors_solver(terminal_steady_state_as_guess_value);
pfwee2 = oo_.endo_simul;

// Now compute the solution by hand to verify the results

perfect_foresight_setup;

// Information arriving in period 1 (temp shock now)
oo_.exo_simul(2,1) = 1.2;
perfect_foresight_solver;

// Information arriving in period 2 (temp shock now + permanent shock in future)
oo_.exo_simul(3,1) = 1.3;
oo_.exo_steady_state = 1.1;
oo_.exo_simul(end, 1) = oo_.exo_steady_state;
oo_.steady_state = evaluate_steady_state(oo_.steady_state, M_, options_, oo_, true);
oo_.endo_simul(:, end) = oo_.steady_state;
options_.periods = 6;
saved_endo = oo_.endo_simul(:, 1);
saved_exo = oo_.exo_simul(1, :);
oo_.endo_simul = oo_.endo_simul(:, 2:end);
oo_.exo_simul = oo_.exo_simul(2:end, :);
perfect_foresight_solver;
oo_.endo_simul = [ saved_endo oo_.endo_simul ];
oo_.exo_simul = [ saved_exo; oo_.exo_simul ];

// Information arriving in period 3 (temp shock now + permanent shock in future)
oo_.exo_simul(4,1) = 1.4;
oo_.exo_steady_state = 1.2;
oo_.exo_simul(end, 1) = oo_.exo_steady_state;
oo_.steady_state = evaluate_steady_state(oo_.steady_state, M_, options_, oo_, true);
oo_.endo_simul(:, end) = oo_.steady_state;
options_.periods = 5;
saved_endo = oo_.endo_simul(:, 1:2);
saved_exo = oo_.exo_simul(1:2, :);
oo_.endo_simul = oo_.endo_simul(:, 3:end);
oo_.exo_simul = oo_.exo_simul(3:end, :);
perfect_foresight_solver;
oo_.endo_simul = [ saved_endo oo_.endo_simul ];
oo_.exo_simul = [ saved_exo; oo_.exo_simul ];

// Information arriving in period 6 (permanent shock arriving now)
oo_.exo_simul(7,1) = 1.1;
oo_.exo_simul(8,1) = 1.1;
oo_.exo_steady_state = 1.1;
oo_.exo_simul(end, 1) = oo_.exo_steady_state;
oo_.steady_state = evaluate_steady_state(oo_.steady_state, M_, options_, oo_, true);
oo_.endo_simul(:, end) = oo_.steady_state;
options_.periods = 2;
saved_endo = oo_.endo_simul(:, 1:5);
saved_exo = oo_.exo_simul(1:5, :);
oo_.endo_simul = oo_.endo_simul(:, 6:end);
oo_.exo_simul = oo_.exo_simul(6:end, :);
perfect_foresight_solver;
oo_.endo_simul = [ saved_endo oo_.endo_simul ];
oo_.exo_simul = [ saved_exo; oo_.exo_simul ];

% We should have strict equality with first pfwee simulation, because algorithm
% and guess values are exactly the same.
if any(any(pfwee1-oo_.endo_simul ~= 0))
    error('Error in perfect_foresight_with_expectation_errors')
end

% For the 2nd simulation, since the guess values are different, there are some
% numerical differences
if max(max(abs(pfwee2-oo_.endo_simul))) > 10*options_.dynatol.f
    error('Error in perfect_foresight_with_expectation_errors + terminal_steady_state_as_guess_value')
end
