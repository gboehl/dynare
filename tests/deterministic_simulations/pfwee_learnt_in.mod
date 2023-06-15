/* Tests perfect_foresight_with_expectation_errors_{setup,solver}
   using the shocks(learnt_in=…) and endval(learnt_in=…) syntax */

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

shocks(learnt_in = 1);
  var x;
  periods 1:2;
  values 1.2;
end;

shocks(learnt_in = 2);
  var x;
  periods 2;
  add 0.1;
end;

endval(learnt_in = 2);
  x = 1.1;
end;

shocks(learnt_in = 3);
  var x;
  periods 3 7;
  values 1.4 1.5;
end;

endval(learnt_in = 3);
  x += 0.1;
end;

// Dummy block, that will be overwritten by the next one
shocks(learnt_in = 6);
  var x;
  periods 6:8;
  values 10;
end;

shocks(learnt_in = 6, overwrite);
  var x;
  periods 6:7;
  multiply 0.8;
end;

endval(learnt_in = 6);
  x *= 0.75;
end;

// Save initial steady state (it will be modified by pfwee)
orig_steady_state = oo_.steady_state;
orig_exo_steady_state = oo_.exo_steady_state;

perfect_foresight_with_expectation_errors_setup(periods = 7);

perfect_foresight_with_expectation_errors_solver;
pfwee_simul = oo_.endo_simul;

// Now compute the solution by hand to verify the results
oo_.steady_state = orig_steady_state;
oo_.exo_steady_state = orig_exo_steady_state;

perfect_foresight_setup;

// Information arriving in period 1 (temp shock now and tomorrow)
oo_.exo_simul(2:3,1) = 1.2;
perfect_foresight_solver;

// Information arriving in period 2 (temp shock now + permanent shock in future)
oo_.exo_simul(3,1) = 1.3;
oo_.exo_steady_state = 1.1;
oo_.exo_simul(end, 1) = oo_.exo_steady_state;
oo_.steady_state = evaluate_steady_state(oo_.steady_state, oo_.exo_steady_state, M_, options_, true);
oo_.endo_simul(:, end) = oo_.steady_state;
options_.periods = 6;
saved_endo = oo_.endo_simul(:, 1);
saved_exo = oo_.exo_simul(1, :);
oo_.endo_simul = oo_.endo_simul(:, 2:end);
oo_.exo_simul = oo_.exo_simul(2:end, :);
perfect_foresight_solver;
oo_.endo_simul = [ saved_endo oo_.endo_simul ];
oo_.exo_simul = [ saved_exo; oo_.exo_simul ];

// Information arriving in period 3 (temp shocks + permanent shock in future)
oo_.exo_simul(4,1) = 1.4;
oo_.exo_simul(8,1) = 1.5;
oo_.exo_steady_state = 1.1+0.1;
oo_.exo_simul(end, 1) = oo_.exo_steady_state;
oo_.steady_state = evaluate_steady_state(oo_.steady_state, oo_.exo_steady_state, M_, options_, true);
oo_.endo_simul(:, end) = oo_.steady_state;
options_.periods = 5;
saved_endo = oo_.endo_simul(:, 1:2);
saved_exo = oo_.exo_simul(1:2, :);
oo_.endo_simul = oo_.endo_simul(:, 3:end);
oo_.exo_simul = oo_.exo_simul(3:end, :);
perfect_foresight_solver;
oo_.endo_simul = [ saved_endo oo_.endo_simul ];
oo_.exo_simul = [ saved_exo; oo_.exo_simul ];

// Information arriving in period 6 (temp shocks + permanent shock)
oo_.exo_simul(7,1) = 1*0.8;
oo_.exo_simul(8,1) = 1.5*0.8;
oo_.exo_steady_state = (1.1+0.1)*0.75;
oo_.exo_simul(end, 1) = oo_.exo_steady_state;
oo_.steady_state = evaluate_steady_state(oo_.steady_state, oo_.exo_steady_state, M_, options_, true);
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
if any(any(pfwee_simul-oo_.endo_simul ~= 0))
    error('Error in perfect_foresight_with_expectation_errors')
end
