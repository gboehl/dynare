function [endo_simul,endo_simul_no_constraint,status] = occbin(M,oo,options)
% function oo=occbin(M,oo,options) solves linear models with occasionally
% binding constraints using OCCBIN by L. Guerrieri

status = 1;
constraint_nbr = sum(strcmp(upper(M.equations_tags(:,2)),'OCCBIN'))/2;

switch(constraint_nbr)
  case 1
    [zdatalinear_ zdatapiecewise_ zdatass_ oobase_ ] = ...
    solve_one_constraint(M,oo,options);
  case 2
    [zdatalinear_ zdatapiecewise_ zdatass_ oobase_ ] = ...
    solve_two_constraints(M,oo,options);
  otherwise
    error('OCCBIN can only handle two constraints in a model')
end
endo_simul = zdatapiecewise_';
endo_simul_no_constraint = zdatalinear_';