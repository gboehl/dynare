% Solve model, generate model IRFs
[zdatalinear, zdatapiecewise, zdatass, oobase_, Mbase_ ] = ...
    solve_one_constraint(modnam,modnamstar,...
                         constraint, constraint_relax,...
                         shockssequence,irfshock,nperiods,maxiter);



% unpack the IRFs
for i=1:Mbase_.endo_nbr
    eval([Mbase_.endo_names{i}, '_uncdifference=zdatalinear(:,i);']);
    eval([Mbase_.endo_names{i}, '_difference=zdatapiecewise(:,i);']);
    eval([Mbase_.endo_names{i}, '_ss=zdatass(i);']);
end


nparams = length(Mbase_.param_names);

for i = 1:nparams
    eval([Mbase_.param_names{i}, '= Mbase_.params(i);']);
end
