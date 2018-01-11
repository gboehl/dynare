% Script that retrieves parameter values once model is solved

nendog = length(Mbase_.endo_names);

for i=1:nendog
    eval([Mbase_.endo_names{i} '_ss = oo_.dr.ys(i); ']);
end

nparams = length(Mbase_.param_names);

for i = 1:nparams
    eval([Mbase_.param_names{i},'= M_.params(i);']);
end
