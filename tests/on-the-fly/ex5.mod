// --+ options: nostrict +--

/*
** It is possible to associate an endogenous variable to an equation were the variable does not appear (if it makes sense).
*/


parameters a0, a1, b0, b1;

a0 = 1.0;
a1 = 0.1;
b0 = 0.5;
b1 = 0.2;



model;

[endogenous='D']
D = a0 - a1*p;

[endogenous='S']
S = b0 + b1*p;

[endogenous='p']
S = D;

end;

if ~isequal(length(intersect(M_.endo_names, {'D'; 'S'; 'p'})), 3)
   error('Endogenous variables are wrong.')
end

if isfield(M_, 'exo_names')
   error('Exogenous variables are wrong.')
end
