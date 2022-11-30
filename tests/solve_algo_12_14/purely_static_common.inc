@#define n = 10

@#for i in 1:n
var x@{i}, t@{i}, s@{i};
varexo e_x@{i}, e_t@{i}, e_s@{i};
@#endfor

model;
@#for i in 1:n
x@{i} = e_x@{i};
t@{i} = x@{i} + s@{i} + e_t@{i};
s@{i} = 3*t@{i} + e_s@{i};
@#endfor
end;

shocks;
@#for i in 1:n
var e_x@{i};
periods 1;
values @{0.1/n};

var e_t@{i};
periods 2;
values @{0.05+0.1/n};

var e_s@{i};
periods 3;
values @{0.1+0.1/n};
@#endfor
end;

perfect_foresight_setup(periods = 100);
