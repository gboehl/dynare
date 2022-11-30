/* A simple purely forward model, consisting of:
  - a forward VAR(3) with one variable in level, one in log and one log-forward-differentiated;
  - two other variables that depend on the contemporaneous values of the VAR(3),
    and which together for a simultaneous block.
  Those five equations are repeated an arbitrary number of times (n). */

@#define n = 10

@#for i in 1:n
var x@{i}, y@{i}, z@{i}, t@{i}, s@{i};
varexo e_x@{i}, e_y@{i}, e_z@{i}, e_t@{i}, e_s@{i};
@#endfor

model;
@#for i in 1:n
x@{i} = 0.1*x@{i}(+1) + 0.2*log(y@{i}(+1)) + 0.3*(log(z@{i})-log(z@{i}(+1))) + e_x@{i};
[name = 'y']
log(y@{i}) = 0.4*x@{i}(+1) + 0.5*log(y@{i}(+1)) + 0.6*(log(z@{i})-log(z@{i}(+1))) + e_y@{i};
[name = 'z']
/* NB: we choose 0.5 as steady state for z, since initial value is drawn from
   [0,1] (in the “simul_backward” case) */
log(z@{i})-log(z@{i}(+1)) = -0.8*(log(z@{i}(+1)) - log(0.5)) + 0.1*x@{i}(+1) + 0.2*log(y@{i}(+1)) - 0.3*(log(z@{i}) - log(z@{i}(+1))) + e_z@{i};

t@{i} = x@{i} + 2*log(y@{i}) + 3*(log(z@{i})-log(z@{i}(+1))) + s@{i} + e_t@{i};
s@{i} = 4*x@{i}(+1) + 3*t@{i} + e_s@{i};
@#endfor
end;

initval;
@#for i in 1:n
y@{i} = 1;
z@{i} = 0.5;
@#endfor
end;

shocks;
@#for i in 1:n
var e_x@{i};
periods 99;
values @{0.1/n};

var e_y@{i};
periods 98;
values @{0.05+0.1/n};

var e_z@{i};
periods 97;
values @{0.1+0.1/n};

var e_t@{i};
periods 96;
values @{0.15+0.1/n};

var e_s@{i};
periods 95;
values @{0.2+0.1/n};
@#endfor
end;

perfect_foresight_setup(periods = 100);
