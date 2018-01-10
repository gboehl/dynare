// --+ options: json=compute +--
% Example from Section 6.1 of
% Ando, Tomohiro and Zellner, Arnold. 2010. Hierarchical Bayesian Analysis of the
% Seemingly Unrelated Regression and Simultaneous Equations Models Using a
% Combination of Direct Monte Carlo and Importance Sampling Techniques.
% Bayesian Analysis Volume 5, Number 1, pp. 65-96.

var qty, price;
varexo res_u, res_v, stormy, mixed;
parameters bq0, bp0, bq1, bp1, bp2;

bp0 =  8.5527;
bp1 = -0.5302;
bp2 = -0.3974;

bq0 =  6.7523;
bq1 = -0.7969;

model(linear);
    price = bp0 + bp1*stormy + bp2*mixed + res_v;
    qty = bq0 + bq1*price + res_u;
end;

surgibbs(dseries('fishdata.csv'), 0.0005.*eye(M_.param_nbr), 10000, 5000);
