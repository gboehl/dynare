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
    qty = bq0 + bq1*price + res_u;
    price = bp0 + bp1*stormy + bp2*mixed + res_v;
end;

% Estimate all parameters
%estparams = M_.param_names;
%estparamsval = M_.params;

% Estimate demand parameters
estparams = {'bq1' 'bq0'};
estparamsval = [bq1 bq0];

A = 0.0005.*eye(length(estparams));
surgibbs(dseries('fishdata.csv'), estparams, estparamsval, A, 20000, 5000, 7);

good = [6.791587808530124
   8.552700000000000
  -0.478275288902356
  -0.530200000000000
  -0.397400000000000];

if sum(abs(M_.params-good)) > 1e-14
    error(['sum of M_.params - good was: ' num2str(sum(abs(M_.params-good)))]);
end
