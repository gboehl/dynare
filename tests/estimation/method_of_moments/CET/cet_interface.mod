% ------------------------------------------------------------------------- 
% Functionality testing of interface for IRF matching
% ------------------------------------------------------------------------- 
 
% Copyright © 2023 Dynare Team 
% 
% This file is part of Dynare. 
% 
% Dynare is free software: you can redistribute it and/or modify 
% it under the terms of the GNU General Public License as published by 
% the Free Software Foundation, either version 3 of the License, or 
% (at your option) any later version. 
% 
% Dynare is distributed in the hope that it will be useful, 
% but WITHOUT ANY WARRANTY; without even the implied warranty of 
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the 
% GNU General Public License for more details. 
% 
% You should have received a copy of the GNU General Public License 
% along with Dynare.  If not, see <https://www.gnu.org/licenses/>. 
 
@#include "cet_model.inc"
xx = [23,24,25];
ww = [51,52];

matched_irfs(overwrite); 
var GDPAGG;  varexo epsR_eps;   periods 5;    values 7;     weights 25;
var GDPAGG;  varexo mupsi_eps;  periods 1,2;  values 17,18;  weights 37,38;
var RAGG;    varexo muz_eps;    periods 3:5;  values (xx);
var ukAGG;   varexo mupsi_eps;  periods 1:2;  values 30;     weights (ww);

varexo epsR_eps;
var wAGG;
periods 1, 13:15, 2:12;
values 2, (xx), 15;
weights 3, (xx), 4;
end;

method_of_moments(mom_method = irf_matching
, cova_compute=0
, irf_matching_file = cet_irf_matching_file
, mcmc_jumping_covariance = prior_variance
, mh_jscale = 0.5
, mh_nblocks = 1
, mh_replic=1000
, mode_compute = 0
, mode_file = cet_original_mode
, plot_priors = 0
, posterior_sampling_method = 'independent_metropolis_hastings'
);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% CHECKS ON INTERFACE AND TRANSFORMATIONS %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
verbatim;

% check on matched_irfs
matched_irfs = [
    {'GDPAGG'} {'epsR_eps'}   {[ {[5] [ 7] [25]} ]};
    {'GDPAGG'} {'mupsi_eps'}  {[ {[1] [17] [37]};
                                 {[2] [18] [38]} ]};
    {'RAGG'}   {'muz_eps'}    {[ {[3 4 5] [23 24 25] [1]} ]};
    {'ukAGG'}  {'mupsi_eps'}  {[ {[1 2] [30] [51 52]} ]};
    {'wAGG'} {'epsR_eps'}     {[ {[1]     [2]        [3]};
                                 {[13:15] [23 24 25] [23 24 25]};
                                 {[2:12]  [15]       [4]}; ]};
];
if ~isequal(M_.matched_irfs,matched_irfs)
    error('Something wrong with the transformation of the matched_irfs block!')
else
    fprintf('matched_irfs transformation was successful!\n\n')
end

% check on data_moments
data_moments = [7 2 15 15 15 15 15 15 15 15 15 15 15 23 24 25 23 24 25 17 18 30 30]';
if ~isequal(oo_.mom.data_moments,data_moments)
    error('Something wrong with the creation of data_moments!')
else
    fprintf('creation of data_moments was successful!\n\n')
end

% check on weighting matrix
weighting_mat = sparse(size(data_moments,1));
weighting_mat(1,1)   = 25;
weighting_mat(2,2)   =  3;
weighting_mat(3,3)   =  4;
weighting_mat(4,4)   =  4;
weighting_mat(5,5)   =  4;
weighting_mat(6,6)   =  4;
weighting_mat(7,7)   =  4;
weighting_mat(8,8)   =  4;
weighting_mat(9,9)   =  4;
weighting_mat(10,10) =  4;
weighting_mat(11,11) =  4;
weighting_mat(12,12) =  4;
weighting_mat(13,13) =  4;
weighting_mat(14,14) = 23;
weighting_mat(15,15) = 24;
weighting_mat(16,16) = 25;
weighting_mat(17,17) =  1;
weighting_mat(18,18) =  1;
weighting_mat(19,19) =  1;
weighting_mat(20,20) = 37;
weighting_mat(21,21) = 38;
weighting_mat(22,22) = 51;
weighting_mat(23,23) = 52;
if ~isequal(oo_.mom.weighting_info.W,weighting_mat)
    error('Something wrong with the creation of weighting_info.W!')
else
    fprintf('creation of weighting_info.W was successful!\n\n')
end

end;//verbatim


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% SIMPLIFIED EXAMPLE TO TEST INTERFACE ON WEIGHTING MATRIX %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

matched_irfs(overwrite); 
var GDPAGG;  varexo epsR_eps;   periods 1 3;  values 11 13; weights 111 222;
var RAGG;    varexo muz_eps;    periods 2 4;  values 22 24; weights 333 444;
end;

matched_irfs_weights(overwrite);
RAGG(2),   muz_eps,  GDPAGG(1), epsR_eps, 555;
GDPAGG(1), epsR_eps, RAGG(4),   muz_eps,  666;
RAGG(2),   muz_eps,  GDPAGG(3), epsR_eps, 777;
GDPAGG(3), epsR_eps, RAGG(4),   muz_eps,  888;
GDPAGG(1), epsR_eps, RAGG(1),   muz_eps,  999;
end;

verbatim;
% check on matched_irfs_weights
matched_irfs_weights = [
    {'GDPAGG'}    {[1]}    {'epsR_eps'}    {'RAGG'  }    {[1]}    {'muz_eps' }    {[999]};
    {'GDPAGG'}    {[1]}    {'epsR_eps'}    {'RAGG'  }    {[4]}    {'muz_eps' }    {[666]};
    {'GDPAGG'}    {[3]}    {'epsR_eps'}    {'RAGG'  }    {[4]}    {'muz_eps' }    {[888]};
    {'RAGG'  }    {[2]}    {'muz_eps' }    {'GDPAGG'}    {[1]}    {'epsR_eps'}    {[555]};
    {'RAGG'  }    {[2]}    {'muz_eps' }    {'GDPAGG'}    {[3]}    {'epsR_eps'}    {[777]};
];
if ~isequal(M_.matched_irfs_weights,matched_irfs_weights)
    error('Something wrong with the transformation of the matched_irfs_weights block!')
else
    fprintf('matched_irfs_weights transformation was successful!\n\n')
end

% test transformation function
[data_moments, W, irfIndex] = mom.matched_irfs_blocks(M_.matched_irfs, M_.matched_irfs_weights, options_mom_.varobs_id, length(options_mom_.varobs_id), M_.exo_nbr, M_.endo_names, M_.exo_names);

if ~isequal(data_moments, [11 13 22 24]')
    error('Something wrong with the creation of data_moments in simple example!');
else
    fprintf('simple example: creation of data_moments was successful!\n\n');
end

if ~isequal(irfIndex, [1 3 58 60]')
    error('Something wrong with the creation of irfIndex in simple example!');
else
    fprintf('simple example: creation of irfIndex was successful!\n\n')
end

W0 = eye(4);
W0(1,1) = 111; % (GDPAGG(1) epsR_eps) x (GDPAGG(1) epsR_eps)
W0(2,2) = 222; % (GDPAGG(3) epsR_eps) x (GDPAGG(3) epsR_eps)
W0(3,3) = 333; % (RAGG(2) muz_eps) x (RAGG(2) muz_eps)
W0(4,4) = 444; % (RAGG(4) muz_eps) x (RAGG(4) muz_eps)
W0(1,3) = 555; % (GDPAGG(1) epsR_eps) x (RAGG(2) muz_eps)
W0(3,1) = 555; % (RAGG(2) muz_eps) x (GDPAGG(1) epsR_eps)
W0(1,4) = 666; % (GDPAGG(1) epsR_eps) x (RAGG(4) muz_eps)
W0(4,1) = 666; % (RAGG(4) muz_eps) x (GDPAGG(1) epsR_eps)
W0(2,3) = 777; % (GDPAGG(3) epsR_eps) x (RAGG(2) muz_eps)
W0(3,2) = 777; % (RAGG(2) muz_eps) x (GDPAGG(3) epsR_eps)
W0(2,4) = 888; % (GDPAGG(3) epsR_eps) x (RAGG(4) muz_eps)
W0(4,2) = 888; % (RAGG(4) muz_eps) x (GDPAGG(3) epsR_eps)

if ~isequal(W, W0)
    error('Something wrong with the creation of W in simple example!')
else
    fprintf('simple example: creation of W was successful!\n\n')
end
end;