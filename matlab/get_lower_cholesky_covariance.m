function chol_sigma=get_lower_cholesky_covariance(Sigma_e,add_tiny_number_to_cholesky)
% function chol_sigma=get_lower_cholesky_covariance(Sigma_e)
% Computes the lower triangular Cholesky decomposition of a covariance matrix,
% working around zero entries on the diagonal and perfect correlation
%
% INPUTS
%   Sigma_e         [double]      covariance matrix
%
% OUTPUTS
%   chol_sigma      [cell]        Cholesky factor 
%
% ALGORITHM
%   Add small value to diagonal to break perfect correlation
%
% SPECIAL REQUIREMENTS.
%   None.
%
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

if nargin<2
    add_tiny_number_to_cholesky=1e-14;
end
std_deviation=sqrt(diag(Sigma_e));
non_zero_indices=find(std_deviation~=0); %find non-zero shocks;
try    
    chol_sigma=zeros(size(Sigma_e));
    chol_sigma(non_zero_indices,non_zero_indices)=chol(Sigma_e(non_zero_indices,non_zero_indices),'lower');
catch
    % cases with perfect correlation
    fprintf('Non-positive definite covariance matrix encountered. Using add_tiny_number_to_cholesky one the diagonal.\n')
    chol_sigma=zeros(size(Sigma_e));
    chol_sigma(non_zero_indices,non_zero_indices)=chol(Sigma_e(non_zero_indices,non_zero_indices)+add_tiny_number_to_cholesky*eye(length(non_zero_indices)),'lower');
    % correlation=diag(std_deviation(non_zero_indices))\Sigma_e(non_zero_indices,non_zero_indices)/diag(std_deviation(non_zero_indices));
end

return % --*-- Unit tests --*--

%@test:1

Sigma_e=diag(4*ones(3,1));
Sigma_e(2,2)=0;
chol_1=get_lower_cholesky_covariance(Sigma_e);
if max(max(abs(chol_1-diag([2,0,2]))))>eps 
    t(1)=false;
else
    t(1)=true;
end

Sigma_e=ones(3,3);
chol_2=get_lower_cholesky_covariance(Sigma_e,1e-14);
chol_3=get_lower_cholesky_covariance(Sigma_e+1e-14*eye(3),1e-14);
if max(max(abs(chol_2-chol_3)))>eps || any(any(triu(chol_3,1)))
    t(2)=false;
else
    t(2)=true;
end

Sigma_e=ones(3,3);
Sigma_e(2,:)=0;
Sigma_e(:,2)=0;
chol_4=get_lower_cholesky_covariance(Sigma_e,1e-14);
if chol_4(2,2)~=0 || any(any(triu(chol_4,1)))
    t(3)=false;
else
    t(3)=true;
end   

Sigma_e=[4 0.5 0; 0.5 9 0; 0 0 16];
chol_5=get_lower_cholesky_covariance(Sigma_e,1e-14);
if any(any(triu(chol_5,1))) %should be lower triangular
    t(4)=false;
else
    t(4)=true;
end

T = all(t);
%@eof:1
