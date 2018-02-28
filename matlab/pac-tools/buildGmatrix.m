function G = buildGmatrix(alpha, beta)
    
% Builds the G matrix needed for PAC.
%
% INPUTS 
% - alpha    [double]    m*1 vector of PAC parameters (lag polynomial parameters).
% - beta     [double]    scalar, discount factor.
%
% OUTPUTS 
% - G         [double]    (m+1)*(m+1) matrix.

if nargin<2
    error('Two input arguments are required (vector of alpha parameters and beta discount parameter)!')
end

% Return an error if the first input is not a real vector.
if ~isnumeric(alpha) || ~isreal(alpha) || ~isvector(alpha)
    error('First input argument has to be a vector of doubles!')
end 

% Return an error if the second input argument is not a discount factor
if ~isnumeric(beta) || ~isreal(beta) || ~isscalar(beta) || beta<eps || beta>1-eps
    error('Second input argument has to be a discount factor!')
end 

% Get the number of parameters
m = length(alpha);

% Return an error if params is too small.
if m<1
    error('First input argument has to be a vector with at least one element.')
end

% Initialize the returned G matrix.
G = zeros(m);

% Fill the returned G matrix.
G(1:m-1,2:m) = eye(m-1);
G(m, :) = -flip(transpose(alpha));
G(m, :) = G(m, :).*flip(cumprod(beta*ones(1,m)));