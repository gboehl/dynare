function [r, J] = jacobian(fun, x, h)

% Returns the Jacobian matrix associated to a function.
%
% INPUTS:
% - fun    [handle]   Function from Rᵖ to Rⁿ.
% - x      [double]   1×p vector, point where the Jacobian has to be evaluated.
% - h      [double]   scalar, perturbation size.
%
% OUTPUTS:
% - r      [double]   n×1 vector, value of the function at point x.
% - J      [double]   n×p matrix of derivatives.

r = fun(x);

if nargout>1
    z = x;
    n = length(r);
    p = length(x);
    J = zeros(n, p);
    for i=1:p
        x = z;
        x(i) = x(i)+h;
        J(:,i) = (r-fun(x))/h;
    end
end