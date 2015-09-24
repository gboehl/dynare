function [fval,DLIK,Hess,exit_flag] = objective_function_penalty(x0,fcn,penalty,varargin)
    [fval,DLIK,Hess,exit_flag,SteadyState,trend_coeff,info] = fcn(x0,varargin{:});

    
    
    if info(1) ~= 0
        fval = penalty + info(2);
    end
end
