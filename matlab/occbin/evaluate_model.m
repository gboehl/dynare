function [r,g1,g2,g3] = evaluate_model(z,x,M,ss)

ll = M.lead_lag_incidence';
y = z(find(ll(:)));

switch nargout
  case 1
    r = feval([M.fname '_dynamic'],y,x, ...
              M.params, ss, 1);
  case 2
    [r,g1] = feval([M.fname '_dynamic'],y,x, ...
                    M.params, ss, 1);
  case 3
    [r,g1,g2] = feval([M.fname '_dynamic'],y,x, ...
                    M.params, ss, 1);
  case 4
    [r,g1,g2,g3] = feval([M.fname '_dynamic'],y,x, ...
                         M.params, ss, 1);
end