function r = get_residuals(ivar,lb,ub,M,oo)

ss = oo.steady_state;
for i = 1:length(ivar)
    % only one is different from zero
    ss(ivar(i)) = lb(i) + ub(i);
end
oo.steady_state = ss;

r = evaluate_model(M,oo);