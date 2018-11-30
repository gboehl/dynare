function [coef_y,coef_u] = get_coef(jacobian,M)

ll = M.lead_lag_incidence;
endo_nbr = M.endo_nbr;
coef_y = zeros(endo_nbr,3*endo_nbr);
coef_u = zeros(endo_nbr,M.exo_nbr);

if M.maximum_lag > 0
    [~,c1,c2] = find(ll(1,:));
    coef_y(:,c1) = jacobian(:,c2);
    [~,c1,c2] = find(ll(2,:));
    coef_y(:,c1+endo_nbr) = jacobian(:,c2);
    if M.maximum_lead > 0
            [~,c1,c2] = find(ll(3,:));
            coef_y(:,c1+2*endo_nbr) = jacobian(:,c2);
    end
else
    [~,c1,c2] = find(ll(1,:));
    coef_y(:,c1+endo_nbr) = jacobian(:,c2);
    if M.maximum_lead > 0
            [~,c1,c2] = find(ll(2,:));
            coef_y(:,c1+2*endo_nbr) = jacobian(:,c2);
    end
end

coef_u = jacobian(:,max(ll(end,:))+1:end);