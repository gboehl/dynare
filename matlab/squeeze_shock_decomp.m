function [oo_,options_] = squeeze_shock_decomp(M_,oo_,options_,sd_vlist)

if nargin==3
    % automatic selection from history of plot_shock_decomp
    sd_vlist = M_.endo_names(options_.plot_shock_decomp.i_var);
end

if isstruct(options_.plot_shock_decomp.q2a)
    
    avname={options_.plot_shock_decomp.q2a.qname};
    sda = options_.plot_shock_decomp.q2a(ismember(avname,sd_vlist));
    for k=1:length(sda)
        if isstruct(sda(k).aux)
            sd_vlist = [sd_vlist; {sda(k).aux.y}];
        end
    end
end
i_var = varlist_indices(sd_vlist,M_.endo_names);

options_.shock_decomp.i_var = i_var;
oo_.shock_decomposition = oo_.shock_decomposition(i_var,:,:);

if isfield (oo_,'realtime_conditional_shock_decomposition')
    oo_.realtime_conditional_shock_decomposition = ...
        my_squeeze(oo_.realtime_conditional_shock_decomposition, i_var);
end
if isfield (oo_,'realtime_forecast_shock_decomposition')
    oo_.realtime_forecast_shock_decomposition = ...
        my_squeeze(oo_.realtime_forecast_shock_decomposition, i_var);
end
if isfield (oo_,'realtime_shock_decomposition')
    oo_.realtime_shock_decomposition = ...
        my_squeeze(oo_.realtime_shock_decomposition, i_var);
end
if isfield (oo_,'conditional_shock_decomposition')
    oo_.conditional_shock_decomposition = ...
        my_squeeze(oo_.conditional_shock_decomposition, i_var);
end
if isfield (oo_,'initval_decomposition')
    oo_.initval_decomposition = oo_.initval_decomposition(i_var,:,:);
end

end

function shock_decomposition = my_squeeze(shock_decomposition, i_var)
fnam = fieldnames(shock_decomposition);
for k=1:length(fnam)
    shock_decomposition.(fnam{k}) =  shock_decomposition.(fnam{k})(i_var,:,:);
end

end
