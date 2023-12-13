function [oo_,options_] = squeeze_shock_decomposition(M_,oo_,options_, sd_vlist)

if isstruct(options_.plot_shock_decomp.q2a)
    
    avname=char({options_.plot_shock_decomp.q2a.qname});
    sda = options_.plot_shock_decomp.q2a(ismember(avname,sd_vlist,'rows'));
    for k=1:length(sda)
        if isstruct(sda(k).aux)
             sd_vlist = [sd_vlist; cellstr(sda(k).aux.y)];
        end
    end
end
i_var = varlist_indices(sd_vlist,M_.endo_names);
sd_vlist = M_.endo_names(i_var);
% first we squeeze usual fields 
oo_ = squeeze_shock_decomposition(M_,oo_,options_,sd_vlist);
i_var = oo_.shock_decomposition_info.i_var;
sd_vlist = M_.endo_names(i_var);

% now we check for occbin SDs
options_.occbin.shock_decomp.i_var = i_var;
if isfield (oo_.occbin.smoother,'decomp')
    oo_.occbin.smoother.decomp = oo_.occbin.smoother.decomp(i_var,:,:);
    oo_.occbin.smoother.wdecomp = oo_.occbin.smoother.wdecomp(i_var,:,:);
end

if isfield(oo_.occbin,'shock_decomp')
    fnames = fieldnames(oo_.occbin.shock_decomp);
    for k=1:length(fnames)
        nendo = numel(oo_.occbin.shock_decomp.(fnames{k}).vname);
        tmp_i_var = varlist_indices(sd_vlist,char(oo_.occbin.shock_decomp.(fnames{k}).vname));
        oo_.occbin.shock_decomp.(fnames{k}).vname = cellstr(sd_vlist);
        tmpnames = fieldnames(oo_.occbin.shock_decomp.(fnames{k}));
        for t=1:length(tmpnames)
            if size(oo_.occbin.shock_decomp.(fnames{k}).(tmpnames{t}),3)==nendo
                oo_.occbin.shock_decomp.(fnames{k}).(tmpnames{t})= oo_.occbin.shock_decomp.(fnames{k}).(tmpnames{t})(:,:,tmp_i_var);
            end
        end
    end
end

end


