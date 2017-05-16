function title=add_filter_subtitle(title,options_)

if ~options_.hp_filter && ~options_.one_sided_hp_filter  && ~options_.bandpass.indicator %do not filter
                                                                                         %nothing to add here
elseif ~options_.hp_filter && ~options_.one_sided_hp_filter && options_.bandpass.indicator
    title = [title ' (Bandpass filter, (' ...
             num2str(options_.bandpass.passband(1)),' ',num2str(options_.bandpass.passband(2)), '))'];
elseif options_.hp_filter && ~options_.one_sided_hp_filter  && ~options_.bandpass.indicator %filter with HP-filter
    title = [title ' (HP filter, lambda = ' ...
             num2str(options_.hp_filter) ')'];
elseif ~options_.hp_filter && options_.one_sided_hp_filter && ~options_.bandpass.indicator
    title = [title ' (One-sided HP filter, lambda = ' ...
             num2str(options_.one_sided_hp_filter) ')'];
end
end