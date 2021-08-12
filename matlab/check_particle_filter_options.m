function [particle_options] = check_particle_filter_options(particle_options)
% function [particle_filter_options, options_] = check_particle_filter_options(particle_filter_options_string, options_)
% initialization of particle filter options
%
% INPUTS
%   particle_filter_options:                structure storing the options

% OUTPUTS
%   particle_filter_options:                checked particle filter options
%
% SPECIAL REQUIREMENTS
%   none

% Copyright (C) 2021 Dynare Team
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
% along with Dynare.  If not, see <http://www.gnu.org/licenses/>.

if strcmpi(particle_options.filter_algorithm, 'sis')
    particle_options.algorithm = 'sequential_importance_particle_filter';
elseif strcmpi(particle_options.filter_algorithm, 'apf')
    particle_options.algorithm = 'auxiliary_particle_filter';
elseif strcmpi(particle_options.filter_algorithm, 'gf')
    particle_options.algorithm = 'gaussian_filter';
elseif strcmpi(particle_options.filter_algorithm,  'gmf')
    particle_options.algorithm = 'gaussian_mixture_filter';
elseif strcmpi(particle_options.filter_algorithm, 'cpf')
    particle_options.algorithm = 'conditional_particle_filter';
elseif strcmpi(particle_options.filter_algorithm, 'nlkf')
    particle_options.algorithm = 'nonlinear_kalman_filter';
else
    error(['Estimation: Unknown filter ' particle_options.filter_algorithm])
end

if ~isempty(particle_options.particle_filter_options)
    % set default options and user defined options
    options_list = read_key_value_string(particle_options.particle_filter_options);
    for i=1:rows(options_list)
        switch options_list{i,1}
            case 'posterior_sampler'
                if ~(strcmpi(options_list{i,2}, 'Herbst_Schorfheide') || ...
                        strcmpi(options_list{i,2}, 'DSMH'))
                    error(['check_particle_filter_options:: the proposal_distribution option to estimation takes either ' ...
                        'Herbst_Schorfheide or Herbst_Schorfheide as options']);
                else
                    particle_options.posterior_sampler=options_list{i,2};
                end
            case 'initial_state_prior_std'
                if options_list{i,2} <= 0
                    error('check_particle_filter_options:: the initial_state_prior_std option takes a positive argument');
                else
                    particle_options.initial_state_prior_std=options_list{i,2};
                end
            case 'pruning'
                if ~islogical(options_list{i,2})
                    error('check_particle_filter_options:: the pruning options takes only true or false');
                else
                    particle_options.pruning=options_list{i,2};
                end
            case 'unscented_alpha'
                if options_list{i,2} <= 0
                    error('check_particle_filter_options:: the unscented_alpha option takes a positive argument');
                else
                    particle_options.unscented.alpha=options_list{i,2};
                end
            case 'unscented_beta'
                if options_list{i,2} <= 0
                    error('check_particle_filter_options:: the unscented_beta option takes a positive argument');
                else
                    particle_options.unscented.beta=options_list{i,2};
                end
            case 'unscented_kappa'
                if options_list{i,2} <= 0
                    error('check_particle_filter_options:: the unscented_kappa option takes a positive argument');
                else
                    particle_options.unscented.kappa=options_list{i,2};
                end
            case 'mixture_state_variables'
                if options_list{i,2} <= 0
                    error('check_particle_filter_options:: the mixture_state_variables option takes a positive integer');
                else
                    particle_options.mixture_state_variables=options_list{i,2};
                end
            case 'mixture_structural_shocks'
                if options_list{i,2} <= 0
                    error('check_particle_filter_options:: the mixture_structural_shocks option takes a positive integer');
                else
                    particle_options.mixture_structural_shocks=options_list{i,2};
                end
            case 'mixture_measurement_shocks'
                if options_list{i,2} <= 0
                    error('check_particle_filter_options:: the mixture_measurement_shocks option takes a positive integer');
                else
                    particle_options.mixture_measurement_shocks=options_list{i,2};
                end
            case 'liu_west_delta'
                if options_list{i,2} <= 0
                    error('check_particle_filter_options:: the liu_west_delta option takes a positive argument');
                else
                    particle_options.liu_west_delta=options_list{i,2};
                end
                
            otherwise
                warning(['check_particle_filter_options: Unknown option (' options_list{i,1} ')!'])
        end
    end
end
