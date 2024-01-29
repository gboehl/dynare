# Copyright © 2018-2021 Dynare Team
#
# This file is part of Dynare.
#
# Dynare is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# Dynare is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with Dynare.  If not, see <https://www.gnu.org/licenses/>.

import re

from pygments.lexer import Lexer, RegexLexer, bygroups, do_insertions, words
from pygments.token import Text, Comment, Operator, Keyword, Name, String, \
    Number, Punctuation, Generic, Whitespace

# Configuration block :BOLD
#conf_blocks = ('hooks','paths','cluster','node')

__all__ = ['DynareLexer']

class DynareLexer(RegexLexer):

    name = 'Dynare'
    aliases = ['dynare']
    filenames = ['*.mod']

    commands = (
        "dynare","var","varexo","varexo_det","parameters","change_type","model_local_variable",
        "predetermined_variables","trend_var","log_trend_var","external_function",
        "write_latex_original_model","write_latex_dynamic_model",
        "write_latex_static_model","write_latex_steady_state_model","resid","initval_file","histval_file","dsample",
        "periods","values","scales","corr","stderr","steady","check","model_diagnostics","model_info",
        "print_bytecode_dynamic_model"," print_bytecode_static_model",
        "perfect_foresight_setup","perfect_foresight_solver","simul","stoch_simul",
        "extended_path","varobs","estimation","unit_root_vars","bvar_density",
        "model_comparison","shock_decomposition","realtime_shock_decomposition",
        "plot_shock_decomposition","calib_smoother","forecast",
        "conditional_forecast","plot_conditional_forecast","bvar_forecast",
        "smoother2histval","osr","osr_params","ramsey_model","ramsey_policy",
        "discretionary_policy","planner_objective","dynare_sensitivity",
        "markov_switching","svar","sbvar","ms_estimation","ms_simulation",
        "ms_compute_mdd","ms_compute_probabilities","ms_irf","ms_forecast",
        "ms_variance_decomposition","rplot","dynatype","dynasave","set_dynare_seed",
        "save_params_and_steady_state","load_params_and_steady_state",
        "dynare_version","write_latex_definitions","write_latex_parameter_table",
        "write_latex_prior_table","collect_latex_files","prior_function",
        "posterior_function","generate_trace_plots","evaluate_planner_objective",
        "occbin_setup","occbin_solver","occbin_write_regimes","occbin_graph","method_of_moments",
        "var_model","trend_component_model","var_expectation_model","pac_model")

    report_commands = ("report","addPage","addSection","addGraph","addTable",
        "addSeries","addParagraph","addVspace","write","compile")

    operators = (
        "STEADY_STATE","EXPECTATION","var_expectation","pac_expectation","pac_target_nonstationary")

    macro_dirs = (
        "@#includepath", "@#include", "@#define", "@#if",
        "@#ifdef", "@#ifndef", "@#else","@#endif",
        "@#for", "@#endfor", "@#echo", "@#error")

    builtin_constants = (
        "inf", "nan")

    tokens = {
        'root': [
            (r'\s*(%|//).*$', Comment.Single),
            (r'/(\\\n)?[*][\w\W]*?[*](\\\n)?/', Comment.Multiline),

            (words((
                'model','steady_state_model','initval','endval','histval','epilogue',
                'shocks','mshocks','homotopy_setup','observation_trends',
                'estimated_params','estimated_params_init','estimated_params_bounds',
                'shock_groups','conditional_forecast_paths','optim_weights',
                'osr_params_bounds','ramsey_constraints','irf_calibration',
                'moment_calibration','identification','svar_identification',
                'matched_moments','occbin_constraints','surprise','overwrite','bind','relax',
                'verbatim','end','node','cluster','paths','hooks','target','pac_target_info','auxname_target_nonstationary',
                'component', 'growth', 'auxname', 'kind'), prefix=r'\b', suffix=r'\s*\b'),Keyword.Reserved),

            # FIXME: Commands following multiline comments are not highlighted properly.
            (words(commands + report_commands,
                   prefix=r'\b', suffix=r'\s*\b'),  Name.Entity),

            (words(operators), Operator.Word),

            (words(macro_dirs,suffix=r'\s*'), Name.Function),

            (words(builtin_constants), Name.Constant),

            (r'\s*[a-zA-Z_]\s*', Name),

            (r'\s*(\d+\.\d+|\d*\.\d+)([eEf][+-]?[0-9]+)?\s*', Number.Float),
            (r'\s*\d+[eEf][+-]?[0-9]+\s*', Number.Float),
            (r'\s*\d+\s*', Number.Integer),

            (r'"[^"]*"', String),
            (r"`[^`]*'", String),
            (r"'[^']*'", String),

            (r'\s*(-|\+|\*|\/|\^)\s*', Operator),
            (r'\s*(==|<=|>=|~=|<|>|&&|!)\s*', Operator),

            (r'\s*[\[\](){}:@.,\|]\s*', Punctuation),
            (r'\s*(=|:|;|>>|#|\$)\s*', Punctuation),
        ]

    }
