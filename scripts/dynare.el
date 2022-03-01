;;; dynare.el --- major mode for editing Dynare mod files

;; Copyright © 2010 Yannick Kalantzis
;; Copyright © 2019-2022 Dynare Team
;;
;; This program is free software; you can redistribute it and/or modify
;; it under the terms of the GNU General Public License as published by
;; the Free Software Foundation, either version 3 of the License, or
;; (at your option) any later version.

;; This program is distributed in the hope that it will be useful,
;; but WITHOUT ANY WARRANTY; without even the implied warranty of
;; MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
;; GNU General Public License for more details.
;;
;; You should have received a copy of the GNU General Public License
;; along with this program.  If not, see
;; <https://www.gnu.org/licenses/>.

;;; Installation:
;;
;;   Put this file somewhere on your load path. The mode will be automatically
;;   loaded and selected when you open a file with the "mod" extension.

;;; TODO:
;;    - font-locking:
;;      + external functions (font-lock-function-name-face),
;;      + change face for macroprocessor?
;;      + handle dates?
;;      + multi-line macro-commands/exprs
;;    - M-a and M-e should skip statements (separated with ;)
;;    - improve indentation
;;      + w.r.t. to statements/equations/macro-commands split on several lines
;;      + for macro-statements (insert space between @# and the keyword, with
;;        an option for controlling this offset)
;;      + add option for controlling whether macrocommands are indented at 0 or
;;        at same level as Dynare code
;;    - basically deactivate the mode within verbatim blocks?
;;    - blocks templates "model/end", "initval/end", etc.
;;    - functions to insert main keywords, with shortcuts in the keymap
;;    - make font-locking case insensitive (as is Dynare lexer)

(defgroup dynare nil
  "Editing Dynare mod files."
  :link '(url-link "https://www.dynare.org")
  :link '(custom-group-link :tag "Font Lock Faces group" font-lock-faces)
  :group 'languages)

(defcustom dynare-block-offset 2
  "Extra indentation applied to statements in Dynare block structures."
  :type 'integer)

;; Those keywords that makes the lexer enter the DYNARE_STATEMENT start
;; condition
;; Also include "end" in this list
(defvar dynare-statements
  '("var" "varexo" "varexo_det" "trend_var" "log_trend_var"
    "predetermined_variables" "parameters" "model_local_variable" "periods"
    "model_info" "estimation" "set_time" "data" "varobs"
    "varexobs" "unit_root_vars" "rplot" "osr_params" "osr" "dynatype"
    "dynasave" "model_comparison" "change_type" "load_params_and_steady_state"
    "save_params_and_steady_state" "write_latex_dynamic_model"
    "write_latex_static_model" "write_latex_original_model"
    "write_latex_steady_state_model" "steady" "check" "simul" "stoch_simul"
    "var_model" "trend_component_model" "var_expectation_model" "pac_model"
    "dsample" "Sigma_e" "planner_objective" "ramsey_model" "ramsey_policy"
    "evaluate_planner_objective" "occbin_setup" "occbin_solver"
    "occbin_write_regimes" "occbin_graph"
    "discretionary_policy" "identification" "bvar_density" "bvar_forecast"
    "dynare_sensitivity" "initval_file" "histval_file" "forecast"
    "shock_decomposition" "realtime_shock_decomposition"
    "plot_shock_decomposition" "initial_condition_decomposition"
    "squeeze_shock_decomposition" "sbvar"
    "ms_estimation" "ms_simulation" "ms_compute_mdd" "ms_compute_probabilities"
    "ms_forecast" "ms_irf" "ms_variance_decomposition" "conditional_forecast"
    "plot_conditional_forecast" "method_of_moments"
    "markov_switching" "svar" "svar_global_identification_check"
    "external_function" "calib_smoother" "model_diagnostics" "extended_path"
    "smoother2histval" "perfect_foresight_setup" "perfect_foresight_solver"
    "perfect_foresight_with_expectation_errors_setup"
    "perfect_foresight_with_expectation_errors_solver"
    "compilation_setup" "resid"
    "std" "corr" "prior_function" "posterior_function" "end")
  "Dynare statement keywords.")

;; Keywords that may appear in blocks, and that begin a statement which will be
;; closed by a semicolon
(defvar dynare-statements-like
  '("stderr" "values" "scales" "restriction" "exclusion" "upper_cholesky" "lower_cholesky"
    "bind" "relax" "error_bind" "error_relax")
  "Dynare statements-like keywords.")

;; Those keywords that makes the lexer enter the DYNARE_BLOCK start condition
;; Also include "verbatim" in this list
;; Needs to be enclosed within eval-when-compile, because this variable is
;; referenced in another eval-when-compile statement in dynare-calculate-indentation
(eval-when-compile
  (defvar dynare-blocks
    '("model" "steady_state_model" "initval" "endval" "histval" "shocks" "heteroskedastic_shocks"
      "shock_groups" "init2shocks" "mshocks" "estimated_params" "epilogue" "priors"
      "estimated_params_init" "estimated_params_bounds" "osr_params_bounds"
      "observation_trends" "deterministic_trends" "optim_weights" "homotopy_setup"
      "conditional_forecast_paths" "svar_identification" "moment_calibration"
      "irf_calibration" "ramsey_constraints" "generate_irfs"
      "matched_moments" "occbin_constraints" "verbatim")
    "Dynare block keywords."))

;; Mathematical functions and operators used in model equations (see "hand_side" in Bison file)
(defvar dynare-functions
  '("expectation" "var_expectation" "pac_expectation" "exp" "diff" "adl" "log"
    "ln" "log10" "sin" "cos" "tan" "asin" "acos" "atan" "sqrt" "cbrt" "abs"
    "sign" "max" "min" "normcdf" "normpdf" "erf" "steady_state")
  "Dynare mathematical functions and operators.")

(defvar dynare-constants
  '("nan" "inf")
  "Dynare constants.")

(defvar dynare-type-attributes
  '("|e" "|x" "|p")
  "Dynare attributes for on-the-fly type declarations.")

(defvar dynare-macro-keywords
  '("line" "include" "includepath" "define" "echo" "error" "if" "ifdef" "ifndef"
    "elseif" "else" "endif" "for" "endfor" "echomacrovars" "in" "when" "save"
    "true" "false" "exp" "log" "ln" "log10" "sin" "cos" "tan" "asin" "acos"
    "atan" "sqrt" "cbrt" "sign" "max" "min" "floor" "ceil" "trunc" "mod" "sum"
    "erf" "erfc" "gamma" "lgamma" "round" "length" "normpdf" "normcdf" "isempty"
    "isboolean" "isreal" "isstring" "istuple" "isarray" "bool" "real" "string"
    "tuple" "array" "defined" "nan" "inf")
  "Dynare macroprocessor keywords.")

(defvar dynare-font-lock-keywords
  `(("@#" . font-lock-variable-name-face) ; Beginning of macro-statement
    ("@#" ,(regexp-opt dynare-macro-keywords 'symbols)
          nil nil (0 font-lock-variable-name-face)) ; Keywords in macro-statements
    ("@{[^}]*}" . font-lock-variable-name-face) ; For macro-substitutions
    ;;; Below is an alternative way of dealing with macro-substitutions
    ;;; Only the delimiters and the keywords are colorized
    ;; ("@{" . font-lock-variable-name-face)
    ;; ("@{" "[^}]*\\(}\\)" nil nil (1 font-lock-variable-name-face))
    ;; ("@{" ,(concat (regexp-opt dynare-macro-keywords 'symbols) "[^}]*}") nil nil (1 font-lock-variable-name-face))
    ("^[ \t]*#" . font-lock-warning-face) ; For model-local variables
    (,(regexp-opt (append dynare-statements dynare-statements-like dynare-blocks) 'symbols) . font-lock-keyword-face)
    (,(regexp-opt dynare-functions 'symbols) . font-lock-builtin-face)
    (,(regexp-opt dynare-constants 'symbols) . font-lock-constant-face)
    (,(concat (regexp-opt dynare-type-attributes) "\\_>") . font-lock-type-face))
  "Keyword highlighting specification for `dynare-mode'.")

(defvar dynare-mode-map
  (let ((map (make-sparse-keymap)))
    ;; TODO: To be filled
    map))

(defvar dynare-mode-syntax-table
  (let ((st (make-syntax-table)))
    ;; mathematical operators are treated as punctuation
    ;; "*" is treated further below
    (modify-syntax-entry ?+ "." st)
    (modify-syntax-entry ?- "." st)
    (modify-syntax-entry ?/ "." st)
    (modify-syntax-entry ?^ "." st)

    ;; symbols for the macrolanguage
    (modify-syntax-entry ?@ "." st)
    (modify-syntax-entry ?# "." st)

    ;; underscores are symbol constituents
    (modify-syntax-entry ?_ "_" st)

    ;; Single-quoted strings
    (modify-syntax-entry ?' "\"" st)

    ;; define C++ style comment  “/* ... */” and “// ...”
    ;; "/" is the 1st and 2nd char of /* and */ (a-style) and the 2nd char of //
    ;; (b-style)
    (modify-syntax-entry ?\/ ". 124" st)
    ;; "*" is the 2nd and 1st char of /* and */ (a-style only)
    (modify-syntax-entry ?* ". 23b" st)
    ;; "%" starts a MATLAB-style comment
    (modify-syntax-entry ?% "<" st)
    ;; newline is the comment-end sequence of b-style and MATLAB-style comments
    (modify-syntax-entry ?\n ">" st)
    st)
  "Syntax table for `dynare-mode'")

(defun dynare-indent-line ()
  "Indent current line of Dynare mod file."
  (interactive)
  (let ((savep (> (current-column) (current-indentation)))
        (indent (max (dynare-calculate-indentation) 0)))
    (if savep
        (save-excursion (indent-line-to indent))
      (indent-line-to indent))))

(defun dynare-calculate-indentation ()
  "Return the column to which the current line should be indented."
  (save-excursion
    (beginning-of-line)
    (cond
      ((bobp)
       0)
      ((looking-at "^[ \t]*end[ \t]*;")
       ;; This is an "end" keyword: decrease the indentation level
       (forward-line -1)
       (max (- (current-indentation) dynare-block-offset)
            0))
      (t
       (let (cur-indent)
         (while (null cur-indent) ; Iterate backwards until we find an indentation hint
           (forward-line -1)
            (cond
              ((looking-at "^[ \t]*end[ \t]*;")
               ;; An "end" was matched: indent at the same level
               (setq cur-indent (current-indentation)))
              ((looking-at (concat "^[ \t]*" (eval-when-compile (regexp-opt
                                                                 dynare-blocks))
                                   "[ \t]*;"))
               ;; A block opening keyword was found: we need to indent an extra level
               (setq cur-indent (+ (current-indentation) dynare-block-offset))) ; Do the actual indenting
              ((bobp)
               ;; No hint was found: indent at 0
               (setq cur-indent 0))))
         cur-indent)))))

;;;###autoload
(define-derived-mode dynare-mode prog-mode "Dynare"
  "Major mode for editing Dynare mod files."
  :syntax-table dynare-mode-syntax-table

  (setq-local comment-start "// ") ; For comment-dwim
  (setq-local font-lock-defaults '(dynare-font-lock-keywords))
  (setq-local indent-line-function 'dynare-indent-line))

;;;###autoload
(add-to-list 'auto-mode-alist '("\\.mod$" . dynare-mode))

(provide 'dynare)

;; The following line is for ELPA compatibility (in particular needed for Debian)
;;; dynare.el ends here
