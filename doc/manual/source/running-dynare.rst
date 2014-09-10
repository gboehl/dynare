.. _running-dynare:

****************
 Running Dynare
****************

^^^^^^^^^^^^^^^^^
Dynare invocation
^^^^^^^^^^^^^^^^^

In order to give instructions to Dynare, the user has to write a model
file whose filename extension must be '.mod'. This file contains the
description of the model and the computing tasks required by the
user. Its contents is described in :ref:`The model file <the-model-file>` section of this reference manual.

Once the model file is written, Dynare is invoked using the ``dynare``
command at the MATLAB or Octave prompt (with the filename of the
'.mod' given as argument).

In practice, the handling of the model file is done in two steps: in
the first one, the model and the processing instructions written by
the user in a model file are interpreted and the proper MATLAB or GNU
Octave instructions are generated; in the second step, the program
actually runs the computations. Both steps are triggered automatically
by the ``dynare`` command.

.. function:: dynare FILENAME[.mod] [OPTIONS…]

   **Description**

   This command launches Dynare and executes the instructions included
   in 'FILENAME.mod'. This user-supplied file contains the model and
   the processing instructions, as described in :ref:`The model file
   <the-model-file>` section. The ``dynare`` command begins by
   launching the preprocessor on the '.mod' file. By default (unless
   ``use_dll`` or ``bytecode`` options have been given to the
   ``model`` block), the preprocessor creates three intermediary
   files:

   * 'FILENAME.m' Contains variable declarations, and computing tasks

   * 'FILENAME_dynamic.m' Contains the dynamic model equations. Note
     that Dynare might introduce auxiliary equations and variables
     (see section Auxiliary variables). Outputs are the residuals of
     the dynamic model equations in the order the equations were
     declared and the Jacobian of the dynamic model equations. For
     higher order approximations also the Hessian and the third-order
     derivatives are provided. When computing the Jacobian of the
     dynamic model, the order of the endogenous variables in the
     columns is stored in ``M_.lead_lag_incidence``. The rows of this
     matrix represent time periods: the first row denotes a lagged
     (time t-1) variable, the second row a contemporaneous (time t)
     variable, and the third row a leaded (time t+1) variable. The
     columns of the matrix represent the endogenous variables in their
     order of declaration. A zero in the matrix means that this
     endogenous does not appear in the model in this time period. The
     value in the ``M_.lead_lag_incidence`` matrix corresponds to the
     column of that variable in the Jacobian of the dynamic
     model. Example: Let the second declared variable be c and the
     (3,2) entry of ``M_.lead_lag_incidence`` be 15. Then the 15th
     column of the Jacobian is the derivative with respect to y(+1).

   * 'FILENAME_static.m' Contains the long run static model
     equations. Note that Dynare might introduce auxiliary equations
     and variables (see section Auxiliary variables). Outputs are the
     residuals of the static model equations in the order the
     equations were declared and the Jacobian of the static
     equations. Entry (i,j) of the Jacobian represents the derivative
     of the ith static model equation with respect to the jth model
     variable in declaration order.

   These files may be looked at to understand errors reported at the
   simulation stage. Dynare will then run the computing tasks by
   executing 'FILENAME.m'.

   A few words of warning is warranted here: the filename of the
   '.mod' file should be chosen in such a way that the generated
   '.m' files described above do not conflict with '.m' files
   provided by MATLAB/Octave or by Dynare. Not respecting this rule
   could cause crashes or unexpected behaviour. In particular, it
   means that the '.mod' file cannot be given the name of a
   MATLAB/Octave or Dynare command. Under Octave, it also means that
   the '.mod' file cannot be named 'test.mod'.

   **Options**

   * ``noclearall``

     By default, dynare will issue a clear all command to MATLAB or
     Octave, thereby deleting all workspace variables; this options
     instructs dynare not to clear the workspace.

   * ``debug``

    Instructs the preprocessor to write some debugging information
    about the scanning and parsing of the '.mod' file.

   * ``notmpterms``

    Instructs the preprocessor to omit temporary terms in the static
    and dynamic files; this generally decreases performance, but is
    used for debugging purposes since it makes the static and dynamic
    files more readable.

   * ``savemacro`` [=FILENAME]

    Instructs dynare to save the intermediary file which is obtained
    after macro-processing (see section :ref:`Macro processing
    language <the-model-file_macro-processing-language>`); the saved
    output will go in the file specified, or if no file is specified
    in 'FILENAME-macroexp.mod'

   * ``onlymacro``

    Instructs the preprocessor to only perform the macro-processing
    step, and stop just after. Mainly useful for debugging purposes or
    for using the macro-processor independently of the rest of Dynare
    toolbox.

   * ``nolinemacro``

    Instructs the macro-preprocessor to omit line numbering
    information in the intermediary '.mod' file created after the
    macro-processing step. Useful in conjunction with savemacro when
    one wants that to reuse the intermediary '.mod' file, without
    having it cluttered by line numbering directives.

   * ``nolog``

    Instructs Dynare to no create a logfile of this run in
    'FILENAME.log'. The default is to create the logfile.

   * ``nowarn``

    Suppresses all warnings.

   * ``warn_uninit``

    Display a warning for each variable or parameter which is not
    initialized. See section :ref:`Parameter
    initialization<the-model-file_parameter-initialization>`, or
    ``load_params_and_steady_state`` for initialization of
    parameters. See section :ref:`Initial and terminal
    conditions<the-model-file_initial-and-terminal-conditions>`,
    or ``load_params_and_steady_state`` for initialization of
    endogenous and exogenous variables.

   * ``console``

    Activate console mode. In addition to the behavior of
    ``nodisplay``, Dynare will not use graphical waitbars for long
    computations.

   * ``nograph``

    Activate the nograph option (see nograph), so that Dynare will not produce any graph.

   * ``nointeractive``

    Instructs Dynare to not request user input.

   * ``cygwin``

    Tells Dynare that your MATLAB is configured for compiling MEX
    files with Cygwin (see section Software requirements). This option
    is only available under Windows, and is used in conjunction with
    the ``use_dll`` option in the model block.

   * ``msvc``

    Tells Dynare that your MATLAB is configured for compiling MEX
    files with Microsoft Visual C++ (see section Software
    requirements). This option is only available under Windows, and is
    used in conjunction with the ``use_dll`` option in the model
    block.

   * ``parallel`` [=CLUSTER_NAME]

    Tells Dynare to perform computations in parallel. If CLUSTER_NAME
    is passed, Dynare will use the specified cluster to perform
    parallel computations. Otherwise, Dynare will use the first
    cluster specified in the configuration file. See section :ref:`The
    Configuration File<the-configuration-file>`, for more information.

   * ``conffile`` =FILENAME

    Specifies the location of the configuration file if it differs
    from the default. See section :ref:`The Configuration
    File<the-configuration-file>`, for more information about the
    configuration file and its default location.

   * ``parallel_slave_open_mode``

    Instructs Dynare to leave the connection to the slave node open
    after computation is complete, closing this connection only when
    Dynare finishes processing.

   * ``parallel_test``

    Tests the parallel setup specified in the configuration file
    without executing the .mod file. See section :ref:`The
    Configuration File<the-configuration-file>`, for more information.

   * ``-DMACRO_VARIABLE`` =MACRO_EXPRESSION

    Defines a macro-variable from the command line (the same effect as
    using the Macro directive ``@#define`` in a model file, see
    section :ref:`Macro processing
    language<the-model-file_macro-processing-language>`).

   * ``nostrict``

     Allows Dynare to issue a warning and continue processing when
     there are more endogenous variables than equations, an undeclared
     symbol is assigned in ``initval`` or ``endval``.


   **Output**

     Depending on the computing tasks requested in the .mod file, executing
     the dynare command will leave variables containing results in the
     workspace available for further processing. More details are given under
     the relevant computing tasks.

     The ``M_``, ``oo_``, and ``options_`` structures are saved in a file called
     FILENAME_results.mat. If they exist, ``estim_params_``, ``bayestopt_``, ``dataset_``,
     and ``estimation_info`` are saved in the same file.


   **Examples**

     .. code-block:: none

        >> dynare ramst
        >> dynare ramst nolog

     Dynare will fail to run a .mod file if this file is not in the current
     directory. For this example, assuming that Dynare is installed under
     'c:\\dynare\\4.x.y' (default Windows path), the ramst.mod is in subfolder
     'c:\\dynare\\4.x.y\\examples'. If the command ``pwd`` on Matlab/Octave command
     prompt returns a different path, you need to change the current path. This
     can be done using the ``cd`` command::

       >> cd c:\dynare\4.x.y\examples

The output of Dynare is left into three main variables in the MATLAB/Octave workspace.

 * ``M_``, a structure containing various information about the model.

 * ``options_``, a structure contains the values of the various options used by Dynare during the computation.

 * ``oo_``, a structure containing the various results of the computations.


^^^^^^^^^^^^
Dynare hooks
^^^^^^^^^^^^

It is possible to call pre and post Dynare preprocessor hooks written as MATLAB
scripts. The script 'MODFILENAME/hooks/priorprocessing.m' is executed before the
call to Dynare’s preprocessor, and can be used to programmatically transform
the mod file that will be read by the preprocessor. The script
'MODFILENAME/hooks/postprocessing.m' is executed just after the call to Dynare’s
preprocessor, and can be used to programmatically transform the files generated
by Dynare’s preprocessor before actual computations start. The pre and/or post
dynare preprocessor hooks are executed if and only if the aforementioned
scripts are detected in the same folder as the the model file, 'FILENAME.mod'.


^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
Understanding Preprocessor Error Messages
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

If the preprocessor runs into an error while processing your '.mod' file, it will
issue an error. Due to the way that a parser works, sometimes these errors can
be misleading. Here, we aim to demystify these error messages.

The preprocessor issues error messages of the form:

 1. ERROR: <<file.mod>>: line A, col B: <<error message>>
 2. ERROR: <<file.mod>>: line A, cols B-C: <<error message>>
 3. ERROR: <<file.mod>>: line A, col B - line C, col D: <<error message>>

The first two errors occur on a single line, with error two spanning multiple
columns. Error three spans multiple rows.

Often, the line and column numbers are precise, leading you directly to the
offending syntax. Infrequently however, because of the way the parser works,
this is not the case. The most common example of misleading line and column
numbers (and error message for that matter) is the case of a missing semicolon,
as seen in the following example::

   varexo a, b
   parameters c, ...;

In this case, the parser doesn’t know a semicolon is missing at the end of the
varexo command until it begins parsing the second line and bumps into the
parameters command. This is because we allow commands to span multiple lines
and, hence, the parser cannot know that the second line will not have a
semicolon on it until it gets there. Once the parser begins parsing the second
line, it realizes that it has encountered a keyword, parameters, which it did
not expect. Hence, it throws an error of the form: ERROR: <<file.mod>>: line 2,
cols 0-9: syntax error, unexpected PARAMETERS. In this case, you would simply
place a semicolon at the end of line one and the parser would continue
processing.
