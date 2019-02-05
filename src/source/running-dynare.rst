.. default-domain:: dynare

##############
Running Dynare
##############

In order to give instructions to Dynare, the user has to write a
*model file* whose filename extension must be ``.mod`` or
``.dyn``. This file contains the description of the model and the
computing tasks required by the user. Its contents are described in
:ref:`model-file`.

.. _dyn-invoc:

Dynare invocation
=================

Once the model file is written, Dynare is invoked using the ``dynare``
command at the MATLAB or Octave prompt (with the filename of the
``.mod`` given as argument).

In practice, the handling of the model file is done in two steps: in
the first one, the model and the processing instructions written by
the user in a *model file* are interpreted and the proper MATLAB or
GNU Octave instructions are generated; in the second step, the program
actually runs the computations. Both steps are triggered automatically
by the ``dynare`` command.

    .. matcomm:: dynare FILENAME[.mod] [OPTIONS…]

    This command launches Dynare and executes the instructions
    included in ``FILENAME.mod``. This user-supplied file contains the
    model and the processing instructions, as described in
    :ref:`model-file`. The options, listed below, can be passed on the
    command line, following the name of the ``.mod`` file or in the
    first line of the ``.mod`` file itself (see below).

    dynare begins by launching the preprocessor on the ``.mod
    file``. By default (unless ``use_dll`` option has been given to
    ``model``), the preprocessor creates three intermediary files:

    - ``filename.m``

        Contains variable declarations, and computing tasks.

    - ``FILENAME_dynamic.m``

        Contains the dynamic model equations. Note that Dynare might
        introduce auxiliary equations and variables (see
        :ref:`aux-variables`). Outputs are the residuals of the
        dynamic model equations in the order the equations were
        declared and the Jacobian of the dynamic model equations. For
        higher order approximations also the Hessian and the
        third-order derivatives are provided. When computing the
        Jacobian of the dynamic model, the order of the endogenous
        variables in the columns is stored in
        ``M_.lead_lag_incidence``. The rows of this matrix represent
        time periods: the first row denotes a lagged (time t-1)
        variable, the second row a contemporaneous (time t) variable,
        and the third row a leaded (time t+1) variable. The columns of
        the matrix represent the endogenous variables in their order
        of declaration. A zero in the matrix means that this
        endogenous does not appear in the model in this time
        period. The value in the ``M_.lead_lag_incidence`` matrix
        corresponds to the column of that variable in the Jacobian of
        the dynamic model. Example: Let the second declared variable
        be ``c`` and the ``(3,2)`` entry of ``M_.lead_lag_incidence``
        be 15. Then the 15th column of the Jacobian is the derivative
        with respect to ``c(+1)``.

    - ``FILENAME_static.m``

        Contains the long run static model equations. Note that Dynare
        might introduce auxiliary equations and variables (see
        :ref:`aux-variables`). Outputs are the residuals of the static
        model equations in the order the equations were declared and
        the Jacobian of the static equations. Entry ``(i,j)`` of the
        Jacobian represents the derivative of the ith static model
        equation with respect to the jth model variable in declaration
        order.


    These files may be looked at to understand errors reported at the
    simulation stage.

    ``dynare`` will then run the computing tasks by executing ``FILENAME.m``.

    A few words of warning are warranted here: the filename of the
    ``.mod`` file should be chosen in such a way that the generated
    ``.m`` files described above do not conflict with ``.m`` files
    provided by MATLAB/Octave or by Dynare. Not respecting this rule
    could cause crashes or unexpected behaviour. In particular, it
    means that the ``.mod`` file cannot be given the name of a
    MATLAB/Octave or Dynare command. Under Octave, it also means that
    the ``.mod`` file cannot be named ``test.mod``.

    *Options*

    .. option:: noclearall

        By default, ``dynare`` will issue a ``clear all`` command to
        MATLAB (<R2015b) or Octave, thereby deleting all workspace
        variables and functions; this option instructs ``dynare`` not
        to clear the workspace. Note that starting with Matlab 2015b
        ``dynare`` only deletes the global variables and the functions
        using persistent variables, in order to benefit from the JIT
        (Just In Time) compilation. In this case the option instructs
        ``dynare`` not to clear the globals and functions.

    .. option:: onlyclearglobals

        By default, ``dynare`` will issue a ``clear all`` command to
        MATLAB versions before 2015b and to Octave, thereby deleting
        all workspace variables; this option instructs ``dynare`` to
        clear only the global variables (i.e. ``M_, options_, oo_,
        estim_params_, bayestopt_``, and ``dataset_``), leaving the
        other variables in the workspace.

    .. option:: debug

        Instructs the preprocessor to write some debugging information
        about the scanning and parsing of the ``.mod`` file.

    .. option:: notmpterms

        Instructs the preprocessor to omit temporary terms in the
        static and dynamic files; this generally decreases
        performance, but is used for debugging purposes since it makes
        the static and dynamic files more readable.

    .. option:: savemacro[=FILENAME]

        Instructs ``dynare`` to save the intermediary file which is
        obtained after macro-processing (see :ref:`macro-proc-lang`);
        the saved output will go in the file specified, or if no file
        is specified in ``FILENAME-macroexp.mod``

    .. option:: onlymacro

        Instructs the preprocessor to only perform the
        macro-processing step, and stop just after. Mainly useful for
        debugging purposes or for using the macro-processor
        independently of the rest of Dynare toolbox.

    .. option:: nolinemacro

        Instructs the macro-preprocessor to omit line numbering
        information in the intermediary ``.mod`` file created after
        the macro-processing step. Useful in conjunction with
        ``savemacro`` when one wants that to reuse the intermediary
        ``.mod`` file, without having it cluttered by line numbering
        directives.

    .. option:: nolog

        Instructs Dynare to no create a logfile of this run in
        ``FILENAME.log.`` The default is to create the logfile.

    .. option:: params_derivs_order=0|1|2

        When :comm:`identification`, :comm:`dynare_sensitivity` (with
        identification), or :ref:`estimation_cmd <estim-comm>` are
        present, this option is used to limit the order of the
        derivatives with respect to the parameters that are calculated
        by the preprocessor. 0 means no derivatives, 1 means first
        derivatives, and 2 means second derivatives. Default: 2

    .. option:: nowarn

        Suppresses all warnings.

    .. option:: json = parse|transform|compute

        Causes the preprocessor to output a version of the ``.mod``
        file in JSON format.

        If ``parse`` is passed, the output will be written after the
        parsing of the ``.mod`` file to a file called
        ``FILENAME.json``.

        If ``transform`` is passed, the JSON output of the transformed
        model (maximum lead of 1, minimum lag of -1, expectation
        operators substituted, etc.) will be written to a file called
        ``FILENAME.json`` and the original, untransformed model will
        be written in ``FILENAME_original.json``.

        And if ``compute`` is passed, the output is written after the
        computing pass. In this case, the transformed model is written
        to ``FILENAME.json``, the original model is written to
        ``FILENAME_original.json``, and the dynamic and static files
        are written to ``FILENAME_dynamic.json`` and
        ``FILENAME_static.json``.

    .. option:: jsonstdout

        Instead of writing output requested by ``json`` to files,
        write to standard out.

    .. option:: onlyjson

        Quit processing once the output requested by ``json`` has been
        written.

    .. option:: jsonderivsimple

        Print a simplified version (excluding variable name(s) and lag
        information) of the static and dynamic files in
        ``FILENAME_static.json`` and ``FILENAME_dynamic.``.

    .. option:: warn_uninit

        Display a warning for each variable or parameter which is not
        initialized. See :ref:`param-init`, or
        :comm:`load_params_and_steady_state
        <load_params_and_steady_state>` for initialization of
        parameters. See :ref:`init-term-cond`, or
        :comm:`load_params_and_steady_state
        <load_params_and_steady_state>` for initialization of
        endogenous and exogenous variables.

    .. option:: console

        Activate console mode. In addition to the behavior of
        ``nodisplay``, Dynare will not use graphical waitbars for long
        computations.

    .. option:: nograph

        Activate the ``nograph`` option (see :opt:`nograph`), so that
        Dynare will not produce any graph.

    .. option:: nointeractive

        Instructs Dynare to not request user input.

    .. option:: nopathchange

        By default Dynare will change Matlab/Octave’s path if
        ``dynare/matlab`` directory is not on top and if Dynare’s
        routines are overriden by routines provided in other
        toolboxes. If one wishes to override Dynare’s routines, the
        ``nopathchange`` options can be used. Alternatively, the path
        can be temporarly modified by the user at the top of the
        ``.mod`` file (using Matlab/Octave’s ``addpath`` command).

    .. option:: nopreprocessoroutput

        Prevent Dynare from printing the output of the steps leading up to the
        preprocessor as well as the preprocessor output itself.

    .. option:: mingw

        Tells Dynare that your MATLAB is configured for compiling MEX
        files with the MinGW compiler from TDM-GCC (see
        :ref:`compil-install`). This option is only available under
        Windows, and is used in conjunction with ``use_dll``.

    .. option:: msvc

        Tells Dynare that your MATLAB is configured for compiling MEX
        files with Microsoft Visual C++ (see
        :ref:`compil-install`). This option is only available under
        Windows, and is used in conjunction with ``use_dll``.

    .. option:: cygwin

        Tells Dynare that your MATLAB is configured for compiling MEX
        files with Cygwin (see :ref:`compil-install`). This option is
        only available under Windows, and is used in conjunction with
        ``use_dll``.

    .. option:: parallel[=CLUSTER_NAME]

        Tells Dynare to perform computations in parallel. If
        CLUSTER_NAME is passed, Dynare will use the specified cluster
        to perform parallel computations. Otherwise, Dynare will use
        the first cluster specified in the configuration file. See
        :ref:`conf-file`, for more information about the configuration
        file.

    .. option:: conffile=FILENAME

        Specifies the location of the configuration file if it differs
        from the default. See :ref:`conf-file`, for more information
        about the configuration file and its default location.

    .. option:: parallel_slave_open_mode

        Instructs Dynare to leave the connection to the slave node
        open after computation is complete, closing this connection
        only when Dynare finishes processing.

    .. option:: parallel_test

        Tests the parallel setup specified in the configuration file
        without executing the ``.mod`` file. See :ref:`conf-file`, for
        more information about the configuration file.

    .. option:: -DMACRO_VARIABLE=MACRO_EXPRESSION

        Defines a macro-variable from the command line (the same
        effect as using the Macro directive ``@#define`` in a model
        file, see :ref:`macro-proc-lang`).

    .. option:: -I<<path>>

        Defines a path to search for files to be included by the
        macroprocessor (using the ``@#include`` command). Multiple
        ``-I`` flags can be passed on the command line. The paths will
        be searched in the order that the ``-I`` flags are passed and
        the first matching file will be used. The flags passed here
        take priority over those passed to ``@#includepath``.

    .. option:: nostrict

        Allows Dynare to issue a warning and continue processing when

        1. there are more endogenous variables than equations.
        2. an undeclared symbol is assigned in ``initval`` or ``endval``.
        3. exogenous variables were declared but not used in the
           ``model`` block.

    .. option:: fast

        Only useful with model option ``use_dll``. Don’t recompile the
        MEX files when running again the same model file and the lists
        of variables and the equations haven’t changed. We use a 32
        bit checksum, stored in ``<model filename>/checksum``. There
        is a very small probability that the preprocessor misses a
        change in the model. In case of doubt, re-run without the fast
        option.

    .. option:: minimal_workspace

        Instructs Dynare not to write parameter assignments to
        parameter names in the .m file produced by the
        preprocessor. This is potentially useful when running
        ``dynare`` on a large ``.mod`` file that runs into workspace
        size limitations imposed by MATLAB.

    .. option:: compute_xrefs

        Tells Dynare to compute the equation cross references, writing
        them to the output ``.m`` file.

    These options can be passed to the preprocessor by listing them
    after the name of the ``.mod`` file. They can alternatively be
    defined in the first line of the ``.mod`` file, this avoids typing
    them on the command line each time a ``.mod`` file is to be
    run. This line must be a Dynare comment (ie must begin with //)
    and the options must be comma separated between ``--+`` options:
    and ``+--``. As in the command line, if an option admits a value
    the equal symbol must not be surrounded by spaces. For instance
    ``json = compute`` is not correct, and should be written
    ``json=compute``.

    *Output*

    Depending on the computing tasks requested in the ``.mod`` file,
    executing the ``dynare`` command will leave variables containing
    results in the workspace available for further processing. More
    details are given under the relevant computing tasks. The
    ``M_``,``oo_``, and ``options_`` structures are saved in a file
    called ``FILENAME_results.mat``. If they exist, ``estim_params_``,
    ``bayestopt_``, ``dataset_``, ``oo_recursive_`` and
    ``estimation_info`` are saved in the same file.

    .. matvar:: M_

        Structure containing various information about the model.

    .. matvar:: options_

        Structure contains the values of the various options used by
        Dynare during the computation.

    .. matvar:: oo_

        Structure containing the various results of the computations.

    .. matvar:: dataset_

        A ``dseries`` object containing the data used for estimation.

    .. matvar:: oo_recursive_

        Cell array containing the ``oo_`` structures obtained when
        estimating the model for the different samples when performing
        recursive estimation and forecasting. The ``oo_`` structure
        obtained for the sample ranging to the `i` -th observation is
        saved in the `i` -th field. The fields for non-estimated
        endpoints are empty.

    *Example*

    Call dynare from the MATLAB or Octave prompt, without or with options:

            .. code-block:: matlab

               >> dynare ramst
               >> dynare ramst.mod savemacro

    Alternatively the options can be passed in the first line of
    ``ramst.mod``:

            .. code-block:: dynare

               // --+ options: savemacro, json=compute +--

    and then dynare called without passing options on the command line:

            .. code-block:: matlab

               >> dynare ramst




Dynare hooks
============

It is possible to call pre and post Dynare preprocessor hooks written
as MATLAB scripts. The script ``MODFILENAME/hooks/priorprocessing.m``
is executed before the call to Dynare’s preprocessor, and can be used
to programmatically transform the mod file that will be read by the
preprocessor. The script ``MODFILENAME/hooks/postprocessing.m`` is
gexecuted just after the call to Dynare’s preprocessor, and can be used
to programmatically transform the files generated by Dynare’s
preprocessor before actual computations start. The pre and/or post
dynare preprocessor hooks are executed if and only if the
aforementioned scripts are detected in the same folder as the the
model file, ``FILENAME.mod``.


Understanding Preprocessor Error Messages
=========================================

If the preprocessor runs into an error while processing your ``.mod``
file, it will issue an error. Due to the way that a parser works,
sometimes these errors can be misleading. Here, we aim to demystify
these error messages.

The preprocessor issues error messages of the form:

   #. ``ERROR: <<file.mod>>: line A, col B: <<error message>>``
   #. ``ERROR: <<file.mod>>: line A, cols B-C: <<error message>>``
   #. ``ERROR: <<file.mod>>: line A, col B - line C, col D: <<error message>>``

The first two errors occur on a single line, with error two spanning
multiple columns. Error three spans multiple rows.

Often, the line and column numbers are precise, leading you directly
to the offending syntax. Infrequently however, because of the way the
parser works, this is not the case. The most common example of
misleading line and column numbers (and error message for that matter)
is the case of a missing semicolon, as seen in the following example::

    varexo a, b
    parameters c, ...;

In this case, the parser doesn’t know a semicolon is missing at the
end of the ``varexo`` command until it begins parsing the second line
and bumps into the ``parameters`` command. This is because we allow
commands to span multiple lines and, hence, the parser cannot know
that the second line will not have a semicolon on it until it gets
there. Once the parser begins parsing the second line, it realizes
that it has encountered a keyword, ``parameters``, which it did not
expect. Hence, it throws an error of the form: ``ERROR: <<file.mod>>:
line 2, cols 0-9: syntax error, unexpected PARAMETERS``. In this case,
you would simply place a semicolon at the end of line one and the
parser would continue processing.
