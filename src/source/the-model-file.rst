.. default-domain:: dynare

.. |br| raw:: html

    <br>

.. _model-file:

##############
The model file
##############

.. _conv:

Conventions
===========

A model file contains a list of commands and of blocks. Each command
and each element of a block is terminated by a semicolon (;). Blocks
are terminated by ``end;``.

Most Dynare commands have arguments and several accept options,
indicated in parentheses after the command keyword. Several options
are separated by commas.

In the description of Dynare commands, the following conventions are
observed:

* Optional arguments or options are indicated between square brackets:
  ‘[]’;
* Repeated arguments are indicated by ellipses: “...”;
* Mutually exclusive arguments are separated by vertical bars: ‘|’;
* INTEGER indicates an integer number;
* INTEGER_VECTOR indicates a vector of integer numbers separated by
  spaces, enclosed by square brackets;
* DOUBLE indicates a double precision number. The following syntaxes
  are valid: ``1.1e3``, ``1.1E3``, ``1.1d3``, ``1.1D3``. In some
  places, infinite Values ``Inf`` and ``-Inf`` are also allowed;
* NUMERICAL_VECTOR indicates a vector of numbers separated by spaces,
  enclosed by square brackets;
* EXPRESSION indicates a mathematical expression valid outside the
  model description (see :ref:`expr`);
* MODEL_EXPRESSION (sometimes MODEL_EXP) indicates a mathematical
  expression valid in the model description (see :ref:`expr` and
  :ref:`model-decl`);
* MACRO_EXPRESSION designates an expression of the macro-processor
  (see :ref:`macro-exp`);
* VARIABLE_NAME (sometimes VAR_NAME) indicates a variable name
  starting with an alphabetical character and can’t contain:
  ‘()+-\*/^=!;:@#.’ or accentuated characters;
* PARAMETER_NAME (sometimes PARAM_NAME) indicates a parameter name
  starting with an alphabetical character and can’t contain:
  ‘()+-\*/^=!;:@#.’ or accentuated characters;
* LATEX_NAME (sometimes TEX_NAME) indicates a valid
  LaTeX expression in math mode (not including the
  dollar signs);
* FUNCTION_NAME indicates a valid MATLAB function name;
* FILENAME indicates a filename valid in the underlying operating
  system; it is necessary to put it between quotes when specifying the
  extension or if the filename contains a non-alphanumeric character;


.. _var-decl:

Variable declarations
=====================

While Dynare allows the user to choose their own variable names, there
are some restrictions to be kept in mind. First, variables and
parameters must not have the same name as Dynare commands or built-in
functions. In this respect, Dynare is not case-sensitive. For example,
do not use ``Ln`` or ``Sigma_e`` to name your variable. Not conforming
to this rule might yield hard-to-debug error messages or
crashes. Second, to minimize interference with MATLAB or Octave
functions that may be called by Dynare or user-defined steady state
files, it is recommended to avoid using the name of MATLAB
functions. In particular when working with steady state files, do not
use correctly-spelled greek names like `alpha`, because there are
Matlab functions of the same name. Rather go for ``alppha`` or
``alph``. Lastly, please do not name a variable or parameter
``i``. This may interfere with the imaginary number i and the index in
many loops. Rather, name investment ``invest``. Using ``inv`` is also
not recommended as it already denotes the inverse operator. Commands
for declaring variables and parameters are described below.

.. command:: var VAR_NAME [$TEX_NAME$] [(long_name=QUOTED_STR|NAME=QUOTED_STR)]...;
    var(deflator=MODEL_EXPR) VAR_NAME (... same options apply)
    var(log_deflator=MODEL_EXPR) VAR_NAME (... same options apply)

    |br| This required command declares the endogenous variables in
    the model. See :ref:`conv` for the syntax of *VAR_NAME* and
    *MODEL_EXPR*. Optionally it is possible to give a
    LaTeX name to the variable or, if it is
    nonstationary, provide information regarding its deflator. ``var``
    commands can appear several times in the file and Dynare will
    concatenate them. Dynare stores the list of declared parameters,
    in the order of declaration, in a column cell array
    ``M_.endo_names``.

    *Options*

    If the model is nonstationary and is to be written as such in the
    ``model`` block, Dynare will need the trend deflator for the
    appropriate endogenous variables in order to stationarize the
    model. The trend deflator must be provided alongside the variables
    that follow this trend.


    .. option:: deflator = MODEL_EXPR

        The expression used to detrend an endogenous variable. All
        trend variables, endogenous variables and parameters
        referenced in MODEL_EXPR must already have been declared by
        the ``trend_var, log_trend_var, var`` and ``parameters``
        commands. The deflator is assumed to be multiplicative; for an
        additive deflator, use ``log_deflator``.

    .. option:: log_deflator = MODEL_EXPR

        Same as ``deflator``, except that the deflator is assumed to
        be additive instead of multiplicative (or, to put it
        otherwise, the declared variable is equal to the log of a
        variable with a multiplicative trend).

    .. _long-name:

    .. option:: long_name = QUOTED_STR

        This is the long version of the variable name. Its value is
        stored in ``M_.endo_names_long`` (a column cell array, in the
        same order as ``M_.endo_names``). In case multiple
        ``long_name`` options are provided, the last one will be
        used. Default: ``VAR_NAME``.

    .. _partitioning:

    .. option:: NAME = QUOTED_STR

        This is used to create a partitioning of variables. It results
        in the direct output in the ``.m`` file analogous to:
        ``M_.endo_partitions.NAME = QUOTED_STR``;.

    *Example (variable partitioning)*

      ::

        var c gnp cva (country=`US', state=`VA')
                  cca (country=`US', state=`CA', long_name=`Consumption CA');
        var(deflator=A) i b;
        var c $C$ (long_name=`Consumption');


.. command :: varexo VAR_NAME [$TEX_NAME$] [(long_name=QUOTED_STR|NAME=QUOTED_STR)...];

    |br| This optional command declares the exogenous variables in the
    model. See :ref:`conv` for the syntax of ``VAR_NAME``. Optionally
    it is possible to give a LaTeX name to the
    variable. Exogenous variables are required if the user wants to be
    able to apply shocks to her model. ``varexo`` commands can appear
    several times in the file and Dynare will concatenate them.

    *Options*

    .. option:: long_name = QUOTED_STRING

        Like :ref:`long_name <long-name>` but value stored in ``M_.exo_names_long``.

    .. option:: NAME = QUOTED_STRING

        Like :ref:`partitioning <partitioning>` but QUOTED_STRING
        stored in ``M_.exo_partitions.NAME``.

    *Example*

        ::

            varexo m gov;


    *Remarks*

    An exogenous variable is an innovation, in the sense
    that this variable cannot be predicted from the knowledge of the
    current state of the economy. For instance, if logged TFP is a first
    order autoregressive process:

    .. math::

       a_t  = \rho a_{t-1} + \varepsilon_t

    then logged TFP :math:`a_t` is an endogenous variable to be
    declared with ``var``, its best prediction is :math:`\rho
    a_{t-1}`, while the innovation :math:`\varepsilon_t` is to be
    declared with ``varexo``.


.. command:: varexo_det VAR_NAME [$TEX_NAME$] [(long_name=QUOTED_STR|NAME=QUOTED_STR)...];

    |br| This optional command declares exogenous deterministic
    variables in a stochastic model. See :ref:`conv` for the syntax of
    VARIABLE_NAME. Optionally it is possible to give a LaTeX
    name to the variable. ``varexo_det`` commands can appear several
    times in the file and Dynare will concatenate them.

    It is possible to mix deterministic and stochastic shocks to build
    models where agents know from the start of the simulation about
    future exogenous changes. In that case ``stoch_simul`` will
    compute the rational expectation solution adding future
    information to the state space (nothing is shown in the output of
    ``stoch_simul``) and forecast will compute a simulation
    conditional on initial conditions and future information.


    *Options*

    .. option:: long_name = QUOTED_STRING

       Like :ref:`long_name <long-name>` but value stored in
       ``M_.exo_det_names_long``.

    .. option:: NAME = QUOTED_STRING

       Like :ref:`partitioning <partitioning>` but QUOTED_STRING stored
       in ``M_.exo_det_partitions.NAME``.

    *Example*

        ::

            varexo m gov;
            varexo_det tau;


.. command :: parameters PARAM_NAME [$TEX_NAME$] [(long_name=QUOTED_STR|NAME=QUOTED_STR)...];

    |br| This command declares parameters used in the model, in variable
    initialization or in shocks declarations. See :ref:`conv` for the
    syntax of ``PARAM_NAME``. Optionally it is possible to give a
    LaTeX name to the parameter.

    The parameters must subsequently be assigned values (see :ref:`param-init`).

    ``parameters`` commands can appear several times in the file and Dynare will concatenate them.

    *Options*

    .. option:: long_name = QUOTED_STRING

        Like :ref:`long_name <long-name>` but value stored in ``M_.param_names_long``.

    .. option:: NAME = QUOTED_STRING

        Like :ref:`partitioning <partitioning>` but QUOTED_STRING stored in ``M_.param_partitions.NAME``.

    *Example*

        ::

            parameters alpha, bet;


.. command :: change_type (var|varexo|varexo_det|parameters) VAR_NAME | PARAM_NAME...;

    Changes the types of the specified variables/parameters to another
    type: endogenous, exogenous, exogenous deterministic or
    parameter. It is important to understand that this command has a
    global effect on the ``.mod`` file: the type change is effective
    after, but also before, the ``change_type`` command. This command
    is typically used when flipping some variables for steady state
    calibration: typically a separate model file is used for
    calibration, which includes the list of variable declarations with
    the macro-processor, and flips some variable.

    *Example*

        ::

            var y, w;
            parameters alpha, beta;
            ...
            change_type(var) alpha, beta;
            change_type(parameters) y, w;

        Here, in the whole model file, ``alpha`` and ``beta`` will be
        endogenous and ``y`` and ``w`` will be parameters.


.. command:: predetermined_variables VAR_NAME...;

    |br| In Dynare, the default convention is that the timing of a variable
    reflects when this variable is decided. The typical example is for
    capital stock: since the capital stock used at current period is
    actually decided at the previous period, then the capital stock
    entering the production function is ``k(-1)``, and the law of
    motion of capital must be written::

        k = i + (1-delta)*k(-1)

    Put another way, for stock variables, the default in Dynare is to
    use a “stock at the end of the period” concept, instead of a
    “stock at the beginning of the period” convention.

    The ``predetermined_variables`` is used to change that
    convention. The endogenous variables declared as predetermined
    variables are supposed to be decided one period ahead of all other
    endogenous variables. For stock variables, they are supposed to
    follow a “stock at the beginning of the period” convention.

    Note that Dynare internally always uses the “stock at the end of
    the period” concept, even when the model has been entered using
    the ``predetermined_variables`` command. Thus, when plotting,
    computing or simulating variables, Dynare will follow the
    convention to use variables that are decided in the current
    period. For example, when generating impulse response functions
    for capital, Dynare will plot ``k``, which is the capital stock
    decided upon by investment today (and which will be used in
    tomorrow’s production function). This is the reason that capital
    is shown to be moving on impact, because it is ``k`` and not the
    predetermined ``k(-1)`` that is displayed. It is important to
    remember that this also affects simulated time series and output
    from smoother routines for predetermined variables. Compared to
    non-predetermined variables they might otherwise appear to be
    falsely shifted to the future by one period.

    *Example*

        The following two program snippets are strictly equivalent.

        Using default Dynare timing convention::

            var y, k, i;
            ...
            model;
            y = k(-1)^alpha;
            k = i + (1-delta)*k(-1);
            ...
            end;

        Using the alternative timing convention::

            var y, k, i;
            predetermined_variables k;
            ...
            model;
            y = k^alpha;
            k(+1) = i + (1-delta)*k;
            ...
            end;


.. command:: trend_var (growth_factor = MODEL_EXPR) VAR_NAME [$LATEX_NAME$]...;

    |br| This optional command declares the trend variables in the
    model. See ref:`conv` for the syntax of MODEL_EXPR and
    VAR_NAME. Optionally it is possible to give a
    LaTeX name to the variable.

    The variable is assumed to have a multiplicative growth trend. For
    an additive growth trend, use ``log_trend_var`` instead.

    Trend variables are required if the user wants to be able to write
    a nonstationary model in the ``model`` block. The ``trend_var``
    command must appear before the var command that references the
    trend variable.

    ``trend_var`` commands can appear several times in the file and
    Dynare will concatenate them.

    If the model is nonstationary and is to be written as such in the
    ``model`` block, Dynare will need the growth factor of every trend
    variable in order to stationarize the model. The growth factor
    must be provided within the declaration of the trend variable,
    using the ``growth_factor`` keyword. All endogenous variables and
    parameters referenced in MODEL_EXPR must already have been
    declared by the var and parameters commands.

    *Example*

        ::

            trend_var (growth_factor=gA) A;


.. command :: log_trend_var (log_growth_factor = MODEL_EXPR) VAR_NAME [$LATEX_NAME$]...;

    |br| Same as ``trend_var``, except that the variable is supposed to
    have an additive trend (or, to put it otherwise, to be equal to
    the log of a variable with a multiplicative trend).


.. command::  model_local_variable VARIABLE_NAME [LATEX_NAME]... ;

    |br| This optional command declares a model local variable. See
    :ref:`conv` for the syntax of VARIABLE_NAME. As you can create
    model local variables on the fly in the model block (see
    :ref:`model-decl`), the interest of this command is primarily to
    assign a LATEX_NAME to the model local variable.

    *Example*

        ::

            model_local_variable GDP_US $GDPUS$;


.. _expr:

Expressions
===========

Dynare distinguishes between two types of mathematical expressions:
those that are used to describe the model, and those that are used
outside the model block (e.g. for initializing parameters or
variables, or as command options). In this manual, those two types of
expressions are respectively denoted by MODEL_EXPRESSION and
EXPRESSION.

Unlike MATLAB or Octave expressions, Dynare expressions are
necessarily scalar ones: they cannot contain matrices or evaluate to
matrices [#f1]_.

Expressions can be constructed using integers (INTEGER), floating
point numbers (DOUBLE), parameter names (PARAMETER_NAME), variable
names (VARIABLE_NAME), operators and functions.

The following special constants are also accepted in some contexts:

.. constant:: inf

    Represents infinity.

.. constant:: nan

    “Not a number”: represents an undefined or unrepresentable value.


Parameters and variables
------------------------

Parameters and variables can be introduced in expressions by simply
typing their names. The semantics of parameters and variables is quite
different whether they are used inside or outside the model block.


Inside the model
^^^^^^^^^^^^^^^^

Parameters used inside the model refer to the value given through
parameter initialization (see :ref:`param-init`) or ``homotopy_setup``
when doing a simulation, or are the estimated variables when doing an
estimation.

Variables used in a MODEL_EXPRESSION denote current period values when
neither a lead or a lag is given. A lead or a lag can be given by
enclosing an integer between parenthesis just after the variable name:
a positive integer means a lead, a negative one means a lag. Leads or
lags of more than one period are allowed. For example, if ``c`` is an
endogenous variable, then ``c(+1)`` is the variable one period ahead,
and ``c(-2)`` is the variable two periods before.

When specifying the leads and lags of endogenous variables, it is
important to respect the following convention: in Dynare, the timing
of a variable reflects when that variable is decided. A control
variable — which by definition is decided in the current period — must
have no lead. A predetermined variable — which by definition has been
decided in a previous period — must have a lag. A consequence of this
is that all stock variables must use the “stock at the end of the
period” convention.

Leads and lags are primarily used for endogenous variables, but can be
used for exogenous variables. They have no effect on parameters and
are forbidden for local model variables (see Model declaration).


Outside the model
^^^^^^^^^^^^^^^^^

When used in an expression outside the model block, a parameter or a
variable simply refers to the last value given to that variable. More
precisely, for a parameter it refers to the value given in the
corresponding parameter initialization (see :ref:`param-init`); for an
endogenous or exogenous variable, it refers to the value given in the
most recent ``initval`` or ``endval`` block.


Operators
---------

The following operators are allowed in both MODEL_EXPRESSION and
EXPRESSION:

* Binary arithmetic operators: ``+``, ``-``, ``*``, ``/``, ``^``
* Unary arithmetic operators: ``+``, ``-``
* Binary comparison operators (which evaluate to either 0 or 1): ``<``,
  ``>``, ``<=``, ``>=``, ``==``, ``!=``

Note the binary comparison operators are differentiable everywhere except on a
line of the 2-dimensional real plane. However for facilitating
convergence of Newton-type methods, Dynare assumes that, at the points
of non-differentiability, the partial derivatives of these operators
with respect to both arguments is equal to 0 (since this is the value
of the partial derivatives everywhere else).

The following special operators are accepted in MODEL_EXPRESSION (but
not in EXPRESSION):

.. operator:: STEADY_STATE (MODEL_EXPRESSION)

    This operator is used to take the value of the enclosed expression
    at the steady state. A typical usage is in the Taylor rule, where
    you may want to use the value of GDP at steady state to compute
    the output gap.

.. operator:: EXPECTATION (INTEGER) (MODEL_EXPRESSION)

    This operator is used to take the expectation of some expression
    using a different information set than the information available
    at current period. For example, ``EXPECTATION(-1)(x(+1))`` is
    equal to the expected value of variable x at next period, using
    the information set available at the previous period. See
    :ref:`aux-variables` for an explanation of how this operator is
    handled internally and how this affects the output.


Functions
---------

Built-in functions
^^^^^^^^^^^^^^^^^^

The following standard functions are supported internally for both
MODEL_EXPRESSION and EXPRESSION:

.. function:: exp(x)

    Natural exponential.

.. function:: log(x)
.. function:: ln(x)

    Natural logarithm.

.. function:: log10(x)

    Base 10 logarithm.

.. function:: sqrt(x)

    Square root.

.. function:: sign(x)

    Signum function, defined as:

        .. math::

           \textrm{sign}(x) =
                  \begin{cases}
                  -1 &\quad\text{if }x<0\\
                  0 &\quad\text{if }x=0\\
                  1 &\quad\text{if }x>0
                  \end{cases}


    Note that this function is not continuous, hence not  differentiable, at
    :math:`x=0`. However, for facilitating convergence of Newton-type
    methods, Dynare assumes that the derivative at :math:`x=0` is
    equal to :math:`0`. This assumption comes from the observation
    that both the right- and left-derivatives at this point exist and
    are equal to :math:`0`, so we can remove the singularity by
    postulating that the derivative at :math:`x=0` is :math:`0`.

.. function:: abs(x)

    Absolute value.

    Note that this continuous function is not differentiable at
    :math:`x=0`. However, for facilitating convergence of Newton-type
    methods, Dynare assumes that the derivative at :math:`x=0` is
    equal to :math:`0` (even if the derivative does not exist). The
    rational for this mathematically unfounded definition, rely on the
    observation that the derivative of :math:`\mathrm{abs}(x)` is equal to
    :math:`\mathrm{sign}(x)` for any :math:`x\neq 0` in :math:`\mathbb R` and
    from the convention for the value of :math:`\mathrm{sign}(x)` at
    :math:`x=0`).

.. function:: sin(x)
.. function:: cos(x)
.. function:: tan(x)
.. function:: asin(x)
.. function:: acos(x)
.. function:: atan(x)

    Trigonometric functions.

.. function:: max(a, b)
.. function:: min(a, b)

    Maximum and minimum of two reals.

    Note that these functions are differentiable everywhere except on
    a line of the 2-dimensional real plane defined by
    :math:`a=b`. However for facilitating convergence of Newton-type
    methods, Dynare assumes that, at the points of
    non-differentiability, the partial derivative of these functions
    with respect to the first (resp. the second) argument is equal to
    :math:`1` (resp. to :math:`0`) (i.e. the derivatives at the kink
    are equal to the derivatives observed on the half-plane where the
    function is equal to its first argument).

.. function:: normcdf(x)
              normcdf(x, mu, sigma)

    Gaussian cumulative density function, with mean *mu* and standard
    deviation *sigma*. Note that ``normcdf(x)`` is equivalent to
    ``normcdf(x,0,1)``.

.. function:: normpdf(x)
              normpdf(x, mu, sigma)

    Gaussian probability density function, with mean *mu* and standard
    deviation *sigma*. Note that ``normpdf(x)`` is equivalent to
    ``normpdf(x,0,1)``.

.. function:: erf(x)

    Gauss error function.


External functions
^^^^^^^^^^^^^^^^^^

Any other user-defined (or built-in) MATLAB or Octave function may be
used in both a MODEL_EXPRESSION and an EXPRESSION, provided that this
function has a scalar argument as a return value.

To use an external function in a MODEL_EXPRESSION, one must declare
the function using the ``external_function`` statement. This is not
required for external functions used in an EXPRESSION outside of a
``model`` block or ``steady_state_model`` block.

.. command:: external_function (OPTIONS...);

    This command declares the external functions used in the model
    block. It is required for every unique function used in the model
    block.

    ``external_function`` commands can appear several times in the
    file and must come before the model block.

    *Options*

    .. option:: name = NAME

        The name of the function, which must also be the name of the
        M-/MEX file implementing it. This option is mandatory.

    .. option:: nargs = INTEGER

        The number of arguments of the function. If this option is not
        provided, Dynare assumes ``nargs = 1``.

    .. option:: first_deriv_provided [= NAME]

        If NAME is provided, this tells Dynare that the Jacobian is
        provided as the only output of the M-/MEX file given as the
        option argument. If NAME is not provided, this tells Dynare
        that the M-/MEX file specified by the argument passed to NAME
        returns the Jacobian as its second output argument.

    .. option:: second_deriv_provided [= NAME]

        If NAME is provided, this tells Dynare that the Hessian is
        provided as the only output of the M-/MEX file given as the
        option argument. If NAME is not provided, this tells Dynare
        that the M-/MEX file specified by the argument passed to NAME
        returns the Hessian as its third output argument. NB: This
        option can only be used if the ``first_deriv_provided`` option
        is used in the same ``external_function`` command.

    *Example*

        ::

           external_function(name = funcname);
           external_function(name = otherfuncname, nargs = 2, first_deriv_provided, second_deriv_provided);
           external_function(name = yetotherfuncname, nargs = 3, first_deriv_provided = funcname_deriv);


A few words of warning in stochastic context
--------------------------------------------

The use of the following functions and operators is strongly
discouraged in a stochastic context: ``max``, ``min``, ``abs``,
``sign``, ``<``, ``>``, ``<=``, ``>=``, ``==``, ``!=``.

The reason is that the local approximation used by ``stoch_simul`` or
``estimation`` will by nature ignore the non-linearities introduced by
these functions if the steady state is away from the kink. And, if the
steady state is exactly at the kink, then the approximation will be
bogus because the derivative of these functions at the kink is bogus
(as explained in the respective documentations of these functions and
operators).

Note that ``extended_path`` is not affected by this problem, because
it does not rely on a local approximation of the mode.


.. _param-init:

Parameter initialization
========================

When using Dynare for computing simulations, it is necessary to
calibrate the parameters of the model. This is done through parameter
initialization.

The syntax is the following::

    PARAMETER_NAME = EXPRESSION;

Here is an example of calibration::

    parameters alpha, beta;

    beta = 0.99;
    alpha = 0.36;
    A = 1-alpha*beta;

Internally, the parameter values are stored in ``M_.params``:

.. matvar:: M_.params

    Contains the values of model parameters. The parameters are in the
    order that was used in the parameters command, hence oredered as
    in ``M_.param_names``.


.. _model-decl:

Model declaration
=================

The model is declared inside a ``model`` block:

.. block:: model ;
   model (OPTIONS...);

    |br| The equations of the model are written in a block delimited by
    ``model`` and ``end`` keywords.

    There must be as many equations as there are endogenous variables
    in the model, except when computing the unconstrained optimal
    policy with ``ramsey_model``, ``ramsey_policy`` or
    ``discretionary_policy``.

    The syntax of equations must follow the conventions for
    MODEL_EXPRESSION as described in :ref:`expr`. Each equation
    must be terminated by a semicolon (‘;’). A normal equation looks
    like:

        MODEL_EXPRESSION = MODEL_EXPRESSION;

    |br| When the equations are written in homogenous form, it is possible
    to omit the ‘=0’ part and write only the left hand side of the
    equation. A homogenous equation looks like:

        MODEL_EXPRESSION;

    |br| Inside the model block, Dynare allows the creation of
    *model-local variables*, which constitute a simple way to share a
    common expression between several equations. The syntax consists
    of a pound sign (#) followed by the name of the new model local
    variable (which must **not** be declared as in :ref:`var-decl`,
    but may have been declared by :comm:`model_local_variable`), an
    equal sign, and the expression for which this new variable will
    stand. Later on, every time this variable appears in the model,
    Dynare will substitute it by the expression assigned to the
    variable. Note that the scope of this variable is restricted to
    the model block; it cannot be used outside. To assign a LaTeX name
    to the model local variable, use the declaration syntax outlined
    by :comm:`model_local_variable`. A model local variable declaration
    looks like:

        #VARIABLE_NAME = MODEL_EXPRESSION;

    |br| It is possible to tag equations written in the model block. A tag
    can serve different purposes by allowing the user to attach
    arbitrary informations to each equation and to recover them at
    runtime. For instance, it is possible to name the equations with a
    ``name``-tag, using a syntax like::

        model;

        [name = 'Budget constraint'];
        c + k = k^theta*A;

        end;

    Here, ``name`` is the keyword indicating that the tag names the
    equation. If an equation of the model is tagged with a name, the
    ``resid`` command will display the name of the equations (which
    may be more informative than the equation numbers) in addition to
    the equation number. Several tags for one equation can be
    separated using a comma::

        model;

        [name='Taylor rule',mcp = 'r > -1.94478']
        r = rho*r(-1) + (1-rho)*(gpi*Infl+gy*YGap) + e;

        end;

    More information on tags is available on the `Dynare wiki`_.

    *Options*

    .. option:: linear

        Declares the model as being linear. It spares oneself from
        having to declare initial values for computing the steady
        state of a stationary linear model. This option can’t be used
        with non-linear models, it will NOT trigger linearization of
        the model.

    .. option:: use_dll

        Instructs the preprocessor to create dynamic loadable
        libraries (DLL) containing the model equations and
        derivatives, instead of writing those in M-files. You need a
        working compilation environment, i.e. a working ``mex``
        command (see :ref:`compil-install` for more details). On
        MATLAB for Windows, you will need to also pass the compiler
        name at the command line. Using this option can result in
        faster simulations or estimations, at the expense of some
        initial compilation time. [#f2]_

    .. option:: block

        Perform the block decomposition of the model, and exploit it
        in computations (steady-state, deterministic simulation,
        stochastic simulation with first order approximation and
        estimation). See `Dynare wiki`_ for details on the algorithms
        used in deterministic simulation and steady-state computation.

    .. option:: bytecode

        Instead of M-files, use a bytecode representation of the
        model, i.e. a binary file containing a compact representation
        of all the equations.

    .. option:: cutoff = DOUBLE

        Threshold under which a jacobian element is considered as null
        during the model normalization. Only available with option
        ``block``. Default: ``1e-15``

    .. option:: mfs = INTEGER

        Controls the handling of minimum feedback set of endogenous
        variables. Only available with option ``block``. Possible
        values:

        ``0``

            All the endogenous variables are considered as feedback
            variables (Default).

        ``1``

            The endogenous variables assigned to equation naturally
            normalized (i.e. of the form :math:`x=f(Y)` where
            :math:`x` does not appear in :math:`Y`) are potentially
            recursive variables. All the other variables are forced to
            belong to the set of feedback variables.

        ``2``

            In addition of variables with ``mfs = 1`` the endogenous
            variables related to linear equations which could be
            normalized are potential recursive variables. All the
            other variables are forced to belong to the set of
            feedback variables.

        ``3``

            In addition of variables with ``mfs = 2`` the endogenous
            variables related to non-linear equations which could be
            normalized are potential recursive variables. All the
            other variables are forced to belong to the set of
            feedback variables.

    .. option:: no_static

        Don’t create the static model file. This can be useful for
        models which don’t have a steady state.

    .. option:: differentiate_forward_vars differentiate_forward_vars = ( VARIABLE_NAME [VARIABLE_NAME ...] )

        Tells Dynare to create a new auxiliary variable for each
        endogenous variable that appears with a lead, such that the
        new variable is the time differentiate of the original
        one. More precisely, if the model contains ``x(+1)``, then a
        variable ``AUX_DIFF_VAR`` will be created such that
        ``AUX_DIFF_VAR=x-x(-1)``, and ``x(+1)`` will be replaced with
        ``x+AUX_DIFF_VAR(+1)``.

        The transformation is applied to all endogenous variables with
        a lead if the option is given without a list of variables. If
        there is a list, the transformation is restricted to
        endogenous with a lead that also appear in the list.

        This option can useful for some deterministic simulations
        where convergence is hard to obtain. Bad values for terminal
        conditions in the case of very persistent dynamics or
        permanent shocks can hinder correct solutions or any
        convergence. The new differentiated variables have obvious
        zero terminal conditions (if the terminal condition is a
        steady state) and this in many cases helps convergence of
        simulations.

    .. option:: parallel_local_files = ( FILENAME [, FILENAME]... )

        Declares a list of extra files that should be transferred to
        slave nodes when doing a parallel computation (see
        :ref:`paral-conf`).

    *Example* (Elementary RBC model)

        ::

            var c k;
            varexo x;
            parameters aa alph bet delt gam;

            model;
            c =  - k + aa*x*k(-1)^alph + (1-delt)*k(-1);
            c^(-gam) = (aa*alph*x(+1)*k^(alph-1) + 1 - delt)*c(+1)^(-gam)/(1+bet);
            end;

    *Example* (Use of model local variables)

        The following program::

            model;
            # gamma = 1 - 1/sigma;
            u1 = c1^gamma/gamma;
            u2 = c2^gamma/gamma;
            end;

        ...is formally equivalent to::

            model;
            u1 = c1^(1-1/sigma)/(1-1/sigma);
            u2 = c2^(1-1/sigma)/(1-1/sigma);
            end;

    *Example* (A linear model)

        ::

         model(linear);
         x = a*x(-1)+b*y(+1)+e_x;
         y = d*y(-1)+e_y;
         end;


Dynare has the ability to output the original list of model equations
to a LaTeX file, using the ``write_latex_original_model``
command, the list of transformed model equations using the
``write_latex_dynamic_model command``, and the list of static model
equations using the ``write_latex_static_model`` command.

.. command:: write_latex_original_model (OPTIONS);

    |br| This command creates two LaTeX files: one
    containing the model as defined in the model block and one
    containing the LaTeX document header information.

    If your ``.mod`` file is ``FILENAME.mod``, then Dynare will create
    a file called ``FILENAME_original.tex``, which includes a file
    called ``FILENAME_original_content.tex`` (also created by Dynare)
    containing the list of all the original model equations.

    If LaTeX names were given for variables and parameters
    (see :ref:`var-decl`), then those will be used; otherwise, the
    plain text names will be used.

    Time subscripts (``t``, ``t+1``, ``t-1``, ...) will be appended to
    the variable names, as LaTeX subscripts.

    Compiling the TeX file requires the following LaTeX
    packages: ``geometry, fullpage, breqn``.

    *Options*

    .. option:: write_equation_tags

        Write the equation tags in the LaTeX output. The
        equation tags will be interpreted with LaTeX markups.

.. command:: write_latex_dynamic_model ;
             write_latex_dynamic_model (OPTIONS);

    |br| This command creates two LaTeX files: one containing
    the dynamic model and one containing the LaTeX document
    header information.

    If your ``.mod`` file is ``FILENAME.mod``, then Dynare will create
    a file called ``FILENAME_dynamic.tex``, which includes a file
    called ``FILENAME_dynamic_content.tex`` (also created by Dynare)
    containing the list of all the dynamic model equations.

    If LaTeX names were given for variables and parameters
    (see :ref:`var-decl`), then those will be used; otherwise, the
    plain text names will be used.

    Time subscripts (``t``, ``t+1``, ``t-1``, ...) will be appended to
    the variable names, as LaTeX subscripts.

    Note that the model written in the TeX file will differ from the
    model declared by the user in the following dimensions:

        * The timing convention of predetermined variables (see
          :comm:`predetermined_variables`) will have been changed to
          the default Dynare timing convention; in other words,
          variables declared as predetermined will be lagged on period
          back,
        * The expectation operators (see :op:`expectation <EXPECTATION
          (INTEGER) (MODEL_EXPRESSION)>`) will have been removed,
          replaced by auxiliary variables and new equations as
          explained in the documentation of the operator,
        * Endogenous variables with leads or lags greater or equal
          than two will have been removed, replaced by new auxiliary
          variables and equations,
        * F_or a stochastic model, exogenous variables with leads or
          lags will also have been replaced by new auxiliary variables
          and equations.

    For the required LaTeX packages, see
    :comm:`write_latex_original_model`.

    *Options*

    .. option:: write_equation_tags

        See :opt:`write_equation_tags`


.. command:: write_latex_static_model (OPTIONS);

    |br| This command creates two LaTeX files: one
    containing the static model and one containing the LaTeX
    document header information.

    If your ``.mod`` file is ``FILENAME.mod``, then Dynare will create
    a file called ``FILENAME_static.tex``, which includes a file
    called ``FILENAME_static_content.tex`` (also created by Dynare)
    containing the list of all the steady state model equations.

    If LaTeX names were given for variables and parameters
    (see :ref:`var-decl`), then those will be used; otherwise, the
    plain text names will be used.

    Note that the model written in the TeX file will differ from the
    model declared by the user in the some dimensions (see
    :comm:`write_latex_dynamic_model` for details).

    Also note that this command will not output the contents of the
    optional ``steady_state_model`` block (see
    :bck:`steady_state_model`); it will rather output a static version
    (i.e. without leads and lags) of the dynamic ``model`` declared in
    the model block. To write the LaTeX contents of the
    ``steady_state_model`` see :comm:`write_latex_steady_state_model`.

    For the required LaTeX packages, see
    :comm:`write_latex_original_model`.

    *Options*

    .. option:: write_equation_tags

        See :opt:`write_equation_tags`.

.. command:: write_latex_steady_state_model

    |br| This command creates two LaTeX files: one containing the steady
    state model and one containing the LaTeX document header
    information.

    If your ``.mod`` file is ``FILENAME.mod``, then Dynare
    will create a file called ``FILENAME_steady_state.tex``,
    which includes a file called
    ``FILENAME_steady_state_content.tex`` (also created by
    Dynare) containing the list of all the steady state model
    equations.

    If LaTeX names were given for variables and parameters
    (see :ref:`var-decl`), then those will be used;
    otherwise, the plain text names will be used.

    Note that the model written in the ``.tex`` file will differ from
    the model declared by the user in some dimensions
    (see :comm:`write_latex_dynamic_model` for details).

    For the required LaTeX packages, see :comm:`write_latex_original_model`.


.. _aux-variables:

Auxiliary variables
===================

The model which is solved internally by Dynare is not exactly the
model declared by the user. In some cases, Dynare will introduce
auxiliary endogenous variables—along with corresponding auxiliary
equations—which will appear in the final output.

The main transformation concerns leads and lags. Dynare will perform a
transformation of the model so that there is only one lead and one lag
on endogenous variables and, in the case of a stochastic model, no
leads/lags on exogenous variables.

This transformation is achieved by the creation of auxiliary variables
and corresponding equations. For example, if ``x(+2)`` exists in the
model, Dynare will create one auxiliary variable ``AUX_ENDO_LEAD =
x(+1)``, and replace ``x(+2)`` by ``AUX_ENDO_LEAD(+1)``.

A similar transformation is done for lags greater than 2 on endogenous
(auxiliary variables will have a name beginning with
``AUX_ENDO_LAG``), and for exogenous with leads and lags (auxiliary
variables will have a name beginning with ``AUX_EXO_LEAD`` or
``AUX_EXO_LAG`` respectively).

Another transformation is done for the ``EXPECTATION`` operator. For
each occurrence of this operator, Dynare creates an auxiliary variable
defined by a new equation, and replaces the expectation operator by a
reference to the new auxiliary variable. For example, the expression
``EXPECTATION(-1)(x(+1))`` is replaced by ``AUX_EXPECT_LAG_1(-1)``,
and the new auxiliary variable is declared as ``AUX_EXPECT_LAG_1 =
x(+2)``.

Auxiliary variables are also introduced by the preprocessor for the
``ramsey_model`` and ``ramsey_policy`` commands. In this case, they
are used to represent the Lagrange multipliers when first order
conditions of the Ramsey problem are computed. The new variables take
the form ``MULT_i``, where *i* represents the constraint with which
the multiplier is associated (counted from the order of declaration in
the model block).

The last type of auxiliary variables is introduced by the
``differentiate_forward_vars`` option of the model block. The new
variables take the form ``AUX_DIFF_FWRD_i``, and are equal to
``x-x(-1)`` for some endogenous variable ``x``.

Once created, all auxiliary variables are included in the set of
endogenous variables. The output of decision rules (see below) is such
that auxiliary variable names are replaced by the original variables
they refer to.

The number of endogenous variables before the creation of auxiliary
variables is stored in ``M_.orig_endo_nbr``, and the number of
endogenous variables after the creation of auxiliary variables is
stored in ``M_.endo_nbr``.

See `Dynare wiki`_ for more technical details on auxiliary variables.


.. _init-term-cond:

Initial and terminal conditions
===============================

For most simulation exercises, it is necessary to provide initial (and
possibly terminal) conditions. It is also necessary to provide initial
guess values for non-linear solvers. This section describes the
statements used for those purposes.

In many contexts (deterministic or stochastic), it is necessary to
compute the steady state of a non-linear model: ``initval`` then
specifies numerical initial values for the non-linear solver. The
command ``resid`` can be used to compute the equation residuals for
the given initial values.

Used in perfect foresight mode, the types of forward-looking models
for which Dynare was designed require both initial and terminal
conditions. Most often these initial and terminal conditions are
static equilibria, but not necessarily.

One typical application is to consider an economy at the equilibrium
at time 0, trigger a shock in first period, and study the trajectory
of return to the initial equilibrium. To do that, one needs
``initval`` and ``shocks`` (see :ref:`shocks-exo`).

Another one is to study how an economy, starting from arbitrary
initial conditions at time 0 converges towards equilibrium. In this
case models, the command ``histval`` permits to specify different
historical initial values for variables with lags for the periods
before the beginning of the simulation. Due to the design of Dynare,
in this case ``initval`` is used to specify the terminal conditions.

.. block:: initval ;
           initval(OPTIONS...);

    |br| The ``initval`` block has two main purposes: providing guess
    values for non-linear solvers in the context of perfect foresight
    simulations and providing guess values for steady state
    computations in both perfect foresight and stochastic
    simulations. Depending on the presence of ``histval`` and
    ``endval`` blocks it is also used for declaring the initial and
    terminal conditions in a perfect foresight simulation
    exercise. Because of this interaction of the meaning of an
    ``initval`` block with the presence of ``histval`` and ``endval``
    blocks in perfect foresight simulations, it is strongly
    recommended to check that the constructed ``oo_.endo_simul`` and
    ``oo_.exo_simul`` variables contain the desired values after
    running ``perfect_foresight_setup`` and before running
    ``perfect_foresight_solver``. In the presence of leads and lags,
    these subfields of the results structure will store the historical
    values for the lags in the first column/row and the terminal
    values for the leads in the last column/row.

    The ``initval`` block is terminated by ``end;`` and contains lines
    of the form:

           VARIABLE_NAME = EXPRESSION;


    |br| *In a deterministic (i.e. perfect foresight) model*

    First, it will fill both the ``oo_.endo_simul`` and
    ``oo_.exo_simul`` variables storing the endogenous and exogenous
    variables with the values provided by this block. If there are no
    other blocks present, it will therefore provide the initial and
    terminal conditions for all the endogenous and exogenous
    variables, because it will also fill the last column/row of these
    matrices. For the intermediate simulation periods it thereby
    provides the starting values for the solver. In the presence of a
    ``histval`` block (and therefore absence of an ``endval`` block),
    this ``histval`` block will provide/overwrite the historical
    values for the state variables (lags) by setting the first
    column/row of ``oo_.endo_simul`` and ``oo_.exo_simul``. This
    implies that the ``initval`` block in the presence of ``histval``
    only sets the terminal values for the variables with leads and
    provides initial values for the perfect foresight solver.

    Because of these various functions of ``initval`` it is often
    necessary to provide values for all the endogenous variables in an
    ``initval`` block. Initial and terminal conditions are strictly
    necessary for lagged/leaded variables, while feasible starting
    values are required for the solver. It is important to be aware
    that if some variables, endogenous or exogenous, are not mentioned
    in the ``initval`` block, a zero value is assumed. It is
    particularly important to keep this in mind when specifying
    exogenous variables using ``varexo`` that are not allowed to take
    on the value of zero, like e.g. TFP.

    Note that if the ``initval`` block is immediately followed by a
    ``steady`` command, its semantics are slightly changed. The
    ``steady`` command will compute the steady state of the model for
    all the endogenous variables, assuming that exogenous variables
    are kept constant at the value declared in the ``initval``
    block. These steady state values conditional on the declared
    exogenous variables are then written into ``oo_.endo_simul`` and
    take up the potential roles as historical and terminal conditions
    as well as starting values for the solver. An ``initval`` block
    followed by ``steady`` is therefore formally equivalent to an
    ``initval`` block with the specified values for the exogenous
    variables, and the endogenous variables set to the associated
    steady state values conditional on the exogenous variables.

    |br| *In a stochastic model*

    The main purpose of ``initval`` is to provide initial guess values
    for the non-linear solver in the steady state computation. Note
    that if the ``initval`` block is not followed by ``steady``, the
    steady state computation will still be triggered by subsequent
    commands (``stoch_simul``, ``estimation``...).

    It is not necessary to declare 0 as initial value for exogenous
    stochastic variables, since it is the only possible value.

    The subsequently computed steady state (not the initial values,
    use histval for this) will be used as the initial condition at all
    the periods preceeding the first simulation period for the three
    possible types of simulations in stochastic mode:

        * :comm:`stoch_simul`, if the ``periods`` option is specified.
        * :comm:`forecast` as the initial point at which the forecasts
          are computed.
        * :comm:`conditional_forecast` as the initial point at which
          the conditional forecasts are computed.

    To start simulations at a particular set of starting values that
    are not a computed steady state, use :bck:`histval`.

    *Options*

    .. option:: all_values_required

        Issues an error and stops processing the .mod file if there is
        at least one endogenous or exogenous variable that has not
        been set in the initval block.

    *Example*
        ::

            initval;
            c = 1.2;
            k = 12;
            x = 1;
            end;

            steady;


.. block:: endval ;
           endval (OPTIONS...);

    |br| This block is terminated by ``end;`` and contains lines of the form:

        VARIABLE_NAME = EXPRESSION;

    |br| The ``endval`` block makes only sense in a deterministic model and
    cannot be used together with ``histval``. Similar to the
    ``initval`` command, it will fill both the ``oo_.endo_simul`` and
    ``oo_.exo_simul`` variables storing the endogenous and exogenous
    variables with the values provided by this block. If no
    ``initval`` block is present, it will fill the whole matrices,
    therefore providing the initial and terminal conditions for all
    the endogenous and exogenous variables, because it will also fill
    the first and last column/row of these matrices. Due to also
    filling the intermediate simulation periods it will provide the
    starting values for the solver as well.

    If an ``initval`` block is present, ``initval`` will provide the
    historical values for the variables (if there are states/lags),
    while ``endval`` will fill the remainder of the matrices, thereby
    still providing *i*) the terminal conditions for variables
    entering the model with a lead and *ii*) the initial guess values
    for all endogenous variables at all the simulation dates for the
    perfect foresight solver.

    Note that if some variables, endogenous or exogenous, are NOT
    mentioned in the ``endval`` block, the value assumed is that of
    the last ``initval`` block or ``steady`` command (if
    present). Therefore, in contrast to ``initval``, omitted variables
    are not automatically assumed to be 0 in this case. Again, it is
    strongly recommended to check the constructed ``oo_.endo_simul``
    and ``oo_.exo_simul`` variables after running
    ``perfect_foresight_setup`` and before running
    ``perfect_foresight_solver`` to see whether the desired outcome
    has been achieved.

    Like ``initval``, if the ``endval`` block is immediately followed
    by a ``steady`` command, its semantics are slightly changed. The
    ``steady`` command will compute the steady state of the model for
    all the endogenous variables, assuming that exogenous variables
    are kept constant to the value declared in the ``endval``
    block. These steady state values conditional on the declared
    exogenous variables are then written into ``oo_.endo_simul`` and
    therefore take up the potential roles as historical and terminal
    conditions as well as starting values for the solver. An
    ``endval`` block followed by ``steady`` is therefore formally
    equivalent to an ``endval`` block with the specified values for
    the exogenous variables, and the endogenous variables set to the
    associated steady state values.

    *Options*

    .. option:: all_values_required

        See :opt:`all_values_required`.

    *Example*

        ::

            var c k;
            varexo x;

            initval;
            c = 1.2;
            k = 12;
            x = 1;
            end;

            steady;

            endval;
            c = 2;
            k = 20;
            x = 2;
            end;

            steady;

        The initial equilibrium is computed by ``steady`` conditional
        on ``x=1``, and the terminal one conditional on ``x=2``. The
        ``initval`` block sets the initial condition for ``k``, while
        the ``endval`` block sets the terminal condition for
        ``c``. The starting values for the perfect foresight solver
        are given by the ``endval`` block. A detailed explanation
        follows below the next example.

    *Example*

        ::

            var c k;
            varexo x;

            model;
            c + k - aa*x*k(-1)^alph - (1-delt)*k(-1);
            c^(-gam) - (1+bet)^(-1)*(aa*alph*x(+1)*k^(alph-1) + 1 - delt)*c(+1)^(-gam);
            end;

            initval;
            k = 12;
            end;

            endval;
            c = 2;
            x = 1.1;
            end;
            simul(periods=200);

        In this example, the problem is finding the optimal path for
        consumption and capital for the periods :math:`t=1` to
        :math:`T=200`, given the path of the exogenous technology
        level ``x``. ``c`` is a forward-looking variable and the
        exogenous variable ``x`` appears with a lead in the expected
        return of physical capital, so we need terminal conditions for
        them, while ``k`` is a purely backward-looking (state)
        variable, so we need an initial condition for it.

        Setting ``x=1.1`` in the ``endval`` block without a ``shocks``
        block implies that technology is at :math:`1.1` in :math:`t=1`
        and stays there forever, because ``endval`` is filling all
        entries of ``oo_.endo_simul`` and ``oo_.exo_simul`` except for
        the very first one, which stores the initial conditions and
        was set to :math:`0` by the ``initval`` block when not
        explicitly specifying a value for it.

        Because the law of motion for capital is backward-looking, we
        need an initial condition for ``k`` at time :math:`0`. Due to
        the presence of ``endval``, this cannot be done via a
        ``histval`` block, but rather must be specified in the
        ``initval`` block. Similarly, because the Euler equation is
        forward-looking, we need a terminal condition for ``c`` at
        :math:`t=201`, which is specified in the ``endval`` block.

        As can be seen, it is not necessary to specify ``c`` and ``x``
        in the ``initval`` block and ``k`` in the ``endval`` block,
        because they have no impact on the results. Due to the
        optimization problem in the first period being to choose
        ``c,k`` at :math:`t=1` given the predetermined capital stock
        ``k`` inherited from :math:`t=0` as well as the current and
        future values for technology ``x``, the values for ``c`` and
        ``x`` at time :math:`t=0` play no role. The same applies to
        the choice of ``c,k`` at time :math:`t=200`, which does not
        depend on ``k`` at :math:`t=201`. As the Euler equation shows,
        that choice only depends on current capital as well as future
        consumption ``c`` and technology ``x``, but not on future
        capital ``k``. The intuitive reason is that those variables
        are the consequence of optimization problems taking place in
        at periods :math:`t=0` and :math:`t=201`, respectively, which
        are not modeled here.

    *Example*

        ::

            initval;
            c = 1.2;
            k = 12;
            x = 1;
            end;

            endval;
            c = 2;
            k = 20;
            x = 1.1;
            end;

        In this example, initial conditions for the forward-looking
        variables ``x`` and ``c`` are provided, together with a
        terminal condition for the backward-looking variable ``k``. As
        shown in the previous example, these values will not affect
        the simulation results. Dynare simply takes them as given and
        basically assumes that there were realizations of exogenous
        variables and states that make those choices equilibrium
        values (basically initial/terminal conditions at the
        unspecified time periods :math:`t<0` and :math:`t>201`).

        The above example suggests another way of looking at the use
        of ``steady`` after ``initval`` and ``endval``. Instead of
        saying that the implicit unspecified conditions before and
        after the simulation range have to fit the initial/terminal
        conditions of the endogenous variables in those blocks, steady
        specifies that those conditions at :math:`t<0` and
        :math:`t>201` are equal to being at the steady state given the
        exogenous variables in the ``initval`` and ``endval``
        blocks. The endogenous variables at :math:`t=0` and
        :math:`t=201` are then set to the corresponding steady state
        equilibrium values.

        The fact that ``c`` at :math:`t=0` and ``k`` at :math:`t=201`
        specified in ``initval`` and ````endval`` are taken as given
        has an important implication for plotting the simulated vector
        for the endogenous variables, i.e. the rows of
        ``oo_.endo_simul``: this vector will also contain the initial
        and terminal conditions and thus is 202 periods long in the
        example. When you specify arbitrary values for the initial and
        terminal conditions for forward- and backward-looking
        variables, respectively, these values can be very far away
        from the endogenously determined values at :math:`t=1` and
        :math:`t=200`. While the values at :math:`t=0` and
        :math:`t=201` are unrelated to the dynamics for
        :math:`0<t<201`, they may result in strange-looking large
        jumps. In the example above, consumption will display a large
        jump from :math:`t=0` to :math:`t=1` and capital will jump
        from :math:`t=200` to :math:`t=201` when using :comm:`rplot`
        or manually plotting ``oo_.endo_val``.


.. block:: histval ;
           histval (OPTIONS...);

    |br| *In a deterministic perfect foresight context*

    In models with lags on more than one period, the ``histval`` block
    permits to specify different historical initial values for
    different periods of the state variables. In this case, the
    ``initval`` block takes over the role of specifying terminal
    conditions and starting values for the solver. Note that the
    ``histval`` block does not take non-state variables.

    This block is terminated by ``end;`` and contains lines of the form:

        VARIABLE_NAME(INTEGER) = EXPRESSION;

    |br| EXPRESSION is any valid expression returning a numerical value
    and can contain already initialized variable names.

    By convention in Dynare, period 1 is the first period of the
    simulation. Going backward in time, the first period before the
    start of the simulation is period 0, then period -1, and so on.

    State variables not initialized in the ``histval`` block are
    assumed to have a value of zero at period 0 and before. Note that
    ``histval`` cannot be followed by ``steady``.

    *Example*

        ::

            model;
            x=1.5*x(-1)-0.6*x(-2)+epsilon;
            log(c)=0.5*x+0.5*log(c(+1));
            end;

            histval;
            x(0)=-1;
            x(-1)=0.2;
            end;

            initval;
            c=1;
            x=1;
            end;

        In this example, ``histval`` is used to set the historical
        conditions for the two lags of the endogenous variable ``x``,
        stored in the first column of ``oo_.endo_simul``. The
        ``initval`` block is used to set the terminal condition for
        the forward looking variable ``c``, stored in the last column
        of ``oo_.endo_simul``. Moreover, the ``initval`` block defines
        the starting values for the perfect foresight solver for both
        endogenous variables ``c`` and ``x``.

    *In a stochastic simulation context*

    In the context of stochastic simulations, ``histval`` allows
    setting the starting point of those simulations in the state
    space. As for the case of perfect foresight simulations, all not
    explicitly specified variables are set to 0. Moreover, as only
    states enter the recursive policy functions, all values specified
    for control variables will be ignored. This can be used

        * In :comm:`stoch_simul`, if the ``periods`` option is
          specified. Note that this only affects the starting point
          for the simulation, but not for the impulse response
          functions. When using the :ref:`loglinear <logl>` option,
          the ``histval`` block nevertheless takes the unlogged
          starting values.
        * In :comm:`forecast` as the initial point at which the
          forecasts are computed. When using the :ref:`loglinear
          <logl>` option, the ``histval`` block nevertheless takes the
          unlogged starting values.
        * In :comm:`conditional_forecast` for a calibrated model as
          the initial point at which the conditional forecasts are
          computed. When using the :ref:`loglinear <logl>` option, the
          histval-block nevertheless takes the unlogged starting
          values.
        * In :comm:`Ramsey policy <ramsey_model>`, where it also
          specifies the values of the endogenous states at which the
          objective function of the planner is computed. Note that the
          initial values of the Lagrange multipliers associated with
          the planner’s problem cannot be set (see
          :ref:`planner_objective_value <plan-obj>`).

    *Options*

    .. option:: all_values_required

        See :opt:`all_values_required`.

    *Example*

        ::

            var x y;
            varexo e;

            model;
            x = y(-1)^alpha*y(-2)^(1-alpha)+e;

            end;

            initval;
            x = 1;
            y = 1;
            e = 0.5;
            end;

            steady;

            histval;
            y(0) = 1.1;
            y(-1) = 0.9;
            end;

            stoch_simul(periods=100);


.. command:: resid ;

    |br| This command will display the residuals of the static
    equations of the model, using the values given for the endogenous
    in the last ``initval`` or ``endval`` block (or the steady state
    file if you provided one, see :ref:`st-st`).

.. command:: initval_file (filename = FILENAME);

    |br| In a deterministic setup, this command is used to specify a
    path for all endogenous and exogenous variables. The length of
    these paths must be equal to the number of simulation periods,
    plus the number of leads and the number of lags of the model (for
    example, with 50 simulation periods, in a model with 2 lags and 1
    lead, the paths must have a length of 53). Note that these paths
    cover two different things:

        * The constraints of the problem, which are given by the path
          for exogenous and the initial and terminal values for
          endogenous
        * The initial guess for the non-linear solver, which is given
          by the path for endogenous variables for the simulation
          periods (excluding initial and terminal conditions)

    The command accepts three file formats:

        * M-file (extension ``.m``): for each endogenous and exogenous
          variable, the file must contain a row or column vector of
          the same name. Their length must be ``periods +
          M_.maximum_lag + M_.maximum_lead``
        * MAT-file (extension ``.mat``): same as for M-files.
        * Excel file (extension ``.xls`` or ``.xlsx``): for each
          endogenous and exogenous, the file must contain a column of
          the same name. NB: Octave only supports the ``.xlsx`` file
          extension and must have the `io`_ package installed (easily
          done via octave by typing ‘``pkg install -forge io``’).

    .. warning:: The extension must be omitted in the command
                 argument. Dynare will automatically figure out the
                 extension and select the appropriate file type.


.. command:: histval_file (filename = FILENAME);

    |br| This command is equivalent to ``histval``, except that it
    reads its input from a file, and is typically used in conjunction
    with ``smoother2histval``.


.. _shocks-exo:

Shocks on exogenous variables
=============================

In a deterministic context, when one wants to study the transition of
one equilibrium position to another, it is equivalent to analyze the
consequences of a permanent shock and this in done in Dynare through
the proper use of ``initval`` and ``endval``.

Another typical experiment is to study the effects of a temporary
shock after which the system goes back to the original equilibrium (if
the model is stable...). A temporary shock is a temporary change of
value of one or several exogenous variables in the model. Temporary
shocks are specified with the command ``shocks``.

In a stochastic framework, the exogenous variables take random values
in each period. In Dynare, these random values follow a normal
distribution with zero mean, but it belongs to the user to specify the
variability of these shocks. The non-zero elements of the matrix of
variance-covariance of the shocks can be entered with the ``shocks``
command. Or, the entire matrix can be directly entered with
``Sigma_e`` (this use is however deprecated).

If the variance of an exogenous variable is set to zero, this variable
will appear in the report on policy and transition functions, but
isn’t used in the computation of moments and of Impulse Response
Functions. Setting a variance to zero is an easy way of removing an
exogenous shock.

Note that, by default, if there are several ``shocks`` or ``mshocks``
blocks in the same ``.mod`` file, then they are cumulative: all the
shocks declared in all the blocks are considered; however, if a
``shocks`` or ``mshocks`` block is declared with the ``overwrite``
option, then it replaces all the previous ``shocks`` and ``mshocks``
blocks.

.. block:: shocks ;
           shocks(overwrite);

    |br| See above for the meaning of the ``overwrite`` option.

    *In deterministic context*

    For deterministic simulations, the ``shocks`` block specifies
    temporary changes in the value of exogenous variables. For
    permanent shocks, use an ``endval`` block.

    The block should contain one or more occurrences of the following
    group of three lines::

      var VARIABLE_NAME;
      periods INTEGER[:INTEGER] [[,] INTEGER[:INTEGER]]...;
      values DOUBLE | (EXPRESSION)  [[,] DOUBLE | (EXPRESSION) ]...;

    It is possible to specify shocks which last several periods and
    which can vary over time. The ``periods`` keyword accepts a list
    of several dates or date ranges, which must be matched by as many
    shock values in the ``values`` keyword. Note that a range in the
    ``periods`` keyword can be matched by only one value in the
    ``values`` keyword. If ``values`` represents a scalar, the same
    value applies to the whole range. If ``values`` represents a
    vector, it must have as many elements as there are periods in the
    range.

    Note that shock values are not restricted to numerical constants:
    arbitrary expressions are also allowed, but you have to enclose
    them inside parentheses.

    *Example* (with scalar values)

    ::

        shocks;

        var e;
        periods 1;
        values 0.5;
        var u;
        periods 4:5;
        values 0;
        var v;
        periods 4:5 6 7:9;
        values 1 1.1 0.9;
        var w;
        periods 1 2;
        values (1+p) (exp(z));

        end;

    *Example* (with vector values)

    ::

        xx = [1.2; 1.3; 1];

        shocks;
        var e;
        periods 1:3;
        values (xx);
        end;

    |br| *In stochastic context*

    For stochastic simulations, the ``shocks`` block specifies the non
    zero elements of the covariance matrix of the shocks of exogenous
    variables.

    You can use the following types of entries in the block:

    * Specification of the standard error of an exogenous variable.

      ::

         var VARIABLE_NAME; stderr EXPRESSION;




    * Specification of the variance of an exogenous variable.

      ::

         var VARIABLE_NAME = EXPRESSION;


    * Specification the covariance of two exogenous variables.

      ::

         var VARIABLE_NAME, VARIABLE_NAME = EXPRESSION;

    * Specification of the correlation of two exogenous variables.

      ::

         corr VARIABLE_NAME, VARIABLE_NAME = EXPRESSION;

    In an estimation context, it is also possible to specify variances
    and covariances on endogenous variables: in that case, these
    values are interpreted as the calibration of the measurement
    errors on these variables. This requires the ``varobs`` command to
    be specified before the ``shocks`` block.

    *Example*

    ::

       shocks;
       var e = 0.000081;
       var u; stderr 0.009;
       corr e, u = 0.8;
       var v, w = 2;
       end;

    *Mixing deterministic and stochastic shocks*

    It is possible to mix deterministic and stochastic shocks to build
    models where agents know from the start of the simulation about
    future exogenous changes. In that case ``stoch_simul`` will
    compute the rational expectation solution adding future
    information to the state space (nothing is shown in the output of
    ``stoch_simul``) and ``forecast`` will compute a simulation
    conditional on initial conditions and future information.

    *Example*

    ::

        varexo_det tau;
        varexo e;
        ...
        shocks;
        var e; stderr 0.01;
        var tau;
        periods 1:9;
        values -0.15;
        end;

        stoch_simul(irf=0);

        forecast;

.. block:: mshocks ;
           mshocks(overwrite);

    |br| The purpose of this block is similar to that of the
    ``shocks`` block for deterministic shocks, except that the numeric
    values given will be interpreted in a multiplicative way. For
    example, if a value of ``1.05`` is given as shock value for some
    exogenous at some date, it means 5% above its steady state value
    (as given by the last ``initval`` or ``endval`` block).

    The syntax is the same than ``shocks`` in a deterministic context.

    This command is only meaningful in two situations:

    * on exogenous variables with a non-zero steady state, in a
      deterministic setup,
    * on deterministic exogenous variables with a non-zero steady
      state, in a stochastic setup.

    See above for the meaning of the ``overwrite`` option.

.. specvar:: Sigma_e

    |br| This special variable specifies directly the covariance
    matrix of the stochastic shocks, as an upper (or lower) triangular
    matrix. Dynare builds the corresponding symmetric matrix. Each row
    of the triangular matrix, except the last one, must be terminated
    by a semi-colon ;. For a given element, an arbitrary *EXPRESSION*
    is allowed (instead of a simple constant), but in that case you
    need to enclose the expression in parentheses. The order of the
    covariances in the matrix is the same as the one used in the
    ``varexo`` declaration.

    *Example*

    ::

        varexo u, e;

        Sigma_e = [ 0.81 (phi*0.9*0.009);
                    0.000081];

    This sets the variance of ``u`` to 0.81, the variance of ``e`` to
    0.000081, and the correlation between ``e`` and ``u`` to ``phi``.

    .. warning:: **The use of this special variable is deprecated and
                 is strongly discouraged**. You should use a
                 ``shocks`` block instead.


Other general declarations
==========================

.. command:: dsample INTEGER [INTEGER];

    |br| Reduces the number of periods considered in subsequent output commands.

.. command:: periods INTEGER

    |br| This command is now deprecated (but will still work for older
    model files). It is not necessary when no simulation is performed
    and is replaced by an option ``periods`` in
    ``perfect_foresight_setup``, ``simul`` and ``stoch_simul``.

    This command sets the number of periods in the simulation. The
    periods are numbered from 1 to INTEGER. In perfect foresight
    simulations, it is assumed that all future events are perfectly
    known at the beginning of period 1.

    *Example*

    ::

       periods 100;


.. _st-st:

Steady state
============

There are two ways of computing the steady state (i.e. the static
equilibrium) of a model. The first way is to let Dynare compute the
steady state using a nonlinear Newton-type solver; this should work
for most models, and is relatively simple to use. The second way is to
give more guidance to Dynare, using your knowledge of the model, by
providing it with a method to compute the steady state, either using a
`steady_state_model` block or writing matlab routine.


Finding the steady state with Dynare nonlinear solver
-----------------------------------------------------

.. command:: steady ;
             steady (OPTIONS...);

    |br| This command computes the steady state of a model using a
    nonlinear Newton-type solver and displays it. When a steady state
    file is used ``steady`` displays the steady state and checks that
    it is a solution of the static model.

    More precisely, it computes the equilibrium value of the
    endogenous variables for the value of the exogenous variables
    specified in the previous ``initval`` or ``endval`` block.

    ``steady`` uses an iterative procedure and takes as initial guess
    the value of the endogenous variables set in the previous
    ``initval`` or ``endval`` block.

    For complicated models, finding good numerical initial values for
    the endogenous variables is the trickiest part of finding the
    equilibrium of that model. Often, it is better to start with a
    smaller model and add new variables one by one.

    *Options*

    .. option:: maxit = INTEGER

       Determines the maximum number of iterations used in the
       non-linear solver. The default value of ``maxit`` is 50.

    .. option:: tolf = DOUBLE

       Convergence criterion for termination based on the function
       value. Iteration will cease when the residuals are smaller than
       ``tolf``. Default: ``eps^(1/3)``

    .. option:: solve_algo = INTEGER

       Determines the non-linear solver to use. Possible values for the
       option are:

           ``0``

                Use ``fsolve`` (under MATLAB, only available if you
                have the Optimization Toolbox; always available under
                Octave).

           ``1``

                Use Dynare’s own nonlinear equation solver (a
                Newton-like algorithm with line-search).

           ``2``

                Splits the model into recursive blocks and solves each
                block in turn using the same solver as value 1.

           ``3``

                Use Chris Sims’ solver.

           ``4``

                Splits the model into recursive blocks and solves each
                block in turn using a trust-region solver with
                autoscaling.

           ``5``

                Newton algorithm with a sparse Gaussian elimination
                (SPE) (requires ``bytecode`` option, see
                :ref:`model-decl`).

           ``6``

                Newton algorithm with a sparse LU solver at each
                iteration (requires ``bytecode`` and/or ``block``
                option, see :ref:`model-decl`).

           ``7``

                Newton algorithm with a Generalized Minimal Residual
                (GMRES) solver at each iteration (requires ``bytecode``
                and/or ``block`` option, see :ref:`model-decl`; not
                available under Octave).

           ``8``

                Newton algorithm with a Stabilized Bi-Conjugate
                Gradient (BICGSTAB) solver at each iteration (requires
                bytecode and/or block option, see :ref:`model-decl`).

           ``9``

                Trust-region algorithm on the entire model.

           ``10``

                Levenberg-Marquardt mixed complementarity problem
                (LMMCP) solver (*Kanzow and Petra (2004)*).

           ``11``

                PATH mixed complementarity problem solver of *Ferris
                and Munson (1999)*. The complementarity conditions are
                specified with an ``mcp`` equation tag, see
                :opt:`lmmcp`. Dynare only provides the interface for
                using the solver. Due to licence restrictions, you have
                to download the solver’s most current version yourself
                from `http://pages.cs.wisc.edu/~ferris/path.html
                <http://pages.cs.wisc.edu/~ferris/path.html>`_ and
                place it in Matlab’s search path.

       |br| Default value is ``4``.

    .. option:: homotopy_mode = INTEGER

       Use a homotopy (or divide-and-conquer) technique to solve for
       the steady state. If you use this option, you must specify a
       ``homotopy_setup`` block. This option can take three possible
       values:

           ``1``

                In this mode, all the parameters are changed
                simultaneously, and the distance between the boundaries
                for each parameter is divided in as many intervals as
                there are steps (as defined by the ``homotopy_steps``
                option); the problem is solves as many times as there
                are steps.

           ``2``

                Same as mode ``1``, except that only one parameter is
                changed at a time; the problem is solved as many times
                as steps times number of parameters.

           ``3``

                Dynare tries first the most extreme values. If it
                fails to compute the steady state, the interval
                between initial and desired values is divided by two
                for all parameters. Every time that it is impossible
                to find a steady state, the previous interval is
                divided by two. When it succeeds to find a steady
                state, the previous interval is multiplied by two. In
                that last case ``homotopy_steps`` contains the maximum
                number of computations attempted before giving up.

    .. option:: homotopy_steps = INTEGER

       Defines the number of steps when performing a homotopy. See
       ``homotopy_mode`` option for more details.

    .. option:: homotopy_force_continue = INTEGER

       This option controls what happens when homotopy fails.

           ``0``

                ``steady`` fails with an error message

           ``1``

                ``steady`` keeps the values of the last homotopy step
                that was successful and continues. **BE CAREFUL**:
                parameters and/or exogenous variables are NOT at the
                value expected by the user

       |br| Default is ``0``.

    .. option:: nocheck

       Don’t check the steady state values when they are provided
       explicitly either by a steady state file or a
       ``steady_state_model`` block. This is useful for models with
       unit roots as, in this case, the steady state is not unique or
       doesn’t exist.

    .. option:: markowitz = DOUBLE

       Value of the Markowitz criterion, used to select the
       pivot. Only used when ``solve_algo = 5``. Default: 0.5.

    *Example*

    See :ref:`init-term-cond`.

After computation, the steady state is available in the following variable:

.. matvar:: oo_.steady_state

    Contains the computed steady state. Endogenous variables are
    ordered in order of declaration used in the ``var`` command (which
    is also the order used in ``M_.endo_names``).


.. block:: homotopy_setup ;

    This block is used to declare initial and final values when using
    a homotopy method. It is used in conjunction with the option
    ``homotopy_mode`` of the steady command.

    The idea of homotopy (also called divide-and-conquer by some
    authors) is to subdivide the problem of finding the steady state
    into smaller problems. It assumes that you know how to compute the
    steady state for a given set of parameters, and it helps you
    finding the steady state for another set of parameters, by
    incrementally moving from one to another set of parameters.

    The purpose of the ``homotopy_setup`` block is to declare the
    final (and possibly also the initial) values for the parameters or
    exogenous that will be changed during the homotopy. It should
    contain lines of the form::

        VARIABLE_NAME, EXPRESSION, EXPRESSION;

    This syntax specifies the initial and final values of a given
    parameter/exogenous.

    There is an alternative syntax::

        VARIABLE_NAME, EXPRESSION;

    Here only the final value is specified for a given
    parameter/exogenous; the initial value is taken from the
    preceeding ``initval`` block.

    A necessary condition for a successful homotopy is that Dynare
    must be able to solve the steady state for the initial
    parameters/exogenous without additional help (using the guess
    values given in the ``initval`` block).

    If the homotopy fails, a possible solution is to increase the
    number of steps (given in ``homotopy_steps`` option of
    ``steady``).

    *Example*

    In the following example, Dynare will first compute the steady
    state for the initial values (``gam=0.5`` and ``x=1``), and then
    subdivide the problem into 50 smaller problems to find the steady
    state for the final values (``gam=2`` and ``x=2``)::

         var c k;
         varexo x;

         parameters alph gam delt bet aa;
         alph=0.5;
         delt=0.02;
         aa=0.5;
         bet=0.05;

         model;
         c + k - aa*x*k(-1)^alph - (1-delt)*k(-1);
         c^(-gam) - (1+bet)^(-1)*(aa*alph*x(+1)*k^(alph-1) + 1 - delt)*c(+1)^(-gam);
         end;

         initval;
         x = 1;
         k = ((delt+bet)/(aa*x*alph))^(1/(alph-1));
         c = aa*x*k^alph-delt*k;
         end;

         homotopy_setup;
         gam, 0.5, 2;
         x, 2;
         end;

         steady(homotopy_mode = 1, homotopy_steps = 50);


Providing the steady state to Dynare
------------------------------------

If you know how to compute the steady state for your model, you can
provide a MATLAB/Octave function doing the computation instead of
using ``steady``. Again, there are two options for doing that:

  * The easiest way is to write a ``steady_state_model`` block, which
    is described below in more details. See also ``fs2000.mod`` in the
    ``examples`` directory for an example. The steady state file
    generated by Dynare will be called ``+FILENAME/steadystate.m.``

  * You can write the corresponding MATLAB function by hand. If your
    MOD-file is called ``FILENAME.mod``, the steady state file must be
    called ``FILENAME_steadystate.m``. See
    ``NK_baseline_steadystate.m`` in the examples directory for an
    example. This option gives a bit more flexibility (loops and
    conditional structures can be used), at the expense of a heavier
    programming burden and a lesser efficiency.

Note that both files allow to update parameters in each call of the
function. This allows for example to calibrate a model to a labor
supply of 0.2 in steady state by setting the labor disutility
parameter to a corresponding value (see ``NK_baseline_steadystate.m``
in the ``examples`` directory). They can also be used in estimation
where some parameter may be a function of an estimated parameter and
needs to be updated for every parameter draw. For example, one might
want to set the capital utilization cost parameter as a function of
the discount rate to ensure that capacity utilization is 1 in steady
state. Treating both parameters as independent or not updating one as
a function of the other would lead to wrong results. But this also
means that care is required. Do not accidentally overwrite your
parameters with new values as it will lead to wrong results.

.. block:: steady_state_model ;

    |br| When the analytical solution of the model is known, this command
    can be used to help Dynare find the steady state in a more
    efficient and reliable way, especially during estimation where the
    steady state has to be recomputed for every point in the parameter
    space.

    Each line of this block consists of a variable (either an
    endogenous, a temporary variable or a parameter) which is assigned
    an expression (which can contain parameters, exogenous at the
    steady state, or any endogenous or temporary variable already
    declared above). Each line therefore looks like::

        VARIABLE_NAME = EXPRESSION;

    Note that it is also possible to assign several variables at the
    same time, if the main function in the right hand side is a
    MATLAB/Octave function returning several arguments::

        [ VARIABLE_NAME, VARIABLE_NAME... ] = EXPRESSION;

    Dynare will automatically generate a steady state file (of the
    form ``+FILENAME/steadystate.m``) using the information provided
    in this block.

    *Steady state file for deterministic models*

    The ``steady_state_model`` block works also with deterministic
    models. An ``initval`` block and, when necessary, an ``endval``
    block, is used to set the value of the exogenous variables. Each
    ``initval`` or ``endval`` block must be followed by ``steady`` to
    execute the function created by ``steady_state_model`` and set the
    initial, respectively terminal, steady state.

    *Example*

        ::

            var m P c e W R k d n l gy_obs gp_obs y dA;
            varexo e_a e_m;

            parameters alp bet gam mst rho psi del;

            ...
            // parameter calibration, (dynamic) model declaration, shock calibration...
            ...

            steady_state_model;
              dA = exp(gam);
              gst = 1/dA; // A temporary variable
              m = mst;

              // Three other temporary variables
              khst = ( (1-gst*bet*(1-del)) / (alp*gst^alp*bet) )^(1/(alp-1));
              xist = ( ((khst*gst)^alp - (1-gst*(1-del))*khst)/mst )^(-1);
              nust = psi*mst^2/( (1-alp)*(1-psi)*bet*gst^alp*khst^alp );

              n  = xist/(nust+xist);
              P  = xist + nust;
              k  = khst*n;

              l  = psi*mst*n/( (1-psi)*(1-n) );
              c  = mst/P;
              d  = l - mst + 1;
              y  = k^alp*n^(1-alp)*gst^alp;
              R  = mst/bet;

              // You can use MATLAB functions which return several arguments
              [W, e] = my_function(l, n);

              gp_obs = m/dA;
              gy_obs = dA;
            end;

            steady;

.. _eq-tag-ss:

Replace some equations during steady state computations
-------------------------------------------------------

When there is no steady state file, Dynare computes the steady state
by solving the static model, i.e. the model from the ``.mod`` file
from which leads and lags have been removed.

In some specific cases, one may want to have more control over the way
this static model is created. Dynare therefore offers the possibility
to explicitly give the form of equations that should be in the static
model.

More precisely, if an equation is prepended by a ``[static]`` tag,
then it will appear in the static model used for steady state
computation, but that equation will not be used for other
computations. For every equation tagged in this way, you must tag
another equation with ``[dynamic]``: that equation will not be used
for steady state computation, but will be used for other computations.

This functionality can be useful on models with a unit root, where
there is an infinity of steady states. An equation (tagged
``[dynamic]``) would give the law of motion of the nonstationary
variable (like a random walk). To pin down one specific steady state,
an equation tagged ``[static]`` would affect a constant value to the
nonstationary variable. Another situation where the ``[static]`` tag
can be useful is when one has only a partial closed form solution for
the steady state.

*Example*

This is a trivial example with two endogenous variables. The second
equation takes a different form in the static model::

    var c k;
    varexo x;
    ...
    model;
    c + k - aa*x*k(-1)^alph - (1-delt)*k(-1);
    [dynamic] c^(-gam) - (1+bet)^(-1)*(aa*alph*x(+1)*k^(alph-1) + 1 - delt)*c(+1)^(-gam);
    [static] k = ((delt+bet)/(x*aa*alph))^(1/(alph-1));
    end;


Getting information about the model
===================================

.. command:: check ;
             check (solve_algo = INTEGER);

    |br| Computes the eigenvalues of the model linearized around the
    values specified by the last ``initval``, ``endval`` or ``steady``
    statement. Generally, the eigenvalues are only meaningful if the
    linearization is done around a steady state of the model. It is a
    device for local analysis in the neighborhood of this steady
    state.

    A necessary condition for the uniqueness of a stable equilibrium
    in the neighborhood of the steady state is that there are as many
    eigenvalues larger than one in modulus as there are forward
    looking variables in the system. An additional rank condition
    requires that the square submatrix of the right Schur vectors
    corresponding to the forward looking variables (jumpers) and to
    the explosive eigenvalues must have full rank.

    *Options*

    .. _solvalg:

    .. option:: solve_algo = INTEGER

        See :ref:`solve_algo <solvalg>`, for the possible values and
        their meaning.

    .. option:: qz_zero_threshold = DOUBLE

        Value used to test if a generalized eigenvalue is :math:`0/0`
        in the generalized Schur decomposition (in which case the
        model does not admit a unique solution). Default: ``1e-6``.

    *Output*

    ``check`` returns the eigenvalues in the global variable ``oo_.dr.eigval``.


.. matvar:: oo_.dr.eigval

    Contains the eigenvalues of the model, as computed by the ``check`` command.

.. command:: model_diagnostics ;

    |br| This command performs various sanity checks on the model, and
    prints a message if a problem is detected (missing variables at
    current period, invalid steady state, singular Jacobian of static
    model).

.. command:: model_info ;
             model_info (OPTIONS...);

    |br| This command provides information about:

    * The normalization of the model: an endogenous variable is
      attributed to each equation of the model;
    * The block structure of the model: for each block ``model_info``
      indicates its type, the equations number and endogenous
      variables belonging to this block.

    This command can only be used in conjunction with the ``block``
    option of the ``model`` block.

    There are five different types of blocks depending on the
    simulation method used:

    * ‘EVALUATE FORWARD’

      In this case the block contains only equations where
      endogenous variable attributed to the equation appears
      currently on the left hand side and where no forward looking
      endogenous variables appear. The block has the form:
      :math:`y_{j,t} = f_j(y_t, y_{t-1}, \ldots, y_{t-k})`.

    * ‘EVALUATE BACKWARD’

      The block contains only equations where endogenous variable
      attributed to the equation appears currently on the left hand
      side and where no backward looking endogenous variables
      appear. The block has the form: :math:`y_{j,t} = f_j(y_t,
      y_{t+1}, \ldots, y_{t+k})`.

    * ‘SOLVE BACKWARD x’

      The block contains only equations where endogenous variable
      attributed to the equation does not appear currently on the
      left hand side and where no forward looking endogenous
      variables appear. The block has the form: :math:`g_j(y_{j,t},
      y_t, y_{t-1}, \ldots, y_{t-k})=0`. x is equal to ‘SIMPLE’
      if the block has only one equation. If several equation
      appears in the block, x is equal to ‘COMPLETE’.

    * ‘SOLVE FORWARD x’

      The block contains only equations where endogenous variable
      attributed to the equation does not appear currently on the
      left hand side and where no backward looking endogenous
      variables appear. The block has the form: :math:`g_j(y_{j,t},
      y_t, y_{t+1}, \ldots, y_{t+k})=0`. x is equal to ‘SIMPLE’
      if the block has only one equation. If several equation
      appears in the block, x is equal to ‘COMPLETE’.

    * ‘SOLVE TWO BOUNDARIES x’

      The block contains equations depending on both forward and
      backward variables. The block looks like: :math:`g_j(y_{j,t},
      y_t, y_{t-1}, \ldots, y_{t-k} ,y_t, y_{t+1}, \ldots,
      y_{t+k})=0`. x is equal to ‘SIMPLE’ if the block has only
      one equation. If several equation appears in the block, x is
      equal to ‘COMPLETE’.

    *Options*

    .. option:: 'static'

       Prints out the block decomposition of the static
       model. Without ’static’ option model_info displays the block
       decomposition of the dynamic model.

    .. option:: 'incidence'

       Displays the gross incidence matrix and the reordered incidence
       matrix of the block decomposed model.


.. command:: print_bytecode_dynamic_model ;

    |br| Prints the equations and the Jacobian matrix of the dynamic
    model stored in the bytecode binary format file. Can only be used
    in conjunction with the ``bytecode`` option of the ``model``
    block.

.. command:: print_bytecode_static_model ;

    |br| Prints the equations and the Jacobian matrix of the static model
    stored in the bytecode binary format file. Can only be used in
    conjunction with the ``bytecode`` option of the ``model`` block.


.. _det-simul:

Deterministic simulation
========================

When the framework is deterministic, Dynare can be used for models
with the assumption of perfect foresight. Typically, the system is
supposed to be in a state of equilibrium before a period ‘1’ when the
news of a contemporaneous or of a future shock is learned by the
agents in the model. The purpose of the simulation is to describe the
reaction in anticipation of, then in reaction to the shock, until the
system returns to the old or to a new state of equilibrium. In most
models, this return to equilibrium is only an asymptotic phenomenon,
which one must approximate by an horizon of simulation far enough in
the future. Another exercise for which Dynare is well suited is to
study the transition path to a new equilibrium following a permanent
shock. For deterministic simulations, the numerical problem consists
of solving a nonlinar system of simultaneous equations in n endogenous
variables in T periods. Dynare offers several algorithms for solving
this problem, which can be chosen via the ``stack_solve_algo``
option. By default (``stack_solve_algo=0``), Dynare uses a Newton-type
method to solve the simultaneous equation system. Because the
resulting Jacobian is in the order of ``n`` by ``T`` and hence will be
very large for long simulations with many variables, Dynare makes use
of the sparse matrix capacities of MATLAB/Octave. A slower but
potentially less memory consuming alternative (``stack_solve_algo=6``)
is based on a Newton-type algorithm first proposed by *Laffargue
(1990)* and *Boucekkine (1995)*, which uses relaxation
techniques. Thereby, the algorithm avoids ever storing the full
Jacobian. The details of the algorithm can be found in *Juillard
(1996)*. The third type of algorithms makes use of block decomposition
techniques (divide-and-conquer methods) that exploit the structure of
the model. The principle is to identify recursive and simultaneous
blocks in the model structure and use this information to aid the
solution process. These solution algorithms can provide a significant
speed-up on large models.

.. command:: perfect_foresight_setup ;
             perfect_foresight_setup (OPTIONS...);

    |br| Prepares a perfect foresight simulation, by extracting the
    information in the ``initval``, ``endval`` and ``shocks`` blocks
    and converting them into simulation paths for exogenous and
    endogenous variables.

    This command must always be called before running the simulation
    with ``perfect_foresight_solver``.

    *Options*

    .. option:: periods = INTEGER

       Number of periods of the simulation.

    .. option:: datafile = FILENAME

       If the variables of the model are not constant over time, their
       initial values, stored in a text file, could be loaded, using
       that option, as initial values before a deterministic
       simulation.

    *Output*

    The paths for the exogenous variables are stored into
    ``oo_.exo_simul``.

    The initial and terminal conditions for the endogenous variables
    and the initial guess for the path of endogenous variables are
    stored into ``oo_.endo_simul``.


.. command:: perfect_foresight_solver ;
             perfect_foresight_solver (OPTIONS...);

    |br| Computes the perfect foresight (or deterministic) simulation
    of the model.

    Note that ``perfect_foresight_setup`` must be called before this
    command, in order to setup the environment for the simulation.

    *Options*

    .. option:: maxit = INTEGER

       Determines the maximum number of iterations used in the
       non-linear solver. The default value of ``maxit`` is ``50``.

    .. option:: tolf = DOUBLE

       Convergence criterion for termination based on the function
       value. Iteration will cease when it proves impossible to
       improve the function value by more than ``tolf``. Default:
       ``1e-5``

    .. option:: tolx = DOUBLE

       Convergence criterion for termination based on the change in
       the function argument. Iteration will cease when the solver
       attempts to take a step that is smaller than ``tolx``. Default:
       ``1e-5``

    .. option:: stack_solve_algo = INTEGER

       Algorithm used for computing the solution. Possible values are:

           ``0``

               Newton method to solve simultaneously all the equations for
               every period, using sparse matrices (Default).

           ``1``

               Use a Newton algorithm with a sparse LU solver at each
               iteration (requires ``bytecode`` and/or ``block``
               option, see :ref:`model-decl`).

           ``2``

               Use a Newton algorithm with a Generalized Minimal
               Residual (GMRES) solver at each iteration (requires
               ``bytecode`` and/or ``block`` option, see
               :ref:`model-decl`; not available under Octave)

           ``3``

               Use a Newton algorithm with a Stabilized Bi-Conjugate
               Gradient (BICGSTAB) solver at each iteration (requires
               ``bytecode`` and/or ``block`` option, see
               :ref:`model-decl`).

           ``4``

               Use a Newton algorithm with a optimal path length at
               each iteration (requires ``bytecode`` and/or ``block``
               option, see :ref:`model-decl`).

           ``5``

               Use a Newton algorithm with a sparse Gaussian
               elimination (SPE) solver at each iteration (requires
               ``bytecode`` option, see :ref:`model-decl`).

           ``6``

               Use the historical algorithm proposed in *Juillard
               (1996)*: it is slower than ``stack_solve_algo=0``, but
               may be less memory consuming on big models (not
               available with ``bytecode`` and/or ``block`` options).

           ``7``

               Allows the user to solve the perfect foresight model
               with the solvers available through option
               ``solve_algo`` (See :ref:`solve_algo <solvalg>` for a
               list of possible values, note that values 5, 6, 7 and
               8, which require ``bytecode`` and/or ``block`` options,
               are not allowed). For instance, the following
               commands::

                    perfect_foresight_setup(periods=400);
                    perfect_foresight_solver(stack_solve_algo=7, solve_algo=9)

               trigger the computation of the solution with a trust
               region algorithm.

    .. option:: robust_lin_solve

       Triggers the use of a robust linear solver for the default
       ``stack_solve_algo=0``.

    .. option:: solve_algo

       See :ref:`solve_algo <solvalg>`. Allows selecting the solver
       used with ``stack_solve_algo=7``.

    .. option:: no_homotopy

       By default, the perfect foresight solver uses a homotopy
       technique if it cannot solve the problem. Concretely, it
       divides the problem into smaller steps by diminishing the size
       of shocks and increasing them progressively until the problem
       converges. This option tells Dynare to disable that
       behavior. Note that the homotopy is not implemented for purely
       forward or backward models.

    .. option:: markowitz = DOUBLE

       Value of the Markowitz criterion, used to select the
       pivot. Only used when ``stack_solve_algo = 5``. Default:
       ``0.5``.

    .. option:: minimal_solving_periods = INTEGER

       Specify the minimal number of periods where the model has to be
       solved, before using a constant set of operations for the
       remaining periods. Only used when ``stack_solve_algo =
       5``. Default: ``1``.

    .. option:: lmmcp

       Solves the perfect foresight model with a Levenberg-Marquardt
       mixed complementarity problem (LMMCP) solver (*Kanzow and Petra
       (2004)*), which allows to consider inequality constraints on
       the endogenous variables (such as a ZLB on the nominal interest
       rate or a model with irreversible investment). This option is
       equivalent to ``stack_solve_algo=7`` **and**
       ``solve_algo=10``. Using the LMMCP solver requires a particular
       model setup as the goal is to get rid of any min/max operators
       and complementary slackness conditions that might introduce a
       singularity into the Jacobian. This is done by attaching an
       equation tag (see :ref:`model-decl`) with the ``mcp`` keyword
       to affected equations. This tag states that the equation to
       which the tag is attached has to hold unless the expression
       within the tag is binding. For instance, a ZLB on the nominal
       interest rate would be specified as follows in the model
       block::

            model;
               ...
               [mcp = 'r > -1.94478']
               r = rho*r(-1) + (1-rho)*(gpi*Infl+gy*YGap) + e;
               ...
            end;

       where ``1.94478`` is the steady state level of the nominal
       interest rate and ``r`` is the nominal interest rate in
       deviation from the steady state. This construct implies that
       the Taylor rule is operative, unless the implied interest rate
       ``r<=-1.94478``, in which case the ``r`` is fixed at
       ``-1.94478`` (thereby being equivalent to a complementary
       slackness condition). By restricting the value of ``r`` coming
       out of this equation, the ``mcp`` tag also avoids using
       ``max(r,-1.94478)`` for other occurrences of ``r`` in the rest
       of the model. It is important to keep in mind that, because the
       ``mcp`` tag effectively replaces a complementary slackness
       condition, it cannot be simply attached to any
       equation. Rather, it must be attached to the correct affected
       equation as otherwise the solver will solve a different problem
       than originally intended.

       Note that in the current implementation, the content of the
       ``mcp`` equation tag is not parsed by the preprocessor. The
       inequalities must therefore be as simple as possible: an
       endogenous variable, followed by a relational operator,
       followed by a number (not a variable, parameter or expression).

    .. option:: endogenous_terminal_period

       The number of periods is not constant across Newton iterations
       when solving the perfect foresight model. The size of the
       nonlinear system of equations is reduced by removing the
       portion of the paths (and associated equations) for which the
       solution has already been identified (up to the tolerance
       parameter). This strategy can be interpreted as a mix of the
       shooting and relaxation approaches. Note that round off errors
       are more important with this mixed strategy (user should check
       the reported value of the maximum absolute error). Only
       available with option ``stack_solve_algo==0``.


    .. option:: linear_approximation

       Solves the linearized version of the perfect foresight
       model. The model must be stationary. Only available with option
       ``stack_solve_algo==0``.

    *Output*

    The simulated endogenous variables are available in global matrix
    ``oo_.endo_simul``.


.. command:: simul ;
             simul (OPTIONS...);

    |br| Short-form command for triggering the computation of a
    deterministic simulation of the model. It is strictly equivalent
    to a call to ``perfect_foresight_setup`` followed by a call to
    ``perfect_foresight_solver``.

    *Options*

    Accepts all the options of ``perfect_foresight_setup`` and
    ``perfect_foresight_solver``.

.. matvar:: oo_.endo_simul

    |br| This variable stores the result of a deterministic simulation
    (computed by ``perfect_foresight_solver`` or ``simul``) or of a
    stochastic simulation (computed by ``stoch_simul`` with the
    periods option or by ``extended_path``). The variables are
    arranged row by row, in order of declaration (as in
    ``M_.endo_names``). Note that this variable also contains initial
    and terminal conditions, so it has more columns than the value of
    ``periods`` option.

.. matvar:: oo_.exo_simul

    |br| This variable stores the path of exogenous variables during a
    simulation (computed by ``perfect_foresight_solver``, ``simul``,
    ``stoch_simul`` or ``extended_path``). The variables are arranged
    in columns, in order of declaration (as in
    ``M_.exo_names``). Periods are in rows. Note that this convention
    regarding columns and rows is the opposite of the convention for
    ``oo_.endo_simul``!


.. _stoch-sol:

Stochastic solution and simulation
==================================

In a stochastic context, Dynare computes one or several simulations
corresponding to a random draw of the shocks.

The main algorithm for solving stochastic models relies on a Taylor
approximation, up to third order, of the expectation functions (see
*Judd (1996)*, *Collard and Juillard (2001a)*, *Collard and Juillard
(2001b)*, and *Schmitt-Grohé and Uríbe (2004)*). The details of the
Dynare implementation of the first order solution are given in
*Villemot (2011)*. Such a solution is computed using the
``stoch_simul`` command.

As an alternative, it is possible to compute a simulation to a
stochastic model using the *extended path* method presented by *Fair
and Taylor (1983)*. This method is especially useful when there are
strong nonlinearities or binding constraints. Such a solution is
computed using the ``extended_path`` command.


Computing the stochastic solution
---------------------------------

.. command:: stoch_simul [VARIABLE_NAME...];
             stoch_simul (OPTIONS...) [VARIABLE_NAME...];

    |br| Solves a stochastic (i.e. rational expectations) model, using
    perturbation techniques.

    More precisely, ``stoch_simul`` computes a Taylor approximation of
    the model around the deterministic steady state and solves of the
    the decision and transition functions for the approximated
    model. Using this, it computes impulse response functions and
    various descriptive statistics (moments, variance decomposition,
    correlation and autocorrelation coefficients). For correlated
    shocks, the variance decomposition is computed as in the VAR
    literature through a Cholesky decomposition of the covariance
    matrix of the exogenous variables. When the shocks are correlated,
    the variance decomposition depends upon the order of the variables
    in the ``varexo`` command.

    The Taylor approximation is computed around the steady state (see
    :ref:`st-st`).

    The IRFs are computed as the difference between the trajectory of
    a variable following a shock at the beginning of period ``1`` and
    its steady state value. More details on the computation of IRFs
    can be found on the `Dynare wiki`_.

    Variance decomposition, correlation, autocorrelation are only
    displayed for variables with strictly positive variance. Impulse
    response functions are only plotted for variables with response
    larger than :math:`10^{-10}`.

    Variance decomposition is computed relative to the sum of the
    contribution of each shock. Normally, this is of course equal to
    aggregate variance, but if a model generates very large variances,
    it may happen that, due to numerical error, the two differ by a
    significant amount. Dynare issues a warning if the maximum
    relative difference between the sum of the contribution of each
    shock and aggregate variance is larger than ``0.01%``.

    The covariance matrix of the shocks is specified with the
    ``shocks`` command (see :ref:`shocks-exo`).

    When a list of ``VARIABLE_NAME`` is specified, results are
    displayed only for these variables.

    The ``stoch_simul`` command with a first order approximation can
    benefit from the block decomposition of the model (see
    :opt:`block`).

    *Options*

    .. option:: ar = INTEGER

       Order of autocorrelation coefficients to compute and to
       print. Default: ``5``.

    .. option:: drop = INTEGER

       Number of points (burnin) dropped at the beginning of
       simulation before computing the summary statistics. Note that
       this option does not affect the simulated series stored in
       ``oo_.endo_simul`` and the workspace. Here, no periods are
       dropped. Default: ``100``.

    .. option:: hp_filter = DOUBLE

       Uses HP filter with :math:`\lambda =` ``DOUBLE`` before
       computing moments. If theoretical moments are requested, the
       spectrum of the model solution is filtered following the
       approach outlined in Uhlig (2001). Default: no filter.

    .. option:: one_sided_hp_filter = DOUBLE

       Uses the one-sided HP filter with :math:`\lambda =` ``DOUBLE``
       described in *Stock and Watson (1999)* before computing
       moments. This option is only available with simulated
       moments. Default: no filter.

    .. option:: hp_ngrid = INTEGER

       Number of points in the grid for the discrete Inverse Fast
       Fourier Transform used in the HP filter computation. It may be
       necessary to increase it for highly autocorrelated
       processes. Default: ``512``.

    .. option:: bandpass_filter

       Uses a bandpass filter with the default passband before
       computing moments. If theoretical moments are requested, the
       spectrum of the model solution is filtered using an ideal
       bandpass filter. If empirical moments are requested, the
       *Baxter and King (1999)* filter is used. Default: no filter.

    .. option:: bandpass_filter = [HIGHEST_PERIODICITY LOWEST_PERIODICITY]

       Uses a bandpass filter before computing moments. The passband
       is set to a periodicity of to LOWEST_PERIODICITY,
       e.g. :math:`6` to :math:`32` quarters if the model frequency is
       quarterly. Default: ``[6,32]``.

    .. option:: irf = INTEGER

       Number of periods on which to compute the IRFs. Setting
       ``irf=0`` suppresses the plotting of IRFs. Default: ``40``.

    .. option:: irf_shocks = ( VARIABLE_NAME [[,] VARIABLE_NAME ...] )

       The exogenous variables for which to compute IRFs. Default:
       all.

    .. option:: relative_irf

       Requests the computation of normalized IRFs. At first order,
       the normal shock vector of size one standard deviation is
       divided by the standard deviation of the current shock and
       multiplied by 100. The impulse responses are hence the
       responses to a unit shock of size 1 (as opposed to the regular
       shock size of one standard deviation), multiplied by 100. Thus,
       for a loglinearized model where the variables are measured in
       percent, the IRFs have the interpretation of the percent
       responses to a 100 percent shock. For example, a response of
       400 of output to a TFP shock shows that output increases by 400
       percent after a 100 percent TFP shock (you will see that TFP
       increases by 100 on impact). Given linearity at ``order=1``, it
       is straightforward to rescale the IRFs stored in ``oo_.irfs``
       to any desired size. At higher order, the interpretation is
       different. The ``relative_irf`` option then triggers the
       generation of IRFs as the response to a 0.01 unit shock
       (corresponding to 1 percent for shocks measured in percent) and
       no multiplication with 100 is performed. That is, the normal
       shock vector of size one standard deviation is divided by the
       standard deviation of the current shock and divided by 100. For
       example, a response of 0.04 of log output (thus measured in
       percent of the steady state output level) to a TFP shock also
       measured in percent then shows that output increases by 4
       percent after a 1 percent TFP shock (you will see that TFP
       increases by 0.01 on impact).

    .. option:: irf_plot_threshold = DOUBLE

       Threshold size for plotting IRFs. All IRFs for a particular
       variable with a maximum absolute deviation from the steady
       state smaller than this value are not displayed. Default:
       ``1e-10``.

    .. option:: nocorr

       Don’t print the correlation matrix (printing them is the
       default).

    .. option:: nodecomposition

       Don’t compute (and don’t print) unconditional variance
       decomposition.

    .. option:: nofunctions

       Don’t print the coefficients of the approximated solution
       (printing them is the default).

    .. option:: nomoments

       Don’t print moments of the endogenous variables (printing them
       is the default).

    .. option:: nograph

       Do not create graphs (which implies that they are not saved to
       the disk nor displayed). If this option is not used, graphs
       will be saved to disk (to the format specified by
       ``graph_format`` option, except if ``graph_format=none``) and
       displayed to screen (unless ``nodisplay`` option is used).

    .. option:: graph

       Re-enables the generation of graphs previously shut off with
       ``nograph``.

    .. option:: nodisplay

       Do not display the graphs, but still save them to disk (unless
       ``nograph`` is used).

    .. option:: graph_format = FORMAT
                graph_format = ( FORMAT, FORMAT... )

       Specify the file format(s) for graphs saved to disk. Possible
       values are ``eps`` (the default), ``pdf``, ``fig`` and ``none``
       (under Octave, only ``eps`` and ``none`` are available). If the
       file format is set equal to ``none``, the graphs are displayed
       but not saved to the disk.

    .. option:: noprint

       Don’t print anything. Useful for loops.

    .. option:: print

       Print results (opposite of ``noprint``).

    .. option:: order = INTEGER

       Order of Taylor approximation. Acceptable values are ``1``,
       ``2`` and ``3``. Note that for third order, ``k_order_solver``
       option is implied and only empirical moments are available (you
       must provide a value for ``periods`` option). Default: ``2``
       (except after an ``estimation`` command, in which case the
       default is the value used for the estimation).

    .. option:: k_order_solver

       Use a k-order solver (implemented in C++) instead of the
       default Dynare solver. This option is not yet compatible with
       the ``bytecode`` option (see :ref:`model-decl`). Default:
       disabled for order 1 and 2, enabled otherwise.

    .. option:: periods = INTEGER

       If different from zero, empirical moments will be computed
       instead of theoretical moments. The value of the option
       specifies the number of periods to use in the
       simulations. Values of the initval block, possibly recomputed
       by ``steady``, will be used as starting point for the
       simulation. The simulated endogenous variables are made
       available to the user in a vector for each variable and in the
       global matrix ``oo_.endo_simul`` (see
       :mvar:`oo_.endo_simul`). The simulated exogenous variables are
       made available in ``oo_.exo_simul`` (see
       :mvar:`oo_.exo_simul`). Default: ``0``.

    .. option:: qz_criterium = DOUBLE

       Value used to split stable from unstable eigenvalues in
       reordering the Generalized Schur decomposition used for solving
       1st-order problems. Default: ``1.000001`` (except when
       estimating with ``lik_init`` option equal to ``1``: the default
       is ``0.999999`` in that case; see :ref:`estim`).

    .. option:: qz_zero_threshold = DOUBLE

       See :opt:`qz_zero_threshold <qz_zero_threshold = DOUBLE>`.

    .. option:: replic = INTEGER

       Number of simulated series used to compute the IRFs. Default:
       ``1`` if ``order=1``, and ``50`` otherwise.

    .. option:: simul_replic = INTEGER

       Number of series to simulate when empirical moments are
       requested (i.e. ``periods`` :math:`>` 0). Note that if this
       option is greater than 1, the additional series will not be
       used for computing the empirical moments but will simply be
       saved in binary form to the file ``FILENAME_simul``. Default:
       ``1``.

    .. option:: solve_algo = INTEGER

       See :ref:`solve_algo <solvalg>`, for the possible values and
       their meaning.

    .. option:: aim_solver

       Use the Anderson-Moore Algorithm (AIM) to compute the decision
       rules, instead of using Dynare’s default method based on a
       generalized Schur decomposition. This option is only valid for
       first order approximation. See `AIM website`_ for more details
       on the algorithm.

    .. option:: conditional_variance_decomposition = INTEGER
                conditional_variance_decomposition = [INTEGER1:INTEGER2]
                conditional_variance_decomposition = [INTEGER1 INTEGER2 ...]

       Computes a conditional variance decomposition for the specified
       period(s). The periods must be strictly positive. Conditional
       variances are given by :math:`var(y_{t+k}\vert t)`. For period
       1, the conditional variance decomposition provides the
       decomposition of the effects of shocks upon impact.

       The results are stored in
       ``oo_.conditional_variance_decomposition`` (see
       :mvar:`oo_.conditional_variance_decomposition`). In the
       presence of measurement error, the
       ``oo_.conditional_variance_decomposition`` field will contain
       the variance contribution after measurement error has been
       taken out, i.e. the decomposition will be conducted of the
       actual as opposed to the measured variables. The variance
       decomposition of the measured variables will be stored in
       ``oo_.conditional_variance_decomposition_ME`` (see
       :mvar:`oo_.conditional_variance_decomposition_ME`).  The
       variance decomposition is only conducted, if theoretical
       moments are requested, *i.e.* using the ``periods=0``-option.
       In case of ``order=2``, Dynare provides a second-order accurate
       approximation to the true second moments based on the linear
       terms of the second-order solution (see *Kim, Kim,
       Schaumburg and Sims (2008)*). Note that the unconditional
       variance decomposition *i.e.* at horizon infinity) is
       automatically conducted if theoretical moments are requested
       and if ``nodecomposition`` is not set (see
       :mvar:`oo_.variance_decomposition`).

    .. option:: pruning

       Discard higher order terms when iteratively computing
       simulations of the solution. At second order, Dynare uses the
       algorithm of *Kim, Kim, Schaumburg and Sims (2008)*, while at
       third order its generalization by *Andreasen,
       Fernández-Villaverde and Rubio-Ramírez (2013)* is used.

    .. option:: partial_information

       Computes the solution of the model under partial information,
       along the lines of *Pearlman, Currie and Levine (1986)*. Agents
       are supposed to observe only some variables of the economy. The
       set of observed variables is declared using the ``varobs``
       command. Note that if ``varobs`` is not present or contains all
       endogenous variables, then this is the full information case
       and this option has no effect. More references can be found
       `here <http://www.dynare.org/DynareWiki/PartialInformation>`_ .

    .. option:: sylvester = OPTION

       Determines the algorithm used to solve the Sylvester equation
       for block decomposed model. Possible values for OPTION are:

           ``default``

                Uses the default solver for Sylvester equations
                (``gensylv``) based on Ondra Kamenik’s algorithm (see
                `here
                <http://www.dynare.org/documentation-and-support/dynarepp/sylvester.pdf/at_download/file>`_
                for more information).

           ``fixed_point``

                Uses a fixed point algorithm to solve the Sylvester
                equation (``gensylv_fp``). This method is faster than
                the default one for large scale models.

       |br| Default value is ``default``.

    .. option:: sylvester_fixed_point_tol = DOUBLE

       The convergence criterion used in the fixed point
       Sylvester solver. Its default value is ``1e-12``.

    .. option:: dr = OPTION

       Determines the method used to compute the decision
       rule. Possible values for OPTION are:

           ``default``

                Uses the default method to compute the decision rule
                based on the generalized Schur decomposition (see
                *Villemot (2011)* for more information).

           ``cycle_reduction``

                Uses the cycle reduction algorithm to solve the
                polynomial equation for retrieving the coefficients
                associated to the endogenous variables in the decision
                rule. This method is faster than the default one for
                large scale models.

           ``logarithmic_reduction``

                Uses the logarithmic reduction algorithm to solve the
                polynomial equation for retrieving the coefficients
                associated to the endogenous variables in the decision
                rule. This method is in general slower than the
                ``cycle_reduction``.

       |br| Default value is ``default``.

    .. option:: dr_cycle_reduction_tol = DOUBLE

       The convergence criterion used in the cycle reduction
       algorithm. Its default value is ``1e-7``.

    .. option:: dr_logarithmic_reduction_tol = DOUBLE

       The convergence criterion used in the logarithmic reduction
       algorithm. Its default value is ``1e-12``.

    .. option:: dr_logarithmic_reduction_maxiter = INTEGER

       The maximum number of iterations used in the logarithmic
       reduction algorithm. Its default value is ``100``.

    .. option:: loglinear

       See :ref:`loglinear <logl>`. Note that ALL variables are
       log-transformed by using the Jacobian transformation, not only
       selected ones. Thus, you have to make sure that your variables
       have strictly positive steady states. ``stoch_simul`` will
       display the moments, decision rules, and impulse responses for
       the log-linearized variables. The decision rules saved in
       ``oo_.dr`` and the simulated variables will also be the ones
       for the log-linear variables.

    .. option:: tex

       Requests the printing of results and graphs in TeX tables and
       graphics that can be later directly included in LaTeX
       files.

    .. option:: dr_display_tol = DOUBLE

       Tolerance for the suppression of small terms in the display of
       decision rules. Rows where all terms are smaller than
       ``dr_display_tol`` are not displayed. Default value: ``1e-6``.

    .. option:: contemporaneous_correlation

       Saves the contemporaneous correlation between the endogenous
       variables in ``oo_.contemporaneous_correlation``. Requires the
       ``nocorr`` option not to be set.

    .. option:: spectral_density

       Triggers the computation and display of the theoretical
       spectral density of the (filtered) model variables. Results are
       stored in ´´oo_.SpectralDensity´´, defined below. Default: do
       not request spectral density estimates.


    *Output*

    This command sets ``oo_.dr``, ``oo_.mean``, ``oo_.var`` and
    ``oo_.autocorr``, which are described below.

    If the ``periods`` option is present, sets ``oo_.skewness``,
    ``oo_.kurtosis``, and ``oo_.endo_simul`` (see
    :mvar:`oo_.endo_simul`), and also saves the simulated variables in
    MATLAB/Octave vectors of the global workspace with the same name
    as the endogenous variables.

    If option ``irf`` is different from zero, sets ``oo_.irfs`` (see
    below) and also saves the IRFs in MATLAB/Octave vectors of the
    global workspace (this latter way of accessing the IRFs is
    deprecated and will disappear in a future version).

    If the option ``contemporaneous_correlation`` is different from
    ``0``, sets ``oo_.contemporaneous_correlation``, which is
    described below.

    *Example*

        ::

            shocks;
            var e;
            stderr 0.0348;
            end;

            stoch_simul;

        Performs the simulation of the 2nd-order approximation of a
        model with a single stochastic shock ``e``, with a standard
        error of ``0.0348``.

    *Example*

        ::

            stoch_simul(irf=60) y k;

        Performs the simulation of a model and displays impulse
        response functions on 60 periods for variables ``y`` and
        ``k``.

.. matvar:: oo_.mean

   |br| After a run of ``stoch_simul``, contains the mean of the
   endogenous variables. Contains theoretical mean if the ``periods``
   option is not present, and simulated mean otherwise. The variables
   are arranged in declaration order.

.. matvar:: oo_.var

    |br| After a run of ``stoch_simul``, contains the
    variance-covariance of the endogenous variables. Contains
    theoretical variance if the ``periods`` option is not present (or
    an approximation thereof for ``order=2``), and simulated variance
    otherwise. The variables are arranged in declaration order.

.. matvar:: oo_.skewness

    |br| After a run of ``stoch_simul`` contains the skewness
    (standardized third moment) of the simulated variables if the
    ``periods`` option is present. The variables are arranged in
    declaration order.

.. matvar:: oo_.kurtosis

   |br| After a run of ``stoch_simul`` contains the kurtosis (standardized
   fourth moment) of the simulated variables if the ``periods`` option
   is present. The variables are arranged in declaration order.

.. matvar:: oo_.autocorr

    |br| After a run of ``stoch_simul``, contains a cell array of the
    autocorrelation matrices of the endogenous variables. The element
    number of the matrix in the cell array corresponds to the order of
    autocorrelation. The option ar specifies the number of
    autocorrelation matrices available. Contains theoretical
    autocorrelations if the ``periods`` option is not present (or an
    approximation thereof for ``order=2``), and simulated
    autocorrelations otherwise. The field is only created if
    stationary variables are present.

    The element ``oo_.autocorr{i}(k,l)`` is equal to the correlation
    between :math:`y^k_t` and :math:`y^l_{t-i}`, where :math:`y^k`
    (resp. :math:`y^l`) is the :math:`k`-th (resp. :math:`l`-th)
    endogenous variable in the declaration order.

    Note that if theoretical moments have been requested,
    ``oo_.autocorr{i}`` is the same than ``oo_.gamma_y{i+1}``.

.. matvar:: oo_.gamma_y

    |br| After a run of ``stoch_simul``, if theoretical moments have been
    requested (i.e. if the ``periods`` option is not present), this
    variable contains a cell array with the following values (where
    ``ar`` is the value of the option of the same name):

        ``oo_.gamma{1}``

            Variance/covariance matrix.

        ``oo_.gamma{i+1}`` (for i=1:ar)

            Autocorrelation function. See :mvar:`oo_.autocorr` for
            more details. **Beware**, this is the autocorrelation
            function, not the autocovariance function.

        ``oo_.gamma{nar+2}``

            Unconditional variance decomposition, see
            :mvar:`oo_.variance_decomposition`.

        ``oo_.gamma{nar+3}``

            If a second order approximation has been requested,
            contains the vector of the mean correction terms.

            In case ``order=2``, the theoretical second moments are a
            second order accurate approximation of the true second
            moments, see conditional_variance_decomposition.

.. matvar:: oo_.variance_decomposition

    |br| After a run of ``stoch_simul`` when requesting theoretical
    moments (``periods=0``), contains a matrix with the result of the
    unconditional variance decomposition (i.e. at horizon
    infinity). The first dimension corresponds to the endogenous
    variables (in the order of declaration after the command or in
    ``M_.endo_names``) and the second dimension corresponds to
    exogenous variables (in the order of declaration). Numbers are in
    percent and sum up to 100 across columns. In the presence of
    measurement error, the field will contain the variance
    contribution after measurement error has been taken out, *i.e.*
    the decomposition will be conducted of the actual as opposed to
    the measured variables.

.. matvar:: oo_.variance_decomposition_ME

    |br| Field set after a run of ``stoch_simul`` when requesting
    theoretical moments (``periods=0``) if measurement error is
    present. It is similar to :mvar:`oo_.variance_decomposition`, but
    the decomposition will be conducted of the measured variables. The
    field contains a matrix with the result of the unconditional
    variance decomposition (*i.e.* at horizon infinity). The first
    dimension corresponds to the observed endoogenous variables (in
    the order of declaration after the command) and the second
    dimension corresponds to exogenous variables (in the order of
    declaration), with the last column corresponding to the
    contribution of measurement error. Numbers are in percent and sum
    up to 100 across columns.

.. matvar:: oo_.conditional_variance_decomposition

    |br| After a run of ``stoch_simul`` with the
    ``conditional_variance_decomposition`` option, contains a
    three-dimensional array with the result of the decomposition. The
    first dimension corresponds to forecast horizons (as declared with
    the option), the second dimension corresponds to endogenous
    variables (in the order of declaration after the command or in
    ``M_.endo_names`` if not specified), the third dimension
    corresponds to exogenous variables (in the order of
    declaration). In the presence of measurement error, the field will
    contain the variance contribution after measurement error has been
    taken out, *i.e.* the decomposition will be conductedof the actual
    as opposed to the measured variables.

.. matvar:: oo_.conditional_variance_decomposition_ME

    |br| Field set after a run of ``stoch_simul`` with the
    ``conditional_variance_decomposition`` option if measurement error
    is present. It is similar to
    :mvar:`oo_.conditional_variance_decomposition`, but the
    decomposition will be conducted of the measured variables.  It
    contains a three-dimensional array with the result of the
    decomposition. The first dimension corresponds to forecast
    horizons (as declared with the option), the second dimension
    corresponds to observed endogenous variables (in the order of
    declaration), the third dimension corresponds to exogenous
    variables (in the order of declaration), with the last column
    corresponding to the contribution of the measurement error.

.. matvar:: oo_.contemporaneous_correlation

    |br| After a run of ``stoch_simul`` with the
    ``contemporaneous_correlation option``, contains theoretical
    contemporaneous correlations if the ``periods`` option is not
    present (or an approximation thereof for ``order=2``), and
    simulated contemporaneous correlations otherwise. The variables
    are arranged in declaration order.

.. matvar:: oo_.SpectralDensity

    |br| After a run of ``stoch_simul`` with option
    ``spectral_density``, contains the spectral density of the model
    variables. There will be a ``nvars`` by ``nfrequencies`` subfield
    ``freqs`` storing the respective frequency grid points ranging
    from :math:`0` to :math:`2\pi` and a same sized subfield
    ``density`` storing the corresponding density.

.. matvar:: oo_.irfs

   |br| After a run of ``stoch_simul`` with option ``irf`` different
   from zero, contains the impulse responses, with the following
   naming convention: `VARIABLE_NAME_SHOCK_NAME`.

    For example, ``oo_.irfs.gnp_ea`` contains the effect on ``gnp`` of
    a one-standard deviation shock on ``ea``.

The approximated solution of a model takes the form of a set of
decision rules or transition equations expressing the current value of
the endogenous variables of the model as function of the previous
state of the model and shocks observed at the beginning of the
period. The decision rules are stored in the structure ``oo_.dr``
which is described below.

.. command:: extended_path ;
             extended_path (OPTIONS...);

    |br| Simulates a stochastic (i.e. rational expectations) model,
    using the extended path method presented by *Fair and Taylor
    (1983)*. Time series for the endogenous variables are generated by
    assuming that the agents believe that there will no more shocks in
    the following periods.

    This function first computes a random path for the exogenous
    variables (stored in ``oo_.exo_simul``, see :mvar:`oo_.exo_simul`)
    and then computes the corresponding path for endogenous variables,
    taking the steady state as starting point. The result of the
    simulation is stored in ``oo_.endo_simul`` (see
    :mvar:`oo_.endo_simul`). Note that this simulation approach does
    not solve for the policy and transition equations but for paths
    for the endogenous variables.

    *Options*

    .. option:: periods = INTEGER

       The number of periods for which the simulation is to be
       computed. No default value, mandatory option.

    .. option:: solver_periods = INTEGER

       The number of periods used to compute the solution of the
       perfect foresight at every iteration of the algorithm. Default:
       ``200``.

    .. option:: order = INTEGER

       If order is greater than ``0`` Dynare uses a gaussian
       quadrature to take into account the effects of future
       uncertainty. If ``order`` :math:`=S` then the time series for
       the endogenous variables are generated by assuming that the
       agents believe that there will no more shocks after period
       :math:`t+S`. This is an experimental feature and can be quite
       slow. Default: ``0``.

    .. option:: hybrid

       Use the constant of the second order perturbation reduced form
       to correct the paths generated by the (stochastic) extended
       path algorithm.


Typology and ordering of variables
----------------------------------

Dynare distinguishes four types of endogenous variables:

*Purely backward (or purely predetermined) variables*

    Those that appear only at current and past period in the model,
    but not at future period (i.e. at :math:`t` and :math:`t-1` but
    not :math:`t+1`). The number of such variables is equal to
    ``M_.npred``.

*Purely forward variables*

    Those that appear only at current and future period in the model,
    but not at past period (i.e. at :math:`t` and :math:`t+1` but not
    :math:`t-1`). The number of such variables is stored in
    ``M_.nfwrd``.

*Mixed variables*

    Those that appear at current, past and future period in the model
    (i.e. at :math:`t`, :math:`t+1` and :math:`t-1`). The number of
    such variables is stored in ``M_.nboth``.

*Static variables*

    Those that appear only at current, not past and future period in
    the model (i.e. only at :math:`t`, not at :math:`t+1` or
    :math:`t-1`). The number of such variables is stored in
    ``M_.nstatic``.

Note that all endogenous variables fall into one of these four
categories, since after the creation of auxiliary variables (see
:ref:`aux-variables`), all endogenous have at most one lead and one
lag. We therefore have the following identity:

    .. code-block:: matlab

       M_.npred + M_.both + M_.nfwrd + M_.nstatic = M_.endo_nbr

Internally, Dynare uses two orderings of the endogenous variables: the
order of declaration (which is reflected in ``M_.endo_names``), and an
order based on the four types described above, which we will call the
DR-order (“DR” stands for decision rules). Most of the time, the
declaration order is used, but for elements of the decision rules, the
DR-order is used.

The DR-order is the following: static variables appear first, then
purely backward variables, then mixed variables, and finally purely
forward variables. Inside each category, variables are arranged
according to the declaration order.

Variable ``oo_.dr.order_var`` maps DR-order to declaration order, and
variable ``oo_.dr.inv_order_var`` contains the inverse map. In other
words, the k-th variable in the DR-order corresponds to the endogenous
variable numbered ``oo_.dr_order_var(k)`` in declaration
order. Conversely, k-th declared variable is numbered
``oo_.dr.inv_order_var(k)`` in DR-order.

Finally, the state variables of the model are the purely backward
variables and the mixed variables. They are ordered in DR-order when
they appear in decision rules elements. There are ``M_.nspred =
M_.npred + M_.nboth`` such variables. Similarly, one has ``M_.nsfwrd =
M_.nfwrd + M_.nboth``, and ``M_.ndynamic = M_.nfwrd + M_.nboth +
M_.npred``.


First-order approximation
-------------------------

The approximation has the stylized form:

.. math::

       y_t = y^s + A y^h_{t-1} + B u_t

where :math:`y^s` is the steady state value of :math:`y` and
:math:`y^h_t=y_t-y^s`.

The coefficients of the decision rules are stored as follows:

* :math:`y^s` is stored in ``oo_.dr.ys``. The vector rows correspond
  to all endogenous in the declaration order.
* :math:`A` is stored in ``oo_.dr.ghx``. The matrix rows correspond to
  all endogenous in DR-order. The matrix columns correspond to state
  variables in DR-order.
* :math:`B` is stored ``oo_.dr.ghu``. The matrix rows correspond to
  all endogenous in DR-order. The matrix columns correspond to
  exogenous variables in declaration order.

Of course, the shown form of the approximation is only stylized,
because it neglects the required different ordering in :math:`y^s` and
:math:`y^h_t`. The precise form of the approximation that shows the
way Dynare deals with differences between declaration and DR-order, is

.. math::

       y_t(\mathrm{oo\_.dr.order\_var}) =
       y^s(\mathrm{oo\_.dr.order\_var}) + A \cdot
       y_{t-1}(\mathrm{oo\_.dr.order\_var(k2)}) -
       y^s(\mathrm{oo\_.dr.order\_var(k2)}) + B\cdot u_t

where :math:`\mathrm{k2}` selects the state variables, :math:`y_t` and
:math:`y^s` are in declaration order and the coefficient matrices are
in DR-order. Effectively, all variables on the right hand side are
brought into DR order for computations and then assigned to
:math:`y_t` in declaration order.


Second-order approximation
--------------------------

The approximation has the form:

.. math::

       y_t = y^s + 0.5 \Delta^2 + A y^h_{t-1} + B u_t + 0.5 C (y^h_{t-1}\otimes y^h_{t-1}) + 0.5 D (u_t \otimes u_t) + E (y^h_{t-1} \otimes u_t)

where :math:`y^s` is the steady state value of :math:`y`,
:math:`y^h_t=y_t-y^s`, and :math:`\Delta^2` is the shift effect of the
variance of future shocks. For the reordering required due to
differences in declaration and DR order, see the first order
approximation.

The coefficients of the decision rules are stored in the variables
described for first order approximation, plus the following variables:

* :math:`\Delta^2` is stored in ``oo_.dr.ghs2``. The vector rows
  correspond to all endogenous in DR-order.
* :math:`C` is stored in ``oo_.dr.ghxx``. The matrix rows correspond
  to all endogenous in DR-order. The matrix columns correspond to the
  Kronecker product of the vector of state variables in DR-order.
* :math:`D` is stored in ``oo_.dr.ghuu``. The matrix rows correspond
  to all endogenous in DR-order. The matrix columns correspond to the
  Kronecker product of exogenous variables in declaration order.
* :math:`E` is stored in ``oo_.dr.ghxu``. The matrix rows correspond
  to all endogenous in DR-order. The matrix columns correspond to the
  Kronecker product of the vector of state variables (in DR-order) by
  the vector of exogenous variables (in declaration order).

Third-order approximation
-------------------------

The approximation has the form:

.. math::

       y_t = y^s + G_0 + G_1 z_t + G_2 (z_t \otimes z_t) + G_3 (z_t \otimes z_t \otimes z_t)

where :math:`y^s` is the steady state value of :math:`y`, and
:math:`z_t` is a vector consisting of the deviation from the steady
state of the state variables (in DR-order) at date :math:`t-1`
followed by the exogenous variables at date :math:`t` (in declaration
order). The vector :math:`z_t` is therefore of size :math:`n_z` =
``M_.nspred`` + ``M_.exo_nbr``.

The coefficients of the decision rules are stored as follows:

* :math:`y^s` is stored in ``oo_.dr.ys``. The vector rows correspond
  to all endogenous in the declaration order.
* :math:`G_0` is stored in ``oo_.dr.g_0``. The vector rows correspond
  to all endogenous in DR-order.
* :math:`G_1` is stored in ``oo_.dr.g_1``. The matrix rows correspond
  to all endogenous in DR-order. The matrix columns correspond to
  state variables in DR-order, followed by exogenous in declaration
  order.
* :math:`G_2` is stored in ``oo_.dr.g_2``. The matrix rows correspond
  to all endogenous in DR-order. The matrix columns correspond to the
  Kronecker product of state variables (in DR-order), followed by
  exogenous (in declaration order). Note that the Kronecker product is
  stored in a folded way, i.e. symmetric elements are stored only
  once, which implies that the matrix has :math:`n_z(n_z+1)/2`
  columns. More precisely, each column of this matrix corresponds to a
  pair :math:`(i_1, i_2)` where each index represents an element of
  :math:`z_t` and is therefore between :math:`1` and :math:`n_z`. Only
  non-decreasing pairs are stored, i.e. those for which :math:`i_1
  \leq i_2`. The columns are arranged in the lexicographical order of
  non-decreasing pairs. Also note that for those pairs where
  :math:`i_1 \neq i_2`, since the element is stored only once but
  appears two times in the unfolded :math:`G_2` matrix, it must be
  multiplied by 2 when computing the decision rules.
* :math:`G_3` is stored in ``oo_.dr.g_3``. The matrix rows correspond
  to all endogenous in DR-order. The matrix columns correspond to the
  third Kronecker power of state variables (in DR-order), followed by
  exogenous (in declaration order). Note that the third Kronecker
  power is stored in a folded way, i.e. symmetric elements are stored
  only once, which implies that the matrix has
  :math:`n_z(n_z+1)(n_z+2)/6` columns. More precisely, each column of
  this matrix corresponds to a tuple :math:`(i_1, i_2, i_3)` where
  each index represents an element of :math:`z_t` and is therefore
  between :math:`1` and :math:`n_z`. Only non-decreasing tuples are
  stored, i.e. those for which :math:`i_1 \leq i_2 \leq i_3`. The
  columns are arranged in the lexicographical order of non-decreasing
  tuples. Also note that for tuples that have three distinct indices
  (i.e. :math:`i_1 \neq i_2` and :math:`i_1 \neq i_3` and :math:`i_2
  \neq i_3`), since these elements are stored only once but appears
  six times in the unfolded :math:`G_3` matrix, they must be
  multiplied by 6 when computing the decision rules. Similarly, for
  those tuples that have two equal indices (i.e. of the form
  :math:`(a,a,b)` or :math:`(a,b,a)` or :math:`(b,a,a)`), since these
  elements are stored only once but appears three times in the
  unfolded :math:`G_3` matrix, they must be multiplied by 3 when
  computing the decision rules.


.. _estim:

Estimation
==========

Provided that you have observations on some endogenous variables, it
is possible to use Dynare to estimate some or all parameters. Both
maximum likelihood (as in *Ireland (2004)*) and Bayesian techniques
(as in *Rabanal and Rubio-Ramirez (2003)*, *Schorfheide (2000)* or
*Smets and Wouters (2003)*) are available. Using Bayesian methods, it
is possible to estimate DSGE models, VAR models, or a combination of
the two techniques called DSGE-VAR.

Note that in order to avoid stochastic singularity, you must have at
least as many shocks or measurement errors in your model as you have
observed variables.

The estimation using a first order approximation can benefit from the
block decomposition of the model (see :opt:`block`).

.. command:: varobs VARIABLE_NAME...;

    |br| This command lists the name of observed endogenous variables
    for the estimation procedure. These variables must be available in
    the data file (see :ref:`estimation_cmd <estim-comm>`).

    Alternatively, this command is also used in conjunction with the
    ``partial_information`` option of ``stoch_simul``, for declaring
    the set of observed variables when solving the model under partial
    information.

    Only one instance of ``varobs`` is allowed in a model file. If one
    needs to declare observed variables in a loop, the macro-processor
    can be used as shown in the second example below.

    *Example*

        ::

           varobs C y rr;

        Declares endogenous variables ``C``, ``y`` and ``rr`` as
        observed variables.

    *Example* (with a macro-processor loop)

        ::

           varobs
           @#for co in countries
           GDP_@{co}
           @#endfor
           ;


.. block:: observation_trends ;

    |br| This block specifies linear trends for observed variables as
    functions of model parameters. In case the ``loglinear`` option is
    used, this corresponds to a linear trend in the logged
    observables, i.e. an exponential trend in the level of the
    observables.

    Each line inside of the block should be of the form::

        VARIABLE_NAME(EXPRESSION);

    In most cases, variables shouldn’t be centered when
    ``observation_trends`` is used.

    *Example*

        ::

            observation_trends;
            Y (eta);
            P (mu/eta);
            end;


.. block:: estimated_params ;

    |br| This block lists all parameters to be estimated and specifies
    bounds and priors as necessary.

    Each line corresponds to an estimated parameter.

    In a maximum likelihood estimation, each line follows this syntax::

        stderr VARIABLE_NAME | corr VARIABLE_NAME_1, VARIABLE_NAME_2 | PARAMETER_NAME
        , INITIAL_VALUE [, LOWER_BOUND, UPPER_BOUND ];

    In a Bayesian estimation, each line follows this syntax::

        stderr VARIABLE_NAME | corr VARIABLE_NAME_1, VARIABLE_NAME_2 | PARAMETER_NAME | DSGE_PRIOR_WEIGHT
        [, INITIAL_VALUE [, LOWER_BOUND, UPPER_BOUND]], PRIOR_SHAPE,
        PRIOR_MEAN, PRIOR_STANDARD_ERROR [, PRIOR_3RD_PARAMETER [,
        PRIOR_4TH_PARAMETER [, SCALE_PARAMETER ] ] ];

    The first part of the line consists of one of the four following
    alternatives:

    * ``stderr VARIABLE_NAME``

      Indicates that the standard error of the exogenous
      variable VARIABLE_NAME, or of the observation
      error/measurement errors associated with endogenous
      observed variable VARIABLE_NAME, is to be estimated.

    * ``corr VARIABLE_NAME1, VARIABLE_NAME2``

      Indicates that the correlation between the exogenous
      variables VARIABLE_NAME1 and VARIABLE_NAME2, or the
      correlation of the observation errors/measurement errors
      associated with endogenous observed variables
      VARIABLE_NAME1 and VARIABLE_NAME2, is to be
      estimated. Note that correlations set by previous
      shocks-blocks or estimation-commands are kept at their
      value set prior to estimation if they are not estimated
      again subsequently. Thus, the treatment is the same as in
      the case of deep parameters set during model calibration
      and not estimated.

    * ``PARAMETER_NAME``

      The name of a model parameter to be estimated

    * ``DSGE_PRIOR_WEIGHT``

      Special name for the weigh of the DSGE model in DSGE-VAR model.


    The rest of the line consists of the following fields, some of
    them being optional:

    .. option:: INITIAL_VALUE

       Specifies a starting value for the posterior mode optimizer or
       the maximum likelihood estimation. If unset, defaults to the
       prior mean.

    .. option:: LOWER_BOUND

       Specifies a lower bound for the parameter value in maximum
       likelihood estimation. In a Bayesian estimation context, sets a
       lower bound only effective while maximizing the posterior
       kernel. This lower bound does not modify the shape of the prior
       density, and is only aimed at helping the optimizer in
       identifying the posterior mode (no consequences for the
       MCMC). For some prior densities (namely inverse gamma, gamma,
       uniform, beta or Weibull) it is possible to shift the support
       of the prior distributions to the left or the right using
       :opt:`prior_3rd_parameter <PRIOR_3RD_PARAMETER>`. In this case
       the prior density is effectively modified (note that the
       truncated Gaussian density is not implemented in Dynare). If
       unset, defaults to minus infinity (ML) or the natural lower
       bound of the prior (Bayesian estimation).

    .. option:: UPPER_BOUND

       Same as ``lower_bound``, but specifying an upper bound instead.

    .. option:: PRIOR_SHAPE

       A keyword specifying the shape of the prior density. The
       possible values are: ``beta_pdf``, ``gamma_pdf``,
       ``normal_pdf``, ``uniform_pdf``, ``inv_gamma_pdf``,
       ``inv_gamma1_pdf``, ``inv_gamma2_pdf`` and
       ``weibull_pdf``. Note that ``inv_gamma_pdf`` is equivalent to
       ``inv_gamma1_pdf``.

    .. option:: PRIOR_MEAN

       The mean of the prior distribution.

    .. option:: PRIOR_STANDARD_ERROR

       The standard error of the prior distribution.

    .. option:: PRIOR_3RD_PARAMETER

       A third parameter of the prior used for generalized beta
       distribution, generalized gamma, generalized Weibull and for
       the uniform distribution. Default: ``0``.

    .. option:: PRIOR_4TH_PARAMETER

        A fourth parameter of the prior used for generalized beta
        distribution and for the uniform distribution. Default: ``1``.

    .. option:: SCALE_PARAMETER

        A parameter specific scale parameter for the jumping
        distribution’s covariance matrix of the Metropolis-Hasting
        algorithm.

    Note that INITIAL_VALUE, LOWER_BOUND, UPPER_BOUND, PRIOR_MEAN,
    PRIOR_STANDARD_ERROR, PRIOR_3RD_PARAMETER, PRIOR_4TH_PARAMETER and
    SCALE_PARAMETER can be any valid EXPRESSION. Some of them can be
    empty, in which Dynare will select a default value depending on
    the context and the prior shape.

    As one uses options more towards the end of the list, all previous
    options must be filled: for example, if you want to specify
    SCALE_PARAMETER, you must specify ``PRIOR_3RD_PARAMETER`` and
    ``PRIOR_4TH_PARAMETER``. Use empty values, if these parameters
    don’t apply.

    *Example*

       ::

          corr eps_1, eps_2, 0.5,  ,  , beta_pdf, 0, 0.3, -1, 1;

       Sets a generalized beta prior for the correlation between
       ``eps_1`` and ``eps_2`` with mean ``0`` and variance
       ``0.3``. By setting ``PRIOR_3RD_PARAMETER`` to ``-1`` and
       ``PRIOR_4TH_PARAMETER`` to ``1`` the standard beta distribution
       with support ``[0,1]`` is changed to a generalized beta with
       support ``[-1,1]``. Note that LOWER_BOUND and UPPER_BOUND are
       left empty and thus default to ``-1`` and ``1``,
       respectively. The initial value is set to ``0.5``.

    *Example*

       ::

          corr eps_1, eps_2, 0.5,  -0.5,  1, beta_pdf, 0, 0.3, -1, 1;

       Sets the same generalized beta distribution as before, but now
       truncates this distribution to ``[-0.5,1]`` through the use of
       LOWER_BOUND and UPPER_BOUND. Hence, the prior does not
       integrate to ``1`` anymore.

    *Parameter transformation*

    Sometimes, it is desirable to estimate a transformation of a
    parameter appearing in the model, rather than the parameter
    itself. It is of course possible to replace the original parameter
    by a function of the estimated parameter everywhere is the model,
    but it is often unpractical.

    In such a case, it is possible to declare the parameter to be
    estimated in the parameters statement and to define the
    transformation, using a pound sign (#) expression (see
    :ref:`model-decl`).

    *Example*

        ::

            parameters bet;

            model;
            # sig = 1/bet;
            c = sig*c(+1)*mpk;
            end;

            estimated_params;
            bet, normal_pdf, 1, 0.05;
            end;


.. block:: estimated_params_init ;
           estimated_params_init (OPTIONS...);

    |br| This block declares numerical initial values for the
    optimizer when these ones are different from the prior mean. It
    should be specified after the ``estimated_params`` block as
    otherwise the specified starting values are overwritten by the
    latter.

    Each line has the following syntax::

        stderr VARIABLE_NAME | corr VARIABLE_NAME_1, VARIABLE_NAME_2 | PARAMETER_NAME, INITIAL_VALUE;

    *Options*

    .. option:: use_calibration

        For not specifically initialized parameters, use the deep
        parameters and the elements of the covariance matrix specified
        in the ``shocks`` block from calibration as starting values
        for estimation. For components of the ``shocks`` block that
        were not explicitly specified during calibration or which
        violate the prior, the prior mean is used.

    See :bck:`estimated_params`, for the meaning and syntax of the
    various components.


.. block:: estimated_params_bounds ;

    |br| This block declares lower and upper bounds for parameters in maximum likelihood estimation.

    Each line has the following syntax::

        stderr VARIABLE_NAME | corr VARIABLE_NAME_1, VARIABLE_NAME_2 | PARAMETER_NAME, LOWER_BOUND, UPPER_BOUND;

    See :bck:`estimated_params`, for the meaning and syntax of the
    various components.


.. _estim-comm:

.. command:: estimation [VARIABLE_NAME...];
             estimation (OPTIONS...) [VARIABLE_NAME...];

    |br| This command runs Bayesian or maximum likelihood estimation.

    The following information will be displayed by the command:

    * Results from posterior optimization (also for maximum likelihood)
    * Marginal log data density
    * Posterior mean and highest posterior density interval (shortest
      credible set) from posterior simulation
    * Convergence diagnostic table when only one MCM chain is used or
      Metropolis-Hastings convergence graphs documented in *Pfeiffer
      (2014)* in case of multiple MCM chains
    * Table with numerical inefficiency factors of the MCMC
    * Graphs with prior, posterior, and mode
    * Graphs of smoothed shocks, smoothed observation errors, smoothed
      and historical variables

    Note that the posterior moments, smoothed variables, k-step ahead
    filtered variables and forecasts (when requested) will only be
    computed on the variables listed after the ``estimation``
    command. Alternatively, one can choose to compute these quantities
    on all endogenous or on all observed variables (see
    ``consider_all_endogenous`` and ``consider_only_observed`` options
    below). If no variable is listed after the estimation command,
    then Dynare will interactively ask which variable set to use.

    Also, during the MCMC (Bayesian estimation with ``mh_replic``
    :math:`>0`) a (graphical or text) waiting bar is displayed showing
    the progress of the Monte-Carlo and the current value of the
    acceptance ratio. Note that if the ``load_mh_file`` option is used
    (see below) the reported acceptance ratio does not take into
    account the draws from the previous MCMC. In the literature there
    is a general agreement for saying that the acceptance ratio should
    be close to one third or one quarter. If this not the case, you
    can stop the MCMC (``Ctrl-C``) and change the value of option
    ``mh_jscale`` (see below).

    Note that by default Dynare generates random numbers using the
    algorithm ``mt199937ar`` (i.e. Mersenne Twister method) with a
    seed set equal to ``0``. Consequently the MCMCs in Dynare are
    deterministic: one will get exactly the same results across
    different Dynare runs (*ceteris paribus*). For instance, the
    posterior moments or posterior densities will be exactly the
    same. This behaviour allows to easily identify the consequences of
    a change on the model, the priors or the estimation options. But
    one may also want to check that across multiple runs, with
    different sequences of proposals, the returned results are almost
    identical. This should be true if the number of iterations
    (i.e. the value of ``mh_replic``) is important enough to ensure
    the convergence of the MCMC to its ergodic distribution. In this
    case the default behaviour of the random number generators in not
    wanted, and the user should set the seed according to the system
    clock before the estimation command using the following command::

        set_dynare_seed('clock');

    so that the sequence of proposals will be different across different runs.

    *Algorithms*

    The Monte Carlo Markov Chain (MCMC) diagnostics are generated by
    the estimation command if :opt:`mh_replic <mh_replic = INTEGER>`
    is larger than 2000 and if option :opt:`nodiagnostic` is not
    used. If :opt:`mh_nblocks <mh_nblocks = INTEGER>` is equal to one,
    the convergence diagnostics of *Geweke (1992,1999)* is
    computed. It uses a chi-square test to compare the means of the
    first and last draws specified by :opt:`geweke_interval
    <geweke_interval = [DOUBLE DOUBLE]>` after discarding the burn-in
    of :opt:`mh_drop <mh_drop = DOUBLE>`. The test is computed using
    variance estimates under the assumption of no serial correlation
    as well as using tapering windows specified in :opt:`taper_steps
    <taper_steps = [INTEGER1 INTEGER2 ...]>`. If :opt:`mh_nblocks
    <mh_nblocks = INTEGER>` is larger than 1, the convergence
    diagnostics of *Brooks and Gelman (1998)* are used instead. As
    described in section 3 of *Brooks and Gelman (1998)* the
    univariate convergence diagnostics are based on comparing pooled
    and within MCMC moments (Dynare displays the second and third
    order moments, and the length of the Highest Probability Density
    interval covering 80% of the posterior distribution). Due to
    computational reasons, the multivariate convergence diagnostic
    does not follow *Brooks and Gelman (1998)* strictly, but rather
    applies their idea for univariate convergence diagnostics to the
    range of the posterior likelihood function instead of the
    individual parameters. The posterior kernel is used to aggregate
    the parameters into a scalar statistic whose convergence is then
    checked using the *Brooks and Gelman (1998)* univariate
    convergence diagnostic.

    The inefficiency factors are computed as in *Giordano et
    al.(2011)* based on Parzen windows as in e.g. *Andrews (1991)*.


    *Options*

    .. _dataf:

    .. option:: datafile = FILENAME

       The datafile: a ``.m`` file, a ``.mat`` file, a ``.csv`` file,
       or a ``.xls/.xlsx`` file (under Octave, the `io
       <http://octave.sourceforge.net/io/>`_ package from Octave-Forge
       is required for the ``.csv`` and ``.xlsx`` formats and the
       ``.xls`` file extension is not supported). Note that the base
       name (i.e. without extension) of the datafile has to be
       different from the base name of the model file. If there are
       several files named FILENAME, but with different file endings,
       the file name must be included in quoted strings and provide
       the file ending like::

            estimation(datafile='../fsdat_simul.mat',...);

    .. option:: dirname = FILENAME

       Directory in which to store ``estimation`` output. To pass a
       subdirectory of a directory, you must quote the
       argument. Default: ``<mod_file>``.

    .. option:: xls_sheet = NAME

       The name of the sheet with the data in an Excel file.

    .. option:: xls_range = RANGE

       The range with the data in an Excel file. For example,
       ``xls_range=B2:D200``.

    .. option:: nobs = INTEGER

       The number of observations following :opt:`first_obs <first_obs
       = [INTEGER1:INTEGER2]>` to be used. Default: all observations
       in the file after ``first_obs``.

    .. option:: nobs = [INTEGER1:INTEGER2]

       Runs a recursive estimation and forecast for samples of size
       ranging of ``INTEGER1`` to ``INTEGER2``. Option ``forecast``
       must also be specified. The forecasts are stored in the
       ``RecursiveForecast`` field of the results structure (see
       :mvar:`RecursiveForecast <oo_.RecursiveForecast>`). The
       respective results structures ``oo_`` are saved in
       ``oo_recursive_`` (see :mvar:`oo_recursive_`) and are indexed
       with the respective sample length.

    .. option:: first_obs = INTEGER

       The number of the first observation to be used. In case of
       estimating a DSGE-VAR, ``first_obs`` needs to be larger than
       the number of lags. Default: ``1``.

    .. option:: first_obs = [INTEGER1:INTEGER2]

       Runs a rolling window estimation and forecast for samples of
       fixed size ``nobs`` starting with the first observation ranging
       from ``INTEGER1`` to ``INTEGER2``. Option ``forecast`` must
       also be specified. This option is incompatible with requesting
       recursive forecasts using an expanding window (see :opt:`nobs
       <nobs = [INTEGER1:INTEGER2]>`). The respective results
       structures ``oo_`` are saved in ``oo_recursive_`` (see
       :mvar:`oo_recursive_`) and are indexed with the respective
       first observation of the rolling window.

    .. option:: prefilter = INTEGER

       A value of 1 means that the estimation procedure will demean
       each data series by its empirical mean. If the :ref:`loglinear
       <logl>` option without the :opt:`logdata` option is requested,
       the data will first be logged and then demeaned. Default:
       ``0``, i.e. no prefiltering.

    .. option:: presample = INTEGER

       The number of observations after :opt:`first_obs <first_obs =
       [INTEGER1:INTEGER2]>` to be skipped before evaluating the
       likelihood. These presample observations do not enter the
       likelihood, but are used as a training sample for starting the
       Kalman filter iterations. This option is incompatible with
       estimating a DSGE-VAR. Default: ``0``.

    .. _logl:

    .. option:: loglinear

       Computes a log-linear approximation of the model instead of a
       linear approximation. As always in the context of estimation,
       the data must correspond to the definition of the variables
       used in the model (see *Pfeifer (2013)* for more details on how
       to correctly specify observation equations linking model
       variables and the data). If you specify the loglinear option,
       Dynare will take the logarithm of both your model variables and
       of your data as it assumes the data to correspond to the
       original non-logged model variables. The displayed posterior
       results like impulse responses, smoothed variables, and moments
       will be for the logged variables, not the original un-logged
       ones. Default: computes a linear approximation.

    .. option:: logdata

       Dynare applies the :math:`log` transformation to the provided
       data if a log-linearization of the model is requested
       (:opt:`loglinear`) unless ``logdata`` option is used. This
       option is necessary if the user provides data already in logs,
       otherwise the :math:`log` transformation will be applied twice
       (this may result in complex data).

    .. option:: plot_priors = INTEGER

       Control the plotting of priors.

           ``0``

                No prior plot.

           ``1``

               Prior density for each estimated parameter is
               plotted. It is important to check that the actual shape
               of prior densities matches what you have in
               mind. Ill-chosen values for the prior standard density
               can result in absurd prior densities.

       |br| Default value is ``1``.

    .. option:: nograph

       See :opt:`nograph`.

    .. option:: posterior_nograph

       Suppresses the generation of graphs associated with Bayesian
       IRFs (:opt:`bayesian_irf`), posterior smoothed objects
       (:opt:`smoother`), and posterior forecasts (:opt:`forecast`).

    .. option:: posterior_graph

       Re-enables the generation of graphs previously shut off with
       :opt:`posterior_nograph`.

    .. option:: nodisplay

       See :opt:`nodisplay`.

    .. option:: graph_format = FORMAT
                graph_format = ( FORMAT, FORMAT... )

       See :opt:`graph_format <graph_format = ( FORMAT, FORMAT... )>`.

    .. option:: lik_init = INTEGER

       Type of initialization of Kalman filter:

           ``1``

               For stationary models, the initial matrix of variance
               of the error of forecast is set equal to the
               unconditional variance of the state variables.

           ``2``

               For nonstationary models: a wide prior is used with an
               initial matrix of variance of the error of forecast
               diagonal with 10 on the diagonal (follows the
               suggestion of *Harvey and Phillips(1979)*).

           ``3``

               For nonstationary models: use a diffuse filter (use
               rather the ``diffuse_filter`` option).

           ``4``

               The filter is initialized with the fixed point of the
               Riccati equation.

           ``5``

               Use i) option 2 for the non-stationary elements by
               setting their initial variance in the forecast error
               matrix to 10 on the diagonal and all covariances to 0
               and ii) option 1 for the stationary elements.

       |br| Default value is 1. For advanced use only.

    .. option:: lik_algo = INTEGER

        For internal use and testing only.

    .. option:: conf_sig = DOUBLE

       Confidence interval used for classical forecasting after
       estimation. See :ref:`conf_sig <confsig>`.

    .. option:: mh_conf_sig = DOUBLE

       Confidence/HPD interval used for the computation of prior and
       posterior statistics like: parameter distributions,
       prior/posterior moments, conditional variance decomposition,
       impulse response functions, Bayesian forecasting. Default:
       ``0.9``.

    .. option:: mh_replic = INTEGER

       Number of replications for Metropolis-Hastings algorithm. For
       the time being, ``mh_replic`` should be larger
       than 1200. Default: ``20000``.

    .. option:: sub_draws = INTEGER

       Number of draws from the MCMC that are used to compute
       posterior distribution of various objects (smoothed variable,
       smoothed shocks, forecast, moments, IRF). The draws used to
       compute these posterior moments are sampled uniformly in the
       estimated empirical posterior distribution (i.e. draws of the
       MCMC). ``sub_draws`` should be smaller than the total number of
       MCMC draws available. Default:
       ``min(posterior_max_subsample_draws, (Total number of
       draws)*(number of chains) )``.

    .. option:: posterior_max_subsample_draws = INTEGER

       Maximum number of draws from the MCMC used to compute posterior
       distribution of various objects (smoothed variable, smoothed
       shocks, forecast, moments, IRF), if not overriden by option
       ``sub_draws``. Default: ``1200``.

    .. option:: mh_nblocks = INTEGER

       Number of parallel chains for Metropolis-Hastings
       algorithm. Default: ``2``.

    .. option:: mh_drop = DOUBLE

       The fraction of initially generated parameter vectors to be
       dropped as a burn-in before using posterior
       simulations. Default: ``0.5``.

    .. option:: mh_jscale = DOUBLE

       The scale parameter of the jumping distribution’s covariance
       matrix (Metropolis-Hastings or TaRB-algorithm). The default
       value is rarely satisfactory. This option must be tuned to
       obtain, ideally, an acceptance ratio of 25%-33%. Basically, the
       idea is to increase the variance of the jumping distribution if
       the acceptance ratio is too high, and decrease the same
       variance if the acceptance ratio is too low. In some situations
       it may help to consider parameter-specific values for this
       scale parameter. This can be done in the
       :bck:`estimated_params` block.

       Note that ``mode_compute=6`` will tune the scale parameter to
       achieve an acceptance rate of
       :ref:`AcceptanceRateTarget<art>`. The resulting scale parameter
       will be saved into a file named
       ``MODEL_FILENAME_mh_scale.mat.`` This file can be loaded in
       subsequent runs via the ``posterior_sampler_options`` option
       :ref:`scale_file <scale-file>`. Both ``mode_compute=6`` and
       ``scale_file`` will overwrite any value specified in
       ``estimated_params`` with the tuned value. Default: ``0.2``.

    .. option:: mh_init_scale = DOUBLE

       The scale to be used for drawing the initial value of the
       Metropolis-Hastings chain. Generally, the starting points
       should be overdispersed for the *Brooks and Gelman (1998)*
       convergence diagnostics to be meaningful. Default:
       ``2*mh_jscale.``

       It is important to keep in mind that ``mh_init_scale`` is set
       at the beginning of Dynare execution, i.e. the default will not
       take into account potential changes in ``mh_jscale`` introduced
       by either ``mode_compute=6`` or the
       ``posterior_sampler_options`` option :ref:`scale_file
       <scale-file>`. If ``mh_init_scale`` is too wide during
       initalization of the posterior sampler so that 100 tested draws
       are inadmissible (e.g. Blanchard-Kahn conditions are always
       violated), Dynare will request user input of a new
       ``mh_init_scale`` value with which the next 100 draws will be
       drawn and tested. If the :opt:`nointeractive` option has been
       invoked, the program will instead automatically decrease
       ``mh_init_scale`` by 10 percent after 100 futile draws and try
       another 100 draws. This iterative procedure will take place at
       most 10 times, at which point Dynare will abort with an error
       message.

    .. option:: mh_recover

       Attempts to recover a Metropolis-Hastings simulation that
       crashed prematurely, starting with the last available saved
       ``mh``-file. Shouldn’t be used together with ``load_mh_file``
       or a different ``mh_replic`` than in the crashed run. Since
       Dynare 4.5 the proposal density from the previous run will
       automatically be loaded. In older versions, to assure a neat
       continuation of the chain with the same proposal density, you
       should provide the ``mode_file`` used in the previous run or
       the same user-defined ``mcmc_jumping_covariance`` when using
       this option. Note that under Octave, a neat continuation of the
       crashed chain with the respective last random number generator
       state is currently not supported.

    .. option:: mh_mode = INTEGER

        ...

    .. option:: mode_file = FILENAME

       Name of the file containing previous value for the mode. When
       computing the mode, Dynare stores the mode (``xparam1``) and
       the hessian (``hh``, only if ``cova_compute=1``) in a file
       called ``MODEL_FILENAME_mode.mat``. After a successful run of
       the estimation command, the ``mode_file`` will be disabled to
       prevent other function calls from implicitly using an updated
       ``mode-file``. Thus, if the mod-file contains subsequent
       ``estimation`` commands, the ``mode_file`` option, if desired,
       needs to be specified again.

    .. option:: mode_compute = INTEGER | FUNCTION_NAME

       Specifies the optimizer for the mode computation:

           ``0``

                The mode isn’t computed. When the ``mode_file`` option
                is specified, the mode is simply read from that file.

                When ``mode_file`` option is not specified, Dynare
                reports the value of the log posterior (log
                likelihood) evaluated at the initial value of the
                parameters.

                When ``mode_file`` is not specified and there is no
                ``estimated_params`` block, but the ``smoother``
                option is used, it is a roundabout way to compute the
                smoothed value of the variables of a model with
                calibrated parameters.

           ``1``

                Uses ``fmincon`` optimization routine (available under
                MATLAB if the Optimization Toolbox is installed; not
                available under Octave).

           ``2``

                Uses the continuous simulated annealing global
                optimization algorithm described in *Corana et
                al.(1987)* and *Goffe et al.(1994)*.

           ``3``

                Uses ``fminunc`` optimization routine (available under
                MATLAB if the Optimization Toolbox is installed;
                available under Octave if the `optim
                <http://octave.sourceforge.net/optim/>`_ package from
                Octave-Forge is installed).

           ``4``

                Uses Chris Sims’s ``csminwel``.

           ``5``

                Uses Marco Ratto’s ``newrat``. This value is not
                compatible with non linear filters or DSGE-VAR
                models. This is a slice optimizer: most iterations are
                a sequence of univariate optimization step, one for
                each estimated parameter or shock. Uses ``csminwel``
                for line search in each step.

           ``6``

                Uses a Monte-Carlo based optimization routine (see
                `Dynare wiki`_ for more details).

           ``7``

                Uses ``fminsearch``, a simplex-based optimization
                routine (available under MATLAB if the Optimization
                Toolbox is installed; available under Octave if the
                optim package from Octave-Forge is installed).

           ``8``

                Uses Dynare implementation of the Nelder-Mead
                simplex-based optimization routine (generally more
                efficient than the MATLAB or Octave implementation
                available with ``mode_compute=7``).

           ``9``

                Uses the CMA-ES (Covariance Matrix Adaptation
                Evolution Strategy) algorithm of *Hansen and Kern
                (2004)*, an evolutionary algorithm for difficult
                non-linear non-convex optimization.

           ``10``

                Uses the ``simpsa`` algorithm, based on the
                combination of the non-linear simplex and simulated
                annealing algorithms as proposed by *Cardoso, Salcedo
                and Feyo de Azevedo (1996)*.

           ``11``

                This is not strictly speaking an optimization
                algorithm. The (estimated) parameters are treated as
                state variables and estimated jointly with the
                original state variables of the model using a
                nonlinear filter. The algorithm implemented in Dynare
                is described in *Liu and West (2001)*.

           ``12``

                Uses the ``particleswarm`` optimization routine
                (available under MATLAB if the Global Optimization
                Toolbox is installed; not available under Octave).

           ``101``

                Uses the SolveOpt algorithm for local nonlinear
                optimization problems proposed by *Kuntsevich and
                Kappel (1997)*.

           ``102``

                Uses ``simulannealbnd`` optimization routine
                (available under MATLAB if the Global Optimization
                Toolbox is installed; not available under Octave)

           ``FUNCTION_NAME``

                It is also possible to give a FUNCTION_NAME to this
                option, instead of an INTEGER. In that case, Dynare
                takes the return value of that function as the
                posterior mode.

       |br| Default value is ``4``.

    .. option:: silent_optimizer

       Instructs Dynare to run mode computing/optimization silently
       without displaying results or saving files in between. Useful
       when running loops.

    .. option:: mcmc_jumping_covariance = OPTION

       Tells Dynare which covariance to use for the proposal density
       of the MCMC sampler. OPTION can be one of the following:

           ``hessian``

               Uses the Hessian matrix computed at the mode.

           ``prior_variance``

               Uses the prior variances. No infinite prior variances
               are allowed in this case.

           ``identity_matrix``

               Uses an identity matrix.

           ``FILENAME``

               Loads an arbitrary user-specified covariance matrix
               from ``FILENAME.mat``. The covariance matrix must be
               saved in a variable named ``jumping_covariance``, must
               be square, positive definite, and have the same
               dimension as the number of estimated parameters.

       Note that the covariance matrices are still scaled with
       :opt:`mh_jscale <mh_jscale = DOUBLE>`. Default value is
       ``hessian``.

    .. option:: mode_check

       Tells Dynare to plot the posterior density for values around
       the computed mode for each estimated parameter in turn. This is
       helpful to diagnose problems with the optimizer. Note that for
       ``order>1`` the likelihood function resulting from the particle
       filter is not differentiable anymore due to random chatter
       introduced by selecting different particles for different
       parameter values. For this reason, the ``mode_check`` plot may
       look wiggly.

    .. option:: mode_check_neighbourhood_size = DOUBLE

       Used in conjunction with option ``mode_check``, gives the width
       of the window around the posterior mode to be displayed on the
       diagnostic plots. This width is expressed in percentage
       deviation. The ``Inf`` value is allowed, and will trigger a
       plot over the entire domain (see also
       ``mode_check_symmetric_plots``). Default:``0.5``.

    .. option:: mode_check_symmetric_plots = INTEGER

       Used in conjunction with option ``mode_check``, if set to
       ``1``, tells Dynare to ensure that the check plots are
       symmetric around the posterior mode. A value of ``0`` allows to
       have asymmetric plots, which can be useful if the posterior
       mode is close to a domain boundary, or in conjunction with
       ``mode_check_neighbourhood_size = Inf`` when the domain in not
       the entire real line. Default: ``1``.

    .. option:: mode_check_number_of_points = INTEGER

       Number of points around the posterior mode where the posterior
       kernel is evaluated (for each parameter). Default is ``20``.

    .. option:: prior_trunc = DOUBLE

       Probability of extreme values of the prior density that is
       ignored when computing bounds for the parameters. Default:
       ``1e-32``.

    .. option:: huge_number = DOUBLE

       Value for replacing infinite values in the definition of
       (prior) bounds when finite values are required for
       computational reasons. Default: ``1e7``.

    .. option:: load_mh_file

       Tells Dynare to add to previous Metropolis-Hastings simulations
       instead of starting from scratch. Since Dynare 4.5 the proposal
       density from the previous run will automatically be loaded. In
       older versions, to assure a neat continuation of the chain with
       the same proposal density, you should provide the ``mode_file``
       used in the previous run or the same user-defined
       ``mcmc_jumping_covariance`` when using this option. Shouldn’t
       be used together with ``mh_recover``. Note that under Octave, a
       neat continuation of the chain with the last random number
       generator state of the already present draws is currently not
       supported.

    .. option:: load_results_after_load_mh

       This option is available when loading a previous MCMC run
       without adding additional draws, i.e. when ``load_mh_file`` is
       specified with ``mh_replic=0``. It tells Dynare to load the
       previously computed convergence diagnostics, marginal data
       density, and posterior statistics from an existing ``_results``
       file instead of recomputing them.

    .. option:: optim = (NAME, VALUE, ...)

       A list of NAME and VALUE pairs. Can be used to set options for
       the optimization routines. The set of available options depends
       on the selected optimization routine (i.e. on the value of
       option :opt:`mode_compute <mode_compute = INTEGER |
       FUNCTION_NAME>`):

           ``1, 3, 7, 12``

               Available options are given in the documentation of the
               MATLAB Optimization Toolbox or in Octave’s
               documentation.

           ``2``

               Available options are:

                   ``'initial_step_length'``

                       Initial step length. Default: ``1``.

                   ``'initial_temperature'``

                       Initial temperature. Default: ``15``.

                   ``'MaxIter'``

                       Maximum number of function
                       evaluations. Default: ``100000``.

                   ``'neps'``

                       Number of final function values used to decide
                       upon termination. Default: ``10``.

                   ``'ns'``

                       Number of cycles. Default: ``10``.

                   ``'nt'``

                       Number of iterations before temperature
                       reduction. Default: ``10``.

                   ``'step_length_c'``

                       Step length adjustment. Default: ``0.1``.

                   ``'TolFun'``

                       Stopping criteria. Default: ``1e-8``.

                   ``'rt'``

                       Temperature reduction factor. Default: ``0.1``.

                   ``'verbosity'``

                       Controls verbosity of display during
                       optimization, ranging from ``0`` (silent) to
                       ``3`` (each function evaluation). Default:
                       ``1``

           ``4``

               Available options are:

                   ``'InitialInverseHessian'``

                       Initial approximation for the inverse of the
                       Hessian matrix of the posterior kernel (or
                       likelihood). Obviously this approximation has
                       to be a square, positive definite and symmetric
                       matrix. Default: ``'1e-4*eye(nx)'``, where nx
                       is the number of parameters to be estimated.

                   ``'MaxIter'``

                       Maximum number of iterations. Default: ``1000``.

                   ``'NumgradAlgorithm'``

                       Possible values are ``2``, ``3`` and ``5``,
                       respectively, corresponding to the two, three
                       and five points formula used to compute the
                       gradient of the objective function (see
                       *Abramowitz and Stegun (1964)*). Values ``13``
                       and ``15`` are more experimental. If
                       perturbations on the right and the left
                       increase the value of the objective function
                       (we minimize this function) then we force the
                       corresponding element of the gradient to be
                       zero. The idea is to temporarily reduce the
                       size of the optimization problem. Default:
                       ``2``.

                   ``'NumgradEpsilon'``

                       Size of the perturbation used to compute
                       numerically the gradient of the objective
                       function. Default: ``1e-6``.

                   ``'TolFun'``

                       Stopping criteria. Default: ``1e-7``.

                   ``'verbosity'``

                       Controls verbosity of display during
                       optimization. Set to ``0`` to set to
                       silent. Default: ``1``.

                   ``'SaveFiles'``

                       Controls saving of intermediate results during
                       optimization. Set to ``0`` to shut off
                       saving. Default: ``1``.

           ``5``

               Available options are:

               ``'Hessian'``

                   Triggers three types of Hessian
                   computations. ``0``: outer product gradient; ``1``:
                   default DYNARE Hessian routine; ``2``: ’mixed’
                   outer product gradient, where diagonal elements are
                   obtained using second order derivation formula and
                   outer product is used for correlation
                   structure. Both {0} and {2} options require
                   univariate filters, to ensure using maximum number
                   of individual densities and a positive definite
                   Hessian. Both {0} and {2} are quicker than default
                   DYNARE numeric Hessian, but provide decent starting
                   values for Metropolis for large models (option {2}
                   being more accurate than {0}). Default: ``1``.

               ``'MaxIter'``

                   Maximum number of iterations. Default: ``1000``.

               ``'TolFun'``

                   Stopping criteria. Default: ``1e-5`` for numerical
                   derivatives, ``1e-7`` for analytic derivatives.

               ``'verbosity'``

                   Controls verbosity of display during
                   optimization. Set to ``0`` to set to
                   silent. Default: ``1``.

               ``'SaveFiles'``

                   Controls saving of intermediate results during
                   optimization. Set to ``0`` to shut off
                   saving. Default: ``1``.

           ``6``

               Available options are:

                   .. _art:

                   ``'AcceptanceRateTarget'``

                       A real number between zero and one. The scale
                       parameter of the jumping distribution is
                       adjusted so that the effective acceptance rate
                       matches the value of option
                       ``'AcceptanceRateTarget'``. Default:
                       ``1.0/3.0``.

                   ``'InitialCovarianceMatrix'``

                       Initial covariance matrix of the jumping
                       distribution. Default is ``'previous'`` if
                       option ``mode_file`` is used, ``'prior'``
                       otherwise.

                   ``'nclimb-mh'``

                       Number of iterations in the last MCMC (climbing
                       mode). Default: ``200000``.

                   ``'ncov-mh'``

                       Number of iterations used for updating the
                       covariance matrix of the jumping
                       distribution. Default: ``20000``.

                   ``'nscale-mh'``

                       Maximum number of iterations used for adjusting
                       the scale parameter of the jumping
                       distribution. Default: ``200000``.

                   ``'NumberOfMh'``

                       Number of MCMC run sequentially. Default: ``3``.

           ``8``

               Available options are:

                   ``'InitialSimplexSize'``

                       Initial size of the simplex, expressed as
                       percentage deviation from the provided initial
                       guess in each direction. Default: ``.05``.

                   ``'MaxIter'``

                       Maximum number of iterations. Default: ``5000``.

                   ``'MaxFunEvals'``

                       Maximum number of objective function
                       evaluations. No default.

                   ``'MaxFunvEvalFactor'``

                       Set ``MaxFunvEvals`` equal to
                       ``MaxFunvEvalFactor`` times the number of
                       estimated parameters. Default: ``500``.

                   ``'TolFun'``

                       Tolerance parameter (w.r.t the objective
                       function). Default: ``1e-4``.

                   ``'TolX'``

                       Tolerance parameter (w.r.t the
                       instruments). Default: ``1e-4``.

                   ``'verbosity'``

                       Controls verbosity of display during
                       optimization. Set to ``0`` to set to
                       silent. Default: ``1``.

           ``9``

               Available options are:

                   ``'CMAESResume'``

                       Resume previous run. Requires the
                       ``variablescmaes.mat`` from the last run. Set
                       to ``1`` to enable. Default: ``0``.

                   ``'MaxIter'``

                       Maximum number of iterations.

                   ``'MaxFunEvals'``

                       Maximum number of objective function
                       evaluations. Default: ``Inf``.

                   ``'TolFun'``

                       Tolerance parameter (w.r.t the objective
                       function). Default: ``1e-7``.

                   ``'TolX'``

                       Tolerance parameter (w.r.t the
                       instruments). Default: ``1e-7``.

                   ``'verbosity'``

                       Controls verbosity of display during
                       optimization. Set to ``0`` to set to
                       silent. Default: ``1``.

                   ``'SaveFiles'``

                       Controls saving of intermediate results during
                       optimization. Set to ``0`` to shut off
                       saving. Default: ``1``.


           ``10``

               Available options are:

                   ``'EndTemperature'``

                       Terminal condition w.r.t the temperature. When
                       the temperature reaches ``EndTemperature``, the
                       temperature is set to zero and the algorithm
                       falls back into a standard simplex
                       algorithm. Default: ``0.1``.

                   ``'MaxIter'``

                       Maximum number of iterations. Default:
                       ``5000``.

                   ``'MaxFunvEvals'``

                       Maximum number of objective function
                       evaluations. No default.

                   ``'TolFun'``

                       Tolerance parameter (w.r.t the objective
                       function). Default: ``1e-4``.

                   ``'TolX'``

                       Tolerance parameter (w.r.t the
                       instruments). Default: ``1e-4``.

                   ``'verbosity'``

                       Controls verbosity of display during
                       optimization. Set to ``0`` to set to
                       silent. Default: ``1``.

           ``101``

               Available options are:

                   ``'LBGradientStep'``

                       Lower bound for the stepsize used for the
                       difference approximation of gradients. Default:
                       ``1e-11``.

                   ``'MaxIter'``

                       Maximum number of iterations. Default: ``15000``

                   ``'SpaceDilation'``

                       Coefficient of space dilation. Default: ``2.5``.

                   ``'TolFun'``

                       Tolerance parameter (w.r.t the objective
                       function). Default: ``1e-6``.

                   ``'TolX'``

                       Tolerance parameter (w.r.t the
                       instruments). Default: ``1e-6``.

                   ``'verbosity'``

                       Controls verbosity of display during
                       optimization. Set to ``0`` to set to
                       silent. Default: ``1``.

           ``102``

               Available options are given in the documentation of the
               MATLAB Global Optimization Toolbox.

       *Example*

           To change the defaults of ``csminwel`` (``mode_compute=4``)::

               estimation(..., mode_compute=4,optim=('NumgradAlgorithm',3,'TolFun',1e-5),...);


    .. option:: nodiagnostic

       Does not compute the convergence diagnostics for
       Metropolis-Hastings. Default: diagnostics are computed and
       displayed.

    .. option:: bayesian_irf

       Triggers the computation of the posterior distribution of
       IRFs. The length of the IRFs are controlled by the ``irf``
       option. Results are stored in ``oo_.PosteriorIRF.dsge`` (see
       below for a description of this variable).

    .. option:: relative_irf

        See :opt:`relative_irf`.

    .. option:: dsge_var = DOUBLE

       Triggers the estimation of a DSGE-VAR model, where the weight
       of the DSGE prior of the VAR model is calibrated to the value
       passed (see *Del Negro and Schorfheide (2004)*). It represents
       the ratio of dummy over actual observations. To assure that the
       prior is proper, the value must be bigger than :math:`(k+n)/T`,
       where :math:`k` is the number of estimated parameters,
       :math:`n` is the number of observables, and :math:`T` is the
       number of observations.

        NB: The previous method of declaring ``dsge_prior_weight`` as
        a parameter and then calibrating it is now deprecated and will
        be removed in a future release of Dynare. Some of objects
        arising during estimation are stored with their values at the
        mode in ``oo_.dsge_var.posterior_mode``.

    .. option:: dsge_var

       Triggers the estimation of a DSGE-VAR model, where the weight
       of the DSGE prior of the VAR model will be estimated (as in
       *Adjemian et al.(2008)*). The prior on the weight of the DSGE
       prior, ``dsge_prior_weight``, must be defined in the
       ``estimated_params`` section.

       NB: The previous method of declaring ``dsge_prior_weight`` as
       a parameter and then placing it in ``estimated_params`` is now
       deprecated and will be removed in a future release of Dynare.

    .. option:: dsge_varlag = INTEGER

       The number of lags used to estimate a DSGE-VAR model. Default:
       ``4``.

    .. option:: posterior_sampling_method = NAME

       Selects the sampler used to sample from the posterior
       distribution during Bayesian
       estimation. Default:``’random_walk_metropolis_hastings’``.

           ``'random_walk_metropolis_hastings'``

               Instructs Dynare to use the Random-Walk
               Metropolis-Hastings. In this algorithm, the proposal
               density is recentered to the previous draw in every
               step.

           ``'tailored_random_block_metropolis_hastings'``

               Instructs Dynare to use the Tailored randomized block
               (TaRB) Metropolis-Hastings algorithm proposed by *Chib
               and Ramamurthy (2010)* instead of the standard
               Random-Walk Metropolis-Hastings. In this algorithm, at
               each iteration the estimated parameters are randomly
               assigned to different blocks. For each of these blocks
               a mode-finding step is conducted. The inverse Hessian
               at this mode is then used as the covariance of the
               proposal density for a Random-Walk Metropolis-Hastings
               step. If the numerical Hessian is not positive
               definite, the generalized Cholesky decomposition of
               *Schnabel and Eskow (1990)* is used, but without
               pivoting. The TaRB-MH algorithm massively reduces the
               autocorrelation in the MH draws and thus reduces the
               number of draws required to representatively sample
               from the posterior. However, this comes at a
               computational cost as the algorithm takes more time to
               run.

           ``'independent_metropolis_hastings'``

               Use the Independent Metropolis-Hastings algorithm where
               the proposal distribution - in contrast to the Random
               Walk Metropolis-Hastings algorithm - does not depend on
               the state of the chain.

           ``'slice'``

               Instructs Dynare to use the Slice sampler of *Planas,
               Ratto, and Rossi (2015)*. Note that ``'slice'`` is
               incompatible with ``prior_trunc=0``.

    .. option:: posterior_sampler_options = (NAME, VALUE, ...)

       A list of NAME and VALUE pairs. Can be used to set options for
       the posterior sampling methods. The set of available options
       depends on the selected posterior sampling routine (i.e. on the
       value of option :opt:`posterior_sampling_method
       <posterior_sampling_method = NAME>`):

           ``'random_walk_metropolis_hastings'``

               Available options are:

           ``'proposal_distribution'``

               Specifies the statistical distribution used for the
               proposal density.

           ``'rand_multivariate_normal'``

               Use a multivariate normal distribution. This is the default.

           ``'rand_multivariate_student'``

               Use a multivariate student distribution.

           ``'student_degrees_of_freedom'``

               Specifies the degrees of freedom to be used with the
               multivariate student distribution. Default: ``3``.

           .. _usemhcov:

           ``'use_mh_covariance_matrix'``

               Indicates to use the covariance matrix of the draws
               from a previous MCMC run to define the covariance of
               the proposal distribution. Requires the
               :opt:`load_mh_file` option to be specified. Default:
               ``0``.

           .. _scale-file:

           ``'scale_file'``

               Provides the name of a ``_mh_scale.mat`` file storing
               the tuned scale factor from a previous run of
               ``mode_compute=6``.

           .. _savetmp:

           ``'save_tmp_file'``

               Save the MCMC draws into a ``_mh_tmp_blck`` file at the
               refresh rate of the status bar instead of just saving
               the draws when the current ``_mh*_blck`` file is
               full. Default: ``0``

           ``'independent_metropolis_hastings'``

               Takes the same options as in the case of
               ``random_walk_metropolis_hastings``.

           ``'slice'``

           ``'rotated'``

               Triggers rotated slice iterations using a covariance
               matrix from initial burn-in iterations. Requires either
               ``use_mh_covariance_matrix`` or
               ``slice_initialize_with_mode``. Default: ``0``.

           ``'mode_files'``

               For multimodal posteriors, provide the name of a file
               containing a ``nparam`` by ``nmodes`` variable called
               ``xparams`` storing the different modes. This array
               must have one column vector per mode and the estimated
               parameters along the row dimension. With this info, the
               code will automatically trigger the ``rotated`` and
               ``mode`` options. Default: ``[]``.

           ``'slice_initialize_with_mode'``

               The default for slice is to set ``mode_compute=0`` and
               start the chain(s) from a random location in the prior
               space. This option first runs the mode-finder and then
               starts the chain from the mode. Together with
               ``rotated``, it will use the inverse Hessian from the
               mode to perform rotated slice iterations. Default:
               ``0``.

           ``'initial_step_size'``

               Sets the initial size of the interval in the
               stepping-out procedure as fraction of the prior
               support, i.e. the size will be ``initial_step_size *
               (UB-LB)``. ``initial_step_size`` must be a real number
               in the interval ``[0,1]``. Default: ``0.8``.

           ``'use_mh_covariance_matrix'``

               See :ref:`use_mh_covariance_matrix <usemhcov>`. Must be
               used with ``'rotated'``. Default: ``0``.

           ``'save_tmp_file'``

               See :ref:`save_tmp_file <savetmp>`. Default: ``1``.

           ``'tailored_random_block_metropolis_hastings'``

           ``new_block_probability = DOUBLE``

               Specifies the probability of the next parameter
               belonging to a new block when the random blocking in
               the TaRB Metropolis-Hastings algorithm is
               conducted. The higher this number, the smaller is the
               average block size and the more random blocks are
               formed during each parameter sweep. Default: ``0.25``.

           ``mode_compute = INTEGER``

               Specifies the mode-finder run in every iteration for
               every block of the TaRB Metropolis-Hastings
               algorithm. See :opt:`mode_compute <mode_compute =
               INTEGER | FUNCTION_NAME>`. Default: ``4``.

           ``optim = (NAME, VALUE,...)``

               Specifies the options for the mode-finder used in the
               TaRB Metropolis-Hastings algorithm. See :opt:`optim
               <optim = (NAME, VALUE, ...)>`.

           ``'scale_file'``

               See :ref:`scale_file <scale-file>`..

           ``'save_tmp_file'``

               See :ref:`save_tmp_file <savetmp>`. Default: ``1``.

    .. option:: moments_varendo

       Triggers the computation of the posterior distribution of the
       theoretical moments of the endogenous variables. Results are
       stored in ``oo_.PosteriorTheoreticalMoments`` (see
       :mvar:`oo_.PosteriorTheoreticalMoments`). The number of lags in
       the autocorrelation function is controlled by the ``ar``
       option.

    .. option:: contemporaneous_correlation

       See :opt:`contemporaneous_correlation`. Results are stored in
       ``oo_.PosteriorTheoreticalMoments``. Note that the ``nocorr``
       option has no effect.

    .. option:: no_posterior_kernel_density

       Shuts off the computation of the kernel density estimator for
       the posterior objects (see :ref:`density <dens>` field).

    .. option:: conditional_variance_decomposition = INTEGER
                conditional_variance_decomposition = [INTEGER1:INTEGER2]
                conditional_variance_decomposition = [INTEGER1 INTEGER2 ...]

       Computes the posterior distribution of the conditional variance
       decomposition for the specified period(s). The periods must be
       strictly positive. Conditional variances are given by
       :math:`var(y_{t+k}\vert t)`. For period 1, the conditional
       variance decomposition provides the decomposition of the
       effects of shocks upon impact. The results are stored in
       ``oo_.PosteriorTheoreticalMoments.dsge.ConditionalVarianceDecomposition``.. Note
       that this option requires the option ``moments_varendo`` to be
       specified. In the presence of measurement error, the field will
       contain the variance contribution after measurement error has
       been taken out, *i.e.* the decomposition will be conducted of the
       actual as opposed to the measured variables. The variance
       decomposition of the measured variables will be stored in
       ``oo_.PosteriorTheoreticalMoments.dsge.ConditionalVarianceDecompositionME``.


    .. option:: filtered_vars

       Triggers the computation of the posterior distribution of
       filtered endogenous variables/one-step ahead forecasts,
       i.e. :math:`E_{t}{y_{t+1}}`. Results are stored in
       ``oo_.FilteredVariables`` (see below for a description of this
       variable)

    .. option:: smoother

       Triggers the computation of the posterior distribution of
       smoothed endogenous variables and shocks, i.e. the expected
       value of variables and shocks given the information available
       in all observations up to the final date
       (:math:`E_{T}{y_t}`). Results are stored in
       ``oo_.SmoothedVariables``, ``oo_.SmoothedShocks`` and
       ``oo_.SmoothedMeasurementErrors``. Also triggers the
       computation of ``oo_.UpdatedVariables``, which contains the
       estimation of the expected value of variables given the
       information available at the current date
       (:math:`E_{t}{y_t}`). See below for a description of all these
       variables.

    .. option:: forecast = INTEGER

       Computes the posterior distribution of a forecast on INTEGER
       periods after the end of the sample used in estimation. If no
       Metropolis-Hastings is computed, the result is stored in
       variable ``oo_.forecast`` and corresponds to the forecast at
       the posterior mode. If a Metropolis-Hastings is computed, the
       distribution of forecasts is stored in variables
       ``oo_.PointForecast`` and ``oo_.MeanForecast``. See
       :ref:`fore`, for a description of these variables.

    .. option:: tex

       See :opt:`tex`.

    .. option:: kalman_algo = INTEGER

           ``0``

               Automatically use the Multivariate Kalman Filter for
               stationary models and the Multivariate Diffuse Kalman
               Filter for non-stationary models.

           ``1``

               Use the Multivariate Kalman Filter.

           ``2``

               Use the Univariate Kalman Filter.

           ``3``

               Use the Multivariate Diffuse Kalman Filter.

           ``4``

               Use the Univariate Diffuse Kalman Filter.

       Default value is ``0``. In case of missing observations of
       single or all series, Dynare treats those missing values as
       unobserved states and uses the Kalman filter to infer their
       value (see e.g. *Durbin and Koopman (2012)*, Ch. 4.10) This
       procedure has the advantage of being capable of dealing with
       observations where the forecast error variance matrix becomes
       singular for some variable(s). If this happens, the respective
       observation enters with a weight of zero in the log-likelihood,
       i.e. this observation for the respective variable(s) is dropped
       from the likelihood computations (for details see *Durbin and
       Koopman (2012)*, Ch. 6.4 and 7.2.5 and *Koopman and Durbin
       (2000)*). If the use of a multivariate Kalman filter is
       specified and a singularity is encountered, Dynare by default
       automatically switches to the univariate Kalman filter for this
       parameter draw. This behavior can be changed via the
       :opt:`use_univariate_filters_if_singularity_is_detected
       <use_univariate_filters_if_singularity_is_detected = INTEGER>`
       option.

    .. option:: fast_kalman_filter

       Select the fast Kalman filter using Chandrasekhar recursions as
       described by ``Herbst (2015)``. This setting is only used with
       ``kalman_algo=1`` or ``kalman_algo=3``. In case of using the
       diffuse Kalman filter (``kalman_algo=3/lik_init=3``), the
       observables must be stationary. This option is not yet
       compatible with :opt:`analytic_derivation`.

    .. option:: kalman_tol = DOUBLE

       Numerical tolerance for determining the singularity of the
       covariance matrix of the prediction errors during the Kalman
       filter (minimum allowed reciprocal of the matrix condition
       number). Default value is ``1e-10``.

    .. option:: diffuse_kalman_tol = DOUBLE

       Numerical tolerance for determining the singularity of the
       covariance matrix of the prediction errors (:math:`F_{\infty}`)
       and the rank of the covariance matrix of the non-stationary
       state variables (:math:`P_{\infty}`) during the Diffuse Kalman
       filter. Default value is ``1e-6``.

    .. option:: filter_covariance

       Saves the series of one step ahead error of forecast covariance
       matrices. With Metropolis, they are saved in
       :mvar:`oo_.FilterCovariance`, otherwise in
       :mvar:`oo_.Smoother.Variance`. Saves also k-step ahead error of
       forecast covariance matrices if ``filter_step_ahead`` is set.

    .. option:: filter_step_ahead = [INTEGER1:INTEGER2]
                filter_step_ahead = [INTEGER1 INTEGER2 ...]

       Triggers the computation k-step ahead filtered values,
       i.e. :math:`E_{t}{y_{t+k}}`. Stores results in
       ``oo_.FilteredVariablesKStepAhead``. Also stores 1-step ahead
       values in
       ``oo_.FilteredVariables``. ``oo_.FilteredVariablesKStepAheadVariances``
       is stored if ``filter_covariance``.

    .. option:: filter_decomposition

       Triggers the computation of the shock decomposition of the
       above k-step ahead filtered values. Stores results in
       ``oo_.FilteredVariablesShockDecomposition``.

    .. option:: smoothed_state_uncertainty

       Triggers the computation of the variance of smoothed estimates,
       i.e. :math:`var_T(y_t)`. Stores results in
       ``oo_.Smoother.State_uncertainty``.

    .. option:: diffuse_filter

       Uses the diffuse Kalman filter (as described in *Durbin and
       Koopman (2012)* and *Koopman and Durbin (2003)* for the
       multivariate and *Koopman and Durbin (2000)* for the univariate
       filter) to estimate models with non-stationary observed
       variables.

       When ``diffuse_filter`` is used the ``lik_init`` option of
       ``estimation`` has no effect.

       When there are nonstationary exogenous variables in a model,
       there is no unique deterministic steady state. For instance, if
       productivity is a pure random walk:

           .. math::

              a_t = a_{t-1} + e_t

       any value of :math:`\bar a` of :math:`a` is a deterministic
       steady state for productivity. Consequently, the model admits
       an infinity of steady states. In this situation, the user must
       help Dynare in selecting one steady state, except if zero is a
       trivial model’s steady state, which happens when the ``linear``
       option is used in the model declaration. The user can either
       provide the steady state to Dynare using a
       ``steady_state_model`` block (or writing a steady state file)
       if a closed form solution is available, see
       :bck:`steady_state_model`, or specify some constraints on the
       steady state, see
       :ref:`equation_tag_for_conditional_steady_state <eq-tag-ss>`,
       so that Dynare computes the steady state conditionally on some
       predefined levels for the non stationary variables. In both
       cases, the idea is to use dummy values for the steady state
       level of the exogenous non stationary variables.

       Note that the nonstationary variables in the model must be
       integrated processes (their first difference or k-difference
       must be stationary).

    .. option:: selected_variables_only

       Only run the classical smoother on the variables listed just
       after the ``estimation`` command. This option is incompatible
       with requesting classical frequentist forecasts and will be
       overridden in this case. When using Bayesian estimation, the
       smoother is by default only run on the declared endogenous
       variables. Default: run the smoother on all the declared
       endogenous variables.

    .. option:: cova_compute = INTEGER

       When ``0``, the covariance matrix of estimated parameters is
       not computed after the computation of posterior mode (or
       maximum likelihood). This increases speed of computation in
       large models during development, when this information is not
       always necessary. Of course, it will break all successive
       computations that would require this covariance
       matrix. Otherwise, if this option is equal to ``1``, the
       covariance matrix is computed and stored in variable ``hh`` of
       ``MODEL_FILENAME_mode.mat``. Default is ``1``.

    .. option:: solve_algo = INTEGER

       See :ref:`solve_algo <solvalg>`.

    .. option:: order = INTEGER

       Order of approximation, either ``1`` or ``2``. When equal to
       ``2``, the likelihood is evaluated with a particle filter based
       on a second order approximation of the model (see
       *Fernandez-Villaverde and Rubio-Ramirez (2005)*). Default is
       ``1``, i.e. the likelihood of the linearized model is evaluated
       using a standard Kalman filter.

    .. option:: irf = INTEGER

       See :opt:`irf <irf = INTEGER>`. Only used if
       :opt:`bayesian_irf` is passed.

    .. option:: irf_shocks = ( VARIABLE_NAME [[,] VARIABLE_NAME ...] )

        See :opt:`irf_shocks <irf_shocks = ( VARIABLE_NAME [[,]
        VARIABLE_NAME ...] )>`. Only used if :opt:`bayesian_irf` is
        passed.

    .. option:: irf_plot_threshold = DOUBLE

       See :opt:`irf_plot_threshold <irf_plot_threshold =
       DOUBLE>`. Only used if :opt:`bayesian_irf` is passed.

    .. option:: aim_solver

       See :opt:`aim_solver`.

    .. option:: sylvester = OPTION

       See :opt:`sylvester <sylvester = OPTION>`.

    .. option:: sylvester_fixed_point_tol = DOUBLE

       See :opt:`sylvester_fixed_point_tol <sylvester_fixed_point_tol
       = DOUBLE>` .

    .. option:: lyapunov = OPTION

       Determines the algorithm used to solve the Lyapunov equation to
       initialized the variance-covariance matrix of the Kalman filter
       using the steady-state value of state variables. Possible
       values for OPTION are:

           ``default``

               Uses the default solver for Lyapunov equations based on
               Bartels-Stewart algorithm.

           ``fixed_point``

               Uses a fixed point algorithm to solve the Lyapunov
               equation. This method is faster than the ``default``
               one for large scale models, but it could require a
               large amount of iterations.

           ``doubling``

               Uses a doubling algorithm to solve the Lyapunov
               equation (``disclyap_fast``). This method is faster
               than the two previous one for large scale models.

           ``square_root_solver``

               Uses a square-root solver for Lyapunov equations
               (``dlyapchol``). This method is fast for large scale
               models (available under MATLAB if the Control System
               Toolbox is installed; available under Octave if the
               `control <http://octave.sourceforge.net/control/>`_
               package from Octave-Forge is installed)

       Default value is ``default``.

    .. option:: lyapunov_fixed_point_tol = DOUBLE

       This is the convergence criterion used in the fixed point
       Lyapunov solver. Its default value is ``1e-10``.

    .. option:: lyapunov_doubling_tol = DOUBLE

       This is the convergence criterion used in the doubling
       algorithm to solve the Lyapunov equation. Its default value is
       ``1e-16``.

    .. option:: use_penalized_objective_for_hessian

       Use the penalized objective instead of the objective function
       to compute numerically the hessian matrix at the mode. The
       penalties decrease the value of the posterior density (or
       likelihood) when, for some perturbations, Dynare is not able to
       solve the model (issues with steady state existence, Blanchard
       and Kahn conditions, ...). In pratice, the penalized and
       original objectives will only differ if the posterior mode is
       found to be near a region where the model is ill-behaved. By
       default the original objective function is used.

    .. option:: analytic_derivation

       Triggers estimation with analytic gradient. The final hessian
       is also computed analytically. Only works for stationary models
       without missing observations, i.e. for ``kalman_algo<3``.

    .. option:: ar = INTEGER

       See :opt:`ar <ar = INTEGER>`. Only useful in conjunction with
       option ``moments_varendo``.

    .. option:: endogenous_prior

       Use endogenous priors as in *Christiano, Trabandt and Walentin
       (2011)*. The procedure is motivated by sequential Bayesian
       learning. Starting from independent initial priors on the
       parameters, specified in the ``estimated_params`` block, the
       standard deviations observed in a "pre-sample", taken to be the
       actual sample, are used to update the initial priors. Thus, the
       product of the initial priors and the pre-sample likelihood of
       the standard deviations of the observables is used as the new
       prior (for more information, see the technical appendix of
       *Christiano, Trabandt and Walentin (2011)*). This procedure
       helps in cases where the regular posterior estimates, which
       minimize in-sample forecast errors, result in a large
       overprediction of model variable variances (a statistic that is
       not explicitly targeted, but often of particular interest to
       researchers).

    .. option:: use_univariate_filters_if_singularity_is_detected = INTEGER

       Decide whether Dynare should automatically switch to univariate
       filter if a singularity is encountered in the likelihood
       computation (this is the behaviour if the option is equal to
       ``1``). Alternatively, if the option is equal to ``0``, Dynare
       will not automatically change the filter, but rather use a
       penalty value for the likelihood when such a singularity is
       encountered. Default: ``1``.

    .. option:: keep_kalman_algo_if_singularity_is_detected

       With the default
       :opt:`use_univariate_filters_if_singularity_is_detected=1
       <use_univariate_filters_if_singularity_is_detected = INTEGER>`,
       Dynare will switch to the univariate Kalman filter when it
       encounters a singular forecast error variance matrix during
       Kalman filtering. Upon encountering such a singularity for the
       first time, all subsequent parameter draws and computations
       will automatically rely on univariate filter, i.e. Dynare will
       never try the multivariate filter again. Use the
       ``keep_kalman_algo_if_singularity_is_detected`` option to have
       the ``use_univariate_filters_if_singularity_is_detected`` only
       affect the behavior for the current draw/computation.

    .. option:: rescale_prediction_error_covariance

       Rescales the prediction error covariance in the Kalman filter
       to avoid badly scaled matrix and reduce the probability of a
       switch to univariate Kalman filters (which are slower). By
       default no rescaling is done.

    .. option:: qz_zero_threshold = DOUBLE

       See :opt:`qz_zero_threshold <qz_zero_threshold = DOUBLE>`.

    .. option:: taper_steps = [INTEGER1 INTEGER2 ...]

       Percent tapering used for the spectral window in the *Geweke
       (1992,1999)* convergence diagnostics (requires
       :opt:`mh_nblocks=1 <mh_nblocks = INTEGER>`). The tapering is
       used to take the serial correlation of the posterior draws into
       account. Default: ``[4 8 15]``.

    .. option:: geweke_interval = [DOUBLE DOUBLE]

       Percentage of MCMC draws at the beginning and end of the MCMC
       chain taken to compute the *Geweke (1992,1999)* convergence
       diagnostics (requires :opt:`mh_nblocks=1 <mh_nblocks =
       INTEGER>`) after discarding the first :opt:`mh_drop = DOUBLE
       <mh_drop>` percent of draws as a burnin. Default: [0.2 0.5].

    .. option:: raftery_lewis_diagnostics

       Triggers the computation of the *Raftery and Lewis (1992)*
       convergence diagnostics. The goal is deliver the number of
       draws required to estimate a particular quantile of the CDF
       ``q`` with precision ``r`` with a probability ``s``. Typically,
       one wants to estimate the ``q=0.025`` percentile (corresponding
       to a 95 percent HPDI) with a precision of 0.5 percent
       (``r=0.005``) with 95 percent certainty (``s=0.95``). The
       defaults can be changed via :opt:`raftery_lewis_qrs
       <raftery_lewis_qrs = [DOUBLE DOUBLE DOUBLE]>`. Based on the
       theory of first order Markov Chains, the diagnostics will
       provide a required burn-in (``M``), the number of draws after
       the burnin (``N``) as well as a thinning factor that would
       deliver a first order chain (``k``). The last line of the table
       will also deliver the maximum over all parameters for the
       respective values.

    .. option:: raftery_lewis_qrs = [DOUBLE DOUBLE DOUBLE]

       Sets the quantile of the CDF ``q`` that is estimated with
       precision ``r`` with a probability ``s`` in the *Raftery and
       Lewis (1992)* convergence diagnostics. Default: ``[0.025 0.005
       0.95]``.

    .. option:: consider_all_endogenous

       Compute the posterior moments, smoothed variables, k-step ahead
       filtered variables and forecasts (when requested) on all the
       endogenous variables. This is equivalent to manually listing
       all the endogenous variables after the ``estimation`` command.

    .. option:: consider_only_observed

       Compute the posterior moments, smoothed variables, k-step ahead
       filtered variables and forecasts (when requested) on all the
       observed variables. This is equivalent to manually listing all
       the observed variables after the ``estimation`` command.

    .. option:: number_of_particles = INTEGER

       Number of particles used when evaluating the likelihood of a
       non linear state space model. Default: ``1000``.

    .. option:: resampling = OPTION

       Determines if resampling of the particles is done. Possible
       values for OPTION are:

           ``none``

               No resampling.

           ``systematic``

               Resampling at each iteration, this is the default value.

           ``generic``

               Resampling if and only if the effective sample size is
               below a certain level defined by
               :opt:`resampling_threshold <resampling_threshold =
               DOUBLE>` * :opt:`number_of_particles
               <number_of_particles = INTEGER>`.

    .. option:: resampling_threshold = DOUBLE

       A real number between zero and one. The resampling step is
       triggered as soon as the effective number of particles is less
       than this number times the total number of particles (as set by
       :opt:`number_of_particles <number_of_particles =
       INTEGER>`). This option is effective if and only if option
       :opt:`resampling <resampling = OPTION>` has value ``generic``.

    .. option:: resampling_method = OPTION

       Sets the resampling method. Possible values for OPTION are:
       ``kitagawa``, ``stratified`` and ``smooth``.

    .. option:: filter_algorithm = OPTION

       Sets the particle filter algorithm. Possible values for OPTION
       are:

           ``sis``

               Sequential importance sampling algorithm, this is the
               default value.

           ``apf``

               Auxiliary particle filter.

           ``gf``

               Gaussian filter.

           ``gmf``

               Gaussian mixture filter.

           ``cpf``

               Conditional particle filter.

           ``nlkf``

               Use a standard (linear) Kalman filter algorithm with
               the nonlinear measurement and state equations.

    .. option:: proposal_approximation = OPTION

       Sets the method for approximating the proposal
       distribution. Possible values for OPTION are: ``cubature``,
       ``montecarlo`` and ``unscented``. Default value is
       ``unscented``.

    .. option:: distribution_approximation = OPTION

       Sets the method for approximating the particle
       distribution. Possible values for OPTION are: ``cubature``,
       ``montecarlo`` and ``unscented``. Default value is
       ``unscented``.

    .. option:: cpf_weights = OPTION

       Controls the method used to update the weights in conditional
       particle filter, possible values are ``amisanotristani``
       (*Amisano et al. (2010)*) or ``murrayjonesparslow`` (*Murray et
       al. (2013)*). Default value is ``amisanotristani``.

    .. option:: nonlinear_filter_initialization = INTEGER

       Sets the initial condition of the nonlinear filters. By default
       the nonlinear filters are initialized with the unconditional
       covariance matrix of the state variables, computed with the
       reduced form solution of the first order approximation of the
       model. If ``nonlinear_filter_initialization=2``, the nonlinear
       filter is instead initialized with a covariance matrix
       estimated with a stochastic simulation of the reduced form
       solution of the second order approximation of the model. Both
       these initializations assume that the model is stationary, and
       cannot be used if the model has unit roots (which can be seen
       with the :comm:`check` command prior to estimation). If the
       model has stochastic trends, user must use
       ``nonlinear_filter_initialization=3``, the filters are then
       initialized with an identity matrix for the covariance matrix
       of the state variables. Default value is
       ``nonlinear_filter_initialization=1`` (initialization based on
       the first order approximation of the model).

    *Note*

    If no ``mh_jscale`` parameter is used for a parameter in
    ``estimated_params``, the procedure uses ``mh_jscale`` for all
    parameters. If ``mh_jscale`` option isn’t set, the procedure uses
    ``0.2`` for all parameters. Note that if ``mode_compute=6`` is
    used or the ``posterior_sampler_option`` called ``scale_file`` is
    specified, the values set in ``estimated_params`` will be
    overwritten.

    *“Endogenous” prior restrictions*

    It is also possible to impose implicit “endogenous” priors about
    IRFs and moments on the model during estimation. For example, one
    can specify that all valid parameter draws for the model must
    generate fiscal multipliers that are bigger than 1 by specifying
    how the IRF to a government spending shock must look like. The
    prior restrictions can be imposed via ``irf_calibration`` and
    ``moment_calibration`` blocks (see :ref:`irf-momcal`). The way it
    works internally is that any parameter draw that is inconsistent
    with the “calibration” provided in these blocks is discarded,
    i.e. assigned a prior density of 0. When specifying these blocks,
    it is important to keep in mind that one won’t be able to easily
    do ``model_comparison`` in this case, because the prior density
    will not integrate to 1.

    *Output*

    After running estimation, the parameters ``M_.params`` and the
    variance matrix ``M_.Sigma_e`` of the shocks are set to the mode
    for maximum likelihood estimation or posterior mode computation
    without Metropolis iterations. After estimation with Metropolis
    iterations (option ``mh_replic > 0`` or option ``load_mh_file``
    set) the parameters ``M_.params`` and the variance matrix
    ``M_.Sigma_e`` of the shocks are set to the posterior mean.

    Depending on the options, ``estimation`` stores results in various
    fields of the ``oo_`` structure, described below. In the following
    variables, we will adopt the following shortcuts for specific
    field names:


        ``MOMENT_NAME``

            This field can take the following values:

            ``HPDinf``

                Lower bound of a 90% HPD interval [#f3]_.

            ``HPDsup``

                Upper bound of a 90% HPD interval.

            ``HPDinf_ME``

                Lower bound of a 90% HPD interval [#f4]_ for
                observables when taking measurement error into account
                (see e.g. *Christoffel et al. (2010*), p.17).

            ``HPDsup_ME``

                Upper bound of a 90% HPD interval for observables when
                taking measurement error into account.

            ``Mean``

                Mean of the posterior distribution.

            ``Median``

                Median of the posterior distribution.

            ``Std``

                Standard deviation of the posterior distribution.

            ``Variance``

                Variance of the posterior distribution.

            ``deciles``

                Deciles of the distribution.

            .. _dens:

            ``density``

                Non parametric estimate of the posterior density
                following the approach outlined in *Skoeld and Roberts
                (2003)*. First and second columns are respectively
                abscissa and ordinate coordinates.

        ``ESTIMATED_OBJECT``

            This field can take the following values:

            ``measurement_errors_corr``

                Correlation between two measurement errors.

            ``measurement_errors_std``

                Standard deviation of measurement errors.

            ``parameters``

                Parameters.

            ``shocks_corr``

                Correlation between two structural shocks.

            ``shocks_std``

                Standard deviation of structural shocks.


    .. matvar:: oo_.MarginalDensity.LaplaceApproximation

        Variable set by the ``estimation`` command. Stores the marginal
        data density based on the Laplace Approximation.


    .. matvar:: oo_.MarginalDensity.ModifiedHarmonicMean

        Variable set by the ``estimation command``, if it is used with
        ``mh_replic > 0`` or ``load_mh_file`` option. Stores the
        marginal data density based on *Geweke (1999)* Modified
        Harmonic Mean estimator.


    .. matvar:: oo_.posterior.optimization

        Variable set by the ``estimation`` command if mode-finding is
        used. Stores the results at the mode. Fields are of the form::

            oo_.posterior.optimization.OBJECT

        where OBJECT is one of the following:

           ``mode``

               Parameter vector at the mode.

           ``Variance``

               Inverse Hessian matrix at the mode or MCMC jumping
               covariance matrix when used with the
               :opt:`MCMC_jumping_covariance <mcmc_jumping_covariance
               = OPTION>` option.

           ``log_density``

               Log likelihood (ML)/log posterior density (Bayesian) at the
               mode when used with ``mode_compute>0``.


    .. matvar:: oo_.posterior.metropolis

        Variable set by the ``estimation`` command if ``mh_replic>0`` is
        used. Fields are of the form::

            oo_.posterior.metropolis.OBJECT

        where OBJECT is one of the following:

            ``mean``

                Mean parameter vector from the MCMC.

            ``Variance``

                Covariance matrix of the parameter draws in the MCMC.


    .. matvar:: oo_.FilteredVariables

        Variable set by the ``estimation`` command, if it is used with the
        ``filtered_vars`` option.

        After an estimation without Metropolis, fields are of the form::

            oo_.FilteredVariables.VARIABLE_NAME

        After an estimation with Metropolis, fields are of the form::

            oo_.FilteredVariables.MOMENT_NAME.VARIABLE_NAME


    .. matvar:: oo_.FilteredVariablesKStepAhead

        Variable set by the ``estimation`` command, if it is used with
        the ``filter_step_ahead`` option. The k-steps are stored along
        the rows while the columns indicate the respective
        variables. The third dimension of the array provides the
        observation for which the forecast has been made. For example,
        if ``filter_step_ahead=[1 2 4]`` and ``nobs=200``, the element
        (3,5,204) stores the four period ahead filtered value of
        variable 5 computed at time t=200 for time t=204. The periods
        at the beginning and end of the sample for which no forecasts
        can be made, e.g. entries (1,5,1) and (1,5,204) in the
        example, are set to zero. Note that in case of Bayesian
        estimation the variables will be ordered in the order of
        declaration after the estimation command (or in general
        declaration order if no variables are specified here). In case
        of running the classical smoother, the variables will always
        be ordered in general declaration order. If the
        :opt:`selected_variables_only` option is specified with the
        classical smoother, non-requested variables will be simply
        left out in this order.


    .. matvar:: oo_.FilteredVariablesKStepAheadVariances

        Variable set by the ``estimation`` command, if it is used with
        the ``filter_step_ahead option``. It is a 4 dimensional array
        where the k-steps are stored along the first dimension, while
        the fourth dimension of the array provides the observation for
        which the forecast has been made. The second and third
        dimension provide the respective variables. For example, if
        ``filter_step_ahead=[1 2 4]`` and ``nobs=200``, the element
        (3,4,5,204) stores the four period ahead forecast error
        covariance between variable 4 and variable 5, computed at time
        t=200 for time t=204. Padding with zeros and variable ordering
        is analogous to ``oo_.FilteredVariablesKStepAhead``.

    .. matvar:: oo_.Filtered_Variables_X_step_ahead

        Variable set by the ``estimation`` command, if it is used with the
        ``filter_step_ahead option`` in the context of Bayesian
        estimation. Fields are of the form::

            oo_.Filtered_Variables_X_step_ahead.VARIABLE_NAME

        The n-th entry stores the k-step ahead filtered variable computed
        at time n for time n+k.


    .. matvar:: oo_.FilteredVariablesShockDecomposition

        Variable set by the ``estimation`` command, if it is used with
        the ``filter_step_ahead`` option. The k-steps are stored along
        the rows while the columns indicate the respective
        variables. The third dimension corresponds to the shocks in
        declaration order. The fourth dimension of the array provides
        the observation for which the forecast has been made. For
        example, if ``filter_step_ahead=[1 2 4]`` and ``nobs=200``,
        the element (3,5,2,204) stores the contribution of the second
        shock to the four period ahead filtered value of variable 5
        (in deviations from the mean) computed at time t=200 for time
        t=204. The periods at the beginning and end of the sample for
        which no forecasts can be made, e.g. entries (1,5,1) and
        (1,5,204) in the example, are set to zero. Padding with zeros
        and variable ordering is analogous to
        ``oo_.FilteredVariablesKStepAhead``.

    .. matvar:: oo_.PosteriorIRF.dsge

        Variable set by the ``estimation`` command, if it is used with the
        ``bayesian_irf`` option. Fields are of the form::

            oo_.PosteriorIRF.dsge.MOMENT_NAME.VARIABLE_NAME_SHOCK_NAME


    .. matvar:: oo_.SmoothedMeasurementErrors

        Variable set by the ``estimation`` command, if it is used with the
        ``smoother`` option. Fields are of the form::

            oo_.SmoothedMeasurementErrors.VARIABLE_NAME


    .. matvar:: oo_.SmoothedShocks

        Variable set by the ``estimation`` command (if used with the
        ``smoother`` option), or by the ``calib_smoother`` command.

        After an estimation without Metropolis, or if computed by
        ``calib_smoother``, fields are of the form::

            oo_.SmoothedShocks.VARIABLE_NAME

        After an estimation with Metropolis, fields are of the form::

            oo_.SmoothedShocks.MOMENT_NAME.VARIABLE_NAME


    .. matvar:: oo_.SmoothedVariables

        Variable set by the ``estimation`` command (if used with the
        ``smoother`` option), or by the ``calib_smoother`` command.

        After an estimation without Metropolis, or if computed by
        ``calib_smoother``, fields are of the form::

            oo_.SmoothedVariables.VARIABLE_NAME

        After an estimation with Metropolis, fields are of the form::

            oo_.SmoothedVariables.MOMENT_NAME.VARIABLE_NAME


    .. matvar:: oo_.UpdatedVariables

        Variable set by the ``estimation`` command (if used with the
        ``smoother`` option), or by the ``calib_smoother``
        command. Contains the estimation of the expected value of
        variables given the information available at the current date.

        After an estimation without Metropolis, or if computed by
        ``calib_smoother``, fields are of the form::

            oo_.UpdatedVariables.VARIABLE_NAME

        After an estimation with Metropolis, fields are of the form::

            oo_.UpdatedVariables.MOMENT_NAME.VARIABLE_NAME


    .. matvar:: oo_.FilterCovariance

        Three-dimensional array set by the ``estimation`` command if
        used with the ``smoother`` and Metropolis, if the
        ``filter_covariance`` option has been requested. Contains the
        series of one-step ahead forecast error covariance matrices
        from the Kalman smoother. The ``M_.endo_nbr`` times
        ``M_.endo_nbr`` times ``T+1`` array contains the variables in
        declaration order along the first two dimensions. The third
        dimension of the array provides the observation for which the
        forecast has been made. Fields are of the form::

            oo_.FilterCovariance.MOMENT_NAME

        Note that density estimation is not supported.


    .. matvar:: oo_.Smoother.Variance

        Three-dimensional array set by the ``estimation`` command (if
        used with the ``smoother``) without Metropolis, or by the
        ``calib_smoother`` command, if the ``filter_covariance``
        option has been requested. Contains the series of one-step
        ahead forecast error covariance matrices from the Kalman
        smoother. The ``M_.endo_nbr`` times ``M_.endo_nbr`` times
        ``T+1`` array contains the variables in declaration order
        along the first two dimensions. The third dimension of the
        array provides the observation for which the forecast has been
        made.


    .. matvar:: oo_.Smoother.State_uncertainty

        Three-dimensional array set by the ``estimation`` command (if
        used with the ``smoother`` option) without Metropolis, or by
        the ``calib_smoother`` command, if the
        ``smoothed_state_uncertainty`` option has been
        requested. Contains the series of covariance matrices for the
        state estimate given the full data from the Kalman
        smoother. The ``M_.endo_nbr`` times ``M_.endo_nbr`` times
        ``T`` array contains the variables in declaration order along
        the first two dimensions. The third dimension of the array
        provides the observation for which the smoothed estimate has
        been made.


    .. matvar:: oo_.Smoother.SteadyState

        Variable set by the ``estimation`` command (if used with the
        ``smoother``) without Metropolis, or by the
        ````calib_smoother`` command. Contains the steady state
        component of the endogenous variables used in the smoother in
        order of variable declaration.


    .. matvar:: oo_.Smoother.TrendCoeffs

        Variable set by the ````estimation`` command (if used with the
        ``smoother``) without Metropolis, or by the ``calib_smoother``
        command. Contains the trend coefficients of the observed
        variables used in the smoother in order of declaration of the
        observed variables.


    .. matvar:: oo_.Smoother.Trend

        Variable set by the ``estimation command`` (if used with the
        ``smoother`` option), or by the ````calib_smoother``
        command. Contains the trend component of the variables used in
        the smoother.

        Fields are of the form::

            oo_.Smoother.Trend.VARIABLE_NAME


    .. matvar:: oo_.Smoother.Constant

        Variable set by the ``estimation`` command (if used with the
        ``smoother`` option), or by the ``calib_smoother``
        command. Contains the constant part of the endogenous
        variables used in the smoother, accounting e.g. for the data
        mean when using the prefilter option.

        Fields are of the form::

            oo_.Smoother.Constant.VARIABLE_NAME


    .. matvar:: oo_.Smoother.loglinear

        Indicator keeping track of whether the smoother was run with
        the :ref:`loglinear <logl>` option and thus whether stored
        smoothed objects are in logs.


    .. matvar:: oo_.PosteriorTheoreticalMoments

        Variable set by the ``estimation`` command, if it is used with the
        ``moments_varendo`` option. Fields are of the form::

            oo_.PosteriorTheoreticalMoments.dsge.THEORETICAL_MOMENT.ESTIMATED_OBJECT.MOMENT_NAME.VARIABLE_NAME

        where *THEORETICAL_MOMENT* is one of the following:

            ``covariance``

                Variance-covariance of endogenous variables.

            ``contemporaneous_correlation``

                Contemporaneous correlation of endogenous variables when the
                :opt:`contemporaneous_correlation` option is specified.

            ``correlation``

                Auto- and cross-correlation of endogenous variables. Fields
                are vectors with correlations from 1 up to order
                ``options_.ar``.

            .. _VarianceDecomposition:

            ``VarianceDecomposition``

                Decomposition of variance (unconditional variance, i.e. at
                horizon infinity). [#f5]_

            ``VarianceDecompositionME``

                Same as `VarianceDecomposition`_, but contains
                theh decomposition of the measured as opposed to the
                actual variable. The joint contribution of the
                measurement error will be saved in a field named
                ``ME``.

            .. _ConditionalVarianceDecomposition:

            ``ConditionalVarianceDecomposition``

                Only if the ``conditional_variance_decomposition``
                option has been specified. In the presence of
                measurement error, the field will contain the variance
                contribution after measurement error has been taken
                out, i.e. the decomposition will be conducted of the
                actual as opposed to the measured variables.

            ``ConditionalVarianceDecompositionME``

                Only if the ``conditional_variance_decomposition``
                option has been specified. Same as
                `ConditionalVarianceDecomposition`_, but contains the
                decomposition of the measured as opposed to the actual
                variable. The joint contribution of the measurement
                error will be saved in a field names ``ME``.


    .. matvar:: oo_.posterior_density

        Variable set by the ``estimation`` command, if it is used with
        ``mh_replic > 0`` or ``load_mh_file`` option. Fields are of
        the form::

            oo_.posterior_density.PARAMETER_NAME


    .. matvar:: oo_.posterior_hpdinf

        Variable set by the ``estimation`` command, if it is used with
        ``mh_replic > 0`` or ``load_mh_file`` option. Fields are of
        the form::

            oo_.posterior_hpdinf.ESTIMATED_OBJECT.VARIABLE_NAME


    .. matvar:: oo_.posterior_hpdsup

        Variable set by the ``estimation`` command, if it is used with
        ``mh_replic > 0`` or ``load_mh_file`` option. Fields are of the
        form::

            oo_.posterior_hpdsup.ESTIMATED_OBJECT.VARIABLE_NAME


    .. matvar:: oo_.posterior_mean

        Variable set by the ``estimation`` command, if it is used with
        ``mh_replic > 0`` or ``load_mh_file`` option. Fields are of the
        form::

            oo_.posterior_mean.ESTIMATED_OBJECT.VARIABLE_NAME


    .. matvar:: oo_.posterior_mode

        Variable set by the ``estimation`` command during
        mode-finding. Fields are of the form::

            oo_.posterior_mode.ESTIMATED_OBJECT.VARIABLE_NAME


    .. matvar:: oo_.posterior_std_at_mode

        Variable set by the ``estimation`` command during mode-finding. It
        is based on the inverse Hessian at ``oo_.posterior_mode``. Fields
        are of the form::

            oo_.posterior_std_at_mode.ESTIMATED_OBJECT.VARIABLE_NAME


    .. matvar:: oo_.posterior_std

        Variable set by the ``estimation`` command, if it is used with
        ``mh_replic > 0`` or ``load_mh_file`` option. Fields are of the
        form::

            oo_.posterior_std.ESTIMATED_OBJECT.VARIABLE_NAME


    .. matvar:: oo_.posterior_var

        Variable set by the ``estimation`` command, if it is used with
        ``mh_replic > 0`` or ``load_mh_file`` option. Fields are of the
        form::

            oo_.posterior_var.ESTIMATED_OBJECT.VARIABLE_NAME


    .. matvar:: oo_.posterior_median

        Variable set by the ``estimation`` command, if it is used with
        ``mh_replic > 0`` or ``load_mh_file`` option. Fields are of the
        form::

            oo_.posterior_median.ESTIMATED_OBJECT.VARIABLE_NAME


    *Example*

    Here are some examples of generated variables::

        oo_.posterior_mode.parameters.alp
        oo_.posterior_mean.shocks_std.ex
        oo_.posterior_hpdsup.measurement_errors_corr.gdp_conso


    .. matvar:: oo_.dsge_var.posterior_mode

        Structure set by the ``dsge_var`` option of the ``estimation``
        command after mode_compute.

        The following fields are saved:

            ``PHI_tilde``

                Stacked posterior DSGE-BVAR autoregressive matrices at the
                mode (equation (28) of *Del Negro and Schorfheide (2004)*).

            ``SIGMA_u_tilde``

                Posterior covariance matrix of the DSGE-BVAR at the mode
                (equation (29) of *Del Negro and Schorfheide (2004)*).

            ``iXX``

                Posterior population moments in the DSGE-BVAR at the mode (
                :math:`inv(\lambda T \Gamma_{XX}^*+ X'X)`).

            ``prior``

                Structure storing the DSGE-BVAR prior.

            ``PHI_star``

                Stacked prior DSGE-BVAR autoregressive matrices at the
                mode (equation (22) of *Del Negro and Schorfheide
                (2004)*).

            ``SIGMA_star``

                Prior covariance matrix of the DSGE-BVAR at the mode
                (equation (23) of *Del Negro and Schorfheide (2004)*).

            ``ArtificialSampleSize``

                Size of the artifical prior sample ( :math:`inv(\lambda T)`).

            ``DF``

                Prior degrees of freedom ( :math:`inv(\lambda T-k-n)`).

            ``iGXX_star``

                Inverse of the theoretical prior “covariance” between
                X and X (:math:`\Gamma_{xx}^*` in *Del Negro and
                Schorfheide (2004)*).


    .. matvar:: oo_.RecursiveForecast

        Variable set by the ``forecast`` option of the ``estimation``
        command when used with the nobs = [INTEGER1:INTEGER2] option (see
        :opt:`nobs <nobs = [INTEGER1:INTEGER2]>`).

        Fields are of the form::

            oo_.RecursiveForecast.FORECAST_OBJECT.VARIABLE_NAME

        where ``FORECAST_OBJECT`` is one of the following [#f6]_ :

        ``Mean``

            Mean of the posterior forecast distribution.

        ``HPDinf/HPDsup``

            Upper/lower bound of the 90% HPD interval taking into account
            only parameter uncertainty (corresponding to
            :mvar:`oo_.MeanForecast`).

        ``HPDTotalinf/HPDTotalsup``.

            Upper/lower bound of the 90% HPD interval taking into account
            both parameter and future shock uncertainty (corresponding to
            :mvar:`oo_.PointForecast`)

        ``VARIABLE_NAME`` contains a matrix of the following size:
        number of time periods for which forecasts are requested using
        the ``nobs = [INTEGER1:INTEGER2]`` option times the number of
        forecast horizons requested by the forecast option. i.e., the
        row indicates the period at which the forecast is performed
        and the column the respective k-step ahead forecast. The
        starting periods are sorted in ascending order, not in
        declaration order.

    .. matvar:: oo_.convergence.geweke

        Variable set by the convergence diagnostics of the ``estimation``
        command when used with ``mh_nblocks=1`` option (see
        :opt:`mh_nblocks <mh_nblocks = INTEGER>`).

        Fields are of the form::

            oo_.convergence.geweke.VARIABLE_NAME.DIAGNOSTIC_OBJECT

        where *DIAGNOSTIC_OBJECT* is one of the following:

        ``posteriormean``

            Mean of the posterior parameter distribution.

        ``posteriorstd``

            Standard deviation of the posterior parameter distribution.

        ``nse_iid``

            Numerical standard error (NSE) under the assumption of iid draws.

        ``rne_iid``

            Relative numerical efficiency (RNE) under the assumption
            of iid draws.

        ``nse_x``

            Numerical standard error (NSE) when using an x% taper.

        ``rne_x``

            Relative numerical efficiency (RNE) when using an x% taper.

        ``pooled_mean``

            Mean of the parameter when pooling the beginning and end parts
            of the chain specified in :opt:`geweke_interval
            <geweke_interval = [DOUBLE DOUBLE]>` and weighting them with
            their relative precision. It is a vector containing the
            results under the iid assumption followed by the ones using
            the ``taper_steps`` option (see :opt:`taper_steps <taper_steps
            = [INTEGER1 INTEGER2 ...]>`).

        ``pooled_nse``

            NSE of the parameter when pooling the beginning and end parts
            of the chain and weighting them with their relative
            precision. See ``pooled_mean``.

        ``prob_chi2_test``

            p-value of a chi-squared test for equality of means in the
            beginning and the end of the MCMC chain. See
            ``pooled_mean``. A value above 0.05 indicates that the null
            hypothesis of equal means and thus convergence cannot be
            rejected at the 5 percent level. Differing values along the
            ``taper_steps`` signal the presence of significant
            autocorrelation in draws. In this case, the estimates using a
            higher tapering are usually more reliable.

.. command:: unit_root_vars VARIABLE_NAME...;

    |br| This command is deprecated. Use ``estimation`` option
    ``diffuse_filter`` instead for estimating a model with
    non-stationary observed variables or ``steady`` option ``nocheck``
    to prevent ``steady`` to check the steady state returned by your
    steady state file.

Dynare also has the ability to estimate Bayesian VARs:

.. command:: bvar_density ;

    |br| Computes the marginal density of an estimated BVAR model, using
    Minnesota priors.

    See ``bvar-a-la-sims.pdf``, which comes with Dynare distribution,
    for more information on this command.



Model Comparison
================

.. command:: model_comparison FILENAME[(DOUBLE)]...;
             model_comparison (marginal_density = ESTIMATOR) FILENAME[(DOUBLE)]...;

    |br| This command computes odds ratios and estimate a posterior density
    over a collection of models (see e.g. *Koop (2003)*, Ch. 1). The
    priors over models can be specified as the *DOUBLE* values,
    otherwise a uniform prior over all models is assumed. In contrast
    to frequentist econometrics, the models to be compared do not need
    to be nested. However, as the computation of posterior odds ratios
    is a Bayesian technique, the comparison of models estimated with
    maximum likelihood is not supported.

    It is important to keep in mind that model comparison of this type
    is only valid with proper priors. If the prior does not integrate
    to one for all compared models, the comparison is not valid. This
    may be the case if part of the prior mass is implicitly truncated
    because Blanchard and Kahn conditions (instability or
    indeterminacy of the model) are not fulfilled, or because for some
    regions of the parameters space the deterministic steady state is
    undefined (or Dynare is unable to find it). The compared marginal
    densities should be renormalized by the effective prior mass, but
    this not done by Dynare: it is the user’s responsibility to make
    sure that model comparison is based on proper priors. Note that,
    for obvious reasons, this is not an issue if the compared marginal
    densities are based on Laplace approximations.

    *Options*

    .. option:: marginal_density = ESTIMATOR

         Specifies the estimator for computing the marginal data
         density. *ESTIMATOR* can take one of the following two values:
         ``laplace`` for the Laplace estimator or
         ``modifiedharmonicmean`` for the *Geweke (1999)* Modified
         Harmonic Mean estimator. Default value: ``laplace``

    *Output*

    The results are stored in ``oo_.Model_Comparison``, which is
    described below.

    *Example*

        ::

            model_comparison my_model(0.7) alt_model(0.3);

        This example attributes a 70% prior over ``my_model`` and 30%
        prior over ``alt_model``.


.. matvar:: oo_.Model_Comparison

    Variable set by the ``model_comparison`` command. Fields are of
    the form::

        oo_.Model_Comparison.FILENAME.VARIABLE_NAME

    where FILENAME is the file name of the model and VARIABLE_NAME is
    one of the following:

        ``Prior``

            (Normalized) prior density over the model.

        ``Log_Marginal_Density``

            Logarithm of the marginal data density.

        ``Bayes_Ratio``

            Ratio of the marginal data density of the model relative
            to the one of the first declared model

        ``Posterior_Model_Probability``

            Posterior probability of the respective model.


Shock Decomposition
===================

.. command:: shock_decomposition [VARIABLE_NAME]...;
             shock_decomposition (OPTIONS...) [VARIABLE_NAME]...;

    |br| This command computes the historical shock decomposition for a
    given sample based on the Kalman smoother, i.e. it decomposes the
    historical deviations of the endogenous variables from their
    respective steady state values into the contribution coming from
    the various shocks. The ``variable_names`` provided govern for
    which variables the decomposition is plotted.

    Note that this command must come after either ``estimation`` (in
    case of an estimated model) or ``stoch_simul`` (in case of a
    calibrated model).

    *Options*

    .. option:: parameter_set = OPTION

        Specify the parameter set to use for running the
        smoother. Possible values for OPTION are:

            * ``calibration``
            * ``prior_mode``
            * ``prior_mean``
            * ``posterior_mode``
            * ``posterior_mean``
            * ``posterior_median``
            * ``mle_mode``

        Note that the parameter set used in subsequent commands like
        ``stoch_simul`` will be set to the specified
        ``parameter_set``. Default value: ``posterior_mean`` if
        Metropolis has been run, ``mle_mode`` if MLE has been run.

    .. option:: datafile = FILENAME

        See :ref:`datafile <dataf>`. Useful when computing the shock
        decomposition on a calibrated model.

    .. option:: first_obs = INTEGER

        See :opt:`first_obs <first_obs = INTEGER>`.

    .. option:: nobs = INTEGER

        See :opt:`nobs <nobs = INTEGER>`.

    .. option:: use_shock_groups [= STRING]

        Uses shock grouping defined by the string instead of
        individual shocks in the decomposition. The groups of shocks
        are defined in the :bck:`shock_groups` block.

    .. option:: colormap = STRING

        Controls the ``colormap`` used for the shocks decomposition
        graphs. See colormap in Matlab/Octave manual for valid
        arguments.

    .. option:: nograph

        See :opt:`nograph`. Suppresses the display and creation only
        within the ``shock_decomposition`` command, but does not
        affect other commands. See :comm:`plot_shock_decomposition`
        for plotting graphs.

    .. option:: init_state = BOOLEAN

        If equal to 0, the shock decomposition is computed conditional
        on the smoothed state variables in period ``0``, i.e. the
        smoothed shocks starting in period 1 are used. If equal to
        ``1``, the shock decomposition is computed conditional on the
        smoothed state variables in period 1. Default: ``0``.

    *Output*

    .. matvar:: oo_.shock_decomposition

        The results are stored in the field
        ``oo_.shock_decomposition``, which is a three dimensional
        array. The first dimension contains the ``M_.endo_nbr``
        endogenous variables. The second dimension stores in the first
        ``M_.exo_nbr`` columns the contribution of the respective
        shocks. Column ``M_.exo_nbr+1`` stores the contribution of the
        initial conditions, while column ``M_.exo_nbr+2`` stores the
        smoothed value of the respective endogenous variable in
        deviations from their steady state, i.e. the mean and trends
        are subtracted. The third dimension stores the time
        periods. Both the variables and shocks are stored in the order
        of declaration, i.e. ``M_.endo_names`` and ``M_.exo_names``,
        respectively.


.. block:: shock_groups ;
           shock_groups(OPTIONS...);

    |br| Shocks can be regrouped for the purpose of shock
    decomposition. The composition of the shock groups is written in a
    block delimited by ``shock_groups`` and ``end``.

    Each line defines a group of shocks as a list of exogenous variables::

        SHOCK_GROUP_NAME   = VARIABLE_1 [[,] VARIABLE_2 [,]...];
        'SHOCK GROUP NAME' = VARIABLE_1 [[,] VARIABLE_2 [,]...];

    *Options*

    .. option:: name = NAME

        Specifies a name for the following definition of shock
        groups. It is possible to use several ``shock_groups`` blocks
        in a model file, each grouping being identified by a different
        name. This name must in turn be used in the
        ``shock_decomposition`` command.

    *Example*

        ::

            varexo e_a, e_b, e_c, e_d;
            ...

            shock_groups(name=group1);
            supply = e_a, e_b;
            'aggregate demand' = e_c, e_d;
            end;

            shock_decomposition(use_shock_groups=group1);

        This example defines a shock grouping with the name
        ``group1``, containing a set of supply and demand shocks and
        conducts the shock decomposition for these two groups.


.. command:: realtime_shock_decomposition [VARIABLE_NAME]...;
             realtime_shock_decomposition (OPTIONS...) [VARIABLE_NAME]...;

    |br| This command computes the realtime historical shock
    decomposition for a given sample based on the Kalman smoother. For
    each period :math:`T=[\texttt{presample},\ldots,\texttt{nobs}]`,
    it recursively computes three objects:

        * Real-time historical shock decomposition :math:`Y(t\vert T)`
          for :math:`t=[1,\ldots,T]`, i.e. without observing data in
          :math:`[T+1,\ldots,\texttt{nobs}]`. This results in a
          standard shock decomposition being computed for each
          additional datapoint becoming available after ``presample``.
        * Forecast shock decomposition :math:`Y(T+k\vert T)` for
          :math:`k=[1,\ldots,forecast]`, i.e. the :math:`k`-step ahead
          forecast made for every :math:`T` is decomposed in its shock
          contributions.
        * Real-time conditional shock decomposition of the difference
          between the real-time historical shock decomposition and the
          forecast shock decomposition. If :opt:`vintage <vintage =
          INTEGER>` is equal to ``0``, it computes the effect of
          shocks realizing in period :math:`T`, i.e. decomposes
          :math:`Y(T\vert T)-Y(T\vert T-1)`. Put differently, it
          conducts a :math:`1`-period ahead shock decomposition from
          :math:`T-1` to :math:`T`, by decomposing the update step of
          the Kalman filter. If ``vintage>0`` and smaller than
          ``nobs``, the decomposition is conducted of the forecast
          revision :math:`Y(T+k\vert T+k)-Y(T+k\vert T)`.

    Like :comm:`shock_decomposition` it decomposes the historical
    deviations of the endogenous variables from their respective
    steady state values into the contribution coming from the various
    shocks. The ``variable_names`` provided govern for which variables
    the decomposition is plotted.

    Note that this command must come after either ``estimation`` (in
    case of an estimated model) or ``stoch_simul`` (in case of a
    calibrated model).

    *Options*

    .. option:: parameter_set = OPTION

        See :opt:`parameter_set <parameter_set = OPTION>` for the
        possible values.

    .. option:: datafile = FILENAME

        See :ref:`datafile <dataf>`.

    .. option:: first_obs = INTEGER

        See :opt:`first_obs <first_obs = INTEGER>`.

    .. option:: nobs = INTEGER

        See :opt:`nobs <nobs = INTEGER>`.

    .. option:: use_shock_groups [= STRING]

        See :opt:`use_shock_groups <use_shock_groups [= STRING]>`.

    .. option:: colormap = STRING

        See :opt:`colormap <colormap = STRING>`.

    .. option:: nograph

        See :opt:`nograph`. Only shock decompositions are computed and
        stored in ``oo_.realtime_shock_decomposition``,
        ``oo_.conditional_shock_decomposition`` and
        ``oo_.realtime_forecast_shock_decomposition`` but no plot is
        made (See :comm:`plot_shock_decomposition`).

    .. option:: presample = INTEGER

        First data point from which recursive realtime shock
        decompositions are computed, i.e. for
        :math:`T=[\texttt{presample} \ldots \texttt{nobs}]`.

    .. option:: forecast = INTEGER

        Compute shock decompositions up to :math:`T+k` periods,
        i.e. get shock contributions to k-step ahead forecasts.

    .. option:: save_realtime = INTEGER_VECTOR

        Choose for which vintages to save the full realtime shock
        decomposition. Default: ``0``..

    *Output*

    .. matvar:: oo_.realtime_shock_decomposition

        Structure storing the results of realtime historical
        decompositions. Fields are three-dimensional arrays with the
        first two dimension equal to the ones of
        :mvar:`oo_.shock_decomposition`. The third dimension stores
        the time periods and is therefore of size
        ``T+forecast``. Fields are of the form::

            oo_.realtime_shock_decomposition.OBJECT

        where OBJECT is one of the following:

            ``pool``

                Stores the pooled decomposition, i.e. for every
                real-time shock decomposition terminal period
                :math:`T=[\texttt{presample},\ldots,\texttt{nobs}]` it
                collects the last period’s decomposition
                :math:`Y(T\vert T)` (see also
                :comm:`plot_shock_decomposition`). The third dimension
                of the array will have size ``nobs+forecast``.

            ``time_*``

                Stores the vintages of realtime historical shock
                decompositions if ``save_realtime`` is used. For
                example, if ``save_realtime=[5]`` and ``forecast=8``,
                the third dimension will be of size ``13``.

    .. matvar:: oo_.realtime_conditional_shock_decomposition

        Structure storing the results of real-time conditional
        decompositions. Fields are of the form::

            oo_.realtime_conditional_shock_decomposition.OBJECT

        where OBJECT is one of the following:

            ``pool``

                Stores the pooled real-time conditional shock
                decomposition, i.e. collects the decompositions of
                :math:`Y(T\vert T)-Y(T\vert T-1)` for the terminal
                periods
                :math:`T=[\texttt{presample},\ldots,\texttt{nobs}]`. The
                third dimension is of size ``nobs``.

            ``time_*``

                Store the vintages of :math:`k`-step conditional
                forecast shock decompositions :math:`Y(t\vert T+k)`,
                for :math:`t=[T \ldots T+k]`. See :opt:`vintage
                <vintage = INTEGER>`. The third dimension is of size
                ``1+forecast``.

    .. matvar:: oo_.realtime_forecast_shock_decomposition

        Structure storing the results of realtime forecast
        decompositions. Fields are of the form::

            oo_.realtime_forecast_shock_decomposition.OBJECT

        where ``OBJECT`` is one of the following:

            ``pool``

                Stores the pooled real-time forecast decomposition of
                the :math:`1`-step ahead effect of shocks on the
                :math:`1`-step ahead prediction, i.e. :math:`Y(T\vert
                T-1)`.

            ``time_*``

                Stores the vintages of :math:`k`-step out-of-sample
                forecast shock decompositions, i.e. :math:`Y(t\vert
                T)`, for :math:`t=[T \ldots T+k]`. See :opt:`vintage
                <vintage = INTEGER>`.


.. command:: plot_shock_decomposition [VARIABLE_NAME]...;
             plot_shock_decomposition (OPTIONS...) [VARIABLE_NAME]...;

    |br| This command plots the historical shock decomposition already
    computed by ``shock_decomposition`` or
    ``realtime_shock_decomposition``. For that reason, it must come
    after one of these commands. The ``variable_names`` provided
    govern which variables the decomposition is plotted for.

    Further note that, unlike the majority of Dynare commands, the
    options specified below are overwritten with their defaults before
    every call to ``plot_shock_decomposition``. Hence, if you want to
    reuse an option in a subsequent call to
    ``plot_shock_decomposition``, you must pass it to the command
    again.

    *Options*

    .. option:: use_shock_groups [= STRING]

        See :opt:`use_shock_groups <use_shock_groups [= STRING]>`.

    .. option:: colormap = STRING

        See :opt:`colormap <colormap = STRING>`.

    .. option:: nodisplay

        See :opt:`nodisplay`.

    .. option:: graph_format = FORMAT
                graph_format = ( FORMAT, FORMAT... )

        See :opt:`graph_format <graph_format = FORMAT>`.

    .. option:: detail_plot

        Plots shock contributions using subplots, one per shock (or
        group of shocks). Default: not activated

    .. option:: interactive

        Under MATLAB, add uimenus for detailed group plots. Default:
        not activated

    .. option:: screen_shocks

        For large models (i.e. for models with more than 16 shocks),
        plots only the shocks that have the largest historical
        contribution for chosen selected
        ``variable_names``. Historical contribution is ranked by the
        mean absolute value of all historical contributions.

    .. option:: steadystate

        If passed, the the :math:`y`-axis value of the zero line in
        the shock decomposition plot is translated to the steady state
        level. Default: not activated

    .. option:: type = qoq | yoy | aoa

        For quarterly data, valid arguments are: ``qoq`` for
        quarter-on-quarter plots, ``yoy`` for year-on-year plots of
        growth rates, ``aoa`` for annualized variables, i.e. the value
        in the last quarter for each year is plotted. Default value:
        empty, i.e. standard period-on-period plots (``qoq`` for
        quarterly data).

    .. option:: fig_name = STRING

        Specifies a user-defined keyword to be appended to the default
        figure name set by ``plot_shock_decomposition``. This can
        avoid to overwrite plots in case of sequential calls to
        ``plot_shock_decomposition``.

    .. option:: write_xls

        Saves shock decompositions to Excel-file in the main
        directory, named
        ``FILENAME_shock_decomposition_TYPE_FIG_NAME.xls``. This
        option requires your system to be configured to be able to
        write Excel files. [#f7]_

    .. option:: realtime = INTEGER

        Which kind of shock decomposition to plot. INTEGER can take
        the following values:

            * ``0``: standard historical shock decomposition. See
              :comm:`shock_decomposition`.
            * ``1``: realtime historical shock decomposition. See
              :comm:`realtime_shock_decomposition`.
            * ``2``: conditional realtime shock decomposition. See
              :comm:`realtime_shock_decomposition`.
            * ``3``: realtime forecast shock decomposition. See
              :comm:`realtime_shock_decomposition`.

        If no vintage is requested, i.e. ``vintage=0`` then the pooled
        objects from :comm:`realtime_shock_decomposition` will be
        plotted and the respective vintage otherwise. Default: ``0``.

    .. option:: vintage = INTEGER

        Selects a particular data vintage in
        :math:`[presample,\ldots,nobs]` for which to plot the results
        from :comm:`realtime_shock_decomposition` selected via the
        :opt:`realtime <realtime = INTEGER>` option. If the standard
        historical shock decomposition is selected (``realtime=0``),
        ``vintage`` will have no effect. If ``vintage=0`` the pooled
        objects from :comm:`realtime_shock_decomposition` will be
        plotted. If ``vintage>0``, it plots the shock decompositions
        for vintage :math:`T=\texttt{vintage}` under the following
        scenarios:

            * ``realtime=1``: the full vintage shock decomposition
              :math:`Y(t\vert T)` for :math:`t=[1,\ldots,T]`
            * ``realtime=2``: the conditional forecast shock
              decomposition from :math:`T`, i.e. plots
              :math:`Y(T+j\vert T+j)` and the shock contributions
              needed to get to the data :math:`Y(T+j)` conditional on
              :math:`T=` vintage, with
              :math:`j=[0,\ldots,\texttt{forecast}]`.
            * ``realtime=3``: plots unconditional forecast shock
              decomposition from :math:`T`, i.e. :math:`Y(T+j\vert
              T)`, where :math:`T=\texttt{vintage}` and
              :math:`j=[0,\ldots,\texttt{forecast}]`.

        Default: ``0``.


Calibrated Smoother
===================

Dynare can also run the smoother on a calibrated model:

.. command:: calib_smoother [VARIABLE_NAME]...;
             calib_smoother (OPTIONS...) [VARIABLE_NAME]...;

    |br| This command computes the smoothed variables (and possible
    the filtered variables) on a calibrated model.

    A datafile must be provided, and the observable variables declared
    with ``varobs``. The smoother is based on a first-order
    approximation of the model.

    By default, the command computes the smoothed variables and shocks
    and stores the results in ``oo_.SmoothedVariables` and
    ``oo_.SmoothedShocks``. It also fills ``oo_.UpdatedVariables``.

    *Options*

    .. option:: datafile = FILENAME

        See :ref:`datafile <dataf>`.

    .. option:: filtered_vars

        Triggers the computation of filtered variables. See
        :opt:`filtered_vars`, for more details.

    .. option:: filter_step_ahead = [INTEGER1:INTEGER2]

        See :opt:`filter_step_ahead <filter_step_ahead = [INTEGER1:INTEGER2]>`.

    .. option:: prefilter = INTEGER

        See :opt:`prefilter <prefilter = INTEGER>`.

    .. option:: parameter_set = OPTION

        See :opt:`parameter_set` for the possible values.

    .. option:: loglinear

        See :ref:`loglinear <logl>`.

    .. option:: first_obs = INTEGER

        See :opt:`first_obs <first_obs = INTEGER>`.

    .. option:: filter_decomposition

        See :opt:`filter_decomposition`.

    .. option:: diffuse_filter = INTEGER

        See :opt:`diffuse_filter`.

    .. option:: diffuse_kalman_tol = DOUBLE

        See :opt:`diffuse_kalman_tol <diffuse_kalman_tol = DOUBLE>`.

.. _fore:

Forecasting
===========

On a calibrated model, forecasting is done using the ``forecast``
command. On an estimated model, use the ``forecast`` option of
``estimation`` command.

It is also possible to compute forecasts on a calibrated or estimated
model for a given constrained path of the future endogenous
variables. This is done, from the reduced form representation of the
DSGE model, by finding the structural shocks that are needed to match
the restricted paths. Use :comm:`conditional_forecast`,
:bck:`conditional_forecast_paths` and :comm:`plot_conditional_forecast` for
that purpose.

Finally, it is possible to do forecasting with a Bayesian VAR using
the :comm:`bvar_forecast` command.

.. command:: forecast [VARIABLE_NAME...];
             forecast (OPTIONS...) [VARIABLE_NAME...];

    |br| This command computes a simulation of a stochastic model from
    an arbitrary initial point.

    When the model also contains deterministic exogenous shocks, the
    simulation is computed conditionally to the agents knowing the
    future values of the deterministic exogenous variables.

    ``forecast`` must be called after ``stoch_simul``.

    ``forecast`` plots the trajectory of endogenous variables. When a
    list of variable names follows the command, only those variables
    are plotted. A 90% confidence interval is plotted around the mean
    trajectory. Use option ``conf_sig`` to change the level of the
    confidence interval.

    *Options*

    .. option:: periods = INTEGER

        Number of periods of the forecast. Default: ``5``.

    .. _confsig:

    .. option:: conf_sig = DOUBLE

        Level of significance for confidence interval. Default: ``0.90``.

    .. option:: nograph

        See :opt:`nograph`.

    .. option:: nodisplay

        See :opt:`nodisplay`.

    .. option:: graph_format = FORMAT
                graph_format = ( FORMAT, FORMAT... )

        See :opt:`graph_format = FORMAT`.

    *Initial Values*

    ``forecast`` computes the forecast taking as initial values the
    values specified in ``histval`` (see :bck:`histval`). When no
    ``histval`` block is present, the initial values are the one
    stated in ``initval``. When ``initval`` is followed by command
    ``steady``, the initial values are the steady state (see
    :comm:`steady`).

    *Output*

    The results are stored in ``oo_.forecast``, which is described below.

    *Example*

        ::

            varexo_det tau;

            varexo e;
            ...
            shocks;
            var e; stderr 0.01;
            var tau;
            periods 1:9;
            values -0.15;
            end;

            stoch_simul(irf=0);

            forecast;


    .. matvar:: oo_.forecast

        Variable set by the ``forecast`` command, or by the
        ``estimation`` command if used with the ``forecast`` option
        and if no Metropolis-Hastings has been computed (in that case,
        the forecast is computed for the posterior mode). Fields are
        of the form::

            oo_.forecast.FORECAST_MOMENT.VARIABLE_NAME

        where ``FORECAST_MOMENT`` is one of the following:

            ``HPDinf``

                Lower bound of a 90% HPD interval [#f8]_ of forecast
                due to parameter uncertainty, but ignoring the effect
                of measurement error on observed variables.

            ``HPDsup``

                Upper bound of a 90% HPD forecast interval due to
                parameter uncertainty, but ignoring the effect of
                measurement error on observed variables.

            ``HPDinf_ME``

                Lower bound of a 90% HPD interval [#f9]_ of forecast
                for observed variables due to parameter uncertainty
                and measurement error.

            ``HPDsup_ME``

                Upper bound of a 90% HPD interval of forecast for
                observed variables due to parameter uncertainty and
                measurement error.

            ``Mean``

                Mean of the posterior distribution of forecasts.

            ``Median``

                Median of the posterior distribution of forecasts.

            ``Std``

                Standard deviation of the posterior distribution of forecasts.

    .. matvar:: oo_.PointForecast

        Set by the ``estimation`` command, if it is used with the
        ``forecast`` option and if either ``mh_replic > 0`` or the
        ``load_mh_file`` option are used.

        Contains the distribution of forecasts taking into account the
        uncertainty about both parameters and shocks.

        Fields are of the form::

            oo_.PointForecast.MOMENT_NAME.VARIABLE_NAME

    .. matvar:: oo_.MeanForecast

        Set by the ``estimation`` command, if it is used with the
        ``forecast`` option and if either ``mh_replic > 0`` or
        ``load_mh_file`` option are used.

        Contains the distribution of forecasts where the uncertainty
        about shocks is averaged out. The distribution of forecasts
        therefore only represents the uncertainty about parameters.

        Fields are of the form::

            oo_.MeanForecast.MOMENT_NAME.VARIABLE_NAME


.. command:: conditional_forecast (OPTIONS...) [VARIABLE_NAME...];

    |br| This command computes forecasts on an estimated or calibrated
    model for a given constrained path of some future endogenous
    variables. This is done using the reduced form first order
    state-space representation of the DSGE model by finding the
    structural shocks that are needed to match the restricted
    paths. Consider the an augmented state space representation that
    stacks both predetermined and non-predetermined variables into a
    vector :math:`y_{t}`:

        .. math::

           y_t=Ty_{t-1}+R\varepsilon_t

    Both :math:`y_t` and :math:`\varepsilon_t` are split up into
    controlled and uncontrolled ones to get:

        .. math::

           y_t(contr\_vars)=Ty_{t-1}(contr\_vars)+R(contr\_vars,uncontr\_shocks)\varepsilon_t(uncontr\_shocks) + R(contr\_vars,contr\_shocks)\varepsilon_t(contr\_shocks)

    which can be solved algebraically for :math:`\varepsilon_t(contr\_shocks)`.

    Using these controlled shocks, the state-space representation can
    be used for forecasting. A few things need to be noted. First, it
    is assumed that controlled exogenous variables are fully under
    control of the policy maker for all forecast periods and not just
    for the periods where the endogenous variables are controlled. For
    all uncontrolled periods, the controlled exogenous variables are
    assumed to be 0. This implies that there is no forecast
    uncertainty arising from these exogenous variables in uncontrolled
    periods. Second, by making use of the first order state space
    solution, even if a higher-order approximation was performed, the
    conditional forecasts will be based on a first order
    approximation. Third, although controlled exogenous variables are
    taken as instruments perfectly under the control of the
    policy-maker, they are nevertheless random and unforeseen shocks
    from the perspective of the households. That is, households are in
    each period surprised by the realization of a shock that keeps the
    controlled endogenous variables at their respective level. Fourth,
    keep in mind that if the structural innovations are correlated,
    because the calibrated or estimated covariance matrix has non zero
    off diagonal elements, the results of the conditional forecasts
    will depend on the ordering of the innovations (as declared after
    ``varexo``). As in VAR models, a Cholesky decomposition is used to
    factorize the covariance matrix and identify orthogonal
    impulses. It is preferable to declare the correlations in the
    model block (explicitly imposing the identification restrictions),
    unless you are satisfied with the implicit identification
    restrictions implied by the Cholesky decomposition.

    This command has to be called after ``estimation`` or ``stoch_simul``.

    Use :bck:`conditional_forecast_paths` block to give the list of
    constrained endogenous, and their constrained future path. Option
    ``controlled_varexo`` is used to specify the structural shocks
    which will be matched to generate the constrained path.

    Use :comm:`plot_conditional_forecast` to graph the results.

    *Options*

    .. option:: parameter_set = OPTION

        Specify the parameter set to use for the forecasting. Possible
        values for OPTION are:

            * ``calibration``
            * ``prior_mode``
            * ``prior_mean``
            * ``posterior_mode``
            * ``posterior_mean``
            * ``posterior_median``

        No default value, mandatory option. Note that in case of
        estimated models, ``conditional_forecast`` does not support
        the ``prefilter`` option.

    .. option:: controlled_varexo = (VARIABLE_NAME...)

        Specify the exogenous variables to use as control
        variables. No default value, mandatory option.

    .. option:: periods = INTEGER

        Number of periods of the forecast. Default:
        ``40``. ``periods`` cannot be smaller than the number of
        constrained periods.

    .. option:: replic = INTEGER

        Number of simulations. Default: ``5000``.

    .. option:: conf_sig = DOUBLE

        Level of significance for confidence interval. Default: ``0.90``.

    *Output*

    The results are not stored in the ``oo_`` structure but in a
    separate structure ``forecasts``, described below, saved to the
    hard disk into a file called ``conditional_forecasts.mat.``

    *Example*

        ::

            var y a;
            varexo e u;
            ...
            estimation(...);

            conditional_forecast_paths;
            var y;
            periods 1:3, 4:5;
            values 2, 5;
            var a;
            periods 1:5;
            values 3;
            end;

            conditional_forecast(parameter_set = calibration, controlled_varexo = (e, u), replic = 3000);

            plot_conditional_forecast(periods = 10) a y;


    .. matvar:: forecasts.cond

        Variable set by the ``conditional_forecast`` command. It
        stores the conditional forecasts. Fields are ``periods+1`` by
        ``1`` vectors storing the steady state (time 0) and the
        subsequent ``periods`` forecasts periods. Fields are of the
        form::

            forecasts.cond.FORECAST_MOMENT.VARIABLE_NAME

        where FORECAST_MOMENT is one of the following:

            ``Mean``

                Mean of the conditional forecast distribution.

            ``ci``

                Confidence interval of the conditional forecast
                distribution. The size corresponds to ``conf_sig``.


    .. matvar:: forecasts.uncond

        Variable set by the ``conditional_forecast`` command. It stores
        the unconditional forecasts. Fields are of the form::

            forecasts.uncond.FORECAST_MOMENT.VARIABLE_NAME


    .. matvar:: forecasts.instruments

        Variable set by the ``conditional_forecast command``. Stores
        the names of the exogenous instruments.


    .. matvar:: forecasts.controlled_variables

        Variable set by the ``conditional_forecast`` command. Stores
        the position of the constrained endogenous variables in
        declaration order.


    .. matvar:: forecasts.controlled_exo_variables

        Variable set by the ``conditional_forecast`` command. Stores
        the values of the controlled exogenous variables underlying
        the conditional forecasts to achieve the constrained
        endogenous variables. Fields are ``[number of constrained
        periods]`` by ``1`` vectors and are of the form::

            forecasts.controlled_exo_variables.FORECAST_MOMENT.SHOCK_NAME

    .. matvar:: forecasts.graphs

        Variable set by the ``conditional_forecast`` command. Stores
        the information for generating the conditional forecast plots.


.. block:: conditional_forecast_paths ;

    |br| Describes the path of constrained endogenous, before calling
    ``conditional_forecast``. The syntax is similar to deterministic
    shocks in ``shocks``, see ``conditional_forecast`` for an example.

    The syntax of the block is the same as for the deterministic
    shocks in the ``shocks`` blocks (see :ref:`shocks-exo`). Note that
    you need to specify the full path for all constrained endogenous
    variables between the first and last specified period. If an
    intermediate period is not specified, a value of 0 is
    assumed. That is, if you specify only values for periods 1 and 3,
    the values for period 2 will be 0. Currently, it is not possible
    to have uncontrolled intermediate periods. In case of the presence
    of ``observation_trends``, the specified controlled path for these
    variables needs to include the trend component. When using the
    :ref:`loglinear <logl>` option, it is necessary to specify the
    logarithm of the controlled variables.


.. command:: plot_conditional_forecast [VARIABLE_NAME...];
             plot_conditional_forecast (periods = INTEGER) [VARIABLE_NAME...];

    |br| Plots the conditional (plain lines) and unconditional (dashed
    lines) forecasts.

    To be used after ``conditional_forecast``.

    *Options*

    .. option:: periods = INTEGER

        Number of periods to be plotted. Default: equal to periods in
        ``conditional_forecast``. The number of periods declared in
        ``plot_conditional_forecast`` cannot be greater than the one
        declared in ``conditional_forecast``.

.. command:: bvar_forecast ;

    |br| This command computes (out-of-sample) forecasts for an
    estimated BVAR model, using Minnesota priors.

    See ``bvar-a-la-sims.pdf``, which comes with Dynare distribution,
    for more information on this command.


If the model contains strong non-linearities or if some perfectly
expected shocks are considered, the forecasts and the conditional
forecasts can be computed using an extended path method. The forecast
scenario describing the shocks and/or the constrained paths on some
endogenous variables should be build. The first step is the forecast
scenario initialization using the function ``init_plan``:

.. matcomm:: HANDLE = init_plan (DATES);

    Creates a new forecast scenario for a forecast period (indicated
    as a dates class, see :ref:`dates class members
    <dates-members>`). This function return a handle on the new
    forecast scenario.

The forecast scenario can contain some simple shocks on the exogenous
variables. This shocks are described using the function
``basic_plan``:

.. matcomm:: HANDLE = basic_plan (HANDLE, 'VAR_NAME', 'SHOCK_TYPE', DATES, MATLAB VECTOR OF DOUBLE | [DOUBLE | EXPR [DOUBLE | | EXPR] ] );

    Adds to the forecast scenario a shock on the exogenous variable
    indicated between quotes in the second argument. The shock type
    has to be specified in the third argument between quotes:
    ’surprise’ in case of an unexpected shock or ’perfect_foresight’
    for a perfectly anticipated shock. The fourth argument indicates
    the period of the shock using a dates class (see :ref:`dates class
    members <dates-members>`). The last argument is the shock path
    indicated as a Matlab vector of double. This function return the
    handle of the updated forecast scenario.

The forecast scenario can also contain a constrained path on an
endogenous variable. The values of the related exogenous variable
compatible with the constrained path are in this case computed. In
other words, a conditional forecast is performed. This kind of shock
is described with the function ``flip_plan``:

.. matcomm:: HANDLE = flip_plan (HANDLE, 'VAR_NAME, 'VAR_NAME', 'SHOCK_TYPE', DATES, MATLAB VECTOR OF DOUBLE | [DOUBLE | EXPR [DOUBLE | | EXPR] ] );

    Adds to the forecast scenario a constrained path on the endogenous
    variable specified between quotes in the second argument. The
    associated exogenous variable provided in the third argument
    between quotes, is considered as an endogenous variable and its
    values compatible with the constrained path on the endogenous
    variable will be computed. The nature of the expectation on the
    constrained path has to be specified in the fourth argument
    between quotes: ’surprise’ in case of an unexpected path or
    ’perfect_foresight’ for a perfectly anticipated path. The fifth
    argument indicates the period where the path of the endogenous
    variable is constrained using a dates class (see :ref:`dates class
    members <dates-members>`). The last argument contains the
    constrained path as a Matlab vector of double. This function
    return the handle of the updated forecast scenario.

Once the forecast scenario if fully described, the forecast is
computed with the command ``det_cond_forecast``:

.. matcomm:: DSERIES = det_cond_forecast (HANDLE[, DSERIES [, DATES]]);

    Computes the forecast or the conditional forecast using an
    extended path method for the given forecast scenario (first
    argument). The past values of the endogenous and exogenous
    variables provided with a dseries class (see :ref:`dseries class
    members <dseries-members>`) can be indicated in the second
    argument. By default, the past values of the variables are equal
    to their steady-state values. The initial date of the forecast can
    be provided in the third argument. By default, the forecast will
    start at the first date indicated in the ``init_plan
    command``. This function returns a dset containing the historical
    and forecast values for the endogenous and exogenous variables.


*Example*

    ::

        % conditional forecast using extended path method
        % with perfect foresight on r path

        var y r;
        varexo e u;
        ...
        smoothed = dseries('smoothed_variables.csv');

        fplan = init_plan(2013Q4:2029Q4);
        fplan = flip_plan(fplan, 'y', 'u', 'surprise', 2013Q4:2014Q4,  [1 1.1 1.2 1.1 ]);
        fplan = flip_plan(fplan, 'r', 'e', 'perfect_foresight', 2013Q4:2014Q4,  [2 1.9 1.9 1.9 ]);

        dset_forecast = det_cond_forecast(fplan, smoothed);

        plot(dset_forecast.{'y','u'});
        plot(dset_forecast.{'r','e'});


.. command:: smoother2histval [(OPTIONS...)]

    The purpose of this command is to construct initial conditions
    (for a subsequent simulation) that are the smoothed values of a
    previous estimation.

    More precisely, after an estimation run with the ``smoother``
    option, ``smoother2histval`` will extract the smoothed values
    (from ``oo_.SmoothedVariables``, and possibly from
    ``oo_.SmoothedShocks`` if there are lagged exogenous), and will
    use these values to construct initial conditions (as if they had
    been manually entered through ``histval``).

    *Options*

    .. option:: period = INTEGER

        Period number to use as the starting point for the subsequent
        simulation. It should be between 1 and the number of
        observations that were used to produce the smoothed
        values. Default: the last observation.

    .. option:: infile = FILENAME

        Load the smoothed values from a ``_results.mat`` file created
        by a previous Dynare run. Default: use the smoothed values
        currently in the global workspace.

    .. option:: invars = ( VARIABLE_NAME [VARIABLE_NAME ...] )

        A list of variables to read from the smoothed values. It can
        contain state endogenous variables, and also exogenous
        variables having a lag. Default: all the state endogenous
        variables, and all the exogenous variables with a lag.

    .. option:: outfile = FILENAME

        Write the initial conditions to a file. Default: write the
        initial conditions in the current workspace, so that a
        simulation can be performed.

    .. option:: outvars = ( VARIABLE_NAME [VARIABLE_NAME ...] )

        A list of variables which will be given the initial
        conditions. This list must have the same length than the list
        given to ``invars``, and there will be a one-to-one mapping
        between the two list. Default: same value as option
        ``invars``.

    *Use cases*

    There are three possible ways of using this command:

        * Everything in a single file: run an estimation with a
          smoother, then run ``smoother2histval`` (without the
          ``infile`` and ``outfile`` options), then run a stochastic
          simulation.
        * In two files: in the first file, run the smoother and then
          run ``smoother2histval`` with the ``outfile`` option; in the
          second file, run ``histval_file`` to load the initial
          conditions, and run a (deterministic or stochastic)
          simulation.
        * In two files: in the first file, run the smoother; in the
          second file, run ``smoother2histval`` with the ``infile``
          option equal to the ``_results.mat`` file created by the
          first file, and then run a (deterministic or stochastic)
          simulation.


Optimal policy
==============

Dynare has tools to compute optimal policies for various types of
objectives. ``ramsey_model`` computes automatically the First Order
Conditions (FOC) of a model, given the ``planner_objective``. You can
then use other standard commands to solve, estimate or simulate this
new, expanded model.

Alternatively, you can either solve for optimal policy under
commitment with ``ramsey_policy``, for optimal policy under discretion
with ``discretionary_policy`` or for optimal simple rule with ``osr``
(also implying commitment).

.. command:: osr [VARIABLE_NAME...];
             osr (OPTIONS...) [VARIABLE_NAME...];

    |br| This command computes optimal simple policy rules for
    linear-quadratic problems of the form:

        .. math::

             \min_\gamma E(y'_tWy_t)

    such that:

        .. math::

             A_1 E_ty_{t+1}+A_2 y_t+ A_3 y_{t-1}+C e_t=0

    where:

        * :math:`E` denotes the unconditional expectations operator;
        * :math:`\gamma` are parameters to be optimized. They must be
          elements of the matrices :math:`A_1`, :math:`A_2`,
          :math:`A_3`, i.e. be specified as parameters in the
          ``params`` command and be entered in the ``model`` block;
        * :math:`y` are the endogenous variables, specified in the
          ``var`` command, whose (co)-variance enters the loss
          function;
        * :math:`e` are the exogenous stochastic shocks, specified in
          the ``varexo``- ommand;
        * :math:`W` is the weighting matrix;

    The linear quadratic problem consists of choosing a subset of
    model parameters to minimize the weighted (co)-variance of a
    specified subset of endogenous variables, subject to a linear law
    of motion implied by the first order conditions of the model. A
    few things are worth mentioning. First, :math:`y` denotes the
    selected endogenous variables’ deviations from their steady state,
    i.e. in case they are not already mean 0 the variables entering
    the loss function are automatically demeaned so that the centered
    second moments are minimized. Second, ``osr`` only solves linear
    quadratic problems of the type resulting from combining the
    specified quadratic loss function with a first order approximation
    to the model’s equilibrium conditions. The reason is that the
    first order state-space representation is used to compute the
    unconditional (co)-variances. Hence, ``osr`` will automatically
    select ``order=1``. Third, because the objective involves
    minimizing a weighted sum of unconditional second moments, those
    second moments must be finite. In particular, unit roots in
    :math:`y` are not allowed.

    The subset of the model parameters over which the optimal simple
    rule is to be optimized, :math:`\gamma`, must be listed with
    ``osr_params``.

    The weighting matrix :math:`W` used for the quadratic objective
    function is specified in the ``optim_weights`` block. By attaching
    weights to endogenous variables, the subset of endogenous
    variables entering the objective function, :math:`y`, is
    implicitly specified.

    The linear quadratic problem is solved using the numerical
    optimizer specified with :opt:`opt_algo <opt_algo = INTEGER>`.

    *Options*

    The ``osr`` command will subsequently run ``stoch_simul`` and
    accepts the same options, including restricting the endogenous
    variables by listing them after the command, as ``stoch_simul``
    (see :ref:`stoch-sol`) plus

    .. option:: opt_algo = INTEGER

        Specifies the optimizer for minimizing the objective
        function. The same solvers as for ``mode_compute`` (see
        :opt:`mode_compute <mode_compute = INTEGER | FUNCTION_NAME>`)
        are available, except for ``5``, ``6``, and ``10``.

    .. option:: optim = (NAME, VALUE, ...)

        A list of NAME`` and VALUE pairs. Can be used to set options
        for the optimization routines. The set of available options
        depends on the selected optimization routine (i.e. on the
        value of option :opt:`opt_algo <opt_algo = INTEGER>`). See
        :opt:`optim <optim = (NAME, VALUE, ...)>`.

    .. option:: maxit = INTEGER

        Determines the maximum number of iterations used in
        ``opt_algo=4``. This option is now deprecated and will be
        removed in a future release of Dynare. Use ``optim`` instead
        to set optimizer-specific values. Default: ``1000``.

    .. option:: tolf = DOUBLE

        Convergence criterion for termination based on the function
        value used in ``opt_algo=4``. Iteration will cease when it
        proves impossible to improve the function value by more than
        tolf. This option is now deprecated and will be removed in a
        future release of Dynare. Use ``optim`` instead to set
        optimizer-specific values. Default: ``e-7``.

    .. option:: silent_optimizer

        See :opt:`silent_optimizer`.

    .. option:: huge_number = DOUBLE

        Value for replacing the infinite bounds on parameters by
        finite numbers. Used by some optimizers for numerical reasons
        (see :opt:`huge_number <huge_number = DOUBLE>`). Users need to
        make sure that the optimal parameters are not larger than this
        value. Default: ``1e7``.

    The value of the objective is stored in the variable
    ``oo_.osr.objective_function`` and the value of parameters at the
    optimum is stored in ``oo_.osr.optim_params``. See below for more
    details.

    After running ``osr`` the parameters entering the simple rule will
    be set to their optimal value so that subsequent runs of
    ``stoch_simul`` will be conducted at these values.


.. command:: osr_params PARAMETER_NAME...;

    |br| This command declares parameters to be optimized by ``osr``.

.. block:: optim_weights ;

    |br| This block specifies quadratic objectives for optimal policy problems.

    More precisely, this block specifies the nonzero elements of the
    weight matrix :math:`W` used in the quadratic form of the
    objective function in ``osr``.

    An element of the diagonal of the weight matrix is given by a line
    of the form::

        VARIABLE_NAME EXPRESSION;

    An off-the-diagonal element of the weight matrix is given
    by a line of the form::

        VARIABLE_NAME,  VARIABLE_NAME EXPRESSION;

*Example*

    ::

        var y inflation r;
        varexo y_ inf_;

        parameters delta sigma alpha kappa gammarr gammax0 gammac0 gamma_y_ gamma_inf_;

        delta =  0.44;
        kappa =  0.18;
        alpha =  0.48;
        sigma = -0.06;

        gammarr = 0;
        gammax0 = 0.2;
        gammac0 = 1.5;
        gamma_y_ = 8;
        gamma_inf_ = 3;

        model(linear);
        y  = delta * y(-1)  + (1-delta)*y(+1)+sigma *(r - inflation(+1)) + y_;
        inflation  =   alpha * inflation(-1) + (1-alpha) * inflation(+1) + kappa*y + inf_;
        r = gammax0*y(-1)+gammac0*inflation(-1)+gamma_y_*y_+gamma_inf_*inf_;
        end;

        shocks;
        var y_; stderr 0.63;
        var inf_; stderr 0.4;
        end;

        optim_weights;
        inflation 1;
        y 1;
        y, inflation 0.5;
        end;

        osr_params gammax0 gammac0 gamma_y_ gamma_inf_;
        osr y;


.. block:: osr_params_bounds ;

    |br| This block declares lower and upper bounds for parameters in
    the optimal simple rule. If not specified the optimization is
    unconstrained.

    Each line has the following syntax::

        PARAMETER_NAME, LOWER_BOUND, UPPER_BOUND;

    Note that the use of this block requires the use of a constrained
    optimizer, i.e. setting :opt:`opt_algo <opt_algo = INTEGER>` to
    ``1``, ``2``, ``5`` or ``9``.

    *Example*

        ::

            osr_params_bounds;
            gamma_inf_, 0, 2.5;
            end;

            osr(solve_algo=9) y;


.. matvar:: oo_.osr.objective_function

    After an execution of the ``osr`` command, this variable contains
    the value of the objective under optimal policy.

.. matvar:: oo_.osr.optim_params

    After an execution of the ``osr`` command, this variable contains
    the value of parameters at the optimum, stored in fields of the
    form ``oo_.osr.optim_params.PARAMETER_NAME``.

.. matvar:: M_.osr.param_names

    After an execution of the ``osr`` command, this cell contains the
    names of the parameters.

.. matvar:: M_.osr.param_indices

    After an execution of the ``osr`` command, this vector contains
    the indices of the OSR parameters in ``M_.params``.

.. matvar:: M_.osr.param_bounds

    After an execution of the ``osr`` command, this two by number of
    OSR parameters matrix contains the lower and upper bounds of the
    parameters in the first and second column, respectively.

.. matvar:: M_.osr.variable_weights

    After an execution of the ``osr`` command, this sparse matrix
    contains the weighting matrix associated with the variables in the
    objective function.

.. matvar:: M_.osr.variable_indices

    After an execution of the ``osr`` command, this vector contains
    the indices of the variables entering the objective function in
    ``M_.endo_names``.


.. command:: ramsey_model (OPTIONS...);

    |br| This command computes the First Order Conditions for maximizing
    the policy maker objective function subject to the constraints
    provided by the equilibrium path of the private economy.

    The planner objective must be declared with the
    ``planner_objective`` command.

    This command only creates the expanded model, it doesn’t perform
    any computations. It needs to be followed by other instructions to
    actually perform desired computations. Note that it is the only
    way to perform perfect foresight simulation of the Ramsey policy
    problem.

    See :ref:`aux-variables`, for an explanation of how Lagrange
    multipliers are automatically created.

    *Options*

    This command accepts the following options:

    .. option:: planner_discount = EXPRESSION

        Declares or reassigns the discount factor of the central
        planner ``optimal_policy_discount_factor``. Default: ``1.0``.

    .. option:: instruments = (VARIABLE_NAME,...)

        Declares instrument variables for the computation of the
        steady state under optimal policy. Requires a
        ``steady_state_model`` block or a ``_steadystate.m`` file. See
        below.

    *Steady state*

    Dynare takes advantage of the fact that the Lagrange multipliers
    appear linearly in the equations of the steady state of the model
    under optimal policy. Nevertheless, it is in general very
    difficult to compute the steady state with simply a numerical
    guess in ``initval`` for the endogenous variables.

    It greatly facilitates the computation, if the user provides an
    analytical solution for the steady state (in
    ``steady_state_model`` block or in a ``_steadystate.m`` file). In
    this case, it is necessary to provide a steady state solution
    CONDITIONAL on the value of the instruments in the optimal policy
    problem and declared with option ``instruments``. Note that
    choosing the instruments is partly a matter of interpretation and
    you can choose instruments that are handy from a mathematical
    point of view but different from the instruments you would refer
    to in the analysis of the paper. A typical example is choosing
    inflation or nominal interest rate as an instrument.


.. block:: ramsey_constraints ;

    |br| This block lets you define constraints on the variables in
    the Ramsey problem. The constraints take the form of a variable,
    an inequality operator (> or <) and a constant.

    *Example*

        ::

            ramsey_constraints;
            i > 0;
            end;


.. command:: ramsey_policy [VARIABLE_NAME...];
             ramsey_policy (OPTIONS...) [VARIABLE_NAME...];

    |br| This command computes the first order approximation of the
    policy that maximizes the policy maker’s objective function
    subject to the constraints provided by the equilibrium path of the
    private economy and under commitment to this optimal policy. The
    Ramsey policy is computed by approximating the equilibrium system
    around the perturbation point where the Lagrange multipliers are
    at their steady state, i.e. where the Ramsey planner acts as if
    the initial multipliers had been set to 0 in the distant past,
    giving them time to converge to their steady state
    value. Consequently, the optimal decision rules are computed
    around this steady state of the endogenous variables and the
    Lagrange multipliers.

    This first order approximation to the optimal policy conducted by
    Dynare is not to be confused with a naive linear quadratic
    approach to optimal policy that can lead to spurious welfare
    rankings (see *Kim and Kim (2003)*). In the latter, the optimal
    policy would be computed subject to the first order approximated
    FOCs of the private economy. In contrast, Dynare first computes
    the FOCs of the Ramsey planner’s problem subject to the nonlinear
    constraints that are the FOCs of the private economy and only then
    approximates these FOCs of planner’s problem to first
    order. Thereby, the second order terms that are required for a
    second-order correct welfare evaluation are preserved.

    Note that the variables in the list after the ``ramsey_policy``
    command can also contain multiplier names. In that case, Dynare
    will for example display the IRFs of the respective multipliers
    when ``irf>0``.

    The planner objective must be declared with the planner_objective command.

    See :ref:`aux-variables`, for an explanation of how this operator
    is handled internally and how this affects the output.

    *Options*

    This command accepts all options of ``stoch_simul``, plus:

    .. option:: planner_discount = EXPRESSION

        See :opt:`planner_discount <planner_discount = EXPRESSION>`.

    .. option:: instruments = (VARIABLE_NAME,...)

        Declares instrument variables for the computation of the
        steady state under optimal policy. Requires a
        ``steady_state_model`` block or a ``_steadystate.m`` file. See
        below.

    Note that only a first order approximation of the optimal Ramsey
    policy is available, leading to a second-order accurate welfare
    ranking (i.e. ``order=1`` must be specified).

    *Output*

    This command generates all the output variables of
    ``stoch_simul``. For specifying the initial values for the
    endogenous state variables (except for the Lagrange multipliers),
    see :bck:`histval`.

    .. _plan-obj:

    In addition, it stores the value of planner objective function
    under Ramsey policy in ``oo_.planner_objective_value``, given the
    initial values of the endogenous state variables. If not specified
    with ``histval``, they are taken to be at their steady state
    values. The result is a 1 by 2 vector, where the first entry
    stores the value of the planner objective when the initial
    Lagrange multipliers associated with the planner’s problem are set
    to their steady state values (see :comm:`ramsey_policy`).

    In contrast, the second entry stores the value of the planner
    objective with initial Lagrange multipliers of the planner’s
    problem set to 0, i.e. it is assumed that the planner exploits its
    ability to surprise private agents in the first period of
    implementing Ramsey policy. This is the value of implementating
    optimal policy for the first time and committing not to
    re-optimize in the future.

    Because it entails computing at least a second order
    approximation, this computation is skipped with a message when the
    model is too large (more than 180 state variables, including
    lagged Lagrange multipliers).

    *Steady state*

    See :comm:`Ramsey steady state <ramsey_model>`.


.. command:: discretionary_policy [VARIABLE_NAME...];
             discretionary_policy (OPTIONS...) [VARIABLE_NAME...];

    |br| This command computes an approximation of the optimal policy
    under discretion. The algorithm implemented is essentially an LQ
    solver, and is described by *Dennis (2007)*.

    You should ensure that your model is linear and your objective is
    quadratic. Also, you should set the ``linear`` option of the
    ``model`` block.

    *Options*

    This command accepts the same options than ``ramsey_policy``, plus:

    .. option:: discretionary_tol = NON-NEGATIVE DOUBLE

        Sets the tolerance level used to assess convergence of the
        solution algorithm. Default: ``1e-7``.

    .. option:: maxit = INTEGER

        Maximum number of iterations. Default: ``3000``.


.. command:: planner_objective MODEL_EXPRESSION ;

    |br| This command declares the policy maker objective, for use
    with ``ramsey_policy`` or ``discretionary_policy``.

    You need to give the one-period objective, not the discounted
    lifetime objective. The discount factor is given by the
    ``planner_discount`` option of ``ramsey_policy`` and
    ``discretionary_policy``. The objective function can only contain
    current endogenous variables and no exogenous ones. This
    limitation is easily circumvented by defining an appropriate
    auxiliary variable in the model.

    With ``ramsey_policy``, you are not limited to quadratic
    objectives: you can give any arbitrary nonlinear expression.

    With ``discretionary_policy``, the objective function must be quadratic.


Sensitivity and identification analysis
=======================================

Dynare provides an interface to the global sensitivity analysis (GSA)
toolbox (developed by the Joint Research Center (JRC) of the European
Commission), which is now part of the official Dynare
distribution. The GSA toolbox can be used to answer the following
questions:

    1. What is the domain of structural coefficients assuring the
       stability and determinacy of a DSGE model?
    2. Which parameters mostly drive the fit of, e.g., GDP and which
       the fit of inflation? Is there any conflict between the optimal
       fit of one observed series versus another?
    3. How to represent in a direct, albeit approximated, form the
       relationship between structural parameters and the reduced form
       of a rational expectations model?

The discussion of the methodologies and their application is described
in *Ratto (2008)*.

With respect to the previous version of the toolbox, in order to work
properly, the GSA toolbox no longer requires that the Dynare
estimation environment is set up.


Performing sensitivity analysis
-------------------------------

.. command:: dynare_sensitivity ;
             dynare_sensitivity(OPTIONS...);

    |br| This command triggers sensitivity analysis on a DSGE model.

    .. _sampl-opt:

    *Sampling Options*

    .. option:: Nsam = INTEGER

        Size of the Monte-Carlo sample. Default: ``2048``.

    .. option:: ilptau = INTEGER

        If equal to ``1``, use :math:`LP_\tau` quasi-Monte-Carlo. If
        equal to ``0``, use LHS Monte-Carlo. Default: ``1``.

    .. option:: pprior = INTEGER

        If equqal to ``1``, sample from the prior distributions. If
        equal to ``0``, sample from the multivariate normal
        :math:`N(\bar{\theta},\Sigma)`, where :math:`\bar{\theta}` is
        the posterior mode and :math:`\Sigma=H^{-1}`, :math:`H` is the
        Hessian at the mode. Default: ``1``.

    .. option:: prior_range = INTEGER

        If equal to ``1``, sample uniformly from prior ranges. If
        equal to ``0``, sample from prior distributions. Default:
        ``1``.

    .. option:: morris = INTEGER

        If equal to ``0``, ANOVA mapping (Type I error) If equal to
        ``1``, Screening analysis (Type II error). If equal to ``2``,
        Analytic derivatives (similar to Type II error, only valid
        when identification=1). Default: ``1`` when
        ``identification=1``, ``0`` otherwise.

    .. option:: morris_nliv = INTEGER

        Number of levels in Morris design. Default: ``6``.

    .. option:: morris_ntra = INTEGER

        Number trajectories in Morris design. Default: ``20``.

    .. option:: ppost = INTEGER

        If equal to ``1``, use Metropolis posterior sample. If equal
        to ``0``, do not use Metropolis posterior sample. Default:
        ``0``.

        NB: This overrides any other sampling option.

    .. option:: neighborhood_width = DOUBLE

        When ``pprior=0`` and ``ppost=0``, allows for the sampling of
        parameters around the value specified in the ``mode_file``, in
        the range :math:`\texttt{xparam1} \pm \left \vert
        \texttt{xparam1} \times \texttt{neighborhood\_width} \right
        \vert`. Default: ``0``.


    *Stability Mapping Options*

    .. option:: stab = INTEGER

        If equal to ``1``, perform stability mapping. If equal to
        ``0``, do not perform stability mapping. Default: ``1``.

    .. option:: load_stab = INTEGER

        If equal to ``1``, load a previously created sample. If equal
        to ``0``, generate a new sample. Default: ``0``.

    .. option:: alpha2_stab = DOUBLE

        Critical value for correlations :math:`\rho` in filtered
        samples: plot couples of parmaters with
        :math:`\left\vert\rho\right\vert>` ``alpha2_stab``. Default:
        ``0``.

    .. option:: pvalue_ks = DOUBLE

        The threshold :math:`pvalue` for significant
        Kolmogorov-Smirnov test (i.e. plot parameters with
        :math:`pvalue<` ``pvalue_ks``). Default: ``0.001``.

    .. option:: pvalue_corr = DOUBLE

        The threshold :math:`pvalue` for significant correlation in
        filtered samples (i.e. plot bivariate samples when
        :math:`pvalue<` ``pvalue_corr``). Default: ``1e-5``.


    *Reduced Form Mapping Options*

    .. option:: redform = INTEGER

        If equal to ``1``, prepare Monte-Carlo sample of reduced form
        matrices. If equal to ``0``, do not prepare Monte-Carlo sample
        of reduced form matrices. Default: ``0``.

    .. option:: load_redform = INTEGER

        If equal to ``1``, load previously estimated mapping. If equal
        to ``0``, estimate the mapping of the reduced form
        model. Default: ``0``.

    .. option:: logtrans_redform = INTEGER

        If equal to ``1``, use log-transformed entries. If equal to
        ``0``, use raw entries. Default: ``0``.

    .. option:: threshold_redform = [DOUBLE DOUBLE]

        The range over which the filtered Monte-Carlo entries of the
        reduced form coefficients should be analyzed. The first number
        is the lower bound and the second is the upper bound. An empty
        vector indicates that these entries will not be
        filtered. Default: empty.

    .. option:: ksstat_redform = DOUBLE

        Critical value for Smirnov statistics :math:`d` when reduced
        form entries are filtered. Default: ``0.001``.

    .. option:: alpha2_redform = DOUBLE

        Critical value for correlations :math:`\rho` when reduced form
        entries are filtered. Default: ``1e-5``.

    .. option:: namendo = (VARIABLE_NAME...)

        List of endogenous variables. ‘:’ indicates all endogenous
        variables. Default: empty.

    .. option:: namlagendo = (VARIABLE_NAME...)

        List of lagged endogenous variables. ‘:’ indicates all lagged
        endogenous variables. Analyze entries [namendo :math:`\times`
        namlagendo] Default: empty.

    .. option:: namexo = (VARIABLE_NAME...)

        List of exogenous variables. ‘:’ indicates all exogenous
        variables. Analyze entries [namendo :math:`\times`
        namexo]. Default: empty.


    *RMSE Options*

    .. option:: rmse = INTEGER

        If equal to ``1``, perform RMSE analysis. If equal to ``0``,
        do not perform RMSE analysis. Default: ``0``.

    .. option:: load_rmse = INTEGER

        If equal to ``1``, load previous RMSE analysis. If equal to
        ``0``, make a new RMSE analysis. Default: ``0``.

    .. option:: lik_only = INTEGER

        If equal to ``1``, compute only likelihood and posterior. If
        equal to ``0``, compute RMSE’s for all observed
        series. Default: ``0``.

    .. option:: var_rmse = (VARIABLE_NAME...)

        List of observed series to be considered. ‘:’ indicates all
        observed variables. Default: ``varobs``.

    .. option:: pfilt_rmse = DOUBLE

        Filtering threshold for RMSE’s. Default: ``0.1``.

    .. option:: istart_rmse = INTEGER

        Value at which to start computing RMSE’s (use ``2`` to avoid
        big intitial error). Default: ``presample+1``.

    .. option:: alpha_rmse = DOUBLE

        Critical value for Smirnov statistics :math:`d`: plot
        parameters with :math:`d>` ``alpha_rmse``. Default: ``0.001``.

    .. option:: alpha2_rmse = DOUBLE

        Critical value for correlation :math:`\rho`: plot couples of
        parmaters with :math:`\left\vert\rho\right\vert=`
        ``alpha2_rmse``. Default: ``1e-5``.

    .. option:: datafile = FILENAME

        See :ref:`datafile <dataf>`.

    .. option:: nobs = INTEGER
                nobs = [INTEGER1:INTEGER2]

        See :opt:`nobs <nobs = INTEGER>`.

    .. option:: first_obs = INTEGER

        See :opt:`first_obs <first_obs = INTEGER>`.

    .. option:: prefilter = INTEGER

        See :opt:`prefilter <prefilter = INTEGER>`.

    .. option:: presample = INTEGER

        See :opt:`presample <presample = INTEGER>`.

    .. option:: nograph

        See :opt:`nograph`.

    .. option:: nodisplay

        See :opt:`nodisplay`.

    .. option:: graph_format = FORMAT
                graph_format = ( FORMAT, FORMAT... )

        See :opt:`graph_format <graph_format = FORMAT>`.

    .. option:: conf_sig = DOUBLE

        See :ref:`conf_sig <confsig>`.

    .. option:: loglinear

        See :ref:`loglinear <logl>`.

    .. option:: mode_file = FILENAME

        See :opt:`mode_file <mode_file = FILENAME>`.

    .. option:: kalman_algo = INTEGER

        See :opt:`kalman_algo <kalman_algo = INTEGER>`.


    *Identification Analysis Options*

    .. option:: identification = INTEGER

        If equal to ``1``, performs identification analysis (forcing
        ``redform=0`` and ``morris=1``) If equal to ``0``, no
        identification analysis. Default: ``0``.

    .. option:: morris = INTEGER

        See :opt:`morris <morris = INTEGER>`.

    .. option:: morris_nliv = INTEGER

        See :opt:`morris_nliv <morris_nliv = INTEGER>`.

    .. option:: morris_ntra = INTEGER

        See :opt:`morris_ntra <morris_ntra = INTEGER>`.

    .. option:: load_ident_files = INTEGER

        Loads previously performed identification analysis. Default: ``0``.

    .. option:: useautocorr = INTEGER

        Use autocorrelation matrices in place of autocovariance
        matrices in moments for identification analysis. Default:
        ``0``.

    .. option:: ar = INTEGER

        Maximum number of lags for moments in identification
        analysis. Default: ``1``.

    .. option:: diffuse_filter = INTEGER

        See :opt:`diffuse_filter`.

.. _irf-momcal:

IRF/Moment calibration
----------------------

The ``irf_calibration`` and ``moment_calibration`` blocks allow
imposing implicit “endogenous” priors about IRFs and moments on the
model. The way it works internally is that any parameter draw that is
inconsistent with the “calibration” provided in these blocks is
discarded, i.e. assigned a prior density of ``0``. In the context of
``dynare_sensitivity``, these restrictions allow tracing out which
parameters are driving the model to satisfy or violate the given
restrictions.

IRF and moment calibration can be defined in ``irf_calibration`` and
``moment_calibration`` blocks:

.. block:: irf_calibration ;
           irf_calibration (OPTIONS...);

    |br| This block allows defining IRF calibration criteria and is
    terminated by ``end;``. To set IRF sign restrictions, the
    following syntax is used::

        VARIABLE_NAME(INTEGER),EXOGENOUS_NAME, -;
        VARIABLE_NAME(INTEGER:INTEGER),EXOGENOUS_NAME, +;

    To set IRF restrictions with specific intervals, the following
    syntax is used::

        VARIABLE_NAME(INTEGER),EXOGENOUS_NAME, [DOUBLE DOUBLE];
        VARIABLE_NAME(INTEGER:INTEGER),EXOGENOUS_NAME, [DOUBLE DOUBLE];

    When ``(INTEGER:INTEGER)`` is used, the restriction is considered
    to be fulfilled by a logical OR. A list of restrictions must
    always be fulfilled with logical AND.

    *Options*

    .. option:: relative_irf

        See :opt:`relative_irf`.

    *Example*

        ::

            irf_calibration;
            y(1:4), e_ys, [ -50 50]; //[first year response with logical OR]
            @#for ilag in 21:40
            R_obs(@{ilag}), e_ys, [0 6]; //[response from 5th to 10th years with logical AND]
            @#endfor
            end;


.. block:: moment_calibration ;
           moment_calibration (OPTIONS...);

    |br| This block allows defining moment calibration criteria. This
    block is terminated by ``end;``, and contains lines of the form::

        VARIABLE_NAME1,VARIABLE_NAME2(+/-INTEGER), [DOUBLE DOUBLE];
        VARIABLE_NAME1,VARIABLE_NAME2(+/-INTEGER), +/-;
        VARIABLE_NAME1,VARIABLE_NAME2(+/-(INTEGER:INTEGER)), [DOUBLE DOUBLE];
        VARIABLE_NAME1,VARIABLE_NAME2((-INTEGER:+INTEGER)), [DOUBLE DOUBLE];

    When ``(INTEGER:INTEGER)`` is used, the restriction is considered
    to be fulfilled by a logical OR. A list of restrictions must
    always be fulfilled with logical AND.

    *Example*

        ::

            moment_calibration;
            y_obs,y_obs, [0.5 1.5]; //[unconditional variance]
            y_obs,y_obs(-(1:4)), +; //[sign restriction for first year acf with logical OR]
            @#for ilag in -2:2
            y_obs,R_obs(@{ilag}), -; //[-2:2 ccf with logical AND]
            @#endfor
            @#for ilag in -4:4
            y_obs,pie_obs(@{ilag}), -; //[-4_4 ccf with logical AND]
            @#endfor
            end;


Performing identification analysis
----------------------------------

.. command:: identification ;
             identification (OPTIONS...);

    |br| This command triggers identification analysis.

    *Options*

    .. option:: ar = INTEGER

        Number of lags of computed autocorrelations (theoretical
        moments). Default: ``1``.

    .. option:: useautocorr = INTEGER

        If equal to ``1``, compute derivatives of autocorrelation. If
        equal to ``0``, compute derivatives of
        autocovariances. Default: ``0``.

    .. option:: load_ident_files = INTEGER

        If equal to ``1``, allow Dynare to load previously computed
        analyzes. Default: ``0``.

    .. option:: prior_mc = INTEGER

        Size of Monte-Carlo sample. Default: ``1``.

    .. option:: prior_range = INTEGER

        Triggers uniform sample within the range implied by the prior
        specifications (when ``prior_mc>1``). Default: ``0``.

    .. option:: advanced = INTEGER

        Shows a more detailed analysis, comprised of an analysis for
        the linearized rational expectation model as well as the
        associated reduced form solution. Further performs a brute
        force search of the groups of parameters best reproducing the
        behavior of each single parameter. The maximum dimension of
        the group searched is triggered by
        ``max_dim_cova_group``. Default: ``0``.

    .. option:: max_dim_cova_group = INTEGER

        In the brute force search (performed when ``advanced=1``) this
        option sets the maximum dimension of groups of parameters that
        best reproduce the behavior of each single model
        parameter. Default: ``2``.

    .. option:: periods = INTEGER

        When the analytic Hessian is not available (i.e. with missing
        values or diffuse Kalman filter or univariate Kalman filter),
        this triggers the length of stochastic simulation to compute
        Simulated Moments Uncertainty. Default: ``300``.

    .. option:: replic = INTEGER

        When the analytic Hessian is not available, this triggers the
        number of replicas to compute Simulated Moments
        Uncertainty. Default: ``100``.

    .. option:: gsa_sample_file = INTEGER

        If equal to ``0``, do not use sample file. If equal to ``1``,
        triggers gsa prior sample. If equal to ``2``, triggers gsa
        Monte-Carlo sample (i.e. loads a sample corresponding to
        ``pprior=0`` and ``ppost=0`` in the ``dynare_sensitivity``
        options). Default: ``0``.

    .. option:: gsa_sample_file = FILENAME

        Uses the provided path to a specific user defined sample
        file. Default: ``0``.

    .. option:: parameter_set = OPTION

        Specify the parameter set to use. Possible values for OPTION are:

            * ``calibration``
            * ``prior_mode``
            * ``prior_mean``
            * ``posterior_mode``
            * ``posterior_mean``
            * ``posterior_median``

        Default: ``prior_mean``.

    .. option:: lik_init = INTEGER

        See :opt:`lik_init <lik_init = INTEGER>`.

    .. option:: kalman_algo = INTEGER

        See :opt:`kalman_algo <kalman_algo = INTEGER>`.

    .. option:: nograph

        See :opt:`nograph`.

    .. option:: nodisplay

        See :opt:`nodisplay`.

    .. option:: graph_format = FORMAT
                graph_format = ( FORMAT, FORMAT... )

        See :opt:`graph_format <graph_format = FORMAT>`.


Types of analysis and output files
----------------------------------

The sensitivity analysis toolbox includes several types of
analyses. Sensitivity analysis results are saved locally in
``<mod_file>/gsa``, where ``<mod_file>.mod`` is the name of the DYNARE
model file.

Sampling
^^^^^^^^

The following binary files are produced:

    * ``<mod_file>_prior.mat``: this file stores information about the
      analyses performed sampling from the prior, i.e. ``pprior=1``
      and ``ppost=0``;
    * ``<mod_file>_mc.mat``: this file stores information about the
      analyses performed sampling from multivariate normal,
      i.e. ``pprior=0`` and ``ppost=0``;
    * ``<mod_file>_post.mat``: this file stores information about
      analyses performed using the Metropolis posterior sample,
      i.e. ``ppost=1``.


Stability Mapping
^^^^^^^^^^^^^^^^^

Figure files produced are of the form ``<mod_file>_prior_*.fig`` and
store results for stability mapping from prior Monte-Carlo samples:

    * ``<mod_file>_prior_stable.fig``: plots of the Smirnov test and
      the correlation analyses confronting the cdf of the sample
      fulfilling Blanchard-Kahn conditions (blue color) with the cdf
      of the rest of the sample (red color), i.e. either instability
      or indeterminacy or the solution could not be found (e.g. the
      steady state solution could not be found by the solver);
    * ``<mod_file>_prior_indeterm.fig``: plots of the Smirnov test and
      the correlation analyses confronting the cdf of the sample
      producing indeterminacy (red color) with the cdf of the rest of
      the sample (blue color);
    * ``<mod_file>_prior_unstable.fig``: plots of the Smirnov test and
      the correlation analyses confronting the cdf of the sample
      producing explosive roots (red color) with the cdf of the rest
      of the sample (blue color);
    * ``<mod_file>_prior_wrong.fig``: plots of the Smirnov test and
      the correlation analyses confronting the cdf of the sample where
      the solution could not be found (e.g. the steady state solution
      could not be found by the solver - red color) with the cdf of
      the rest of the sample (blue color);
    * ``<mod_file>_prior_calib.fig``: plots of the Smirnov test and
      the correlation analyses splitting the sample fulfilling
      Blanchard-Kahn conditions, by confronting the cdf of the sample
      where IRF/moment restrictions are matched (blue color) with the
      cdf where IRF/moment restrictions are NOT matched (red color);

Similar conventions apply for ``<mod_file>_mc_*.fig`` files, obtained
when samples from multivariate normal are used.


IRF/Moment restrictions
^^^^^^^^^^^^^^^^^^^^^^^

The following binary files are produced:

    * ``<mod_file>_prior_restrictions.mat``: this file stores
      information about the IRF/moment restriction analysis performed
      sampling from the prior ranges, i.e. ``pprior=1`` and
      ``ppost=0``;
    * ``<mod_file>_mc_restrictions.mat``: this file stores information
      about the IRF/moment restriction analysis performed sampling
      from multivariate normal, i.e. ``pprior=0`` and ``ppost=0``;
    * ``<mod_file>_post_restrictions.mat``: this file stores
      information about IRF/moment restriction analysis performed
      using the Metropolis posterior sample, i.e. ``ppost=1``.

Figure files produced are of the form
``<mod_file>_prior_irf_calib_*.fig`` and
``<mod_file>_prior_moment_calib_*.fig`` and store results for mapping
restrictions from prior Monte-Carlo samples:

    * ``<mod_file>_prior_irf_calib_<ENDO_NAME>_vs_<EXO_NAME>_<PERIOD>.fig``:
      plots of the Smirnov test and the correlation analyses splitting
      the sample fulfilling Blanchard-Kahn conditions, by confronting
      the cdf of the sample where the individual IRF restriction
      ``<ENDO_NAME>`` vs. ``<EXO_NAME>`` at period(s) ``<PERIOD>`` is
      matched (blue color) with the cdf where the IRF restriction is
      NOT matched (red color)
    * ``<mod_file>_prior_irf_calib_<ENDO_NAME>_vs_<EXO_NAME>_ALL.fig``:
      plots of the Smirnov test and the correlation analyses splitting
      the sample fulfilling Blanchard-Kahn conditions, by confronting
      the cdf of the sample where ALL the individual IRF restrictions
      for the same couple ``<ENDO_NAME>`` vs. ``<EXO_NAME>`` are
      matched (blue color) with the cdf where the IRF restriction is
      NOT matched (red color)
    * ``<mod_file>_prior_irf_restrictions.fig``: plots visual
      information on the IRF restrictions compared to the actual Monte
      Carlo realization from prior sample.
    * ``<mod_file>_prior_moment_calib_<ENDO_NAME1>_vs_<ENDO_NAME2>_<LAG>.fig``:
      plots of the Smirnov test and the correlation analyses splitting
      the sample fulfilling Blanchard-Kahn conditions, by confronting
      the cdf of the sample where the individual acf/ccf moment
      restriction ``<ENDO_NAME1>`` vs. ``<ENDO_NAME2>`` at lag(s)
      ``<LAG>`` is matched (blue color) with the cdf where the IRF
      restriction is NOT matched (red color)
    * ``<mod_file>_prior_moment_calib_<ENDO_NAME>_vs_<EXO_NAME>_ALL.fig``:
      plots of the Smirnov test and the correlation analyses splitting
      the sample fulfilling Blanchard-Kahn conditions, by confronting
      the cdf of the sample where ALL the individual acf/ccf moment
      restrictions for the same couple ``<ENDO_NAME1>``
      vs. ``<ENDO_NAME2>`` are matched (blue color) with the cdf where
      the IRF restriction is NOT matched (red color)
    * ``<mod_file>_prior_moment_restrictions.fig``: plots visual
      information on the moment restrictions compared to the actual
      Monte Carlo realization from prior sample.

Similar conventions apply for ``<mod_file>_mc_*.fig`` and
``<mod_file>_post_*.fig`` files, obtained when samples from
multivariate normal or from posterior are used.


Reduced Form Mapping
^^^^^^^^^^^^^^^^^^^^

When the option ``threshold_redform`` is not set, or it is empty (the
default), this analysis estimates a multivariate smoothing spline
ANOVA model (the ’mapping’) for the selected entries in the transition
matrix of the shock matrix of the reduce form first order solution of
the model. This mapping is done either with prior samples or with MC
samples with ``neighborhood_width``. Unless ``neighborhood_width`` is
set with MC samples, the mapping of the reduced form solution forces
the use of samples from prior ranges or prior distributions, i.e.:
``pprior=1`` and ``ppost=0``. It uses 250 samples to optimize
smoothing parameters and 1000 samples to compute the fit. The rest of
the sample is used for out-of-sample validation. One can also load a
previously estimated mapping with a new Monte-Carlo sample, to look at
the forecast for the new Monte-Carlo sample.

The following synthetic figures are produced:

    * ``<mod_file>_redform_<endo name>_vs_lags_*.fig``: shows bar
      charts of the sensitivity indices for the ten most important
      parameters driving the reduced form coefficients of the selected
      endogenous variables (``namendo``) versus lagged endogenous
      variables (``namlagendo``); suffix ``log`` indicates the results
      for log-transformed entries;
    * ``<mod_file>_redform_<endo name>_vs_shocks_*.fig``: shows bar
      charts of the sensitivity indices for the ten most important
      parameters driving the reduced form coefficients of the selected
      endogenous variables (``namendo``) versus exogenous variables
      (``namexo``); suffix ``log`` indicates the results for
      log-transformed entries;
    * ``<mod_file>_redform_gsa(_log).fig``: shows bar chart of all
      sensitivity indices for each parameter: this allows one to
      notice parameters that have a minor effect for any of the
      reduced form coefficients.

Detailed results of the analyses are shown in the subfolder
``<mod_file>/gsa/redform_prior`` for prior samples and in
``<mod_file>/gsa/redform_mc`` for MC samples with option
``neighborhood_width``, where the detailed results of the estimation
of the single functional relationships between parameters
:math:`\theta` and reduced form coefficient (denoted as :math:`y`
hereafter) are stored in separate directories named as:

    * ``<namendo>_vs_<namlagendo>``, for the entries of the transition matrix;
    * ``<namendo>_vs_<namexo>``, for entries of the matrix of the shocks.

The following files are stored in each directory (we stick with prior
sample but similar conventions are used for MC samples):

    * ``<mod_file>_prior_<namendo>_vs_<namexo>.fig``: histogram and
      CDF plot of the MC sample of the individual entry of the shock
      matrix, in sample and out of sample fit of the ANOVA model;
    * ``<mod_file>_prior_<namendo>_vs_<namexo>_map_SE.fig``: for
      entries of the shock matrix it shows graphs of the estimated
      first order ANOVA terms :math:`y = f(\theta_i)` for each deep
      parameter :math:`\theta_i`;
    * ``<mod_file>_prior_<namendo>_vs_<namlagendo>.fig``: histogram
      and CDF plot of the MC sample of the individual entry of the
      transition matrix, in sample and out of sample fit of the ANOVA
      model;
    * ``<mod_file>_prior_<namendo>_vs_<namlagendo>_map_SE.fig``: for
      entries of the transition matrix it shows graphs of the
      estimated first order ANOVA terms :math:`y = f(\theta_i)` for
      each deep parameter :math:`\theta_i`;
    * ``<mod_file>_prior_<namendo>_vs_<namexo>_map.mat``,
      ``<mod_file>_<namendo>_vs_<namlagendo>_map.mat``: these files
      store info in the estimation;

When option ``logtrans_redform`` is set, the ANOVA estimation is
performed using a log-transformation of each y. The ANOVA mapping is
then transformed back onto the original scale, to allow comparability
with the baseline estimation. Graphs for this log-transformed case,
are stored in the same folder in files denoted with the ``_log``
suffix.

When the option ``threshold_redform`` is set, the analysis is
performed via Monte Carlo filtering, by displaying parameters that
drive the individual entry ``y`` inside the range specified in
``threshold_redform``. If no entry is found (or all entries are in the
range), the MCF algorithm ignores the range specified in
``threshold_redform`` and performs the analysis splitting the MC
sample of ``y`` into deciles. Setting ``threshold_redform=[-inf inf]``
triggers this approach for all ``y``’s.

Results are stored in subdirectories of ``<mod_file>/gsa/redform_prior`` named

    * ``<mod_file>_prior_<namendo>_vs_<namlagendo>_threshold``, for
      the entries of the transition matrix;
    * ``<mod_file>_prior_<namendo>_vs_<namexo>_threshold``, for
      entries of the matrix of the shocks.

The files saved are named:

    * ``<mod_file>_prior_<namendo>_vs_<namexo>_threshold.fig``, ``<mod_file>_<namendo>_vs_<namlagendo>_threshold.fig``: graphical outputs;
    * ``<mod_file>_prior_<namendo>_vs_<namexo>_threshold.mat``, ``<mod_file>_<namendo>_vs_<namlagendo>_threshold.mat``: info on the analysis;


RMSE
^^^^

The RMSE analysis can be performed with different types of sampling options:

    1. When ``pprior=1`` and ``ppost=0``, the toolbox analyzes the
       RMSEs for the Monte-Carlo sample obtained by sampling
       parameters from their prior distributions (or prior ranges):
       this analysis provides some hints about what parameter drives
       the fit of which observed series, prior to the full estimation;
    2. When ``pprior=0`` and ``ppost=0``, the toolbox analyzes the
       RMSEs for a multivariate normal Monte-Carlo sample, with
       covariance matrix based on the inverse Hessian at the optimum:
       this analysis is useful when maximum likelihood estimation is
       done (i.e. no Bayesian estimation);
    3. When ``ppost=1`` the toolbox analyzes the RMSEs for the
       posterior sample obtained by Dynare’s Metropolis procedure.

The use of cases 2 and 3 requires an estimation step beforehand. To
facilitate the sensitivity analysis after estimation, the
``dynare_sensitivity`` command also allows you to indicate some
options of the ``estimation command``. These are:

    * ``datafile``
    * ``nobs``
    * ``first_obs``
    * ``prefilter``
    * ``presample``
    * ``nograph``
    * ``nodisplay``
    * ``graph_format``
    * ``conf_sig``
    * ``loglinear``
    * ``mode_file``

Binary files produced my RMSE analysis are:

    * ``<mod_file>_prior_*.mat``: these files store the filtered and
      smoothed variables for the prior Monte-Carlo sample, generated
      when doing RMSE analysis (``pprior=1`` and ``ppost=0``);
    * ``<mode_file>_mc_*.mat``: these files store the filtered and
      smoothed variables for the multivariate normal Monte-Carlo
      sample, generated when doing RMSE analysis (``pprior=0`` and
      ``ppost=0``).

Figure files <mod_file>_rmse_*.fig store results for the RMSE analysis.

    * ``<mod_file>_rmse_prior*.fig``: save results for the analysis
      using prior Monte-Carlo samples;
    * ``<mod_file>_rmse_mc*.fig``: save results for the analysis using
      multivariate normal Monte-Carlo samples;
    * ``<mod_file>_rmse_post*.fig``: save results for the analysis
      using Metropolis posterior samples.

The following types of figures are saved (we show prior sample to fix
ideas, but the same conventions are used for multivariate normal and
posterior):

    * ``<mod_file>_rmse_prior_params_*.fig``: for each parameter,
      plots the cdfs corresponding to the best 10% RMSEs of each
      observed series (only those cdfs below the significance
      threshold ``alpha_rmse``);
    * ``<mod_file>_rmse_prior_<var_obs>_*.fig``: if a parameter
      significantly affects the fit of ``var_obs``, all possible
      trade-off’s with other observables for same parameter are
      plotted;
    * ``<mod_file>_rmse_prior_<var_obs>_map.fig``: plots the MCF
      analysis of parameters significantly driving the fit the
      observed series ``var_obs``;
    * ``<mod_file>_rmse_prior_lnlik*.fig``: for each observed series,
      plots in BLUE the cdf of the log-likelihood corresponding to the
      best 10% RMSEs, in RED the cdf of the rest of the sample and in
      BLACK the cdf of the full sample; this allows one to see the
      presence of some idiosyncratic behavior;
    * ``<mod_file>_rmse_prior_lnpost*.fig``: for each observed series,
      plots in BLUE the cdf of the log-posterior corresponding to the
      best 10% RMSEs, in RED the cdf of the rest of the sample and in
      BLACK the cdf of the full sample; this allows one to see
      idiosyncratic behavior;
    * ``<mod_file>_rmse_prior_lnprior*.fig``: for each observed
      series, plots in BLUE the cdf of the log-prior corresponding to
      the best 10% RMSEs, in RED the cdf of the rest of the sample and
      in BLACK the cdf of the full sample; this allows one to see
      idiosyncratic behavior;
    * ``<mod_file>_rmse_prior_lik.fig``: when ``lik_only=1``, this
      shows the MCF tests for the filtering of the best 10%
      log-likelihood values;
    * ``<mod_file>_rmse_prior_post.fig``: when ``lik_only=1``, this
      shows the MCF tests for the filtering of the best 10%
      log-posterior values.


Screening Analysis
^^^^^^^^^^^^^^^^^^

Screening analysis does not require any additional options with
respect to those listed in :ref:`Sampling Options <sampl-opt>`. The
toolbox performs all the analyses required and displays results.

The results of the screening analysis with Morris sampling design are
stored in the subfolder ``<mod_file>/gsa/screen``. The data file
``<mod_file>_prior`` stores all the information of the analysis
(Morris sample, reduced form coefficients, etc.).

Screening analysis merely concerns reduced form coefficients. Similar
synthetic bar charts as for the reduced form analysis with Monte-Carlo
samples are saved:

    * ``<mod_file>_redform_<endo name>_vs_lags_*.fig``: shows bar
      charts of the elementary effect tests for the ten most important
      parameters driving the reduced form coefficients of the selected
      endogenous variables (``namendo``) versus lagged endogenous
      variables (``namlagendo``);
    * ``<mod_file>_redform_<endo name>_vs_shocks_*.fig``: shows bar
      charts of the elementary effect tests for the ten most important
      parameters driving the reduced form coefficients of the selected
      endogenous variables (``namendo``) versus exogenous variables
      (``namexo``);
    * ``<mod_file>_redform_screen.fig``: shows bar chart of all
      elementary effect tests for each parameter: this allows one to
      identify parameters that have a minor effect for any of the
      reduced form coefficients.


Identification Analysis
^^^^^^^^^^^^^^^^^^^^^^^

Setting the option ``identification=1``, an identification analysis
based on theoretical moments is performed. Sensitivity plots are
provided that allow to infer which parameters are most likely to be
less identifiable.

Prerequisite for properly running all the identification routines, is
the keyword ``identification``; in the Dynare model file. This keyword
triggers the computation of analytic derivatives of the model with
respect to estimated parameters and shocks. This is required for
option ``morris=2``, which implements *Iskrev (2010)* identification
analysis.

For example, the placing::

    identification;
    dynare_sensitivity(identification=1, morris=2);

in the Dynare model file triggers identification analysis using
analytic derivatives *Iskrev (2010)*, jointly with the mapping of the
acceptable region.

The identification analysis with derivatives can also be triggered by
the commands ``identification;`` This does not do the mapping of
acceptable regions for the model and uses the standard random sampler
of Dynare. It completely offsets any use of the sensitivity analysis
toolbox.



Markov-switching SBVAR
======================

Given a list of variables, observed variables and a data file, Dynare
can be used to solve a Markov-switching SBVAR model according to
*Sims, Waggoner and Zha (2008)* [#f10]_ . Having done this, you can
create forecasts and compute the marginal data density, regime
probabilities, IRFs, and variance decomposition of the model.

The commands have been modularized, allowing for multiple calls to the
same command within a ``<mod_file>.mod`` file. The default is to use
``<mod_file>`` to tag the input (output) files used (produced) by the
program. Thus, to call any command more than once within a
``<mod_file>.mod`` file, you must use the ``*_tag`` options described
below.


.. command:: markov_switching (OPTIONS...);

    |br| Declares the Markov state variable information of a
    Markov-switching SBVAR model.

    *Options*

    .. option:: chain = INTEGER

        The Markov chain considered. Default: ``none``.

    .. option:: number_of_regimes = INTEGER

        Specifies the total number of regimes in the Markov
        Chain. This is a required option.

    .. option:: duration = DOUBLE | [ROW VECTOR OF DOUBLES]

        The duration of the regimes or regimes. This is a required
        option. When passed a scalar real number, it specifies the
        average duration for all regimes in this chain. When passed a
        vector of size equal ``number_of_regimes``, it specifies the
        average duration of the associated regimes
        (``1:number_of_regimes``) in this chain. An absorbing state
        can be specified through the :opt:`restrictions <restrictions
        = [[ROW VECTOR OF 3 DOUBLES],[ROW VECTOR OF 3 DOUBLES],...]>`
        option.

    .. option:: restrictions = [[ROW VECTOR OF 3 DOUBLES],[ROW VECTOR OF 3 DOUBLES],...]

        Provides restrictions on this chain’s regime transition
        matrix. Its vector argument takes three inputs of the form:
        ``[current_period_regime, next_period_regime,
        transition_probability]``.

        The first two entries are positive integers, and the third is
        a non-negative real in the set [0,1]. If restrictions are
        specified for every transition for a regime, the sum of the
        probabilities must be 1. Otherwise, if restrictions are not
        provided for every transition for a given regime the sum of
        the provided transition probabilities msut be <1. Regardless
        of the number of lags, the restrictions are specified for
        parameters at time ``t`` since the transition probability for
        a parameter at t is equal to that of the parameter at ``t-1``.

    In case of estimating a MS-DSGE model, [#f11]_ in addition the
    following options are allowed:

    .. option:: parameters = [LIST OF PARAMETERS]

        This option specifies which parameters are controlled by this
        Markov Chain.

    .. option:: number_of_lags = DOUBLE

        Provides the number of lags that each parameter can take
        within each regime in this chain.

    *Example*

        ::

            markov_switching(chain=1, duration=2.5, restrictions=[[1,3,0],[3,1,0]]);

        Specifies a Markov-switching BVAR with a first chain with 3
        regimes that all have a duration of 2.5 periods. The
        probability of directly going from regime 1 to regime 3 and
        vice versa is 0.

    *Example*

        ::

            markov_switching(chain=2, number_of_regimes=3, duration=[0.5, 2.5, 2.5],
            parameter=[alpha, rho], number_of_lags=2, restrictions=[[1,3,0],[3,3,1]]);

        Specifies a Markov-switching DSGE model with a second chain
        with 3 regimes that have durations of 0.5, 2.5, and 2.5
        periods, respectively. The switching parameters are ``alpha``
        and ``rho``. The probability of directly going from regime 1
        to regime 3 is 0, while regime 3 is an absorbing state.


.. command:: svar (OPTIONS...);

    |br| Each Markov chain can control the switching of a set of
    parameters. We allow the parameters to be divided equation by
    equation and by variance or slope and intercept.

    *Options*

    .. option:: coefficients

        Specifies that only the slope and intercept in the given
        equations are controlled by the given chain. One, but not
        both, of ``coefficients`` or ``variances`` must
        appear. Default: ``none``.

    .. option:: variances

        Specifies that only variances in the given equations are
        controlled by the given chain. One, but not both, of
        ``coefficients`` or ``variances`` must appear. Default:
        ``none``.

    .. option:: equations

        Defines the equation controlled by the given chain. If not
        specified, then all equations are controlled by
        ``chain``. Default: ``none``.

    .. option:: chain = INTEGER

        Specifies a Markov chain defined by
        :comm:`markov_switching`. Default: ``none``.


.. command:: sbvar (OPTIONS...);

    |br| To be documented. For now, see the wiki:
    `<https://www.dynare.org/DynareWiki/SbvarOptions>`_

    *Options*

    ``datafile``,
    ``freq``,
    ``initial_year``,
    ``initial_subperiod``,
    ``final_year``,
    ``final_subperiod``,
    ``data``,
    ``vlist``,
    ``vlistlog``,
    ``vlistper``,
    ``restriction_fname``,
    ``nlags``,
    ``cross_restrictions``,
    ``contemp_reduced_form``,
    ``real_pseudo_forecast``,
    ``no_bayesian_prior``,
    ``dummy_obs``,
    ``nstates``,
    ``indxscalesstates``,
    ``alpha``,
    ``beta``,
    ``gsig2_lmdm``,
    ``q_diag``,
    ``flat_prior``,
    ``ncsk``,
    ``nstd``,
    ``ninv``,
    ``indxparr``,
    ``indxovr``,
    ``aband``,
    ``indxap``,
    ``apband``,
    ``indximf``,
    ``indxfore``,
    ``foreband``,
    ``indxgforhat``,
    ``indxgimfhat``,
    ``indxestima``,
    ``indxgdls``,
    ``eq_ms``,
    ``cms``,
    ``ncms``,
    ``eq_cms``,
    ``tlindx``,
    ``tlnumber``,
    ``cnum``,
    ``forecast``,
    ``coefficients_prior_hyperparameters``

.. block:: svar_identification ;

    |br| This block is terminated by ``end;`` and contains lines of the form::

        UPPER_CHOLESKY;
        LOWER_CHOLESKY;
        EXCLUSION CONSTANTS;
        EXCLUSION LAG INTEGER; VARIABLE_NAME [,VARIABLE_NAME...];
        EXCLUSION LAG INTEGER; EQUATION INTEGER, VARIABLE_NAME [,VARIABLE_NAME...];
        RESTRICTION EQUATION INTEGER, EXPRESSION = EXPRESSION;

    To be documented. For now, see the wiki:
    `<http://www.dynare.org/DynareWiki/MarkovSwitchingInterface>`_


.. command:: ms_estimation (OPTIONS...);

    |br| Triggers the creation of an initialization file for, and the
    estimation of, a Markov-switching SBVAR model. At the end of the
    run, the :math:`A^0`, :math:`A^+`, :math:`Q` and :math:`\zeta`
    matrices are contained in the ``oo_.ms`` structure.

    *General Options*

    .. option:: file_tag = FILENAME

        The portion of the filename associated with this run. This
        will create the model initialization file,
        ``init_<file_tag>.dat``. Default: ``<mod_file>``.

    .. option:: output_file_tag = FILENAME

        The portion of the output filename that will be assigned to
        this run. This will create, among other files,
        ``est_final_<output_file_tag>.out``,
        ``est_intermediate_<output_file_tag>.out``. Default:
        ``<file_tag>``.

    .. option:: no_create_init

        Do not create an initialization file for the model. Passing
        this option will cause the *Initialization Options* to be
        ignored. Further, the model will be generated from the output
        files associated with the previous estimation run
        (i.e. ``est_final_<file_tag>.out``,
        ``est_intermediate_<file_tag>.out`` or
        ``init_<file_tag>.dat``, searched for in sequential
        order). This functionality can be useful for continuing a
        previous estimation run to ensure convergence was reached or
        for reusing an initialization file. NB: If this option is not
        passed, the files from the previous estimation run will be
        overwritten. Default: off (i.e. create initialization file)

    *Initialization Options*

    .. option:: coefficients_prior_hyperparameters = [DOUBLE1 DOUBLE2 ... DOUBLE6]

        Sets the hyper parameters for the model. The six elements of
        the argument vector have the following interpretations:

        ``1``

            Overall tightness for :math:`A^0` and :math:`A^+`.

        ``2``

            Relative tightness for :math:`A^+`.

        ``3``

            Relative tightness for the constant term.

        ``4``

            Tightness on lag decay (range: 1.2 - 1.5); a faster decay
            produces better inflation process.

        ``5``

            Weight on nvar sums of coeffs dummy observations (unit roots).

        ``6``

            Weight on single dummy initial observation including constant.

        Default: ``[1.0 1.0 0.1 1.2 1.0 1.0]``

    .. option:: freq = INTEGER | monthly | quarterly | yearly

        Frequency of the data (e.g. ``monthly, 12``). Default: ``4``.

    .. option:: initial_year = INTEGER

        The first year of data. Default: ``none``.

    .. option:: initial_subperiod = INTEGER

        The first period of data (i.e. for quarterly data, an integer
        in ``[1,4]``). Default: ``1``.

    .. option:: final_year = INTEGER

        The last year of data. Default: Set to encompass entire dataset.

    .. option:: final_subperiod = INTEGER

        The final period of data (i.e. for monthly data, an integer in
        ``[1,12]``. Default: When final_year is also missing, set to
        encompass entire dataset; when ``final_year`` is indicated,
        set to the maximum number of subperiods given the frequency
        (i.e. 4 for quarterly data, 12 for monthly,...).

    .. option:: datafile = FILENAME

        See :ref:`datafile <dataf>`.

    .. option:: xls_sheet = NAME

        See :opt:`xls_sheet <xls_sheet = NAME>`.

    .. option:: xls_range = RANGE

        See :opt:`xls_range <xls_range = RANGE>`.

    .. option:: nlags = INTEGER

        The number of lags in the model. Default: ``1``.

    .. option:: cross_restrictions

        Use cross :math:`A^0` and :math:`A^+` restrictions. Default: ``off``.

    .. option:: contemp_reduced_form

        Use contemporaneous recursive reduced form. Default: ``off``.

    .. option:: no_bayesian_prior

        Do not use Bayesian prior. Default: ``off`` (i.e. use Bayesian prior).

    .. option:: alpha = INTEGER

        Alpha value for squared time-varying structural shock
        lambda. Default: ``1``.

    .. option:: beta = INTEGER

        Beta value for squared time-varying structural shock
        lambda. Default: ``1``.

    .. option:: gsig2_lmdm = INTEGER

        The variance for each independent :math:`\lambda` parameter
        under ``SimsZha`` restrictions. Default: ``50^2``.

    .. option:: specification = sims_zha | none

        This controls how restrictions are imposed to reduce the
        number of parameters. Default: ``Random Walk``.

    *Estimation Options*

    .. option:: convergence_starting_value = DOUBLE

        This is the tolerance criterion for convergence and refers to
        changes in the objective function value. It should be rather
        loose since it will gradually be tightened during
        estimation. Default: ``1e-3``.

    .. option:: convergence_ending_value = DOUBLE

        The convergence criterion ending value. Values much smaller
        than square root machine epsilon are probably
        overkill. Default: ``1e-6``.

    .. option:: convergence_increment_value = DOUBLE

        Determines how quickly the convergence criterion moves from
        the starting value to the ending value. Default: ``0.1``.

    .. option:: max_iterations_starting_value = INTEGER

        This is the maximum number of iterations allowed in the
        hill-climbing optimization routine and should be rather small
        since it will gradually be increased during
        estimation. Default: ``50``.

    .. option:: max_iterations_increment_value = DOUBLE

        Determines how quickly the maximum number of iterations is
        increased. Default: ``2``.

    .. option:: max_block_iterations = INTEGER

        The parameters are divided into blocks and optimization
        proceeds over each block. After a set of blockwise
        optimizations are performed, the convergence criterion is
        checked and the blockwise optimizations are repeated if the
        criterion is violated. This controls the maximum number of
        times the blockwise optimization can be performed. Note that
        after the blockwise optimizations have converged, a single
        optimization over all the parameters is performed before
        updating the convergence value and maximum number of
        iterations. Default: ``100``.

    .. option:: max_repeated_optimization_runs = INTEGER

        The entire process described by :opt:`max_block_iterations
        <max_block_iterations = INTEGER>` is repeated until
        improvement has stopped. This is the maximum number of times
        the process is allowed to repeat. Set this to ``0`` to not
        allow repetitions. Default: ``10``.

    .. option:: function_convergence_criterion = DOUBLE

        The convergence criterion for the objective function when
        ``max_repeated_optimizations_runs`` is positive. Default:
        ``0.1``.

    .. option:: parameter_convergence_criterion = DOUBLE

        The convergence criterion for parameter values when
        ``max_repeated_optimizations_runs`` is positive. Default:
        ``0.1``.

    .. option:: number_of_large_perturbations = INTEGER

        The entire process described by :opt:`max_block_iterations
        <max_block_iterations = INTEGER>` is repeated with random
        starting values drawn from the posterior. This specifies the
        number of random starting values used. Set this to ``0`` to
        not use random starting values. A larger number should be
        specified to ensure that the entire parameter space has been
        covered. Default: ``5``.

    .. option:: number_of_small_perturbations = INTEGER

        The number of small perturbations to make after the large
        perturbations have stopped improving. Setting this number much
        above ``10`` is probably overkill. Default: ``5``.

    .. option:: number_of_posterior_draws_after_perturbation = INTEGER

        The number of consecutive posterior draws to make when
        producing a small perturbation. Because the posterior draws
        are serially correlated, a small number will result in a small
        perturbation. Default: ``1``.

    .. option:: max_number_of_stages = INTEGER

        The small and large perturbation are repeated until
        improvement has stopped. This specifies the maximum number of
        stages allowed. Default: ``20``.

    .. option:: random_function_convergence_criterion = DOUBLE

        The convergence criterion for the objective function when
        ``number_of_large_perturbations`` is positive. Default:
        ``0.1``.

    .. option:: random_parameter_convergence_criterion = DOUBLE

        The convergence criterion for parameter values when
        ``number_of_large_perturbations`` is positive. Default:
        ``0.1``.

    *Example*

        ::

            ms_estimation(datafile=data, initial_year=1959, final_year=2005,
            nlags=4, max_repeated_optimization_runs=1, max_number_of_stages=0);

            ms_estimation(file_tag=second_run, datafile=data, initial_year=1959,
            final_year=2005, nlags=4, max_repeated_optimization_runs=1,
            max_number_of_stages=0);

            ms_estimation(file_tag=second_run, output_file_tag=third_run,
            no_create_init, max_repeated_optimization_runs=5,
            number_of_large_perturbations=10);


.. command:: ms_simulation ;
             ms_simulation (OPTIONS...);

    |br| Simulates a Markov-switching SBVAR model.

    *Options*

    .. option:: file_tag = FILENAME

        The portion of the filename associated with the
        ``ms_estimation`` run. Default: ``<mod_file>``.

    .. option:: output_file_tag = FILENAME

        The portion of the output filename that will be assigned to
        this run. Default: ``<file_tag>``.

    .. option:: mh_replic = INTEGER

        The number of draws to save. Default: ``10,000``.

    .. option:: drop = INTEGER

        The number of burn-in draws. Default: ``0.1*mh_replic*thinning_factor``.

    .. option:: thinning_factor = INTEGER

        The total number of draws is equal to
        ``thinning_factor*mh_replic+drop``. Default: ``1``.

    .. option:: adaptive_mh_draws = INTEGER

        Tuning period for Metropolis-Hastings draws. Default: ``30,000``.

    .. option:: save_draws

        Save all elements of :math:`A^0`, :math:`A^+`, :math:`Q`, and
        :math:`\zeta`, to a file named ``draws_<<file_tag>>.out`` with
        each draw on a separate line. A file that describes how these
        matrices are laid out is contained in
        ``draws_header_<<file_tag>>.out``. A file called
        ``load_flat_file.m`` is provided to simplify loading the saved
        files into the corresponding variables ``A0``, ``Aplus``,
        ``Q``, and ``Zeta`` in your MATLAB/Octave workspace. Default:
        ``off``.

    *Example*

        ::

            ms_simulation(file_tag=second_run);
            ms_simulation(file_tag=third_run, mh_replic=5000, thinning_factor=3);


.. command:: ms_compute_mdd ;
             ms_compute_mdd (OPTIONS...);

    |br| Computes the marginal data density of a Markov-switching
    SBVAR model from the posterior draws. At the end of the run, the
    Muller and Bridged log marginal densities are contained in the
    ``oo_.ms`` structure.

    *Options*

    .. option:: file_tag = FILENAME

        See :opt:`file_tag <file_tag = FILENAME>`.

    .. option:: output_file_tag = FILENAME

        See :opt:`output_file_tag <output_file_tag = FILENAME>`.

    .. option:: simulation_file_tag = FILENAME

        The portion of the filename associated with the simulation
        run. Default: ``<file_tag>``.

    .. option:: proposal_type = INTEGER

        The proposal type:

        ``1``

            Gaussian.

        ``2``

            Power.

        ``3``

            Truncated Power.

        ``4``

            Step.

        ``5``

            Truncated Gaussian.

        Default: ``3``

    .. option:: proposal_lower_bound = DOUBLE

        The lower cutoff in terms of probability. Not used for
        ``proposal_type`` in ``[1,2]``. Required for all other
        proposal types. Default: ``0.1``.

    .. option:: proposal_upper_bound = DOUBLE

        The upper cutoff in terms of probability. Not used for
        ``proposal_type`` equal to ``1``. Required for all other
        proposal types. Default: ``0.9``.

    .. option:: mdd_proposal_draws = INTEGER

        The number of proposal draws. Default: ``100,000``.

    .. option:: mdd_use_mean_center

        Use the posterior mean as center. Default: ``off``.


.. command:: ms_compute_probabilities ;
             ms_compute_probabilities (OPTIONS...);

    |br| Computes smoothed regime probabilities of a Markov-switching SBVAR
    model. Output ``.eps`` files are contained in
    ``<output_file_tag/Output/Probabilities>``.

    *Options*

    .. option:: file_tag = FILENAME

        See :opt:`file_tag <file_tag = FILENAME>`.

    .. option:: output_file_tag = FILENAME

        See :opt:`output_file_tag <output_file_tag = FILENAME>`.

    .. option:: filtered_probabilities

        Filtered probabilities are computed instead of
        smoothed. Default: ``off``.

    .. option:: real_time_smoothed

        Smoothed probabilities are computed based on time ``t``
        information for :math:`0\le t\le nobs`. Default: ``off``


.. command:: ms_irf ;
             ms_irf (OPTIONS...);

    |br| Computes impulse response functions for a Markov-switching SBVAR
    model. Output ``.eps`` files are contained in
    ``<output_file_tag/Output/IRF>``, while data files are contained
    in ``<output_file_tag/IRF>``.

    *Options*

    .. option:: file_tag = FILENAME

        See :opt:`file_tag <file_tag = FILENAME>`.

    .. option:: output_file_tag = FILENAME

        See :opt:`output_file_tag <output_file_tag = FILENAME>`.

    .. option:: simulation_file_tag = FILENAME

        See :opt:`simulation_file_tag <simulation_file_tag = FILENAME>`.

    .. option:: horizon = INTEGER

        The forecast horizon. Default: ``12``.

    .. option:: filtered_probabilities

        Uses filtered probabilities at the end of the sample as
        initial conditions for regime probabilities. Only one of
        ``filtered_probabilities``, ``regime`` and ``regimes`` may be
        passed. Default: ``off``.

    .. option:: error_band_percentiles = [DOUBLE1 ...]

        The percentiles to compute. Default: ``[0.16 0.50 0.84]``. If
        ``median`` is passed, the default is ``[0.5]``.

    .. option:: shock_draws = INTEGER

        The number of regime paths to draw. Default: ``10,000``.

    .. option:: shocks_per_parameter = INTEGER

        The number of regime paths to draw under parameter
        uncertainty. Default: ``10``.

    .. option:: thinning_factor = INTEGER

        Only :math:`1/ \texttt{thinning\_factor}` of the draws in
        posterior draws file are used. Default: ``1``.

    .. option:: free_parameters = NUMERICAL_VECTOR

        A vector of free parameters to initialize theta of the
        model. Default: use estimated parameters

    .. option:: parameter_uncertainty

        Calculate IRFs under parameter uncertainty. Requires that
        ``ms_simulation`` has been run. Default: ``off``.

    .. option:: regime = INTEGER

        Given the data and model parameters, what is the ergodic
        probability of being in the specified regime. Only one of
        ``filtered_probabilities``, ``regime`` and ``regimes`` may be
        passed. Default: ``off``.

    .. option:: regimes

        Describes the evolution of regimes. Only one of
        ``filtered_probabilities``, ``regime`` and ``regimes`` may be
        passed. Default: ``off``.

    .. option:: median

        A shortcut to setting
        ``error_band_percentiles=[0.5]``. Default: ``off``.


.. command:: ms_forecast ;
             ms_forecast (OPTIONS...);

    |br| Generates forecasts for a Markov-switching SBVAR
    model. Output ``.eps`` files are contained in
    ``<output_file_tag/Output/Forecast>``, while data files are
    contained in ``<output_file_tag/Forecast>``.

    *Options*

    .. option:: file_tag = FILENAME

        See :opt:`file_tag <file_tag = FILENAME>`.

    .. option:: output_file_tag = FILENAME

        See :opt:`output_file_tag <output_file_tag = FILENAME>`.

    .. option:: simulation_file_tag = FILENAME

        See :opt:`simulation_file_tag <simulation_file_tag = FILENAME>`.

    .. option:: data_obs_nbr = INTEGER

        The number of data points included in the output. Default: ``0``.

    .. option:: error_band_percentiles = [DOUBLE1 ...]

        See :opt:`error_band_percentiles <error_band_percentiles =
        [DOUBLE1 ...]>`.

    .. option:: shock_draws = INTEGER

        See :opt:`shock_draws <shock_draws = INTEGER>`.

    .. option:: shocks_per_parameter = INTEGER

        See :opt:`shocks_per_parameter <shocks_per_parameter = INTEGER>`.

    .. option:: thinning_factor = INTEGER

        See :opt:`thinning_factor <thinning_factor = INTEGER>`.

    .. option:: free_parameters = NUMERICAL_VECTOR

        See :opt:`free_parameters <free_parameters = NUMERICAL_VECTOR>`.

    .. option:: parameter_uncertainty

        See :opt:`parameter_uncertainty`.

    .. option:: regime = INTEGER

        See :opt:`regime <regime = INTEGER>`.

    .. option:: regimes

        See :opt:`regimes`.

    .. option:: median

        See :opt:`median`.

    .. option:: horizon = INTEGER

        See :opt:`horizon <horizon = INTEGER>`.


.. command:: ms_variance_decomposition ;
             ms_variance_decomposition (OPTIONS...);

    |br| Computes the variance decomposition for a Markov-switching
    SBVAR model. Output ``.eps`` files are contained in
    ``<output_file_tag/Output/Variance_Decomposition>``, while data
    files are contained in
    ``<output_file_tag/Variance_Decomposition>``.

    *Options*

    .. option:: file_tag = FILENAME

        See :opt:`file_tag <file_tag = FILENAME>`.

    .. option:: output_file_tag = FILENAME

        See :opt:`output_file_tag <output_file_tag = FILENAME>`.

    .. option:: simulation_file_tag = FILENAME

        See :opt:`simulation_file_tag <simulation_file_tag = FILENAME>`.

    .. option:: horizon = INTEGER

        See :opt:`horizon <horizon = INTEGER>`.

    .. option:: filtered_probabilities

        See :opt:`filtered_probabilities`.

    .. option:: no_error_bands

        Do not output percentile error bands (i.e. compute
        mean). Default: ``off`` (i.e. output error bands)

    .. option:: error_band_percentiles = [DOUBLE1 ...]

        See :opt:`error_band_percentiles <error_band_percentiles =
        [DOUBLE1 ...]>`.

    .. option:: shock_draws = INTEGER

        See :opt:`shock_draws <shock_draws = INTEGER>`.

    .. option:: shocks_per_parameter = INTEGER

        See :opt:`shocks_per_parameter <shocks_per_parameter = INTEGER>`.

    .. option:: thinning_factor = INTEGER

        See :opt:`thinning_factor <thinning_factor = INTEGER>`.

    .. option:: free_parameters = NUMERICAL_VECTOR

        See :opt:`free_parameters <free_parameters = NUMERICAL_VECTOR>`.

    .. option:: parameter_uncertainty

        See :opt:`parameter_uncertainty`.

    .. option:: regime = INTEGER

        See :opt:`regime <regime = INTEGER>`.

    .. option:: regimes

        See :opt:`regimes`.


Displaying and saving results
=============================

Dynare has comments to plot the results of a simulation and to save the results.

.. command:: rplot VARIABLE_NAME...;

    |br| Plots the simulated path of one or several variables, as
    stored in ``oo_.endo_simul`` by either
    ``perfect_foresight_solver``, ``simul`` (see :ref:`det-simul`) or
    ``stoch_simul`` with option ``periods`` (see
    :ref:`stoch-sol`). The variables are plotted in levels.

.. command:: dynatype (FILENAME) [VARIABLE_NAME...];

    |br| This command prints the listed variables in a text file named
    FILENAME. If no VARIABLE_NAME is listed, all endogenous variables
    are printed.

.. command:: dynasave (FILENAME) [VARIABLE_NAME...];

    |br| This command saves the listed variables in a binary file
    named FILENAME. If no VARIABLE_NAME are listed, all endogenous
    variables are saved.

    In MATLAB or Octave, variables saved with the ``dynasave command``
    can be retrieved by the command::

        load -mat FILENAME


.. _macro-proc-lang:

Macro-processing language
=========================

It is possible to use “macro” commands in the ``.mod`` file for doing
the following tasks: including modular source files, replicating
blocks of equations through loops, conditionally executing some code,
writing indexed sums or products inside equations...

The Dynare macro-language provides a new set of *macro-commands* which
can be inserted inside ``.mod`` files. It features:

    * File inclusion
    * Loops (for structure)
    * Conditional inclusion (``if/then/else`` structures)
    * Expression substitution

Technically, this macro language is totally independent of the basic
Dynare language, and is processed by a separate component of the
Dynare pre-processor. The macro processor transforms a ``.mod`` file
with macros into a ``.mod`` file without macros (doing
expansions/inclusions), and then feeds it to the Dynare parser. The
key point to understand is that the macro-processor only does text
substitution (like the C preprocessor or the PHP language). Note that
it is possible to see the output of the macro-processor by using the
``savemacro`` option of the ``dynare`` command (see :ref:`dyn-invoc`).

The macro-processor is invoked by placing *macro directives* in the
``.mod`` file. Directives begin with an at-sign followed by a pound
sign (``@#``). They produce no output, but give instructions to the
macro-processor. In most cases, directives occupy exactly one line of
text. In case of need, two backslashes (``\\``) at the end of the line
indicate that the directive is continued on the next line. The main
directives are:

    * ``@#includepath``, paths to search for files that are to be included,
    * ``@#include``, for file inclusion,
    * ``@#define``, for defining a macro-processor variable,
    * ``@#if, @#ifdef, @#ifndef, @#else, @#endif`` for conditional statements,
    * ``@#for, @#endfor`` for constructing loops.

The macro-processor maintains its own list of variables (distinct of
model variables and of MATLAB/Octave variables). These macro-variables
are assigned using the ``@#define`` directive, and can be of four
types: integer, character string, array of integers, array of strings.


.. _macro-exp:

Macro expressions
-----------------

It is possible to construct macro-expressions which can be assigned to
macro-variables or used within a macro-directive. The expressions are
constructed using literals of the four basic types (integers, strings,
arrays of strings, arrays of integers), macro-variables names and
standard operators.

String literals have to be enclosed between **double** quotes (like
``"name"``). Arrays are enclosed within brackets, and their elements
are separated by commas (like ``[1,2,3]`` or ``["US", "EA"]``).

Note that there is no boolean type: *false* is represented by integer
zero and *true* is any non-null integer. Further note that, as the
macro-processor cannot handle non-integer real numbers, integer
division results in the quotient with the fractional part truncated
(hence, :math:`5/3=3/3=1`).

The following operators can be used on integers:

    * Arithmetic operators: ``+, -, *, /``
    * Comparison operators: ``<, >, <=, >=, ==, !=``
    * Logical operators: ``&&, ||, !``
    * Integer ranges, using the following syntax:
      ``INTEGER1:INTEGER2`` (for example, ``1:4`` is equivalent to
      integer array ``[1,2,3,4]``)

The following operators can be used on strings:

    * Comparison operators: ``==, !=``
    * Concatenation of two strings: ``+``
    * Extraction of substrings: if ``s`` is a string, then ``s[3]`` is
      a string containing only the third character of ``s``, and
      ``s[4:6]`` contains the characters from 4th to 6th

The following operators can be used on arrays:

    * Dereferencing: if ``v`` is an array, then ``v[2]`` is its 2nd element.
    * Concatenation of two arrays: ``+``.
    * Difference ``-``: returns the first operand from which the
      elements of the second operand have been removed.
    * Extraction of sub-arrays: e.g. ``v[4:6]``.
    * Testing membership of an array: ``in`` operator (for example:
      ``"b"`` in ``["a", "b", "c"]`` returns ``1``)
    * Getting the length of an array: ``length`` operator (for
      example: ``length(["a", "b", "c"])`` returns ``3`` and, hence,
      ``1:length(["a", "b", "c"])`` is equivalent to integer array
      ``[1,2,3]``)

Macro-expressions can be used at two places:

    * Inside macro directives, directly;
    * In the body of the ``.mod`` file, between an at-sign and curly
      braces (like ``@{expr}``): the macro processor will substitute
      the expression with its value.

In the following, MACRO_EXPRESSION designates an expression
constructed as explained above.


Macro directives
----------------

.. macrodir:: @#includepath "PATH"
              @#includepath MACRO_VARIABLE

    |br| This directive adds the colon-separated paths contained in PATH to
    the list of those to search when looking for a ``.mod`` file
    specified by ``@#include``. Note that these paths are added
    *after* any paths passed using :opt:`-I <-I\<\<path\>\>>`.

    *Example*

        ::

            @#includepath "/path/to/folder/containing/modfiles:/path/to/another/folder"
            @#includepath folders_containing_mod_files


.. macrodir:: @#include "FILENAME"
              @#include MACRO_VARIABLE

    |br| This directive simply includes the content of another file at the
    place where it is inserted. It is exactly equivalent to a
    copy/paste of the content of the included file. Note that it is
    possible to nest includes (i.e. to include a file from an included
    file). The file will be searched for in the current directory. If
    it is not found, the file will be searched for in the folders
    provided by :opt:`-I <-I\<\<path\>\>>` and ``@#includepath``.

    *Example*

        ::

            @#include "modelcomponent.mod"
            @#include location_of_modfile


.. macrodir:: @#define MACRO_VARIABLE = MACRO_EXPRESSION

    |br| Defines a macro-variable.

    *Example*

        ::

            @#define x = 5              // Integer
            @#define y = "US"           // String
            @#define v = [ 1, 2, 4 ]    // Integer array
            @#define w = [ "US", "EA" ] // String array
            @#define z = 3 + v[2]       // Equals 5
            @#define t = ("US" in w)    // Equals 1 (true)

    *Example*

        ::

            @#define x = [ "B", "C" ]
            @#define i = 2

            model;
              A = @{x[i]};
            end;

        The latter is strictly equivalent to::

            model;
              A = C;
            end;


.. macrodir:: @#if MACRO_EXPRESSION
.. macrodir:: @#ifdef MACRO_VARIABLE
.. macrodir:: @#ifndef MACRO_VARIABLE
.. macrodir:: @#else
.. macrodir:: @#endif

    |br| Conditional inclusion of some part of the ``.mod`` file. The
    lines between ``@#if``, ``@#ifdef`` or ``@#ifndef`` and the next
    ``@#else`` or ``@#endif`` is executed only if the condition
    evaluates to a non-null integer. The ``@#else`` branch is optional
    and, if present, is only evaluated if the condition evaluates to
    ``0``.

    *Example*

        Choose between two alternative monetary policy rules using a
        macro-variable::

            @#define linear_mon_pol = 0 // or 1
            ...
            model;
            @#if linear_mon_pol
              i = w*i(-1) + (1-w)*i_ss + w2*(pie-piestar);
            @#else
              i = i(-1)^w * i_ss^(1-w) * (pie/piestar)^w2;
            @#endif
            ...
            end;

    *Example*

        Choose between two alternative monetary policy rules using a
        macro-variable. As ``linear_mon_pol`` was not previously defined
        in this example, the second equation will be chosen::

            model;
            @#ifdef linear_mon_pol
              i = w*i(-1) + (1-w)*i_ss + w2*(pie-piestar);
            @#else
              i = i(-1)^w * i_ss^(1-w) * (pie/piestar)^w2;
            @#endif
            ...
            end;


.. macrodir:: @#for MACRO_VARIABLE in MACRO_EXPRESSION
.. macrodir:: @#endfor

    |br| Loop construction for replicating portions of the ``.mod``
    file. Note that this construct can enclose variable/parameters
    declaration, computational tasks, but not a model declaration.

    *Example*

        ::

            model;
            @#for country in [ "home", "foreign" ]
              GDP_@{country} = A * K_@{country}^a * L_@{country}^(1-a);
            @#endfor
            end;

        The latter is equivalent to::

            model;
              GDP_home = A * K_home^a * L_home^(1-a);
              GDP_foreign = A * K_foreign^a * L_foreign^(1-a);
            end;


.. macrodir:: @#echo MACRO_EXPRESSION

    |br| Asks the preprocessor to display some message on standard
    output. The argument must evaluate to a string.


.. macrodir:: @#error MACRO_EXPRESSION

    |br| Asks the preprocessor to display some error message on standard
    output and to abort. The argument must evaluate to a string.

.. macrodir:: @#echomacrovars MACRO_EXPRESSION
              @#echomacrovars(save) MACRO_EXPRESSION

    |br| Asks the preprocessor to display the value of all macro
    variables up until this point. If the ``save`` option is passed,
    theh values of the macro variables are saved to
    ``options_.macrovars``.

Typical usages
--------------

Modularization
^^^^^^^^^^^^^^

The ``@#include`` directive can be used to split ``.mod`` files into
several modular components.

Example setup:

``modeldesc.mod``

    Contains variable declarations, model equations and shocks declarations.

``simul.mod``

    Includes ``modeldesc.mod``, calibrates parameters and runs
    stochastic simulations.

``estim.mod``

    Includes ``modeldesc.mod``, declares priors on parameters and runs
    Bayesian estimation.

Dynare can be called on ``simul.mod`` and ``estim.mod``, but it makes
no sense to run it on ``modeldesc.mod``.

The main advantage is that it is no longer needed to manually
copy/paste the whole model (at the beginning) or changes to the model
(during development).


Indexed sums of products
^^^^^^^^^^^^^^^^^^^^^^^^

The following example shows how to construct a moving average::

    @#define window = 2

    var x MA_x;
    ...
    model;
    ...
    MA_x = 1/@{2*window+1}*(
    @#for i in -window:window
            +x(@{i})
    @#endfor
           );
    ...
    end;

After macro-processing, this is equivalent to::

    var x MA_x;
    ...
    model;
    ...
    MA_x = 1/5*(
            +x(-2)
            +x(-1)
            +x(0)
            +x(1)
            +x(2)
           );
    ...
    end;


Multi-country models
^^^^^^^^^^^^^^^^^^^^

Here is a skeleton example for a multi-country model::

    @#define countries = [ "US", "EA", "AS", "JP", "RC" ]
    @#define nth_co = "US"

    @#for co in countries
    var Y_@{co} K_@{co} L_@{co} i_@{co} E_@{co} ...;
    parameters a_@{co} ...;
    varexo ...;
    @#endfor

    model;
    @#for co in countries
     Y_@{co} = K_@{co}^a_@{co} * L_@{co}^(1-a_@{co});
    ...
    @#if co != nth_co
     (1+i_@{co}) = (1+i_@{nth_co}) * E_@{co}(+1) / E_@{co}; // UIP relation
    @#else
     E_@{co} = 1;
    @#endif
    @#endfor
    end;


Endogeneizing parameters
^^^^^^^^^^^^^^^^^^^^^^^^

When doing the steady state calibration of the model, it may be useful
to consider a parameter as an endogenous (and vice-versa).

For example, suppose production is defined by a CES function:

    .. math::

           y = \left(\alpha^{1/\xi} \ell^{1-1/\xi}+(1-\alpha)^{1/\xi}k^{1-1/\xi}\right)^{\xi/(\xi-1)}

The labor share in GDP is defined as:

    .. math::

        \textrm{lab\_rat} = (w \ell)/(p y)

In the model, :math:`\alpha` is a (share) parameter, and ``lab_rat``
is an endogenous variable.

It is clear that calibrating :math:`\alpha` is not straightforward;
but on the contrary, we have real world data for ``lab_rat``, and it
is clear that these two variables are economically linked.

The solution is to use a method called *variable flipping*, which
consists in changing the way of computing the steady state. During
this computation, :math:`\alpha` will be made an endogenous variable
and ``lab_rat`` will be made a parameter. An economically relevant
value will be calibrated for ``lab_rat``, and the solution algorithm
will deduce the implied value for :math:`\alpha`.

An implementation could consist of the following files:

``modeqs.mod``

    This file contains variable declarations and model equations. The
    code for the declaration of :math:`\alpha` and ``lab_rat`` would
    look like::

        @#if steady
        var alpha;
         parameter lab_rat;
        @#else
         parameter alpha;
         var lab_rat;
        @#endif

``steady.mod``

    This file computes the steady state. It begins with::

        @#define steady = 1
        @#include "modeqs.mod"

    Then it initializes parameters (including ``lab_rat``, excluding
    :math:`\alpha`), computes the steady state (using guess values for
    endogenous, including :math:`\alpha`), then saves values of
    parameters and endogenous at steady state in a file, using the
    ``save_params_and_steady_state`` command.

``simul.mod``

    This file computes the simulation. It begins with::

        @#define steady = 0
        @#include "modeqs.mod"

    Then it loads values of parameters and endogenous at steady state
    from file, using the ``load_params_and_steady_state`` command, and
    computes the simulations.


MATLAB/Octave loops versus macro-processor loops
------------------------------------------------

Suppose you have a model with a parameter :math:`\rho`, and you want
to make simulations for three values: :math:`\rho = 0.8, 0.9,
1`. There are several ways of doing this:

*With a MATLAB/Octave loop*

    ::

        rhos = [ 0.8, 0.9, 1];
        for i = 1:length(rhos)
          rho = rhos(i);
          stoch_simul(order=1);
        end

    Here the loop is not unrolled, MATLAB/Octave manages the
    iterations. This is interesting when there are a lot of iterations.

*With a macro-processor loop (case 1)*

    ::

        rhos = [ 0.8, 0.9, 1];
        @#for i in 1:3
          rho = rhos(@{i});
          stoch_simul(order=1);
        @#endfor

    This is very similar to the previous example, except that the loop
    is unrolled. The macro-processor manages the loop index but not
    the data array (``rhos``).

*With a macro-processor loop (case 2)*

    ::

        @#for rho_val in [ "0.8", "0.9", "1"]
          rho = @{rho_val};
          stoch_simul(order=1);
        @#endfor

    The advantage of this method is that it uses a shorter syntax,
    since list of values directly given in the loop construct. Note
    that values are given as character strings (the macro-processor
    does not know floating point values). The inconvenience is that
    you can not reuse an array stored in a MATLAB/Octave variable.


Verbatim inclusion
==================

Pass everything contained within the verbatim block to the
``<mod_file>.m`` file.

.. block:: verbatim ;

    |br| By default, whenever Dynare encounters code that is not
    understood by the parser, it is directly passed to the
    preprocessor output.

    In order to force this behavior you can use the ``verbatim``
    block. This is useful when the code you want passed to the
    ``<mod_file>.m`` file contains tokens recognized by the Dynare
    preprocessor.

    *Example*

        ::

            verbatim;
            % Anything contained in this block will be passed
            % directly to the <modfile>.m file, including comments
            var = 1;
            end;


Misc commands
=============

.. command:: set_dynare_seed (INTEGER)
             set_dynare_seed (`default')
             set_dynare_seed (`clock')
             set_dynare_seed (`reset')
             set_dynare_seed (`ALGORITHM', INTEGER)

    |br| Sets the seed used for random number generation. It is
    possible to set a given integer value, to use a default value, or
    to use the clock (by using the latter, one will therefore get
    different results across different Dynare runs). The ``reset``
    option serves to reset the seed to the value set by the last
    ``set_dynare_seed`` command. On MATLAB 7.8 or above, it is also
    possible to choose a specific algorithm for random number
    generation; accepted values are ``mcg16807``, ``mlfg6331_64``,
    ``mrg32k3a``, ``mt19937ar`` (the default), ``shr3cong`` and
    ``swb2712``.

.. command:: save_params_and_steady_state (FILENAME);

    |br| For all parameters, endogenous and exogenous variables,
    stores their value in a text file, using a simple name/value
    associative table.

        * for parameters, the value is taken from the last parameter
          initialization.
        * for exogenous, the value is taken from the last ``initval``
          block.
        * for endogenous, the value is taken from the last steady
          state computation (or, if no steady state has been computed,
          from the last ``initval`` block).

    Note that no variable type is stored in the file, so that the
    values can be reloaded with ``load_params_and_steady_state`` in a
    setup where the variable types are different.

    The typical usage of this function is to compute the steady-state
    of a model by calibrating the steady-state value of some
    endogenous variables (which implies that some parameters must be
    endogeneized during the steady-state computation).

    You would then write a first ``.mod`` file which computes the
    steady state and saves the result of the computation at the end of
    the file, using ``save_params_and_steady_state``.

    In a second file designed to perform the actual simulations, you
    would use ``load_params_and_steady_state`` just after your
    variable declarations, in order to load the steady state
    previously computed (including the parameters which had been
    endogeneized during the steady state computation).

    The need for two separate ``.mod`` files arises from the fact that
    the variable declarations differ between the files for steady
    state calibration and for simulation (the set of endogenous and
    parameters differ between the two); this leads to different
    ``var`` and ``parameters`` statements.

    Also note that you can take advantage of the ``@#include``
    directive to share the model equations between the two files (see
    :ref:`macro-proc-lang`).

.. command:: load_params_and_steady_state (FILENAME);

    |br| For all parameters, endogenous and exogenous variables, loads
    their value from a file created with
    ``save_params_and_steady_state``.

        * for parameters, their value will be initialized as if they
          had been calibrated in the ``.mod`` file.
        * for endogenous and exogenous variables, their value will be
          initialized as they would have been from an ``initval``
          block .

    This function is used in conjunction with
    ``save_params_and_steady_state``; see the documentation of that
    function for more information.

.. matcomm:: dynare_version ;

    |br| Output the version of Dynare that is currently being used
    (i.e. the one that is highest on the MATLAB/Octave path).

.. matcomm:: write_latex_definitions ;

    |br| Writes the names, LaTeX names and long names of
    model variables to tables in a file named
    ``<<M_.fname>>_latex_definitions.tex``. Requires the following
    LaTeX packages: ``longtable``.

.. matcomm:: write_latex_parameter_table ;

    |br| Writes the LaTeX names, parameter names, and long
    names of model parameters to a table in a file named
    ``<<M_.fname>>_latex_parameters.tex.`` The command writes the
    values of the parameters currently stored. Thus, if parameters are
    set or changed in the steady state computation, the command should
    be called after a steady-command to make sure the parameters were
    correctly updated. The long names can be used to add parameter
    descriptions. Requires the following LaTeX packages:
    ``longtable, booktabs``.

.. matcomm:: write_latex_prior_table ;

    |br| Writes descriptive statistics about the prior distribution to
    a LaTeX table in a file named
    ``<<M_.fname>>_latex_priors_table.tex``. The command writes the
    prior definitions currently stored. Thus, this command must be
    invoked after the ``estimated_params`` block. If priors are
    defined over the measurement errors, the command must also be
    preceeded by the declaration of the observed variables (with
    ``varobs``). The command displays a warning if no prior densities
    are defined (ML estimation) or if the declaration of the observed
    variables is missing. Requires the following LaTeX
    packages: ``longtable, booktabs``.

.. matcomm:: collect_latex_files ;

    |br| Writes a LaTeX file named
    ``<<M_.fname>>_TeX_binder.tex`` that collects all TeX output
    generated by Dynare into a file. This file can be compiled using
    ``pdflatex`` and automatically tries to load all required
    packages. Requires the following LaTeX packages:
    ``breqn``, ``psfrag``, ``graphicx``, ``epstopdf``, ``longtable``,
    ``booktabs``, ``caption``, ``float,`` ``amsmath``, ``amsfonts``,
    and ``morefloats``.


.. _Dynare wiki: http://www.dynare.org/DynareWiki/EquationsTags
.. _io: http://octave.sourceforge.net/io/
.. _AIM website: http://www.federalreserve.gov/Pubs/oss/oss4/aimindex.html

.. rubric:: Footnotes

.. [#f1] Note that arbitrary MATLAB or Octave expressions can be put
         in a ``.mod`` file, but those expressions have to be on
         separate lines, generally at the end of the file for
         post-processing purposes. They are not interpreted by Dynare,
         and are simply passed on unmodified to MATLAB or
         Octave. Those constructions are not addresses in this
         section.

.. [#f2] In particular, for big models, the compilation step can be
         very time-consuming, and use of this option may be
         counter-productive in those cases.

.. [#f3] See option :ref:`conf_sig <confsig>` to change the size of
         the HPD interval.

.. [#f4] See option :ref:`conf_sig <confsig>` to change the size of
         the HPD interval.

.. [#f5] When the shocks are correlated, it is the decomposition of
         orthogonalized shocks via Cholesky decomposition according to
         the order of declaration of shocks (see :ref:`var-decl`)

.. [#f6] See :opt:`forecast <forecast = INTEGER>` for more information.

.. [#f7] In case of Excel not being installed,
         `<https://mathworks.com/matlabcentral/fileexchange/38591-xlwrite--generate-xls-x--files-without-excel-on-mac-linux-win>`_
         may be helpful.

.. [#f8] See option :ref:`conf_sig <confsig>` to change the size of
         the HPD interval.

.. [#f9] See option :ref:`conf_sig <confsig>` to change the size of
         the HPD interval.

.. [#f10] If you want to align the paper with the description herein,
          please note that :math:`A` is :math:`A^0` and :math:`F` is
          :math:`A^+`.

.. [#f11] An example can be found at
          `<https://github.com/DynareTeam/dynare/blob/master/tests/ms-dsge/test_ms_dsge.mod>`_.
