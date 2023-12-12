.. default-domain:: dynare

.. |br| raw:: html

    <br>

.. _conf-file:

######################
The configuration file
######################

The configuration file is used to provide Dynare with information not
related to the model (and hence not placed in the model file). At the
moment, it is only used when using Dynare to run parallel
computations.

On Linux and macOS, the configuration file is searched by default under
``dynare/dynare.ini`` in the configuration directories defined by the XDG
specification (typically ``$HOME/.config/dynare/dynare.ini`` for the
user-specific configuration and ``/etc/xdg/dynare/dynare.ini`` for the
system-wide configuration, the former having precedence over the latter). Under
Windows, the configuration file is searched by default in
``%APPDATA%\dynare\dynare.ini`` (typically
``c:\Users\USERNAME\AppData\Roaming\dynare\dynare.ini``). You can specify a non
standard location using the ``conffile`` option of the ``dynare`` command (see
:ref:`dyn-invoc`).

The parsing of the configuration file is case-sensitive and it should
take the following form, with each option/choice pair placed on a
newline::

    [command0]
    option0 = choice0
    option1 = choice1

    [command1]
    option0 = choice0
    option1 = choice1

The configuration file follows a few conventions (self-explanatory
conventions such as ``USER_NAME`` have been excluded for concision):

``COMPUTER_NAME``

    Indicates the valid name of a server (e.g. ``localhost``,
    ``server.cepremap.org``) or an IP address.

``DRIVE_NAME``

    Indicates a valid drive name in Windows, without the trailing
    colon (e.g. ``C``).

``PATH``

    Indicates a valid path in the underlying operating system
    (e.g. ``/home/user/dynare/matlab/``).

``PATH_AND_FILE``

    Indicates a valid path to a file in the underlying operating
    system (e.g. ``/usr/local/MATLAB/R2023b/bin/matlab``).

``BOOLEAN``

    Is ``true`` or ``false``.


Dynare Configuration
====================

This section explains how to configure Dynare for general
processing. Currently, there is only one option available.

.. confblock:: [hooks]

    |br| This block can be used to specify configuration options that will
    be used when running Dynare.

    *Options*

    .. option:: GlobalInitFile = PATH_AND_FILE

        The location of a global initialization file that can be used to
        customize some Dynare internals (typically default option values). This
        is a MATLAB/Octave script.

        If this option is not specified, Dynare will look for a
        ``global_init.m`` file in its configuration directory (typically
        ``$HOME/.config/dynare/global_init.m`` under Linux and macOS, and
        ``c:\Users\USERNAME\AppData\Roaming\dynare\global_init.m`` under
        Windows).

    *Example*

        ::

            [hooks]
            GlobalInitFile = /home/usern/dynare/myInitFile.m


.. confblock:: [paths]

    |br| This block can be used to specify paths that will be used
    when running dynare.

    *Options*

    .. option:: Include = PATH

        A colon-separated path to use when searching for files to
        include via ``@#include``. Paths specified via :opt:`-I
        <-I\<\<path\>\>>` take priority over paths specified here,
        while these paths take priority over those specified by
        ``@#includepath``.

    *Example*

        ::

            [paths]
            Include = /path/to/folder/containing/modfiles:/path/to/another/folder

.. _paral-conf:

Parallel Configuration
======================

This section explains how to configure Dynare for parallelizing some
tasks which require very little inter-process communication.

The parallelization is done by running several MATLAB or Octave
processes, either on local or on remote machines. Communication
between leader and follower processes are done through SMB on Windows and
SSH on UNIX. Input and output data, and also some short status
messages, are exchanged through network filesystems. Currently the
system works only with homogenous grids: only Windows or only Unix
machines.

The following routines are currently parallelized:

    * the posterior sampling algorithms when using multiple chains;
    * the Metropolis-Hastings diagnostics;
    * the posterior IRFs;
    * the prior and posterior statistics;
    * some plotting routines.

Note that creating the configuration file is not enough in order to
trigger parallelization of the computations: you also need to specify
the ``parallel`` option to the ``dynare`` command. For more details,
and for other options related to the parallelization engine, see
:ref:`dyn-invoc`.

You also need to verify that the following requirements are met by
your cluster (which is composed of a leader and of one or more
followers):

For a Windows grid:

        * a standard Windows network (SMB) must be in place;
        * the `PsTools`_ suite must be installed in the path of the
          leader Windows machine;
        * the Windows user on the leader machine has to be user of any
          other follower machine in the cluster, and that user will be
          used for the remote computations.
        * detailed step-by-step setup instructions can be found in
          :ref:`win-ssg`.

For a UNIX grid:

        * SSH must be installed on the leader and on the follower machines;
        * SSH keys must be installed so that the SSH connection from
          the leader to the follower can be done without passwords, or
          using an SSH agent.

.. warning:: Compatibility considerations between leader and follower

    It is highly recommended to use the same version of Dynare on both the
    leader and all followers. Different versions regularly cause problems like
    zero acceptance rates during estimation. When upgrading to a newer Dynare
    version do not forget to adjust the ``DynarePath``.

We now turn to the description of the configuration directives. Note
that comments in the configuration file can be provided by separate
lines starting with a hashtag (#).

.. confblock:: [cluster]

    |br| When working in parallel, ``[cluster]`` is required to specify the
    group of computers that will be used. It is required even if you
    are only invoking multiple processes on one computer.

    *Options*

    .. option:: Name = CLUSTER_NAME

        The reference name of this cluster.

    .. option:: Members = NODE_NAME[(WEIGHT)] NODE_NAME[(WEIGHT)] ...

        A list of nodes that comprise the cluster with an optional
        computing weight specified for that node. The computing weight
        indicates how much more powerful one node is with respect to
        the others (e.g. ``n1(2) n2(1) n3(3)`` means that ``n1`` is
        two times more powerful than ``n2`` whereas ``n3`` is three
        times more powerful than ``n2``). Each node is separated by at
        least one space and the weights are in parenthesis with no
        spaces separating them from their node.

    *Example*

        ::

            [cluster]
            Name = c1
            Members = n1 n2 n3

            [cluster]
            Name = c2
            Members = n1(4) n2 n3


.. confblock:: [node]

    |br| When working in parallel, ``[node]`` is required for every
    computer that will be used. The options that are required differ,
    depending on the underlying operating system and whether you are
    working locally or remotely.

    *Options*

    .. option:: Name = NODE_NAME

        The reference name of this node.

    .. option:: CPUnbr = INTEGER | [INTEGER:INTEGER]

        If just one integer is passed, the number of processors to
        use. If a range of integers is passed, the specific processors
        to use (processor counting is defined to begin at one as
        opposed to zero). Note that using specific processors is only
        possible under Windows; under Linux and macOS, if a range is
        passed the same number of processors will be used but the
        range will be adjusted to begin at one.

    .. option:: ComputerName = COMPUTER_NAME

        The name or IP address of the node. If you want to run
        locally, use ``localhost`` (case-sensitive).

    .. option:: Port = INTEGER

        The port number to connect to on the node. The default is
        empty, meaning that the connection will be made to the default
        SSH port (22).

    .. option:: UserName = USER_NAME

        The username used to log into a remote system. Required for
        remote runs on all platforms.

    .. option:: Password = PASSWORD

        The password used to log into the remote system. Required for
        remote runs originating from Windows.

    .. option:: RemoteDrive = DRIVE_NAME

        The drive to be used for remote computation. Required for
        remote runs originating from Windows.

    .. option:: RemoteDirectory = PATH

        The directory to be used for remote computation. Required for
        remote runs on all platforms.

    .. option:: DynarePath = PATH

        The path to the matlab subdirectory within the Dynare
        installation directory. The default is the empty string.

    .. option:: MatlabOctavePath = PATH_AND_FILE

        The path to the MATLAB or Octave executable. The default value
        is ``matlab`` as MATLAB’s executable is typically in the %PATH% environment 
        variable. When using full paths on Windows, you may need to enclose the path
        in quoted strings, e.g. ``MatlabOctavePath="C:\Program Files\MATLAB\R2023b\bin\matlab.exe"``

    .. option:: NumberOfThreadsPerJob = INTEGER

        This option controls the distribution of jobs (e.g. MCMC chains) across additional MATLAB instances that are run in parallel.
        Needs to be an exact divisor of the number of cores.
        The formula :opt:`CPUnbr <CPUnbr = INTEGER | [INTEGER:INTEGER]>` divided by :opt:`NumberOfThreadsPerJob <NumberOfThreadsPerJob = INTEGER>`
        calculates the number of MATLAB/Octave instances that will be launched in parallel,
        where each instance will then execute a certain number of jobs sequentially.
        For example, if you run a MCMC estimation with 24 chains on a 12 core machine, setting ``CPUnbr = 12`` and ``NumberOfThreadsPerJob = 4``
        will launch 3 MATLAB instances in parallel, each of which will compute 8 chains sequentially.
        Note that this option does not dictate the number of maximum threads utilized by each MATLAB/Octave instance,
        see related option :opt:`SingleCompThread <SingleCompThread = BOOLEAN>` for this.
        Particularly for very large models, setting this option to 2 might distribute the workload in a 
        more efficient manner, depending on your hardware and task specifics.
        It’s advisable to experiment with different values to achieve optimal performance.
        The default value is ``1``.        
        

    .. option:: SingleCompThread = BOOLEAN

        This option allows you to enable or disable MATLAB’s native multithreading capability. When set to ``true``, 
        the additional MATLAB instances are initiated in single thread mode utilizing the ``-singleCompThread`` startup option, 
        thereby disabling MATLAB’s native multithreading. When set to ``false``, MATLAB’s native multithreading 
        is enabled, e.g. the actual number of threads utilized by each MATLAB instance is usually determined by the number of CPU cores
        (you can check this by running ``maxNumCompThreads`` in MATLAB’s command window).
        Note: While MATLAB aims to accelerate calculations by distributing them across your computer’s threads, 
        certain tasks, like MCMC estimations, may exhibit slowdowns with MATLAB’s multitasking especially when Dynare’s parallel computing is turned on
        as we do not use MATLAB’s parallel toolbox.
        So in many cases it is advisable to set this setting to ``true``.
        If you want to have more control, you can manually add the MATLAB command `maxNumCompThreads(N)` at the beginning of `fParallel.m`.
        The default value is ``false``. This option is ineffective under Octave.

        
    .. option:: OperatingSystem = OPERATING_SYSTEM

        The operating system associated with a node. Only necessary
        when creating a cluster with nodes from different operating
        systems. Possible values are ``unix`` or ``windows``. There is
        no default value.

    *Example*

        ::

            [node]
            Name = n1
            ComputerName = localhost
            CPUnbr = 1

            [node]
            Name = n2
            ComputerName = dynserv.cepremap.org
            CPUnbr = 5
            UserName = usern
            RemoteDirectory = /home/usern/Remote
            DynarePath = /home/usern/dynare/matlab
            MatlabOctavePath = matlab

            [node]
            Name = n3
            ComputerName = dynserv.dynare.org
            Port = 3333
            CPUnbr = [2:4]
            UserName = usern
            RemoteDirectory = /home/usern/Remote
            DynarePath = /home/usern/dynare/matlab
            MatlabOctavePath = matlab

.. _win-ssg:

Windows Step-by-Step Guide
==========================

This section outlines the steps necessary on most Windows systems to
set up Dynare for parallel execution. Note that the steps 3 to 6 are 
required unless parallel execution is confined to a local pool 
with the ``parallel_use_psexec=false`` option. 

    1. Write a configuration file containing the options you want. A
       mimimum working example setting up a cluster consisting of two
       local CPU cores that allows for e.g. running two Monte Carlo
       Markov Chains in parallel is shown below.
    2. Save the configuration file somwhere. The name and file ending
       do not matter if you are providing it with the ``conffile``
       command line option. The only restrictions are that the path
       must be a valid filename, not contain non-alpha-numeric
       characters, and not contain any whitespaces. For the
       configuration file to be accessible without providing an
       explicit path at the command line, you must save it under the
       name ``dynare.ini`` into your user account’s ``Application
       Data`` folder.
    3. Install `PSTools`_ to your system, e.g. into ``C:\PSTools.``
    4. Set the Windows System Path to the ``PSTools`` folder
       (e.g. using something along the line of pressing Windows
       Key+Pause to open the System Configuration, then go to Advanced
       -> Environment Variables -> Path).
    5. Restart your computer to make the path change effective.
    6. Open MATLAB and type into the command window::

           !psexec

       This executes the ``psexec.exe`` from PSTools on your system
       and shows whether Dynare will be able to locate it. If MATLAB
       complains at this stage, you did not correctly set your Windows
       system path for the ``PSTools`` folder.
    7. If ``psexec.exe`` was located in the previous step, a popup
       will show up, asking for confirmation of the license
       agreement. Confirm this copyright notice of ``psexec`` (this
       needs to be done only once). After this, Dynare should be ready
       for parallel execution.
    8. Call Dynare on your mod-file invoking the ``parallel`` option
       and providing the path to your configuration file with the
       ``conffile`` option (if you did not save it as
       ``%APPDATA%\dynare.ini`` in step 2 where it should be detected
       automatically)::

            dynare ls2003 parallel conffile='C:\Users\Dynare~1\parallel\conf_file.ini'

    Please keep in mind that no white spaces or names longer than 8
    characters are allowed in the ``conffile`` path. The 8-character
    restriction can be circumvented by using the tilde Windows path
    notation as in the above example.

*Example*::

    #cluster needs to always be defined first
    [cluster]
    #Provide a name for the cluster
    Name=Local
    #declare the nodes being member of the cluster
    Members=n1

    #declare nodes (they need not all be part of a cluster)
    [node]
    #name of the node
    Name=n1
    #name of the computer (localhost for the current machine)
    ComputerName=localhost
    #cores to be included from this node
    CPUnbr=[1:2]
    #path to matlab.exe; on Windows, the MATLAB bin folder is in the system path
    #so we only need to provide the name of the exe file
    MatlabOctavePath=matlab
    #Dynare path you are using
    DynarePath=C:/dynare/4.7.0/matlab

.. _PsTools: https://technet.microsoft.com/sysinternals/pstools.aspx
