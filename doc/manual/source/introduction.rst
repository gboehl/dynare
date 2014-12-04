.. _introduction:

**************
 Introduction
**************

.. introduction_what-is_dynare:

What is Dynare
--------------


Dynare is a software platform for handling a wide class of economic
models, in particular dynamic stochastic general equilibrium (DSGE)
and overlapping generations (OLG) models. The models solved by Dynare
include those relying on the rational expectations hypothesis, wherein
agents form their expectations about the future in a way consistent
with the model. But Dynare is also able to handle models where
expectations are formed differently: on one extreme, models where
agents perfectly anticipate the future; on the other extreme, models
where agents have limited rationality or imperfect knowledge of the
state of the economy and, hence, form their expectations through a
learning process. In terms of types of agents, models solved by Dynare
can incorporate consumers, productive firms, governments, monetary
authorities, investors and financial intermediaries. Some degree of
heterogeneity can be achieved by including several distinct classes of
agents in each of the aforementioned agent categories.

Dynare offers a user-friendly and intuitive way of describing these
models. It is able to perform simulations of the model given a
calibration of the model parameters and is also able to estimate these
parameters given a dataset. In practice, the user will write a text
file containing the list of model variables, the dynamic equations
linking these variables together, the computing tasks to be performed
and the desired graphical or numerical outputs.

A large panel of applied mathematics and computer science techniques
are internally employed by Dynare: multivariate nonlinear solving and
optimization, matrix factorizations, local functional approximation,
Kalman filters and smoothers, MCMC techniques for Bayesian estimation,
graph algorithms, optimal control, …

Various public bodies (central banks, ministries of economy and
finance, international organisations) and some private financial
institutions use Dynare for performing policy analysis exercises and
as a support tool for forecasting exercises. In the academic world,
Dynare is used for research and teaching purposes in postgraduate
macroeconomics courses.

Dynare is a free software, which means that it can be downloaded free
of charge, that its source code is freely available, and that it can
be used for both non-profit and for-profit purposes. Most of the
source files are covered by the GNU General Public Licence (GPL)
version 3 or later (there are some exceptions to this, see the file
license.txt in Dynare distribution). It is available for the Windows,
Mac and Linux platforms and is fully documented through a user guide
and a reference manual. Part of Dynare is programmed in C++, while the
rest is written using the `MATLAB`_ programming language. The latter
implies that commercially-available `MATLAB`_ software is required in
order to run Dynare. However, as an alternative to `MATLAB`_, Dynare
is also able to run on top of `GNU Octave`_ (basically a free clone of
`MATLAB`_): this possibility is particularly interesting for students
or institutions who cannot afford, or do not want to pay for,
`MATLAB`_ and are willing to bear the concomitant performance loss.

The development of Dynare is mainly done at `Cepremap`_ by a core team
of researchers who devote part of their time to software
development. Currently the development team of Dynare is composed of
Stéphane Adjemian (`Université du Maine`_, Gains), Houtan Bastani
(`Cepremap`_), Michel Juillard (`Banque de France`_), Frédéric Karamé
(`Université du Maine`_, Gains and `Cepremap`_), Junior Maih (Norges
Bank), Ferhat Mihoubi (Université Paris-Est Créteil, Epee and
Cepremap), George Perendia, Johannes Pfeifer (University of Mannheim),
Marco Ratto (JRC) and Sébastien Villemot (`OFCE – Sciences Po`_). Increasingly,
the developer base is expanding, as tools developed by researchers
outside of Cepremap are integrated into Dynare. Financial support is
provided by `Cepremap`_, `Banque de France`_ and `DSGE-net`_ (an
international research network for DSGE modeling). The Dynare project
also received funding through the Seventh Framework Programme for
Research (FP7) of the European Commission’s Socio-economic Sciences
and Humanities (SSH) Program from October 2008 to September 2011 under
grant agreement SSH-CT-2009-225149.

Interaction between developers and users of Dynare is central to the
project. A web forum is available for users who have questions about
the usage of Dynare or who want to report bugs. Training sessions are
given through the Dynare Summer School, which is organized every year
and is attended by about 40 people. Finally, priorities in terms of
future developments and features to be added are decided in
cooperation with the institutions providing financial support.

.. introduction_documentation-sources:

Documentation sources
---------------------

The present document is the reference manual for Dynare. It documents
all commands and features in a systematic fashion. New users may
rather begin with Dynare User Guide (Mancini (2007)), distributed with
Dynare and also available from the official Dynare web site.

Other useful sources of information include the `Dynare wiki`_ and the `Dynare forums`_.

.. introduction_citation:

Citing Dynare in your research
------------------------------

If you would like to refer to Dynare in a research article, the
recommended way is to cite the present manual, as follows:


    Stéphane Adjemian, Houtan Bastani, Michel Juillard, Frédéric
    Karamé, Ferhat Mihoubi, George Perendia, Johannes Pfeifer, Marco
    Ratto and Sébastien Villemot (2011), “*Dynare: Reference Manual,
    Version 4*,” Dynare Working Papers, 1, CEPREMAP

Note that citing the Dynare Reference Manual in your research is a
good way to help the Dynare project.

If you want to give a URL, use the address of the Dynare website:
http://www.dynare.org.




.. _Dynare wiki: http://www.dynare.org/DynareWiki
.. _Dynare forums: http://www.dynare.org/phpBB3
.. _GNU Octave: http://www.octave.org/
.. _MATLAB: http://www.mathworks.com/products/matlab/
.. _Cepremap: http://www.cepremap.fr
.. _Université du Maine: http://www.univ-lemans.fr
.. _Banque de France: http://www.banque-france.fr
.. _DSGE-net: http://www.dsge.net
.. _OFCE – Sciences Po: http://www.ofce.fr
