![](https://travis-ci.org/DynareTeam/reporting.svg?branch=master)

# Dynare Reporting

Originally designed to provide reporting functionality for
[Dynare](http://www.dynare.org), it has been moved to its own
repository in the hopes that it can be used without obliging the user
to download Dynare in its entirety.

# License

Dynare Reporting is covered by the GNU General Public Licence version 3 or later
(see [license.txt](license.txt) in the Dynare Reporting distribution for
specifics).

# Obtain the code that Dynare Reporting depends on

Dynare ```reporting``` depends on the Dynare
[```dates```](https://github.com/DynareTeam/dates) and
[```dseries```](https://github.com/DynareTeam/dseries) repositories.

# Use the Dynare Reporting code

- Open Matlab/Octave
- Change into the ```reporting``` directory
- Ensure that the [```dates```](https://github.com/DynareTeam/dates)
  and [```dseries```](https://github.com/DynareTeam/dseries)
  directories are in your path and initialized correctly (follow the
  directions on the respective sites)
- Run ```initialize_reporting_toolbox```
- Use the reporting code as outlined in the Dynare documentation

# Run the example code

- Follow the steps above
- Change into the ```reporting/examples``` directory
- Run the ```run_report_example``` matlab script
