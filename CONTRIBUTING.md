# Instructions for Contributing to Dynare

## Introduction

Hello from the Dynare Team! We're happy you're on this page because hopefully that means you're thinking of getting directly involved with the Dynare project. Herein, we outline how you can contribute to Dynare. Please read this document all the way through before contributing.

Please follow the steps in the sections below in order. Note that, though we'd love for you to contribute code, you don't need to be a programmer to contribute to Dynare. You can report bugs, ask for enhancements, fix typos in the manual, contribute tests to the test suite, or do something we haven't thought of yet!

If something is not clear, don't hesitate to ask if you can't find the answer online. You can contact us directly at [dev@dynare.org](mailto:dev@dynare.org).

Please note that the repositories under the purview of this page are:

* [Dynare](https://git.dynare.org/Dynare/dynare)
* [Particles](https://git.dynare.org/Dynare/particles)
* [Dates](https://git.dynare.org/Dynare/dates)
* [Dseries](https://git.dynare.org/Dynare/dseries)
* [Reporting](https://git.dynare.org/Dynare/reporting)
* [Testsuite](https://git.dynare.org/Dynare/testsuite)
* [M-unit-tests](https://git.dynare.org/Dynare/m-unit-tests)

## Making your Intentions Known

Before making changes to the codebase, it'd be helpful if you communicated your intentions with us. This will avoid double work and ensure that you don't contribute code that won't be included in Dynare for one reason or another. The way to communicate with us is via [GitHub Issues](https://guides.github.com/features/issues/).

### Issues: Reporting Bugs

You can report bugs in both the stable and unstable versions of Dynare. Before reporting a bug in the stable version of Dynare, please check the [Known Bugs](http://www.dynare.org/DynareWiki/KnownBugs) page to ensure it has not already been encountered/fixed. If reporting a bug in the unstable version of Dynare, please ensure the bug exists in the latest [unstable Dynare snapshot](http://www.dynare.org/download/dynare-unstable).

To report a bug in Dynare, simply open a GitHub issue in the repository where the bug resides. For example, to report a bug in Dynare itself, go to the [Dynare repository issue page](https://git.dynare.org/Dynare/dynare/issues) and click on "New Issue."

The minimal information to add is a subject and a description of the steps needed to reproduce the bug. However, the most helpful description would also provide the code to reproduce the bug (often times a `.mod` file). The most helpful `.mod` file is a minimal, quick-running example that reproduces the bug, but we'll take anything that will help us squash a bug.

To include short amounts of code, please paste it into the description box, using the appropriate [GitHub markdown](https://help.github.com/articles/github-flavored-markdown/) code. For larger amounds of code like `.mod` files, please create a new [GitHub Gist](https://gist.github.com) and provide the link in the description box.

### Issues: Enhancements

Issues are not only used to report bugs. They are also used to ask for improvements to the codebase or new features to Dynare in general. Please be descriptive when asking for improvements or new features. Links or references to papers or detailed examples are helpful.

Though our development priorities lay with those who finance Dynare and with what we think may most benefit the Dynare community, this does not mean we are closed to outside ideas for enhancements. On the contrary: we invite them! Moreover, if you are willing to program the enhancement you want, the odds of it being included in Dynare are much higher than if you needed us to do it. That said, it is best to create an issue with an enhancement idea **before** beginning the work. As stated above, this is important to avoid duplication of work and also because we wouldn't want you to take the time to work on something that would not end up being included in Dynare.

## Get to Coding!

So, now you've reported the bug or asked for an enhancemnt by creating a GitHub issue. That's already a great help. Thank you!

Now, if you want to go the extra mile, you'll volunteer to contribute code to fix the GitHub issue you created above. Once we've agreed that you'll do it, please do the following:

1. Clone the Dynare repository:
   * `git clone https://git.dynare.org/Dynare/dynare.git`
1. [Fork the Dynare repository](https://help.github.com/articles/fork-a-repo)
1. Change into the `dynare` folder and add the forked repository as a remote:
   * `cd dynare`
   * `git remote add personal https://github.com/<<GitHub username>>/dynare.git`
1. Create a branch to work on
   * `git checkout -b <<descriptive branch name>>`
1. Do your work, all the while respecting the [Dynare Coding Standards](http://www.dynare.org/DynareWiki/CodingStandards)

As you work, your forked repository will likely fall out of sync with the main Dynare repository as we'll be working in parallel. No matter. Follow these steps to ensure your changes will be merge-able when they're done:

1. Get the changes from the main Dynare repository:
   * `git checkout master`
   * `git fetch`
   * `git rebase origin master`
1. Move your changes on top of the current master branch of Dynare
   * `git checkout <<descriptive branch name>>`
   * `git rebase origin/master`
   * This last command may cause a conflict. It is up to you to resolve this conflict.

Once you've made the changes necessary to fix the bug or add an enhancement, ensure that it has been rebased on the master branch (following the steps above), commit it, push it to your forked Dynare repository, and create a pull request:

1. Get the latest changes from Dynare and rebase your branch on top of them (see above)
1. Commit your changes:
   * `git add <<files to commit>>`
   * `git commit -m "<<descriptive commit message.>> Closes #<<Ref. to GitHub issue number fixed by this commit>>"`
1. Push to your forked Dynare repository
   * `git push personal <<descriptive branch name>>`
1. Create a [Pull Request](https://help.github.com/articles/creating-a-pull-request/) from the branch in your forked Dynare repository

## Tests

The Dynare Test Suite runs nightly. It's how we quickly catch bugs that may have been introduced by changes made during the day. The output from the test suite can be found here: [http://www.dynare.org/testsuite/master/](http://www.dynare.org/testsuite/master/). This is also a good place to start fixing bugs. If you see a `.mod` file that doesn't run in the test suite and think you can fix it, create an issue and once you have the go ahead, go for it!

### Test `.mod` File

It's useful to contribute `.mod` files that test some aspect of Dynare that is not currently tested. A `.mod` file that runs into a bug is perfect. As the test suite currently takes several hours to run, we prefer you modify a current test to also create the bug you've found. If you can't do that, please add a new test that runs as quickly as possible. It will contain only those commands necessary to create the bug, nothing more. To contribute a test, after having made an issue and cloned and forked the repository as described above, do the following:

1. Modify the `MODFILES` variable in `tests/Makefile.am` with a line containing your test file name
1. If any ancillary files are needed to run your test, please include them in the `EXTRA_DIST` variable in `tests/Makefile.am`
1. Add and commit your test file and `tests/Makefile.am` as described above
1. Push and create a pull request as described above

NB: Please do not contribute non-text files (e.g. `.xls`). If you absolutely need such a file for your test to run, please contact us to see how to go about contributing it.

### Unitary Tests
