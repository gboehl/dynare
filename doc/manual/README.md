# Dynare reference manual

Source for the Dynare's reference manual which is available at
<http://www.dynare.org/manual-unstable>. The documentation is obtained
using [Sphinx](http://www.sphinx-doc.org/) Pyhthon Documentation
Generator. To build the `html` version of the reference manual just type:

```bash
~$ make html
```

on the command line. The reference manual will be available under
`src/build/html`, provided all the following (debian package)
requirements are met:

 - python3,
 - python3-sphinx,
 - python3-recommonmark, and
 - libjs-mathjax
