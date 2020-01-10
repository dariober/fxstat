[![Build Status](https://travis-ci.com/dariober/fxstat.svg?branch=master)](https://travis-ci.com/dariober/fxstat)


<!-- vim-markdown-toc GFM -->

* [Compile and setup](#compile-and-setup)
* [Usage](#usage)
* [Unit tests and code coverage](#unit-tests-and-code-coverage)

<!-- vim-markdown-toc -->

Fast and easy stats for FASTA and FASTQ files

Compile and setup
=================

Assuming you are in the root directory of this project (i.e. where this README
file is), compile the fxstat executable with:

```
gcc -O3 -Wall -g src/fxstat.c src/utils.c -o fxstat
```

Then move the `fxstat` executable to a directory of your choice like `$HOME/bin/`.

If you are on a Unix system, all the requirements should be satisfied as fxstat
depends only on gcc for compilation and standard C libraries.

Usage
=====

Usage should be straightforward. Assuming `fxstat` is on your `PATH` (if not, use `path/to/fxstat`):

```
fxstat -h
fxstat reads.fastq
cat reads.fastq | fxstat
```

Unit tests and code coverage
============================

To run tests:

```
cd test
python test.py
```

Code coverage report will be in [test/code_coverage](https://htmlpreview.github.io/?https://github.com/dariober/fxstat/blob/master/test/code_coverage/index.html).

Requirements for running the test suite:

* [python](https://www.python.org/) should be available on virtually every system

* [valgrind](https://valgrind.org/) for checking memory leaks. If not
  available, install it with any package manager like `apt` or
  [conda](https://anaconda.org/conda-forge/valgrind)

* [lcov](https://github.com/linux-test-project/lcov) for code coverage report.
  As above, install with a package manager or
  [conda](https://anaconda.org/conda-forge/lcov)

There are tests to check for memory leaks using valgrind. Please check 

To compute code coverage you need
[lcov](https://github.com/linux-test-project/lcov). It can be installed
with your package manager (e.g., for Ubuntu `sudo apt-get install lcov`) or with
[conda](https://anaconda.org/conda-forge/lcov).
