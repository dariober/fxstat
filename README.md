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
gcc -g src/fxstat.c src/utils.c -o fxstat
```

Then move the `fxstat` executable to a directory of your choice like `$HOME/bin/`.

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

Run tests:

```
cd test
python test.py
```

Code coverage report will be in [code_coverage](https://htmlpreview.github.io/?https://github.com/dariober/fxstat/blob/master/code_coverage/index.html).

To compute code coverage you need
[lcov](https://github.com/linux-test-project/lcov). It can be installed
with your package manager (e.g., for Ubuntu `sudo apt-get install lcov`) or with
[conda](https://anaconda.org/conda-forge/lcov).
