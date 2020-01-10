<!-- vim-markdown-toc GFM -->

* [Compile](#compile)
* [Usage](#usage)
* [Run unit tests and compute test code coverage](#run-unit-tests-and-compute-test-code-coverage)

<!-- vim-markdown-toc -->

Fast and easy stats for FASTA and FASTQ files

Compile
=======

```
gcc -g src/fxstat.c src/utils.c -o fxstat
```

Usage
=====

Usage should be straightforward:

```
fxstat -h
fxstat reads.fastq
cat reads.fastq | fxstat
```

Run unit tests and compute test code coverage
=============================================

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
