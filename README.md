![Language](https://img.shields.io/badge/language-C-blue.svg)
[![Build Status](https://travis-ci.com/dariober/fxstat.svg?branch=master)](https://travis-ci.com/dariober/fxstat)
[![Test coverage](https://img.shields.io/badge/coverage-gcov-grey.svg)](https://htmlpreview.github.io/?https://github.com/dariober/fxstat/blob/master/test/code_coverage/index.html)

<!-- vim-markdown-toc GFM -->

* [Compile and setup](#compile-and-setup)
* [Usage & Output](#usage--output)
* [Performance](#performance)
* [Unit tests and code coverage](#unit-tests-and-code-coverage)
* [Known issues](#known-issues)

<!-- vim-markdown-toc -->

Fast and easy stats for FASTA and FASTQ files

Compile and setup
=================

Assuming you are in the root directory of this project (i.e. where this README
file is), compile the fxstat executable with:

```
gcc -O3 -Wall -g src/fxstat.c -lz -lm -o fxstat
```

Then move the `fxstat` executable to a directory of your choice like `$HOME/bin/`.

If you are on a Unix system, all the requirements should be satisfied as fxstat
depends only on gcc for compilation and on standard C libraries.

Usage & Output
==============

Usage should be straightforward. Assuming `fxstat` is on your `PATH` (if not,
use `path/to/fxstat`):

```
fxstat -h
fxstat r1.fastq.gz r2.fastq.gz
fxstat sequence_dir/
cat reads.fastq | fxstat
```

Stat all `*.fastq` files in directory `fastq_pass` and in its subdirectories,
pull all statistics in a single summary (default is one summary per file):

```
fxstat -p fastq_pass
```

For more control over files to include, consider using the `find` command. *E.g.*,

```
fxstat -p `find fastq_pass -name '*.fastq'
```

----

The output should be self explanatory:

```
[53oxa1_fast_nano.fastq.gz]
n_seq              35028
n_bases            365454001
mean_length        10433.20
median_length      8224.50
mean_read_quality  10.40
N0                 68226  # MAX length
N25                23612
N50                17594
N75                11909
N100               57  # MIN length > 0
A                  24.87%
T                  24.76%
C                  25.05%
G                  25.32%
N                  0.00%
GC                 50.37%
n_files: 1
# Proc time 00:00:07
```

The [N50](https://en.wikipedia.org/wiki/N50,_L50,_and_related_statistics) is
the sequence length of the shortest read at 50% of the total read length. Other
Nx statistics are defined similarly for other thresholds.

Performance
===========

Speed and ease of use should be the main features of `fxstat`. Quick & dirty
benchmark on v0.1.0 commit 8c5290c. Test file is 370 Mbytes (compressed), 365
million bases:

```
time zcat 53oxa1_fast_nano.fastq.gz | wc > /dev/null

real    0m11.941s
user    0m20.965s
sys     0m0.395s

time fxstat 53oxa1_fast_nano.fastq.gz > /dev/null

real    0m7.572s
user    0m7.486s
sys     0m0.085s
```

Test file from:

```
curl -s https://sra-download.ncbi.nlm.nih.gov/traces/sra54/SRZ/011038/SRR11038963/53oxa1_fast_nano.fastq \
| gzip > 53oxa1_fast_nano.fastq.gz
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

Known issues
============

* There is no check the input files are correctly formatted
