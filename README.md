<!-- vim-markdown-toc GFM -->

* [Compile](#compile)
* [Run unit tests and compute test code coverage](#run-unit-tests-and-compute-test-code-coverage)

<!-- vim-markdown-toc -->

Fast and easy stats for FASTA and FASTQ files

Compile
=======

```
gcc -g src/fxstat.c src/utils.c -o fxstat
```

Run unit tests and compute test code coverage
=============================================

Run test:

```
cd test
python test.py
```

Code coverage report will be in `code_coverage`.

To compute code coverage you need
[lcov](https://github.com/linux-test-project/lcov). It can be installed
with your package manager (e.g. `sudo apt-get install lcov`) or with
[conda](https://anaconda.org/conda-forge/lcov).
