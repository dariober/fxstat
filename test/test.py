#!/usr/bin/env python3

import unittest
import shutil
import os
import re
import sys
import subprocess as sp

class Fxstat(unittest.TestCase):

    def setUp(self):
        sys.stderr.write('\n' + self.id().split('.')[-1] + ' ') # Print test name

    def valgrind(self, stderr):
        self.assertTrue('All heap blocks were freed' in stderr)
        self.assertTrue('ERROR SUMMARY: 0 errors' in stderr)
        

    def testHelp(self):
        p= sp.Popen('./fxstat -h',
                shell=True, stdout= sp.PIPE, stderr= sp.PIPE)
        stdout, stderr= p.communicate()
        self.assertEqual(0, p.returncode)
        self.assertTrue(len(stdout) > 20)

    def testVersion(self):
        p= sp.Popen('./fxstat -V',
                shell=True, stdout= sp.PIPE, stderr= sp.PIPE)
        stdout, stderr= p.communicate()
        self.assertEqual(0, p.returncode)
        self.assertTrue(len(stdout) >=5)

    def testBasicFastq(self):
        p= sp.Popen('./fxstat ../data/basic.fq',
                shell=True, stdout= sp.PIPE, stderr= sp.PIPE)
        stdout, stderr= p.communicate()
        self.assertEqual(0, p.returncode)
        self.assertTrue(re.search('n_sequences +9\n', stdout.decode()))
        self.assertTrue(re.search('mean_length +28.33\n', stdout.decode()))
        self.assertTrue(re.search('n_bases +255\n', stdout.decode()))

    def testCarriageReturn(self):
        p= sp.Popen('./fxstat ../data/basic_crlf.fq',
                shell=True, stdout= sp.PIPE, stderr= sp.PIPE)
        stdout, stderr= p.communicate()
        self.assertEqual(0, p.returncode)
        self.assertTrue(re.search('n_bases +255\n', stdout.decode()))

    def testMedian(self):
        p= sp.Popen('./fxstat ../data/basic.fq',
                shell=True, stdout= sp.PIPE, stderr= sp.PIPE)
        stdout, stderr= p.communicate()
        self.assertEqual(0, p.returncode)
        self.assertTrue(re.search('median_length +21.00\n', stdout.decode()))

        p= sp.Popen('cat ../data/basic.fq ../data/basic.fq | ./fxstat',
                shell=True, stdout= sp.PIPE, stderr= sp.PIPE)
        stdout, stderr= p.communicate()
        self.assertEqual(0, p.returncode)
        self.assertTrue(re.search('median_length +21.00\n', stdout.decode()))

        p= sp.Popen('head -n 16 ../data/basic.fq | ./fxstat',
                shell=True, stdout= sp.PIPE, stderr= sp.PIPE)
        stdout, stderr= p.communicate()
        self.assertEqual(0, p.returncode)
        self.assertTrue(re.search('median_length +31.50\n', stdout.decode()))

    def testBasicFasta(self):
        p= sp.Popen('./fxstat ../data/basic.fa',
                shell=True, stdout= sp.PIPE, stderr= sp.PIPE)
        stdout, stderr= p.communicate()
        self.assertEqual(0, p.returncode)
        self.assertTrue(re.search('n_sequences +7\n', stdout.decode()))
        self.assertTrue(re.search('n_bases +181\n', stdout.decode()))

    def testStopAfter(self):
        p= sp.Popen('./fxstat -s 5 ../data/basic.fq',
                shell=True, stdout= sp.PIPE, stderr= sp.PIPE)
        stdout, stderr= p.communicate()
        self.assertEqual(0, p.returncode)
        self.assertTrue(re.search('n_sequences +5\n', stdout.decode()))

    def testBaseQualityStats(self):
        p= sp.Popen('head -n 4 ../data/quality.fq | ./fxstat',
                shell=True, stdout= sp.PIPE, stderr= sp.PIPE)
        stdout, stderr= p.communicate()
        self.assertEqual(0, p.returncode)
        self.assertTrue(re.search('mean_read_quality +3.59\n', stdout.decode()))
    
    def testInvalidArg(self):
        # Unknown arg
        p= sp.Popen('./fxstat -Z',
                shell=True, stdout= sp.PIPE, stderr= sp.PIPE)
        stdout, stderr= p.communicate()
        self.assertEqual(1, p.returncode)

        # Option required
        p= sp.Popen('./fxstat -o',
                shell=True, stdout= sp.PIPE, stderr= sp.PIPE)
        stdout, stderr= p.communicate()
        self.assertEqual(1, p.returncode)

    def testNx(self):
        p= sp.Popen('./fxstat ../data/basic.fq',
                shell=True, stdout= sp.PIPE, stderr= sp.PIPE)
        stdout, stderr= p.communicate()
        self.assertEqual(0, p.returncode)
        self.assertTrue(re.search('N0 +60', stdout.decode()))
        self.assertTrue(re.search('N50 +42', stdout.decode()))
        self.assertTrue(re.search('N100 +9', stdout.decode()))

        p= sp.Popen('./fxstat -N 33,66,33 ../data/basic.fq',
                shell=True, stdout= sp.PIPE, stderr= sp.PIPE)
        stdout, stderr= p.communicate()
        self.assertEqual(0, p.returncode)
        self.assertTrue(re.search('N33 ', stdout.decode()))
        self.assertTrue(re.search('N66 ', stdout.decode()))
        # Check dups removed
        out = stdout.decode().split('\n')
        self.assertEqual(1, len([x for x in out if x.startswith('N33')]))

        p= sp.Popen('./fxstat -N 101 ../data/basic.fq',
                shell=True, stdout= sp.PIPE, stderr= sp.PIPE)
        stdout, stderr= p.communicate()
        self.assertEqual(1, p.returncode)

        p= sp.Popen('./fxstat -N -1 ../data/basic.fq',
                shell=True, stdout= sp.PIPE, stderr= sp.PIPE)
        stdout, stderr= p.communicate()
        self.assertEqual(1, p.returncode)

    def testBasicFastqNucPercent(self):
        p= sp.Popen('./fxstat ../data/basic.fq',
                shell=True, stdout= sp.PIPE, stderr= sp.PIPE)
        stdout, stderr= p.communicate()
        self.assertEqual(0, p.returncode)
        self.assertTrue(re.search('A +35.69', stdout.decode()))
        self.assertTrue(re.search('GC +37.25', stdout.decode()))
        self.assertTrue(re.search('N +1.18', stdout.decode()))

    def testLongRead(self):
        p= sp.Popen('./fxstat ../data/long.fq',
                shell=True, stdout= sp.PIPE, stderr= sp.PIPE)
        stdout, stderr= p.communicate()
        self.assertEqual(0, p.returncode)

    def testReadStdin(self):
        p= sp.Popen('cat ../data/basic.fq | ./fxstat', shell=True, stdout= sp.PIPE, stderr= sp.PIPE)
        stdout, stderr= p.communicate()
        self.assertEqual(0, p.returncode)
        
        p= sp.Popen('cat ../data/basic.fq | ./fxstat -', shell=True, stdout= sp.PIPE, stderr= sp.PIPE)
        stdout, stderr= p.communicate()
        self.assertEqual(0, p.returncode)

        p= sp.Popen('./fxstat <(cat ../data/basic.fq) <(gzip -c ../data/basic.fq)', shell=True, stdout= sp.PIPE, stderr= sp.PIPE, executable='/bin/bash')
        stdout, stderr= p.communicate()
        self.assertEqual(0, p.returncode)
        self.assertTrue('n_files: 2', stdout.decode())

    def testGzipStdin(self):
        p= sp.Popen('gzip -c ../data/long.fq | ./fxstat',
                shell=True, stdout= sp.PIPE, stderr= sp.PIPE)
        stdout, stderr= p.communicate()
        self.assertEqual(0, p.returncode)
        self.assertTrue(re.search('n_bases +40000\n', stdout.decode()))

    def testGzipUnexpectedEOF(self):
        p= sp.Popen('gzip -c ../data/long.fq | head -c 100 | ./fxstat',
                shell=True, stdout= sp.PIPE, stderr= sp.PIPE)
        stdout, stderr= p.communicate()
        self.assertTrue(p.returncode != 0);
        self.assertEqual("", stdout.decode());

    def testWriteToFile(self):
        p= sp.Popen('./fxstat ../data/basic.fq -o out.txt',
                shell=True, stdout= sp.PIPE, stderr= sp.PIPE)
        stdout, stderr= p.communicate()
        self.assertEqual(0, p.returncode)
        self.assertTrue(os.path.getsize('out.txt') > 20);

        p= sp.Popen('./fxstat ../data/basic.fq -o -',
                shell=True, stdout= sp.PIPE, stderr= sp.PIPE)
        stdout, stderr= p.communicate()
        self.assertEqual(0, p.returncode)
        self.assertTrue(len(stdout.decode()) > 20)

        p= sp.Popen('./fxstat ../data/basic.fq -o foobar/out.txt',
                shell=True, stdout= sp.PIPE, stderr= sp.PIPE)
        stdout, stderr= p.communicate()
        self.assertTrue(p.returncode != 0)

    def testMemoryLeakFastQ(self):
        p= sp.Popen('valgrind --leak-check=full ./fxstat ../data/basic.fq',
                shell=True, stdout= sp.PIPE, stderr= sp.PIPE)
        stdout, stderr= p.communicate()
        self.assertEqual(0, p.returncode)
        self.valgrind(stderr.decode())

        p= sp.Popen('head -n 4 ../data/basic.fq | valgrind ./fxstat',
                shell=True, stdout= sp.PIPE, stderr= sp.PIPE)
        stdout, stderr= p.communicate()
        self.assertEqual(0, p.returncode)
        self.valgrind(stderr.decode())

        p= sp.Popen('cat ../data/long.fq ../data/long.fq ../data/long.fq | valgrind ./fxstat',
                shell=True, stdout= sp.PIPE, stderr= sp.PIPE)
        stdout, stderr= p.communicate()
        self.assertEqual(0, p.returncode)
        self.valgrind(stderr.decode())

    def testMemoryLeakFastA(self):
        p= sp.Popen('valgrind --leak-check=full ./fxstat ../data/basic.fa',
                shell=True, stdout= sp.PIPE, stderr= sp.PIPE)
        stdout, stderr= p.communicate()
        self.assertEqual(0, p.returncode)
        self.valgrind(stderr.decode())

        p= sp.Popen('head -n 4 ../data/basic.fq | valgrind ./fxstat',
                shell=True, stdout= sp.PIPE, stderr= sp.PIPE)
        stdout, stderr= p.communicate()
        self.assertEqual(0, p.returncode)
        self.valgrind(stderr.decode())

    def testCompileWithoutWarnings(self):
        p= sp.Popen('gcc -O3 -Wall ../../src/fxstat.c -lz -lm -o test_fxstat',
                shell=True, stdout= sp.PIPE, stderr= sp.PIPE)
        stdout, stderr= p.communicate()
        self.assertEqual('', stderr.decode())

    def testZeroLengthRead(self):
        p= sp.Popen('./fxstat -N 0,100 ../data/zerolength.fq',
                shell=True, stdout= sp.PIPE, stderr= sp.PIPE)
        stdout, stderr= p.communicate()
        self.assertEqual(0, p.returncode)
        self.assertTrue(re.search('n_sequences +3\n', stdout.decode()))
        self.assertTrue(re.search('mean_length +10.00\n', stdout.decode()))
        # NB: N100 reports the shortest read with length >0
        self.assertTrue(re.search('N100 +10', stdout.decode()))
        self.assertTrue(re.search('N0 +20', stdout.decode()))

    def testMultipleFiles(self):
        p= sp.Popen('./fxstat ../data/basic.fq ../data/basic.fa ../data/basic.fq',
                shell=True, stdout= sp.PIPE, stderr= sp.PIPE)
        stdout, stderr= p.communicate()
        self.assertEqual(0, p.returncode)
        self.assertTrue("[../data/basic.fq]" in stdout.decode())
        self.assertTrue("[../data/basic.fa]" in stdout.decode())
        self.assertEqual(1, stdout.decode().count('basic.fq')) # Dup file removed
        self.assertTrue(re.search('n_bases +181', stdout.decode()))
        self.assertTrue(re.search('n_bases +255', stdout.decode()))

    def testMultipleFilesPulled(self):
        p= sp.Popen('valgrind --leak-check=full ./fxstat -p ../data/quality.fq ../data/basic.fq',
                shell=True, stdout= sp.PIPE, stderr= sp.PIPE)
        stdout, stderr= p.communicate()
        self.assertEqual(0, p.returncode)
        self.valgrind(stderr.decode())

        self.assertTrue(re.search('n_sequences +17\n', stdout.decode()))
        self.assertTrue(re.search('n_bases +375', stdout.decode()))
        self.assertTrue("[../data/basic.fq]\n[../data/quality.fq]" in stdout.decode())
        
        # These two executions are equivalent:        
        p= sp.Popen("./fxstat -p ../data/long.fq ../data/quality.fq ../data/basic.fq | grep -v -F '[' | grep -v n_files",
                shell=True, stdout= sp.PIPE, stderr= sp.PIPE)
        stdout, stderr= p.communicate()

        pexp= sp.Popen("cat ../data/basic.fq ../data/quality.fq ../data/long.fq | ./fxstat | grep -v -F '[' | grep -v n_files",
                shell=True, stdout= sp.PIPE, stderr= sp.PIPE)
        stdoutExp, stderr2= pexp.communicate()
        self.assertEqual(stdoutExp.decode(), stdout.decode())


    def testMultipleFilesPulledStopAfter(self):
        p= sp.Popen('valgrind --leak-check=full ./fxstat -s 5 -p ../data/quality.fq ../data/basic.fq',
                shell=True, stdout= sp.PIPE, stderr= sp.PIPE)
        stdout, stderr= p.communicate()
        self.assertEqual(0, p.returncode)
        self.valgrind(stderr.decode())

        self.assertTrue(re.search('n_sequences +5\n', stdout.decode()))
        self.assertTrue("[../data/basic.fq]" in stdout.decode()) # We don't read the 2nd file at all

        p= sp.Popen('valgrind --leak-check=full ./fxstat -s 10 -p ../data/quality.fq ../data/basic.fq',
                shell=True, stdout= sp.PIPE, stderr= sp.PIPE)
        stdout, stderr= p.communicate()
        self.assertEqual(0, p.returncode)
        self.valgrind(stderr.decode())

        self.assertTrue(re.search('n_sequences +10\n', stdout.decode()))
        self.assertTrue("[../data/basic.fq]\n[../data/quality.fq]" in stdout.decode())

    def testInputDir(self):
        p= sp.Popen('valgrind --leak-check=full ./fxstat ../data ../data',
                shell=True, stdout= sp.PIPE, stderr= sp.PIPE)
        stdout, stderr= p.communicate()
        self.assertEqual(0, p.returncode)
        self.valgrind(stderr.decode())
        self.assertEqual(9, stdout.decode().count('['));
        self.assertTrue('n_files: 9' in stdout.decode());
        self.assertTrue('.txt' not in stdout.decode());

    def testFileNotfound(self):
        p= sp.Popen('./fxstat foo.fq',
                shell=True, stdout= sp.PIPE, stderr= sp.PIPE)
        stdout, stderr= p.communicate()
        self.assertEqual(1, p.returncode)
        self.assertTrue("Invalid input file" in stderr.decode());

    def testOntStats(self):
        p= sp.Popen('valgrind ./fxstat ../data/ont.fq',
                shell=True, stdout= sp.PIPE, stderr= sp.PIPE)
        stdout, stderr= p.communicate()
        self.valgrind(stderr.decode())

        self.assertTrue(re.search('ont_time_span +203', stdout.decode()))
        self.assertTrue('00:03:23' in stdout.decode())

        self.assertTrue(re.search('ont_n_channels +10', stdout.decode()))
        self.assertTrue(re.search('ont_ch_mean_reads +1.80', stdout.decode()))
        self.assertTrue(re.search('ont_ch_sd_reads +1.03', stdout.decode()))
        self.assertTrue(re.search('ont_ch_dispersion +0.59', stdout.decode()))

        p= sp.Popen('valgrind ./fxstat -p ../data/53oxa1_fast_nano.fq ../data/ont.fq',
                shell=True, stdout= sp.PIPE, stderr= sp.PIPE)
        stdout, stderr= p.communicate()
        self.valgrind(stderr.decode())
        self.assertTrue(re.search('ont_n_channels +187', stdout.decode()))

    def testUnrecognizedInputError(self):
        pass

if __name__ == '__main__':
    # Set up
    # First off: compile code and check ok. -g is to improve valgrind report
    # compiled code and test output will be in 'out'
    if os.path.exists('out'):
        shutil.rmtree('out')
    os.makedirs('out')
    os.chdir('out') # NB we are in 'out'

    p= sp.Popen('gcc -fprofile-arcs -ftest-coverage -g ../../src/fxstat.c -lz -lm -o fxstat',
            shell=True, stdout= sp.PIPE, stderr= sp.PIPE)
    stdout, stderr= p.communicate()
    if p.returncode != 0:
        print(stderr.decode())
        sys.exit(1)
    print("Compiled OK")
    
    unittest.main(exit= False)
    
    # code coverage
    p = sp.Popen('lcov --help && genhtml --help', shell= True, stdout= sp.PIPE, stderr= sp.PIPE)
    p.communicate()
    if p.returncode == 0:
        try:
            os.symlink('../../src/fxstat.c', './fxstat.c')
            p= sp.Popen('gcov fxstat.c', shell= True, stdout= sp.PIPE, stderr= sp.PIPE)
            p.communicate()
            p= sp.Popen('lcov --capture --directory . --output-file gcov.info', shell= True, stdout= sp.PIPE, stderr= sp.PIPE)
            p.communicate()
            p= sp.Popen('genhtml gcov.info --output-directory ../code_coverage', shell= True, stdout= sp.PIPE, stderr= sp.PIPE)
            p.communicate()
        except:
            print('Error performing code coverage')
    else:
        print('lcov and/or genhtml not available. Code coverage report skipped')
    shutil.rmtree('../out')
