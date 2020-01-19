#!/usr/bin/env python3

import unittest
import shutil
import os
import re
import sys
import subprocess as sp

class Fqstat(unittest.TestCase):

    def setUp(self):
        sys.stderr.write('\n' + self.id().split('.')[-1] + ' ') # Print test name

    #def tearDown(self):
    #    if os.path.exists('test_out'):
    #        shutil.rmtree('test_out')

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
        self.assertTrue(re.search('n_seq +9\n', stdout.decode()))
        self.assertTrue(re.search('mean_length +28.33\n', stdout.decode()))
        self.assertTrue(re.search('n_bases +255\n', stdout.decode()))

    def testBasicFasta(self):
        p= sp.Popen('./fxstat ../data/basic.fa',
                shell=True, stdout= sp.PIPE, stderr= sp.PIPE)
        stdout, stderr= p.communicate()
        self.assertEqual(0, p.returncode)
        self.assertTrue(re.search('n_seq +7\n', stdout.decode()))
        self.assertTrue(re.search('n_bases +181\n', stdout.decode()))

    def testStopAfter(self):
        p= sp.Popen('./fxstat -s 5 ../data/basic.fq',
                shell=True, stdout= sp.PIPE, stderr= sp.PIPE)
        stdout, stderr= p.communicate()
        self.assertEqual(0, p.returncode)
        self.assertTrue(re.search('n_seq +5\n', stdout.decode()))

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
        self.assertTrue(re.search('N90 +9', stdout.decode()))

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
        self.assertTrue(re.search('CG +37.25', stdout.decode()))
        self.assertTrue(re.search('N +1.18', stdout.decode()))

    def testLongRead(self):
        p= sp.Popen('./fxstat ../data/long.fq',
                shell=True, stdout= sp.PIPE, stderr= sp.PIPE)
        stdout, stderr= p.communicate()
        self.assertEqual(0, p.returncode)

    def testGzipStdin(self):
        p= sp.Popen('gzip -c ../data/long.fq | ./fxstat',
                shell=True, stdout= sp.PIPE, stderr= sp.PIPE)
        stdout, stderr= p.communicate()
        self.assertEqual(0, p.returncode)
        self.assertTrue(re.search('n_bases +40000\n', stdout.decode()))

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
        self.assertTrue('All heap blocks were freed' in stderr.decode())
        self.assertTrue('ERROR SUMMARY: 0 errors' in stderr.decode())

        p= sp.Popen('head -n 4 ../data/basic.fq | valgrind ./fxstat',
                shell=True, stdout= sp.PIPE, stderr= sp.PIPE)
        stdout, stderr= p.communicate()
        self.assertEqual(0, p.returncode)
        self.assertTrue('All heap blocks were freed' in stderr.decode())
        self.assertTrue('ERROR SUMMARY: 0 errors' in stderr.decode())

        p= sp.Popen('cat ../data/long.fq ../data/long.fq ../data/long.fq | valgrind ./fxstat',
                shell=True, stdout= sp.PIPE, stderr= sp.PIPE)
        stdout, stderr= p.communicate()
        self.assertEqual(0, p.returncode)
        self.assertTrue('All heap blocks were freed' in stderr.decode())
        self.assertTrue('ERROR SUMMARY: 0 errors' in stderr.decode())

    def testMemoryLeakFastA(self):
        p= sp.Popen('valgrind --leak-check=full ./fxstat ../data/basic.fa',
                shell=True, stdout= sp.PIPE, stderr= sp.PIPE)
        stdout, stderr= p.communicate()
        self.assertEqual(0, p.returncode)
        self.assertTrue('All heap blocks were freed' in stderr.decode())
        self.assertTrue('ERROR SUMMARY: 0 errors' in stderr.decode())

        p= sp.Popen('head -n 4 ../data/basic.fq | valgrind ./fxstat',
                shell=True, stdout= sp.PIPE, stderr= sp.PIPE)
        stdout, stderr= p.communicate()
        self.assertEqual(0, p.returncode)
        self.assertTrue('All heap blocks were freed' in stderr.decode())
        self.assertTrue('ERROR SUMMARY: 0 errors' in stderr.decode())

    def testCompileWithoutWarnings(self):
        p= sp.Popen('gcc -O3 -Wall ../../src/fxstat.c -lz -o test_fxstat',
                shell=True, stdout= sp.PIPE, stderr= sp.PIPE)
        stdout, stderr= p.communicate()
        self.assertEqual('', stderr.decode())

    def testZeroLengthRead(self):
        pass 
    def testTruncatedGZip(self):
        pass
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

    p= sp.Popen('gcc -fprofile-arcs -ftest-coverage -g ../../src/fxstat.c -lz -o fxstat',
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
