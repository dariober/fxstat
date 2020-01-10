#!/usr/bin/env python3

import unittest
import shutil
import os
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
        self.assertTrue('n_seq\t9\n' in stdout.decode())
        self.assertTrue('mean_length\t28.33\n' in stdout.decode())
        self.assertTrue('n_bases\t255\n' in stdout.decode())

    def testBaseQualityStats(self):
        p= sp.Popen('head -n 4 ../data/quality.fq | ./fxstat',
                shell=True, stdout= sp.PIPE, stderr= sp.PIPE)
        stdout, stderr= p.communicate()
        self.assertEqual(0, p.returncode)
        self.assertTrue('mean_read_quality\t3.59\n' in stdout.decode())

    def testBasicFastqNucPercent(self):
        p= sp.Popen('./fxstat ../data/basic.fq',
                shell=True, stdout= sp.PIPE, stderr= sp.PIPE)
        stdout, stderr= p.communicate()
        self.assertEqual(0, p.returncode)
        self.assertTrue('A\t35.69' in stdout.decode())
        self.assertTrue('CG\t37.25' in stdout.decode())
        self.assertTrue('N\t1.18' in stdout.decode())

    def testLongRead(self):
        p= sp.Popen('./fxstat ../data/long.fq',
                shell=True, stdout= sp.PIPE, stderr= sp.PIPE)
        stdout, stderr= p.communicate()
        self.assertEqual(0, p.returncode)

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
        p= sp.Popen('gcc -O3 -Wall ../../src/fxstat.c ../../src/utils.c -o test_fxstat',
                shell=True, stdout= sp.PIPE, stderr= sp.PIPE)
        stdout, stderr= p.communicate()
        self.assertEqual('', stderr.decode())

    def testZeroLengthRead(self):
        pass 

if __name__ == '__main__':
    # Set up
    # First off: compile code and check ok. -g is to improve valgrind report
    # compiled code and test output will be in 'out'
    if os.path.exists('out'):
        shutil.rmtree('out')
    os.makedirs('out')
    os.chdir('out') # NB we are in 'out'

    p= sp.Popen('gcc -fprofile-arcs -ftest-coverage -g ../../src/fxstat.c ../../src/utils.c -o fxstat',
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
