"""Test rejection outlier using built-in grcollect executable."""

from os import path
import sys
import subprocess
from subprocess import Popen, PIPE
from functools import partial

import unittest
import random
import statistics
from pandas import read_csv
import pandas

#TODO need for fix import grcollect path for windows; points to c:/local
#from astrowisp import grcollect_path
from astrowisp.tests.utilities import FloatTestCase

#_test_data_dir = path.join(path.dirname(path.abspath(__file__)),
#                           'test_data')

#TODO remove manual overwrite of directory paths here
_test_data_dir='B:/Github/AstroWISP/PythonPackage/astrowisp/tests/test_data/'
grcollect_path='B:/fitsh-0.9.4-minimum/builddir/src/grcollect.exe'

class TestGRCollect(FloatTestCase):
    """Test cases for the grcollect executable."""

    def test_dataset(self):
        """Create dataset to test rejection outlier"""
        # parse_result = partial(read_csv, sep=r'\s+', comment='#', header=None)

        #setup subprocess for grcollect
        #TODO do we need to use median or mean rejection for grcollect
        grcollect = Popen(
                [
                grcollect_path,
                '-',
                '-V',
                '--col-base', '1',
                '--col-stat', '2,3,4,5,6,7,8,9,10,11,12,13,14',
                '--output', str(path.join(_test_data_dir, 'test_stat_grcollect.txt')),
                '--stat', 'count,rcount,mean,rmean,median,rmedian,mediandev,rmediandev',
                '--rejection', 'column=4,iterations=10,mean,stddev=4.0',
                # '--rejection', 'column=4,,iterations=10,median,stddev=4.0',
                '--rejection', 'column=7,iterations=10,mean,stddev=4.0',
                # '--rejection', 'column=7,,iterations=10,median,stddev=4.0',
                '--rejection', 'column=10,iterations=10,mean,meddev=15.0',
                # '--rejection', 'column=10,,iterations=10,median,stddev=4.0',
                '--rejection', 'column=13,iterations=10,mean,meddev=15.0',
                # '--rejection', 'column=13,,iterations=10,median,stddev=4.0',
                '--max-memory', '4g',
                #TODO add tmpdir and remove after usage to get rid of all the grcollect temp files
                #temp file module python
                #'--tmpdir', ''
                ],
                stdin=PIPE,
                stdout=PIPE
                )
        
        col_gaussian=[]
        mu = 100
        sigma = 50
        #TODO maybe do 1000000, although it takes like an hour
        niter=200000
        #TODO add gaussian to inner iteration
        for i in range(niter):  
            temp = random.gauss(mu, sigma)  
            col_gaussian.append(temp)

        stdev_gaus = statistics.stdev(col_gaussian)
        mean_gaus = statistics.mean(col_gaussian)

        #setup outlier counters
        one_outliers=0
        onefloat_outliers=0
        mil_outliers=0
        gaus_outliers=0

        #array values to add random statistics to
        #every 3rd version of values is being altered i.e. values[3n]
        values = [0, 
                  1, 1, 1, 
                  1.0, 1.0, 1.0, 
                  1.0, 1.0, 1.0, 
                  1.0, 1.0, 1.0, 
                  0]

        #TODO maybe do 1000000, although it takes like an hour
        for i in range(niter):
            values[0] = random.randint(0, 10)

            #TODO need to find how to exclude the outliers so we can match rmedian and rmean to 
            # the regular mean and median excluded columns
            if i%20==0:
                values[2] = 0 #100
                values[3] = 100
                one_outliers += 1
            else:
                values[2] = 1
                values[3] = 1
            
            if i%20==0:
                values[5] = 0 #100.0
                values[6] = 100.0
                onefloat_outliers += 1
            else:
                values[5] = 1.0
                values[6] = 1.0

            values[7] += 1

            if i%20==0:
                values[8] = 0 #values[7] + 100000000
                values[9] = values[7] + 100000000
                mil_outliers += 1
            elif i%50==0:
                values[8] = 0 #values[7] - 100000000
                values[9] = values[7] - 100000000
                mil_outliers += 1
            else:
                values[8] = values[7]
                values[9] = values[7]

            values[10] = col_gaussian[i]  

            if i%20==0:
                values[11] = 0 #values[10] + 3*100000000*stdev_gaus
                values[12] = values[10] + 3*100000000*stdev_gaus
                gaus_outliers += 1
            elif i%50==0:
                values[11] = 0 #values[10] - 3*100000000*stdev_gaus
                values[12] = values[10] - 3*100000000*stdev_gaus
                gaus_outliers += 1
            else:
                values[11] = values[10]
                values[12] = values[10]
            
            #writes each value list as a line in subprocess POPEN for grcollect
            grcollect.stdin.write(((str(' '.join(str(v) for v in values))+'\n')).encode('ascii'))

        #makes subprocess wait for grcollect to finish
        grcollect.communicate()
        #closes subprocess
        grcollect.stdout.close()

        #TODO remove print statements for everything from here forward
        print('outliers:')
        print('one_outliers='+repr(one_outliers))
        print('onefloat_outliers='+repr(onefloat_outliers))
        print('mil_outliers='+repr(mil_outliers))
        print('gaus_outliers='+repr(gaus_outliers))

        #write grcollect stat file to pandas dataframe
        df=pandas.read_csv(path.join(_test_data_dir, 'test_stat_grcollect.txt'), 
                           header=None, 
                           sep=r'\s+', 
                           index_col=0)
        #TODO disable to csv after grcollect is tested
        df.to_csv(path.join(_test_data_dir, 'test_stat_grcollect.csv'))

        #every 24 columns is the columns we perform rejection to in this case
        #columns are 'count,rcount,mean,rmean,median,rmedian,mediandev,medianmeddev' 
        #for every column that has statistics added to grcollect
        #i.e columns '2,3,4,5,6,7,8,9,10,11,12,13,14'

        print((df[17]-df[18]).sum())
        print((df[17+24]-df[18+24]).sum())
        print((df[17+24+24]-df[18+24+24]).sum())
        print((df[17+24+24+24]-df[18+24+24+24]).sum())



        print((df[17]-df[18]).sum())
        #self assert true if every rejection outlier column from grcollect matches the number of outliers we produced
        #TODO cleanup this +24 junk to something nicer maybe do something with range(24,3,24) and iterate over those
        self.assertTrue((df[17]-df[18]).sum()==one_outliers,
                        f'The one outliers mismatch: {(df[17+24+24]-df[18+24+24]).sum()!r} not equal to '
                        f'{one_outliers!r}')
        self.assertTrue((df[17+24]-df[18+24]).sum()==onefloat_outliers,
                        f'The floating point one outliers mismatch: {(df[17+24+24]-df[18+24+24]).sum()!r} not equal to '
                        f'{onefloat_outliers!r}')
        self.assertTrue((df[17+24+24]-df[18+24+24]).sum()==mil_outliers,
                        f'The counted outliers mismatch: {(df[17+24+24]-df[18+24+24]).sum()!r} not equal to '
                        f'{mil_outliers!r}')

        #gaussian outlier rejection should be within 2 percent of its total outliers
        self.assertTrue(gaus_outliers-0.02*gaus_outliers<=
                        (df[17+24+24+24]-df[18+24+24+24]).sum()
                        <=gaus_outliers+0.02*gaus_outliers,
                        f'Gaussian outliers mismatch: {(df[17+24+24+24]-df[18+24+24+24]).sum()!r} not within 2% of '
                        f'{gaus_outliers!r}')
        
        #now check if the rmeans and rmedians should be what they should if the outliers werent present
        #columns are 'count,rcount,mean,rmean,median,rmedian,mediandev,medianmeddev' 

        #TODO need to find how to exclude the outliers so we can match rmedian and rmean to 
        # the regular mean and median excluded columns
        #maybe doing a antimean then applying a new mean with the outlier number but having 0s for the 
        #outliers for the dataset that dont get rejected
        # print(df[18-8+2].sum()*(niter)/(niter-one_outliers))
        # print(df[18+2].sum())
        print('means with excluded')
        print(df[18-8+2].sum()*(niter)/(niter-one_outliers))
        print(df[18+2].sum())
        print(df[18-8+2+24].sum()*(niter)/(niter-onefloat_outliers))
        print(df[18+2+24].sum())
        print(df[18-8+2+24+24].sum()*(niter)/(niter-mil_outliers))
        print(df[18+2+24+24].sum())
        print(df[18-8+2+24+24+24].sum()*(niter)/(niter-gaus_outliers))
        print(df[18+2+24+24+24].sum())

        print('means')
        print(df[18-8+2].sum())
        print(df[18+2].sum())
        print(df[18-8+2+24].sum())
        print(df[18+2+24].sum())
        print(df[18-8+2+24+24].sum())
        print(df[18+2+24+24].sum())
        print(df[18-8+2+24+24+24].sum())
        print(df[18+2+24+24+24].sum())

        # print('means')
        # print(df[18-8+2])
        # print(df[18+2])
        # print(df[18-8+2+24])
        # print(df[18+2+24])
        # print(df[18-8+2+24+24])
        # print(df[18+2+24+24])
        # print(df[18-8+2+24+24+24])
        # print(df[18+2+24+24+24])

        print('medians')
        print(df[18-8+4].sum())
        print(df[18+4].sum())
        print(df[18-8+4+24].sum())
        print(df[18+4+24].sum())
        print(df[18-8+4+24+24].sum())
        print(df[18+4+24+24].sum())
        print(df[18-8+4+24+24+24].sum())
        print(df[18+4+24+24+24].sum())

if __name__ == '__main__':
    unittest.main()
