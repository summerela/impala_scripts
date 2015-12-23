# impala connection strings
impala_host = 'glados18'
impala_port = '21050'
impala_user = 'selasady'

# Database to test
input_db = "p7_product"
input_table = "ens_partitioned"

# enter rsID to use as test
test_rsid = 'rs334'

####################
## import modules ##
####################
import pandas as pd
from impala.dbapi import connect
import unittest
import urllib, json

# disable extraneous pandas warning
pd.options.mode.chained_assignment = None

#################
# Unit Testing ##
#################

class table_test(unittest.TestCase):

        # establish connection to database
        def setUp(self):
            db =  None
            self.db = connect(host=impala_host, port=impala_port, timeout=10000, user=impala_user)

        # close database connection when done
        def tearDown(self):
            if self.db: self.db.close()


        def test_chromosomes(self):
            '''
            Test that all the chromosomes are loaded and that we
            have a reasonable number of rows for each
            '''
            print "Testing that all chromosomes were uploaded... \n"
            expect = map( str, range(1,23) ) + ['X','Y','MT', 'M']
            found  = []
            with self.db.cursor() as c:
                chrom_cmd = 'SELECT chrom, COUNT(1) AS cnt FROM {}.{} GROUP BY chrom'.format(input_db, input_table)
                c.execute(chrom_cmd)
                for r in c:
                    found.append( r[0] )
                    if r[0] == 'MT':
                        self.assertGreater( r[1], 1000, 'Chromosome has at least # variants' )
                    else:
                        self.assertGreater( r[1], 100000, 'Chromosome has at least # variants' )
            self.assertItemsEqual( found, expect, 'Chromosomes all loaded'  )


        def test_ref_alt(self):

            '''
            Test that the ref and alt columns all have values
            '''
            print "Testing that there are no null values in ref and alt columns... \n"
            colvals_cmd = 'SELECT COUNT(1) FROM {}.{} WHERE ref IS NULL OR alt IS NULL'.format(input_db, input_table)
            with self.db.cursor() as c:
                c.execute( colvals_cmd )
                self.assertEqual( c.fetchone()[0], 0,
                                  'One or more null values found in ref/alt columns' )


        # def test_one_based(self):
        #     '''
        #     Ensure that the data was loaded as 1-based not 0-based by checking the
        #     position of a common SNP, rs334 (Sickle Cell Anemia).
        #     '''
        #     print "Testing that coordinates are 1-based... \n"
        #     pos_cmd = 'SELECT chrom, pos FROM {}.{} WHERE rs_id = "rs334"'.format(input_db, input_table)
        #     with self.db.cursor() as c:
        #         c.execute( pos_cmd )
        #         # Position will be dependent on reference build
        #         # for build GRCh37.p13
        #         self.assertEqual( c.fetchone()[1], 5248232, "1-based coordinates" )
        #         self.assertEqual( c.fetchone()[0], "11", "1-based chromosome" )

        def test_one_based(self):
            '''
            Ensure that the data was loaded as 1-based not 0-based by checking the
            position of a common SNP
            '''
            print "Testing that coordinates are 1-based... \n"
            test_url = 'http://myvariant.info/v1/query?q={}'.format(test_rsid)
            response = urllib.urlopen(test_url)
            data = json.loads(response.read())
            test_pos = data["hits"][0]["dbnsfp"]["hg19"]["start"]
            test_chrom = data["hits"][0]["dbsnp"]["chrom"]
            pos_cmd = 'SELECT chrom, pos FROM {}.{} WHERE rs_id = "{}"'.format(input_db, input_table, test_rsid)
            with self.db.cursor() as c:
                c.execute( pos_cmd )
                # Position will be dependent on reference build
                # for build GRCh37.p13
                self.assertEqual( c.fetchone()[1], int(test_pos), "1-based coordinates" )
                self.assertEqual( c.fetchone()[0], test_chrom, "1-based chromosome" )


if __name__ == '__main__':
    unittest.main()