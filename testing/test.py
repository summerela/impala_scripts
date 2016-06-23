import unittest



############################
### Unit Test with Spark ###
############################

class TestETL(unittest.TestCase):

    def __init__(self):
        self.somatic_chr = sorted(map(str, range(1, 23) + ["M", "X", "Y"]))
        self.test_chroms = sorted(map(str, range(1, 20) + ["M", "X", "Y"]))

    def check_chrom_set(self):
        assert len(set(self.somatic_chr).difference(self.test_chroms)) == 0,"The following chromosomes are not loaded: \n {}".format(set(self.somatic_chr).difference(self.test_chroms))

    def chrom_check2(self):
        self.assertListEqual(self.somatic_chr, self.test_chroms, "The following chromosomes are not loaded: \n {}".format(set(self.somatic_chr).difference(self.test_chroms)))

###############
if __name__ == '__main__':

    test = TestETL()

    # test.check_chrom_set()

    test.chrom_check2()

