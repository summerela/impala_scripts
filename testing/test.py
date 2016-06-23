#!/usr/bin/env pyspark

'''
Documemtation here
'''

# import the python packages you will need
import os

# create a class to hold the functions we'll use
class test_class():

    # define variables that won't change here
    days_in_year = 365

    # initialize these variables every time when class is instantiated
    def __init__(self, in_file, out_file):
        self.in_file = in_file
        self.out_file = out_file

        # define functions here
        def test_function(self):
            print ("We will process this file {}").format(self.in_file)

###############
# tell python what to do when this program is run
if __name__ == '__main__':

    my_file = 'get/from/here.txt'
    my_out_file = 'put/this/here.txt'

    # instantiate the class, leave at default or input user args
    test = test_class(my_file, my_out_file)


