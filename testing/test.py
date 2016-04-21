#!/usr/bin/env python

# read in file
in_file = "/Users/selasady/Downloads/demo.txt"

import csv

f = open(in_file)
csv_f = csv.reader(f)

for row in csv_f:
    rowtotal = 0
    vals = [row[3].split('|')[1], row[4:]]
    row_sum = 0
    for val in vals:
        sum += val
    print row_sum
        # print "{}:Average:{}".format(row[0], sum(vals))



