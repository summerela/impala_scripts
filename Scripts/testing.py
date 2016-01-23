#!/usr/bin/env python

import sys
import collections as coll

# specify input file
in_file = '/Users/summerrae/Documents/dreamworld.txt'

# parse command line if no input file given
if in_file:
    inf = in_file
else:
    inf = sys.stdin

with open(inf, 'r') as f:
    p = f.read()
    words = p.split()
    word_count = len(words)
    print "The total word count is {}. ".format(word_count)
    sentences = p.split('.')
    sentence_count = len(sentences)
    print "There are {} sentences.".format(sentence_count)
    avg_len = sum(len(x.split()) for x in sentences) / sentence_count
    print "The average length of each sentance is {}.".format(avg_len)
    distinct_words = coll.Counter(p.split())
    print "The count of distinct words = ", distinct_words
    unique_words = [word for (word, count) in distinct_words.iteritems() if count==1]
    print "There are {} words that occur only once.".format(len(unique_words))




