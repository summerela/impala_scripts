#!/usr/bin/env python

import argparse
import sys
import os
from impala.dbapi import connect
from impala.util import as_pandas

# connect to impala db
conn=connect(host='local', port=21050)

#################
### User Args ###
#################
# parse input family id's either in text or as a list
def parse_family_arg(parser, arg):
    if os.path.exists(arg):
        with open(arg) as f:
            for line in f:
                fam_list = (line.strip('\n').split(','))
                return fam_list
    elif "," in arg:
        print("Input is a list!")
    else:
        raise SystemExit("Enter comma separated list of family id's or the path to a .txt file containing them")

def get_parser():

    try:
        parser = argparse.ArgumentParser(description="Find variants located in ACMG gene regions")
        parser.add_argument("--kav_freq", metavar='-K', help="enter max kaviar freq percent",
                            default=.03, type=float)
        parser.add_argument("--hdfs_path", metavar='-H', help="enter HDFS path to output result table",
                            default='default', type=str)
        parser.add_argument("--platform", metavar='-P', help="specify either illumina or cgi platform",
                            default='illumina', type=str)
        parser.add_argument("fam", metavar='F', help="family id of trios to analyze as comma separated list or path to a .txt file containing this information as a comma-separated list",
                            type=lambda x: parse_family_arg(parser, x))
        args = parser.parse_args()
        return args
    except:
        e = sys.exc_info()[0]
        print (e)

# get variants that match family id's
    # subject_id, family_id,  var_id, gt, chrom, pos
    # return only trios and add m/f/n annotation (sequential if more than one nb)
    # return only het parents
    # exclude X,Y,MT chromosomes

# match with acmg genes
    # by chrom, start, stop
    # return phenotype, mim_disorder, age_onset, inheritance, variants_to_report
    # return only variants in ACMG gene regions
    # examine inheritance for newborns (add label in query with case statement)
        # dominant = (inheritance = AD)
        # hom_recessive (inheritance = AR & newborn gt = 1/1)

# annotate with global vars table
    # join by var_id
    # return kav_freq, clin_sig, cadd_raw, dann_score, interpro_domain, gene_name, transcript_name, transcript_id,
    #    exon_name, exon_number, effect, impact, feature, hgvs_c, hgvs_p, ppc_rating
    # filter by must be kaviar rare and (impact = high or ppc=synonymous)


# label compound heterozygotes using python
    # label possible de novo or toss out?

# return dominant, hom_recessive and comp het variants that pass kav and pathogenic filters

# any other reports/visualizations?










if __name__ == "__main__":
    args = get_parser()
    print (args)

