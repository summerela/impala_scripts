#!/usr/bin/env python

import sys
import os.path
import datetime



################################################################################
#   Nov 17, 2011
#   Authors: Vlad Makarov, Chris Yoon
#   Language: Python
#   OS: UNIX/Linux, MAC OSX
#   Copyright (c) 2011, The Mount Sinai School of Medicine

#   Available under BSD  licence

#  Modified by Summer Elasady: accept direct command line arg, remove uppercase
#  requirement for column names
################################################################################

#  Purpose: Converts tabular format to VCF.
#  #Column names CHROM, POS, REF, ALT  are required for parsing
#  Input: Text file in the VCF format
#  Output: Tabular format with "parsed" extension


#  Arguments:
#  infile - name of tabular file

#  To run:
#  python tab2vcf.py
#
################################################################################
# script.py
import argparse

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument(
        'input_tsv', metavar='i', help='Enter the full file path to the tsv file')
    args = parser.parse_args()
 
ACCEPTED_CHR = ["1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", "16", "17", "18", "19", "20","21","22", "X", "Y", "MT"]

""" Linear search, returns first index """
def find_first_index(lst, elem):
    ind = 0
    for l in lst:
        if str(elem).strip() == str(l).strip():
            return ind
        ind=ind+1
    return -1

def vcfheader(filename):
    """ Generates VCF header """
    filename = os.path.basename(filename)
    filename = os.path.splitext(filename)[0]
    now = datetime.datetime.now()
    curdate=str(now.year)+'-'+str(now.month)+'-'+str(now.day)
    lines=[]
    lines.append('##fileformat=VCFv4.0')
    lines.append('##fileDate='+curdate)
    lines.append('##reference=1000Genomes-NCBI37')
    lines.append('#CHROM'+'\t'+'POS'+'\t'+'ID'+'\t'+'REF'+'\t'+'ALT'+'\t'+'QUAL'+'\t'+'FILTER'+'\t'+'INFO'+'\t'+'FORMAT' + '\t'+ 'Genotype')
    return '\n'.join(lines)

def tab2vcf(filename, sep='\t'):

    outfile=filename+'.vcf'
    fh_out = open(outfile, "w")
    fh_out.write(vcfheader(filename)+'\n')
    fh = open(filename)
    linenum = 0
    chromind = -1
    posind = -1
    refind = -1
    altind = -1
    rsind = -1

    for line in fh:
        line=line.strip()
        fields=line.split(sep)
        if linenum==0:
            chromind=find_first_index(fields, 'chrom')
            posind=find_first_index(fields, 'pos')
            refind=find_first_index(fields, 'ref')
            altind=find_first_index(fields, 'alt')
            rsind=find_first_index(fields, 'id')
            sampleind = find_first_index(fields, 'sample_id')
            geneind = find_first_index(fields, 'gene')
            qualind = find_first_index(fields, 'qual')
            filterind = find_first_index(fields, 'filter')
            gtind = find_first_index(fields, 'gt')

            print(str(chromind) + ' ' +  str(posind) + ' ' +  str(refind) + ' ' +  str(altind) + ' ' +  str(rsind) + ' ' +  str(sampleind) + ' ' +  str(geneind)  + ' ' +  str(qualind)  + ' ' +  str(filterind) + ' ' +  str(gtind))

            if(chromind < 0 or posind < 0 or refind < 0 or altind < 0):
                print("Column names chrom, pos, ref and alt are mandatory")
                break;
        else:

            chr=fields[chromind].strip().replace('chr', '')
            pos=str(fields[posind]).strip()
            ref=str(fields[refind]).strip()
            alt=str(fields[altind]).strip()
            id=str(fields[rsind]).strip()
            sample_id=str(fields[sampleind]).strip()
            gene= str(fields[geneind]).strip()
            qual = str(fields[qualind]).strip()
            filter = str(fields[filterind]).strip()
            gt = str(fields[gtind]).strip()
            
            if len(fields)>4:
                info=gene
            if (alt != ref) and (find_first_index(ACCEPTED_CHR, chr.strip()) > -1):
                l= (chr + sep + pos + sep + id +sep + ref + sep + alt + sep + qual + sep + filter + sep + "gene={},sample_id={}".format(info, sample_id) + sep + "GT" + sep + gt ).strip()
                #print l
            fh_out.write(l+'\n')

        linenum=linenum+1

    fh.close()
    fh_out.close()


def run(filename):
    if os.path.exists(filename) and os.path.isfile(filename):
        tab2vcf(filename)
    else:
        print("File does not exist " + filename)



filename = ''
while True:
    filename= args.input_tsv
    if os.path.exists(filename) and os.path.isfile(filename):
        break
    else:
        print 'File ' + filename + ' does not exists'

print 'You entered ' + filename

run(filename=filename)