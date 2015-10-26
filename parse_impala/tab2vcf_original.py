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

#   Redistribution and use in source and binary forms, with or without modification,
#   are permitted provided that the following conditions are met:
#
#   Redistributions of source code must retain the above copyright notice,
#   this list of conditions and the following disclaimer.
#
#   Redistributions in binary form must reproduce the above copyright notice,
#   this list of conditions and the following disclaimer in the documentation and/or
#   other materials provided with the distribution.
#
#   THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
#   ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
#   WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED.
#   IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT,
#   INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING,
#   BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
#   DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY
#   OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
#   NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE,
#   EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
#
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
    lines.append('#CHROM'+'\t'+'POS'+'\t'+'ID'+'\t'+'REF'+'\t'+'ALT'+'\t'+'QUAL'+'\t'+'FILTER'+'\t'+'INFO'+'\t'+'FORMAT' + '\t'+ filename)
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

    for line in fh:
        line=line.strip()
        fields=line.split(sep)
        if linenum==0:
            chromind=find_first_index(fields, 'CHROM')
            posind=find_first_index(fields, 'POS')
            refind=find_first_index(fields, 'REF')
            altind=find_first_index(fields, 'ALT')

            print(str(chromind) + ' ' +  str(posind) + ' ' +  str(refind) + ' ' +  str(altind))

            if(chromind < 0 or posind < 0 or refind < 0 or altind < 0):
                print("Column names CHROM, POS, REF and ALT are mandatory")
                break;
        else:

            chr=fields[chromind].strip().replace('chr', '')
            pos=str(fields[posind]).strip()
            ref=str(fields[refind]).strip()
            alt=str(fields[altind]).strip()
            info='.'
            if len(fields)>4:
                info=';'.join( map (str, fields[4:len(fields)]) )
            if (alt != ref) and (find_first_index(ACCEPTED_CHR, chr.strip()) > -1):
                l= (chr + sep + pos + sep + '.' +sep + ref + sep + alt + sep + '.' + sep + 'PASS' + sep + info + sep + 'GT' + sep + './.' ).strip()
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
    filename=raw_input("Enter pileup file name: ")
    if os.path.exists(filename) and os.path.isfile(filename):
        break
    else:
        print 'File ' + filename + ' does not exists'

print 'You entered ' + filename

run(filename=filename)
