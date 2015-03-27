# -*- coding: utf-8 -*-

from __future__ import division


from   collections          import defaultdict, namedtuple

import pysam
import sys

from   glu.lib.recordtype   import recordtype


VarInfo    = namedtuple('VarInfo',   'names')
VarRecord  = recordtype('VarRecord', 'chromosome start stop names')
VariantKey = namedtuple('VariantKey', 'chromosome start stop')

goFlyAKite = VarInfo( ())
emptyMap = VarRecord('','','','')


class MIRBase(object):
  def __init__(self, mirbase):
    self.vars      = pysam.Tabixfile(mirbase)

  def query_variants(self, chromosome, start, stop):

    chrmap = {'X':23,'Y':24,'MT':25,'M':25}
    try:
     score  = self.vars.fetch(chrmap.get(chromosome,chromosome), start, stop)
     for s in score:
      chrom,vstart,vstop,names = s.split('\t')
      yield VarRecord(chrom, int(vstart), int(vstop), names)
    except KeyError as e:
     yield emptyMap

  def build_variants_lookup(self, chromosome, start, stop):
    vdict = defaultdict(list)
    vdata = self.query_variants(chromosome, start, stop)

    for d in vdata:
        vdict[VariantKey(d.chromosome,d.start,d.stop)].append(d)
    return vdict

  def find_mirna(self, chromosome, start, stop):
    vdata          = self.build_variants_lookup(chromosome, start, stop)
    names = []
    if not vdata:
      return goFlyAKite
    for v,vrecs in vdata.iteritems():
      for vrec in vrecs:
       names.append(vrec.names)
    names   = sorted(set(names))
    return VarInfo(names)
