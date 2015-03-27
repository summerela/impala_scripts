# -*- coding: utf-8 -*-

#**********************************************************************************
# 01/25/13 (Prachi): - iterate through all ALT alleles
#                    - display annotations per allele/transcript
#                    - change annotation format to alleleIndex:ref->var:gene:geneid:location:function:mRNA:protein:exon:strand:aaPos:impact
#                    - display transcript ID for intronic changes
#                    - evaluate PPC status based on all allele/transcript combinations
#                    - update vannotator.py for AA position derivation from chrpos
#                    - show Indel in FILTER if one of the variants is an indel
# 06/03/14 (Prachi): - add annotation for miRNA from mirBase
#***********************************************************************************

from __future__ import division

__gluindex__  = True
__abstract__  = 'Annotate genomic coordinates and nucleotide variation'
__copyright__ = 'Copyright (c) 2010, BioInformed LLC and the U.S. Department of Health & Human Services. Funded by NCI under Contract N01-CO-12400.'
__license__   = 'See GLU license for terms by running: glu license'
__revision__  = '$Id$'

import sys

from   operator                     import itemgetter

import pysam
import re

from   glu.lib.utils                import unique
from   glu.lib.fileutils            import list_reader, autofile, hyphen, table_reader, table_writer
from   glu.lib.progressbar          import progress_loop
from   glu.lib.recordtype           import recordtype

from   glu.lib.genedb.queries       import query_segdups, query_repeats

from   glu.lib.seqlib.vcf           import VCFReader, VCFWriter
from   glu.lib.seqlib.cga           import cga_reader
from   vannotator                   import VariantAnnotator
from   mirbase                      import MIRBase
from   cgfvariants                  import CGFVariants
from   glu.lib.seqlib.kaviar        import kaviar_reader
from   glu.lib.seqlib.refvariants   import ReferenceVariants


class OrderedReader(object):
  def __init__(self, data, references, get_contig=itemgetter(0), get_loc=itemgetter(1)):
    self.data       = data
    self.refmap     = dict( (r,i) for i,r in enumerate(references) )
    self.current    = next(self.data,None)
    self.get_contig = get_contig
    self.get_loc    = get_loc

  def get(self, chrom, loc):
    current = self.current
    if current is None:
      return None

    try:
      data = self.data
      get_contig = self.get_contig

      current_chrom = get_contig(current)

      if chrom!=current_chrom:
        refmap   = self.refmap
        chromidx = refmap[chrom]

        # Stream is past current contig, do nothing...
        if chromidx<refmap[current_chrom]:
          return None

        while chromidx>refmap[current_chrom]:
          last_chrom = current_chrom
          while last_chrom==current_chrom:
            current       = next(data)
            current_chrom = get_contig(current)

      assert chrom==current_chrom

      get_loc     = self.get_loc
      current_loc = get_loc(current)

      while current_loc<loc:
        current       = next(data)
        current_chrom = get_contig(current)
        current_loc   = get_loc(current)

        if current_chrom!=chrom or current_loc>loc:
          self.current = current
          return None

      if current_loc>loc:
        self.current = current
        return None

      assert current_loc==loc
      results = []

      while 1:
        results.append(current)

        self.current = current = next(data,None)

        if current is None:
          break

        current_chrom = get_contig(current)
        current_loc   = get_loc(current)

        if current_chrom!=chrom or current_loc!=loc:
          break

      return results

    except StopIteration:
      self.current = None
      return None


def make_infomap(info):
  infomap = {}
  for inf in info:
    if '=' in inf:
      key,value = inf.split('=',1)
    else:
      key,value = inf,''

    infomap[key] = value
  return infomap


def polyphen2_code(c):
  if c=='?':
    return '.'
  elif len(c)==1:
    return c
  elif 'D' in c:
    return'D'
  elif 'd' in c:
    return'd'
  elif '?' in c:
    return'.'
  elif 'b' in c:
    return'b'
  else:
    return'.'


def aa_change(e):
  return '%s->%s' % (e.ref_aa,e.var_aa) if e.ref_aa or e.var_aa else ''


def update_vcf_annotation(v, vs, cv, mirbase, esp, kaviar, refvars, polyphen2, options):
  new_info = []
  ppc_flag = 0
  indel_flag = 0

  if vs:
    # change annotation format and iterate over all ALT alleles - Prachi
    # FIXME: Order genes and evidence consistently
    anno_string = []
    cytoband = []
    gene_flag = 0
    genes = []

    for alleleIndex,allele in enumerate(v.var):
        allele = '' if allele=='.' else allele
	evidence   = list(vs.annotate(v.chrom, v.start, v.end, allele, nsonly=False))

	for e in evidence:
		ppc_flag_anno = ''
                gene = ''
		if e.func:
			ppc_flag = 1
			ppc_flag_anno = 'PPC'
                if e.gene and e.gene.symbol:
                        gene = e.gene.symbol
                        genes.append(gene)
                if e.gene:
			gene_flag = 1
                if not v.ref or allele=='' or len(v.ref)!=len(allele) or v.ref=='.' or allele=='.':
                        indel_flag = 1
		geneid	 = e.gene.geneid if e.gene and e.gene.geneid else ''
		location = e.intersect if e.intersect else ''
		function = e.func_type or e.func_class if e.func else ''
		mRNA = re.findall('(?=mRNA=(.*?):)', e.details)[0] if re.findall('(?=mRNA=(.*?):)', e.details) else ''
		protein = re.findall('(?=protein=(.*?):)', e.details)[0] if re.findall('(?=protein=(.*?):)', e.details) else ''
		exon = re.findall('(?=exon=(.*?):)', e.details)[0] if re.findall('(?=exon=(.*?):)', e.details) else ''
		strand = re.findall('(?=strand=([\+\-]))', e.details)[0] if re.findall('(?=strand=([\+\-]))', e.details) else ''
		aaPos = re.findall('(?=aa=(.*):)', e.details)[0] if re.findall('(?=aa=(.*):)', e.details) else ''
                cDNA = ('c.' + re.findall('(?=cdna=(.*))', e.details)[0]) if re.findall('(?=cdna=(.*))', e.details) else ''
                impact	 = aa_change(e) if aa_change(e) else ''
		anno_string.append(('%s|%s->%s|%s|%s|%s|%s|%s|%s|%s|%s|%s|%s|%s|%s' % (alleleIndex+1,v.ref,allele,ppc_flag_anno,gene,geneid,location,function,mRNA,protein,exon,strand,aaPos,impact,cDNA)))
		if e.cytoband:
			cytoband.append(e.cytoband)
 
    anno_string = list(unique(anno_string))
    cytoband = list(unique(cytoband))
    genes = list(unique(genes))
    
    segdups    = query_segdups(vs.con, v.chrom, v.start, v.end)
    repeats    = query_repeats(vs.con, v.chrom, v.start, v.end)

    while 'PASS' in v.filter:
      v.filter.remove('PASS')

    if not gene_flag:
      v.filter.append('Intergenic')

    if ppc_flag:
      v.filter.append('PPC')

    if cytoband:
      new_info.append('CYTOBAND=%s' % (','.join(cytoband)))

    if genes:
      new_info.append('GENE_NAMES=%s' % (','.join(genes)))

    if anno_string:
      new_info.append('ITMI_FI=%s' % ('&'.join(anno_string  )))

    if segdups:
      segdup_info = ('chr%s:%s-%s:%.2f' % (s.other_chrom,s.other_start,s.other_stop,s.matchFraction*100) for s in segdups)
      segdup_info = ','.join(segdup_info)
      v.filter.append('SegDup')
      new_info.append('SEGDUP_COUNT=%d' % len(segdups))
      new_info.append('SEGDUP_INFO=%s' % segdup_info)

    if repeats:
      repeat_info = ('%s:%s:%s' % (r.repeatName,r.repeatClass,r.repeatFamily) for r in repeats)
      repeat_info = ','.join(repeat_info)
      v.filter.append('Repeat')
      new_info.append('REPEAT_COUNT=%d' % len(repeats))
      new_info.append('REPEAT_INFO=%s' % repeat_info)

  if cv:
    cvinfo  = cv.score_and_classify(v.chrom,v.start,v.end,[v.ref,v.var[0]])
    if cvinfo.exact_vars:
      v.names = sorted(set(v.replace('dbsnp:','rs') for v in cvinfo.exact_vars)|set(v.names))
      #Prachi
      if 'esp' in v.names:
        v.names.remove('esp')

    if cvinfo.exact_vars or cvinfo.common_score:
      #new_info.append('COMMON_SCORE=%.2f' % cvinfo.common_score)
      common_score_list = []
      for cs_allele,cs_score in cvinfo.common_score.items():
          common_score_list.append('%s:%.2f' % (cs_allele,cs_score))
      new_info.append('COMMON_SCORE=%s' % (','.join(common_score_list)))

    if cvinfo.function_info:
      function_info = ','.join(cvinfo.function_info)
      new_info.append('FUNCTION_INFO=%s' % function_info)

    if cvinfo.inexact_vars:
      inexact = ','.join(sorted(set(cvinfo.inexact_vars)))
      new_info.append('INEXACT_VARIANTS=%s' % inexact)

    #if cvinfo.common_score>options.commonscore:
    #  v.filter.append('Common')

  if mirbase:
    mirbaseinfo = mirbase.find_mirna(v.chrom,v.start,v.end)
    if mirbaseinfo.names:
     mirnames = ','.join(mirbaseinfo.names)
     new_info.append('MIRNA=%s' % mirnames)

  if esp:
    chrom = v.chrom
    vars  = set(v.var)
    if chrom.startswith('chr'):
      chrom = chrom[3:]

    try:
      esp_info = esp.fetch(chrom,v.start,v.end)
      esp_info = [ e for e in esp_info if e.start==v.start and e.end==v.end and vars.intersection(e.var) ]
    except ValueError:
      esp_info = None
    except KeyError:
      esp_info = None

    if esp_info:
      gtc = maf_ea = maf_aa = 0

      for einfo in esp_info:
        infomap  = make_infomap(einfo.info)
        if 'MAF' in infomap:
          mafs     = map(float,infomap['MAF'].split(','))
          maf_ea   = max(maf_ea,mafs[0]/100)
          maf_aa   = max(maf_aa,mafs[1]/100)

        if 'GTC' in infomap:
          gtcs     = map(int,infomap['GTC'].split(','))
          gtc      = max(gtc,gtcs[0]+gtcs[1])

      maf = max(maf_ea,maf_aa)

      v.filter.append('ESP')
      if maf>options.commonscore:
        v.filter.append('ESPCommon')

      new_info.append('ESP_COUNT=%d' % gtc)
      new_info.append('ESP_MAF_EA=%.4f' % maf_ea)
      new_info.append('ESP_MAF_AA=%.4f' % maf_aa)

  if kaviar:
    kinfo = kaviar.get(v.chrom,v.start) or []
    ktext = [ stuff.replace(';',',').replace(' ','_') for chrom,loc,mallele,maf,allele,stuff in kinfo if allele in v.var ]

    if ktext:
      v.filter.append('Kaviar')

      mallele = kinfo[0][2]
      maf     = kinfo[0][3]

      if mallele and mallele in v.var:
        new_info.append('KAVIAR_MAF=%.2f' % maf)

        if maf>options.commonscore:
          v.filter.append('KaviarCommon')

      new_info.append('KAVIAR_NAMES=%s' % ','.join(ktext))


  if refvars:
    ingroup,outgroup = refvars.get(v.chrom,v.start,v.end,v.var) if refvars else ([],[])

    new_info.append('REFVAR_INGROUP_COUNT=%d'  % len(ingroup))
    new_info.append('REFVAR_OUTGROUP_COUNT=%d' % len(outgroup))

    if ingroup:
      ingroup  = ','.join(ingroup)
      new_info.append('REFVAR_INGROUP_NAMES=%s'  % ingroup)

    if outgroup:
      outgroup = ','.join(outgroup)
      new_info.append('REFVAR_OUTGROUP_NAMES=%s' % outgroup)
      v.filter.append('RefVar')

  if vs and polyphen2 and v.end-v.start==1 and 'CDS' in location:
    pmap  = {'A':0,'C':1,'G':2,'T':3}

    try:
      pvars = [ p.rstrip().split('\t') for p in polyphen2.fetch(v.chrom,v.start,v.end) ]
    except ValueError:
      pvars = None
    except KeyError:
      pvars = None

    if pvars:
      hdivs = []
      hvars = []

      for a in v.var:
        if a in pmap:
          i = pmap[a]

          hdiv = hvar = ''
          for p in pvars:
            hdiv += p[3][i]
            hvar += p[4][i]

          hdiv = polyphen2_code(hdiv)
          hvar = polyphen2_code(hvar)
        else:
          hdiv = '.'
          hvar = '.'

        assert hdiv!='r' and hvar!='r'

        hdivs.append(hdiv)
        hvars.append(hvar)

      if hdivs.count('.')!=len(hdivs):
        new_info.append('POLYPHEN2_HDIV=%s' % (','.join(hdivs)))
      if hvars.count('.')!=len(hvars):
        new_info.append('POLYPHEN2_HVAR=%s' % (','.join(hvars)))

  if 'tgp' in v.names:
    v.names.remove('tgp')
    v.filter.append('1000G')

  if indel_flag:
    v.filter.append('Indel')

  # Remove any old fields that have been replaced by a new field
  # Remove old ITMI annotation fields if VCF is being reannotated - Prachi
  new_info_fields = set(f.split('=',1)[0] for f in new_info)
  old_fields_to_remove = ['GENE_NAME', 'GENE_ID', 'GENE_LOCATION', 'GENE_FUNCTION', 'GENE_FUNCTION_DETAILS']
  old_info = v.info
  v.info          = new_info+[ f for f in v.info if f.split('=',1)[0] not in new_info_fields ]
  v.info          = [ f for f in v.info if f.split('=',1)[0] not in old_fields_to_remove ]
 
  v.filter = sorted(set(v.filter)) or ['PASS']

  return v


def annotate_vcf(options):
  vs       = VariantAnnotator(options.genedb, options.reference)
  vcf      = VCFReader(options.variants,sys.stdin)
  cv       = CGFVariants(options.cgfvariants, options.reference) if options.cgfvariants else None
  mirbase  = MIRBase(options.mirbase) if options.mirbase else None
  esp      = VCFReader(options.esp) if options.esp else None
  polyphen2= pysam.Tabixfile(options.polyphen2) if options.polyphen2 else None

  if options.kaviar:
    references = list_reader(options.reference+'.fai')
    kaviar     = kaviar_reader(options.kaviar)
    kaviar     = OrderedReader(kaviar, references)
  else:
    kaviar     = None

  refvars  = ReferenceVariants(options.refvars,options.refingroup) if options.refvars else None

  metadata = vcf.metadata

  metadata.setdefault('FILTER',[])
  metadata.setdefault('INFO',[])

  #Prachi - prefix contig names with "chr"
  metadata['contig'] = [ctg.replace( 'contig=<ID=','contig=<ID=chr' ) for ctg in metadata['contig']]

  metadata['FILTER'].append('##FILTER=<ID=Indel,Description="Variant is an insertion or deletion">')
  metadata['FILTER'].append('##FILTER=<ID=Intergenic,Description="Variant not in or near a gene">')
  #metadata['FILTER'].append('##FILTER=<ID=NPF,Description="Variant is not predicted to alter a protein">')
  metadata['FILTER'].append('##FILTER=<ID=PPC,Description="Variant is predicted to alter a protein">')
  metadata['FILTER'].append('##FILTER=<ID=SegDup,Description="Variant occurs in a segmentally duplicated region">')
  metadata['FILTER'].append('##FILTER=<ID=Repeat,Description="Variant occurs in a repetitive or low-complexity region">')

  metadata['INFO'].append('##INFO=<ID=CYTOBAND,Number=.,Type=String,Description="Name of cytoband(s) containing variant">')
  metadata['INFO'].append('##INFO=<ID=GENE_NAMES,Number=.,Type=String,Description="Name of gene(s) containing variant">')
  #metadata['INFO'].append('##INFO=<ID=GENE_ID,Number=.,Type=String,Description="Entrez/LocusLink gene identifiers of genes containing variant">')
  #metadata['INFO'].append('##INFO=<ID=GENE_LOCATION,Number=.,Type=String,Description="Location of variant in gene(s)">')
  #metadata['INFO'].append('##INFO=<ID=GENE_FUNCTION,Number=.,Type=String,Description="Functional classification of variant for each gene and transcript">')
  metadata['INFO'].append('##INFO=<ID=ITMI_FI,Number=.,Type=String,Description="Annotation of variant for each allele/gene/transcript in the order AlleleIndex:Ref->Alt:PPCStatus:GeneName:GeneID:Location:Function:TranscriptID:ProteinID:ExonNum:Strand:AAPosition:AAImpact">')
  #metadata['INFO'].append('##INFO=<ID=GENE_FUNCTION_DETAILS,Number=.,Type=String,Description="Functional details of variant for each gene and transcript">')
  metadata['INFO'].append('##INFO=<ID=SEGDUP_COUNT,Number=1,Type=Integer,Description="Number of segmental duplications that overlap variant locus">')
  metadata['INFO'].append('##INFO=<ID=SEGDUP_INFO,Number=.,Type=String,Description="Details of segmental duplications that overlap variant locus">')
  metadata['INFO'].append('##INFO=<ID=REPEAT_COUNT,Number=1,Type=Integer,Description="Number of repetitive or low complexity elements that overlap variant locus">')
  metadata['INFO'].append('##INFO=<ID=REPEAT_INFO,Number=.,Type=String,Description="Details of repetitive or low complexity elements that overlap variant locus">')

  if cv:
    metadata['FILTER'].append('##FILTER=<ID=Common,Description="Variant is likely common with common score>%f">' % options.commonscore)
    metadata['FILTER'].append('##FILTER=<ID=1000G,Description="Variant was reported by 1000 Genomes project">')
    metadata['INFO'  ].append('##INFO=<ID=COMMON_SCORE,Number=1,Type=Float,Description="Common score: maximum allele frequency in any population for rarest allele">')
    metadata['INFO'  ].append('##INFO=<ID=FUNCTION_INFO,Number=.,Type=String,Description="Annotated as function by OMIM, dbSNP, or COSMIC">')
    metadata['INFO'  ].append('##INFO=<ID=INEXACT_VARIANTS,Number=.,Type=String,Description="Observed variant matches inexactly: it has different alleles or overlaps observed">')
  if mirbase:
    metadata['INFO'  ].append('##INFO=<ID=MIRNA,Number=.,Type=String,Description="mirBase IDs overlapping with the variant">')
  if esp:
    metadata['FILTER'].append('##FILTER=<ID=ESP,Description="Variant appears in the UW ESP database">')
    metadata['FILTER'].append('##FILTER=<ID=ESPCommon,Description="Variant appears in the UW ESP database and appears to be common with MAF>%f">' % options.commonscore)
    metadata['INFO'  ].append('##INFO=<ID=ESP_COUNT,Number=1,Type=Integer,Description="Count of individuals with one or more variant alleles in UW ESP database">')
    metadata['INFO'  ].append('##INFO=<ID=ESP_MAF_EA,Number=1,Type=Float,Description="Minor allele frequency in European-Americans according to UW ESP database">')
    metadata['INFO'  ].append('##INFO=<ID=ESP_MAF_AA,Number=1,Type=Float,Description="Minor allele frequency in African-Americans according to UW ESP database">')

  if kaviar:
    metadata['FILTER'].append('##FILTER=<ID=Kaviar,Description="Variant appears in the Kaviar database">')
    metadata['FILTER'].append('##FILTER=<ID=KaviarCommon,Description="Variant appears in the Kaviar database and appears to be common with MAF>%f">' % options.commonscore)
    metadata['INFO'  ].append('##INFO=<ID=KAVIAR_MAF,Number=1,Type=Float,Description="Minor allele frequency according to Kaviar database">')
    metadata['INFO'  ].append('##INFO=<ID=KAVIAR_NAMES,Number=.,Type=String,Description="Samples or datasets from Kaviar in which variant was found">')

  if refvars:
    metadata['FILTER'].append('##FILTER=<ID=RefVar,Description="Variant appears in the local reference variant list">')
    metadata['INFO'  ].append('##INFO=<ID=REFVAR_INGROUP_COUNT,Number=1,Type=Integer,Description="Count of times variant is present in intra-group samples">')
    metadata['INFO'  ].append('##INFO=<ID=REFVAR_OUTGROUP_COUNT,Number=1,Type=Integer,Description="Count of times variant is present in extra-group samples">')
    metadata['INFO'  ].append('##INFO=<ID=REFVAR_INGROUP_NAMES,Number=.,Type=String,Description="Intra-group samples in which Variant is present">')
    metadata['INFO'  ].append('##INFO=<ID=REFVAR_OUTGROUP_NAMES,Number=.,Type=String,Description="Extra-group samples in which Variant is present">')

  if polyphen2:
    metadata['INFO'  ].append('##INFO=<ID=POLYPHEN2_HDIV,Number=.,Type=String,Description="Polyphen2 HDIV prediction code for each variant SNV allele (b=benign, d=possibly damaging, D=probably damaging, .=unknown)">')
    metadata['INFO'  ].append('##INFO=<ID=POLYPHEN2_HVAR,Number=.,Type=String,Description="Polyphen2 HVAR prediction code for each variant SNV allele (b=benign, d=possibly damaging, D=probably damaging, .=unknown)">')

  # Prachi - for writing to anno file
  out = VCFWriter(options.output, metadata, vcf.samples, options.reference)

  for v in vcf:
    update_vcf_annotation(v, vs, cv, mirbase, esp, kaviar, refvars, polyphen2, options)

    out.write_locus(v)
    


def update_ion_annotation(v, vs, cv, esp, kaviar, refvars, polyphen2, options):
  new_info = []

  if vs:
    # FIXME: Order genes and evidence consistently
    evidence   = list(vs.annotate(v.chrom, v.start, v.end, v.var, nsonly=False))
    #v.names   = sorted(set(str(v) for e in evidence for v in e.varid_exact)|set(v.names))
    cytoband   = sorted(set(e.cytoband    for e in evidence if e.cytoband))
    genes      = sorted(set(e.gene.symbol for e in evidence if e.gene and e.gene.symbol))
    geneids    = sorted(set(e.gene.geneid for e in evidence if e.gene and e.gene.geneid))
    location   = sorted(set(e.intersect   for e in evidence if e.intersect))
    function   = sorted(set(e.func_type or e.func_class  for e in evidence if e.func))
    nsevidence = [ e for e in evidence if e.func ]
    nsinfo     = ( '%s:%s:%s:%s' % (e.gene.symbol,e.func_type,e.details,aa_change(e))
                   for e in nsevidence )
    nsinfo     = list(unique(nsinfo))

    segdups    = query_segdups(vs.con, v.chrom, v.start, v.end)
    repeats    = query_repeats(vs.con, v.chrom, v.start, v.end)

    if not genes:
      v.Intergenic = 'Y'

    if not nsevidence:
      v.NPF = 'Y'

    if cytoband:
      v.CYTOBAND = ','.join(cytoband)

    if genes:
      v.GENE_NAME = ','.join(genes)

    if geneids:
      v.GENE_ID = ','.join(str(g) for g in geneids)

    if location:
      v.GENE_LOCATION = ','.join(location)

    if function:
      v.GENE_FUNCTION = ','.join(function)

    if nsinfo:
      v.GENE_FUNCTION_DETAILS = ','.join(nsinfo)

    if segdups:
      segdup_info = ('chr%s:%s-%s:%.2f' % (s.other_chrom,s.other_start,s.other_stop,s.matchFraction*100) for s in segdups)

      v.SegDup = 'Y'
      v.SEGDUP_COUNT = len(segdups)
      v.SEGDUP_INFO  = ','.join(segdup_info)

    if repeats:
      repeat_info = ('%s:%s:%s' % (r.repeatName,r.repeatClass,r.repeatFamily) for r in repeats)

      v.Repeat = 'Y'
      v.REPEAT_COUNT = len(repeats)
      v.REPEAT_INFO  = ','.join(repeat_info)

  if cv:
    cvinfo  = cv.score_and_classify(v.chrom,v.start,v.end,[v.ref,v.var])

    if cvinfo.exact_vars:
      if 'tgp' in cvinfo.exact_vars:
        v.TGP = 'Y'

      vars = set(v.replace('dbsnp:','rs') for v in cvinfo.exact_vars)
      vars.discard('tgp')

      v.VAR_NAME = ','.join(sorted(vars))

    if cvinfo.exact_vars or cvinfo.common_score>0:
      v.COMMON_SCORE = '%.2f' % cvinfo.common_score

    if cvinfo.function_info:
      v.FUNCTION_INFO = ','.join(cvinfo.function_info)

    if cvinfo.inexact_vars:
      v.INEXACT_VARIANTS = ','.join(sorted(set(cvinfo.inexact_vars)))

    if cvinfo.common_score>options.commonscore:
      v.Common = 'Y'

  if esp:
    chrom = v.chrom
    vars  = set(v.var)
    if chrom.startswith('chr'):
      chrom = chrom[3:]

    try:
      esp_info = esp.fetch(chrom,v.start,v.end)
      esp_info = [ e for e in esp_info if e.start==v.start and e.end==v.end and vars.intersection(e.var) ]
    except ValueError:
      esp_info = None

    if esp_info:
      gtc = maf_ea = maf_aa = 0

      for einfo in esp_info:
        infomap  = make_infomap(einfo.info)
        if 'MAF' in infomap:
          mafs     = map(float,infomap['MAF'].split(','))
          maf_ea   = max(maf_ea,mafs[0]/100)
          maf_aa   = max(maf_aa,mafs[1]/100)

        if 'GTC' in infomap:
          gtcs     = map(int,infomap['GTC'].split(','))
          gtc      = max(gtc,gtcs[0]+gtcs[1])

      maf = max(maf_ea,maf_aa)

      v.ESP = 'Y'
      if maf>options.commonscore:
        v.ESPCommon = 'Y'

      v.ESP_COUNT  = gtc
      v.ESP_MAF_EA = '%.4f' % maf_ea
      v.ESP_MAF_AA = '%.4f' % maf_aa

  if kaviar:
    kinfo = kaviar.get(v.chrom,v.start) or []
    ktext = [ stuff.replace(';',',').replace(' ','_') for chrom,loc,mallele,maf,allele,stuff in kinfo if allele in v.var ]

    if ktext:
      v.Kaviar = 'Y'

      mallele = kinfo[0][2]
      maf     = kinfo[0][3]

      if mallele and mallele in v.var:
        v.KAVIAR_MAF = '%.2f' % maf

        if maf>options.commonscore:
          v.KaviarCommon = 'Y'

      v.KAVIAR_NAMES = ','.join(ktext)


  if refvars:
    ingroup,outgroup = refvars.get(v.chrom,v.start,v.end,v.var) if refvars else ([],[])

    v.REFVAR_INGROUP_COUNT  = len(ingroup)
    v.REFVAR_OUTGROUP_COUNT = len(outgroup)

    if ingroup:
      v.REFVAR_INGROUP_NAMES = ','.join(ingroup)

    if outgroup:
      v.RefVar = 'Y'
      v.REFVAR_OUTGROUP_NAMES = ','.join(outgroup)

  if vs and polyphen2 and v.end-v.start==1 and 'CDS' in location:
    pmap  = {'A':0,'C':1,'G':2,'T':3}

    try:
      pvars = [ p.rstrip().split('\t') for p in polyphen2.fetch(v.chrom,v.start,v.end) ]
    except ValueError:
      pvars = None

    if pvars:
      hdivs = []
      hvars = []

      for a in [v.var]:
        if a in pmap:
          i = pmap[a]

          hdiv = hvar = ''
          for p in pvars:
            hdiv += p[3][i]
            hvar += p[4][i]

          hdiv = polyphen2_code(hdiv)
          hvar = polyphen2_code(hvar)
        else:
          hdiv = '.'
          hvar = '.'

        #assert hdiv!='r' and hvar!='r'

        hdivs.append(hdiv)
        hvars.append(hvar)

      if hdivs.count('.')!=len(hdivs):
        v.POLYPHEN2_HDIV = ','.join(hdivs)
      if hvars.count('.')!=len(hvars):
        v.POLYPHEN2_HVAR = ','.join(hvars)


  if not v.ref or not v.var:
    v.Indel = 'Y'

  return v


def annotate_ion(options):
  vars     = table_reader(options.variants,sys.stdin)
  header   = next(vars)

  fields  = ['Chrom','Position','Ref','Variant']

  try:
    indices = [ header.index(f) for f in fields ]
  except ValueError:
    missing = ', '.join(f for f in fields if f not in header)
    raise ValueError('Invalid IonTorrent Variant file. Missing headers: %s' % missing)

  var_fields = itemgetter(*indices)

  vs       = VariantAnnotator(options.genedb, options.reference)
  cv       = CGFVariants(options.cgfvariants, options.reference) if options.cgfvariants else None
  esp      = VCFReader(options.esp) if options.esp else None
  polyphen2= pysam.Tabixfile(options.polyphen2) if options.polyphen2 else None

  if options.kaviar:
    references = list_reader(options.reference+'.fai')
    kaviar     = kaviar_reader(options.kaviar)
    kaviar     = OrderedReader(kaviar, references)
  else:
    kaviar     = None

  refvars  = ReferenceVariants(options.refvars,options.refingroup) if options.refvars else None

  filters = ['Indel','Intergenic','NPF','SegDup','Repeat']
  info    = ['VAR_NAME','CYTOBAND','GENE_NAME','GENE_ID','GENE_LOCATION','GENE_FUNCTION',
             'GENE_FUNCTION_DETAILS','SEGDUP_COUNT','SEGDUP_INFO','REPEAT_COUNT','REPEAT_INFO']

  if cv:
    filters += ['Common','TGP']
    info    += ['COMMON_SCORE','FUNCTION_INFO','INEXACT_VARIANTS']


  if esp:
    filters += ['ESP', 'ESPCommon']
    info    += ['ESP_COUNT','ESP_MAF_EA','ESP_MAF_AA']

  if kaviar:
    filters += ['Kaviar','KaviarCommon']
    info    += ['KAVIAR_MAF','KAVIAR_NAMES']

  if refvars:
    filters += ['RefVar']
    info    += ['REFVAR_INGROUP_COUNT','REFVAR_OUTGROUP_COUNT','REFVAR_INGROUP_NAMES','REFVAR_OUTGROUP_NAMES']

  if polyphen2:
    info    += ['POLYPHEN2_HDIV','POLYPHEN2_HVAR']

  new_fields  = set(filters+info)
  keep_fields = itemgetter( *(i for i,h in enumerate(header) if h and h not in new_fields ) )

  new_fields  = filters+info
  blanks      = ['']*len(new_fields)
  keep_header = list(keep_fields(header))
  keep_len    = len(keep_header)
  new_header  = keep_header+filters+info
  Var         = recordtype('Var', ['chrom','start','end','ref','var']+filters+info)
  out         = table_writer(options.output,hyphen=sys.stdout)

  out.writerow(new_header)


  for rec in vars:
    chrom,pos,ref,var = var_fields(rec)

    pos   = int(pos)
    start = pos-1
    end   = pos

    # Fix indels encoded in the absurd VCF-like manner
    if ref and var and ref[0]==var[0]:
      ref    = ref[1:]
      var    = var[1:]
      start += 1
      end   += 1

      rec[indices[1]] = pos+1
      rec[indices[2]] = ref
      rec[indices[3]] = var

    v = Var( *([chrom,start,end,ref,var]+blanks) )
    update_ion_annotation(v, vs, cv, esp, kaviar, refvars, polyphen2, options)

    new_rec = list(keep_fields(rec))
    if len(new_rec)<keep_len:
      new_rec += ['']*(keep_len-len(new_rec))

    new_rec += list(v[5:])

    out.writerow(new_rec)


def valid_allele(a):
  return a is not None and a!='=' and 'N' not in a and '?' not in a


def annotate_mastervar(options):
  vs       = VariantAnnotator(options.genedb, options.reference)
  cv       = CGFVariants(options.cgfvariants, options.reference) if options.cgfvariants else None
  out      = autofile(hyphen(options.output,sys.stdout))

  if options.kaviar:
    references = list_reader(options.reference+'.fai')
    kaviar     = kaviar_reader(options.kaviar)
    kaviar     = OrderedReader(kaviar, references)
  else:
    kaviar     = None

  refvars    = ReferenceVariants(options.refvars,options.refingroup) if options.refvars else None

  attrs,header,vars = cga_reader(options.variants,sys.stdin)

  for var in vars:
    varid_exact     = []
    varid_inexact   = []
    common_score    = 0
    function_score  = 0
    kaviar_maf      = 0
    kaviar_details  = []
    geneinfo        = []

    allele1         = var.allele1Seq if valid_allele(var.allele1Seq) else None
    allele2         = var.allele2Seq if valid_allele(var.allele2Seq) else None

    alleles = []
    if allele1 is not None:
      alleles.append(allele1)
    if allele2 is not None:
      alleles.append(allele2)

    if len(alleles)!=2:
      continue

    if vs:
      evidence1     = list(vs.classify(var.chromosome, var.begin, var.end, allele1)) if allele1 is not None else []
      evidence2     = list(vs.classify(var.chromosome, var.begin, var.end, allele2)) if allele2 is not None else []
      evidence      = evidence1+evidence2

      geneinfo      = [ (e.gene.symbol, e.gene.geneid,
                         e.intersect,e.func_class,e.func_type,e.details,e.ref_aa,e.var_aa) for e in evidence if e.gene ]

      geneinfo      = list(unique(geneinfo))

    if cv:
      cvinfo        = cv.score_and_classify(var.chromosome,var.begin,var.end,alleles)
      if cvinfo.exact_vars:
        varid_exact = sorted(set(v.replace('dbsnp:','rs') for v in cvinfo.exact_vars)|set(varid_exact))
      if cvinfo.inexact_vars:
        varid_inexact = sorted(set(v.replace('dbsnp:','rs') for v in cvinfo.inexact_vars)|set(varid_inexact))

      #cvinfo.common_score)
      #cvinfo.function_score)

    if kaviar and var.end-var.begin==1:
      kinfo = kaviar.get(var.chromosome,var.begin) or []
      kaviar_details += [ stuff.replace(';',',') for chrom,loc,mallele,maf,allele,stuff in kinfo if allele in (allele1,allele2) ]

      if ktext:
        mallele = kinfo[0][2]
        maf     = kinfo[0][3]

        if mallele and mallele in (allele1,allele2):
          kaviar_maf = maf

    #ingroup,outgroup = refvars.get(var.chrom,var.start,var.end,var.var) if refvars else ([],[])

    #ingroup  = ','.join(ingroup)  if  ingroup else ''
    #outgroup = ','.join(outgroup) if outgroup else ''

    print var
    print '    varid_exact=',varid_exact
    print '  varid_inexact=',varid_inexact
    if geneinfo:
      print '      gene info=',geneinfo[0]
      for g in geneinfo[1:]:
        print '                ',g
    print


def option_parser():
  from glu.lib.glu_argparse import GLUArgumentParser

  parser = GLUArgumentParser(description=__abstract__)

  parser.add_argument('variants', help='Input variant file')

  parser.add_argument('-f', '--format',   metavar='NAME', default='VCF',
                      help='File format (VCF, ION)')
  parser.add_argument('-g', '--genedb',   metavar='NAME',
                      help='Genedb genome annotation database name or file')
  parser.add_argument('-r', '--reference',   metavar='NAME', required=True,
                      help='Reference genome sequence (FASTA + FAI files)')
  parser.add_argument('--cgfvariants',   metavar='NAME',
                      help='CGFvariant database annotation')
  parser.add_argument('--mirbase',   metavar='NAME',
                      help='mirBase database annotation')
  parser.add_argument('--commonscore', metavar='T', type=float, default=0.05,
                      help='Annotate all variants with common score > T')
  parser.add_argument('--esp',   metavar='NAME',
                        help='UW Exome Sequencing Project (ESP) annotation (optional)')
  parser.add_argument('--kaviar',   metavar='NAME',
                        help='Kaviar annotation (optional)')
  parser.add_argument('--refvars',   metavar='NAME',
                        help='Reference variant list')
  parser.add_argument('--polyphen2',   metavar='NAME',
                        help='Polyphen2 exome annotation (optional)')
  parser.add_argument('--refingroup',   metavar='NAME',
                        help='List of subjects defined to be intra-group reference variants')
  parser.add_argument('-o', '--output', metavar='FILE', default='-',
                    help='Output variant file')
  return parser


def main():
  parser  = option_parser()
  options = parser.parse_args()
  format  = options.format.upper()

  if format=='VCF':  # *** VCF ***
    annotate_vcf(options)
  elif format=='ION':
    annotate_ion(options)
  elif format=='MASTERVAR':
    annotate_mastervar(options)
  else:
    raise ValueError('Unknown or Unsupported format specified: %s' % options.format)


if __name__=='__main__':
  main()

