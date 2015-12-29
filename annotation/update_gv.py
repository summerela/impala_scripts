# setup variables
input_db = "p7_platform"
input_table =  "wgs_illumina_variant"

# download variants from input table that are not in global_vars
select chrom, pos, id, ref, alt
from p7_platform.wgs_illumina_variant v
where v.chrom = '1'
and not exists (
  (select chrom, pos, rs_id as id, ref, alt
   from global_variants gv
   where gv.chrom = '1'
   and v.chrom = gv.chrom
   and v.pos = gv.pos
   and v.id = gv.rs_id
   and v.ref  = gv.ref
   and v.alt = gv.alt)
  )

# if len(pandas table) > 0:
# annoate with sub-tables
# run through snpeff
# add to table