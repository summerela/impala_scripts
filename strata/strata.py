#!/usr/bin/env python

from impala.dbapi import connect
from impala.util import as_pandas
import subprocess as sp

impala_host = 'glados14'
impala_port = 21050
impala_name = 'selasady'

conn = connect(host=impala_host, port=impala_port, timeout=10000, user=impala_name)
cur = conn.cursor()

# create reusable query function
def run_query(input_query):
    try:
        cur.execute(input_query)
    except Exception as e:
        print e

# add mac id to irg data
create_irg_mac = '''
create table strata.irg_mac
(
  pos int,
  subject_id string,
  zygosity int,
  mac int,
  adj_mac double
  )
partitioned by (chrom string, blk_pos int)
'''

run_query(create_irg_mac)

insert_irg_mac = '''
insert into table strata.irg_mac partition (chrom, blk_pos)
select i.strata_id as pos, i.subject_id, i.zygosity,
  s.mac, (((s.mac + RAND()*(1-0)+0)/5008) * 100) as adj_mac,
  i.chrom, cast(i.strata_id/1000000 as int) as blk_pos
  from strata.irg_chrom i, strata.sts_partitioned s
  where i.chrom = s.chrom
  and i.strata_id = s.pos
'''

run_query(insert_irg_mac)

# join sts and irg information, with 'all' for sts subject_id
union_sts_irg = '''

'''



