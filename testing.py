import ibis
import os

ic = ibis.impala.connect(host='glados19', port=21050)

con = ibis.make_client(ic)


table = con.table('clinvar', database='public_hg19')

print table.limit(10)