import sys

from glu.lib.fileutils import table_reader, table_writer

seen = set()

out = table_writer(sys.stdout)

for row in table_reader(sys.stdin):
  name = row[0]
  if name in seen:
    continue
  seen.add(name)
  out.writerow(row)
