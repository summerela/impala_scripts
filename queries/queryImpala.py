__author__ = 'summerrae'

from impala.dbapi import connect
conn = connect(host='my.impala.host', port=21050)
cur = conn.cursor()