from impala.dbapi import connect

#connect to impala
conn=connect(host='glados19', port=21050)
#create a cursor object to interact with db
cur = conn.cursor()

# #view databases
# cur.execute('SHOW DATABASES')
# #fetch results of cur.execute
# print cur.fetchall()
#
# #view tables
# cur.execute('SHOW TABLES')
# print cur.fetchall()
#
# #execute a query
# cur.execute('SELECT * FROM feature_fm LIMIT 5')
# #view schema of result set
# print cur.description
# #get results and view
# for row in cur:
#     print row

##grab a specific column
# cur.execute('SELECT feature_id FROM feature_fm LIMIT 5')
# for row in cur:
#     print row

#grab multiple columns
# cur.execute('SELECT feature_id, pos FROM feature_fm LIMIT 5')
# for row in cur:
#     print row

#grab unique rows
#cur.execute('SELECT DISTINCT feature_id FROM feature_fm LIMIT 5')

#show columns
# cur.execute('DESCRIBE feature_fm')
# for row in cur:
#     print row

#select chr 1
# cur.execute('SELECT * FROM feature_fm WHERE chr IN ("1")')
# for row in cur:
#     print row

#select chr 1,2 and 3
# cur.execute('SELECT * FROM feature_fm WHERE chr BETWEEN "1" AND "3"')
# for row in cur:
#     print row

#find all entries where feature ID contains meth and order by chr
# cur.execute('SELECT * FROM feature_fm WHERE feature_id LIKE "%METH%" ORDER BY chr DESC')
# for row in cur:
#     print row

##count entries for each chromosome
# cur.execute('SELECT COUNT(*) AS "Chromsome Count" FROM feature_fm GROUP BY chr')
# for row in cur:
#     print row


# #now lets store results as a pandas table
from impala.util import  as_pandas
import numpy as np
import matplotlib.pyplot as plt
cur.execute('SELECT * FROM feature_fm LIMIT 5')
df = as_pandas(cur)
type(df)
print df