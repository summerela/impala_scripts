#################
# load modules  #
#################
import pandas as pd
from impala.dbapi import connect
from impala.udf import udf, ship_udf
from string import split

#create database connection
conn = connect(host='glados19', port=21050)
cur = conn.cursor()




@udf
def hour_from_weird_data_format(context, date):
    return split(date, '-')[1]

ship_udf(cur, hour_from_weird_data_format, '/users_selasady/test_udf.ll', 'glados19')

cur.execute('SELECT hour_from_weird_data_format(date) AS hour FROM public_hg19.clinvar LIMIT 100')
cur.fetchall()