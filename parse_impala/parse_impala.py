__author__ = 'selasady'

from impala.dbapi import connect
import os
import subprocess as subp
import pandas as pd

def label_member(tbl_name, trio_arg):
    """
    function to create a sql statement from user trio argument for which
    trio members to include in analysis
    :param tbl_name: alias of tablename used in query as string, ex. 'bv'
    :param trio_arg: list of trio members to include, ex. 'M','F','NB'
    :return: member_arg
    """
    member_list = []
    for member in trio_arg:
        if member == 'NB':
            member_list.append("{}.sample_id LIKE '%03'".format(tbl_name))
        if member == 'M':
            member_list.append("{}.sample_id LIKE '%01'".format(tbl_name))
        if member == 'F':
            member_list.append("{}.sample_id LIKE '%02'".format(tbl_name))
        if member_list == 'all':
            member_list =''
    if len(member_list) > 0:
        member_arg = 'AND (' + ' OR '.join(member_list) + ')'
    # otherwise statment is empty
    else:
        member_arg = ''
    return member_arg


def select_samples(tbl_name, sample_list):
    """
    format sample id argument from user arg for which sample id's to include
    :param tbl_name: tbl_name: alias of tablename used in query as string, ex. 'bv'
    :param sample_list: list of sample id's to include, ex. '101-555-02','101-555-01'
    :return: subject_list
    """
    sample_arg = []
    if sample_list != 'all':
        sample_arg.append("AND {}.sample_id IN ".format(tbl_name) + str(sample_list))
        subject_list = ", ".join(str(i) for i in sample_arg)
    else:
        subject_list = ''
    return subject_list

def join_args(arg_list):
    """
    joins together sql statements, adding 'WITH' to
    first argument, and AND to each following arg
    :param arg_list: list of sql statements to join
    :return: subject_statement
    """
    if len(arg_list) > 0:
        subject_statement = ' '.join(arg_list)
    # otherwise return empty string
    else:
        subject_statement = ''
    return subject_statement

def run_query(query_name, remove_name, out_db, out_name):
    """
    opens odbc connection to impala, drops table if exists, runs query
     saves results as table, closes connection
    :param query_name: name of query statement to run
    :param remove_name: name of the table to drop if it already exists
    :param out_db: name of output database
    :param out_name: name of output table name
    :return: table of results saved on impala in specified output dir
    """
    # create connection object
    conn=connect(host='glados19', port=21050, timeout=120)
    # drop table if it exists
    cur = conn.cursor()
    print 'Removing table if it already exists...'
    cur.execute('DROP TABLE IF EXISTS {}.{}_{}'.format(out_db,out_name, remove_name))
    # run query
    print 'Running the following query on impala: \n' + query_name
    cur.execute(query_name)
    cur.execute('COMPUTE STATS {}.{}_{}'.format(out_db,out_name, remove_name))
    print 'Query finished. Closing connection.'
    cur.close()
    conn.close()

def df_to_snpeff(input_df, snpeff_path, out_name):
    """
    outputs variants to vcf format and run through snpeff
    :param input_df: pandas dataframe to annotate with snpeff
    :param snpeff_path: path to snpeff jar file
    :param out_name: table name to use for output
    :return: vcf file of snpeff coding consequence predictions
    to be matched back to original pandas dataframe
    """
    if os.path.exists(snpeff_path):
        # setup file names
        tsv_outname = '{}.tsv'.format(os.path.join(os.getcwd(), out_name))
        vcf_outname = '{}.vcf'.format(os.path.join(os.getcwd(), out_name))
        # these columns are output to vcf file
        df = input_df[['chrom', 'pos', 'ref', 'alt', 'qual', 'filter', 'gt']]
        # write to file for conversion to vcf
        df.to_csv(tsv_outname, header=True, encoding='utf-8', sep="\t", index=False)
        # run tab2vcf and upon success, run snpeff
        vcf_process = subp.Popen(['python', './nbs_genes/tab2vcf.py', tsv_outname])
        vcf_process.wait()
        # command to run snpeff
        snpeff_cmd = 'java -Xmx4g -jar {}  -t -v -noStats GRCh37.74 {} > {}_snpeff.vcf'.format(
            snpeff_path, vcf_outname, out_name)
        # run snpeff on vcf file
        snpeff_process = subp.Popen(snpeff_cmd, shell=True)
        snpeff_process.wait()
    else:
        print "Make sure you entered the correct path to snpEff.jar"

def drop_y(df):
    """
    drop extra columns ending in _y from df merge
    :param df: pandas dataframe to remove _y from column names
    :return:same pandas dataframe, with column names parsed
    """
    for col in df:
        if col.endswith('_y'):
            df.drop(col, axis=1, inplace=True)

def parse_snpeff(tx_df, notx_df, input_vcf):
    """
    match snpeff annotations back to data frame
    :param tx_df: dataframe of variants that contain transcript id's
    :param notx_df: dataframe of variants that have no transcript id,
    will be matched by genoimic position instead
    :param input_vcf:snpeff vcf file of annotations to match with dataframe
    :return: df_functional a pandas dataframe with snpeff annotations
    """
    # read in snpeff vcf file
    annot_vcf = pd.read_csv(input_vcf, sep='\t', skiprows=8)
    # split info field into separate rows for each transcript
    info_df = pd.Series([j for i in annot_vcf['INFO'].str.split(',') for j in i])
    # split each rown into separate columns by the pipe
    info_df = pd.DataFrame(list(info_df.str.split('|')))
    # drop emtpy/unnecessary columns
    info_df = info_df[list(info_df.columns[0:11])]
    info_df.columns = ['alt', 'effect', 'impact', 'gene_name', 'gene_id', 'feature_type', 'transcript_id',
                     'tx_biotype', 'rank', 'hgvs_c', 'hgvs_p']
    # remove the annotation header 'ANN=' from the alt field
    info_df['alt'] = info_df['alt'].str.replace('ANN=', '')
    # keep only transcript level feature types
    info_df = info_df[(info_df['feature_type'] == 'transcript')]
    # merge annotations with variant table
    tx_functional = pd.merge(tx_df, info_df, on=['transcript_id', 'alt'], how='left')
    notx_functional = pd.merge(notx_df, info_df, on=['chrom', 'pos', 'ref', 'alt'], how='left')
    # drop y columns
    drop_y(tx_functional)
    drop_y(notx_functional)
    # rename columns ending in _x from merge (TODO combine first then parse)
    tx_functional.rename(columns=lambda x: x.replace('_x', ''), inplace=True)
    tx_functional.drop_duplicates(inplace=True)
    notx_functional.rename(columns=lambda x: x.replace('_x', ''), inplace=True)
    notx_functional.drop_duplicates(inplace=True)
    # merge data frames back together
    functional_annot = pd.concat([tx_functional, notx_functional]).drop_duplicates().reset_index(drop=True)
    return functional_annot

# # function to match snpeff annotations back to data frame
# # when no tx id available by genomic pos
# def snpeff_notx(input_df, input_vcf):
#     # read in snpeff vcf file
#     annot_vcf = pd.read_csv(input_vcf, sep='\t', skiprows=8)
#     # split info field into separate rows for each transcript
#     info_df = pd.Series([j for i in annot_vcf['INFO'].str.split(',') for j in i])
#     # split each rown into separate columns by the pipe
#     info_df = pd.DataFrame(list(info_df.str.split('|')))
#     # drop emtpy/unnecessary columns
#     info_df = info_df[list(info_df.columns[0:11])]
#     info_df.columns = ['alt', 'effect', 'impact', 'gene_name', 'gene_id', 'feature_type', 'transcript_id',
#                      'tx_biotype', 'rank', 'hgvs_c', 'hgvs_p']
#     # remove the annotation header 'ANN=' from the alt field
#     info_df['alt'] = info_df['alt'].str.replace('ANN=', '')
#     # keep only transcript level feature types
#     info_df = info_df[(info_df['feature_type'] == 'transcript')]
#     #merge annotations with variant table
#     df_functional =
#     drop_y(df_functional)
#     # rename columns ending in _x from merge
#     df_functional.rename(columns=lambda x: x.replace('_x', ''), inplace=True)
#     df_functional.drop_duplicates(inplace=True)
#     return df_functional