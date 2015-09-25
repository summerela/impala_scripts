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

def parse_snpeff(input_df, input_vcf):
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
    # split info field into separate rows for each consequence
    temp = pd.concat([pd.Series(row['ID'], row['INFO'].split(','))
                    for _, row in annot_vcf.iterrows()]).reset_index()
    # split each row into separate columns by the pipe
    info_df = pd.DataFrame(list(temp['index'].str.split('|')))
    # add variant id to beginning of data frame
    info_df.insert(0, 'var_id', temp[0])
    # # drop emtpy/unnecessary columns
    info_df = info_df[list(info_df.columns[0:12])]
    info_df.columns = ['var_id', 'alt', 'effect', 'impact', 'gene_name', 'gene_id', 'feature_type', 'transcript_id',
                       'tx_biotype', 'rank', 'hgvs_c', 'hgvs_p']
    # remove the annotation header 'ANN=' from the alt field
    info_df['alt'] = info_df['alt'].str.replace('ANN=', '')
    # keep only transcript level feature types
    info_df = info_df[(info_df['feature_type'] == 'transcript')]
    # drop extraneous columns
    info_df.drop(info_df.columns[[1,4,6]], axis=1, inplace=True)
    # recombine info_df with annot_df
    snp_df = pd.merge(annot_vcf, info_df, left_on=['ID'], right_on=['var_id'])
    # drop extraneous snp_df columns
    snp_df.drop(snp_df.columns[[2, 7]], axis=1, inplace=True)
    # edit col names to lower case for matching
    snp_df.columns = map(str.lower, snp_df.columns)
    snp_df.columns = snp_df.columns.str.replace('#','')
    # drop duplicates from each df before merge
    snp_df.drop_duplicates(inplace=True)
    info_df.drop_duplicates(inplace=True)
    # merge annotations with variant table
    annot_vars = pd.merge(input_df, snp_df, on=['var_id'], how='left')
    # # remove all duplicated columns from merge ending in _y
    cols = [c for c in annot_vars.columns if not c.endswith('_y')]
    annotated_df =annot_vars[cols]
    # rename columns ending in _x from merge
    annotated_df.rename(columns=lambda x: x.replace('_x', ''), inplace=True)
    # drop duplicate rows
    annotated_df.drop_duplicates(inplace=True)
    # add family id
    annotated_df['family_id'] = annotated_df['sample_id'].apply(lambda x: x.split('-')[1])
    return annotated_df

def label_predictive(input_df):
    """
    create a column to label variants as predictive if
    they are marked as pathogenic in clinvar or
    are rare in kaviar and are high impact or predicted to
    be determintal by dbNSFP
    :param input_df: dataframe of variants to label
    :return: same database with new column 'predictive'
    """
    input_df['predictive'] = ((input_df['clin_patho'] == 'Y') |
                           ((input_df['kaviar_rare'] == 'Y') &
                           ((input_df['impact'] == 'HIGH')
                           | (input_df['dbnfsp_predicted'] == 'Y'))))

def find_AR_cands(input_df):
    """
    subset dataframe for predictive heterozygous AR variants, used in comp-het analysis
    :param input_df: dataframe of variants to subset
    :return: five data frames, one for mom, dad, nb1, nb2, and nb3
    :use: momAR, dadAR, nb1Ar,nb2Ar,nb2AR = find_AR_cands(input_df)
    """
    #subset for het AR parent variants
    mom_AR = input_df[((input_df['member'] == 'M') & (input_df['gt'] == '0/1')\
                      & (input_df['predictive'] == True) & (input_df['inheritance'] == 'AR'))]
    dad_AR = input_df[((input_df['member'] == 'F') & (input_df['gt'] == '0/1')\
                      & (input_df['predictive'] == True) & (input_df['inheritance'] == 'AR'))]
    #subset newborn het AR variants
    nb1_AR = input_df[((input_df['member'] == 'NB') & (input_df['gt'] == '0/1')\
                     & (input_df['predictive'] == True) & (input_df['inheritance'] == 'AR'))]
    nb2_AR = input_df[((input_df['member'] == 'NB2') & (input_df['gt'] == '0/1')\
                     & (input_df['predictive'] == True)& (input_df['inheritance'] == 'AR'))]
    nb3_AR = input_df[((input_df['member'] == 'NB3') & (input_df['gt'] == '0/1')\
                     & (input_df['predictive'] == True)& (input_df['inheritance'] == 'AR'))]
    return mom_AR, dad_AR, nb1_AR, nb2_AR, nb3_AR

def find_dominant_nb(input_df):
    """
    subset dataframe for predictive newborn dominant disorders
    :param input_df: data frame to subset
    :return: nb_dominant dataframe of predictive newborn variants
    """
    # subset newborn variants by variant MOI and/or zygosity
    nb_dominant = input_df[((input_df['member'] == 'NB') & (input_df['inheritance'].isin(['AD'])) \
                            & (input_df['predictive'] == True))]
    nb_dominant.name = 'dominant'
    return nb_dominant

def find_hom_rec(input_df):
    """
    subset dataframe for predictive newborn hom_rec disorders
    :param input_df: data frame to subset
    :return: nb_hom_red dataframe of homozygous alt, predictive newborn variants
    """
    # subset newborn variants by variant MOI and/or zygosity
    nb_hom_rec = input_df[((input_df['member'] == 'NB') & (input_df['inheritance'].isin(['AR'])) \
                           & (input_df['gt'] == '1/1') & (input_df['predictive'] == True))]
    nb_hom_rec.name = 'hom_recessive'
    return nb_hom_rec

# function to find matching parent variants
def find_parent_vars(nb_df, parent_df):
    # merge dataframes on by variant position
    merged_df = pd.merge(nb_df, parent_df, on=['chrom', 'pos', 'ref', 'alt'], how='inner')
    # rename parent sample_id column to avoid dropping when removing '_y' cols
    merged_df.rename(columns = {'member_y':'from_parent'}, inplace=True)
    # drop extra y columns from merge with fathers
    drop_y(merged_df)
    #remove _x from colnames
    merged_df.rename(columns=lambda x: x.replace('_x', ''), inplace=True)
    return merged_df

# run function for each group of newborns
def match_parents(nb_df):
    if (len(mom_hets) > 0) and (len(dad_hets) > 0):
        nb_and_mom = find_parent_vars(nb_df, mom_hets)
        nb_and_dad = find_parent_vars(nb_df, dad_hets)
        # merge variants found in either mother or father
        het_cands = pd.concat([nb_and_mom,nb_and_dad]).drop_duplicates().reset_index(drop=True)
        # group variants by gene name
        by_gene = het_cands.groupby(['gene_name', 'family_id'])
        return by_gene
    else:
        print "No compound het variants"