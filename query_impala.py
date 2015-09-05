#! /usr/bin/env python
# Filename : query_impala.py

''' Add documentation on overall purpose here.
'''

#create lists of regular and wildcard args
def find_wildcards(user_arg):
    '''separate wildcards from regular args
    in comma-separated list of user args
    to create statements with either LIKE or = '''
    user_args = []
    wildcards = []
    # seperate user args by comma
    arg_list = user_arg.replace("'", "").split(',')
    for arg in arg_list:
        if arg.endswith('%'):
            wildcards.append(arg)
        else:
            user_args.append(arg)
    return user_args, wildcards

# function to turn arg_list and wild_list into statements
def process_args(arg_list, wild_list, name_arg, arg_db):
    '''
    :param arg_list: input list of non-wildcard user arguments
    :param wild_list: input list of wildcard arguments
    :param name_arg: argument name, such as 'genes'
    :param arg_db: name of annotation or search database
    :return: user arguments converted into query statements with LIKE for wildcards and = for non-wildcards
    '''
    # create statement for one user arg and no wildcards
    if (len(arg_list) == 1 and len(wild_list) == 0):
        query_args = "AND {}.{} = '{}'".format(arg_db, name_arg, ', '.join(str(e) for e in arg_list))
    # create statement for more than one user argument and no wildcards
    elif (len(arg_list) > 1 and len(wild_list) == 0):
        query_args = "AND ({})".format(" OR ".join("{}.{} = '{}'".format(arg_db, name_arg, s) for s in arg_list))
    # create statement for one wildcard and no regular args
    elif (len(wild_list) == 1 and len(arg_list) == 0):
         query_args = "AND {}.{} LIKE '{}'".format(arg_db, name_arg, ', '.join(str(e) for e in wild_list))
    # create statement for more than one wildcard and no regular args
    elif (len(wild_list) > 1 and len(arg_list) == 0):
        query_args = "AND ({})".format(" OR ".join("{}.{} LIKE '{}'".format(arg_db, name_arg, s) for s in wild_list))
    # create statement for more than one wildcard and more than one user arg
    elif (len(arg_list) > 1 and len(wild_list) > 1):
        reg_args = "AND ({})".format(" OR ".join("{}.{} = '{}'".format(arg_db, name_arg, s) for s in arg_list))
        wild_args = "AND ({})".format(" OR ".join("{}.{} LIKE '{}'".format(arg_db, name_arg, s) for s in wild_list))
        query_args = (reg_args + " " + wild_args).replace(' AND ', ' OR ', 2) + ")"
    # create statement for more than one of each
    elif (len(arg_list) == 1 and len(wild_list) == 1):
        reg_args = "AND {}.{} = '{}'".format(arg_db, name_arg, ', '.join(str(e) for e in arg_list))
        wild_args = "AND {}.{} LIKE '{}'".format(arg_db, name_arg, ', '.join(str(e) for e in wild_list))
        query_args = (reg_args + " " + wild_args).replace(') AND (', ' OR ', 1)
    #if none of the above match, something is wrong
    else:
        print "Length of arg_list = " + str(len(arg_list))
        print "Length of wild_list = " + str(len(wild_list))
        print "Check your {0} argument and try again.".format(name_arg)
    return query_args

# merge wildcard and argument list for each user arg
def create_query(argument, arg_name, arg_db):
    '''
    :param argument: specify user argument to run through function, such as 'genes'
    :param arg_name: the actual column name of the field to input into a query, such as 'gt' for genotype
    :param arg_db: specify annot_db or search_db, wherever this value will be coming from
    :return: returns a combined list of regular and wildcard statements for each user argument
    '''
    if argument != 'all':
        arg_list, wild_list = find_wildcards(argument)
        arg_query = process_args(arg_list, wild_list, arg_name, arg_db)
    else:
        #return an empty string if the user arg is 'all'
        arg_query = ''
    return arg_query

# merge statements for all user arguments
def merge_query(query_arg, merge_list):
    '''
    :param query_arg: list of combined user/wildcard args object created from create_query() to be merged into a final
    query statement
    :param merge_list: name of empty list to hold final merged query
    :return: final merged query composed of all regular and wildcard statemetns for each user argument
    '''
    if query_arg:
        merge_list.append(query_arg)
    impala_query= ' '.join(merge_list).replace('AND', 'WHERE', 1)
    return impala_query

# create list of columns to return from the query
def get_cols(db, col_arg):
    '''
    :param db: name of database to return columns from (search_db or annot_db)
    :param col_arg: name of variable containing a list of column names to return
    :return: list of columns to add to a SELECT columns FROM statement
    '''
    if col_arg == 'all':
        col_query = '{0}.*'.format(db)
    else:
        col_list = col_arg.replace("'", "").split(',')
        if col_list:
            for col in col_list:
                cols = ["{0}.".format(db) + s for s in col_list]
                col_query = ",".join(map(str, cols))
        elif len(col_list) == 1:
            col_query = "{0}.".format(db) + ",".join(map(str, col_list))
        else:
            print "Check your search_db_cols and annot_db_cols arguments and try again please."
            col_query = ''
    return col_query

# check if tables contain start/stop or pos, and ref columns
def check_colnames(cursor, table_name, print_out=True):
    '''
    :param cursor: impala database connection object
    :param table_name: name of table to check
    :param print_out: print out the list to screen
    :return: a list of columns in the specified impala table
    '''
    cur.execute('DESCRIBE {0}'.format(table_name))
    info = cur.fetchall()
    cols = []
    for col in info:
         cols.append(col[0])
    return cols

# function to match cgi with annot_db
def join_annot_to_cgi(annot_db):
    '''
    :param annot_db: name of annotation source, example 'ensemble_genes'
    :return: returns query statements used to match tables on position/ref/alt etc.
    '''
    # check for ref columns
    if 'ref' in annot_cols:
        # if there is a start column
        if 'start' in annot_cols:
            # if there is a stop position reported
            if 'stop' in annot_cols:
                using_statement = "USING (start,stop,ref,alt)"
                match_statement = " AND {0}.start = {1}.start AND {0}.stop = {1}.stop AND {0}.ref = {1}.ref AND { \
                                  0}.alt = {1}.alt ".format(search_db, annot_db)
            # there is start but no stop position
            else:
                using_statement = "USING (start,ref,alt)"
                match_statement = " AND {0}.start = {1}.start AND {0}.ref = {1}.ref AND { \
                                  0}.alt = {1}.alt ".format(search_db, annot_db)
        # there is pos column instead of start
        elif 'pos' in annot_cols:
            using_statement = "USING (pos,ref,alt)"
            match_statement = " AND ({1}.pos BETWEEN {0}.start AND {0}.stop) AND {0}.ref = {1}.ref \
            AND {0}.alt = {1}.alt ".format(search_db, annot_db)
    elif 'ref' not in annot_cols:
        # if there is a start column
        if 'start' in annot_cols:
            # if there is a stop position reported
            if 'stop' in annot_cols:
                using_statement = "USING (start,stop)"
                match_statement = " AND ({0}.start >= {1}.start AND {0}.stop <= {1}.stop) ".format(search_db, annot_db)
            # there is start but no stop position
            else:
                using_statement = "USING (start)"
                match_statement = " AND {0}.start = {1}.start ".format(search_db, annot_db)
        # there is pos column instead of start
        elif 'pos' in annot_cols:
            using_statement = "USING (pos)"
            match_statement = " AND ({1}.pos BETWEEN {0}.start AND {0}.stop) ".format(search_db, annot_db)
    # if no start, stop or pos, match on gene name
    elif 'gene_name' in annot_cols:
        using_statement = "USING (gene_name)"
        match_statement = " AND {0}.gene_name = {1}.gene_name ".format(search_db, annot_db)
    else:
        print "Cannot find method for joining tables."
    return using_statement,match_statement

# function to match illumina with annot_db
def join_annot_to_illumina(annot_db):
    '''
    :param annot_db: name of annotation source, example 'ensemble_genes'
    :return: returns query statements used to match tables on position/ref/alt etc.
    '''
    # check for ref columns
    if 'ref' in annot_cols:
        # if there is a start column
        if 'start' in annot_cols:
            # if there is a stop position reported
            if 'stop' in annot_cols:
                using_statement = "USING (start,stop,ref,alt)"
                match_statement = " AND ({0}.pos BETWEEN {1}.start AND {1}.stop) AND {0}.ref = {1}.ref \
                                  AND {0}.alt = {1}.alt ".format(search_db, annot_db)
            # there is start but no stop position
            else:
                using_statement = "USING (start,ref,alt)"
                match_statement = " AND {0}.pos = {1}.start AND {0}.ref = {1}.ref AND { \
                                  0}.alt = {1}.alt ".format(search_db, annot_db)
        # there is pos column instead of start
        elif 'pos' in annot_cols:
            using_statement = "USING (pos,ref,alt)"
            match_statement = " AND {0}.pos = {1}.pos AND {0}.ref = {1}.ref \
            AND {0}.alt = {1}.alt ".format(search_db, annot_db)
    elif 'ref' not in annot_cols:
        # if there is a start column
        if 'start' in annot_cols:
            # if there is a stop position reported
            if 'stop' in annot_cols:
                using_statement = "USING (start,stop)"
                match_statement = " AND ({0}.pos BETWEEN {1}.start AND {1}.stop) ".format(search_db,annot_db)
            # there is start but no stop position
            else:
                using_statement = ""
                match_statement = " AND {0}.pos = {1}.start ".format(search_db, annot_db)
        # there is pos column instead of start
        elif 'pos' in annot_cols:
            using_statement = "USING (pos)"
            match_statement = " AND {0}.pos = {1}.pos ".format(search_db, annot_db)
    # if no start, stop or pos, match on gene name
    elif 'gene_name' in annot_cols:
        using_statement = "USING (gene_name)"
        match_statement = " AND {0}.gene_name = {1}.gene_name ".format(search_db, annot_db)
    else:
        print "Cannot find method for joining tables."
        using_statement = ""
        match_statement = ""
    return using_statement, match_statement

