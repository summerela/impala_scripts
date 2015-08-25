import pandas as pd

cand_script = pd.DataFrame(open('/Users/selasady/impala_annot/candidate_patterns_for_itmi.103-00001_at_0713190239.txt', 'rw').readlines())

print cand_script.head(15)