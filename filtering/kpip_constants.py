# constants based of PCA/metrics from ur_metrics.py
HSS_P = 0.074  # average proportion of HSS found in the test genomes
GC_LB_P = 0.35  # lower bound of noramlised GC proportion where HSS_1s are found
GC_UB_P = 0.85
OLAP_SIZE_LB_P = 0.2  # lower bound of normalised average StORF overlap size where HSS_1s are found
MIN_NORM_SIZE = 0.1  # majority of HSS_1 StORFs are larger than this average normalised StORF size
ARB_MAX_STORF_SIZE = 50000  # based of original StORF_finder.py
ARB_MAX_OLAP_SIZE = 50000
