from pulp import LpVariable, LpProblem, LpMaximize, value, LpAffineExpression, LpConstraint
import re
import argparse

FILTER_ARGS = None


def is_allowed_storf_range(storf: list) -> bool:
    storf_len = len(re.findall('[AGTC]', storf[1]))
    return FILTER_ARGS.storf_range[0] <= storf_len <= FILTER_ARGS.storf_range[0]


def get_overlaps(uf_value_storfs: list) -> list:
    """
    Assign any StORFs that overlap on the genome into
    a distinct list. Each group of overlapping StORFs
    will be used to create constraints in the integer
    program.
    """
    overlapping_groups = []
    group_i = -1
    for storf in uf_value_storfs:
        colon_delimiters = [s.start() for s in re.finditer(r":",storf[0])] 
        pipe_delimiters = [s.start() for s in re.finditer(r"\|",storf[0])]

        chromosome_UR_loci = (storf[0][colon_delimiters[1] + 1:pipe_delimiters[1]])
        overlap_num = int(chromosome_UR_loci[len(chromosome_UR_loci) - 1])
        
        # TODO implement StORF overlap min max filtering
        # if StORF is the start of a new overlapping group
        if overlap_num == 0:
            if FILTER_ARGS.storf_range:
                if is_allowed_storf_range(storf):
                    group_i += 1
                    overlapping_groups.append([])
                    overlapping_groups[group_i].append(storf)
            else:
                group_i += 1
                overlapping_groups.append([])
                overlapping_groups[group_i].append(storf)
        else:
            if FILTER_ARGS.storf_range:
                if is_allowed_storf_range(storf):
                    overlapping_groups[group_i].append(storf)
            else:
                overlapping_groups[group_i].append(storf)

    print(overlapping_groups)
    return overlapping_groups


def filter(uf_value_storfs: list) -> list:
    '''
    Filter StORFs favouring higher valued ones using integer programming to
    exclude overlapping StORFs.
    The IP model is informally defined as follows:
    Maximise Sum(all identified StORF values)
    S.t:
        sum(s ∈ S) = 1, where S is a group of overlapping StORFs
        # only one StORF in overlapping group can be selected
    with bounds:
        s ∈ {0,1}/(0 <= s <= 1, where s ∈ N+)
        # each StORF is either selected or unselected

    return: a list of filtered StORFs
    '''
    # convert to list for fast storf lookups later
    storf_dict = {i: uf_value_storfs[i] for i in range(0, len(uf_value_storfs))}
    # initialise ip instance
    prob = LpProblem("StORF_Knapsack_Problem", LpMaximize)
    # initialise ip variables
    storf_ids = [f"x_{i}" for i in range(0, len(uf_value_storfs))]
    ip_vars = [LpVariable(storf_ids[i], lowBound=0, upBound=1, cat='Integer') for i in range(0, len(uf_value_storfs))]
    # get overlapping StORFs
    overlapping_storf_groups = get_overlaps(uf_value_storfs)

    ip_var_id = 0
    obj_expression = []
    # create objective function
    for storf_group in overlapping_storf_groups:
        for storf in storf_group:
            # storf id and storf value
            obj_expression.append((ip_vars[ip_var_id], storf[2]))
            ip_var_id += 1
    # initialise objective function
    prob += LpAffineExpression(obj_expression)

    ip_var_id = 0
    # for each group of overlapping storfs, create a constraint
    for group_i, storf_group in enumerate(overlapping_storf_groups):
        c_g = []  # constraint group
        for storf_i, storf in enumerate(storf_group):
            # storf id as a coefficient, e.g. selected storfs = 1*1
            c_g.append((ip_vars[ip_var_id], 1))
            ip_var_id += 1
        # expression construction
        e = LpAffineExpression(c_g)
        # initialise constraint e = 1
        c = LpConstraint(e, 0, f"StORF_group_constraint_{group_i}", 1)
        # assign constraint to IP
        prob += c
    # solve the IP model
    prob.solve()

    # get selected StORFs
    selected = []
    for var in prob.variables():
        if value(var) == 1:
            selected.append(storf_dict[int(var.name[2:])])
    return selected


def read_fasta() -> list:
    unfiltered_storfs = []
    with open("../../testout/smallstorftest.fasta") as storf_file:
        for line in storf_file:
            if line[0] == ">":
                unfiltered_storfs.append([line, next(storf_file)])
    #with open("../../testout/E-coli_output_no_filt.fasta") as storf_file:
    #    for line in storf_file:
    #        if line[0] == ">":
    #            unfiltered_storfs.append([line, next(storf_file)])
    return unfiltered_storfs


def get_values(unfiltered_storfs: list):
    """
    Define values for each StORF according to set run-time parameters
    """
    # TODO rewrite lenght attribute already exists with the storf meta data
    largest_storf = max(unfiltered_storfs, key=lambda x: len(x[1]))
    ls_len = len(largest_storf[1])  # largest storf length
    # define value of each StORF
    for storf in unfiltered_storfs:
        value = 0
        if FILTER_ARGS.storf_len:
            value += int(len(storf[1]))/ls_len
        if FILTER_ARGS.simple_gc:
            # GC-content of each StORF
            value += len(re.findall('[GC]', storf[1])) / len(storf[1])
        # TODO add all filtering option's logic
        # GC content, overlap nt size, etc.
        storf.append(value)

    
    return unfiltered_storfs


def write_fasta(filtered_storfs: list) -> None:
    f = open("../../testout/filtered_StORFs.fasta", "w")
    for storf in filtered_storfs:
        f.write(f"{storf[0]}{storf[1]}")


def init_cl_args() -> argparse.ArgumentParser:
    """
    Run-time parameters dictate the type(s) of filtering to take place
    """
    # TODO add argument exceptions
    args = argparse.ArgumentParser(description='StORF filter options')
    args.add_argument('-olr', '--overlap_range', dest='overlap_range', 
                      nargs=2, metavar=('MIN', 'MAX'), type=int,
                      help="the min and max StORF overlap range (nt)"
                     )
    # default ~ 30, 100000
    args.add_argument('-sr' ,'--storf_range', dest='storf_range', nargs=2,
                      metavar=("MIN","MAX"), type=int,
                      help="(default: 30, 100000) the min and max StORF size (nt)"
                     )
    args.add_argument('-gc', '--gc_profile', dest='gc', nargs=2,
                      type=int, metavar=("VAR","TYPE"),
                      help='VAR: (default: 10) acceptable percentage range ' \
                      "TYPE: (default: 0) 0=mean, 1=median, 2=mode"
                     )
    args.add_argument('-sgc' ,'--simple_gc', dest='simple_gc', action='store_true',
                      help='filter favouring highest GC content StORFs when choosing between overlaps'
                     )
    args.add_argument('-sct', '--stop_codon_type', dest='stop_codon_type', nargs='?',
                      choices=('UAG', 'UGA', 'UAA'), help='filter StORFs delimited by ' \
                      'one or more of the selected stop codons'
                     ) 
    args.add_argument('-asct', '--average_stop_codon', dest='stop_codon_av',
                      metavar=("TYPE"), type=int, 
                      help=" TYPE: (default: 0) 0=mean, 1=median, 2=mode"
                     )
    args.add_argument('-slen', '--storf_length', dest="storf_len", action='store_true',
                     help="favour selection of largest StORF in overlapping StORF groups"
                     )


    global FILTER_ARGS
    FILTER_ARGS = args.parse_args()
    if not (FILTER_ARGS.overlap_range or FILTER_ARGS.storf_range or FILTER_ARGS.simple_gc
            or FILTER_ARGS.gc or FILTER_ARGS.stop_codon_type or FILTER_ARGS.stop_codon_av
            or FILTER_ARGS.storf_len):
        args.error("No filtering options provided")


def main():
    init_cl_args()
    uf_storfs = read_fasta()
    uf_value_storfs = get_values(uf_storfs)
    filtered_storfs = filter(uf_value_storfs)
    write_fasta(filtered_storfs)
    get_stats(filtered_storfs, uf_storfs)


if __name__ == '__main__':
    main()
