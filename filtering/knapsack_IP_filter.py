from pulp import LpVariable, LpProblem, LpMaximize, value, LpAffineExpression, LpConstraint
import re
import argparse
import itertools

FILTER_ARGS = None


def get_storf_delim(storf: list) -> list:
    colon_delimiters = [s.start() for s in re.finditer(r":",storf[0])] 
    pipe_delimiters = [s.start() for s in re.finditer(r"\|",storf[0])]
    return colon_delimiters, pipe_delimiters


def get_overlaps(uf_value_storfs: list) -> list:
    """
    Assign any StORFs that overlap on the genome into
    a distinct list. Each group of overlapping StORFs
    will be used to create constraints in the integer
    program.
    """
    overlapping_groups = []
    group_i = -1
    for index, storf in enumerate(uf_value_storfs):
        colon_delimiters, pipe_delimiters = get_storf_delim(storf)
        # the location of the StORF relative to its overlaps
        chromosome_rel_loci = (storf[0][colon_delimiters[1] + 1:pipe_delimiters[1]])
        # get _x at end of StORF location indicating its position in overlapping group
        overlap_num = int(chromosome_rel_loci[len(chromosome_rel_loci) - 1])
        # iff StORF is the start of a new overlapping group
        if overlap_num == 0:
            overlapping_groups.append([])
            group_i += 1
        overlapping_groups[group_i].append(storf)
    return overlapping_groups


def ip_set_storf_range(prob: LpProblem, storf_group: list, ip_vars: list, ip_var_id: int) -> None:
    # for each StORF, constrain length to set range
    storf_min = FILTER_ARGS.storf_range[0]
    storf_max = FILTER_ARGS.storf_range[1]
    for storf in storf_group:
        storf_len = len(re.findall('[AGTC]', storf[1]))
        prob += ip_vars[ip_var_id] * storf_len  <= storf_max
        prob += ip_vars[ip_var_id] * storf_len >= storf_min
        ip_var_id += 1
        

def ip_set_overlap_range(prob: LpProblem, storf_group: list, ip_vars: list, ip_var_id: int, group_id: int) -> None:
    overlap_min = FILTER_ARGS.overlap_range[0]
    overlap_max = FILTER_ARGS.overlap_range[1]
    # get each StORF pair combination, c, from the group
    s_pair_combination = [i for i in itertools.combinations(storf_group, 2)]  # combinations in order
    s_pair_ids = [i for i in itertools.combinations([j for j in range(0, len(storf_group))], 2)] 
    # for each c, get c[0].stop, c[1].start
    for index, s_pair in enumerate(s_pair_combination):
        # get stop location of StORF 1
        s1_colon_delimiters, s1_pipe_delimiters = get_storf_delim(s_pair[0])
        storf_1_loci = s_pair[0][0][s1_colon_delimiters[0] + 1:s1_pipe_delimiters[0]]
        storf_1_stop = int(storf_1_loci[storf_1_loci.index('-') + 1:])

        # get start location of StORF 2
        s2_colon_delimiters, s2_pipe_delimiters = get_storf_delim(s_pair[1])
        storf_2_loci = s_pair[1][0][s2_colon_delimiters[0] + 1:s2_pipe_delimiters[0]]
        storf_2_start = int(storf_2_loci[:storf_2_loci.index('-')])
        
        # if overlap, apply constraints
        if storf_1_stop - storf_2_start > 0:
            id_1 = s_pair_ids[index][0]
            id_2 = s_pair_ids[index][1]
            # overlap must be within range
            prob += (ip_vars[id_1] * storf_1_stop) - (ip_vars[id_2] * storf_2_start) <= overlap_max
            prob += (ip_vars[id_1] * storf_1_stop) - (ip_vars[id_2] * storf_2_start) >= overlap_min
            # all or none of each StORF in pair must be selected
            prob += (ip_vars[id_1] * 1) + (ip_vars[id_2] * 1) <= 2
        ip_var_id += 1


def ip_set_storf_value(prob: LpProblem, storf_group: list, ip_vars: list) -> None:
    # TODO
    return


def ip_set_no_overlaps(prob: LpProblem, storf_group: list, ip_vars: list, ip_var_id: int, group_id: int) -> None:
    # constraint group
    cg = [(i,1) for i in ip_vars[ip_var_id:ip_var_id + len(storf_group)]]
    e = LpAffineExpression(cg)
    # ensure only one StORF selected in group therefore RHS: <= 1
    c = LpConstraint(e, -1, f"no_overlap_constraint(group {group_id})", 1)
    prob += c
    

def ip_set_obj(prob: LpProblem, storf_groups: list, ip_vars: list) -> None:
    ip_var_id = 0
    obj_expression = []
    for storf_group in storf_groups:
        for storf in storf_group:
            # storf id and storf value
            obj_expression.append((ip_vars[ip_var_id], storf[2]))
            ip_var_id += 1
    # initialise objective function
    prob += LpAffineExpression(obj_expression)


def ip_set_constraints(prob: LpProblem, storf_groups: list, ip_vars: list) -> None:
    ip_var_id = 0  # track the current IP variable (StORF) being used
    # create constraints for each storf group
    for group_id, storf_group in enumerate(storf_groups):
        # TODO apply all run-time parameter constraints
        
        if FILTER_ARGS.storf_range:
            ip_set_storf_range(prob, storf_group, ip_vars, ip_var_id)
        
        if FILTER_ARGS.overlap_range:
            # TODO
            ip_set_overlap_range(prob, storf_group, ip_vars, ip_var_id, group_id)  
        else:
            ip_set_no_overlaps(prob, storf_group, ip_vars, ip_var_id, group_id)

        if FILTER_ARGS.stop_codon_type:
            return  # TODO

        ip_var_id += len(storf_group)
    

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
    # convert to list for fast lookup of selected storfs from IP solution
    storf_dict = {i: uf_value_storfs[i] for i in range(0, len(uf_value_storfs))}
    # initialise ip instance
    prob = LpProblem("StORF_Knapsack_Problem", LpMaximize)
    # initialise IP variables
    storf_ids = [f"x_{i}" for i in range(0, len(uf_value_storfs))]
    ip_vars = [LpVariable(storf_ids[i], lowBound=0, upBound=1, cat='Integer') for i in range(0, len(uf_value_storfs))]
    # get overlapping StORFs
    overlapping_storf_groups = get_overlaps(uf_value_storfs)

    # create objective function
    ip_set_obj(prob, overlapping_storf_groups, ip_vars)

    # create all constraints
    ip_set_constraints(prob, overlapping_storf_groups, ip_vars)
    
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
        value = 1

        # TODO make this the default behaviour
        if FILTER_ARGS.storf_len:
            value += int(len(storf[1]))/ls_len


        if FILTER_ARGS.simple_gc:
            # GC-content of each StORF
            value += len(re.findall('[GC]', storf[1])) / len(storf[1])
        

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

if __name__ == '__main__':
    main()
