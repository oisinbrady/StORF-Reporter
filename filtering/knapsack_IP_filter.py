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
        # the location of the StORF relative to its overlapping region
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
        prob += ip_vars[ip_var_id] * storf_len <= storf_max
        prob += ip_vars[ip_var_id] * storf_len >= storf_min
        ip_var_id += 1
        

def ip_set_overlap_range(prob: LpProblem, storf_group: list, ip_vars: list, ip_var_id: int, group_id: int) -> None:
    overlap_min = FILTER_ARGS.overlap_range[0]
    overlap_max = FILTER_ARGS.overlap_range[1]
    # get each StORF pair combination, c, from the group
    s_pair_combination = [i for i in itertools.combinations(storf_group, 2)]  # combinations in order
    s_pair_ids = [i for i in itertools.combinations([j for j in range(0, len(storf_group))], 2)] 

    for index, s_pair in enumerate(s_pair_combination):
        # for each s_pair, get s_pair[0].stop, s_pair[1].start
        # get stop location of StORF 1
        s1_colon_delimiters, s1_pipe_delimiters = get_storf_delim(s_pair[0])
        storf_1_loci = s_pair[0][0][s1_colon_delimiters[0] + 1:s1_pipe_delimiters[0]]
        storf_1_stop = int(storf_1_loci[storf_1_loci.index('-') + 1:])

        # get start location of StORF 2
        s2_colon_delimiters, s2_pipe_delimiters = get_storf_delim(s_pair[1])
        storf_2_loci = s_pair[1][0][s2_colon_delimiters[0] + 1:s2_pipe_delimiters[0]]
        storf_2_start = int(storf_2_loci[:storf_2_loci.index('-')])
        
        # if overlap, apply constraints
        # TODO REFACTOR - MIN CONSTRAINT NOT WORKING PROPERLY
        if storf_1_stop - storf_2_start > 0:
            # determine the IP variable IDs
            id_1 = s_pair_ids[index][0] + ip_var_id
            id_2 = s_pair_ids[index][1] + ip_var_id
            # CONSTRAINTS (1) & (2): overlap must be within range
            prob += ip_vars[id_1] * storf_1_stop - ip_vars[id_2] * storf_2_start <= overlap_max
            prob += (ip_vars[id_1] * storf_1_stop) - (ip_vars[id_2] * storf_2_start) >= overlap_min
            
    
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


def non_ip_filter(storf_list) -> list:
    """
    Currently includes filtering by stop codon type, and via average gc profile
    """
    if FILTER_ARGS.gc_profile:
        ip_var_id = 0
        gc_average = get_gc_average(storf_list)
        percentage_variance = FILTER_ARGS.gc_profile[0]/100
        gc_average_min = gc_average - ((percentage_variance/2) * gc_average)
        gc_average_max = gc_average + ((percentage_variance/2) * gc_average)
    
    for index, storf in enumerate(storf_list):        
        if FILTER_ARGS.gc_profile: # filter according to average of gc_profile
            # get gc_content of StORF
            storf_gc_len = len(re.findall('[GC]', storf[1]))
            # remove any StORF not within the deviation range of the input file's average gc value
            if not gc_average_min <= storf_gc_len <= gc_average_max:
                del storf_list[index]
                continue
        if FILTER_ARGS.stop_codon_type:  # filter according to allowed stop codon types
            stop_codon_1 = storf[1][0:3] 
            stop_codon_2 = storf[1][len(storf[1]) - 4:len(storf[1])]
            remove = True
            for sct in FILTER_ARGS.stop_codon_type:
                if sct == stop_codon_1 or sct == stop_codon_2:
                    remove = False
            if remove:
                del storf_list[index]
    return storf_list

def ip_set_constraints(prob: LpProblem, storf_groups: list, ip_vars: list) -> None:
    # create constraints for each storf
    ip_var_id = 0  # track the current IP variable (StORF) being used
    
    for group_id, storf_group in enumerate(storf_groups):
        
        if FILTER_ARGS.storf_range:
            ip_set_storf_range(prob, storf_group, ip_vars, ip_var_id)
        
        if FILTER_ARGS.overlap_range:
            ip_set_overlap_range(prob, storf_group, ip_vars, ip_var_id, group_id)  
        else:
            ip_set_no_overlaps(prob, storf_group, ip_vars, ip_var_id, group_id)

        ip_var_id += len(storf_group)
    

def ip_filter(uf_value_storfs: list) -> list:
    '''
    Filter StORFs using integer programming
    The IP model is informally defined as follows:
    Maximise Sum(all identified StORF values)
    S.t:  # !!!TODO!!! add each potential model in documentation
        (1)
        (2)
    with bounds:
        (a)

    return: a list of filtered StORFs
    '''
    # convert to list for lookup of selected storfs from IP solution
    storf_dict = {i: uf_value_storfs[i] for i in range(0, len(uf_value_storfs))}
    # initialise IP instance
    prob = LpProblem("StORF_Knapsack_Problem", LpMaximize)
    # initialise IP variables with bounds

    # TODO make bound low=0, upper=1 if filtering for specific stop codons and storf stop codon is not allowed
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

    # get selected StORFs from IP solution
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


def get_gc_average(storf_list: list) -> float: 
    average_type = FILTER_ARGS.gc_profile[1]  # 0=mean, 1=median, 2=mode
    if average_type == 0: # mean
        total_gc = 0
        num_storfs = len(storf_list)
        for storf in storf_list:
            total_gc += re.findall('[GC]', storf[1])
        return total_gc/num_storfs
    elif average_type == 1: # median
        # order StORFs by length
        ord_storfs = sorted(storf_list, key=lambda x: len(re.findall('[GC]', x[1])))
        mid = (len(ord_storfs) - 1) // 2
        if len(ord_storfs) % 2:  # odd
            mid_gc_value = len(re.findall('[GC]', ord_storfs[mid][1]))
            return mid_gc_value
        else:  # even
            # return interpolated median, average of both middle values
            mid_1_gc = len(re.findall('[GC]', ord_storfs[mid][1]))
            mid_2_gc = len(re.findall('[GC]', ord_storfs[mid + 1][1]))
            mid_gc_value = (mid_1_gc + mid_2_gc) / 2.0
            return mid_gc_value
    elif average_type == 2:  # mode
        # get gc count of each StORF
        s_gc_len = [len(re.findall('[GC]', s[1])) for s in storf_list]
        # get the most occuring gc count
        mode_largest = max(set(s_gc_len), key=s_gc_len.count)

        # TODO currently only selects the largest mode if multiple exist
        # add some functionality for determining which mode to take?
        return mode_largest


def get_values(unfiltered_storfs: list):
    """
    Define values for each StORF according to set run-time parameters. 
    These values will be added to the coefficients of each StORF variable 
    in the maximisation objective function of the integer program.
    """
    # TODO rewrite lenght attribute already exists with the storf meta data
    largest_storf = max(unfiltered_storfs, key=lambda x: len(x[1]))
    ls_len = len(largest_storf[1])  # largest storf length

    # define value of each StORF
    for storf in unfiltered_storfs:
        value = 1

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
    args.add_argument('-gc', '--gc_profile', dest='gc_profile', nargs=2,
                      type=int, metavar=("VAR","TYPE"),
                      help='VAR: (default: 10) acceptable percentage range ' \
                      'TYPE: (default: 0) 0=mean, 1=median, 2=mode. Filter ' \
                      'favouring StORFs closest to average gc content'
                     )
    args.add_argument('-sgc' ,'--simple_gc', dest='simple_gc', action='store_true',
                      help='filter favouring highest GC content StORFs when choosing between overlaps'
                     )
    args.add_argument('-sct', '--stop_codon_type', dest='stop_codon_type', nargs='+',
                      choices=('TAG', 'TAA', 'TGA'), help='filter StORFs delimited by ' \
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
            or FILTER_ARGS.gc_profile or FILTER_ARGS.stop_codon_type or FILTER_ARGS.stop_codon_av
            or FILTER_ARGS.storf_len):
        args.error("No filtering options provided")


def main():
    init_cl_args()
    uf_storfs = read_fasta()
    uf_value_storfs = get_values(uf_storfs)  # unfiltered StORFs with heuristic values
    ip_filtered_storfs = ip_filter(uf_value_storfs)
    fully_filtered_storfs = non_ip_filter(ip_filtered_storfs)
    write_fasta(fully_filtered_storfs)

if __name__ == '__main__':
    main()
