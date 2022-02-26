from enum import Enum
import pulp
import re
import itertools

class HARD_FILTER(Enum):
    overlap_range = [0, 1000]  #67, 1000
    size_range = [0, 100]
    gc_range = [0, 2]  # 10=percentage variance, 0=mean, 1=median, 2=mode


class SOFT_FILTER(Enum):
    storf_length = False  # True=favour the largest StORFs (more value) in groups
    gc_length = False  # True=favour the StORFs with highest gc% (more value) in groups


class WEIGHT_CONSTRAINTS(Enum):
    # knapsack constraints, -1=infinity, AMAP
    sum_total = -1  # of all StORF group total weights
    storf_group = -1  


def get_storf_delim(storf: list) -> list:
    colon_delimiters = [s.start() for s in re.finditer(r":",storf[0])] 
    pipe_delimiters = [s.start() for s in re.finditer(r"\|",storf[0])]
    return colon_delimiters, pipe_delimiters


def filter_by_overlap(storf_group_values: list, storf_group: list) -> list:
    # Adjust objective variable coefficient to 0 iff StORF not within overlap constraint
    o_min = HARD_FILTER.overlap_range.value[0]
    o_max = HARD_FILTER.overlap_range.value[1]
    s_pair_combinations = [i for i in itertools.combinations(storf_group, 2)]  # combinations in order
    banned_storfs = []
    for i in s_pair_combinations:
        x_storf_id = int(i[0][0])
        y_storf_id = int(i[1][0])
        if x_storf_id in banned_storfs or y_storf_id in banned_storfs:
            banned_storfs.append(x_storf_id)
            banned_storfs.append(y_storf_id)
            continue
        x = i[0][1][0]
        storf_x_locus = x[x.index(":")+1:x.index("|")] 
        stop_x = int(storf_x_locus[storf_x_locus.index("-") + 1:])
        y = i[1][1][0]
        storf_y_locus = y[y.index(":")+1:y.index("|")] 
        start_y = int(storf_y_locus[:storf_x_locus.index("-")])
        if stop_x - start_y > 0:
            if stop_x - start_y < o_min or stop_x - start_y > o_max:
                banned_storfs.append(y_storf_id)
                banned_storfs.append(x_storf_id)
    for b in banned_storfs:
        storf_group_values[b] = (b, 0) 
    return storf_group_values


def filter_by_size_range(storf_group_values: list, storf_group: list) -> list:
    return storf_group_values


def filter_by_gc_range(storf_group_values: list, storf_group: list) -> list:
    return storf_group_values


def filter_favour_length(storf_group_values: list, storf_group: list) -> list:
    return storf_group_values


def filter_favour_gc_percentage(storf_group_values: list, storf_group: list) -> list:
    return storf_group_values


def set_group_values(storf_group: list) -> list:
    # default value is of StORF is 1 
    # If restrictions clash then overwrite offending StORF value with 0
    storf_group_values = [(s[0],1) for s in storf_group]

    #print("before")
    #print(storf_group_values)
    if len(storf_group_values) != 1:
        if HARD_FILTER.overlap_range.value is not None:
            storf_group_values = filter_by_overlap(storf_group_values, storf_group)
            
        #TODO# if HARD_FILTER.size_range.value is not None:
            # storf_group_values = filter_by_size_range(storf_group_values, storf_group)

        #TODO# if HARD_FILTER.gc_range.value is not None:
            # storf_group_values = filter_by_gc_range(storf_group_values, storf_group)
        
        #TODO# if SOFT_FILTER.storf_length.value is not None:
            # storf_group_values = filter_favour_length(storf_group_values, storf_group)

        #TODO# if SOFT_FILTER.gc_length.value is not None:
            # storf_group_values = filter_favour_gc_percentage(storf_group_values, storf_group)

    #print("after")
    #print(storf_group_values)
    return storf_group_values


def is_new_group(storf: list, storf_id, total_num_storfs: int) -> bool:
    # Does the StORF belong to a new collection of overlapping StORFs
    colon_delims, pipe_delims = get_storf_delim(storf)
    # the location of the StORF relative to its overlapping region
    chromosome_rel_loci = (storf[0][colon_delims[1] + 1:pipe_delims[1]])
    # get _x at end of StORF location indicating its position in overlapping group
    overlap_num = int(chromosome_rel_loci[len(chromosome_rel_loci) - 1])
    if storf_id + 1 == total_num_storfs:
        return True
    return overlap_num == 0 and storf_id != 0


def read_fasta() -> list:
    unfiltered_storfs = []
    with open("../../testin/smallstorftest.fasta") as storf_file:
        for line in storf_file:
            if line[0] == ">":
                unfiltered_storfs.append([line, next(storf_file)])
    # TMP test files
    # "../../testout/smallstorftest.fasta"
    # "../../testout/E-coli_output_no_filt.fasta"
    return unfiltered_storfs


def write_fasta(filtered_storfs: list) -> None:
    f = open("../../testout/oisin_output/output.fasta", "w")
    for storf in filtered_storfs:
        f.write(f"{storf[0]}{storf[1]}")


def ip_set_obj_func(prob: pulp.LpProblem, obj_variables: list, ip_vars: list) -> None:
    obj_expression = []
    for variable in obj_variables:
        obj_expression.append((ip_vars[variable[0]], variable[1]))
    e = pulp.LpAffineExpression(obj_expression)
    # add objective function (in turn, adds variable bounds)
    prob += e


def ip_filter(storfs: list) -> list:

    # N.B., Weight of each StORF = 1
    # C_t = cap of all StORFs e.g., how many StORFs can be selected total
    # C_k = cap selection for each StORF group (sub-knapsack)   

    # initialise IP instance
    prob = pulp.LpProblem("StORF_IP_Filter", pulp.LpMaximize)
    s_total = len(storfs)
    g = 0  # group id
    s = 0  # storf id
    # initialise ILP variables
    storf_ids = [f"x_{s}" for s in range(0, s_total)]
    ip_vars = [pulp.LpVariable(storf_ids[i], lowBound=0, upBound=1, cat='Integer') for i in range(0, s_total)]
    group = []  # overlapping group of StORFs
    same_group = True
    obj_variables = []
    # determine values of future IP variables
    for storf in storfs:
        # many StORF values are dependent on their group
        if is_new_group(storf, s, s_total):
            obj_variables += set_group_values(group)
            group = [(s, storf)]  # init new group
        else:
            group.append((s, storf))
        s+=1 
    # add values to last group of StORFs...  
    obj_variables += set_group_values(group)

    # construct objective function
    ip_set_obj_func(prob, obj_variables, ip_vars)
    

    # TODO # Add knapsack constraints C_t, C_k

    # get selected StORFs from IP solution
    prob.solve()
    selected = []
    for i, var in enumerate(prob.variables()):
        print(f"{i}={pulp.value(var)}")

        if pulp.value(var) == 1.0:
            selected.append(storfs[i])
    return selected


def main():
    storfs = read_fasta()
    storfs = ip_filter(storfs)
    write_fasta(storfs)

if __name__ == '__main__':
    main()