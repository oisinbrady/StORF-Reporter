from enum import Enum
import pulp
import re
import itertools
import collections

class HARD_FILTER(Enum):
    overlap_range = [0, 50] # default [0, 50]  #67, 1000
    size_range = [100, 50000] # default [100, 50000]
    gc_range = None# [1000, 0]  # 10=percentage variance, 0=mean, 1=median, 2=mode


class SOFT_FILTER(Enum):
    # usefull if hard filters return more candidates for a group than what 
    # the storf_group constraint allows, this will then try to favour selection
    # based of largest length, gc-content, etc.
    # (by increasing the value of the obj var coefficient)

    # TODO make soft filters for smallest storf/gc len
    storf_length = None  # True=favour the largest StORFs (more value) in groups
    gc_length = None  # True=favour the StORFs with highest gc% (more value) in groups


class WEIGHT_CONSTRAINTS(Enum):
    # knapsack constraints, default should be -1=infinity, i.e., AMAP
    sum_total = -1  # of all StORF group total weights
    storf_group = -1


def get_storf_delim(storf: list) -> list:
    colon_delimiters = [s.start() for s in re.finditer(r":",storf[0])] 
    pipe_delimiters = [s.start() for s in re.finditer(r"\|",storf[0])]
    return colon_delimiters, pipe_delimiters


def get_ave_gc(average_type: int, storfs: list) -> float:
    if average_type == 0: # mean
        total_gc = 0
        num_storfs = len(storfs)
        for storf in storfs:
            total_gc += len(re.findall('[GC]', storf[1]))
        return float(total_gc/num_storfs)
    elif average_type == 1: # median
        # order StORFs by length
        ord_storfs = sorted(storfs, key=lambda x: len(re.findall('[GC]', x[1])))
        mid = (len(ord_storfs) - 1) // 2 
        if len(ord_storfs) % 2:  # odd
            mid_gc_value = len(re.findall('[GC]', ord_storfs[mid][1]))
            return float(mid_gc_value)
        else:  # even
            # return interpolated median, average of both middle values
            mid_1_gc = len(re.findall('[GC]', ord_storfs[mid][1]))
            mid_2_gc = len(re.findall('[GC]', ord_storfs[mid + 1][1]))
            mid_gc_value = (mid_1_gc + mid_2_gc) / 2
            return mid_gc_value
    else: # mode
        # get gc count of each StORF
        s_gc_len = [len(re.findall('[GC]', s[1])) for s in storfs]
        # get the most occuring gc count
        mode_largest = max(set(s_gc_len), key=s_gc_len.count)
        # TODO currently only selects the largest mode if multiple exist
        # add some functionality for determining which mode to take?
        return float(mode_largest)



def filter_by_overlap(storf_group_values: list, storf_group: list) -> list:
    # TODO #this is not reproducing the same results as Nicks
    # Adjust objective variable coefficient to 0 iff StORF not within overlap constraint
    o_min = HARD_FILTER.overlap_range.value[0]
    o_max = HARD_FILTER.overlap_range.value[1]
    s_pair_combinations = [i for i in itertools.combinations(storf_group, 2)]  # combinations in order
    banned_storfs = []
    for i in s_pair_combinations:
        x = i[0][1][0]  # StORF x meta info
        storf_x_locus = x[x.index(":")+1:x.index("|")] 
        stop_x = int(storf_x_locus[storf_x_locus.index("-") + 1:])
        y = i[1][1][0]  # StORF y meta info
        storf_y_locus = y[y.index(":")+1:y.index("|")] 
        start_y = int(storf_y_locus[:storf_x_locus.index("-")])

        if not stop_x - start_y > 0:  # if no overlap
            continue

        x_storf_id = int(i[0][0])
        y_storf_id = int(i[1][0])
        
        if x_storf_id in banned_storfs or y_storf_id in banned_storfs:
            if y_storf_id not in banned_storfs:
                banned_storfs.append(y_storf_id)
            if x_storf_id not in banned_storfs:
                banned_storfs.append(x_storf_id)
            continue

        if stop_x - start_y <= o_min or stop_x - start_y >= o_max:
            banned_storfs.append(y_storf_id)
            banned_storfs.append(x_storf_id)

    if len(banned_storfs) > 0:
        ord_banned = sorted(banned_storfs)
        offset = ord_banned[0]
        for storf_id in banned_storfs:
            # overwrite selection value for banned StORF(s)
            storf_group_values[storf_id - offset][1] = 0
    return storf_group_values


def filter_by_size_range(storf_group_values: list, storf_group: list) -> list:
    for i, storf in enumerate(storf_group):
        minl, maxl = HARD_FILTER.size_range.value[0], HARD_FILTER.size_range.value[1]
        if not minl <= len(re.findall('[AGCT]', storf[1][1])) <= maxl:
            storf_group_values[i][1] = 0
    return storf_group_values


def filter_by_gc_range(storf_group_values: list, storf_group: list, ave_gc: float) -> list:
    mingc = ave_gc - (ave_gc * 0.05)
    maxgc = ave_gc + (ave_gc * 0.05)
    for i, storf in enumerate(storf_group):
        if not mingc <= len(re.findall('[GC]', storf[1][1])) <= maxgc:
            storf_group_values[i][1] = 0
    return storf_group_values


def filter_favour_length(storf_group_values: list, storf_group: list) -> list:
    for i, storf in enumerate(storf_group):
        # must be multiplication so that hard filters still affect coefficient values
        storf_group_values[i][1] *= len(re.findall('[AGCT]', storf[1][1]))
    return storf_group_values


def filter_favour_most_gc(storf_group_values: list, storf_group: list) -> list:

    for i, storf in enumerate(storf_group):
        storf_group_values[i][1] *= len(re.findall('[GC]', storf[1][1]))
    return storf_group_values


def set_group_values(storf_group: list, ave_gc: None | int) -> list:
    # default value is of StORF is 1 
    # If restrictions clash then overwrite offending StORF value with 0
    storf_group_values = [[s[0],1] for s in storf_group]  # storf id and value
    #print("before")
    #print(storf_group_values)
    
    # HARD FILTERS, final subset discretely defined using ranges
    if HARD_FILTER.overlap_range.value is not None:
        storf_group_values = filter_by_overlap(storf_group_values, storf_group)

    if HARD_FILTER.size_range.value is not None:
        storf_group_values = filter_by_size_range(storf_group_values, storf_group)

    if HARD_FILTER.gc_range.value is not None:
        storf_group_values = filter_by_gc_range(storf_group_values, storf_group, ave_gc)

    # TODO HARD FILTER VIA STOP CODON SELECTION
    
    # SOFT FILTERS, final subset effected by less defined StORF characteristics
    # e.g. FAVOUR larger sized StORFs(soft), only StORFs in range x-y (HARD)
    if SOFT_FILTER.storf_length.value is not None:
        storf_group_values = filter_favour_length(storf_group_values, storf_group)

    if SOFT_FILTER.gc_length.value is not None:
        storf_group_values = filter_favour_most_gc(storf_group_values, storf_group)

    #print("after")
    #print(storf_group_values)
    #b = storf_group[0][1][0].find("4970-5284")
    #if b != -1:
    #    print(storf_group[0][1][0])
    #    print(storf_group_values)
    return storf_group_values


def is_new_group(storf: list, storf_id: int, total_num_storfs: int) -> bool:
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
    with open("../../testin/E-coli_output_no_filt.fasta") as storf_file:
        for line in storf_file:
            if line[0] == ">":
                unfiltered_storfs.append([line, next(storf_file)])
    # TMP test files
    # "../../testin/smallstorftest.fasta"
    # "../../testin/E-coli_output_no_filt.fasta"
    return unfiltered_storfs


def write_fasta(filtered_storfs: list) -> None:
    f = open("../../testout/oisin_output/output.fasta", "w")
    for storf in filtered_storfs:
        f.write(f"{storf[0]}{storf[1]}")


def ip_set_obj_func(prob: pulp.LpProblem, obj_variables: list, ip_vars: list) -> None:
    obj_expression = []
    #for tmp in range(0,4):
        #print(obj_variables[tmp])
    
    for variable in obj_variables:
        obj_expression.append((ip_vars[variable[0]], variable[1]))

    #for tmp in range(0,4):
    #    print(obj_expression[tmp])
    e = pulp.LpAffineExpression(obj_expression)
    # add objective function (in turn, adds variable bounds)
    prob += e


def ip_set_group_constraint(prob: pulp.LpProblem, ip_vars: list, group: list, group_id: int) -> None:
    if WEIGHT_CONSTRAINTS.storf_group.value == -1:
        # iff infinity (i.e. AMAP StORFs in each group), no need for constraint
        return
    else:
        # constraint group, where weight of each var = 1
        C_g = WEIGHT_CONSTRAINTS.storf_group.value
        g = []
        for i, storf in enumerate(group):
            g.append((ip_vars[storf[0]], 1))
        e = pulp.LpAffineExpression(g)
        # RHS: <= C_g
        c = pulp.LpConstraint(e, -1, f"internal_knapsack_constraint_{group_id}", C_g)
        prob += c


def ip_set_total_constraint(prob: pulp.LpProblem, ip_vars: list) -> None:
    # constrain the sum of all knapsack weights 
    # I.e, of all StORFs in file to be selected
    if WEIGHT_CONSTRAINTS.sum_total.value == -1:
        return
    else:
        C_t = WEIGHT_CONSTRAINTS.sum_total.value
        g = [(i,1) for i in ip_vars]
        e = pulp.LpAffineExpression(g)
        # RHS: <= C_t
        c = pulp.LpConstraint(e, -1, f"external_knapsack_constraint", C_t)
        prob += c


def ip_filter(storfs: list) -> list:

    # N.B.: Weight of each StORF = 1,
    # C_t = cap of all StORFs e.g., how many StORFs can be selected total,
    # C_k = cap selection for each StORF group (sub-knapsack). 

    # initialise IP instance
    prob = pulp.LpProblem("StORF_IP_Filter", pulp.LpMaximize)
    s_total = len(storfs)
    g = 0  # group id
    s = 0  # storf id
    # Create IP variables
    storf_ids = [f"x_{s}" for s in range(0, s_total)]
    ip_vars = [pulp.LpVariable(storf_ids[i], lowBound=0, upBound=1, cat='Integer') for i in range(0, s_total)]
    group = []  # overlapping group of StORFs
    same_group = True
    obj_variables = []


    # pre-compute average gc content
    if HARD_FILTER.gc_range.value is not None: 
        typ = HARD_FILTER.gc_range.value[1]  # 0=mean, 1=median, 2=mode
        ave_gc = get_ave_gc(typ, storfs)
    else:
        ave_gc = None
    

    # determine coefficient values of future IP variables
    for storf in storfs:
        # StORF values are dependent on their group
        if is_new_group(storf, s, s_total):
            obj_variables += set_group_values(group, ave_gc)
            # add weight constraint to each group
            ip_set_group_constraint(prob, ip_vars, group, g)
            group = [(s, storf)]  # init new group
            g += 1
        else:
            group.append((s, storf))
        s+=1
    # add values to last group of StORFs...  
    obj_variables += set_group_values(group, ave_gc)
    # set constraint for last group (auto-adds IP variable bounds)
    ip_set_group_constraint(prob, ip_vars, group, g)

    # construct objective function (auto-adds IP variable bounds)
    ip_set_obj_func(prob, obj_variables, ip_vars)
    
    # Add knapsack sum constraint
    ip_set_total_constraint(prob, ip_vars)

    # get selected StORFs from IP solution
    prob.solve()

    selected = {}
    for var in prob.variables():
        print(f"{var.name}={pulp.value(var)}")
        selected[var.name] = pulp.value(var)

    ordered_selected = collections.OrderedDict(sorted(selected.items(), key=lambda t: int(t[0][2:]) ))
    final_filter = []
    for i, v in enumerate(ordered_selected.values()):
        if v == 1.0:
            final_filter.append(storfs[i])
        
    return final_filter


def main():
    storfs = read_fasta()
    storfs = ip_filter(storfs)
    write_fasta(storfs)

if __name__ == '__main__':
    main()