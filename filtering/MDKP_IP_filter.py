from enum import Enum
import pulp
import re
import itertools
import collections

# filters to match default behaviour of Nick's program is:
# storf_length=True
# overlap_range = [0, 50]
# size_range = [100, 50000]

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
    storf_length = True  # True=favour the largest StORFs (more value) in groups
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
    # modified version of Nick's tile_filtering function
    o_min, o_max = HARD_FILTER.overlap_range.value[0], HARD_FILTER.overlap_range.value[1]
    ################ - Order greatest value first filtering
    ordered_by_value = []
    sorted_values = sorted(storf_group_values, key=lambda storf: storf[1], reverse=True)
    offset = storf_group_values[0][0]
    for s in sorted_values:
        ordered_by_value.append(storf_group[s[0]-offset])
    ############## - For each StORF, remove all smaller overlapping STORFs according to filtering rules
    length = len(ordered_by_value)
    i = 0
    while i < length:
        x_index = ordered_by_value[i][0] - offset
        storf_x_meta = storf_group[x_index][1][0]
        storf_x_id = storf_group[i][0]
        x_locus = storf_x_meta[storf_x_meta.find(":")+1:storf_x_meta.find("|")]
        start_x = int(x_locus[:x_locus.find("-")])
        stop_x = int(x_locus[x_locus.find("-")+1:])
        j = i+1
        while j < length:
            y_index = ordered_by_value[j][0] - offset
            storf_y_meta = storf_group[y_index][1][0]
            y_locus = storf_y_meta[storf_y_meta.find(":")+1:storf_y_meta.find("|")]
            start_y = int(y_locus[:y_locus.find("-")])
            stop_y = int(y_locus[y_locus.find("-")+1:])
            if start_y >= stop_x or stop_y <= start_x:
                j+=1
                continue  # Not caught up yet / too far
            elif start_y >= start_x and stop_y <= stop_x:
                ordered_by_value.pop(j)
                length = len(ordered_by_value)
            else: # +1 needed for stop codon
                x = set(range(start_x,stop_x+1))
                y = set(range(start_y,stop_y+1))
                overlap = len(x.intersection(y))
                if overlap >= o_max:
                    ordered_by_value.pop(j)
                    length = len(ordered_by_value)
                else:
                    j += 1
        length = len(ordered_by_value)
        i+=1

    selected = [i[0] for i in ordered_by_value]
    
    for storf in storf_group_values:
        if storf[0] not in selected:
            storf[1] = 0
        
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


def storf_get_value(storf: list) -> int:
    value = 0
    if SOFT_FILTER.storf_length.value is not None:
        value += len(re.findall('[AGCT]', storf[1]))

    if SOFT_FILTER.gc_length.value is not None:
        value += len(re.findall('[GC]', storf[1]))
    return value


def set_group_values(storf_group: list, ave_gc: None | int) -> list:
    # default value is of StORF is 1 
    # If restrictions clash then overwrite offending StORF value with 0
    storf_group_values = [[s[0],1] for s in storf_group]  # storf id and value
    #print("before")
    #print(storf_group_values)
    
    # SOFT FILTERS, final subset effected by less defined StORF characteristics
    # e.g. FAVOUR larger sized StORFs(soft), only StORFs in range x-y (HARD)
    if SOFT_FILTER.storf_length.value is not None:
        storf_group_values = filter_favour_length(storf_group_values, storf_group)

    if SOFT_FILTER.gc_length.value is not None:
        storf_group_values = filter_favour_most_gc(storf_group_values, storf_group)



    # TODO convert to constraints
    # HARD FILTERS, final subset discretely defined using ranges
    #if HARD_FILTER.overlap_range.value is not None:
    #    storf_group_values = filter_by_overlap(storf_group_values, storf_group)

    #if HARD_FILTER.size_range.value is not None:
    #    storf_group_values = filter_by_size_range(storf_group_values, storf_group)

    #if HARD_FILTER.gc_range.value is not None:
    #    storf_group_values = filter_by_gc_range(storf_group_values, storf_group, ave_gc)

    # TODO HARD FILTER VIA STOP CODON SELECTION
    

    #print("after")
    #print(storf_group_values)
    return storf_group_values


def ip_set_constraints(prob: pulp.LpProblem, storfs: list, total_storfs: int, ip_vars: list) -> list:
    group = []
    g = 0
    s = 0
    group_tracker = [[]]
    # get StORF "overlap" groups
    for storf in storfs:
        if is_next_new_group(storf, s, total_storfs):
            # ip_set_group_constraint(prob, ip_vars, group, g)
            ip_set_overlap_constraint(prob, group, total_storfs, ip_vars)

            group_tracker.append([s])
            group = [(s, storf)]  # init new group
            g += 1
        else:
            group_tracker[g].append(s)
            group.append((s, storf))

        s+=1
    return group_tracker
    # ip_set_group_constraint(prob, ip_vars, group, g)


def ip_set_overlap_constraint(prob: pulp.LpProblem, group: list, total_storfs: int, ip_vars: list):
    max_overlap = HARD_FILTER.overlap_range.value[1]
    # get each pair combination from StORF group
    pairs = [i for i in itertools.combinations(group, 2)]  # combinations in order
    if max_overlap == -1 and len(group) > 1:
        # iff no max overlap set, no need for constraint
        return
    elif len(group) > 1:
        for i, pair in enumerate(pairs):
            # get x stop location
            storf_x_meta = pair[0][1][0]
            storf_x_id = pair[0][0]
            x_locus = storf_x_meta[storf_x_meta.find(":")+1:storf_x_meta.find("|")]
            stop_x = int(x_locus[x_locus.find("-")+1:])
            # get y start location
            storf_y_meta = pair[1][1][0]
            storf_y_id = pair[1][0]
            y_locus = storf_y_meta[storf_y_meta.find(":")+1:storf_y_meta.find("|")]
            start_y = int(y_locus[:y_locus.find("-")])

            # TODO # this is not correct
            # 16956-17539 compare results for this StORF group
            # not the same! the constraint must be wrong

            # in the case that no storfs in group are selected:
            # need to select StORF of highest value in group
            # this is done at end of ip_filter (where it looks at max of all stORF values in group)

            prob += ip_vars[storf_x_id]*stop_x - ip_vars[storf_y_id]*start_y <= max_overlap
            # prob += ip_vars[storf_x_id]*stop_x - ip_vars[storf_y_id]*start_y >= 0


def is_next_new_group(storf: list, storf_id: int, total_num_storfs: int) -> bool:
    # determine if the StORF after this current StORF is a member of a new group

    # Does the StORF belong to a new collection of overlapping StORFs
    colon_delims, pipe_delims = get_storf_delim(storf)
    # the location of the StORF relative to its overlapping region
    chromosome_rel_loci = (storf[0][colon_delims[1] + 1:pipe_delims[1]])

    # get _x at end of StORF location indicating its position in overlapping group
    overlap_num = int(chromosome_rel_loci[chromosome_rel_loci.find("_")+1:len(chromosome_rel_loci)])


    if storf_id + 1 == total_num_storfs:
        return True
    if overlap_num == 0:
        if storf_id != 0:
            return True
        else: return False
    else: return False
    


def read_fasta() -> list:
    unfiltered_storfs = []
    with open("../../testin/E-coli_output_no_filt.fasta") as storf_file:
        for line in storf_file:
            if line[0] == ">":
                unfiltered_storfs.append([line, next(storf_file)])
    # TMP test files
    # "../../testout/nick_output/E-Coli_storfs_no_filter"
    # "../../testin/smallstorftest.fasta"
    # "../../testin/E-coli_output_no_filt.fasta"
    return unfiltered_storfs


def write_fasta(filtered_storfs: list) -> None:
    f = open("../../testout/oisin_output/output.fasta", "w")
    for storf in filtered_storfs:
        f.write(f"{storf[0]}{storf[1]}")


def ip_set_obj_func(prob: pulp.LpProblem, obj_values: list, ip_vars: list) -> None:
    obj_expression = []
    for storf, value in enumerate(obj_values):
        obj_expression.append((ip_vars[storf], value))
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
    s_total = len(storfs)
    # initialise IP instance
    prob = pulp.LpProblem("StORF_IP_Filter", pulp.LpMaximize)
    # Create IP StORF variables
    storf_ids = [s for s in range(0, s_total)]
    ip_vars = [pulp.LpVariable(f"x_{i}", lowBound=0, upBound=1) for i in range(0, s_total)]
    # create IP y variables (for conditional constraints determined by x variables)
    # ip_vars += [pulp.LpVariable(f"y{i}_", lowBound=0, upBound=1, cat='Integer') for i in range(0, s_total)]
    group = []  # overlapping group of StORFs
    same_group = True
    obj_variables = []


    # pre-compute average gc content
    if HARD_FILTER.gc_range.value is not None: 
        typ = HARD_FILTER.gc_range.value[1]  # 0=mean, 1=median, 2=mode
        ave_gc = get_ave_gc(typ, storfs)
    else:
        ave_gc = None


    # determine coefficient values of IP variables (StORFs)
    values = []
    for storf in storfs:
        values.append(storf_get_value(storf))
    # values += [1 for i in range(0, len(storfs))]  # for y (constraint conditional) variables
    original_values = values

    # initialise objective function
    ip_set_obj_func(prob, values, ip_vars)
    # total_value = sum(x * storf for x,storf in zip(ip_vars, values))
    
    # initialise constraints
    group_tracker = ip_set_constraints(prob, storfs, s_total, ip_vars)

    # Add knapsack sum constraint
    ip_set_total_constraint(prob, ip_vars)


    prob.solve()


    

    solved_variables = sorted(prob.variables(), key=lambda x:int(x.name[2:]))
    filtered_storfs = []
    #21
    s = 0
    for g, group in enumerate(group_tracker):
        selected = False
        for storf in group:
            
            if s == 24:
                print(group)
                for x in group:
                    print(solved_variables[x].name)
                    print(storfs[int(solved_variables[x].name[2:])])
                    print(solved_variables[x].value())
                
            s += 1
            if solved_variables[storf].value() == 1.0:
                filtered_storfs.append(storfs[storf])
                selected = True
        if selected == False:
            # if no permitted overlaps in a group
            group_start = group_tracker[g][0]
            group_end = group_tracker[g][len(group)-1]
            # get highest value in group
            highest_index = original_values.index(max(original_values[group_start:group_end+1]))
            filtered_storfs.append(storfs[highest_index])

    return filtered_storfs


def main():
    storfs = read_fasta()
    storfs = ip_filter(storfs)
    write_fasta(storfs)

if __name__ == '__main__':
    main()