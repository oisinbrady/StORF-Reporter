from enum import Enum

import pulp
import re
import collections

""" Filters to match default behaviour of StORF_finder.py program is:
        - storf_length=True
        - overlap_range=[0, 50]
        - size_range=[100, 50000] 
"""


class HardFilter(Enum):
    """
    Filters that bound the sub-set of final selected StORFs
    """
    overlap_range = [0, 50]  # default [0, 50]  #67, 1000
    size_range = [100, 50000]  # default [100, 50000]
    gc_range = None  # [1000, 0]  # 10=percentage variance, 0=mean, 1=median, 2=mode


class SoftFilter(Enum):
    """
    Directly affect values of each StORF in objective function based of selected booleans.
    storf_length: favour larger StORFs
    gc_length: favour higher gc count in StORFs
    N.B.:
        IF |"overlapping" group| > WeightConstraints.storf_group:
            THEN soft filters will decide which StORFs are taken as to maximise the objective function
    """
    storf_length = True  # True=favour the largest StORFs (more value) in groups
    gc_length = None  # True=favour the StORFs with highest gc% (more value) in groups


class WeightConstraints(Enum):
    """
    The knapsack constraints.
    sum_total: total number of StORFs allowed
    storf_group: total number of StORFs allowed in each group
    A value of -1 indicates that there is no restriction. I.e., constraints don't apply; get as many as possible
    """
    sum_total = -1
    storf_group = -1


def get_storf_delim(storf: list) -> tuple[list[int], list[int]]:
    """
    The indexes of ':' and '|' characters in a StORF.
    Useful for finding specific StORF information E.g., locus
    :param storf: (list) the StORF
    :return: (tuple(listA, listB)) listA=colon indexes, ListB=pipe indexes
    """
    colon_delimiters = [s.start() for s in re.finditer(r":", storf[0])]
    pipe_delimiters = [s.start() for s in re.finditer(r"\|", storf[0])]
    return colon_delimiters, pipe_delimiters


def get_ave_gc(average_type: int, storfs: list) -> float:
    """
    The average percentage of G and C bases of all StORFs in unannotated genome
    :param average_type: (int) 0=mean,1=median,2=mode
    :param storfs: (list) all StORFs in unannotated genome
    :return: (float) the average gc percentage of all StORFs
    """
    if average_type == 0:  # mean
        total_gc = 0
        num_storfs = len(storfs)
        for storf in storfs:
            total_gc += len(re.findall('[GC]', storf[1]))
        return float(total_gc / num_storfs)
    elif average_type == 1:  # median
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
    else:  # mode
        # get gc count of each StORF
        s_gc_len = [len(re.findall('[GC]', s[1])) for s in storfs]
        # get the most occurring gc count
        mode_largest = max(set(s_gc_len), key=s_gc_len.count)
        # TODO currently only selects the largest mode if multiple exist
        # add some functionality for determining which mode to take?
        return float(mode_largest)


def get_storf_loc(storf: list) -> tuple[int, int]:
    """
    The stop and start location of a StORF relative to the unannotated genome
    :param storf: (list) the StORF
    :return: (tuple(intA, intB)) intA=start location, intB=stop location
    """
    storf_x_meta = storf[1][0]
    locus = storf_x_meta[storf_x_meta.find(":") + 1:storf_x_meta.find("|")]
    start = int(locus[:locus.find("-")])
    stop = int(locus[locus.find("-") + 1:])
    return start, stop


def filter_by_overlap(storf_group_values: list, storf_group: list) -> list:
    """
    Assign a coefficient value of 0 to all StORFs in group that don't satisfy the overlapping constraints.
    I.e., Get the highest valued StORF(s) obeying knapsack constraints.
    Logic is as follows:
        - Return at least one StORF of the highest value
        Exclude:
            - any overlapping StORFs outside defined overlap range
            - fully embedded StORFs (where y is found completely in x start and stop)
        Multiple StORFs can be returned as long as they do not overlap, or are within
        the allowed overlap constraints. If every StORF has an overlap outside the
        allowed range, then return only the highest value StORF.
    :param storf_group_values: (list) current values of each StORF.
    :param storf_group: (list) a contiguous group of StORFs containing some overlaps.
    :return: (list) An adjusted list of StORF values according to constraints
    """
    # modified version of "StORF_finder.py" "tile_filtering" function
    o_min, o_max = HardFilter.overlap_range.value[0], HardFilter.overlap_range.value[1]
    # Order according to the greatest length first
    ordered_by_value = []
    sorted_values = sorted(storf_group_values, key=lambda k: k[1], reverse=True)
    offset = storf_group_values[0][0]
    for s in sorted_values:
        ordered_by_value.append(storf_group[s[0] - offset])
    # For each StORF, remove all smaller overlapping STORFs according to filtering rules
    length = len(ordered_by_value)
    i = 0
    while i < length:
        x_index = ordered_by_value[i][0] - offset
        start_x, stop_x = get_storf_loc(storf_group[x_index])
        j = i + 1
        while j < length:
            y_index = ordered_by_value[j][0] - offset
            start_y, stop_y = get_storf_loc(storf_group[y_index])
            if start_y >= stop_x or stop_y <= start_x:
                j += 1
                continue  # Not caught up yet / too far
            elif start_y >= start_x and stop_y <= stop_x:
                ordered_by_value.pop(j)
                length = len(ordered_by_value)
            else:  # +1 needed for stop codon
                x = set(range(start_x, stop_x + 1))
                y = set(range(start_y, stop_y + 1))
                overlap = len(x.intersection(y))
                if o_min <= overlap >= o_max:
                    ordered_by_value.pop(j)
                    length = len(ordered_by_value)
                else:
                    j += 1
        length = len(ordered_by_value)
        i += 1
    selected = [i[0] for i in ordered_by_value]
    for storf in storf_group_values:
        if storf[0] not in selected:
            storf[1] = 0
    return storf_group_values


def filter_by_size_range(storf_group_values: list, storf_group: list) -> list:
    """
    Assign a coefficient value of 0 to all StORFs in group that don't satisfy the size constraint.
    I.e., filter out StORFs that are too large/small
    :param storf_group_values: (list) current coefficient values of StORFs in group
    :param storf_group: (list) a contiguous region of StORFs containing some overlaps
    :return: (list) An adjusted list of StORF values according to size constraint
    """
    for i, storf in enumerate(storf_group):
        min_l, max_l = HardFilter.size_range.value[0], HardFilter.size_range.value[1]
        if not min_l <= len(re.findall('[AGCT]', storf[1][1])) <= max_l:
            storf_group_values[i][1] = 0
    return storf_group_values


def filter_by_gc_range(storf_group_values: list, storf_group: list, ave_gc: float) -> list:
    """
    Assign a coefficient value of 0 to all StORFs in group that don't satisfy the gc constraint.
    I.e., filter out StORFs that are not within a percentage deviation of the average gc
    :param storf_group_values: (list) Current coefficient values of StORFs in group
    :param storf_group: (list) A contiguous region of StORFs containing some overlaps
    :param ave_gc: (float) The average gc percentage of all StORFs in unannotated genome (see get_ave_gc())
    :return: (list) An adjusted list of StORF values according to gc constraint
    """
    min_gc = ave_gc - (ave_gc * 0.05)
    max_gc = ave_gc + (ave_gc * 0.05)
    for i, storf in enumerate(storf_group):
        if not min_gc <= len(re.findall('[GC]', storf[1][1])) <= max_gc:
            storf_group_values[i][1] = 0
    return storf_group_values


def filter_favour_length(storf_group_values: list, storf_group: list) -> list:
    """
    Adjust the coefficient value of all StORFs in group according to their length.
    :param storf_group_values: (list) Current coefficient values of StORFs in group
    :param storf_group: (list) A contiguous region of StORFs containing some overlaps
    :return: (list) An adjusted list of StORF values according to their length
    """
    for i, storf in enumerate(storf_group):
        # must be multiplication so that hard filters still affect coefficient values
        storf_group_values[i][1] *= len(re.findall('[AGCT]', storf[1][1]))
    return storf_group_values


def filter_favour_most_gc(storf_group_values: list, storf_group: list) -> list:
    """
    Adjust the coefficient value of all StORFs in group according to their gc-content.
    :param storf_group_values: (list) Current coefficient values of StORFs in group.
    :param storf_group: (list) A contiguous region of StORFs containing some overlaps.
    :return: (list) An adjusted list of StORF values according to their gc-content.
    """
    for i, storf in enumerate(storf_group):
        storf_group_values[i][1] *= len(re.findall('[GC]', storf[1][1]))
    return storf_group_values


def set_group_values(storf_group: list, ave_gc: None | int) -> list:
    """
    Adjust the coefficient value of all StORFs in group according to the user defined filters.
    :param storf_group: (list) A contiguous region of StORFs containing some overlap
    :param ave_gc: (float) The average gc percentage of each StORF
    :return: (list) An adjusted list of StORF values according to the selected filters.
    """
    storf_group_values = [[s[0], 1] for s in storf_group]  # storf id and value
    # Default value of StORF is 1, if StORF doesn't obey restrictions then overwrite to 0
    """
    SOFT FILTERS, final subset effected by less defined StORF characteristics
    # e.g. FAVOUR larger sized StORFs(soft), only StORFs in range x-y (HARD)
    """
    if SoftFilter.storf_length.value is not None:
        storf_group_values = filter_favour_length(storf_group_values, storf_group)
    if SoftFilter.gc_length.value is not None:
        storf_group_values = filter_favour_most_gc(storf_group_values, storf_group)
    # HARD FILTERS, final subset discretely defined using ranges
    if HardFilter.overlap_range.value is not None:
        storf_group_values = filter_by_overlap(storf_group_values, storf_group)
    if HardFilter.size_range.value is not None:
        storf_group_values = filter_by_size_range(storf_group_values, storf_group)
    if HardFilter.gc_range.value is not None:
        storf_group_values = filter_by_gc_range(storf_group_values, storf_group, ave_gc)
    # TODO HardFilter by stop codon type
    return storf_group_values


def is_next_new_group(storf: list, storf_id: int, total_num_storfs: int) -> bool:
    """
    Determine if the current StORF belongs to a new group from the previous stream of StORFs.
    A group is defined as: StORFs within a contiguous region of the unannotated genome containing some overlap
    :param storf: (list) The StORF
    :param storf_id: (int) The StORF's id relative to all StORFs
    :param total_num_storfs: (int) The total number of StORFs in unannotated genome
    :return: (bool) True if StORF belongs to new group, else False
    """
    colon_delim, pipe_delim = get_storf_delim(storf)
    # the location of the StORF relative to its overlapping region
    chromosome_rel_loci = (storf[0][colon_delim[1] + 1:pipe_delim[1]])
    # get _x at end of StORF location indicating its position in overlapping group
    overlap_num = int(chromosome_rel_loci[chromosome_rel_loci.find("_") + 1:len(chromosome_rel_loci)])
    if storf_id + 1 == total_num_storfs:
        return True
    if overlap_num == 0:
        if storf_id != 0:
            return True
        else:
            return False
    else:
        return False


def read_fasta() -> list:
    """
    Convert unannotated genome ".fasta" file to a list of StORFs.
    :return: (list) A list of all StORFs in fasta file input
    """
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
    """
    Convert the final filtered subset of StORFs into a fasta file
    :param filtered_storfs: (list) StORFs from optimal knapsack solution
    :return: (None)
    """
    f = open("../../testout/oisin_output/output.fasta", "w")
    for storf in filtered_storfs:
        f.write(f"{storf[0]}{storf[1]}")


def ip_set_obj_func(prob: pulp.LpProblem, obj_values: list, ip_vars: list) -> None:
    """
    Initialise the objective function of the integer program.
    :param prob: (pulp.LpProblem) The integer program
    :param obj_values: (list) The IP variable (StORF) reference and its value
    :param ip_vars: (list of pulp.LpVariable) the IP variables
    :return: (None)
    """
    obj_expression = []
    for obj_id, value in obj_values:
        obj_expression.append((ip_vars[obj_id], value))
    e = pulp.LpAffineExpression(obj_expression)
    prob += e  # add objective function (in turn, Pulp adds variable bounds)


def ip_set_group_constraint(prob: pulp.LpProblem, ip_vars: list, group: list, group_id: int) -> None:
    """
    Set the knapsack GROUP capacity of selected StORFs in the unannotated genome
    :param prob: (pulp.LpProblem) The integer program
    :param ip_vars: (list of pulp.LpVariable) the IP variables
    :param group: (list) A contiguous region of StORFs containing some overlap
    :param group_id: (int) The current group's ID
    :return: (None)
    """
    if WeightConstraints.storf_group.value == -1:
        # iff infinity (i.e. AMAP StORFs in each group), no need for constraint
        return
    else:
        # constrain group, where weight of each var = 1
        gc = WeightConstraints.storf_group.value  # group constraint
        g = []
        for i, storf in enumerate(group):
            g.append((ip_vars[storf[0]], 1))
        e = pulp.LpAffineExpression(g)
        # LHS: sum of selected StORFs in group, RHS: <= gc
        c = pulp.LpConstraint(e, -1, f"internal_knapsack_constraint_{group_id}", gc)
        prob += c


def ip_set_total_constraint(prob: pulp.LpProblem, ip_vars: list) -> None:
    """
    Set the knapsack capacity of selected StORFs in unannotated genome
    :param prob: (pulp.LpProblem) The integer program
    :param ip_vars: (list of pulp.LpVariable) the IP variables
    :return: (None)
    """
    # constrain the sum of all knapsack weights 
    # I.e, of all StORFs in file to be selected
    if WeightConstraints.sum_total.value == -1:
        return
    else:
        tc = WeightConstraints.sum_total.value  # total selected StORF constraint
        g = [(i, 1) for i in ip_vars]
        e = pulp.LpAffineExpression(g)
        # RHS: <= C_t
        c = pulp.LpConstraint(e, -1, f"external_knapsack_constraint", tc)
        prob += c


def ip_filter(storfs: list) -> list:
    """
    Filter StORFs in unannotated genome file using integer programming
    :param storfs: (list) all StORFs in unannotated genome
    :return: (list) optimal IP solution of selected StORFs
    """
    """
    Knapsack problem: 
    - Weight of each StORF = 1
    - The following default to infinity:
        - gc = capacity on number of selected StORFs in each group (sub-knapsack)
        - tc = capacity of total selected StORFs
    - All other constraints are selected by the user and so the model is defined accordingly
    """
    # initialise IP instance
    prob = pulp.LpProblem("StORF_IP_Filter", pulp.LpMaximize)
    s_total = len(storfs)
    g = 0  # group id
    s = 0  # storf id
    # Create IP variables
    storf_ids = [f"x_{s}" for s in range(0, s_total)]
    ip_vars = [pulp.LpVariable(storf_ids[i], lowBound=0, upBound=1, cat='Integer') for i in range(0, s_total)]
    group = []  # overlapping group of StORFs
    obj_values = []
    # pre-compute average gc content
    if HardFilter.gc_range.value is not None:
        typ = HardFilter.gc_range.value[1]  # 0=mean, 1=median, 2=mode
        ave_gc = get_ave_gc(typ, storfs)
    else:
        ave_gc = None
    # determine coefficient values of future IP variables
    for storf in storfs:
        # StORF values are dependent on their group
        if is_next_new_group(storf, s, s_total):
            # create list of to be objective variables with coefficients for each StORF
            obj_values += set_group_values(group, ave_gc)
            # add weight constraint to each group
            ip_set_group_constraint(prob, ip_vars, group, g)
            group = [(s, storf)]  # init new group
            g += 1
        else:
            group.append((s, storf))
        s += 1
    # add values to last group of StORFs...  
    obj_values += set_group_values(group, ave_gc)
    # set constraint for last group (auto-adds IP variable bounds)
    ip_set_group_constraint(prob, ip_vars, group, g)
    # construct objective function (auto-adds IP variable bounds)
    ip_set_obj_func(prob, obj_values, ip_vars)
    # Add knapsack sum constraint
    ip_set_total_constraint(prob, ip_vars)
    # get selected StORFs from IP solution
    prob.solve()
    selected = {}
    for var in prob.variables():
        # print(f"{var.name}={pulp.value(var)}")
        selected[var.name] = pulp.value(var)
    ordered_selected = collections.OrderedDict(sorted(selected.items(), key=lambda t: int(t[0][2:])))
    final_filter = []
    for i, v in enumerate(ordered_selected.values()):
        if v == 1.0:
            final_filter.append(storfs[i])
    return final_filter


def main():
    """
    Filter the unannotated genome according to the user's selections
    :return: (None)
    """
    storfs = read_fasta()
    storfs = ip_filter(storfs)
    write_fasta(storfs)


if __name__ == '__main__':
    main()
