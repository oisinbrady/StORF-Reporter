import collections  # for ordered dictionary operations
import re  # for StORF sequence identifications (e.g. GC-content)
import pulp  # for IP modelling and solving
import argparse  # for command line argument parsing
import time  # for run time measurements
from math import ceil, floor  # for rounding in "plateau point" calculations

# constants used as filtering bounds
from filter_constants import STORF_CAP_PERCENTAGE, GC_LB, GC_UB, OLAP_LB, SIZE_LB, OLAP_COUNT_UB, ARB_MAX_STORF_SIZE, \
    ARB_MAX_OLAP_SIZE

# set by argparse at run-time
FILTERS = {
    "min_orf": None,
    "max_orf": None,
    "min_olap": None,
    "max_olap": None,
    "max_olap_count": None,
    "min_gc": None,
    "max_gc": None,
    "plateau_value": None,
    "stop_codons": None,
    "mode_stop_codons": None,
    "test_disable_o_count": False,  # for unit testing
}


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


def get_start_stop_bases(storf: list) -> tuple[str, str]:
    start = storf[0].find("Start_Stop=") + len("Start_Stop=")
    start_codon = storf[0][start:start + 3]
    stop = storf[0].find("End_Stop=") + len("End_Stop=")
    stop_codon = storf[0][stop:stop + 3]
    return start_codon, stop_codon


def get_storf_loc(storf: list) -> tuple[int, int]:
    """
    The stop and start location of a StORF relative to the unannotated genome
    :param storf: (list) the StORF
    :return: (tuple(intA, intB)) intA=start location, intB=stop location
    """
    storf_x_meta = storf[0]
    locus = storf_x_meta[storf_x_meta.find(":") + 1:storf_x_meta.find("|")]
    start = int(locus[:locus.find("-")])
    stop = int(locus[locus.find("-") + 1:])
    return start, stop


def is_plateau(percentile_bound: list, initial_m: float, storfs: list) -> bool:
    # assume plateau as a difference in rate of change (m) that is less than theta (0.01)
    theta = 0.01  # threshold
    m = rate_of_change(percentile_bound, storfs)
    return m / initial_m < theta


def rate_of_change(percentile_bound: list, storfs: list) -> float:
    lb = percentile_bound[0]  # bounds of the ECDF percentile group
    ub = percentile_bound[1]
    lb_index = floor(lb * len(storfs))  # get the StORFs at bounds
    ub_index = ceil(ub * len(storfs))
    # return delta y (y-axis: percentile range) / delta x (x-axis: StORF size range)
    delta_x = len(storfs[ub_index][1]) - len(storfs[lb_index][1])
    delta_y = ub - lb
    return delta_y / delta_x


def get_len_ecdf_plateau(unfiltered_storfs: list) -> int:
    # get the StORF length (nt) at the plateau point of an ECDF distribution of StORF lengths
    # none of the selected filtering options will apply to any StORFs greater than this size
    sorted_storfs = sorted(unfiltered_storfs, key=lambda s: len(s[1]))  # order StORFs by len
    percentile_bound = [0.0, 0.5]
    initial_m = rate_of_change(percentile_bound, sorted_storfs)
    i = 0
    lim = 100000  # hard limit in-case no plateau region is detected
    while not is_plateau(percentile_bound, initial_m, sorted_storfs) and i < lim:
        # if difference in rate of change is not enough to be a plateau
        lb = percentile_bound[1]  # lower bound becomes previous upper bound
        ub = (1 - lb) / 2 + lb  # new upper bound definition
        percentile_bound = [lb, ub]  # move to the next percentile slot of StORFs
        i += 1
    # mid-point of percentile bound
    mp = ((percentile_bound[1] - percentile_bound[0]) / 2) + percentile_bound[0]
    mp_index = floor(mp * len(sorted_storfs))
    # determine max StORF length where filtering will be applied
    return len(sorted_storfs[mp_index][1])


def filter_by_overlap(storf_group_values: list, storf_group: list) -> list:
    """
    Assign a value of 0 to all StORFs in group that don't satisfy the overlapping constraints.
    Logic is as follows, for each contig:
        - Return at least one StORF of the highest value
        Exclude:
            - any overlapping StORFs outside defined overlap range
            - fully embedded StORFs (where y is found completely in x start and stop)
        Multiple StORFs can be returned as long as they do not overlap, or are within
        the allowed overlap constraints. If every StORF has an overlap outside the
        allowed range, then return only the highest value StORF.
    :param storf_group_values: (list) current values of each StORF in contig group.
    :param storf_group: (list) a contiguous group of StORFs containing some overlaps.
    :return: (list) An adjusted list of StORF values according to constraints
    """
    if len(storf_group) == 1:
        return storf_group_values
    # modified version of "StORF_finder.py" "tile_filtering" function
    o_min, o_max = FILTERS.get('min_olap'), FILTERS.get('max_olap')
    # Order according to the greatest length first
    sorted_storfs = sorted(storf_group, key=lambda s: len(s[1]), reverse=True)
    # For each StORF, remove all smaller overlapping STORFs according to filtering rules
    length = len(sorted_storfs)
    overlap_count = [0 for i in range(len(storf_group))]
    i = 0  # first StORF of comparison pair
    while i < length:
        start_x, stop_x = get_storf_loc(sorted_storfs[i])
        j = i + 1  # second StORF of comparison pair
        overlaps_x = 0  # number of overlaps of StORF x
        while j < length:
            start_y, stop_y = get_storf_loc(sorted_storfs[j])
            if start_y >= stop_x or stop_y <= start_x:
                # iff no overlap b/w pair
                j += 1
                continue
            elif start_y >= start_x and stop_y <= stop_x:
                # iif StORF j fully embedded in i
                sorted_storfs.pop(j)
                length = len(sorted_storfs)
                overlaps_x += 1
            else:  # +1 needed for stop codon
                overlaps_x += 1
                overlap = stop_x - start_y + 1 if start_x <= start_y else stop_y - start_x + 1
                if not o_min < overlap < o_max:
                    # if j overlap is out of allowed bound
                    sorted_storfs.pop(j)
                    length = len(sorted_storfs)
                else:  # StORF j remains selected
                    j += 1
        overlap_count[i] = overlaps_x
        length = len(sorted_storfs)
        i += 1
    selected = [i[0] for i in sorted_storfs]
    for i, storf in enumerate(storf_group_values):
        if storf[0] not in selected:
            storf[1] = 0  # change value to filter out disallowed StORFs
        # inelegant work-around for unit-testing - b/c unable to mock OPTIONS env variable
        try:
            if not FILTERS.get('test_disable_o_count') or not OPTIONS.disable_olap_count:
                if overlap_count[i] > FILTERS.get("max_olap_count"):
                    storf[1] = 0  # filter by overlap *count* (number of overlaps per StORF)
        except NameError:  # if function called by unit test
            if not FILTERS.get('test_disable_o_count'):
                if overlap_count[i] > FILTERS.get("max_olap_count"):
                    storf[1] = 0
    return storf_group_values


def filter_by_size_range(storf_group_values: list, storf_group: list) -> list:
    # remove StORFs not in the allowed size range
    min_l, max_l = FILTERS.get('min_orf'), FILTERS.get('max_orf')
    for i, storf in enumerate(storf_group):
        if not min_l <= len(re.findall('[AGCT]', storf[1])) <= max_l:
            storf_group_values[i][1] = 0
    return storf_group_values


def filter_by_gc_range(storf_group_values: list, storf_group: list) -> list:
    # remove StORFs not in the allowed GC range
    for i, storf in enumerate(storf_group):
        total_len = len(re.findall('[AGCT]', storf[1]))
        if not FILTERS.get('min_gc') <= len(re.findall('[GC]', storf[1])) / total_len <= FILTERS.get('max_gc'):
            storf_group_values[i][1] = 0
    return storf_group_values


def filter_by_stop_codons(storf_group_values: list, storf_group: list) -> list:
    # remove StORFs that don't start AND end in the allowed stop-codons
    allowed_codons = FILTERS.get('stop_codons')
    for i, storf in enumerate(storf_group):
        sl = len(storf[1])
        start_codon = storf[1][0:3]
        stop_codon = storf[1][sl - 4:sl - 1]
        if start_codon not in allowed_codons or stop_codon not in allowed_codons:
            storf_group_values[i][1] = 0
    return storf_group_values


def filter_favour_length(storf_group_values: list, storf_group: list) -> list:
    # DEPRECATED IP FUNCTION
    """
    Adjust the coefficient value of all StORFs in group according to their length.
    :param storf_group_values: (list) Current coefficient values of StORFs in group
    :param storf_group: (list) A contiguous region of StORFs containing some overlaps
    :return: (list) An adjusted list of StORF values according to their length
    """
    for i, storf in enumerate(storf_group):
        storf_group_values[i][1] *= len(re.findall('[AGCT]', storf[1]))
    return storf_group_values


def filter_favour_most_gc(storf_group_values: list, storf_group: list) -> list:
    # DEPRECATED IP FUNCTION
    """
    Adjust the coefficient value of all StORFs in group according to their gc-content.
    :param storf_group_values: (list) Current coefficient values of StORFs in group.
    :param storf_group: (list) A contiguous region of StORFs containing some overlaps.
    :return: (list) An adjusted list of StORF values according to their gc-content.
    """
    for i, storf in enumerate(storf_group):
        storf_group_values[i][1] *= len(re.findall('[GC]', storf[1]))
    return storf_group_values


def filter_by_mode_stop_codons(storf_group_values: list, storf_group: list, mode_codon: list) -> list:
    for i, storf in enumerate(storf_group):
        start = storf[1][0].find("Start_Stop=") + len("Start_Stop=")
        start_codon = storf[1][0][start:start + 3]
        stop = storf[1][0].find("End_Stop=") + len("End_Stop=")
        stop_codon = storf[1][0][stop:stop + 3]
        codon_id = (start_codon, stop_codon)
        if codon_id != mode_codon:
            storf_group_values[i][1] = 0
    return storf_group_values


def filter_favour_gc_by_len(storf_group_values: list, storf_group: list) -> list:
    # DEPRECATED IP FUNCTION
    for i, storf in enumerate(storf_group):
        storf_len = len(re.findall('[AGCT]', storf[1]))
        storf_gc = len(re.findall('[GC]', storf[1]))
        # based of cumulative PCA gc vs size
        gc_by_size = storf_gc / storf_len * storf_len
        storf_group_values[i][1] = gc_by_size
    return storf_group_values


def ip_set_weighted_values(storf_contig_values: list, contig_group: list) -> list:
    # DEPRECATE IP FUNCTION
    # WEIGHTED VALUE FILTERS: coefficient values equivalent to attributes of StORF, e.g. length
    # Only makes a difference if knapsack capacity constraints are a factor
    if OPTIONS.len_weighted_value:
        storf_contig_values = filter_favour_length(storf_contig_values, contig_group)
    if OPTIONS.gc_weighted_value:
        storf_contig_values = filter_favour_most_gc(storf_contig_values, contig_group)
    if OPTIONS.len_gc_weighted_value:
        storf_contig_values = filter_favour_gc_by_len(storf_contig_values, contig_group)
    return storf_contig_values


def filter_by_bounds(storf_contig_values: list, contig_group: list) -> list:
    # value is 0,1 iff isn't or is in the following bounds
    if not OPTIONS.disable_olap:
        storf_contig_values = filter_by_overlap(storf_contig_values, contig_group)
    if not OPTIONS.disable_size_filter:
        storf_contig_values = filter_by_size_range(storf_contig_values, contig_group)
    if not OPTIONS.disable_gc:
        if not OPTIONS.original_filter:
            storf_contig_values = filter_by_gc_range(storf_contig_values, contig_group)
    if FILTERS.get('stop_codons') is not None:
        storf_contig_values = filter_by_stop_codons(storf_contig_values, contig_group)
    return storf_contig_values


def set_group_values(contig_group: list) -> list:
    # DEPRECATED IP FUNCTION
    """
    Adjust the coefficient value of all StORFs in group according to the user defined filters.
    :param contig_group: (list) A contiguous region of StORFs containing some overlap
    :return: (list) An adjusted list of StORF values according to the selected filters.
    """
    # Default value of StORF is 1, if StORF doesn't obey restrictions then overwrite to 0
    storf_contig_values = [[s[0], 1] for s in contig_group]  # storf id and value
    storf_contig_values = filter_by_bounds(storf_contig_values, contig_group)
    storf_contig_values = ip_set_weighted_values(storf_contig_values, contig_group)
    # apply less strict filtering for sub-set S
    # S has a much higher distribution of HSS (High Significance StORFs)
    if not OPTIONS.disable_ecdf_relaxation:
        if not OPTIONS.original_filter:
            p = FILTERS.get("plateau_value")
            # overwrite filter result for S
            for i in range(0, len(storf_contig_values)):
                if len(contig_group[i][1]) >= p:
                    storf_contig_values[i][1] = 1
    return storf_contig_values


def is_next_new_group(storfs: list, storf_id: int) -> bool:
    if storf_id == 0:
        return False
    colon_delim, pipe_delim = get_storf_delim(storfs[storf_id])
    # StORF location relative to contig region
    chromosome_rel_loci = (storfs[storf_id][0][colon_delim[1] + 1:pipe_delim[1]])
    # get _x at end of StORF location indicating its position in overlapping group
    overlap_num = int(chromosome_rel_loci[chromosome_rel_loci.find("_") + 1:len(chromosome_rel_loci)])
    return overlap_num == 0


def read_fasta(file_name=None) -> list:
    """
    Convert unannotated genome ".fasta" file to a list of StORFs.
    :return: (list) A list of all StORFs in fasta file input
    """
    unfiltered_storfs = []
    if file_name == None:
        file = OPTIONS.fasta
    else:
        file = file_name
    with open(file) as storf_file:
        for line in storf_file:
            if line[0] == ">":
                unfiltered_storfs.append([line, next(storf_file)])
    return unfiltered_storfs


def write_fasta(filtered_storfs: list) -> None:
    """
    Convert the final filtered subset of StORFs into a fasta file
    :param filtered_storfs: (list) StORFs from optimal knapsack solution
    :return: (None)
    """
    f = open(f"{OPTIONS.output}", "w")
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
    obj_i = 0
    for _, value in obj_values:
        obj_expression.append((ip_vars[obj_i], value))
        obj_i += 1
    e = pulp.LpAffineExpression(obj_expression)
    prob += e  # add objective function (in turn, Pulp adds variable bounds = [0,1])


def ip_set_total_constraint(prob: pulp.LpProblem, ip_vars: list) -> None:
    """
    Set the knapsack capacity of selected StORFs in unannotated genome
    :param prob: (pulp.LpProblem) The integer program
    :param ip_vars: (list of pulp.LpVariable) the IP variables
    :return: (None)
    """
    # constrain the sum of all knapsack weights 
    # I.e, of all StORFs in file to be selected
    if OPTIONS.total_cap == -1:
        # average HSS percentage is ~7.4% per test genomes
        tc = round(len(ip_vars) * STORF_CAP_PERCENTAGE)
    else:
        tc = OPTIONS.total_cap  # total selected StORF constraint
    g = [(i, 1) for i in ip_vars]
    e = pulp.LpAffineExpression(g)
    # RHS: <= C_t
    c = pulp.LpConstraint(e, -1, f"external_knapsack_constraint", tc)
    prob += c


def propability_distribution(x_sorted_storfs: list, percentiles: list) -> list:
    # I.E. box and whiskers plot for StORF attribute ('x')
    total_storfs = len(x_sorted_storfs)
    distributions = []
    for percentile in percentiles:
        distributions.append(x_sorted_storfs[round(total_storfs * percentile)])
    return distributions


def get_con_groups(storfs: list) -> list:
    con_groups = []
    con_group = []
    for i, storf in enumerate(storfs):
        if is_next_new_group(storfs, i):
            con_groups.append(con_group)
            con_group = [storf]
        else:
            con_group.append(storf)
    if len(con_group) != 0:
        con_groups.append(con_group)
    return con_groups


def set_overlap_bounds(storfs) -> None:
    if FILTERS.get("min_olap") is None:
        FILTERS["min_olap"] = OLAP_LB
    if FILTERS.get("max_olap") is None:
        FILTERS["max_olap"] = ARB_MAX_OLAP_SIZE
    if FILTERS.get("max_olap") < FILTERS.get("min_olap"):
        FILTERS["min_olap"] = 0


def set_orf_bounds(storfs) -> None:
    if FILTERS.get("min_orf") is None:
        # size based of cumulative distribution PCA size vs gc-content scatter-plot
        FILTERS['min_orf'] = SIZE_LB
    if FILTERS.get("max_orf") is None:
        FILTERS['max_orf'] = ARB_MAX_STORF_SIZE


def set_gc_bounds(storfs) -> None:
    if FILTERS.get('min_gc') is None:
        FILTERS['min_gc'] = GC_LB
    if FILTERS.get('max_gc') is None:
        FILTERS['max_gc'] = GC_UB


def def_user_filter_args(storfs):
    if not OPTIONS.disable_size_filter:
        FILTERS["min_orf"] = OPTIONS.min_orf
        FILTERS["max_orf"] = OPTIONS.max_orf
        if OPTIONS.min_orf is None or OPTIONS.max_orf is None:
            set_orf_bounds(storfs)
    if not OPTIONS.disable_olap:
        FILTERS["min_olap"] = OPTIONS.min_olap
        FILTERS["max_olap"] = OPTIONS.max_olap
        if OPTIONS.min_olap is None or OPTIONS.max_olap is None:
            set_overlap_bounds(storfs)
    if not OPTIONS.disable_olap_count:
        FILTERS["max_olap_count"] = OPTIONS.max_olap_count
        if OPTIONS.max_olap_count is None:
            FILTERS["max_olap_count"] = OLAP_COUNT_UB
    if not OPTIONS.disable_gc:
        FILTERS["min_gc"] = OPTIONS.min_gc
        FILTERS["max_gc"] = OPTIONS.max_gc
        if OPTIONS.min_gc is None or OPTIONS.max_gc is None:
            set_gc_bounds(storfs)
    if OPTIONS.stop_codons is not None and OPTIONS.stop_codons is not ["TAG", "TGA", "TAA"]:
        FILTERS['stop_codons'] = OPTIONS.stop_codons


def print_all_filter_params(storfs: list) -> None:
    print(FILTERS)
    if not OPTIONS.basic_filter:
        if OPTIONS.total_cap == -1:
            print(f"Total weight capacity: {round(len(storfs) * STORF_CAP_PERCENTAGE)}")
        else:
            print(f"Total weight capacity: {OPTIONS.total_cap}")
    exit()


def set_filters(storfs) -> None:
    if OPTIONS.original_filter:
        # original filter parameters of StORF_finder.py
        FILTERS['min_olap'] = 0
        FILTERS['max_olap'] = 50
        FILTERS['min_orf'] = 100
        FILTERS['max_orf'] = 50000
        OPTIONS.disable_ecdf_relaxation = True
        OPTIONS.disable_olap_count = True
        OPTIONS.total_cap = len(storfs)
        if OPTIONS.print_filter_params:
            print_all_filter_params(storfs)
        return None
    def_user_filter_args(storfs)
    FILTERS['plateau_value'] = get_len_ecdf_plateau(storfs) if not OPTIONS.disable_ecdf_relaxation else None
    if OPTIONS.print_filter_params:
        print_all_filter_params(storfs)


def ip_process_output(prob: pulp.LpProblem, storfs: list) -> list:
    selected = {}
    for var in prob.variables():  # get selected StORFs from IP solution
        selected[var.name] = pulp.value(var)
    ordered_selected = collections.OrderedDict(sorted(selected.items(), key=lambda t: int(t[0][2:])))
    output = []
    for i, v in enumerate(ordered_selected.values()):
        if v >= 1.0:
            output.append(storfs[i])
    return output


def basic_filter(storfs: list) -> list:
    selected = []  # 1=selected, 0=not-selected
    contigs = get_con_groups(storfs)
    for g, contig in enumerate(contigs):
        contig_values = [[s[0], 1] for s in contig]  # default all StORFs to selected
        selected += filter_by_bounds(contig_values, contig)  # filter according to bounds
    filtered_storfs = []
    # prepare output
    for i, storf_data in enumerate(selected):
        is_selected = storf_data[1]
        storf_fasta = storfs[i]  # StORF meta data and sequence
        if is_selected:
            filtered_storfs.append(storf_fasta)
    return filtered_storfs


def ip_filter(storfs: list) -> list:
    # DEPRECATED IP FUNCTION
    """
    Filter StORFs in unannotated genome file using integer programming.
    Knapsack problem:
    - Weight of each StORF = 1
    - The following default to infinity:
        - gc = capacity on number of selected StORFs in each group (sub-knapsack)
        - tc = capacity of total selected StORFs
    - All other constraints are selected by the user and so the model is defined accordingly
    :param storfs: (list) all StORFs in unannotated genome
    :return: (list) optimal IP solution of selected StORFs
    """
    # initialise IP instance
    prob = pulp.LpProblem("StORF_IP_Filter", pulp.LpMaximize)
    s_total = len(storfs)
    # Create IP variables
    storf_ids = [f"x_{s}" for s in range(0, s_total)]
    ip_vars = [pulp.LpVariable(storf_ids[i], lowBound=0, upBound=1, cat='Integer') for i in range(0, s_total)]
    obj_values = []
    # value-setting of IP vars to filter unwanted StORFs
    contigs = get_con_groups(storfs)
    for g, contig in enumerate(contigs):
        obj_values += set_group_values(contig)
    # construct objective function (auto-adds IP variable bounds)
    ip_set_obj_func(prob, obj_values, ip_vars)
    ip_set_total_constraint(prob, ip_vars)  # add knapsack sum capacity constraint
    prob.solve()  # solve the knapsack problem
    return ip_process_output(prob, storfs)


def init_argparse():
    parser = argparse.ArgumentParser(description='StORF filtering parameters')
    parser.add_argument('-f', action="store", dest='fasta', required=True,
                        help='Input FASTA File')
    parser.add_argument('-o', action="store", dest='output', required=True,
                        help='Output name as FASTA File')
    parser.add_argument('-min_orf', action="store", dest='min_orf', type=int,
                        help='Minimum StORF size (nt)')
    parser.add_argument('-max_orf', action="store", dest='max_orf', type=int,
                        help='Maximum StORF size (nt)')
    parser.add_argument('-min_olap', action="store", dest='min_olap', type=int,
                        help='Maximum StORF overlap size (nt)')
    parser.add_argument('-max_olap', action="store", dest='max_olap', type=int,
                        help='Maximum StORF overlap size (nt)')
    parser.add_argument('-max_olap_count', action="store", dest='max_olap_count', type=int,
                        help='Maximum number of overlaps per StORF')
    parser.add_argument('-min_gc', action="store", dest='min_gc', type=float,
                        help='Minimum gc percentage of StORFs as float')
    parser.add_argument('-max_gc', action="store", dest='max_gc', type=float,
                        help='Maximum gc percentage of StORFs as float')
    parser.add_argument('-codons', action="store", dest='stop_codons', default="TAG,TGA,TAA",
                        help='Default - (\'TAG,TGA,TAA\'): List Stop Codons to use')
    # filtering using IP method: DEPRECATED
    parser.add_argument('-ip', action="store_true", dest='kpip_filter',
                        help='filter by basic bounds for StORF attributes')
    parser.add_argument('-ip_v_l', action="store_true", dest='len_weighted_value',
                        help='Let IP favour StORF selection by size')
    parser.add_argument('-ip_v_gc', action="store_true", dest='gc_weighted_value',
                        help='Let IP favour StORF selection by gc content'),
    parser.add_argument('-ip_v_gcl', action="store_true", dest='len_gc_weighted_value',
                        help='Let IP favour StORF length * GC percentage'),
    # ip weight capacity
    parser.add_argument('-cap', action="store", dest='total_cap', type=int, default=-1,
                        help='Default - infinity/No capacity: set IP maximum StORF selection capacity')
    # mimic original StORF filtering ("StORF-Reporter")
    parser.add_argument('-original_filter', action="store_true", dest='original_filter',
                        help='Set all filtering parameters equal to StORF_finder.py')
    # the following commands allow certain filtering methods to be disabled at run-time
    parser.add_argument('-d_gc', action="store_true", dest='disable_gc',
                        help='Disable filtering by GC content')
    parser.add_argument('-d_olap', action="store_true", dest='disable_olap',
                        help='Disable filtering by overlap size (nt)')
    parser.add_argument('-d_olap_count', action="store_true", dest='disable_olap_count',
                        help='')
    parser.add_argument('-d_size', action="store_true", dest='disable_size_filter',
                        help='Disable filtering by StORF size (nt)'),
    parser.add_argument('-d_ecdf', action="store_true", dest='disable_ecdf_relaxation',
                        help='Disable relaxation of filtering for the largest StORFs (using ECDF plateau region)'),
    # Utility/Debug arguments
    parser.add_argument('-print_params', action="store_true", dest='print_filter_params',
                        help='Print all filters after pre-calculations')
    parser.add_argument('-run_time', action="store_true", dest='run_time',
                        help='Print total program execution time')

    # TODO ADD post filter "twin" StORF labelling with blastx
    global OPTIONS
    OPTIONS = parser.parse_args()


def main():
    """
    Filter the unannotated genome according to the user's selections
    """
    init_argparse()
    storfs = read_fasta()
    set_filters(storfs)
    if not OPTIONS.kpip_filter:
        storfs = basic_filter(storfs)
    else:
        storfs = ip_filter(storfs)
    write_fasta(storfs)


if __name__ == '__main__':
    start_time = time.time()
    main()
    if OPTIONS.run_time:
        print(f"{time.time() - start_time}")
