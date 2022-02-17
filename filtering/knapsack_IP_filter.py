from pulp import LpVariable, LpProblem, LpMaximize, value, LpAffineExpression, LpConstraint
import re
import argparse

CONSOLE_ARGS = None


def get_overlaps(uf_value_storfs: list) -> list:
    """
    Assign any StORFs that overlap on the genome into
    a distinct list. Each group of overlapping StORFs
    will be used to create constraints in the integer
    program.
    """
    overlapping_groups = [[]]
    current_stop = None
    group_i = 0
    for storf in uf_value_storfs:
        loci = (storf[0][15:(storf[0].index('|'))])
        stop = int(loci[loci.index('-') + 1:])
        if current_stop is None:
            current_stop = stop
            overlapping_groups[group_i].append(storf)
            continue
        # if overlap exists
        elif stop <= current_stop:
            # !TODO allow argparse to determine range of overlap allowed
            # add to current overlapping group
            overlapping_groups[group_i].append(storf)
        # no overlap
        else:
            group_i += 1
            overlapping_groups.append([storf])
            current_stop = stop
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
    ip_vars = [LpVariable(storf_ids[i], lowBound=0, upBound=1) for i in range(0, len(uf_value_storfs))]
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
    with open("../../testout/E-coli_output_no_filt.fasta") as storf_file:
        for line in storf_file:
            if line[0] == ">":
                unfiltered_storfs.append([line, next(storf_file)])
    return unfiltered_storfs


def get_values(unfiltered_storfs: list):
    """
    Heuristically defined value for each storf within the range 0-1
    """
    # !TODO rewrite lenght attribute already exists with the storf meta data
    largest_storf = max(unfiltered_storfs, key=lambda x: len(x[1]))
    ls_len = len(largest_storf[1])  # largest storf length
    # assess the value of each StORF
    for storf in unfiltered_storfs:
        # !TODO implement additional heuristics for value of StORF

        # This is currently very naive
        # should the program only value according to one
        # how are we to weigh each characteristic?
        # should certain scores take more precedent or are all to
        # be treated equally?

        # value of storf as percentage of ls_len
        length_score = int(len(storf[1]))/ls_len
        # GC percentage
        gc_score = len(re.findall('[GC]', storf[1])) / len(storf[1])

        value = length_score + gc_score
        storf.append(value)

    # GC content, overlap nt size, etc.
    return unfiltered_storfs


def write_fasta(filtered_storfs: list) -> None:
    f = open("../../testout/filtered_StORFs.fasta", "w")
    for storf in filtered_storfs:
        f.write(f"{storf[0]}{storf[1]}")


def get_stats(filtered_storfs: list, uf_storfs: list) -> None:
    """
    Derive statistics from the filtered StORFs such as coverage
    optimality run-time performance, graphs(?) etc.
    # !TODO what other stats could be useful?
    """
    # percentage of StORFs selected (optimality?)
    print(f"{len(filtered_storfs)/len(uf_storfs) * 100} of "
          f"total unfiltered StORFs selected")


def init_cl_args() -> argparse.ArgumentParser:
    args = argparse.ArgumentParser(description='StORF filter options')
    args.add_argument('-mino', '--minoverlap', action="store",
                      dest='min_o', default=0, type=int,
                      help='default=0: Minimum StORF overlap in nt')
    args.add_argument('-maxol', '--maxoverlap', action="store",
                      dest='maxo', default=0, type=int,
                      help='default=0: Maximum StORF overlap in nt')
    args.parse_args()
    global CONSOLE_ARGS
    CONSOLE_ARGS = args


def main():
    init_cl_args()
    uf_storfs = read_fasta()
    uf_value_storfs = get_values(uf_storfs)
    filtered_storfs = filter(uf_value_storfs)
    write_fasta(filtered_storfs)
    get_stats(filtered_storfs, uf_storfs)


if __name__ == '__main__':
    main()
