from pulp import LpVariable, LpProblem, LpMaximize, value, LpAffineExpression, LpConstraint, LpStatus

# read StORF files - all_StORFs Data type dict of lists?
    # each storf holds:
        # its meta-data (e.g. start-stop, UR id, etc.)
        # selected {0,1} - dictates if the StORF has been selected via the IP
        # value - the heuristic value of the StORF
            # represent as a normalised value 0-1 e.g. 0.95 is a very 'good' StORF

# initialise IP instance

# establish constraints
    # for each COMBINATION from overlapping StORFs group, e.g. x_1, x_2, x_3: (find_overlaps(storfs) -> overlapping_storfs)

    # For all storfs in a combination group: sum.selected = {0,1}
        # x_1.selected + x_2.selected + x_3.selected >=0
        # x_1.selected + x_2.selected + x_3.selected <=1

        # i.e. ONLY 1 storf of the combination group can be selected

# establish bound(s)
    # x.selected C- {0,1} x = 1,2,...n
        # each StORF is either selected or not for filtering

# establish objective function
    # MAXIMISE(SIGMA(S)), where S = x C- all_StORFs
    # e.g. MAX ({0,1}x_1.value + {0,1}x_2.value + {0,1}x_3.value)


def get_overlaps(uf_value_storfs: list) -> list:
    '''
    assign any storfs that overlap on the genome into
    a distinct list in overlapping_groups
    '''

    # all storfs in overlapping_groups will be variables
    # in the IP and these groups will be constraints
    # (0 <= sum <= 1) in the IP
    # the objective function will be (max of) the sum of
    # all variables
    overlapping_groups = [[]]
    current_stop = None
    group_i = 0
    for storf in uf_value_storfs:
        loci = (storf[0][15:(storf[0].index('|'))])
        start = int(loci[:loci.index('-')])
        stop = int(loci[loci.index('-') + 1:])

        if current_stop is None:
            current_stop = stop
            overlapping_groups[group_i].append(storf)
            continue

        # if overlap exists
        elif stop <= current_stop:
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
    exclude overlapping StORFs
    '''
    # convert to list for fast storf value lookups
    storf_dict = {i: uf_value_storfs[i] for i in range(0, len(uf_value_storfs))}
    # initialise lp instance
    prob = LpProblem("StORF_Knapsack_Problem", LpMaximize)
    # initialise lp variables
    storf_ids = [f"x_{i}" for i in range(0, len(uf_value_storfs))]
    ip_vars = [LpVariable(storf_ids[i], lowBound=0, upBound=1) for i in range(0, len(uf_value_storfs))]
    # get overlapping StORFs
    overlapping_storf_groups = get_overlaps(uf_value_storfs)
    # for each group of overlapping storfs, create a constraint
    ip_var_id = 0
    obj_expression = []
    for group_i, storf_group in enumerate(overlapping_storf_groups):
        c_g = []  # constraint group
        for storf_i, storf in enumerate(storf_group):
            # storf id and its heuristic value
            c_g.append((ip_vars[ip_var_id], storf[2]))
            # get all distinct variable expressions for the objective function
            if storf_i == 1:
                obj_expression.append((ip_vars[ip_var_id], storf[2]))
            ip_var_id += 1
        # expression construction
        e = LpAffineExpression(c_g)
        # initialise constraint e <= 1
        c = LpConstraint(e, -1, f"StORF_group_constraint_{group_i}", 1)
        # assign constraint to IP
        prob += c

    # initialise objective function
    prob += LpAffineExpression(obj_expression)
    # print(prob.objective)
    # print(prob.variables())
    # print(prob) # prints all constraints (and bounds)
    status = prob.solve()
    LpStatus[status]
    cnt = 0
    for i in range(0, len(ip_vars)):
        if value(ip_vars[i]) != 0.0:
            print(value(ip_vars[i]))
            cnt += 1
            print(f"storf:{storf_dict[i]}")
    print(cnt)



def read_fasta() -> list:
    '''
    convert the fasta file into a dictionary of StORFs
    return: dict of storfs
    '''
    unfiltered_storfs = []
    with open("../../testout/E-coli_output_no_filt.fasta") as storf_file:
        for line in storf_file:
            if line[0] == ">":
                unfiltered_storfs.append([line, next(storf_file)])
    return unfiltered_storfs


def get_values(unfiltered_storfs: dict):
    '''
    Heuristically defined value for each storf within the range 0-1
    '''
    # initially use length until IP function is implemented properly
    # todo rewrite lenght attribute already exists with the storf meta data
    largest_storf = max(unfiltered_storfs, key=lambda x: len(x[1]))
    ls_len = len(largest_storf[1])  # largest storf length
    for storf in unfiltered_storfs:
        # value of storf as percentage of ls_len
        storf.append(int(len(storf[1]))/ls_len)

    # TODO implement additional heuristics for value of StORF
    # GC content, overlap nt size, etc.
    return unfiltered_storfs


def main():
    uf_storfs = read_fasta()
    uf_value_storfs = get_values(uf_storfs)
    filtered_storfs = filter(uf_value_storfs)

    # output filtered StORFs to fasta file


if __name__ == '__main__':
    main()
