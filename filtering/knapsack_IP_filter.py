from pulp import LpVariable, LpProblem, LpMaximize, value

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
    overlapping_groups = []
    current_stop = None
    for storf in uf_value_storfs:
        if current_stop is None:
            # current_stop =
            return


def filter(uf_value_storfs: list) -> list:
    '''
    Filter StORFs favouring higher valued ones using integer programming to
    exclude overlapping StORFs
    '''
    # initialise lp instance
    prob = LpProblem("StORF Knapsack Problem", LpMaximize)
    # get overlapping StORFs
    overlapping_storfs = get_overlaps(uf_value_storfs)


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
    print(uf_storfs)
    uf_value_storfs = get_values(uf_storfs)
    filtered_storfs = filter(uf_value_storfs)

    # output filtered StORFs to fasta file


if __name__ == '__main__':
    main()
