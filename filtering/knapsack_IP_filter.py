# TODO read StORF files - all_StORFs Data type dict of lists?
    # each storf holds:
        # its meta-data (e.g. start-stop, UR id, etc.)
        # selected {0,1} - dictates if the StORF has been selected via the IP
        # value - the heuristic value of the StORF
            # represent as a normalised value 0-1 e.g. 0.95 is a very 'good' StORF

# TODO initialise IP instance

# TODO establish constraints
    # for each COMBINATION from overlapping StORFs group, e.g. x_1, x_2, x_3: (find_overlaps(storfs) -> overlapping_storfs)

    # For all storfs in a combination group: sum.selected = {0,1}
        # x_1.selected + x_2.selected + x_3.selected >=0
        # x_1.selected + x_2.selected + x_3.selected <=1

        # i.e. ONLY 1 storf of the combination group can be selected

# TODO establish bound(s)
    # x.selected C- {0,1} x = 1,2,...n
        # each StORF is either selected or not for filtering

# TODO establish objective function
    # MAXIMISE(SIGMA(S)), where S = x C- all_StORFs
    # e.g. MAX ({0,1}x_1.value + {0,1}x_2.value + {0,1}x_3.value)