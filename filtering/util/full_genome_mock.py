def get_locus(storf_meta):
    offset = len("UR_Stop_Locations=")
    i = storf_meta.find("UR_Stop_Locations=") + offset
    j = storf_meta.find(";Length") 
    return storf_meta[i:j].split("-")


def update_storf_meta(meta, storf_loc, con_iter_number):
    # initialise dummy storf meta information
    return f">Chromosome-UR:{storf_loc[0]}-{storf_loc[1]}|Chromosome-UR:0-460_{con_iter_number}|UR_Stop_Locations:18-177|Length:159|Frame:1|UR_Frame:1|Start_Stop=TAG|End_Stop=TAA|StORF_Type:Stop-ORF"


def write_fasta(uf_storfs):
    f = open(f"../storf_finder_output/ecoli_FULL_GENOME_mock.fasta", "w")
    for storf in uf_storfs:
        f.write(f"{storf[0]}\n{storf[1]}")


def main():
    # trying to mock input to be similar in format to:
    # >Chromosome-UR:18-177|Chromosome-UR:0-460_0|UR_Stop_Locations:18-177|
        # Length:159|Frame:1|UR_Frame:1|Start_Stop=TAG|End_Stop=TAA|StORF_Type:Stop-ORF
    uf_storfs = []
    with open("../../../ecoli_full_genome.fasta") as storf_file:
        for line in storf_file:
            if line[0] == ">":
                uf_storfs.append([line, next(storf_file)])
    print(len(uf_storfs))

    con_number = 0  # the iteration of a StORF in a contig
    previous_locus = None
    i = 0
    for meta, seq in uf_storfs:
        current_locus = get_locus(meta)
        if previous_locus == None:
            previous_locus = get_locus(meta)
            uf_storfs[i][0] = update_storf_meta(meta, current_locus, con_number)
        else:
            start_x = int(current_locus[0])
            stop_x = int(current_locus[1])
            start_y = int(previous_locus[0])
            stop_y = int(previous_locus[1])
            if start_y > stop_x or stop_y < start_x:  # if no overlap
                con_number = 0
                uf_storfs[i][0] = update_storf_meta(meta, current_locus, con_number)
            else:
                # if there is an overlap then in the same contig
                con_number += 1
                uf_storfs[i][0] = update_storf_meta(meta, current_locus, con_number)
            previous_locus = get_locus(meta)
        i += 1
    write_fasta(uf_storfs)
    

if __name__ == '__main__':
    main()