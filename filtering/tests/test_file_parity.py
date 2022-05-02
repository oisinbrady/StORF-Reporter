def read_file(file_path: str) -> list:
    file_storfs = []
    bases = ["A", "G", "C", "T", "N"]
    with open(file_path) as storf_file:
        for line in storf_file:
            if line[0] in bases:
                file_storfs.append(line)
    return file_storfs


def get_matches(file_1_name, file_2_name):
    file_1_ids = []
    with open(file_1_name) as storf_file:
        for line in storf_file:
            if line[0] == ">":
                storf_id = line[0:line.find("|")]
                file_1_ids.append(storf_id)
    file_2_ids = []
    with open(file_2_name) as storf_file:
        for line in storf_file:
            if line[0] == ">":
                storf_id = line[0:line.find("|")]
                file_2_ids.append(storf_id)
    matches = len([s for s in file_1_ids if s in file_2_ids])
    print(f"file 1 length = {len(file_1_ids)}")
    print(f"file 2 length = {len(file_2_ids)}")
    print(f"matches:{matches}/{max((len(file_1_ids), len(file_2_ids)))}")


def main():
    file_1_name = "../kpip_output/test/test"
    file_2_name = "../storf_finder_output/storf_finder_ecoli_out.fasta"
    file_1_storfs = read_file(file_1_name)
    file_2_storfs = read_file(file_2_name)
    # print(len([s for s in file_1_storfs if s in file_s_storfs]) # what if duplicates? - use id instead
    get_matches(file_1_name, file_2_name)
    for line, a in enumerate(file_1_storfs):
        try:
            b = file_2_storfs[line]
            assert a == b
        except AssertionError:
            print(f"Mismatch on line:{line} ({line*2} on input fasta) between:\n{a}\nand\n{b}")
            exit()
                
 

if __name__ == '__main__':
    main()