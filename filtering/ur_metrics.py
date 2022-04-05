import csv
import os
import re
from itertools import permutations, combinations
from math import ceil

import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import matplotlib.style as mplstyle
import pandas as pd
import seaborn as sns
import random as rand

from kpip_filter import is_next_new_group, get_start_stop_bases

ABSPATH = os.path.abspath(__file__)
ABSPATH = os.path.dirname(ABSPATH)  # absolute path to working directory
PLOT = False  # if graphics are to be produced


def read_filtered(file_name: str) -> list:
    unfiltered_storfs = []
    with open(f"kpip_output/{file_name}") as storf_file:
        for line in storf_file:
            if line[0] == ">":
                unfiltered_storfs.append([line, next(storf_file)])
    return unfiltered_storfs


def read_unfiltered(file_name: str) -> list:
    uf_storfs = []
    with open(f"input/{file_name}") as storf_file:
        for line in storf_file:
            if line[0] == ">":
                uf_storfs.append([line, next(storf_file)])
    return uf_storfs


def read_blastx(file_name: str) -> list:
    matches = []
    file_name = file_name[:file_name.find(".")] + ".out"
    with open(f"metrics/input/blastx/output/{file_name}") as storf_file:
        for line in storf_file:
            matches.append(line.split())
    return matches


def get_graph_path(hss: bool, genome_name: str) -> str:
    if hss:
        graph_path = ABSPATH + f"/metrics/output/{genome_name}/hss_metrics"
    else:
        graph_path = ABSPATH + f"/metrics/output/{genome_name}/unfiltered_storf_metrics"
    return graph_path


def remove_embedded_storfs(blastx_matches: list) -> list:
    # remove any *fully* embedded StORFs
    print("\t\tFiltering embedded StORFs...")
    bm = blastx_matches
    embedded_storfs = []
    for j in range(0, len(bm)):
        locus_j = bm[j][0][bm[j][0].find(":") + 1: bm[j][0].find("|")]
        for k in range(j, len(bm)):
            locus_k = bm[k][0][bm[k][0].find(":") + 1: bm[k][0].find("|")]
            embedded_storf = get_embedded_storf(locus_j, locus_k)
            if embedded_storf == locus_k:
                embedded_storfs.append(k)
            elif embedded_storf == locus_j:
                embedded_storfs.append(j)
            else:
                continue
    return [storf for i, storf in enumerate(bm) if i not in embedded_storfs]


def get_hss(file_name: str, unfiltered_storfs: list) -> list[list, list]:
    print("\tDetermining HSS datasets...")
    # TODO general refactor
    # TODO plot HSS stats histogram
    # TODO tag each category in returned hss for distinction in PCA
    # TODO need a lot of research on these metrics for HSS definition
    blastx_matches = read_blastx(file_name)
    #blastx_matches = remove_embedded_storfs(blastx_matches)
    #blastx_matches = remove_repeat_matches(blastx_matches)
    #blastx_matches = remove_repeat_storfs(blastx_matches)
    # matched_subset = monte_carlo_subsambple(blastx_matches, 0.1)  # 10% random subsample
    # https://www.biostars.org/p/187230/
    print("\t\tFiltering duplicate StORFs")
    high_bitscore_bms = remove_repeat_storfs([m for m in blastx_matches if float(m[11]) >= 50])
    medium_bitscore_bms = remove_repeat_storfs([m for m in blastx_matches if 40 < float(m[11]) < 50])
    high_evalue_bms = remove_repeat_storfs([m for m in blastx_matches if float(m[10]) < float("1e-50")])
    medium_evalue_bms = remove_repeat_storfs([m for m in blastx_matches if float("1e-50") < float(m[10]) < 0.01])
   
    h_evalue_h_bitscore = [s for s in high_evalue_bms if s in high_bitscore_bms] # HSS group 1
    h_evalue_m_bitscore = [s for s in high_evalue_bms if s in medium_bitscore_bms] # " " 2
    m_evalue_h_bitscore = [s for s in medium_evalue_bms if s in high_bitscore_bms]  # " " 3
    m_evalue_m_bitscore = [s for s in medium_evalue_bms if s in medium_bitscore_bms]  # " " 4

    hss = h_evalue_h_bitscore+h_evalue_m_bitscore+m_evalue_h_bitscore+m_evalue_m_bitscore
    # start, stop for each HSS group
    hss_2_start = len(h_evalue_h_bitscore)
    hss_2_stop = hss_2_start + len(h_evalue_m_bitscore)
    hss_3_start = hss_2_stop
    hss_3_stop = hss_3_start + len(m_evalue_h_bitscore)
    hss_4_start = hss_3_stop
    group_identity_matrix = [
        [0, len(h_evalue_h_bitscore)], 
        [hss_2_start, hss_2_stop],
        [hss_3_start, hss_3_stop],
        [hss_4_start, len(hss)-1]
    ]

    hss = get_blastx_sequences(hss, unfiltered_storfs)
    return group_identity_matrix, hss


def get_blastx_sequences(blastx_matches: list, storf_list: list) -> list:
    matching_storfs = []
    for m in blastx_matches:
        # get gene/storf sequence id via its locus relative to UR
        seq_id = m[0][m[0].find(":") + 1: m[0].find("|")]
        for storf in storf_list:
            if seq_id in storf[0]:
                matching_storfs.append(storf)
                break
    return matching_storfs


def get_start_stop(locus: list) -> tuple[int, int]:
    start = int(locus[:locus.find("-")])
    stop = int(locus[locus.find("-") + 1:])
    return start, stop


def get_embedded_storf(locus_j: list, locus_k: list) -> list or None:
    # determine which StORF, if any, is embedded in the other
    start_j, stop_j = get_start_stop(locus_j)
    start_k, stop_k = get_start_stop(locus_k)
    if start_j > start_k and stop_j < stop_k:
        return locus_j
    elif start_k > start_j and stop_k < stop_j:
        return locus_k
    else:
        return None


def remove_repeat_storfs(blastx_matches: list) -> list:
    seen_storfs = []
    duplicates = []
    for i, storf in enumerate(blastx_matches):
        storf_id = storf[0]
        if storf_id not in seen_storfs:
            seen_storfs.append(storf_id)
        else:
            duplicates.append(i)
    return [storf for i, storf in enumerate(blastx_matches) if i not in duplicates]


def remove_repeat_matches(blastx_matches: list) -> list:
    print("\t\tFiltering repeat blastx matches")
    sequence_refs = []
    duplicates = []
    for i, storf in enumerate(blastx_matches):
        sequence_ref = storf[1]
        if sequence_ref not in sequence_refs:
            sequence_refs.append(sequence_ref)
        else:
            duplicates.append(i)
    return [storf for i, storf in enumerate(blastx_matches) if i not in duplicates]


def get_overlap(locus_j: list, locus_k: list) -> int or None:
    start_j, stop_j = get_start_stop(locus_j)
    start_k, stop_k = get_start_stop(locus_k)
    # if no overlap
    if start_j >= stop_k or stop_j <= start_k:
        return None
    else:
        # +1 needed for stop codon
        return stop_k - start_j + 1 if start_k <= start_j else stop_j - start_k + 1


def get_overlaps(con_group: list) -> list:
    overlaps = []
    for j in range(0, len(con_group) - 1):
        storf_j_id = con_group[j][1][0]
        locus_j = storf_j_id[storf_j_id.find(":") + 1: storf_j_id.find("|")]
        for k in range(j, len(con_group) - 1):
            storf_k_id = con_group[k][1][0]
            locus_k = storf_k_id[storf_k_id.find(":") + 1: storf_k_id.find("|")]
            overlap = get_overlap(locus_j, locus_k)
            if overlap is not None:
                overlaps.append(overlap)
    return overlaps


def gc_metric(contig_group: list) -> list:
    gc_percentages = []
    for _,storf in contig_group:
        sequence = storf[1]
        gc_count = len(re.findall('[GC]', sequence))
        nt_count = len(re.findall('[AGTC]', sequence))
        gc_percentages.append(gc_count / nt_count)
    return gc_percentages


def kmer_metric(storfs: list, k: int):
    kmers = [''.join(i) for i in list(permutations('AGTC', k))]
    kmer_count = {k: 0 for k in kmers}
    for storf in storfs:
        seq = storf[1]
        for kmer in kmers:
            kmer_count[kmer] = kmer_count.get(kmer) + seq.count(kmer)
    print(f"k_mer frequencies: {kmer_count}")
    df = pd.DataFrame(kmer_count.items())
    ax = sns.barplot(x=0, y=1, data=df)
    ax.set(xlabel=f'{k}-mer', ylabel='Frequency', title=f'{k}-mer distribution')
    ax.set_xticklabels(ax.get_xticklabels(), rotation=90)
    plt.show()


def codon_metric(storfs: list, genome_name: str, hss=False):
    graph_path = get_graph_path(hss, genome_name)
    # get combinations of stop-stop codons in StORF
    codon_pair_count = {}
    start_codons = {}
    stop_codons = {}
    for storf in storfs:
        # count the stop-start codon pair
        start_codon, stop_codon = get_start_stop_bases(storf[0])

        if start_codon in start_codons.keys():
            start_codons[start_codon] += 1
        else:
            start_codons[start_codon] = 0

        if stop_codon in stop_codons.keys():
            stop_codons[stop_codon] += 1
        else:
            stop_codons[stop_codon] = 0

        codon_pair = start_codon, stop_codon
        if codon_pair in codon_pair_count.keys():
            codon_pair_count[codon_pair] += 1
        else:
            codon_pair_count[codon_pair] = 0
    sorted_codon_pairs = {k: v for k, v in sorted(codon_pair_count.items(), key=lambda item: item[1], reverse=True)}
    sorted_start = {k: v for k, v in sorted(start_codons.items(), key=lambda item: item[1], reverse=True)}
    sorted_stop = {k: v for k, v in sorted(stop_codons.items(), key=lambda item: item[1], reverse=True)}
    print(f"Start codons: {sorted_start}")
    print(f"Stop codons: {sorted_stop}")
    print(f"Start-stop pair count: {sorted_codon_pairs}")
    print(f"Mode StORF start-end: {max(sorted_codon_pairs, key=sorted_codon_pairs.get)}\n")
    # plot pair frequencies
    pair_df = pd.DataFrame(sorted_codon_pairs.items())
    ax = sns.barplot(x=0, y=1, data=pair_df)
    ax.set(xlabel='Codon pair', ylabel='Frequency', title='Start-Stop codons')
    fig = ax.get_figure()
    fig.savefig(graph_path + f"/{genome_name}_pair_frequencies.png")
    # plot start codon frequencies
    # TODO this needs to read the second occurrence of a stop codon
    start_df = pd.DataFrame(sorted_start.items())
    ax = sns.barplot(x=0, y=1, data=start_df)
    ax.set(xlabel='Codon type', ylabel='Frequency', title='Start codons')
    fig = ax.get_figure()
    fig.savefig(graph_path + f"/{genome_name}_stop_codon_frequencies.png")
    # plot stop codon frequencies
    stop_df = pd.DataFrame(sorted_stop.items())
    ax = sns.barplot(x=0, y=1, data=stop_df)
    ax.set(xlabel='Codon type', ylabel='Frequency', title='Stop codons')
    fig = ax.get_figure()
    fig.savefig(graph_path + f"/{genome_name}_start_codon_frequencies.png")


def get_summary(unfilt_storfs: list, hss: list, filt_storfs: list) -> list:
    f_hss_intersect = [s for s in filt_storfs if s in hss]
    hss_in_f = len(f_hss_intersect)
    hss_lost = len(hss) - len(f_hss_intersect)
    return [len(unfilt_storfs), len(hss), len(filt_storfs), hss_in_f, hss_lost]


def write_csv_metrics(genome_name: str, data: list) -> None:
    file_path = ABSPATH + f"/metrics/output/{genome_name}/{genome_name}_metrics.csv"
    header = [
        "storf_set_type",
        "total_storfs",
        "mean_storf_size",
        "smallest_storf",
        "largest_storf",
        "total_overlaps",
        "mean_overlap_size"
    ]
    with open(file_path, "w") as f:
        writer = csv.writer(f)
        writer.writerow(header)
        for r in data:
            writer.writerow(r)


def write_csv_summary(data: list, genome_name: str) -> None:
    file_path = ABSPATH + f"/metrics/output/{genome_name}/{genome_name}_summary_metrics.csv"
    header = [
        "unfiltered_storfs",
        "hss",
        "filtered_storfs",
        "hss_in_filtered",
        "hss_lost_by_filtering",
    ]
    with open(file_path, "w") as f:
        writer = csv.writer(f)
        writer.writerow(header)
        writer.writerow(data)


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


def get_metrics_dataframe(storf_set: list) -> pd.DataFrame:
    # TODO this needs a refactor
    # ALSO, metrics for:
        # contig sizes
        # sense/anti-sense strand
        # overlap pair strands (same or different?)
    con_group = []  # StORF contig group
    overlaps = []
    set_length = len(storf_set)
    all_storf_metrics = []
    for s, storf in enumerate(storf_set):
        # get metrics for each contig StORF group
        if is_next_new_group(storf, s, set_length):
            gc_perc = gc_metric(con_group)
            lengths = [len(s[1][1]) for s in con_group]
            overlaps = get_overlaps(con_group)  # list of (1 to 1/Many)
            storf_ave_overlap = sum(overlaps)/len(overlaps) if len(overlaps) != 0 else "N.a"
            con_storf_metrics = [[] for i in con_group]
            for i in range(0, len(con_storf_metrics)):
                if len(overlaps) > 0:
                    storf_total_overlaps = len(overlaps)
                    storf_ave_overlap = sum(overlaps)/storf_total_overlaps
                else:
                    storf_ave_overlap = 0
                    storf_total_overlaps = 0
                con_storf_metrics[i].append([gc_perc[i], lengths[i], storf_ave_overlap, storf_total_overlaps]) 
            for storf_metrics in con_storf_metrics:
                for metric in storf_metrics:
                    all_storf_metrics.append(metric)
            con_group = [(s, storf)]  # start finding new contig-group
        else:
            con_group.append((s, storf))
    data = pd.DataFrame(all_storf_metrics)
    for storf_metric in data.columns:  # normalise the values
        data[storf_metric] = data[storf_metric] / data[storf_metric].max()
    data.columns = ["gc_percentage", "size", "average_overlap_size", "total_overlaps"]
    return data


def monte_carlo_subsambple(data: list, subsample_percentage: float) -> pd.DataFrame:
    samples = 0
    subsamples = []
    subsample_chosen_indices = []
    limit = round(len(data) * subsample_percentage)
    while samples <= limit:
        i = rand.randrange(0, len(data))
        if i in subsample_chosen_indices:  # repeat index
            continue
        subsample_chosen_indices.append(i)
        subsamples.append(data[i])
        samples += 1
    return subsamples


def plot(data: pd.DataFrame, genome_name: str) -> None:
    summary = genome_name == "all" 
    if summary:
        print("Plotting cumulative genome 2D PCAs...")
    else:
        print("\t\tplotting 2D PCAs...")
    graph_path = ABSPATH + f"/metrics/output/" if summary else ABSPATH + f"/metrics/output/{genome_name}/pca"
    pca_2d_perm = list(combinations(data.columns, 2))
    for perm in pca_2d_perm:
        pc1 = perm[0]
        pc2 = perm[1]
        print(f"\t\twriting to: metrics/output/{genome_name}_{pc1}_{pc2}_pca.png")
        title = f"All Genomes: {pc1}/{pc2} PCA" if summary else f"{genome_name}: {pc1}/{pc2} PCA"
        #data = monte_carlo_subsambple(data, 0.2)
        mplstyle.use('fast')  # automatically simplify plot for speed
        ax = sns.scatterplot(x=pc1, y=pc2, hue="set_type", data=data)
        ax.xaxis.set_major_locator(ticker.MultipleLocator(20))  # x axis intervals of 10(/100) 
        ax.xaxis.set_major_formatter(ticker.PercentFormatter())
        ax.set_title(title)
        plt.savefig(graph_path + f"/{genome_name}_{pc1}_{pc2}_pca.png")
        plt.close()
    

def set_dataframe_group_identities(hss_data: list, group_identities: list) -> None:
    hss_data['set_type'] = 'N.a'
    for group_id, locus in enumerate(group_identities):
        start = locus[0]
        stop = locus[1]
        for i in range(start, stop):
            hss_data.at[i, 'set_type'] = f"HSS_{group_id+1}"


def metrics() -> None:
    # N.b. see ../Genomes/x/x_metrics.csv for full x genome metrics
    unfiltered_ur = [
    "bascillus_no_filt.fasta",
    "E-coli_no_filt.fasta",
    "caul_no_filt.fasta",
    "staph_no_filt.fasta",
    "myco_no_filt.fasta"
    ]
    genomes_dataframe = pd.DataFrame()
    # For each genome, produce metrics
    for i in range(0, len(unfiltered_ur)):
        genome_name = unfiltered_ur[i][:unfiltered_ur[i].find("_")]
        print(f"\nProcessing {genome_name} UR...")
        unfiltered_storfs = read_unfiltered(unfiltered_ur[i])
        # High Significance StORF groups and their index location for DataFrame referencing
        group_identities, hss = get_hss(unfiltered_ur[i], unfiltered_storfs)
        print(f"\tMetricising {len(unfiltered_storfs)} unfiltered StORFs & {len(hss)} HSS")
        unfiltered_storfs = [s for s in unfiltered_storfs if s not in hss]
        uf_data = get_metrics_dataframe(unfiltered_storfs)
        uf_data["set_type"] = "unfiltered StORFs"
        hss_data = get_metrics_dataframe(hss)
        set_dataframe_group_identities(hss_data, group_identities)
        uf_hss_data = pd.concat([uf_data, hss_data])
        plot(uf_hss_data, genome_name)
        genomes_dataframe = pd.concat([genomes_dataframe, uf_hss_data])

    genomes_dataframe = pd.DataFrame(genomes_dataframe)
    plot(genomes_dataframe, "all")


def main():
    metrics()


if __name__ == '__main__':
    main()
