import csv
import os
import re
from itertools import permutations
from math import ceil

import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns

from kpip_filter import is_next_new_group, get_start_stop_bases

ABSPATH = os.path.abspath(__file__)
ABSPATH = os.path.dirname(ABSPATH)  # absolute path to working directory


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
    with open(f"metrics/blastx/output/{file_name}") as storf_file:
        for line in storf_file:
            matches.append(line)
    return matches


def get_graph_path(hss: bool, genome_name: str) -> str:
    if hss:
        graph_path = ABSPATH + f"/metrics/output/{genome_name}/hss_metrics"
    else:
        graph_path = ABSPATH + f"/metrics/output/{genome_name}/unfiltered_storf_metrics"
    return graph_path


def get_hss(file_name: str, unfiltered_storfs: list) -> list:
    blastx_matches = read_blastx(file_name)
    blastx_matches = remove_repeat_matches(blastx_matches)
    blastx_matches = remove_repeat_storfs(blastx_matches)
    # Filter by pident, e-value, bit-score
    # N.b. e-value and bit-score heavy reliance on database in use
    bm_sorted = sorted(blastx_matches, key=lambda s: (s[2], s[10], s[11]), reverse=True)
    hss = [s for s in bm_sorted[:round(0.01 * len(blastx_matches))]]  # get top 1%
    hss = get_blastx_sequences(hss, unfiltered_storfs)
    hss = filter_embedded_storfs(hss)
    return hss


def get_blastx_sequences(blastx_matches: list, storf_list: list) -> list:
    matching_storfs = []
    for m in blastx_matches:
        # get gene/storf sequence id via its locus relative to UR
        seq_id = m[m.find(":") + 1: m.find("|")]
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


def remove_repeat_storfs(hss: list) -> list:
    seen_storfs = []
    duplicates = []
    for i, storf in enumerate(hss):
        storf_id = storf.split()[0]
        if storf_id not in seen_storfs:
            seen_storfs.append(storf_id)
        else:
            duplicates.append(i)
    return [storf for i, storf in enumerate(hss) if i not in duplicates]


def remove_repeat_matches(hss: list) -> list:
    sequence_refs = []
    duplicates = []
    for i, storf in enumerate(hss):
        sequence_ref = storf.split()[1]
        if sequence_ref not in sequence_refs:
            sequence_refs.append(sequence_ref)
        else:
            duplicates.append(i)
    return [storf for i, storf in enumerate(hss) if i not in duplicates]


def get_overlap(locus_j: list, locus_k: list) -> int or None:
    start_j, stop_j = get_start_stop(locus_j)
    start_k, stop_k = get_start_stop(locus_k)
    # if no overlap
    if start_j >= stop_k or stop_j <= start_k:
        return None
    else:
        # +1 needed for stop codon
        return stop_k - start_j + 1 if start_k <= start_j else stop_j - start_k + 1


def filter_embedded_storfs(storfs) -> list:
    # remove any *fully* embedded StORFs
    embedded_storfs = []
    for j in range(0, len(storfs)):
        locus_j = storfs[j][0][storfs[j][0].find(":") + 1: storfs[j][0].find("|")]
        for k in range(j, len(storfs)):
            locus_k = storfs[k][0][storfs[k][0].find(":") + 1: storfs[k][0].find("|")]
            embedded_storf = get_embedded_storf(locus_j, locus_k)
            if embedded_storf == locus_k:
                embedded_storfs.append(k)
            elif embedded_storf == locus_j:
                embedded_storfs.append(j)
            else:
                continue
    return [storf for i, storf in enumerate(storfs) if i not in embedded_storfs]


def overlap_metric(storfs: list, genome_name: str, hss=False) -> list:
    graph_path = get_graph_path(hss, genome_name)
    con_group = []  # StORF contig group
    s_total = len(storfs)
    overlaps = []
    for s, storf in enumerate(storfs):
        if is_next_new_group(storf, s, s_total):
            # calculate all overlaps in previous contiguous group
            for j in range(0, len(con_group) - 1):
                storf_j_id = con_group[j][1][0]
                locus_j = storf_j_id[storf_j_id.find(":") + 1: storf_j_id.find("|")]
                for k in range(j, len(con_group) - 1):
                    storf_k_id = con_group[k][1][0]
                    locus_k = storf_k_id[storf_k_id.find(":") + 1: storf_k_id.find("|")]
                    overlap = get_overlap(locus_j, locus_k)
                    if overlap is not None:
                        overlaps.append(overlap)
            con_group = [(s, storf)]  # start finding new contig-group
        else:
            con_group.append((s, storf))
    set_type = "hss" if hss else "unfiltered"
    if len(overlaps) == 0:
        print(f"{set_type} set has no overlaps, skipping distribution plotting...")
        return [0, "N.A."]
    mean_overlap_len = int(ceil(sum(overlaps) / len(overlaps)))  # round up
    lm = [s for s in overlaps]
    print(f"writing {set_type} ecdf overlap plot ({len(overlaps)} overlaps)...")
    g = sns.displot(data=lm, kind="ecdf")  # computationally expensive job
    g.set_axis_labels("overlap length (nt)", "Proportion")
    g.figure.savefig(graph_path + f"/{genome_name}_overlap_len_ecdf.png")
    # for csv data printing
    return [len(overlaps), mean_overlap_len]


def len_metric(storfs: list, genome_name: str, hss=False) -> list:
    graph_path = get_graph_path(hss, genome_name)
    sl = [storf for storf in sorted(storfs, key=lambda x: len(x[1]))]
    total = sum(len(s[1]) for s in sl)
    mean = round(total / len(sl))
    largest_storf = sl[len(sl) - 1]
    smallest_storf = sl[0]
    # plot empirical cumulative distribution function, i.e., proportional StORF lengths
    lm = [len(s[1]) for s in sl]
    set_type = "hss" if hss else "unfiltered"
    print(f"writing {set_type} ecdf length plot ({total} storfs)...")
    g = sns.displot(data=lm, kind="ecdf")
    g.set_axis_labels("StORF length (nt)", "Proportion")
    g.figure.savefig(graph_path + f"/{genome_name}_len_ecdf.png")
    # for csv data printing
    return [total, mean, len(smallest_storf[1]), len(largest_storf[1])]


def gc_metric(storfs: list) -> None:
    gc = {}
    for i, storf in enumerate(storfs):
        sequence = storf[1]
        gc_count = len(re.findall('[GC]', sequence))
        nt_count = len(re.findall('[AGTC]', sequence))
        storf_gc_percentage = gc_count / nt_count
        gc[f"{i}"] = storf_gc_percentage
    gc_sum = sum(s for s in gc.values())
    average_storf_gc_percentage = gc_sum / len(storfs) * 100
    print(f"StORF average GC percentage: {average_storf_gc_percentage}")
    # get variance from genome gc
    genome_gc_total = 0
    nt_count = 0
    with open("../Genomes/E-coli/E-coli.fa") as genome:
        for line in genome:
            if line[0] != ">":
                nt_count += len(re.findall('[AGTC]', line))  # exclude undetermined nucleotides
                genome_gc_total += len(re.findall('[GC]', line))
    genome_gc = genome_gc_total / nt_count * 100
    print(f"Original genome GC percentage: {genome_gc}\n")
    # plot
    lm = [s for s in list(gc.values())]
    g = sns.displot(data=lm, bins=30)
    g.set_axis_labels("gc percentage", "Count")
    plt.show()


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


def metrics() -> None:
    # N.b. see ../Genomes/x/x_metrics.csv for full x genome metrics
    unfiltered_ur = [
        "bascillus_no_filt.fasta",
        "E-coli_no_filt.fasta",
        "caul_no_filt.fasta",
        "staph_no_filt.fasta",
        "myco_no_filt.fasta",
    ]
    kpip_output = [
        "bascillus_output.fasta",
        "E-coli_output.fasta",
        "caul_output.fasta",
        "staph_output.fasta",
        "myco_output.fasta"
    ]
    # For each genome, produce metrics
    for i in range(0, len(unfiltered_ur)):
        genome_name = unfiltered_ur[i][:unfiltered_ur[i].find("_")]
        print(f"\nprocessing {genome_name} metrics...")
        # Non-filtered StORFs
        unfiltered_storfs = read_unfiltered(unfiltered_ur[i])

        uf_data = len_metric(unfiltered_storfs, genome_name)
        uf_data += overlap_metric(unfiltered_storfs, genome_name)
        uf_data.insert(0, "unfiltered_storfs")
        # Get high significance StORFs (HSS)
        hss = get_hss(unfiltered_ur[i], unfiltered_storfs)
        hss_data = len_metric(hss, f"{genome_name}", hss=True)
        hss_data += overlap_metric(hss, f"{genome_name}", hss=True)
        hss_data.insert(0, "hss")
        csv_data = [uf_data, hss_data]
        write_csv_metrics(genome_name, csv_data)
        filtered_storfs = read_filtered(kpip_output[i])
        summary_csv_data = get_summary(unfiltered_storfs, hss, filtered_storfs)
        write_csv_summary(summary_csv_data, genome_name)


def main():
    metrics()


if __name__ == '__main__':
    main()
