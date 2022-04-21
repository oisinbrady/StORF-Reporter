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
import argparse

from kpip_filter import get_con_groups

ABSPATH = os.path.abspath(__file__)
ABSPATH = os.path.dirname(ABSPATH)  # absolute path to working directory
OPTIONS = None  # argparse run-time commands

# N.b. see ../Genomes/x/x_metrics.csv for full x genome metrics
UNFILTERED_UR = [
    "input/bascillus_no_filt.fasta",
    "input/E-coli_no_filt.fasta",
    "input/caul_no_filt.fasta",
    "input/staph_no_filt.fasta",
    "input/myco_no_filt.fasta",
    "input/pseudo_no_filt.fasta",
]
KPIP_OUTPUT = [
    "kpip_output/bascillus_output.fasta",
    "kpip_output/E-coli_output.fasta",
    "kpip_output/caul_output.fasta",
    "kpip_output/staph_output.fasta",
    "kpip_output/myco_output.fasta",
    "kpip_output/pseudo_output.fasta",
]
BLASTX_FILES = [
    "metrics/input/blastx/blastx_output/bascillus_no_filt.out",
    "metrics/input/blastx/blastx_output/E-coli_no_filt.out",
    "metrics/input/blastx/blastx_output/caul_no_filt.out",
    "metrics/input/blastx/blastx_output/staph_no_filt.out",
    "metrics/input/blastx/blastx_output/myco_no_filt.out",
    "metrics/input/blastx/blastx_output/pseudo_no_filt.out",
]


def read_unfiltered(file_path: str) -> list:
    uf_storfs = []
    with open(file_path) as storf_file:
        for line in storf_file:
            if line[0] == ">":
                uf_storfs.append([line, next(storf_file)])
    return uf_storfs


def read_filtered(file_path: str) -> list:
    unfiltered_storfs = []
    with open(file_path) as storf_file:
        for line in storf_file:
            if line[0] == ">":
                unfiltered_storfs.append([line, next(storf_file)])
    return unfiltered_storfs


def read_blastx(file_path: str) -> list:
    matches = []
    with open(file_path) as storf_file:
        for line in storf_file:
            matches.append(line.split())
    return matches


def remove_embedded_storfs(blastx_matches: list) -> list:
    # remove any *fully* embedded StORFs
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


def get_hss_1(blastx_matches: list) -> list:
    # bests StORFs (most likely CDS)
    if not OPTIONS.filt:
        print("\tDetermining HSS_1 datasets...")
    bm = blastx_matches
    hss_1 = [m for m in bm 
             if float(m[2]) >= 50  # percentage identity >= 50
             and float(m[3]) > 100  # alignment length > 100
             and float(m[10]) <= float("1e-50")  # e-value <= 1e-50
             and float(m[11]) >= 50  # bit-score >= 50
            ]  
    if not OPTIONS.filt:
        print("\t\tFiltering duplicate & fully-embedded StORFs")
    hss_1 = remove_repeat_storfs(hss_1)
    hss_1 = remove_embedded_storfs(hss_1)
    return hss_1


def get_hss_2(blastx_matches: list) -> list:
    # StORFs less likely than hss_1 to be CDS but still with some significance
    if not OPTIONS.filt:
        print("\tDetermining HSS_2 datasets...")
    bm = blastx_matches
    hss_2 = [m for m in bm 
             if 40 <= float(m[2]) < 50  # percentage identity = [30, 50)
             and float(m[3]) > 100  # alignment length > 100
             and float("1e-50") < float(m[10]) < 0.01  # e-value = (1e-50, 0.01)
             and 40 <= float(m[11]) < 50  #  bit-score = [40, 50]
            ]
    if not OPTIONS.filt:
        print("\t\tFiltering duplicate & fully-embedded StORFs")
    hss_2 = remove_repeat_storfs(hss_2)
    hss_2 = remove_embedded_storfs(hss_2)
    return hss_2


def get_hss_groups(unfiltered_storfs: list, blastx_matches: list) -> list:
    # Reference links on statistics:
    # https://www.biostars.org/p/187230/
    # https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3820096/
    # matched_subset = monte_carlo_subsambple(blastx_matches, 0.1)  # 10% random subsample
    hss_1_matches = get_hss_1(blastx_matches)
    hss_2_matches = get_hss_2(blastx_matches)
    hss_1_sequences = get_blastx_sequences(hss_1_matches, unfiltered_storfs)
    hss_2_sequences = get_blastx_sequences(hss_2_matches, unfiltered_storfs)
    hss_groups = {"HSS_1": hss_1_sequences, "HSS_2": hss_2_sequences}
    return hss_groups


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


def get_overlaps(con_group: list) -> list:
    overlaps = []
    for j in range(0, len(con_group) - 1):
        storf_j_id = con_group[j][0]
        locus_j = storf_j_id[storf_j_id.find(":") + 1: storf_j_id.find("|")]
        for k in range(j, len(con_group) - 1):
            storf_k_id = con_group[k][0]
            locus_k = storf_k_id[storf_k_id.find(":") + 1: storf_k_id.find("|")]
            start_j, stop_j = get_start_stop(locus_j)
            start_k, stop_k = get_start_stop(locus_k)
            # if no overlap
            if start_j >= stop_k or stop_j <= start_k:
                continue
            else:
                # +1 needed for stop codon
                overlap = stop_k - start_j + 1 if start_k <= start_j else stop_j - start_k + 1
                overlaps.append(overlap)
    return overlaps


def gc_metric(contig_group: list) -> list:
    gc_percentages = []
    for storf in contig_group:
        sequence = storf[1]
        gc_count = len(re.findall('[GC]', sequence))
        nt_count = len(re.findall('[AGTC]', sequence))
        gc_percentages.append(gc_count / nt_count)
    return gc_percentages


def kmer_metric(storfs: list, k: int):
    # TODO
    return


def codon_metric(storfs: list, genome_name: str, hss=False):
    # TODO
    return


def write_csv_summary(cumulative: pd.DataFrame, filtered: list, unfiltered: list, all_hss_groups: list) -> None:
    file_path = ABSPATH + f"/metrics/output/summary_metrics.csv"
    header = [
        "mean StORF size(normalised)", "mean overlap size(normalised)", "unfiltered StORFs",
        "filtered StORFs", "total unfiltered HSS", "total filtered HSS", "total HSS lost",
        "HSS_1 unfiltered", "HSS_1 filtered", "HSS_1 lost", "HSS_2 unfiltered", "HSS_2 filtered",
        "HSS_2 lost"
    ]
    _,cumulative = normalise_dataframe(cumulative)
    data = []
    total_size = 0
    total_overlap_size = 0
    for storf_values in cumulative.values:
        total_size += storf_values[1]
        # N.b: overlap_size is an average if StORF has multiple overlaps
        total_overlap_size += storf_values[2]
    mean_storf_size = total_size / len(cumulative.values)
    mean_overlap_size = total_overlap_size / len(cumulative.values)
    data.extend([mean_storf_size, mean_overlap_size, len(unfiltered), len(filtered)])
    unfilt_hss_groups = [[0,0] for i in range(0, len(UNFILTERED_UR))]
    filtered_hss_groups = [[0,0] for i in range(0, len(UNFILTERED_UR))]
    for genome in range(0, len(UNFILTERED_UR)):
        for g, hss_group in enumerate(all_hss_groups[genome]):
            # get all HSS StORF groups in filtered output for each genome
            filtered_hss_groups[genome][g] = [s for s in all_hss_groups[genome][hss_group] if s in filtered]
            # get all HSS StORF groups in unfiltered input for each genome
            unfilt_hss_groups[genome][g] = [s for s in all_hss_groups[genome][hss_group] if s in unfiltered]
    # get cumulative HSS_1 and HSS_2 stats
    total_filtered_hss_1 = total_unfiltered_hss_1 = total_filtered_hss_2 = total_unfiltered_hss_2 = 0
    for genome in range(0, len(UNFILTERED_UR)):
        total_filtered_hss_1 += len(filtered_hss_groups[genome][0])
        total_filtered_hss_2 += len(filtered_hss_groups[genome][1])
        total_unfiltered_hss_1 += len(unfilt_hss_groups[genome][0])
        total_unfiltered_hss_2 += len(unfilt_hss_groups[genome][1])
    t_uf_hss = total_unfiltered_hss_1 + total_unfiltered_hss_2
    t_f_hss = total_filtered_hss_1 + total_filtered_hss_2
    data.extend([t_uf_hss, t_f_hss, t_uf_hss - t_f_hss])
    # HSS_1 unfiltered, "" filtered, "" lost
    data.extend([total_unfiltered_hss_1, total_filtered_hss_1,
                 total_unfiltered_hss_1 - total_filtered_hss_1])
    # HSS_2 unfiltered, "" filtered, "" lost
    data.extend([total_unfiltered_hss_2, total_filtered_hss_2,
                 total_unfiltered_hss_2 - total_filtered_hss_2])
    # write all stats to csv file
    with open(file_path, "w") as f:
        writer = csv.writer(f)
        writer.writerow(header)
        writer.writerow(data)


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
    con_groups = get_con_groups(storf_set)
    for group in con_groups:
        # get metrics for each contig StORF group
        gc_perc = gc_metric(group)
        lengths = [len(s[1]) for s in group]
        overlaps = get_overlaps(group)  # list of (1 to 1/Many)
        storf_ave_overlap = sum(overlaps) / len(overlaps) if len(overlaps) != 0 else 0
        group_metrics = [[] for i in group]
        for i in range(0, len(group_metrics)):
            if len(overlaps) > 0:
                storf_total_overlaps = len(overlaps)
                storf_ave_overlap = sum(overlaps) / storf_total_overlaps
            else:
                storf_ave_overlap = 0
                storf_total_overlaps = 0
            group_metrics[i].append([gc_perc[i], lengths[i], storf_ave_overlap, storf_total_overlaps])
        for storf_metrics in group_metrics:
            for metric in storf_metrics:
                all_storf_metrics.append(metric)
    data = pd.DataFrame(all_storf_metrics)
    if len(data) > 0:
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


def plot_ecdf(data: list, genome_name: str, x_title: str) -> None:
    file_name = x_title.replace(" ", "_")
    if genome_name == "all":
        graph_path = ABSPATH + f"/metrics/output/"
    else:    
        graph_path = ABSPATH + f"/metrics/output/{genome_name}/ecdf"
    g = sns.displot(data=data, kind="ecdf")
    g.set_axis_labels(f"{x_title} (max: {max(data)})", "Proportion")
    print(f"\t\twriting to: /metrics/output/{genome_name}/ecdf/{file_name.lower()}ecdf.png")
    g.figure.savefig(graph_path + f"/{file_name.lower()}ecdf.png")
    plt.close()


def normalise_dataframe(data_frame: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    original_data = data_frame.copy()
    for storf_metric in data_frame.columns[0:len(data_frame.columns) - 2]:
        # normalise the values, excluding set_type column
        data_frame[storf_metric] = data_frame[storf_metric] / data_frame[storf_metric].max()
    return original_data, data_frame


def plot_pca(data: pd.DataFrame, genome_name: str) -> None:
    summary = genome_name == "all"
    if summary:
        print("\nPlotting cumulative genome 2D PCAs...")
    else:
        print("\t\tplotting 2D PCAs...")
    graph_path = ABSPATH + f"/metrics/output/cumulative_pca" if summary else ABSPATH + f"/metrics/output/{genome_name}/pca"
    original_data, data = normalise_dataframe(data)
    pca_2d_perm = list(combinations(data.columns, 2))
    for perm in pca_2d_perm:
        pc1 = perm[0]
        pc2 = perm[1]
        print(f"\t\twriting to: metrics/output/{genome_name}_{pc2}_{pc1}_pca.png")
        #title = f"All Genomes: {pc2}_{pc1} PCA" if summary else f"{genome_name}: {pc2}_{pc1} PCA"
        palette = {
            "unfiltered StORFs": "tab:grey",
            "HSS_1": "tab:blue",
            "HSS_2": "tab:green",
            "HSS_3": "tab:orange",
            "HSS_4": "tab:red"
        }
        if pc1 == "set_type" or pc2 == "set_type":
            ax = sns.boxplot(x=pc2, y=pc1, data=data, palette=palette)
        else:
            ax = sns.scatterplot(x=pc2, y=pc1, hue="set_type", data=data, palette=palette, alpha=0.8)
            ax.set(xlabel=f'{pc2.replace("_", " ")} (max:{original_data[pc2].max()})')
        ax.set(ylabel=f'{pc1.replace("_", " ")} (max:{original_data[pc1].max()})')
        #ax.set_title(title)
        plt.savefig(graph_path + f"/{genome_name}_{pc2}_{pc1}_pca.png")
        plt.close()


def set_dataframe_group_identities(hss_data: list, group_identities: list) -> None:
    hss_data['set_type'] = 'N.a'
    for group_id, locus in enumerate(group_identities):
        start = locus[0]
        stop = locus[1] + 1
        for i in range(start, stop):
            hss_data.at[i, 'set_type'] = f"HSS_{group_id + 1}"


def get_hss_in_storfs(storfs: list, hss: list, group_locus: list) -> list:
    # TODO seriously need a refactor if i need a function like this...
    hss_groups = [[] for i in range(0, 4)]  # for the 4 hss categories
    for g, group in enumerate(group_locus):
        start = group[0]
        stop = group[1]
        for hss_storf in hss[start:stop]:
            for storf in storfs:
                if hss_storf[0] == storf[0]:
                    hss_groups[g].append(storf)
                    break
    return hss_groups


def print_hss_group_info(hss_groups) -> None:
    print_str = []
    print_str.append(f"\tHSS_1: {len(hss_groups.get('HSS_1'))}\n")
    print_str.append(f"\tHSS_2: {len(hss_groups.get('HSS_2'))}")
    sum_hss = len(hss_groups.get('HSS_1')) + len(hss_groups.get('HSS_2'))
    print(*print_str)
    print(f"\tTotal HSS: {sum_hss}\n")


def metric_one_genome() -> None:
    unfiltered_storfs = read_unfiltered(OPTIONS.unfilt)
    filt = read_filtered(OPTIONS.filt)
    blastx_matches = read_blastx(OPTIONS.blastx)
    uf_hss_groups = get_hss_groups(unfiltered_storfs, blastx_matches)
    print(f"unfiltered StORFs: {len(unfiltered_storfs)}")
    print_hss_group_info(uf_hss_groups)
    print(f"filtered StORFs: {len(filt)}")
    filt_hss_1 = [s for s in filt if s in uf_hss_groups["HSS_1"]]
    filt_hss_2 = [s for s in filt if s in uf_hss_groups["HSS_2"]]
    hss_groups = {"HSS_1": filt_hss_1, "HSS_2": filt_hss_2}
    print_hss_group_info(hss_groups)


def metrics() -> None:
    # for manual file simple metric
    if OPTIONS.unfilt is not None:
        metric_one_genome()
        exit()
    cumulative_dataframe = pd.DataFrame()
    # list of each genome containing a list of each HSS group
    all_hss_groups = [{"HSS_1": None, "HSS_2": None} for i in range(0, len(UNFILTERED_UR))]
    all_unfiltered = [[] for i in range(0, len(UNFILTERED_UR))]
    # For each genome, produce metrics
    for i in range(0, len(UNFILTERED_UR)):
        genome_name = UNFILTERED_UR[i][UNFILTERED_UR[i].find("/")+1:UNFILTERED_UR[i].find("_")]
        print(f"\nProcessing {genome_name} UR...")
        unfiltered_storfs = read_unfiltered(UNFILTERED_UR[i])
        all_unfiltered += unfiltered_storfs
        blastx_matches = read_blastx(BLASTX_FILES[i])
        hss_group_sequences = get_hss_groups(unfiltered_storfs, blastx_matches)
        hss_sequences = hss_group_sequences.get("HSS_1") + hss_group_sequences.get("HSS_2")
        print(f"\tMetricising {len(unfiltered_storfs)} unfiltered StORFs & {len(hss_sequences)} HSS")
        unfiltered_storfs = [s for s in unfiltered_storfs if s not in hss_sequences]
        uf_data = get_metrics_dataframe(unfiltered_storfs)
        uf_data["set_type"] = "unfiltered StORFs"
        hss_1_data = get_metrics_dataframe(hss_group_sequences.get("HSS_1"))
        hss_1_data["set_type"] = "HSS_1"
        hss_2_data = get_metrics_dataframe(hss_group_sequences.get("HSS_2"))
        hss_2_data["set_type"] = "HSS_2"
        hss_data = pd.concat([hss_2_data, hss_1_data])
        # add HSS seqeunces to cumulative HSS lists
        all_hss_groups[i]["HSS_1"] = hss_group_sequences.get("HSS_1")
        all_hss_groups[i]["HSS_2"] = hss_group_sequences.get("HSS_2")
        # combine HSS_1, HSS_2, and non-HSS (low CDS probability) StORFs
        all_storf_data = pd.concat([uf_data, hss_data])
        # add to cumulative StORF dataframe
        cumulative_dataframe = pd.concat([cumulative_dataframe, all_storf_data])
        plot_pca(all_storf_data, genome_name)
        print("\t\tPlotting StORF lengths ECDF...")
        plot_ecdf([len(s[1]) for s in unfiltered_storfs], genome_name, "All StORF lengths")
        print("\t\tPlotting HSS lengths ECDF...")
        plot_ecdf([len(s[1]) for s in hss_sequences], genome_name, "HSS lengths")
    # plot cumulative PCA of all test genomes
    plot_pca(cumulative_dataframe, "all")
    plot_ecdf([len(s[1]) for s in all_unfiltered], "all", "Cumulative StORF Size ECDF")
    all_filtered = []
    for file_name in KPIP_OUTPUT:
        all_filtered += read_filtered(file_name)
    write_csv_summary(cumulative_dataframe, 
                      all_filtered, 
                      all_unfiltered, 
                      all_hss_groups)


def init_argparse() -> None:
    parser = argparse.ArgumentParser(description='StORF filtering parameters')
    parser.add_argument('-filtered', action="store", dest='filt', help='Input FASTA File')
    parser.add_argument('-unfiltered', action="store", dest='unfilt', help='Input FASTA File')
    parser.add_argument('-blastx', action="store", dest='blastx', required=False,
                        help='Input FASTA File')
    global OPTIONS
    OPTIONS = parser.parse_args()
    if OPTIONS.filt is not None or OPTIONS.unfilt is not None or OPTIONS.blastx is not None:
        if OPTIONS.filt is None or OPTIONS.unfilt is None or OPTIONS.blastx is None:
            parser.error("The following must be specified together: -filt, -unfilt, -blastx")


def main():
    # TODO general refactor
    init_argparse()
    metrics()


if __name__ == '__main__':
    main()
