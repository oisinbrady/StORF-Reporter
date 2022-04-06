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


def read_unfiltered(file_name: str) -> list:
    uf_storfs = []
    with open(f"input/{file_name}") as storf_file:
        for line in storf_file:
            if line[0] == ">":
                uf_storfs.append([line, next(storf_file)])
    return uf_storfs


def read_filtered(file_name: str) -> list:
    unfiltered_storfs = []
    with open(f"kpip_output/{file_name}") as storf_file:
        for line in storf_file:
            if line[0] == ">":
                unfiltered_storfs.append([line, next(storf_file)])
    return unfiltered_storfs


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


def get_con_groups(blastx_matches: list) -> list:
    con_groups = [] 
    con_group = []
    total_matches = len(blastx_matches)
    for i, match in enumerate(blastx_matches):
        con_group.append(match)
        if is_next_new_group(match, i, total_matches):
            con_groups.append(con_group)
            con_group = []
        else:
            continue
    return con_groups


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


def get_hss(file_name: str, unfiltered_storfs: list) -> list[list, list]:
    print("\tDetermining HSS datasets...")
    # TODO need a lot of research on these metrics for HSS definition
    blastx_matches = read_blastx(file_name)
    # matched_subset = monte_carlo_subsambple(blastx_matches, 0.1)  # 10% random subsample
    # https://www.biostars.org/p/187230/
    print("\t\tFiltering duplicate & fully-embedded StORFs")
    high_bitscore_bms = remove_embedded_storfs(remove_repeat_storfs([m for m in blastx_matches if float(m[11]) >= 50]))
    medium_bitscore_bms = remove_embedded_storfs(remove_repeat_storfs([m for m in blastx_matches if 40 <= float(m[11]) < 50]))
    high_evalue_bms = remove_embedded_storfs(remove_repeat_storfs([m for m in blastx_matches if float(m[10]) < float("1e-50")]))
    medium_evalue_bms = remove_embedded_storfs(remove_repeat_storfs([m for m in blastx_matches if float("1e-50") < float(m[10]) < 0.01]))
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


def get_overlaps(con_group: list) -> list:
    overlaps = []
    for j in range(0, len(con_group) - 1):
        storf_j_id = con_group[j][1][0]
        locus_j = storf_j_id[storf_j_id.find(":") + 1: storf_j_id.find("|")]
        for k in range(j, len(con_group) - 1):
            storf_k_id = con_group[k][1][0]
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
    for _,storf in contig_group:
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


def write_csv_summary(cumulative: pd.DataFrame, unfiltered: list, filtered: list, genomes_hss_group_size: list, genomes_filtered_hss_group_size: list) -> None:
    # ["gc_percentage", "size", "average_overlap_size", "total_overlaps"]

    file_path = ABSPATH + f"/metrics/output/summary_metrics.csv"
    header = [
        "mean StORF size(normalised)", "mean overlap size(normalised)", "unfiltered StORFs",
        "filtered StORFs", "total unfiltered HSS", "total filtered HSS", "total HSS lost",
        "HSS_1 unfiltered", "HSS_1 filtered", "HSS_1 lost", "HSS_2 unfiltered", "HSS_2 filtered",
        "HSS_2 lost", "HSS_3 unfiltered", "HSS_3 filtered", "HSS_3 lost", "HSS_4 unfiltered",
        "HSS_4 filtered", "HSS_4 lost"
    ]
    data = []

    total_size = 0
    total_overlap_size = 0
    for storf_values in cumulative.values:
        total_size += storf_values[1]
        # N.b: overlap_size is an average if StORF has multiple overlaps
        total_overlap_size += storf_values[2]
    mean_storf_size = total_size/len(cumulative.values)
    mean_overlap_size = total_overlap_size/len(cumulative.values)
    data.extend([mean_storf_size, mean_overlap_size, len(unfiltered), len(filtered)])
    
    total_unfiltered_hss = 0
    total_filtered_hss = 0 
    total_filtered_hss_groups = [0,0,0,0]
    total_unfiltered_hss_groups = [0,0,0,0]
    for genome in range(0,5):
        total_unfiltered_hss += sum(genomes_hss_group_size[genome])
        total_filtered_hss += sum(genomes_filtered_hss_group_size[genome])
        for group in range(0,4):
            total_unfiltered_hss_groups[group] += genomes_hss_group_size[genome][group]
            total_filtered_hss_groups[group] += genomes_filtered_hss_group_size[genome][group]

    data.extend([total_unfiltered_hss, total_filtered_hss, total_unfiltered_hss - total_filtered_hss])
    data.extend([total_unfiltered_hss_groups[0], total_filtered_hss_groups[0] , total_unfiltered_hss_groups[0]- total_filtered_hss_groups[0]])
    data.extend([total_unfiltered_hss_groups[1], total_filtered_hss_groups[1], total_unfiltered_hss_groups[1]- total_filtered_hss_groups[1]])
    data.extend([total_unfiltered_hss_groups[2], total_filtered_hss_groups[2], total_unfiltered_hss_groups[2]- total_filtered_hss_groups[2]])
    data.extend([total_unfiltered_hss_groups[3], total_filtered_hss_groups[3], total_unfiltered_hss_groups[3]- total_filtered_hss_groups[3]])


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

def plot_ecdf(data: list, genome_name: str, x_title: str) -> None:
    file_name = x_title.replace(" ", "_")
    graph_path = ABSPATH + f"/metrics/output/{genome_name}/ecdf"
    g = sns.displot(data=data, kind="ecdf")
    g.set_axis_labels(x_title, "Proportion")
    print(f"\t\twriting to: /metrics/output/{genome_name}/ecdf/{file_name.lower()}ecdf.png")
    g.figure.savefig(graph_path + f"/{file_name.lower()}ecdf.png")
    plt.close()


def plot_pca(data: pd.DataFrame, genome_name: str) -> None:
    summary = genome_name == "all" 
    if summary:
        print("\nPlotting cumulative genome 2D PCAs...")
    else:
        print("\t\tplotting 2D PCAs...")
    graph_path = ABSPATH + f"/metrics/output/cumulative_pca" if summary else ABSPATH + f"/metrics/output/{genome_name}/pca"
    pca_2d_perm = list(combinations(data.columns, 2))
    for perm in pca_2d_perm:
        pc1 = perm[0]
        pc2 = perm[1]
        print(f"\t\twriting to: metrics/output/{genome_name}_{pc2}_{pc1}_pca.png")
        title = f"All Genomes: {pc2}_{pc1} PCA" if summary else f"{genome_name}: {pc2}_{pc1} PCA"
        #data = monte_carlo_subsambple(data, 0.2)
        #mplstyle.use('fast')  # automatically simplify plot for speed
        palette = {
            "unfiltered StORFs":"tab:grey",
            "HSS_1":"tab:blue",
            "HSS_2":"tab:green",
            "HSS_3":"tab:orange",
            "HSS_4":"tab:red"
        }
        if pc1 == "set_type" or pc2 == "set_type":
            ax = sns.boxplot(x=pc2, y=pc1, data=data, palette=palette)
        else:
            ax = sns.scatterplot(x=pc2, y=pc1, hue="set_type", data=data, palette=palette)
        ax.set_title(title)
        plt.savefig(graph_path + f"/{genome_name}_{pc2}_{pc1}_pca.png")
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
    kpip_output = [
        "bascillus_output.fasta",
        "E-coli_output.fasta",
        "caul_output.fasta",
        "staph_output.fasta",
        "myco_output.fasta"
    ]
    genomes_dataframe = pd.DataFrame()
    # list of each genome containing a list of each HSS group
    genomes_hss_group_size = [[0 for i in range(0,4)] for i in range(0,5)]
    genomes_filtered_hss_group_size = [[0 for i in range(0,4)] for i in range(0,5)]
    total_unfiltered = [[] for i in range(0,len(unfiltered_ur))]
    total_hss = [[] for i in range(0,len(unfiltered_ur))]
    # For each genome, produce metrics
    for i in range(0, len(unfiltered_ur)):
        genome_name = unfiltered_ur[i][:unfiltered_ur[i].find("_")]
        print(f"\nProcessing {genome_name} UR...")
        unfiltered_storfs = read_unfiltered(unfiltered_ur[i])
        total_unfiltered += unfiltered_storfs
        # High Significance StORF groups and their index location for DataFrame referencing
        hss_group_identities, hss = get_hss(unfiltered_ur[i], unfiltered_storfs)
        print(f"\tMetricising {len(unfiltered_storfs)} unfiltered StORFs & {len(hss)} HSS")
        total_hss += hss
        for group_id, locus in enumerate(hss_group_identities):
            start = locus[0]
            stop = locus[1]  if locus[0] != 0 else locus[1] + 1
            genomes_hss_group_size[i][group_id] += stop - start
            filtered_storfs = read_filtered(kpip_output[i])
            for hs_storf in hss[start:stop]:
                # get number of HSS_x that are in the filtered output
                for storf in filtered_storfs:
                    if hs_storf[0] == storf[0]:
                        genomes_filtered_hss_group_size[i][group_id] += 1
                        break
                    else:
                        continue
        unfiltered_storfs = [s for s in unfiltered_storfs if s not in hss]
        uf_data = get_metrics_dataframe(unfiltered_storfs)
        uf_data["set_type"] = "unfiltered StORFs"
        hss_data = get_metrics_dataframe(hss)
        set_dataframe_group_identities(hss_data, hss_group_identities)
        uf_hss_data = pd.concat([uf_data, hss_data])
        plot_pca(uf_hss_data, genome_name)
        print("\t\tPlotting StORF lengths ECDF...")
        plot_ecdf([len(s[1]) for s in unfiltered_storfs], genome_name, "All StORF lengths")
        print("\t\tPlotting HSS lengths ECDF...")
        plot_ecdf([len(s[1]) for s in hss], genome_name, "HSS lengths")
        genomes_dataframe = pd.concat([genomes_dataframe, uf_hss_data])
    genomes_dataframe = pd.DataFrame(genomes_dataframe)
    # print quartile, mean, metrics
    plot_pca(genomes_dataframe, "all")
    total_filtered = []
    for file_name in kpip_output:
        total_filtered += read_filtered(file_name)
    write_csv_summary(genomes_dataframe, total_unfiltered, total_filtered, genomes_hss_group_size, genomes_filtered_hss_group_size)


def main():
    # TODO general refactor
    metrics()


if __name__ == '__main__':
    main()
