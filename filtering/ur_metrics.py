import re
from itertools import product


def read_ur_fasta() -> list:
    unfiltered_storfs = []
    with open("../../testin/E-coli_output_no_filt.fasta") as storf_file:
        for line in storf_file:
            if line[0] == ">":
                unfiltered_storfs.append([line, next(storf_file)])
    return unfiltered_storfs


def read_prod_genes() -> list:
	# get a list of the highest bit-score match for each distinct StORF
	prod_storfs = []
	with open("../../testout/blastx/E-coli_no_filt.out") as storf_file:
		blastx_matches = [line.split() for line in storf_file]
	current_storf = None
	best_score_id = None
	best_score = None  # highest bit score of current StORF
	for i, match in enumerate(blastx_matches):
		storf_id = match[0][:match[0].find("|")]
		bit_score = float(match[11])
		# initial StORF
		if current_storf is None:
			current_storf = storf_id
			best_score = bit_score
			best_score_id = i
		# next Blastx match for same StORF instance
		elif current_storf == storf_id:  
			if bit_score > best_score:
				best_score = bit_score
				best_score_id = i
		# new StORF instance
		else:  
			# add the highest previous StORFs highest blasx bit score match
			prod_storfs.append(blastx_matches[best_score_id])
			current_storf = storf_id
			best_score = bit_score
			best_score_id = i
	return prod_storfs


def read_kpip_storfs() -> list:
	storfs = []
	with open("/home/oisin/Projects/mp/testout/oisin_output/output.fasta") as storf_file:
		for line in storf_file:
			if line[0] == ">":
				storfs.append([line, next(storf_file)])
	return storfs


def get_sequences(genes:list, storf_list:list) -> list:
	matching_storfs = []
	for gene in genes:
		# get gene/storf sequence id
		seq_id = gene[0][gene[0].find("|Chromosome")+1:gene[0].find("|UR_Stop_Locations")]
		for storf in storf_list:
			if seq_id in storf[0]:
				matching_storfs.append(storf)
				break

	return matching_storfs


def get_embedded_storf(locus_j: list, locus_k: list) -> list or None:
	# determine which StORF, if any, is embedded in the other
	start_j = int(locus_j[:locus_j.find("-")])
	stop_j = int(locus_j[locus_j.find("-") + 1:])
	start_k = int(locus_k[:locus_k.find("-")])
	stop_k = int(locus_k[locus_k.find("-") + 1:])
	if start_j > start_k and stop_j < stop_k:
		return locus_j
	elif start_k > start_j and stop_k < stop_j:
		return locus_k
	else: 
		return None


def remove_repeat_matches(kps_meta: list) -> list:
	sequence_ref = []
	duplicates = []
	for i, storf in enumerate(kps_meta):
		if storf[1] not in sequence_ref:
			sequence_ref.append(storf[1])
		else:
			duplicates.append(i)
	return [storf for i, storf in enumerate(kps_meta) if i not in duplicates]


def get_overlap(locus_j: list, locus_k: list) -> int or None:
	start_j = int(locus_j[:locus_j.find("-")])
	stop_j = int(locus_j[locus_j.find("-") + 1:])
	start_k = int(locus_k[:locus_k.find("-")])
	stop_k = int(locus_k[locus_k.find("-") + 1:])
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


def overlap_metric(storfs: list) -> None:
	# display stats on the overlapping behaviour of StORFs
	# get all overlaps
	overlaps = []
	for j in range(0, len(storfs) - 1):
		locus_j = storfs[j][0][storfs[j][0].find(":") + 1: storfs[j][0].find("|")]
		for k in range(j, len(storfs) - 1):
			locus_k = storfs[k][0][storfs[k][0].find(":") + 1: storfs[k][0].find("|")]
			overlap = get_overlap(locus_j, locus_k)
			if overlap is not None:
				overlaps.append(overlap)
	print(f"Mean overlap length: {round(sum(overlaps)/len(overlaps))}")


def len_metric(storfs: list) -> None:
	osl = [storf for storf in sorted(storfs, key=lambda x: len(x[1]))]
	total = sum(len(s[1]) for s in osl)
	mean = round(total/len(osl))
	largest_storf = osl[len(osl)-1]
	smallest_storf = osl[0]
	print(f"Smallest StORF: {len(smallest_storf[1])}\nLargest StORF: {len(largest_storf[1])}\nMean StORF length: {mean}\n")


def gc_metric(storfs: list) -> None:
	gc = {}
	for i, storf in enumerate(storfs):
		sequence = storf[1]
		gc_count = len(re.findall('[GC]', sequence))
		nt_count = len(re.findall('[AGTC]', sequence))
		storf_gc_percentage = gc_count/nt_count
		gc[f"{i}"] = storf_gc_percentage
	gc_sum = sum(s for s in gc.values())
	average_storf_gc_percentage = gc_sum/len(storfs)*100
	print(f"StORF average GC percentage: {average_storf_gc_percentage}")
	# get variance from genome gc
	genome_gc_total = 0
	nt_count = 0
	with open("../Genomes/E-coli/E-coli.fa") as genome:
		for line in genome:
			if line[0] != ">":
				nt_count += len(re.findall('[AGTC]', line))  # exclude undetermined nucleotides
				genome_gc_total += len(re.findall('[GC]', line))
	genome_gc = genome_gc_total/nt_count*100
	print(f"Original genome GC percentage: {genome_gc}")


def codon_metric(storfs: list) -> None:
	# get combinations of stop-stop codons in StORF
	codon_pair_count = {}
	start_codons = {}
	stop_codons = {}
	for storf in storfs:
		# count the stop-start codon pair
		start = storf[0].find("Start_Stop=")+len("Start_Stop=")
		start_codon = storf[0][start:start+3]
		stop = storf[0].find("End_Stop=")+len("End_Stop=")
		stop_codon = storf[0][stop:stop+3]

		if start_codon in start_codons.keys():
			start_codons[start_codon] += 1
		else:
			start_codons[start_codon] = 0

		if stop_codon in stop_codons.keys():
			stop_codons[stop_codon] += 1
		else:
			stop_codons[stop_codon] = 0

		codon_pair = (start_codon,stop_codon)
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
	print(f"Mode StORF start-end: {max(sorted_codon_pairs, key=sorted_codon_pairs.get)}")


def metrics(storfs: list) -> None:
	# see StORF-Reporter/Genomes/E-coli/E-coli_metrics.csv for genome metrics

	print("\n===========All UR E-coli StORF Metrics===========\n")
	# metrics for all StORFs
	len_metric(storfs)
	gc_metric(storfs)
	codon_metric(storfs)
	# overlap_metric(storfs)
	print("\n===========KPIP Filtered StORF Metrics===========\n")
	kpip_storfs = read_kpip_storfs()
	len_metric(kpip_storfs)
	gc_metric(kpip_storfs)
	codon_metric(kpip_storfs)
	# overlap_metric(kpip_storfs)
	
	print("\n===========Key Prodigal UR StORF Metrics===========\n")
	# get the highest bit-score match of each StORF
	prod_storfs = read_prod_genes()
	# get key prodigal StORFs above a certain bit-score (and E-value?)]
	kps_meta = [s for s in prod_storfs if float(s[11]) > 50.0 and float(f"{s[10]}") < float("1e-50")]
	kps_meta = remove_repeat_matches(kps_meta)
	# get relative StORF sequences of each kps
	kps_sequences = get_sequences(kps_meta, storfs)
	# filter out kps sequences that are fully embedded i.e. that produce the same matches
	kps_sequences = filter_embedded_storfs(kps_sequences)
	# metrics of all key prodigal StORFs
	len_metric(kps_sequences)
	gc_metric(kps_sequences)
	codon_metric(kps_sequences)
	overlap_metric(kps_sequences)

	print("\n===========Summary Stats===========\n")
	print(f"total StORFs in prod_storfs: {len(prod_storfs)}")
	print(f"total (non-embedded) StORFs with matching bit scores > 50.0 and e-values < 1e-50: {len(kps_sequences)}")
	print(f"\ntotal StORFs in kpip_storfs: {len(kpip_storfs)}")
	print(f"total StORFs in unfiltered_storfs: {len(storfs)}")


	# get similarities between kpip filtered and prodigal StORFs 
	keyprod_kpip_matches = get_sequences(kps_sequences, kpip_storfs)
	print(f"\nkey prodigal StORFs in kpip filtered StORFs: {len(keyprod_kpip_matches)}/{len(kpip_storfs)} ({len(keyprod_kpip_matches)/len(kpip_storfs) * 100}%)")

	print(f"\nkey prodigal StORFs in all StORFs: {len(kps_sequences)}/{len(storfs)} ({len(kps_sequences)/len(storfs) * 100}%)")

	print(f"\nkey prodigal StORFs lost using kpip: {(len(kps_sequences) - len(keyprod_kpip_matches))} ({(len(kps_sequences) - len(keyprod_kpip_matches)) / len(kps_sequences) * 100 }%)")

def main():
	storfs = read_ur_fasta()
	metrics(storfs)


if __name__ == '__main__':
	main()