import re
from itertools import product


def read_ur_fasta() -> list:
    unfiltered_storfs = []
    with open("../../testin/E-coli_output_no_filt.fasta") as storf_file:
        for line in storf_file:
            if line[0] == ">":
                unfiltered_storfs.append([line, next(storf_file)])
    return unfiltered_storfs


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
	# get usage of each stop codon
	taa = tag = tga = 0
	stop_codons = ["TAA", "TAG", "TGA"]
	# get combinations of start-stop, stop-stop codons in StORF
	stpperm = product(stop_codons, repeat=2)
	codon_delimit_count = dict()
	for perm in stpperm:
		codon_delimit_count[perm] = 0
	for storf in storfs:
		# count the stop-start codon pair
		start = storf[0].find("Start_Stop=")+len("Start_Stop=")
		start_codon = storf[0][start:start+3]
		if start_codon == "TAA":
			taa += 1
		elif start_codon == "TAG":
			tag += 1
		else:
			tga += 1
		stop = storf[0].find("End_Stop=")+len("End_Stop=")
		stop_codon = storf[0][stop:stop+3]
		if stop_codon == "TAA":
			taa += 1
		elif stop_codon == "TAG":
			tag += 1
		else:
			tga += 1
		codon_id = (start_codon,stop_codon)
		codon_delimit_count[codon_id] += 1
	print(f"TAA count: {taa}\nTAG count: {tag}\nTGA count: {tga}")
	print(f"Start-stop pair count: {codon_delimit_count}")
	print(f"Mode StORF start-end: {max(codon_delimit_count, key=codon_delimit_count.get)}")


def metrics(storfs: list) -> None:
	len_metric(storfs)
	gc_metric(storfs)
	codon_metric(storfs)


def main():
	storfs = read_ur_fasta()
	metrics(storfs)


if __name__ == '__main__':
	main()