Directory:Description
blastx/output/: blastx match output for each StORF in unfiltered UR
output/: Graphs and CSV files processed from ../ur_metrics.py for each genome UR in ../input/


To reduce the size of filtered StORFs and still get as mant HSS as possible: 
Given that overlap and length distributions for HSS correlate with unfiltered StORFs
- Derive upper and lower quartiles average length & overlap metrics of unfiltered StORFs
- Bound kpip-filtering to [q1, q2] for size and overlap length
- IFF positive output:
    - Refactor and import the functions for deriving overlap/length metrics into kpip_filter.py

- Could try alternative of implementing elastic constraints using these [q1, q2] bounds
    - penalises StORF values in a contig not in [q1, q2] bound. 
    - requires a lot of considerations... not sure if I have time left to do that


Additional metric ideas:
	- complex motif/k-mer distributions
	- https://en.wikipedia.org/wiki/Alignment-free_sequence_analysis
