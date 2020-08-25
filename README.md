# Sequence_Data_Analysis


Step 1. Discard Non-viral DNA reads(Readsequence)

During low molecular DNA harvesting, cellular DNA such as mitochondrial DNA were also harvested and sequenced. In this step, cellular DNA reads were omitted by comparing reads with BK polyomavirus (BKPyV) architype genome template(AB211369.1 BK polyomavirus DNA, complete genome, isolate: Dik). File path is passed to this script as argument. BKPyV architype genome sequence is required in FAS format. The control plasmid sequence is also included in this script.
1. All read that matches template will be omitted because no recombination exist.
2. Each read (251bp in length) was break into ten 25-mers.If none of the ten 25-mers matches to BKPyV architype genome sequence, the original read will be omitted.
3.All reads that were not omitted in the first two steps will be collected to a file in FASTA format for BLAST later. FASTA file will be export to the working directory.

Step 2. BLAST
All sequence collected in the previous step will be blasted against BKPyV architype genome sequence with local BLAST+ with output format "10 qseqid qstart qend sstart send frames sseq qseq stitle".

Step 3. Format BLAST result
Blasted results are further analyzed to get recombination result. Because BKPyV have a circular genome, extra code was added to resolve errors because of circular genome vs linear genome. 2bp Gaps are allowed.

Step 4. Format BLAST result
To generate recombination map. A R script is used for NCCR recombination.
