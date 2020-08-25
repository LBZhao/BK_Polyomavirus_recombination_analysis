# Sequence_Data_Analysis

##**Step 1. Discard Non-viral DNA reads(Readsequence)**

During low molecular DNA harvesting, cellular DNA such as mitochondrial DNA was also harvested and sequenced. In this step, cellular DNA reads were omitted by comparing reads with BK polyomavirus (BKPyV) archetype genome template(AB211369.1 BK polyomavirus DNA, complete genome, isolate: Dik). FASTQ.GZ File path is passed to this script as an argument. BKPyV archetype genome sequence is required in FAS format. The control plasmid sequence is also included in this script.

  1. All read that matches template will be omitted because no recombination exists.
  2. Each read (251bp in length) was broken into ten 25-mers. If none of the ten 25-mers matches to BKPyV archetype genome sequence, the original read will be omitted.
  3. All reads that were not omitted in the first two steps will be collected to a file in FASTA format for BLAST later. FASTA file will be export to the working directory.

##**Step 2. BLAST (BLAST+)**
All sequences collected in the previous step will be blasted against BKPyV archetype genome sequence with local BLAST+ with output format "10 qseqid qstart qend sstart send frames sseq qseq stitle".

##**Step 3. Format BLAST result (Format_blast_result)**
Blasted results are further analyzed to get recombination results. Because BKPyV has a circular genome, extra code was added to resolve errors because of circular genome.
FASTA file used for BLAST is passed to this script as an argument.


##**Step 4. Circular diagram**
To generate a recombination map. An R script is used for NCCR recombination.
