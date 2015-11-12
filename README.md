# Detect_NDRs_and_flanking_nucleosomes
The basic units of DNA packaging are called nucleosomes. In the past decade, genome-wide nucleosome mapping experiments suggested a conserved nucleosome organization near gene promoters, consisting of regular nucleosome arrays on the gene bodies and a nucleosome-depleted region (NDR a.k.a. NFR) immediately upstream of the transcription start site (TSS). 

Using this MATLAB script we detect the locations of the nucleosome depleted regions and the flanking nucleosomes (+1/-1) corresponding to each gene promoter.

Example run:
    Detect_NDRs_and_flanking_nucleosomes('Dyads_AGH01-1_120_160.mat', 'AGH01-1')
