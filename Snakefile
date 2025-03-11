##############################
# Compare alignments by Maaft
# and Nextclade/Nextalign
# Using different parameters
###############################

# Note, using the cds_coords to only compare coding regions is simplified for viruses with a polyprotein
# (like enteroviruses) and may not work well for other viruses - may need to adjust for other viruses

# Quickstart:
# Download the comparison script and copy into 'scripts' folder from here: https://github.com/hodcroftlab/poliovirus_recombination/blob/main/scripts/evaluate_codon_alignments.py
# **Recommended to use Emma's modified version on the branch: https://github.com/hodcroftlab/poliovirus_recombination/blob/improve-align-compare/scripts/evaluate_codon_alignments.py

# Fill in the files you need in 'files' - include a maaft alignment if you want
#   (if not, remove this from the evaluate_codon_alignments input rule)
# Adjust parameters & output file names for each alignment rule, and ensure they match what's in 'alignments'
# Tip: generate GFF3 files using https://github.com/nextstrain/nextclade_data/blob/master/docs/example-workflow/scripts/generate_from_genbank.py


rule files:
    input:
        mafft_align = "results/aligned_maaft.fasta", #maaft alignment
        gff_reference = "input/annotation.gff", # GFF3 version of the reference seq
        gb_reference = "input/ev_d68_reference_genome.gb", # Genbank file of the reference seq
        unaligned = "input/nextstrain_filtered.fasta", # sequences you want to align
files = rules.files.input

cds_coords = "733,7299" # coordinates of the CDS in the reference sequence (see genbank file of reference)

rule all:
    input:
        "comparison-results/alignment_evaluation.csv",
        "comparison-results/alignment_evaluation.png"

# add the output for each alignment you want to compare, here
# so you can easily run them all at once by them being pulled into the comparison step
rule alignments:
    input:
        "results/aligned_hiv.fasta",
        "results/aligned_high-diversity.fasta",
        "results/aligned_a71.fasta"

# This just converts reference GB to fasta seq - if you already have this
# you can skip this rule by copying in your fasta reference file to 'results/reference_sequence.fasta'
rule reference_gb_to_fasta:
    message:
        """
        Converting reference sequence from genbank to fasta format
        """
    input:
        reference = files.gb_reference
    output:
        reference = "results/reference_sequence.fasta"
    run:
        from Bio import SeqIO 
        SeqIO.convert(input.reference, "genbank", output.reference, "fasta")

# Align with parameters designed for divergent HIV seqs
rule align_hiv_params: 
    message:
        """
        Aligning sequences to {input.reference} using Nextalign.
        """
    input:
        gff_reference = files.gff_reference,
        sequences = files.unaligned,
        reference = rules.reference_gb_to_fasta.output.reference
    output:
        alignment = "results/aligned_hiv.fasta"
    params:
            #HIV-nextclade
            penalty_gap_extend = 0, #make longer gaps more costly - default is 0
            penalty_gap_open = 8,  #make gaps more expensive relative to mismatches - default is 6
            penalty_gap_open_in_frame = 12, #make gaps more expensive relative to mismatches - default is 7
            penalty_gap_open_out_of_frame = 14, #make out of frame gaps more expensive - default is 8
            kmer_length = 10, #reduce to find more matches - default is 10
            kmer_distance = 50, #reduce to try more seeds - default is 50
            min_match_length = 40, #reduce to keep more seeds - default is 40
            allowed_mismatches = 8, #increase to keep more seeds - default is 8
            min_length = 30, # min_length - default is 100
            #cost of a mutation is 4
    shell:
        """
        nextclade run \
        {input.sequences} \
        --input-ref {input.reference} \
        --input-annotation {input.gff_reference} \
        --penalty-gap-open {params.penalty_gap_open} \
        --penalty-gap-extend {params.penalty_gap_extend} \
        --penalty-gap-open-in-frame {params.penalty_gap_open_in_frame} \
        --penalty-gap-open-out-of-frame {params.penalty_gap_open_out_of_frame} \
        --kmer-length {params.kmer_length} \
        --kmer-distance {params.kmer_distance} \
        --min-match-length {params.min_match_length} \
        --allowed-mismatches {params.allowed_mismatches} \
        --min-length {params.min_length} \
        --include-reference false \
        --output-fasta {output.alignment} 
        """

# Align with parameters designed for high diversity seqs
rule align_high_diversity: 
    message:
        """
        Aligning sequences to {input.reference} using Nextalign.
        """
    input:
        gff_reference = files.gff_reference,
        sequences = files.unaligned,
        reference = rules.reference_gb_to_fasta.output.reference
    output:
        alignment = "results/aligned_high-diversity.fasta"
    params:
            #high-diversity 
            penalty_gap_extend = 1, #make longer gaps more costly - default is 0
            penalty_gap_open = 13,  #make gaps more expensive relative to mismatches - default is 13
            penalty_gap_open_in_frame = 18, #make gaps more expensive relative to mismatches - default is 7
            penalty_gap_open_out_of_frame = 23, #make out of frame gaps more expensive - default is 8 # prev was 19
            kmer_length = 6, #reduce to find more matches - default is 10
            kmer_distance = 25, #reduce to try more seeds - default is 50
            min_match_length = 30, #reduce to keep more seeds - default is 40
            allowed_mismatches = 15, #increase to keep more seeds - default is 8
            min_length = 30, # min_length - default is 100
            #cost of a mutation is 4
    shell:
        """
        nextclade run \
        {input.sequences} \
        --input-ref {input.reference} \
        --input-annotation {input.gff_reference} \
        --penalty-gap-open {params.penalty_gap_open} \
        --penalty-gap-extend {params.penalty_gap_extend} \
        --penalty-gap-open-in-frame {params.penalty_gap_open_in_frame} \
        --penalty-gap-open-out-of-frame {params.penalty_gap_open_out_of_frame} \
        --kmer-length {params.kmer_length} \
        --kmer-distance {params.kmer_distance} \
        --min-match-length {params.min_match_length} \
        --allowed-mismatches {params.allowed_mismatches} \
        --min-length {params.min_length} \
        --include-reference false \
        --output-fasta {output.alignment} 
        """

# Align with parameters taken from the current A71 build (as of 11 Mar 2025)
rule align_a71: 
    message:
        """
        Aligning sequences to {input.reference} using Nextalign.
        """
    input:
        gff_reference = files.gff_reference,
        sequences = files.unaligned,
        reference = rules.reference_gb_to_fasta.output.reference
    output:
        alignment = "results/aligned_a71.fasta"
    params:
            #a71
            allowed_mismatches = 10, #allowed_mismatches - default is 8
            min_length = 30, # min_length - default is 100
            #cost of a mutation is 4
    shell:
        """
        nextclade run \
        {input.sequences} \
        --input-ref {input.reference} \
        --input-annotation {input.gff_reference} \
        --allowed-mismatches {params.allowed_mismatches} \
        --min-length {params.min_length} \
        --include-reference false \
        --output-fasta {output.alignment} 
        """


##############################

# Compare the results of the different codon-aware alignment tools by computing number of incomplete codons and stop codons
rule evaluate_codon_alignments:
    input:
        alignments = rules.alignments.input + [files.mafft_align],
    params:
        cdscoords = cds_coords
    output:
        csv = "comparison-results/alignment_evaluation.csv",
        plot = "comparison-results/alignment_evaluation.png",
    shell:
        """
        python scripts/evaluate_codon_alignments.py \
            --alignments {input.alignments} \
            --cds-coords {params.cdscoords} \
            --output {output.csv} \
            --plot {output.plot} 
        """
