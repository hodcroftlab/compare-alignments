# Compare Alignments

A simple snakemake setup that allows Nextclade/Nextalign to be run with a variety of parameters (more can easily be added by adding new rules) and then runs a comparison script between them (and an additional Maaft alignment if wished) to compare alignments.

## Quickstart

### Download comparison script
Download the [alignment comparison script](https://github.com/hodcroftlab/poliovirus_recombination/blob/main/scripts/evaluate_codon_alignments.py) from the polio repository. 
**As of 11 March 2025** it's recommended to use [the script from the branch with Emma's modifications](https://github.com/hodcroftlab/poliovirus_recombination/blob/improve-align-compare/scripts/evaluate_codon_alignments.py), not yet merged into the main repo.
Add this scripts to a folder within the repo called `scripts/`

### Add starting files
Ensure you add a starting unaligned fasta file, fasta reference, and reference as GFF3 format (can generate using [this script](https://github.com/nextstrain/nextclade_data/blob/master/docs/example-workflow/scripts/generate_from_genbank.py)), to an `input/` folder, and adjust the `files` rule to the names, if needed.

If you want to also include a Maaft alignment, add this to a `results/` folder and adjust the name (if needed) in the `files` rule. If you don't want to include a Maaft alignment, remove it from the input of the `evaluate_codon_alignments` rule.

### Adjust the CDS coordinates
This specifies which part of the whole-genome alignment should be compared - if using a virus that has a polyprotein, you can use this to remove the UTRs from the comparison. If you want to compare the entire alignment, remove the `--cds-coords` line from the `evaluate_codon_alignments` rule. 

Note this ability is intended for virus with a polyprotein (one concatinated coding region - like enteroviruses) and may not work well for other viruses. For example, non-coding regions between genes will remain (and thus may increase your 'stop-codon' counts). This may be ok as long as you consider this when interpreting results. However, if the coding frame changes between genes you will definitely get bad results and should probably use another method.

### Adjust the alignment parameters to be compared
If you wish to change the names of the alignments or add additional alignment parameter sets, just ensure you add the output file names to the `alignments` rule so that they are all included in the `evaluate_codon_alignments` rule at the end.
