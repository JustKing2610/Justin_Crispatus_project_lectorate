# Phylogenetic analysis workflow

For this analyis, first 700 _L. crispatus_ genomes were downloaded from NCBI as see in *list*

next, Prokka was run on all genomes using: 
```
for file in /700_genomes_crispatus/*.fasta; do prokka --outdir "/output_prokka/$(basename "$file" .fasta)" "$file"; done
```

then, a python [script](Python/extract_SlpH_sequences.py) was used to extract the positions of annotated SlpH from the prokka GFF, find the sequence between those positions in the assembly file and write them to a new file named after the basename of the assembly file.
