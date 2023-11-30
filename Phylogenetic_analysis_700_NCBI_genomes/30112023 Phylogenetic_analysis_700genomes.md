# 30112023 Phylogenetic_analysis_700genomes
# introduction
For this analysis, 700 genomes, contigs and scaffolds were downloaded from NCBI as .fasta files. the question at hand was: Is the SlpH gene a suitable target to characterize and identify strains of _Lactobacillus crispatus_.
To answer this question, the 700 genomes were to be quality checked, after which the sequence of SlpH would be extracted from the assemblies containing the gene, a good quality MSA would be made before a phylogenetic tree would be constructed, giving a idea of the support values for each branch and indicating how much SlpH differs between strains.

# sequence assemblies processing
The 700 assemblies on NCBI for _Lactobacillus crispatus_ were downloaded as fasta's and transferred to our Asgard server for analysis and processing. 

# extraction of SlpH gene sequence in fasta format from the genomes
Before anything could be said and done with the SlpH genes of as many strains as possible, they first had to be extracted from the genomes downloaded from NCBI. For this, Prokka was used to annotate each assembly. The following bash command was used to run prokka on each genome assembly.

```
asfada

```
