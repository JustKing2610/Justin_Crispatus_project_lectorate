# Phylogenetic analysis workflow

For this analyis, first 700 _L. crispatus_ genomes were downloaded from NCBI as see in *list*

next, Prokka was run on all genomes using: 
```
for file in /input/folder/genomes/*.fasta; do prokka --outdir "/output_prokka/$(basename "$file" .fasta)" "$file"; done
```
Then, the gff files in the prokka output were renamed after the folder name:
```
find . -type f -name "*.gff" -exec sh -c 'dir=$(dirname "{}"); cp "{}" "prokka_gff/$(basename "$dir").gff"' \;
```
then, a python [script](Python/extract_SlpH_sequences.py) was used to extract the positions of annotated SlpH from the prokka GFF, find the sequence between those positions in the assembly file and write them to a new file named after the basename of the prokka folder (named after GCA accession).

after these fasta files containing SlpH were combined to a single file:
```
cat *.fasta > combined_slph_for_alignment.fasta
```
This same strategy was utelized for the extraction of the SlpH locus of _L. helveticus_ to serve as an outgroup

next, a chewbbaca cgMLST scheme was created for all 700 genomes of _l. crispatus_.
ChewBBACA (version 3.3.1) was used. First, Busco was run on the genomes, determining the completeness of the genomes. The highest number of complete genes was determined to be 4020. It was decided to not include assemblies with a complete gene number lower then 3900, to increase the quality and completeness of the cgMLST scheme.
after running Busco, chewBBACA was run according to the documentation found here:https://github.com/B-UMMI/chewBBACA/blob/master/CHEWBBACA/docs/user/tutorials/chewie_step_by_step.rst this resulted in a cgMLST scheme consisting of the remaining genomes.

next, allele schema was made into ST types and assigned to the sequences in the fasta file. Since both files had the same order of genomes, the following [script](Python/assign_sequences_types.py) was used to assign sequence types to the GCA numbers. After, this [script](Python/change_identifiers_sequence_types.py) was used to add the ST types to the headers of the fasta file so they would show up in the subsequent analysis.

In an initial alignment was manually inspected and multiple sequences were discarded due to extreme divergence to the refseq. 

SlpH amplicons were sequenced and a consensus was made using samtools consensus 
```
input_folder="input_bam_sorted"; for bam_file in "$input_folder"/*.bam; do [ -e "$bam_file" ] && bam_basename=$(basename "$bam_file" .bam) && samtools consensus "$bam_file" > "${bam_basename}_consensus.fasta"; done 
```
The files somehow received the header of the refseq, probably happend during mapping. So the headers of these 11 files were manually changed to sample names for alignment. after these files were concatenated and added to the fasta file with all other sequences with ST types in the headers. only sample A2/A3 and all RL strains (RL09/10/11/17) were copied to this file, due to other samples either containing no crispatus, or multiple strains.

a mafft alignment using the G_INS_I algorithm was made with options with the following command
```
"/usr/bin/mafft"  --globalpair --maxiterate 16 --inputorder "amplicon_fasta_file.fasta"

```

This alignment file was then feeded into iqtree using the following command
```
iqtree -s alignment_file_slph_sequences -B 1000 -alrt 1000
```

and was visualised in ITOL
