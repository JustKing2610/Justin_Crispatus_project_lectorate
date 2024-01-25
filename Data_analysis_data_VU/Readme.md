# Mapping reads on SlpH from raw fasta (illumina paired-end) data
Then create an index of the SlpH gene from the reference genome (GCA_009769205.1) fasta using BWA index:
```
bwa index L_crispatus_SlpH_fasta.fasta
```

continue with a loop to map every read to the index using BWA mem, make sure to map r1 and r2 to 1 sam file: 
```
for R1_file in *_R1.fastq; do R2_file="${R1_file/_R1/_R2}"; bwa mem L_crispatus_SlpH_fasta.fasta "$R1_file" "$R2_file" > "${R1_file%%_R1.fastq}_aligned.sam"; done
```
After, move all .sam files to a new folder and convert all the sam files to Bam files using samtools view:
```
mv *.sam sam_files_folder/
for f in sam_files_folder/*.sam; do n=$(basename "$f" .sam) echo "$n" samtools view -bS sam_files_folder/$n.sam | samtools sort - -o bam_files_folder/$n.bam done
```
then, use samtools consensus to create a concensus strand from each bam file: 
```
for f in bam_files_folder/*.bam; do n=$(basename "$f" .bam); echo "${n}"; samtools consensus -f fasta "$f" -o consensus_strands_folder/"${n}_consensus.fasta"; done
```
concatenate the consensus strands using cat
```
cat *.fasta > concatenated_consensus.fasta
```
create an MSA using MAFFT
```
mafft concatenated_consensus.fasta > MSA_consensus.fasta
```
Visualise in a MSA viewer or run iqtree to construct a phylogenetic analysis.
