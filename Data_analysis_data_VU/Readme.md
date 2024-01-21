# Mapping reads on SlpH from raw fasta (illumina paired-end) data
rename files, in this case:
```
for file in *.filt.fastq; do
    new_name=$(echo "$file" | awk -F'_' '{print $1"_"$2"_"$6".fastq"}')
    mv "$file" "$new_name"
done 
```
Then create an index of the SlpH gene fasta using BWA index:
```
bwa index L_crispatus_SlpH_fasta.fasta
```

continue with a loop to map every read to the index using BWA mem, make sure to map r1 and r2 to 1 sam file: 
```
for R1_file in *_R1.fastq; do R2_file="${R1_file/_R1/_R2}"; bwa mem L_crispatus_SlpH_fasta.fasta "$R1_file" "$R2_file" > "${R1_file%%_R1.fastq}_aligned.sam"; done
```
After, move all .sam files to a new folder and convert all the sam files to Bam files using samtools view:
```
mv *.sam sam_files_rosanne/
for f in sam_files_rosanne/*.sam; do n=$(basename "$f" .sam) echo "$n" samtools view -bS sam_files_rosanne/$n.sam | samtools sort - -o bam_files_rosanne/$n.bam done
```
then, use samtools consensus to create a concensus strand from each bam file: 
```
for f in bam_files_rosanne/*.bam; do n=$(basename "$f" .bam); echo "${n}"; samtools consensus -f fasta "$f" -o consensus_strands_data_rosanne/"${n}_consensus.fastq"; done
```
concatenate the consensus strands using cat
```
cat *.fasta > concatenated_consensus.fasta
```
Make a MSA using MAFFT
```
mafft concatenated_consensus.fasta > MSA_consensus.fasta
```
