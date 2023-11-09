# Justin_Crispatus_project_lectorate
In this repository you will find analysis, scripts, plots, commands and more created during my internship at the lectorate analysis techniques in life sciences at avans university of applied sciences for the research project on lactobacillus crispatus and adaptive sampling


# Notes
## scripts and commands used

## preprocessing 

## concatenating reads
```bash

zcat 2324-012_Eefje_JI_LC_RBK_BC01-04_1/no_sample/20231004_1231_MN29431_FAV04503_fb5bc7e2/fastq_pass/barcode01/*.fastq.gz > fastq/ENR_01.fastq
zcat 2324-012_Eefje_JI_LC_RBK_BC01-04_1/no_sample/20231004_1231_MN29431_FAV04503_fb5bc7e2/fastq_pass/barcode02/*.fastq.gz > fastq/ENR_02.fastq
zcat 2324-012_Eefje_JI_LC_RBK_BC01-04_1/no_sample/20231004_1231_MN29431_FAV04503_fb5bc7e2/fastq_pass/barcode03/*.fastq.gz > fastq/ENR_03.fastq
zcat 2324-012_Eefje_JI_LC_RBK_BC01-04_1/no_sample/20231004_1231_MN29431_FAV04503_fb5bc7e2/fastq_pass/barcode04/*.fastq.gz > fastq/ENR_04.fastq
zcat 2324-012_Eefje_JI_LC_RBK_BC01-04_2/no_sample/20231004_1844_MN29431_FAV04503_1fb3955d/fastq_pass/barcode01/*.fastq.gz > fastq/DEP_01.fastq
zcat 2324-012_Eefje_JI_LC_RBK_BC01-04_2/no_sample/20231004_1844_MN29431_FAV04503_1fb3955d/fastq_pass/barcode02/*.fastq.gz > fastq/DEP_02.fastq
zcat 2324-012_Eefje_JI_LC_RBK_BC01-04_2/no_sample/20231004_1844_MN29431_FAV04503_1fb3955d/fastq_pass/barcode03/*.fastq.gz > fastq/DEP_03.fastq
zcat 2324-012_Eefje_JI_LC_RBK_BC01-04_2/no_sample/20231004_1844_MN29431_FAV04503_1fb3955d/fastq_pass/barcode04/*.fastq.gz > fastq/DEP_04.fastq
zcat 2324-012_Eefje_JI_LC_RBK_BC01-04_3/no_sample/20231005_0110_MN29431_FAV04503_c3668816/fastq_pass/barcode01/*.fastq.gz > fastq/NORM_01.fastq
zcat 2324-012_Eefje_JI_LC_RBK_BC01-04_3/no_sample/20231005_0110_MN29431_FAV04503_c3668816/fastq_pass/barcode02/*.fastq.gz > fastq/NORM_02.fastq
zcat 2324-012_Eefje_JI_LC_RBK_BC01-04_3/no_sample/20231005_0110_MN29431_FAV04503_c3668816/fastq_pass/barcode03/*.fastq.gz > fastq/NORM_03.fastq
zcat 2324-012_Eefje_JI_LC_RBK_BC01-04_3/no_sample/20231005_0110_MN29431_FAV04503_c3668816/fastq_pass/barcode04/*.fastq.gz > fastq/NORM_04.fastq
```
## concatenating reads Asgard server
```
zcat /mnt/StudentFiles/2023/Justin/Data_Usable/2324-012_Eefje_JI_LC_RBK_BC01-04_1/2324-012_Eefje_JI_LC_RBK_BC01-04_1/no_sample/20231004_1231_MN29431_FAV04503_fb5bc7e2/fastq_pass/barcode01/*.fastq.gz > fastq/ENR_01.fastq && \
zcat /mnt/StudentFiles/2023/Justin/Data_Usable/2324-012_Eefje_JI_LC_RBK_BC01-04_1/2324-012_Eefje_JI_LC_RBK_BC01-04_1/no_sample/20231004_1231_MN29431_FAV04503_fb5bc7e2/fastq_pass/barcode02/*.fastq.gz > fastq/ENR_02.fastq && \
zcat /mnt/StudentFiles/2023/Justin/Data_Usable/2324-012_Eefje_JI_LC_RBK_BC01-04_1/2324-012_Eefje_JI_LC_RBK_BC01-04_1/no_sample/20231004_1231_MN29431_FAV04503_fb5bc7e2/fastq_pass/barcode03/*.fastq.gz > fastq/ENR_03.fastq && \
zcat /mnt/StudentFiles/2023/Justin/Data_Usable/2324-012_Eefje_JI_LC_RBK_BC01-04_1/2324-012_Eefje_JI_LC_RBK_BC01-04_1/no_sample/20231004_1231_MN29431_FAV04503_fb5bc7e2/fastq_pass/barcode04/*.fastq.gz > fastq/ENR_04.fastq && \
zcat /mnt/StudentFiles/2023/Justin/Data_Usable/2324-012_Eefje_JI_LC_RBK_BC01-04_1/2324-012_Eefje_JI_LC_RBK_BC01-04_1/no_sample/20231004_1231_MN29431_FAV04503_fb5bc7e2/fastq_pass/barcode01/*.fastq.gz > fastq/DEP_01.fastq && \
zcat /mnt/StudentFiles/2023/Justin/Data_Usable/2324-012_Eefje_JI_LC_RBK_BC01-04_1/2324-012_Eefje_JI_LC_RBK_BC01-04_1/no_sample/20231004_1231_MN29431_FAV04503_fb5bc7e2/fastq_pass/barcode02/*.fastq.gz > fastq/DEP_02.fastq && \
zcat /mnt/StudentFiles/2023/Justin/Data_Usable/2324-012_Eefje_JI_LC_RBK_BC01-04_1/2324-012_Eefje_JI_LC_RBK_BC01-04_1/no_sample/20231004_1231_MN29431_FAV04503_fb5bc7e2/fastq_pass/barcode03/*.fastq.gz > fastq/DEP_03.fastq && \
zcat /mnt/StudentFiles/2023/Justin/Data_Usable/2324-012_Eefje_JI_LC_RBK_BC01-04_1/2324-012_Eefje_JI_LC_RBK_BC01-04_1/no_sample/20231004_1231_MN29431_FAV04503_fb5bc7e2/fastq_pass/barcode04/*.fastq.gz > fastq/DEP_04.fastq && \
zcat /mnt/StudentFiles/2023/Justin/Data_Usable/2324-012_Eefje_JI_LC_RBK_BC01-04_1/2324-012_Eefje_JI_LC_RBK_BC01-04_1/no_sample/20231004_1231_MN29431_FAV04503_fb5bc7e2/fastq_pass/barcode01/*.fastq.gz > fastq/NORM_01.fastq && \
zcat /mnt/StudentFiles/2023/Justin/Data_Usable/2324-012_Eefje_JI_LC_RBK_BC01-04_1/2324-012_Eefje_JI_LC_RBK_BC01-04_1/no_sample/20231004_1231_MN29431_FAV04503_fb5bc7e2/fastq_pass/barcode02/*.fastq.gz > fastq/NORM_02.fastq && \
zcat /mnt/StudentFiles/2023/Justin/Data_Usable/2324-012_Eefje_JI_LC_RBK_BC01-04_1/2324-012_Eefje_JI_LC_RBK_BC01-04_1/no_sample/20231004_1231_MN29431_FAV04503_fb5bc7e2/fastq_pass/barcode03/*.fastq.gz > fastq/NORM_03.fastq && \
zcat /mnt/StudentFiles/2023/Justin/Data_Usable/2324-012_Eefje_JI_LC_RBK_BC01-04_1/2324-012_Eefje_JI_LC_RBK_BC01-04_1/no_sample/20231004_1231_MN29431_FAV04503_fb5bc7e2/fastq_pass/barcode04/*.fastq.gz > fastq/NORM_04.fastq

```
## concatenating reads for ENR / DEP / NORM comparison 
```
cat /mnt/StudentFiles/2023/Justin/Data_Usable/fastq/ENR_*.fastq > fastq/fastq_combined_barcodes/ENR_BC01_04.fastq
cat /mnt/StudentFiles/2023/Justin/Data_Usable/fastq/DEP_*.fastq > fastq/fastq_combined_barcodes/DEP_BC01_04.fastq
cat /mnt/StudentFiles/2023/Justin/Data_Usable/fastq/NORM_*.fastq > fastq/fastq_combined_barcodes/NORM_BC01_04.fastq
```
## Running QC & porechop
```bash
for file in ./*.fastq;
  do n=$(basename $file .fastq);
  NanoPlot -t 8 -c pink -o nanoplot/${n}_rawreads --fastq $file;
  mkdir fastqc/${n}_rawreads;
  fastqc -t 8 -o fastqc/${n}_rawreads $file;
  porechop -v 2 -i $file > chopped/${n}.chopped.fastq;
  NanoPlot -t 8 -c pink -o nanoplot/${n}_chopped --fastq chopped/${n}.chopped.fastq;
  mkdir fastqc/${n}_chopped;
  fastqc -t 8 -o fastqc/${n}_chopped chopped/${n}.chopped.fastq;
done
```
## Kraken2 for chopped and filt fastq's 
```
for f in concatenate_chopped_NAS/input_filt/*.filt.fastq; do n=$(basename "$f" _BC01_BC04_chopped.filt.fastq) echo "$n" kraken2 -db ../../../kraken_standard/ 
             --threads 10 
             --unclassified-out kraken2/kraken_genus/unclassified/"${n}_unclass.fastq"
             --confidence 0.0 
             --report kraken2/kraken_genus/reports/"${n}.tsv"
             --out kraken_genus/"${n}_kraken.txt"
 done
```

## Running EPI2MElabs on server with docker containers
To install nextflow using conda, first create an environment where you want to install nextflow:
```
conda create -n nexflow_env
```
To install nextflow in the environment, run the following command:
```
conda activate nextflow_env
conda install nextflow
```
For this project, nextflow version 23.10.0 was used.

To run epi2me, the account on the server needs to recieve permission to run docker commands, for this you can contact your systems administrator. 
to check if docker works, run:
```
docker run hello-world
```
if it works, docker is succesfully set up, you may now proceed to running the epi2me labs workflow on a server!

to to this, use the following command and change filepaths and settings as desired or neccesary:
```bash
nextflow run epi2me-labs/wf-metagenomics --fastq /mnt/StudentFiles/2023/Justin/Data_Usable/2324-012_Eefje_JI_LC_RBK_BC01-04/analsis/Input_minimap_metagenomic_epi2me/ -with-docker 'epi2me/image:tag' -work-dir /mnt/StudentFiles/2023/Justin/Data_Usable/2324-012_Eefje_JI_LC_RBK_BC01-04/analsis/epi2me_output/Metagenomics_pipeline_RBK/ --threads 8 --database_set PlusPF-8
```

## Mapping reads on SlpH from raw fasta (illumina paired-end) data
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

Now, use IQTREE to make a .phy.iqtree file to use in the visualisation of the phylogenetic tree (in this case, bootstrapping was set using -B for ultrafast bootstrapping (1000 times) same with -altr for SH approximate likelihood ratio test 1000 times) 
```
iqtree -s MSA_consensus.fasta -B 1000 -alrt 1000
```
after gathering the iqtree outputs, visualize it using itol tree

# using prokka for annotation of 700 genomes crispatus for phylogenetic analysis
use the following loop to iterate over all the genomes in 1 folder and annotate them for further phylogenetic analysis
```
for f in ../700_genomes_crispatus/*.fasta; do base_name=$(basename "$f" .fasta); prokka "$f" --cpus 1 --outdir output_prokka/"$base_name" "$f";   echo "Prokka processing done for "$f"; done
```

# Count mapped reads
## to count reads mapped to a reference from bam files
use the following for loop:
```
for bam_file in *.bam; do     file_name=$(basename "$bam_file" .bam);      counts=$(samtools view -F260 "$bam_file" | cut -f3 | datamash -g 1 count 1);      echo "$file_name     $counts" >> "$output_file"; done
```
This command uses samtools view to filter and extract mapped reads from a BAM file, excluding unmapped and secondary alignments by excluding reads with specific flag bits (0x100 and 0x020). It then uses cut to extract the reference sequence names from the filtered output and datamash to count the occurrences of each unique name. The result is a summary of the number of mapped reads per reference sequence in the BAM file.
The output looks like the following when displayed on the command line where the first column is the basename (name of bam file), the second column is the refernce the reads were mapped to, and the third column includes the read count for each file:
```
Bella1_1_aligned     NZ_CP039266.1:164865-166220	6536
Bonita1_1_aligned     NZ_CP039266.1:164865-166220	4925
Campanarius1_1_aligned     NZ_CP039266.1:164865-166220	7650
Campanarius1_5_aligned     NZ_CP039266.1:164865-166220	6365
Chodsko1_1_aligned     NZ_CP039266.1:164865-166220	7393
Crisje1_1_aligned     NZ_CP039266.1:164865-166220	7224
Gaia1_2_aligned     NZ_CP039266.1:164865-166220	6063
Haifa1_1_aligned     NZ_CP039266.1:164865-166220	7089
Mechelea1_1_aligned     NZ_CP039266.1:164865-166220	5830
Nefesh1_10_aligned     NZ_CP039266.1:164865-166220	6530
Nefesh1_2_aligned     NZ_CP039266.1:164865-166220	6505
Nefesh1_3_aligned     NZ_CP039266.1:164865-166220	7893
Nefesh1_4_aligned     NZ_CP039266.1:164865-166220	6127
Nefesh1_5_aligned     NZ_CP039266.1:164865-166220	6962
Nefesh1_6_aligned     NZ_CP039266.1:164865-166220	7620
Nefesh1_7_aligned     NZ_CP039266.1:164865-166220	8748
Nefesh1_8_aligned     NZ_CP039266.1:164865-166220	4
Nefesh1_9_aligned     NZ_CP039266.1:164865-166220	6972
Paulie1_1_aligned     NZ_CP039266.1:164865-166220	6203
Reina1_2_aligned     NZ_CP039266.1:164865-166220	8083
Rhea1_1_aligned     NZ_CP039266.1:164865-166220	6043
Rivka2_1_aligned     NZ_CP039266.1:164865-166220	6429
Rivka2_2_aligned     NZ_CP039266.1:164865-166220	7194
Rivka2_3_aligned     NZ_CP039266.1:164865-166220	4903
Rivka2_4_aligned     NZ_CP039266.1:164865-166220	7581
Rivka2_7_aligned     NZ_CP039266.1:164865-166220	7949
Sucaria1_1_aligned     NZ_CP039266.1:164865-166220	9359
Tonna1_2_aligned     NZ_CP039266.1:164865-166220	6834
```
#Extracting SlpH from NCBI genomes
Download fasta files from all available genomes on NCBI for crispatus.

Run prokka using the following bash loop:
```
for f in ../700_genomes_crispatus/*.fasta; do base_name=$(basename "$f" .fasta); prokka "$f" --cpus 1 --outdir output_prokka/"$base_name" "$f";   echo "Prokka processing done for "$f"; done
```

Then run the "extract_positions_SlpH_GFF.py" script using python 3.10.

after, concatenate all the fasta file into 1 for alignment using:
```
cat *.fasta > cmbined.fasta
```
making sure to specify file paths
