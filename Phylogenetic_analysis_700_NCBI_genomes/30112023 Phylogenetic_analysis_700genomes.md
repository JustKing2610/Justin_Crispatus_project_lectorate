# 30112023 Phylogenetic_analysis_700genomes
# introduction
For this analysis, 700 genomes, contigs and scaffolds were downloaded from NCBI as .fasta files. the question at hand was: Is the SlpH gene a suitable target to characterize and identify strains of _Lactobacillus crispatus_.
To answer this question, the 700 genomes were to be quality checked, after which the sequence of SlpH would be extracted from the assemblies containing the gene, a good quality MSA would be made before a phylogenetic tree would be constructed, giving a idea of the support values for each branch and indicating how much SlpH differs between strains.

# sequence assemblies processing
The 700 assemblies on NCBI for _Lactobacillus crispatus_ were downloaded as fasta's and transferred to our Asgard server for analysis and processing. 

# extraction of SlpH gene sequence in fasta format from the genomes
Before anything could be said and done with the SlpH genes of as many strains as possible, they first had to be extracted from the genomes downloaded from NCBI. For this, Prokka was used to annotate each assembly. The following bash command was used to run prokka on each genome assembly.

```
for file in /mnt/StudentFiles/2023/Justin/Phylogenetic_analysis/700_genomes_crispatus/*.fasta; do prokka --outdir "/mnt/StudentFiles/2023/Justin/Phylogenetic_analysis/analysis/output_prokka/$(basename "$file" .fasta)" "$file"; done
```
This bash command created directories for each assembly in wich the prokka result was saved. The next step was to rename the prokka-output GFF files to the basename of the fasta file, which was the GCA identifier. 
```
find . -type f -name "*.gff" -exec sh -c 'dir=$(dirname "{}"); cp "{}" "prokka_gff/$(basename "$dir").gff"' \;
```

The newly created folder contained the GFF files made by prokka named after the GCA identifier to tell them apart in further analysis. The next step was creating a python script to extract the positions of SlpH annotated by prokka and write the corresponding sequence in the assembly file with the same GCA identifier to a new file, named after the GCA identifier. the following script was written for this purpose:
```
from Bio.Seq import Seq
from Bio import SeqIO
from BCBio import GFF
from BCBio.GFF import GFFExaminer
import os
import glob


def extract_slpH_gene_positions(gff_file):
    target_feature_type = "CDS"
    for gff_file in gff_file:
        folder_name = os.path.basename(os.path.dirname(gff_file))
        output = f"{folder_name}_SlpH.fasta"
        with open(gff_file) as in_handle, open(output, 'w') as out_handle:
            for rec in GFF.parse(in_handle):
                for feature in rec.features:
                    if feature.type == target_feature_type and "Name" in feature.qualifiers and feature.qualifiers["Name"][0] == "slpH":
                        start = feature.location.start
                        end = feature.location.end
                        seq = str(Seq(feature.extract(rec.seq)))
                        out_handle.write(f">SlpH_{folder_name} \n{seq}\n")
                        print(f"Gene Name: slpH, Start: {start}, End: {end}")
                        

    gff_file = glob.glob("/mnt/StudentFiles/2023/Justin/Phylogenetic_analysis/analysis/output_prokka/GCA*/PROKKA_*.gff")
 
    extract_slpH_gene_positions(gff_file)
```

When this script was ran, a folder named extracted_SlpH_fastas was created containing individual fasta sequences of each assembly containing the SlpH sequence as annotated by prokka. After this step, 297 sequences were discarded as they did not contain SlpH.

```
cat *.fasta > combined_slph_for_alignment.fasta
```

#quality control of the extracted SlpH sequences
After aligning the combined SlpH sequences using MUSCLE in the MEGA11 application, a lot of big gaps were witnessed, and certain sequences started only at position 1300 or later. To eliminate the gaps and mismatching in the alignment, a script was written to calculate the average length, and write sequences above 1300 bases (refseq was 1356bases) to a new file called seq_right_length_slph.fasta. it also printed the amount of sequences written to this file, which were 311 sequences. 
```
from Bio import SeqIO

def calculate_average_seqlength(input_file):
    sequences = list(SeqIO.parse(input_file, "fasta"))
    total_length = sum(len(seq_record.seq) for seq_record in sequences)
    return total_length / len(sequences)

def filter_by_length(input_file, minimum_length):
    good_seqs = []
    for seq_record in SeqIO.parse(input_file, "fasta"):
        if len(seq_record.seq) > minimum_length:
            good_seqs.append(seq_record)
    return good_seqs

input_file = "combined_slph_for_alignment.fasta"
output_file = "seq_right_length.fasta"
min_length = 1300

filtered_sequences = filter_by_length(input_file, min_length)

average_length = calculate_average_seqlength(input_file)
print(f"Average sequence length: {average_length:.2f} bases")

SeqIO.write(filtered_sequences, output_file, "fasta")
print(f"{len(filtered_sequences)} sequences above {min_length} bases written to {output_file}")
```

When the resulting fasta was used in a MSA, all sequences aligned well, however, the some were longer than others, which called for some consideration as to which sequences should be included in the phylogenetic analysis

