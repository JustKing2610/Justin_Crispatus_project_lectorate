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

Prokka was also run on all 700 crispatus genomes collected from NCBI. When after this same script was run on all the gff files produced. a folder named extracted_SlpH_fastas was created containing individual fasta sequences of each assembly containing the SlpH sequence as annotated by prokka. After this step, 297 sequences were discarded as they did not contain SlpH.

To combine all the seperate gene .fasta files, they were concatenated into 1 file for alignment using the following command.

```
cat *.fasta > combined_slph_for_alignment.fasta
```

# quality control of the extracted SlpH sequences
After aligning the combined SlpH sequences using MUSCLE in the MEGA11 application, a lot of big gaps were witnessed, and certain sequences started only at position 1300 or later. To eliminate the gaps and mismatching in the alignment, a script was written to calculate the average length, and write sequences above 1300 bases (refseq was 1356bases) to a new file called seq_right_length_slph.fasta. it also printed the amount of sequences written to this file, which were 311 sequences. meaning 92 sequences were considered too short.
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
next, the newly written sequences longer than 1300bp were concatonated into 1 fasta containing all SlpH sequences from the remaining genomes.
```
cat *.fasta > combined_fastas.fasta
```
# creating an outgroup from _Lactobacillus helveticus_
To continue in the creation of a solid phylogenetic analysis of all remaining SlpH sequences from _L. crisptatus_, a good outgroup needed to be selected. For this, lactobacillus helveticus was chosen, based on an earlier made alignment which showed similarities, yet enough differences.

To start, the assembly of _L. helveticus_ GCA_003053085.1 was downloaded from NCBI. Prokka was run on this genome to extract the SlpH gene in the same way as for _L. crisptatus_. Since in previous alignments, the helveticus SlpH gene extracted directly from NCBI was noticably longer and did not contribute to a good multiple sequence alignment.

prokka was run using the following command
```
prokka --outdir prokka_out/ genome/ncbi_dataset/data/GCA_003053085.1/GCA_003053085.1_ASM305308v1_genomic.fna
```

to extract the SlpH sequence, the script used for _L. crispatus_ was slightly modified as seen below:
```

from Bio.Seq import Seq
from Bio import SeqIO
from BCBio import GFF
import os
import glob

def extract_slpH_gene_positions(gff_files):
    target_feature_type = "CDS"
    target_gene_name = "slpH"

    for gff_file in gff_files:
        folder_name = os.path.basename(os.path.dirname(gff_file))
        output = f"{folder_name}_SlpH.fasta"
        with open(gff_file) as in_handle, open(output, 'w') as out_handle:
            for rec in GFF.parse(in_handle):
                for feature in rec.features:
                    gene_name = feature.qualifiers.get('Name', [''])[0]
                    if feature.type.lower() == target_feature_type.lower() and gene_name.lower() == target_gene_name.lower():
                        start = feature.location.start
                        end = feature.location.end
                        seq = str(Seq(feature.extract(rec.seq)))
                        out_handle.write(f">SlpH_{folder_name} \n{seq}\n")
                        print(f"Gene Name: SlpH, Start: {start}, End: {end}")
                    else:
                        print("Skipped feature:", feature)

gff_files = glob.glob("/mnt/StudentFiles/2023/Justin/Phylogenetic_analysis/helveticus_refseq/prokka_out/helveticus_refseq/PROKKA_12072023.gff")

extract_slpH_gene_positions(gff_files)
```

# Multiple sequence alignment using MAFFT
when doing an initial alignment, it was witnessed that the helveticus gene started 162 bases earlier than annotated SlpH genes from crispatus strains. This anomaly was caused by a an earlier ATG site in the genomic assembly of helveticus, which prokka picked up to be the start codon. this was trimmed until the proper ATG which matched with the crispatus SlpH genes.

the E-INS-i strategy was chosen in mafft due to it's accuracy on similar types of alignments like on this dataset. as seen in the image below. This strategy seemed optimal since initial alignments showed a very similar profile.
![image](https://github.com/JustKing2610/Justin_Crispatus_project_lectorate/assets/127951903/4b64dfcf-0e44-48f1-8a52-a4f375e33b45)

Mafft was run on the remaining genomes longer than 1300bp. the command used was, this was created by following the menu of mafft by typing "mafft" in the command line. after specifying input, output, output format, alignment method.
```
mafft --genafpair  --maxiterate 16 --inputorder "../only_good_seqs_slph.fasta" > "only_good_seqs_slph_E_INS_I_alignment_final" 

```
However, after performing this initial alignment, it was seen through NCBI MSA viewer version 1.25.0, there were a lot of big gaps, created by sequences that had big insertions in random places, as seen in a small portion this alignment below:
![image](https://github.com/JustKing2610/Justin_Crispatus_project_lectorate/assets/127951903/7edfb074-76c7-467c-b2bb-94567a399100)

In this alignment, the very top sequence is the consensus, below that is a slph gene from _L. helveticus_ and below _L. helveticus_ is the _L. crispatus_ reference sequence. It is clear there are some very different sequences in this alignment. to improve the quality of the alignment and simultaneously create a more supported phylogenetic tree, all highly variable sequences, ande sequences with an insert that was not present in the majority of the sequences, where discarded and saved in a fasta file for later analysis. The discarded sequences and the reason are shown in the flowchart below:

![image](https://github.com/JustKing2610/Justin_Crispatus_project_lectorate/assets/127951903/d21a82e6-9728-4866-91c0-fca6551190ae)

After inspecting the alignment after discarding all possible sequences which did not match the set criteria, the following alignment was generated with the same parameters as before (only a part is shown): 
![image](https://github.com/JustKing2610/Justin_Crispatus_project_lectorate/assets/127951903/651cc557-92d7-4070-8cfe-58873f064f2d)
in this alignment, it is clear that all remaining sequences align much better, the only big gaps are created by the helveticus refseq (sequence 1) and crispatus refseq (sequence 2)

# Initial Phylogenetic analysis

After using mafft to align the sequences, the output file from mafft was used in iqtree to create a phylogenetic tree to analyse, using the following command.
```
iqtree -s only_good_seqs_slph_E_INS_I_alignment_final -B 1000 -alrt 1000

```
where -s specified the input
-B 1000 set the number of bootstrap replicates to 1000 to create a higher reliability

-alrt 1000 set the number of ultra fast bootstrap replicates to 1000 and stands for approximate likelihood ratio test which is a second method of assesing bootstrap support.

This created the following phylogenetic tree:

![image](https://github.com/JustKing2610/Justin_Crispatus_project_lectorate/assets/127951903/80fc5645-9a4c-4418-a4e1-0f37101bfd6e)

# adding Sequence types to the analysis with a self made scheme using chewbbaca
## running chewBBACA
ChewBBACA (version 3.3.1) was used to create a cgMLST schema from the crispatus genomes gathered from NCBI. First, Busco was run on the genomes, determining the completeness of the genomes. The highest number of complete genes was determined to be 4020. It was decided to not include assemblies with a complete gene number lower then 3900, to increase the quality and completeness of the cgMLST scheme. 

after running Busco, chewBBACA was run according to the documentation found here:https://github.com/B-UMMI/chewBBACA/blob/master/CHEWBBACA/docs/user/tutorials/chewie_step_by_step.rst 
this resulted in a cgMLST scheme consisting of the remaining genomes.

## Assigning sequence types using python
To assign the sequence Types to each genome, the cgMLST scheme (excel format) was ran through the following script, extracting the cgMLST profiles for each identifier, saving it to a dictonary. After which every unique dictonary got assigned a ST starting from 1. 

```
import pandas as pd

tsv_file = "/mnt/StudentFiles/2023/Justin/mlst_crispatus/tsv_MLST_scheme/cgMLST_profiles.tsv"

df = pd.read_csv(tsv_file, delimiter='\t')

dict_GCA = {}

#for each row, take the GCA ID from the column named "file" and join the values of the row together with "-"

for index, row in df.iterrows():
    GCA_ID = row["FILE"]

    row_values = list(row[1:])

    row_string = '-'.join(map(str, row_values))

    dict_GCA[GCA_ID] = row_string

results_df = pd.DataFrame(list(dict_GCA.items()), columns=["GCA", 'Allele_values'])


results_df.to_excel("/mnt/StudentFiles/2023/Justin/mlst_crispatus/output_python_mlst_scheme/output_dict_Alleleschema.xlsx", index=False)



values = list(dict_GCA.values())

labels = pd.factorize(values)[0] + 1

ST_types_label = [f'ST{label}' for label in labels]

final_dict = dict(zip(dict_GCA.keys(), ST_types_label))


print("Check final dictonary for GCA numbers and ST types:", final_dict)

ST_dataframe = pd.DataFrame(list(final_dict.items()), columns=["Identifier", "Sequence_types"])

ST_dataframe.to_excel("/mnt/StudentFiles/2023/Justin/mlst_crispatus/output_python_mlst_scheme/ST_types_Crispatus_chewbacca.xlsx", index=False)
```

This created an excel file with the following format: | Identifier | Sequence Type |

![image](https://github.com/JustKing2610/Justin_Crispatus_project_lectorate/assets/127951903/48b4a673-2a32-4a21-b44e-0d14084d4060)


Next, these assigned sequence types should be added to the identifiers in the fasta file, this way, they will show up in the alignment and subsequent phylogenetic analysis, which allows for better visualisation of clusters. For this, another python script was written, this script Took a file of all previous selected and quality controlled sequences in fasta format and the previously described excel file containing Identifiers and sequence types. 

The read the excel file as a panda DataFrame, then created a dictonary of the identifiers and sequence types and matched the identifiers in the excel file to the identifiers in the fasta file. When it found a matching GCA identifier, the sequence type assigned in the excel file would be parsed into a modified header in the Fasta file and written to a new fasta file. 
```
from Bio import SeqIO
import pandas as pd
import re

# Path to the input FASTA file
fasta_file_path = "/mnt/StudentFiles/2023/Justin/mlst_crispatus/ST_all_700_Genomes/combined_slph_for_alignment.fasta"

# Path to the Excel file with identifiers and sequence types
excel_file_path = "/mnt/StudentFiles/2023/Justin/mlst_crispatus/test_new_identifiers/output_python_mlst_scheme/ST_types_Crispatus_chewbacca.xlsx"

# Path to the output FASTA file
output_fasta_path = "/mnt/StudentFiles/2023/Justin/mlst_crispatus/ST_all_700_Genomes/Assigned_st_all_extracted_slph.fasta "

# Read Excel file into a pandas DataFrame
df = pd.read_excel(excel_file_path, names=['Identifier', 'Sequence_types'])

# Create a dictionary to map identifiers to sequence types
identifier_map = dict(zip(df['Identifier'], df['Sequence_types']))
 
# Open the output FASTA file for writing
with open(output_fasta_path, 'w') as output_fasta:
    # Iterate over sequences in the input FASTA file
    for record in SeqIO.parse(fasta_file_path, 'fasta'):
        # Use regular expression to extract everything starting from 'GCA_' until the first period
        match = re.match(r"^.*?(GCA_\d+)\..*$", record.id)
        
        if match:
            base_identifier = match.group(1)
            
            # Check if the base identifier is in the mapping dictionary
            if base_identifier in identifier_map:
                # Concatenate the sequence type to the modified identifier
                new_identifier = f'>{identifier_map[base_identifier]}_{base_identifier}'
                
                # Create a new SeqRecord with the updated identifier and the original sequence
                updated_record = record
                updated_record.id = new_identifier
                updated_record.description = ''
                
                
                SeqIO.write(updated_record, output_fasta, 'fasta')
                
                
                print(f"Updated: {new_identifier}")
        else:
            print(f"Skipping record with invalid identifier format: {record.id}")

print(f"Output written to {output_fasta_path}")
```
The running of this script lead to the loss of 81 sequences. Which was caused by these sequences being lost in the cgMLST scheme creation. So this alignment continued with 127 sequences This lead to an updated flowchart with discarded sequences: 
![Genome _ SlpH quality assesment  (1)](https://github.com/JustKing2610/Justin_Crispatus_project_lectorate/assets/127951903/ba303010-d812-4dbb-a171-6f033a6e60b8)

# Final alignment and phylogenetic analysis of SlpH 
First, Sequences with the same sequence type were discarded randomly and 1 was left. this was done by hand. 

The alignment with 1 of each sequence type was achieved, simultaneously, our sequenced SlpH amplicons were added to this alignment.
## alignment consensus sequenced amplicons
For the addition of the consensus, the raw sequenced reads were first filterered, with the read length criteria between 900 and 1200, since the amplicons were 1150bp long, this was done using chopper:
```
for i in {1..11}; do   barcode="barcode$(printf "%02d" $i)";   input_file="../../input/$barcode/$barcode*.fastq";   output_file="${barcode}_filtered.fastq";   cat $input_file | chopper -l 900 --maxlength 1200 -t 4 > $output_file; done
```

After this, a concensus was made using EPI2MElabs wf-amplicons using the following command for the filtered raw data:
```

```


