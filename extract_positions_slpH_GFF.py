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

