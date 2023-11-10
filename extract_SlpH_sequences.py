from Bio import SeqIO
from BCBio.GFF import GFFExaminer
from BCBio import GFF

positions = []

gff = "/mnt/StudentFiles/2023/Justin/Phylogenetic_analysis/analysis/test_python_GFF_parse/PROKKA_11012023.gff"

with open(gff) as handle:
    for rec in GFF.parse(handle):
        for feature in rec.features:
            if feature.type == "CDS" and feature.qualifiers.get("Name") == ["slpH"]: 
                start = feature.location.start
                end = feature.location.end
                positions.append((start, end))

print(positions)

fasta = "/mnt/StudentFiles/2023/Justin/Phylogenetic_analysis/700_genomes_crispatus/GCA_000091765.1_ASM9176v1_genomic.fasta"
sequence = []

record = SeqIO.read(fasta, "fasta")

with open("/mnt/StudentFiles/2023/Justin/Phylogenetic_analysis/analysis/test_python_GFF_parse/GCA_000091765_slph", "w") as out:
    SeqIO.write(positions, out, "fasta") 


