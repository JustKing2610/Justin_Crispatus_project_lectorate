

# Goals
- Build cgMLST schema of *Lactobacillus crispatus*
- Evaluate schema


# Notes:


- There is no prodigal training file yet for L crispatus within ChewBBACA, therefore I'll make it from the reference genome on [ncbi](https://www.ncbi.nlm.nih.gov/nuccore/CP039266.1?report=fasta):
```bash
prodigal -i genomes/CP039266.fasta -t crispatus.trn -p single
```



- Starting point is the 700 genomes from NCBI. I ran quast & busco and selected the 30 most complete & contiguous assemblies to base the initial scheme on.

```bash
genomes/
├── allelecall
│   ├── 670x genome.fna
└── schema
    ├── 30x genome.fna
```


Creating the schema:
```bash
chewBBACA.py CreateSchema -i genomes/schema/ -o schema --cpu 20
# populate schema with alleles from schema genomes
chewBBACA.py AlleleCall -i genomes/schema/ -g schema/schema_seed/ -o schema_wgmlst --cpu 20
# remove paralogous genes
chewBBACA.py RemoveGenes -i schema_wgmlst/results_alleles.tsv -g schema_wgmlst/paralogous_counts.tsv -o schema_wgmlst/results_alleles_NoParalogs.tsv
```

Determination of cgMLST schema.
```bash
chewBBACA.py ExtractCgMLST -i schema_wgmlst/results_alleles_NoParalogs.tsv -o schema_wgmlst/cgMLST
```
![[Pasted image 20231206160317.png]]

A little under 1200 loci will be kept for the 95% cgMLST. For 99 and 100%, this will be a little less than 1000.


- Allele calling the 670 genomes
```bash
chewBBACA.py AlleleCall -i genomes/allelecall/ -g schema/schema_seed/ --gl schema_wgmlst/cgMLST/cgMLSTschema95.txt -o schema_cgMLST --cpu 20
# concatenate allelic profile with extracted cgMLST
chewBBACA.py JoinProfiles -p schema_wgmlst/cgMLST/cgMLST95.tsv schema_cgMLST/results_alleles.tsv -o cgMLST.tsv
# redetermine cgMLST scheme
chewBBACA.py ExtractCgMLST -i cgMLST.tsv -o final_cgMLST
```

![[Pasted image 20231206162927.png]]


* I filtered out low quality genomes from the cgMLST based on the criteria in the tutorial:
![[Pasted image 20231206162804.png]]

```bash
chewBBACA.py ExtractCgMLST -i cgMLST.tsv -o cgMLST_qualityFiltered --g to_exclude.txt
```


- Lastly, I annotated the schema loci and generated a schema report
```bash
chewBBACA.py UniprotFinder -i schema/schema_seed/ -o uniprotfinder -t schema_cgMLST/cds_coordinates.tsv --taxa "Lactobacillus crispatus" --cpu 20

chewBBACA.py SchemaEvaluator -g schema/schema_seed/ -a uniprotfinder/schema_seed_annotations.tsv --cpu 32 --loci-reports -o schema_eval_full
```
