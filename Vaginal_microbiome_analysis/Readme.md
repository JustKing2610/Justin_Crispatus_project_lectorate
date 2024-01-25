# Analysis of the vaginal microbiome using oxford nanopore sequencing
In this folder, you will find the results of each run as analysed by EPI2ME labs. Within the folders, you can find the final HTML summary, the report html file containing the command / parameters and the versions of used tools used.

Data used is stored on the Avans university of applied sciences servers

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
nextflow run epi2me-labs/wf-metagenomics --fastq /file/path/to/raw/input/data -with-docker 'epi2me/image:tag' -work-dir /path/to/desired/working/directory/ [other settings]
```
