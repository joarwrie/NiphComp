# NiphComp
Pipeline used to compare results from groundwater eDNA metabarcoding to DNA barcoding on individual amphipod collected via a citizen science campaign. It analyses a COI library from groundwater eDNA amplified with primers modified for Niphargus amplification

This pipeline was used in the article “Integrating citizen science and environmental DNA metabarcoding to study biodiversity of groundwater amphipods in Switzerland.” by Marjorie Couton, Angela Studer, Samuel Hürlemann, Nadine Locher, Mara Knüsel, Roman Alther & Florian Altermatt, submitted for publication in Scientific Reports.

The pipeline can be performed by running the file named "file_to_run.sh". All softwares or R packages and their version number are indicated as comments.

An additional step of phylogenetic placement was performed using the software DARN. As this tool is presented as a docker, it was not included in the script but the results are provided in the file named DARN_res_per_query.tsv.

Raw sequencing data on which this pipeline was used are stored as SRA in ncbi under the BioProject ID: PRJNA905821.
