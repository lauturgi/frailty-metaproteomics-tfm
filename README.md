# Identification of bacterial protein markers associated with frailty

## Metadata descriptive analysis
1. Create **data** and **plots** directories
2. Copy the following file inside **data**:
  - Tabla metadatos ENRICA.csv
3. Execute *metadata_descriptive_analysis.R*

## Protein-centric analysis
1. Create **data**, **output** and **plots** directories in **./prdea/**
2. Copy the following files inside **data**:
 - combined_protein.tsv
 - experiment_annotation.tsv
3. Execute R scripts in the following order:
 - *01_preprocessing_lfq.R*
 - *02_limma_prde_analysis.R*

## Peptide-centric analysis

1. Match UHGP proteins to UniRef90 (in Masterbio3)
- Download uniref90.fasta
```bash
wget https://ftp.uniprot.org/pub/databases/uniprot/uniref/uniref90/uniref90.fasta.gz
```
- Download UHGP protein catalogue
```bash
wget https://ftp.ebi.ac.uk/pub/databases/metagenomics/mgnify_genomes/human-gut/v2.0.2/protein_catalogue/uhgp-90.tar.gz
```
- Unzip UniRef90 and UHGP files
```bash
gunzip uniref90.fasta.gz
tar -xvzf uhgp-90.tar.gz
```
- Download diamond-aligner
```bash
wget https://github.com/bbuchfink/diamond/releases/download/v2.1.11/diamond-linux64.tar.gz
tar -xvzf diamond-linux64.tar.gz
```
- Build DIAMOND database from FASTA file
```bash
./diamond makedb --in uniref90.fasta -d uniref90
```
- Align amino acid query sequences from UHGP-90 against UniRef90 protein reference database
```bash
nohup ./diamond blastp --db uniref90.dmnd -q uhgp-90/uhgp-90.faa -o diamond_uhgp_out.tsv --sensitive -e 0.1 --top 5 -f 6 qseqid qlen sseqid slen evalue length nident > diamond.log 2>&1 &
```
- Keep only best hits
```bash
sort -k1,1 -k5,5g diamond_uhgp_out.tsv | awk '!seen[$1]++' > diamond_uhgp_best_hits_out.tsv
```

2. In silico tryptic protein digestion
- Download trypsin script and install required libraries
```bash
wget https://github.com/northomics/bin/blob/master/trypsin.py
pip install biopython
pip install tqdm
```
- Create a file of peptides with their corresponding protein: 
-i: fasta file of proteins to digest
-c: missed trypsin cleavages
-m: minimum length of analyzed tryptic peptides
-a: maximum length of analyzed tryptic peptides
-o: fasta file of peptide sequences
```bash
nohup python3 trypsin.py -i uhgp-90/uhgp-90.faa -c 2 -m 5 -a 50 -o uhgp_pep.tsv > output.log 2>&1 &
```

3. Create **data**, **output** and **plots** directories in **./psva/**
4. Copy the following files inside **data**:
  - combined_peptide.tsv 
  - experiment_annotation.tsv
  - ko00001.keg
  - uhgp-90_eggNOG.tsv
  - diamond_uhgp_best_hits_out.tsv
  - uhgp_pep.tsv

5. Execute R scripts in the following order:
  - *01_get_kegg_orthology_l3_4.R*
  - *02_create_peptide_kegg_db_2.R* (in Masterbio3)
  - *03_create_peptide_ko_count_db_2.R* (in Masterbio3)
  - *04_create_kegg_peptidesets.R* (in Masterbio3)
  - *05_preprocessing_peptide_lfq.R*
  - *06_match_core_pep_kegg.R* (in Masterbio3)
  - *07_peptide_set_variation_analysis.R* (GSVA and extract peptides from significant pathways in Masterbio3)
