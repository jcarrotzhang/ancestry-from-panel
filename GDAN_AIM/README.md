# Local ancestry analyses used in the TCGA ancestry paper. 

## Usage example

#### Local ancestry enrichment by gene. 

```
python find_local_ancestry_by_gene.py PIK3CA 3:178864311-178959881
```

#### Genome-wide local ancestry enricment (generate zscores). 

#### Step 1: read in local ancestry breakpoints identified in the cohort.
```
python local_ancestry_count_blocks.py samplelist.txt
```
#####    Input: samplelist.txt   (prepare a list of local ancestry output file per sample in bed format.)
```
TCGA-D8-A1JP_A.bed
TCGA-D8-A1JP_B.bed
TCGA-E2-A1L9_A.bed
TCGA-E2-A1L9_B.bed
```
#### Step 2:
```
python local_ancestry_enrichment.py
```

