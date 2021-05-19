# Local ancestry analyses used in the TCGA ancestry paper. 

## Usage example

#### Local ancestry enrichment by gene. 

```
python find_local_ancestry_by_gene.py PIK3CA 3:178864311-178959881
```

#### Genome-wide local ancestry enricment (generate zscores). 

#### Step 1:
```
python find_local_ancestry_count_blocks.py samplelist.txt
```
#####    Input: samplelist.txt   (prepare a list of local ancestry output file per sample in bed format.)
```
TCGA-35-A46O_A.bed
TCGA-35-A46O_B.bed
TCGA-36-A47O_A.bed
TCGA-36-A47O_B.bed
```
#### Step 2:
```
python local_ancestry_enrichment.py
```

