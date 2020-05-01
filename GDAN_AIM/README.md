# Local ancestry analyses used in the TCGA ancestry paper. 

## Usage example

#### local ancestry enrichment by gene. 

```
python TCGA_find_local_ancestry_by_gene.py PIK3CA 3:178864311-178959881
```

#### genome-wide local ancestry enricment (generate zscores). 

```
python TCGA_geno_scripts/TCGA_find_local_ancestry_count_blocks.py
python local_ancestry_enrichment.py
```

