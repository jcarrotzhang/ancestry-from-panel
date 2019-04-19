# ancestry-from-panel
Ancestry identification from panel sequencing data using off-target reads
## Requirements
* Samtools
* LASER
## Usage 
A BED file of 1000 genome SNP sites and a list of BAM files with full path are required. 

```
sh ancestry_from_panel.sh SNP_sites.bed Bamfiles_list.txt output_prefix
```

** ** HGDP.extract.site for LASER is provided. Please prepare the HGDP.extract.geno file with SNPs matching the .site file.
