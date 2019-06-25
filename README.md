# ancestry-from-panel
Ancestry identification from panel sequencing data using off-target reads
## Requirements
* Samtools
* LASER
## Usage 

#### Run LASER for ancestry identification. 
A BED file of 1000 genome SNP sites and a list of BAM files with full path are required. 

```
sh ancestry_from_panel.sh SNP_sites.bed Bamfiles_list.txt output_prefix
```
** ** HGDP.extract.site for LASER is provided. Please prepare the HGDP.extract.geno file with SNPs matching the .site file.


#### Prepare RFMIX input.
```
do for i in {1..22}; do python preRFMIX_imputation.py 
\ ${prefix}.chr${i}.beagle_phased.filtered.vcf 
\ reference/chr${i}.1kg.phase3.v5a.vcf.gz 
\ reference/genetic_map_chr${i}_combined_b37.txt 
\ ${prefix}.chr${i}.alleles 
\ ${prefix}.chr${i} 
\ ${prefix}.chr${i}.classes 
\ ${i}"; done
```
