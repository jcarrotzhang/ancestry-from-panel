# ancestry-from-panel
Ancestry identification from panel sequencing data using off-target reads
## Requirements
* Samtools
* LASER
* BEAGLE
* RFMIX (https://github.com/slowkoni/rfmix)

## Usage 

#### Run LASER for ancestry identification. 
A BED file of 1000 genome SNP sites and a list of BAM files with full path are required. 

```
sh ancestry_from_panel.sh SNP_sites.bed Bamfiles_list.txt output_prefix
```
* HGDP.extract.site for LASER is provided in reference folder. Please prepare the HGDP.extract.geno file with SNPs matching the .site file.


#### Run local ancestry identification.
Please prepare a BED file of selected SNPs to assign ancestry for.

```
sh offtarget_local_ancestry.sh Bamfiles_list.txt output_prefix
```

* BEAGLE reference vcf : http://bochet.gcc.biostat.washington.edu/beagle/genetic_maps/ 
* BEAGLE reference map : http://bochet.gcc.biostat.washington.edu/beagle/1000_Genomes_phase3_v5a/ 

#### Run local ancestry risk score.
This is a simple python script that is an example of how I calculated the cross-validated local ancstry risk score. (k=10 for ten-fold cross-validation)
```
python local_ancestry_risk_score_kfold.py k
```
