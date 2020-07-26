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
** ** HGDP.extract.site for LASER is provided. Please prepare the HGDP.extract.geno file with SNPs matching the .site file.


#### Run local ancestry identification.
```
sh offtarget_local_ancestry.sh Bamfiles_list.txt output_prefix
```

** ** BEAGLE reference vcf : http://bochet.gcc.biostat.washington.edu/beagle/genetic_maps/ ** **
** ** BEAGLE reference map : http://bochet.gcc.biostat.washington.edu/beagle/1000_Genomes_phase3_v5a/ ** **

