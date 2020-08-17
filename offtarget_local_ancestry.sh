## low coverage calling local ancestry pipeline
## -----------------------------
## Input data ---
# This file must contain a list of your bam files (full paths)
bamlist=$1
# List of polymorphisms from 1000 genomes
BED="reference/ALL.wgs.phase3_shapeit2_mvncall_integrated_v5b.20130502.sites.bed"
# Reference genome (from correct build)
REF="/seq/references/Homo_sapiens_assembly19/v1/Homo_sapiens_assembly19.fasta"
# Output prefix
OUT=$2
## ---

# 1. Generate low-pass genotype calls
cat $bamlist | while read bam; do
       samtools mpileup -q 20 -Q 20 -l $BED -uf $REF $bam | bcftools call -mv -Oz > $bam.vcf.gz
       tabix $bam.vcf.gz
       done

# 2. Merge into a single vcf file per chromosome
for chr in `seq 1 22`; do
       bcftools merge -l $OUT.vcflist -r ${chr} -Oz > $OUT.chr${chr}.vcf.gz
done

bcftools merge -l $OUT.vcflist -r X -Oz > $OUT.chrX.vcf.gz

# 3. Imputation
# can replace the chunklist file defining regions from 1MB-5MB across the genome. 
cat reference/chunklist.txt | while read line; do
        chr=`echo $line | awk '{ print $1 }'`
        p0=`echo $line | awk '{ print $2 }'`
        p1=`echo $line | awk '{ print $3 }'`

        java -Xmx10g -jar ~bin/beagle.08Jun17.d8b.jar \
        gt=$OUT.chr${chr}.vcf.gz \
        ref=reference/beagle/chr${chr}.1kg.phase3.v5a.vcf.gz \
        out=$OUT.$chr.$p0.$p1.beagle_phased \
        map=reference/beagle/plink.chr${chr}.GRCh37.map \
        chrom=${chr}:${p0}-${p1} \
        impute=true

done;

# 4. Local ancestry identification
# 4.1 merge chunk for each chromosome (list beagle_phased.vcf.gz files generated from step 3 for each chromosome) in $OUT.chr${chr}.vcflist file.
# Select alleles with different population frequencies (alf_dif_snp.txt). 
for chr in {1..22}; do 
  bcftools concat -f $OUT.chr${chr}.vcflist | bcftools sort -Oz -o $OUT.chr${chr}.beagle_phased.vcf.gz; 
done

for chr in {1..22}; do 
  tabix $OUT.chr${chr}.beagle_phased.vcf.gz 
  bcftools view -R alf_dif_snp.txt $OUT.chr${chr}.beagle_phased.vcf.gz >$OUT.chr${chr}.beagle_phased.filtered.vcf
done

# 4.2 run RFMix for each chromosome. 
for chr in `seq 1 22`; do
  ~bin/rfmix-master/rfmix -f $OUT.chr${chr}.beagle_phased.filtered.vcf -r ~reference/beagle/chr${chr}.1kg.phase3.v5a.vcf.gz \
  -g ~reference/1000GP_Phase3/genetic_map_chr${chr}_combined_b37.rfmix.txt -m ref_samplelist.txt -o $OUT.chr${chr}.results --chromosome=${chr} \
done

## repeat step 4 for chrX. 
## ---
