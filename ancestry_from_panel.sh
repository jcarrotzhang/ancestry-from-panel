# This file must contain a list of your bam files (full paths)
bamlist=$1
# List of polymorphisms from 1000 genomes (from correct build)
BED="reference/ALL.wgs.phase3_shapeit2_mvncall_integrated_v5b.20130502.sites.bed"
# Reference genome (from correct build)
REF="/seq/references/Homo_sapiens_assembly19/v1/Homo_sapiens_assembly19.fasta"
# Output prefix
OUT=$2


# 1. Generate low-pass pileup files
# Call this once for each bam file
cat $bamlist | while read bam; do
         echo "samtools mpileup -q 20 -Q 20 -l $BED -f $REF $bam > $bam.pileup"
    done

# 2. Convert pileup to Laser input (.seq file)
cat $bamlist | while read bam; do
        echo "python ~bin/LASER-2.04/pileup2seq/pileup2seq.py -f $REF -m reference/HGDP.extract.site -o $bam $bam.pileup"
    done

# 3. prepare for LASER
cat $bam.seq >$OUT.seq
head -n 1 $bamlist | while read bam; do
        cp $bam.site $OUT.site
    done

#4. run LASER for both reference and cases
# Note:
# -- change SEQ_FILE to $OUT.seq
~bin/LASER-2.04/laser -p Laser.example.conf -k 10
