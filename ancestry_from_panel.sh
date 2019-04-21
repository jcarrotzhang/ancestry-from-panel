# List of polymorphisms from 1000 genomes
BED=$1
# List of your bam files with full paths
bamlist=$2
# Reference genome (change it to your own path)
REF="/seq/references/Homo_sapiens_assembly19/v1/Homo_sapiens_assembly19.fasta"
# Output prefix
OUT=$2

# 1. Generate low-pass pileup files
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
