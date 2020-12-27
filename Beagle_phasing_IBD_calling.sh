#!/bin/bash
#SBATCH --time=5-3:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=2
#SBATCH --error=error_ibd%j
#SBATCH --job-name=vcftools
#SBATCH --partition=ABGC_Std
#SBATCH --mem=40000
#SBATCH --output=output_ibd%j
module load java
date
#Beagle 5.0
#for chr in 4;do \
# java -Xmx40g -jar /home/WUR/wu090/Install/refined-ibd.26Feb19.29e.jar \
# gt=Galgal5_traditional_adiv_cgn.fd.depth.rm.chr4.phased.vcf.gz.vcf.gz \
# out=Refined_Galgal5_traditional_adiv_cgn.fd.depth.rm.chr4.phased.vcf
#window=0.02
#nthreads=8 \
#length=0.0002 \
#lod=3;
#done
# Explaination
#If no genetic map is specified, Refined IBD will assume a constant recombination rate of 1 cM per Mb
#window=[positive number] specifies the cM length of the sliding marker window (default: window=40.0cM) The overlap between adjacent windows will be twice the length of the ibd parameter(?).
#length=[positive number] specifies the minimum cM length for reported IBD segments (default: length=1.5)
#----------------------------------
###Phase genotype with Beagle v4.1
## input gt= [unphased vcf]
indir="/lustre/nobackup/WUR/ABGC/wu090/Bantam/adiv_cgn_fli/data_origin"
#Galgal5_traditional_adiv_cgn.fd.depth.rm.chr1.vcf.gz
for chr in 1;do \
        java -Xmx40g -jar /home/WUR/wu090/Install/beagle.27Jan16.aae.jar \
        gt=Galgal5_traditional_adiv_cgn.fd.depth.rm.chr${chr}.vcf.gz \
        out=LOD2_Galgal5_traditional_adiv_cgn.fd.depth.rm.chr${chr}.phased.vcf.gz \
        window=20000 \
        overlap=1000 \
        nthreads=16 \
        niterations=10 \
        ne=10000 \
        ibd=true \
        ibdlod=2.0
done
date
