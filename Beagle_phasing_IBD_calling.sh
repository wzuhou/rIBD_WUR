#!/bin/bash
#SBATCH --time=2-00:10:00
#SBATCH --mail-user=zhou.wu@wur.nl
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --error=error.IBD%j
#SBATCH --job-name=vcftools
#SBATCH --qos=std
#SBATCH --mem=20000
#SBATCH --output=output.IBD%j
module load vcftools
module load bcftools
module load samtools
module load plink/160607
module load java
###
#Author:ZhouWu
#Usage: raw vcf file; step1 filtering based on few criteria;step2 change multi-snp name to seperate nm; step3 converting to plink format.
###
beagle="/home/WUR/wu090/Install/beagle/beagle.16May19.351.jar"
for chr in `cat candidate_chr` ;do \
#or chr=""
mkdir -p ${chr} && cd $_
vcf="/lustre/backup/WUR/ABGC/wu090/Backup/Bantam_gal6_vcf/Multi_vcf/gal6_CHR${chr}.rnm.vcf.gz"
ref_fa="/lustre/nobackup/WUR/ABGC/shared/public_data_store/genomes/chicken/Ensembl95/Gallus_gallus.GRCg6a.dna.toplevel.fa"
beagle="/home/WUR/wu090/Install/beagle/beagle.16May19.351.jar"
refine="/home/WUR/wu090/Install/beagle/refined-ibd.16May19.ad5.jar"

###-----------------------------###
### STEP1: REMEMBER TO MODIFY THE MISSING GENOTYPE NOTATION FROM . TO ./.
zcat ${vcf}| perl -pe "s/\s\.:/\t.\/.:/g"| bgzip -c >./gal6_CHR${chr}.rnm.beagle.vcf.gz

### STEP2: REMOVE THE LOW QUALITY SAMPLES AND GALLUS GALLUS SAMPLES >> 210
#vcftools --gzvcf gal6_CHR${chr}.rnm.beagle.vcf.gz --remove exclude.sample --recode --out gal6_CHR${chr}.rnm.beagle.210
#bgzip -c gal6_CHR${chr}.rnm.beagle.210.recode.vcf>gal6_CHR${chr}.rnm.beagle.210.vcf.gz

### STEP3: FILTER SNP SITE FOR COVERAGE AND DEPTH
vcf="gal6_CHR${chr}.rnm.beagle.vcf.gz"
#vcf="gal6_CHR4.rnm.beagle.210.vcf.gz"
#--remove exclude.sample
vcftools --gzvcf ${vcf} --min-meanDP 3 --max-meanDP 100 --minQ 10 --minDP 3 --maxDP  180 --recode --recode-INFO-all  --out ${vcf}.filtered
#bgzip -c ${vcf}.filtered.recode.vcf>${vcf}.filtered.recode.vcf.gz
#tabix -p vcf ${vcf}.filtered.recode.vcf.gz

### STEP4:MAKING MAP FILE
plink --vcf ${vcf}.filtered.recode.vcf --keep-allele-order --allow-extra-chr --allow-no-sex --cow --recode --out ${vcf}.filtered
k=`cat ../recombination.csv |awk -F, -v chr=${chr} '$1==chr {print $4}'`
cat ${vcf}.filtered.map|awk -F" " -v k=${k} '{$3=$4*k*0.000001 ; print}'> gal6_CHR${chr}.map

### STEP5:Phase genotype with Beagle v5
#java -Xmx32768m -jar ${beagle}  gt=${vcf}.filtered.recode.vcf  out=${vcf}.filtered.recode.phased_B5 map=gal6_CHR${chr}.map nthreads=12 window=0.02 overlap=0.01

### STEP6:REFINED IBD
refine="/home/WUR/wu090/Install/beagle/refined-ibd.16May19.ad5.jar"
#The "window" parameter must be at least 3 times the "ibd" parameter
java -Xmx36409m -jar ${refine} gt=${vcf}.filtered.recode.phased_B5.vcf.gz map=gal6_CHR${chr}.New.map out=LOD3_${vcf}.10k.filtered.recode.phased_B5.vcf.gz \
length=0.02 window=0.06 lod=3 nthreads=12 trim=0.001

cd ../..
done
date

### STEP7:CHANGE INDEL NAME
#tabix -p vcf ${vcf}.filtered.recode.phased_B5.vcf.gz
#bcftools norm -Ou -m -any ${vcf}.filtered.recode.phased_B5.vcf.gz|
#  bcftools norm -Ou -f ${ref_fa}|
#  bcftools annotate -Oz -x ID \
#  -I +'%CHROM:%POS:%REF:%ALT' >${vcf}.filtered.recode.phased_B5.vcf.gz.change_indel.vcf.gz
cd ..
done
