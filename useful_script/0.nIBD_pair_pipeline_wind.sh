#!/bin/bash
#SBATCH --time=10:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=4
#SBATCH --error=error.count
#SBATCH --job-name=count
#SBATCH --partition=ABGC_Std
#SBATCH --mem=10000
#SBATCH --output=output.count
#### USAGE
#### ./sh 0.nIBD_pair_pipeline_wind.sh DB DB_cl.neoBantam
### ./sh 0.nIBD_pair_pipeline_wind.sh DB_cl.large DB_cl.neoBantam
#### Author: ZhouWu
mkdir -p ibd_breed
mkdir -p pairwise
## STEP1:CREATE THE IBD FILE OF TWO BREEDS
#for the breed related ibd files
breed1=$1
breed2=$2
chr="Merg"
#ibd="gal6_CHR1.ibd.merg"
#ibd="LOD2_gal6_CHR1.rnm.beagle.vcf.gz.10k.filtered.recode.phased_B5.vcf.gz.ibd.lod2"
#ibd="LOD1_gal6_CHR1.rnm.beagle.vcf.gz.10k.filtered.recode.phased_B5.vcf.gz.ibd.lod1"
ibd="LOD3_gal6_CHR1.rnm.beagle.vcf.gz.10k.filtered.recode.phased_B5.vcf.gz.ibd.lod3"
ibd="New_Map.LOD3_gal6_CHR1.rnm.beagle.vcf.gz.10k.filtered.recode.phased_B5.vcf.gz.ibd"
ibd="LOD3_gal6_Merged.rnm.beagle.vcf.gz.10k.filtered.recode.phased_B5.vcf.gz.ibd"
#ibd="LOD3_gal6_CHR1.rnm.beagle.vcf.gz.10k.filtered.recode.phased_B5.vcf.gz.hbd.lod3"
#ibd="gal6_CHR1.rnm.beagle.210.vcf.gz.filtered.recode.vcf.phased_B5.vcf.gz.vcf.gz.ibd.lod1"
#ibd="LOD2_Galgal5_traditional_adiv_cgn.fd.depth.rm.chr${chr}.phased.vcf.gz.ibd"
#gunzip -c 
#echo "OPEN IBD FILE OF ./ibd_breed/ibd.${breed1}"
##DB;DB_cl.large or DB_cl.neoBantam
#for the id of breed2 ind
n_breed1=`wc -l /lustre/nobackup/WUR/ABGC/wu090/Bantam/List_animal_212/List.${breed1}|cut -f1 -d " "`
m_breed2=`wc -l /lustre/nobackup/WUR/ABGC/wu090/Bantam/List_animal_212/List.${breed2}|cut -f1 -d " "`
#echo $n_breed1
#echo $m_breed2
mkdir -p ./ibd_breed/${chr}
grep -f /lustre/nobackup/WUR/ABGC/wu090/Bantam/List_animal_212/List.$breed1 \
	${ibd} \
	>./ibd_breed/${chr}/${ibd}.$breed1

#---------------------------------------------------------------------------------
## STEP2 from the breed1 ibd file, extract the pairwise comparsion with the breed2
#echo "EXTRACT INDVIDUALS OF List.${breed2}"
mkdir -p ./pairwise/${chr}
grep -f /lustre/nobackup/WUR/ABGC/wu090/Bantam/List_animal_212/List.${breed2} \
	./ibd_breed/${chr}/${ibd}.${breed1} \
	>./pairwise/${chr}/ibdpair.${ibd}.${breed1}_${breed2}

echo "WRITE FILE TO ./pairwise/${chr}/ibdpair.${breed1}_${breed2}"

#---------------------------------------------------------------------------------
##STEP3 :bedtools
## make bed file out of ibd and calculate based
module load bedtools
prefix="pairwise/${chr}/ibdpair.${ibd}.${breed1}_${breed2}" ##prefix= pairwise ibd file
## 3.1 Input the pairwise IBD file
sort -n -k 6 ${prefix} > $prefix.sort
awk '{print $5"\t"$6"\t"$7}' $prefix.sort >${prefix}.sort.bed
bedtools coverage -a ./pairwise/New_nIBD/split_chr/${chr}_split \
	-b ${prefix}.sort.bed \
	>${prefix}.summary.wind.ibd

#---------------------------------------------------------------------------------
## STEP4: post  bedtools conversion
tIBD=$[$n_breed1*$m_breed2*2*2]
echo  ${prefix} $tIBD
#awk '$2>=34080000 && $2<=34370000 {print $8=$4*$7/12}' summary.ibdpair.DB_FriFwB.sort.wind.ibd
awk -v  tIBD="$tIBD" '{print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$4/tIBD"\t"$4*$7/tIBD}' \
	${prefix}.summary.wind.ibd \
	>${prefix}.summary.sort.wind.ibd2

#---------------------------------------------------------------------------------
#input_file="/pairwise/ibdpair.${breed1}_${breed2}"
#sort -n -k 6 ${input_file}  > ${input_file}_sort
#awk '{print $7-$6}' ${input_file}_sort>${input_file}_size
#sh count_ibd.sh ./pairwise/ibdpair.${breed1}_${breed2}
echo "FINISH IBD SUMMARY"

