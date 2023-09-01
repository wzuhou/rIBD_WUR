
#!/bin/bash
#### USAGE
#### ./sh 0.nIBD_pair_pipeline_wind.sh breed1 breed2 ibd widow step output
#### Author: ZhouWu
#echo 'WARNING: Make sure bedtools is intalled and in your PATH
#Otherwise, ignore this message'
#echo "
#File ./CHR_coords.bed is needed for chromosome info, e.g. CHR START END
#6       0       36365699
#7       0       1640102
#"

## STEP1:CREATE THE IBD FILE OF TWO BREEDS
#for the breed related ibd files
B1=$1 #name of breed list
B2=$2 #name of breed2 list
IBD=$3 #test ibd=LOD3_gal6_Merged.rnm.beagle.vcf.gz.10k.filtered.recode.phased_B5.vcf.gz.ibd
breed1=$(basename ${B1})
breed2=$(basename ${B2})
ibd=$(basename ${IBD})
window=$4
step=$5
prefix=$6
#chr="Merg"
#ibd=LOD3_gal6_Merged.rnm.beagle.vcf.gz.10k.filtered.recode.phased_B5.vcf.gz.ibd
#################################################
n_breed1=`wc -l ${B1}|cut -f1 -d " "`
m_breed2=`wc -l ${B2}|cut -f1 -d " "`
#echo $n_breed1
#echo $m_breed2
mkdir -p ./ibd_breed
grep -f $B1 \
        ${IBD} \
        >./ibd_breed/${ibd}.$breed1

#---------------------------------------------------------------------------------
## STEP2 from the breed1 ibd file, extract the pairwise comparsion with the breed2
#echo "EXTRACT INDVIDUALS OF List.${breed2}"
mkdir -p ./pairwise/
grep -f ${B2} \
        ./ibd_breed/${ibd}.${breed1} \
        >./pairwise/ibdpair.${ibd}.${breed1}_${breed2}

if [[ -f ./pairwise/ibdpair.${ibd}.${breed1}_${breed2} ]] ; then
echo "
WROTE FILE TO ./pairwise/ibdpair.${ibd}.${breed1}_${breed2} "
fi
#---------------------------------------------------------------------------------
##STEP3 :bedtools
## make bed file out of ibd and calculate based
#module load bedtools
#prefix="pairwise/ibdpair.${ibd}.${breed1}_${breed2}" ##prefix= pairwise ibd file
prefix=${6}.${breed1}_${breed2}
## 3.1 Input the pairwise IBD file
sort -n -k 6 ./pairwise/ibdpair.${ibd}.${breed1}_${breed2} > $prefix.sort
awk '{print $5"\t"$6"\t"$7}' $prefix.sort >${prefix}.sort.bed
# make sliding win bed file
#if [[ ! -f ./pairwise/Chr_split]] ; then
bedtools makewindows -b CHR_coords.bed -w $window -s $step > ./pairwise/Chr_split
#fi
bedtools coverage -a ./pairwise/Chr_split \
        -b ${prefix}.sort.bed \
        >${prefix}.summary.wind.ibd

#---------------------------------------------------------------------------------
## STEP4: post  bedtools conversion
tIBD=$[$n_breed1*$m_breed2*2*2]
#echo  ${prefix} $tIBD
#awk '$2>=34080000 && $2<=34370000 {print $8=$4*$7/12}' summary.ibdpair.DB_FriFwB.sort.wind.ibd
awk -v  tIBD="$tIBD" '{print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$4/tIBD"\t"$4*$7/tIBD}' \
        ${prefix}.summary.wind.ibd \
        >${prefix}.nibd

#---------------------------------------------------------------------------------
ibd2size=`wc -l ${prefix}.nibd | cut -f 1 -d' ' `
if [ $ibd2size -gt 0 ] ; then
echo "CLEAN UP
Clean up unimportant intermediate files"
rm -r ibd_breed
rm ${prefix}.summary.wind.ibd  ${prefix}.sort.bed ${prefix}.sort pairwise/Chr_split
fi

echo "FINISH IBD SUMMARY
Final output: ${prefix}.nibd"

