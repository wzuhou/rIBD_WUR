#!/bin/bash
#To use this script
#cd Example/
#usage: sh run_nIBD_pair.sh -i Example/test_chr6.ibd -A Example/List.DrFwB -B Example/List.DB -C Example/List.DrFw -O rIBD_DrFwB_DB_DrFw -W 20000 -S 10000

usage() { echo "Usage: $0 [-A <file>] [-B <file>] [-C <file>] [-i <ibd_input>] [-W <window in bp>] [-S <step in bp>] [-O <output prefix>]" 1>&2; exit 1; }

while getopts ":A:B:C:i:W:S:O:" o; do
    case "${o}" in
        A)
            A=${OPTARG}
            ;;
        B)
            B=${OPTARG}
            ;;

        C)
            C=${OPTARG}
            ;;
        i)
            i=${OPTARG}
            ;;
        W)
            W=${OPTARG}
            ;;
        S)
            S=${OPTARG}
            ;;

        O)
            O=${OPTARG}
            ;;

        *)
            usage
            ;;
    esac
done
shift $((OPTIND-1))

if [ -z "${i}" ] || [ -z "${A}" ] ||[ -z "${B}" ]|| [ -z "${C}" ] ; then
    usage
fi

if [ -z "${W}" ] ; then
S=$W #=step
fi

if [ -z "${O}" ] ; then
O="Output_prefix"
fi


echo "A = ${A}"
echo "B = ${B}"
echo "C = ${C}"
echo "i = ${i}"
echo "W = ${W}"
echo "S = ${S}"
echo "O = ${O}"


echo '
File ./CHR_coords.bed is needed for chromosome info, try to use .fai to make this file. e.g. CHR START END
awk "{print $1,0,$2}" *.fasta.fai > CHR_coords.bed
6       0       36365699
7       0       1640102
'

sh 0.nIBD_pair_pipeline_wind_v2.sh $A $B $i $W $S $O
sh 0.nIBD_pair_pipeline_wind_v2.sh $A $C $i $W $S $O
#Inside this 0.nIBD code
#1=list1
#2=list2 or list3
#3=ibd
#4=window
#5=step
#6=output
