#! /bin/bash
#for the breed related ibd files
breed1=$1
echo "OPEN IBD FILE OF ./ibd_breed/ibd.${breed1}"
#for the id of breed2 ind
breed2=$2
echo "EXTRACT INDVIDUALS OF ./List_animal_ID/List.${breed2}"
#from the breed1 ibd file, extract the pairwise comparsion with the breed2
grep -f ./List_animal_ID/List.${breed2} ./ibd_breed/ibd.${breed1}>./pairwise/ibdpair.${breed1}_${breed2}
echo "WRITE FILE TO ./pairwise/ibdpair.${breed1}_${breed2}"
