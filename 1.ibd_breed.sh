#! /bin/bash
breed=$1
grep -f /lustre/nobackup/WUR/ABGC/wu090/Bantam/List_animal_212/List.$breed Galgal5_traditional_adiv_cgn.fd.depth.rm.chr1.phased.vcf.gz.ibd >./ibd_breed/ibd.$breed

#for breed in `cat list_breedname`; 
#do grep -f /lustre/nobackup/WUR/ABGC/wu090/Bantam/List_animal_212/List.$breed Galgal5_traditional_adiv_cgn.fd.depth.rm.chr1.phased.vcf.gz.ibd >./ibd_breed/ibd.$breed;done

