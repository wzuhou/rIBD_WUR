# rIBD_WUR

**Authors: Zhou Wu & Langqing Liu**
**Update: 2022/04/20**

## What does it do?
- Calcuate the relative haplotype sharing (rIBD) between three populations(A, B and C)

> Normalized IBD for one  group: nIBD=cIBD/tIBD (where cIBD=count of all haplotypes IBD between A and one group (B or C) and tIBD=total pairwise comparisons between A and one  group).
>
> rIBD between two pig groups: rIBD=nIBD<sub>AB</sub> – nIBD<sub>AC</sub>
>
> Adopted from *[Bosse, M. et al.(2014)](https://www.nature.com/articles/ncomms5392)*

- Codes in this repository

`Beagle_phasing_IBD_calling.sh` This code shows an example of using Beagle to perform the **phasing** and IBD calling

`rIBD_pd.py` This code shows an example of using the python code to compute the **rIBD**

`plot_rIBD.R` This code shows an example of **plotting** the outcome of the pyton codes 


## Citation

If you use codes from this repository, please cite the following papers:

Wu Z, Bortoluzzi C, Derks MFL, Liu L, Bosse M, Hiemstra SJ, Groenen MAM, Crooijmans RPMA. Heterogeneity of a dwarf phenotype in Dutch traditional chicken breeds revealed by genomic analyses.2020.Evol. Appl.eva.13183. [https://doi.org/10.1111/eva.13183](https://doi.org/10.1111/eva.13183)

Wu, Z., Bosse, M., Rochus, C.M. et al. Genomic insight into the influence of selection, crossbreeding, and geography on population structure in poultry. Genet Sel Evol 55, 5 (2023). [https://doi.org/10.1186/s12711-022-00775-x](https://doi.org/10.1186/s12711-022-00775-x)

## rIBD between population A,B, and C
<p align="center">
  <img src="https://github.com/wzuhou/rIBD_WUR/blob/main/Github_rIBD.png">
</p>

## How to use it? Prepare your input and env.
### Dependence: 
- python >= v3 (Tested in v3.8.5)
- python modules: pandas, pybedtools
- R 
- R package: ggplot2, gridExtra

### Example
```bash
python rIBD_pd.py -i Example/test_chr6.ibd -A Example/List.DrFwB -B Example/List.DB -C Example/List.DrFw -o rIBD_DrFwB_DB_DrFw -W 20000 -S 10000 -M 1
#expected output, see: Example/rIBD_DrFwB_DB_DrFw
```

Example breeds:DrFw, DrFwB, DB
- Input files:
1. A ibd file generated by Beagle: test_chr6.ibd
2. Three list files(A,B,C): A: List.DrFwB, B: List.DB , C: List.DrFw
3. For using the bash alternative, you also need a chromosome size file, see below (CHR_coords.bed)

- Parameters
1. Window size (-W): 20,000 bp; 
2. Step size (-S): 10,000 bp

- Method options

> Method `1` corresponds to *[Bosse, M. et al.(2014)](https://www.nature.com/articles/ncomms5392)*
> 
> Method `2` computed the weighted cIBD in each window 

- Output format
- 
1.Chr 2.Start 3.End 4.rIBD 5. weighted rIBD

## Bash alternative(untested but fast)

### Dependence: 
- bedtools (tested in v2.30.0)

**NB: this is an untested version to compute nIBD using bash scripts.**

**rIBD will need to be manually computated by simply rIBD=nIBD</sub>AB</sub> – nIBD<sub>AC</sub>**

**A chromosome size file (CHR_coords.bed) is required for this run**

You can try to use .fai to make this file. e.g. CHR START END

```sh
awk '{print $1,0,$2}' *.fasta.fai > CHR_coords.bed
6       0       36365699
7       0       1640102
```

### Example

Download scripts `run_nIBD_pair.sh` and `0.nIBD_pair_pipeline_wind_v2.sh` into same directory.
Similar usage and parameters as the python version.

```sh
sh run_nIBD_pair.sh -i Example/test_chr6.ibd -A Example/List.DrFwB -B Example/List.DB -C Example/List.DrFw -O rIBD_DrFwB_DB_DrFw -W 20000 -S 10000
#NB: note the difference in this script when giving output name, here is with uppercase -O ; And no -M function supported.
```

- Output

You will have 2 **nIBD** output files
1. rIBD_DrFwB_DB_DrFw.List.DrFwB_List.DB.nibd
2. rIBD_DrFwB_DB_DrFw.List.DrFwB_List.DrFw.nibd

To get rIBD simply follow this rIBD=nIBD</sub>AB</sub> – nIBD<sub>AC</sub>  
(We are considering to include this in the future, if you need help with this, do not hesitate to get in touch)

- Output format

"CHR","START","END","COUNT","BASE","WIN","COV_PERCENT","nIBD","WEIGHTED_nIBD"

---

## Visualization of the result
### Example to plot

Here we highlight the candidate region on Chromosome 6: 10000000-11000000

```bash
Rscript plot_rIBD.R rIBD_DrFwB_DB_DrFw 6 10000000 11000000
```

**Arguments:**
1. Input rIBD file (Example: rIBD_DrFwB_DB_DrFw. **NB**: Whole genome output or a single chromosome output? Both are possible.)
2. Chromosome of interest
3. Start position for region to be highlighted (in bp)
4. End position for region to be highlighted (in bp)

Output files is plot: Chr_RIBD.pdf   (*Optional*) + WG_RIBD.pdf (Only when the input contains more than two chromosomes) 
