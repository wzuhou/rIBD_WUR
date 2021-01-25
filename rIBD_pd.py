#!/usr/bin/env/ python3
"""
Author: Langqing Liu & Zhou Wu

usage: python rIBD_pd.py -i example/Tibetan-Zebu-Yak.Chr28.beagle.ibd -A example/Tibetan.list -B example/Zebu.list -C example/Yak.list -o Tibetan-Zebu-Yak.Chr28 -W 20000 -S 10000 -M Mirte (or Weighted)
"""
from __future__ import division
import sys, getopt
from sys import argv
import subprocess
import os.path
import pandas as pd
import pybedtools


def main(argv):
    inputfile = ''
    outputfile = ''
    A_pop_list = ''
    B_pop_list = ''
    C_pop_list = ''
    Window = ''
    Step = ''
    Method = ''
    try:
        opts, args = getopt.getopt(argv,"hi:o:A:B:C:W:S:M:",["ifile=","ofile=","A_pop_list=","B_pop_list=","C_pop_list=","windowsize=","stepsize=","method="])
    except getopt.GetoptError:
        print('USAGE:rIBD.py -i <inputfile> -o <outputfile> -A <A_pop_list> -B <B_pop_list> -C <C_pop_list> -W <windowsize> -S <stepsize> -M <method>')
        sys.exit(2)
    for opt, arg in opts:
        if opt == '-h':
            print('USAGE:rIBD.py -i <inputfile> -o <outputfile> -A <A_pop_list> -B <B_pop_list> -C <C_pop_list> -W <windowsize> -S <stepsize> -M <method>')
            sys.exit()
        elif opt in ("-i", "--ifile"):
            inputfile = arg
        elif opt in ("-o", "--ofile"):
            outputfile = arg
        elif opt in ("-A", "--A_pop_list"):
            A_pop_list = arg
        elif opt in ("-B", "--B_pop_list"):
            B_pop_list = arg
        elif opt in ("-C", "--C_pop_list"):
            O_pop_list = arg
        elif opt in ("-W", "--windowsize"):
            Window = arg
        elif opt in ("-S", "--stepsize"):
            Step = arg
        elif opt in ("-M", "--method"):
            Method = arg
    print('Input file is :', inputfile)
    print('Output file is :', outputfile)
    print('Method file is :', Method)
    return inputfile,outputfile,A_pop_list,B_pop_list,C_pop_list,Window,Step,Method

def read_ibd(inputfile,A_pop_list,B_pop_list,C_pop_list):
    A_pop=pd.read_csv(A_pop_list,sep='\t',header=None)
    A_pop_list=A_pop[0].values.tolist()
    B_pop=pd.read_csv(B_pop_list,sep='\t',header=None)
    B_pop_list=B_pop[0].values.tolist()
    C_pop=pd.read_csv(C_pop_list,sep='\t',header=None)
    C_pop_list=C_pop[0].values.tolist()
    ibd=pd.read_csv(inputfile,sep='\t',header=None)
    ibd_AB=ibd.loc[((ibd[0].isin(A_pop_list)) & (ibd[2].isin(B_pop_list)))|((ibd[2].isin(A_pop_list)) & (ibd[1].isin(B_pop_list))),[4,5,6]]
    ibd_AC=ibd.loc[((ibd[0].isin(A_pop_list)) & (ibd[2].isin(C_pop_list)))|((ibd[2].isin(A_pop_list)) & (ibd[1].isin(C_pop_list))),[4,5,6]]
    Chrsize = max(ibd[6])
    Chrname = ibd[4][1]
    nA = len(A_pop_list)
    nB = len(B_pop_list)
    nC = len(C_pop_list)
    ibd_AB = pybedtools.BedTool(ibd_AB.values.tolist())
    ibd_AC = pybedtools.BedTool(ibd_AC.values.tolist())
    return ibd_AB,ibd_AC,Chrname,Chrsize,nA,nB,nC

def genome_windows(Chrname,Chrsize,Window,Step):
    Start = 1
    End = int(Window)
    nWindow=int(Chrsize)//int(Step)
    #print(nWindow)
    chr_windows = []
    for i in range(1,nWindow):
        chr_windows.append((Chrname,Start,End))
        Start = int(Start) + int(Step)
        End = int(End) + int(Step)
    chr_windows = pybedtools.BedTool(chr_windows)
    chr_windows.saveas('w20ks10k.bed')
    return chr_windows

def bedtools_coverage(chr_windows,ibd_AB,ibd_AC):
    nibd_AB = chr_windows.coverage(ibd_AB)
    nibd_AO = chr_windows.coverage(ibd_AC)
    #print(nibd_AB.head(3))
    nibd_AB.saveas('AB.bed')
    nibd_AC.saveas('AC.bed')
    return nibd_AB,nibd_AC

def calculate_ribd(nibd_AB,nibd_AC,nA,nB,nC,outputfile):
    #print(len(nibd_AB))
    #print(len(nibd_AC))
    f = open(outputfile,"a")
    for i in range(0,len(nibd_AB)):
        Chr=str(nibd_AB[i][0])
        Start=str(nibd_AB[i][1])
        End=str(nibd_AB[i][2])
        sAB=float(nibd_AB[i][3])
        fAB=float(nibd_AB[i][6])
        sAC=float(nibd_AC[i][3])
        fAC=float(nibd_AC[i][6])
        ribd=(sAB/(nA*nB*4)) - (sAC/(nA*nC*4))
        w_ribd=(sAB*fAB/(nA*nB*4)) - (sAC*fAC/(nA*nC*4))
        #print("\t".join([Chr,Start,End,str(ribd),str(w_ribd)]))
        f.write("\t".join([Chr,Start,End,str(ribd),str(w_ribd)]))
    f.close()

def bedtools_genomecov(chr_windows,ibd_AB,ibd_AC,Chrname,Chrsize):
    ibd_AB_cov=ibd_AB.genome_coverage(bga=True, genome={Chrname:(1,Chrsize)})
    ibd_AC_cov=ibd_AC.genome_coverage(bga=True, genome={Chrname:(1,Chrsize)})
    ibd_AB_cov.saveas('AB.bed')
    ibd_AC_cov.saveas('AC.bed')
    ibd_AB_cov_pd=pybedtools.BedTool.to_dataframe(ibd_AB_cov)
    ibd_AC_cov_pd=pybedtools.BedTool.to_dataframe(ibd_AC_cov)
    sum_ibd=[]
    for i in chr_windows:
        AB=ibd_AB_cov_pd.loc[((ibd_AB_cov_pd['start']>int(i[1])) & (ibd_AB_cov_pd['start']<int(i[2])))|((ibd_AB_cov_pd['end']>int(i[1])) & (ibd_AB_cov_pd['end']<int(i[2])))]
        AC=ibd_AC_cov_pd.loc[((ibd_AC_cov_pd['start']>int(i[1])) & (ibd_AC_cov_pd['start']<int(i[2])))|((ibd_AC_cov_pd['end']>int(i[1])) & (ibd_AC_cov_pd['end']<int(i[2])))]
        #data["x1"]=data[["a","b"]].apply(lambda x:x["a"]+x["b"],axis=1)
        #print(AB)
        #print(AC)
        try:
            AB['w']=AB[['start','end','name']].apply(lambda x:x['name']*(x['end']-x['start']),axis=1)
            wAB=sum(AB['w'])/(max(AB['end'])-min(AB['start']))
        except:
            wAB=0
        try:
            AC['w']=AC[['start','end','name']].apply(lambda x:x['name']*(x['end']-x['start']),axis=1)
            wAC=sum(AC['w'])/(max(AC['end'])-min(AC['start']))
        except:
            wAC=0
        try:
            sum_ibd.append([i[0],i[1],i[2],wAB,wAC])
        except:
            sum_ibd.append([i[0],i[1],i[2],0,0])
    return sum_ibd

def calculate_ribd_gcov(sum_ibd,nA,nB,nC,outputfile):
    f = open(outputfile,"a")
    for i in sum_ibd:
        ribd=(i[3]/(nA*nB*4)) - (i[4]/(nA*nC*4))
        #print("\t".join([i[0],i[1],i[2],str(ribd)]))
        f.write("\t".join([i[0],i[1],i[2],str(ribd)]))
    f.close()
    return


if __name__ == "__main__":
    inputfile,outputfile,A_pop_list,B_pop_list,C_pop_list,Window,Step,Method = main(sys.argv[1:])
    #print(main(sys.argv[1:]))
    ibd_AB,ibd_AC,Chrname,Chrsize,nA,nB,nC = read_ibd(inputfile,A_pop_list,B_pop_list,C_pop_list)
    #print(nA)
    #print(nB)
    #print(nC)
    #print(Chrsize)
    #print(Window)
    #print(Step)
    chr_windows = genome_windows(Chrname,Chrsize,Window,Step)
    if Method == "Mirte":
        nibd_AB,nibd_AC=bedtools_coverage(chr_windows,ibd_AB,ibd_AC)
        calculate_ribd(nibd_AB,nibd_AC,nA,nB,nC,outputfile)
    else:
        sum_ibd=bedtools_genomecov(chr_windows,ibd_AB,ibd_AC,Chrname,Chrsize)
        calculate_ribd_gcov(sum_ibd,nA,nB,nC,outputfile)
