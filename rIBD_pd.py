#!/usr/bin/env/ python3
"""
Author: Langqing Liu & Zhou Wu

usage: python rIBD.py -i example/Tibetan-Zebu-Yak.Chr28.beagle.ibd -A example/Tibetan.list -B example/Zebu.list -O example/Yak.list -o Tibetan-Zebu-Yak.Chr28 -W 20000 -S 10000
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
    O_pop_list = ''
    Window = ''
    Step = ''
    try:
        opts, args = getopt.getopt(argv,"hi:o:A:B:O:W:S:",["ifile=","ofile=","A_pop_list=","B_pop_list=","O_pop_list=","windowsize=","stepsize="])
    except getopt.GetoptError:
        print('USAGE:rIBD.py -i <inputfile> -o <outputfile> -A <A_pop_list> -B <B_pop_list> -O <outgroup_pop_list> -W <windowsize> -S <stepsize>')
        sys.exit(2)
    for opt, arg in opts:
        if opt == '-h':
            print('USAGE:rIBD.py -i <inputfile> -o <outputfile> -A <A_pop_list> -B <B_pop_list> -O <outgroup_pop_list> -W <windowsize> -S <stepsize>')
            sys.exit()
        elif opt in ("-i", "--ifile"):
            inputfile = arg
        elif opt in ("-o", "--ofile"):
            outputfile = arg
        elif opt in ("-A", "--A_pop_list"):
            A_pop_list = arg
        elif opt in ("-B", "--B_pop_list"):
            B_pop_list = arg
        elif opt in ("-O", "--O_pop_list"):
            O_pop_list = arg
        elif opt in ("-W", "--windowsize"):
            Window = arg
        elif opt in ("-S", "--stepsize"):
            Step = arg
    print('Input file is :', inputfile)
    print('Output file is :', outputfile)
    return inputfile,outputfile,A_pop_list,B_pop_list,O_pop_list,Window,Step

def read_ibd(inputfile,A_pop_list,B_pop_list,O_pop_list):
    A_pop=pd.read_csv(A_pop_list,sep='\t',header=None)
    A_pop_list=A_pop[0].values.tolist()
    B_pop=pd.read_csv(B_pop_list,sep='\t',header=None)
    B_pop_list=B_pop[0].values.tolist()
    O_pop=pd.read_csv(O_pop_list,sep='\t',header=None)
    O_pop_list=O_pop[0].values.tolist()
    ibd=pd.read_csv(inputfile,sep='\t',header=None)
    ibd_AB=ibd.loc[((ibd[0].isin(A_pop_list)) & (ibd[2].isin(B_pop_list)))|((ibd[2].isin(A_pop_list)) & (ibd[1].isin(B_pop_list))),[4,5,6]]
    ibd_AO=ibd.loc[((ibd[0].isin(A_pop_list)) & (ibd[2].isin(O_pop_list)))|((ibd[2].isin(A_pop_list)) & (ibd[1].isin(O_pop_list))),[4,5,6]]
    Chrsize = max(ibd[6])
    Chrname = ibd[4][1]
    nA = len(A_pop_list)
    nB = len(B_pop_list)
    nO = len(O_pop_list)
    ibd_AB = pybedtools.BedTool(ibd_AB.values.tolist())
    ibd_AO = pybedtools.BedTool(ibd_AO.values.tolist())
    return ibd_AB,ibd_AO,Chrname,Chrsize,nA,nB,nO

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

def bedtools_coverage(chr_windows,ibd_AB,ibd_AO):
    nibd_AB = chr_windows.coverage(ibd_AB)
    nibd_AO = chr_windows.coverage(ibd_AO)
    #print(nibd_AB.head(3))
    nibd_AB.saveas('AB.bed')
    nibd_AO.saveas('AO.bed')
    return nibd_AB,nibd_AO

def calculate_ribd(nibd_AB,nibd_AO,nA,nB,nO):
    #print(len(nibd_AB))
    #print(len(nibd_AO))
    for i in range(0,len(nibd_AB)):
        Chr=str(nibd_AB[i][0])
        Start=str(nibd_AB[i][1])
        End=str(nibd_AB[i][2])
        sAB=float(nibd_AB[i][3])
        fAB=float(nibd_AB[i][6])
        sAO=float(nibd_AO[i][3])
        fAO=float(nibd_AO[i][6])
        ribd=(sAB/(nA*nB*4)) - (sAO/(nA*nO*4))
        w_ribd=(sAB*fAB/(nA*nB*4)) - (sAO*fAO/(nA*nO*4))
        print("\t".join([Chr,Start,End,str(ribd),str(w_ribd)]))

def bedtools_genomecov(chr_windows,ibd_AB,ibd_AO,Chrname,Chrsize):
    ibd_AB_cov=ibd_AB.genome_coverage(bga=True, genome={Chrname:(1,Chrsize)})
    ibd_AO_cov=ibd_AO.genome_coverage(bga=True, genome={Chrname:(1,Chrsize)})
    ibd_AB_cov_pd=pybedtools.BedTool.to_dataframe(ibd_AB_cov)
    ibd_AO_cov_pd=pybedtools.BedTool.to_dataframe(ibd_AO_cov)
    sum_ibd=[]
    for i in chr_windows:
        AB=ibd_AB_cov_pd.loc[((ibd_AB_cov_pd['start']>int(i[1])) & (ibd_AB_cov_pd['start']<int(i[2])))|((ibd_AB_cov_pd['end']>int(i[1])) & (ibd_AB_cov_pd['end']<int(i[2])))]
        AO=ibd_AO_cov_pd.loc[((ibd_AO_cov_pd['start']>int(i[1])) & (ibd_AO_cov_pd['start']<int(i[2])))|((ibd_AO_cov_pd['end']>int(i[1])) & (ibd_AO_cov_pd['end']<int(i[2])))]
        #data["x1"]=data[["a","b"]].apply(lambda x:x["a"]+x["b"],axis=1)
        try:
            AB['w']=AB[['start','end','name']].apply(lambda x:x['name']*(x['end']-x['start']),axis=1)
            AO['w']=AO[['start','end','name']].apply(lambda x:x['name']*(x['end']-x['start']),axis=1)
            wAB=sum(AB['w'])/(max(AB['end'])-min(AB['start']))
            wAO=sum(AO['w'])/(max(AO['end'])-min(AO['start']))
            sum_ibd.append([i[0],i[1],i[2],wAB,wAO])
        except:
            pass
    return sum_ibd

def calculate_ribd_gcov(sum_ibd,nA,nB,nO):
    for i in sum_ibd:
        ribd=(i[3]/(nA*nB*4)) - (i[4]/(nA*nO*4))
        print("\t".join([i[0],i[1],i[2],str(ribd)]))
    return


if __name__ == "__main__":
    inputfile,outputfile,A_pop_list,B_pop_list,O_pop_list,Window,Step = main(sys.argv[1:])
    #print(main(sys.argv[1:]))
    ibd_AB,ibd_AO,Chrname,Chrsize,nA,nB,nO = read_ibd(inputfile,A_pop_list,B_pop_list,O_pop_list)
    #print(nA)
    #print(nB)
    #print(nO)
    #print(Chrsize)
    #print(Window)
    #print(Step)
    chr_windows = genome_windows(Chrname,Chrsize,Window,Step)
    #nibd_AB,nibd_AO=bedtools_coverage(chr_windows,ibd_AB,ibd_AO)
    #calculate_ribd(nibd_AB,nibd_AO,nA,nB,nO)

    sum_ibd=bedtools_genomecov(chr_windows,ibd_AB,ibd_AO,Chrname,Chrsize)
    calculate_ribd_gcov(sum_ibd,nA,nB,nO)
