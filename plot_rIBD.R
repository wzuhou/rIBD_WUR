#!/bin/Rscript
###Author: wu090(2020)
library(ggplot2)
require(gridExtra)
#library(data.table) #Optional,to speedup the read in file process
args <- commandArgs(T)
setwd(system("pwd",intern=T))
## arguments:
## [1] Input rIBD file (example: rIBD_DrFwB_DB_DrFw)
## [2] Chromosome of interest
## [3] Start position for candidate region(in bp)
## [4] End position for candidate region(in bp)
## Output files are plots: Chr_RIBD.pdf+ WG_RIBD.pdf (Only when the input contains more than two chromosomes) 

#############
##Plot rIBD #
#############
rIBD<- read.table(args[1],sep="\t",header = F)
if (ncol(rIBD) == 5){
	names(rIBD) <- c("CHR","START","END","R_IBD","WR_IBD")
	method=1
} else {
	names(rIBD) <- c("CHR","START","END","WR_IBD")
	method=2
}
#########################
#Zoom in Chr of interest#
#########################
Chr=args[2]
sline=as.numeric(args[3])/1000000
eline=as.numeric(args[4])/1000000
rIBD <- subset(rIBD,CHR==Chr)
if (method==1){
gg2 <- ggplot(rIBD ) +  
  geom_col(size=5,fill="black",aes(x=START/1000000, y=R_IBD))+
  scale_x_continuous(breaks =waiver(),expand = c(0, 0))+
  #scale_y_continuous(breaks =waiver(),limits = c(-1,1),expand = c(0, 0))+
  labs(title= NULL,x=paste0("Position on Chr ",Chr," (Mb)"),y="rIBD")+
  geom_vline(xintercept =c(sline,eline),colour="#990000",linetype='dashed',alpha=0.7)+
  theme_light()+theme(axis.title.y = element_text(size=10))
ggsave(paste0("Chr",Chr,"_RIBD.pdf"),gg2,width = 5,height = 2.5)
# plot a second figure for weighted rIBD
gg3 <- ggplot(rIBD ) +  
  geom_col(size=5,fill="black",aes(x=START/1000000, y=WR_IBD))+
  scale_x_continuous(breaks =waiver(),expand = c(0, 0))+
  #scale_y_continuous(breaks =waiver(),limits = c(-1,1),expand = c(0, 0))+
  labs(title= NULL,x=paste0("Position on Chr ",Chr," (Mb)"),y="rIBD")+
  geom_vline(xintercept =c(sline,eline),colour="#990000",linetype='dashed',alpha=0.7)+
  theme_light()+theme(axis.title.y = element_text(size=10))
ggsave(paste0("Chr",Chr,"_WRIBD.pdf"),gg3,width = 5,height = 2.5)
} else {
  gg3 <- ggplot(rIBD ) +
  geom_col(size=5,fill="black",aes(x=START/1000000, y=WR_IBD))+
  scale_x_continuous(breaks =waiver(),expand = c(0, 0))+
  #scale_y_continuous(breaks =waiver(),limits = c(-1,1),expand = c(0, 0))+
  labs(title= NULL,x=paste0("Position on Chr ",Chr," (Mb)"),y="rIBD")+
  geom_vline(xintercept =c(sline,eline),colour="#990000",linetype='dashed',alpha=0.7)+
  theme_light()+theme(axis.title.y = element_text(size=10))
ggsave(paste0("Chr",Chr,"_WSRIBD.pdf"),gg3,width = 5,height = 2.5)

}
########################
# whole genome overview#
########################
if (length(unique(rIBD$CHR))>1){
	if (method==1){
		gg1 <- ggplot(rIBD, aes(x=START/1000000, y=R_IBD) ) + 
		geom_col(size=6,fill="black",alpha=1)+
		scale_x_continuous(breaks =waiver(),expand = c(0, 0))+
		#scale_y_continuous(breaks =c(0),expand = c(0, 0))+
		theme_light()+
		facet_wrap(~ CHR,ncol=7,scales="free_x",drop=T)+ #wrap ncol=7 or grid,space="free_x"
		labs(title= NULL,x=paste0("Position on Chr","(Mb)"),y="rIBD")
		ggsave("WG_RIBD.pdf",gg1,width = 15,height = 7.5)
	}
	else{
		gg1 <- ggplot(rIBD, aes(x=START/1000000, y=WR_IBD) ) + 
		geom_col(size=6,fill="black",alpha=1)+
		scale_x_continuous(breaks =waiver(),expand = c(0, 0))+
		#scale_y_continuous(breaks =c(0),expand = c(0, 0))+
		theme_light()+
		facet_wrap(~ CHR,ncol=7,scales="free_x",drop=T)+ #wrap ncol=7 or grid,space="free_x"
		labs(title= NULL,x=paste0("Position on Chr","(Mb)"),y="rIBD")
		ggsave("WG_RIBD.pdf",gg1,width = 15,height = 7.5)
	}

}
