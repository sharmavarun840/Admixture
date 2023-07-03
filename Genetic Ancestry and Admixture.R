#for admixture run#
#same binary files used for PCA are used for finding the admixture proportions##
#three type of files are generated .P, .Q and .log files##
#log files are used for finding out the optimum CV Value for admixture components##

## making a bash file to run admixture###

#!/bin/bash
#SBATCH -c 128
#SBATCH -a 2-25
#SBATCH -t 7-00:00:00

##this batch file is prepared if using server for running admixture##
#NAME THIS FILE WITH EXTENSION filename.sh
source `/.bashrc
conda activate myenv

admixture -s time --cv=5 input.bed ${SLURM_ARRAY_TASKID} -j128 | tee ADMIXTURE_log{SLURM_ARRAY_TASK_ID}.out

#use command###
sbatch filename.sh

#it will tak days to b done# after it is completed# plotting can be done#


####for CV error plot####
# make the cv error file in SHELL##
grep "CV error" ADMIXTURE_log{2..25}.out >> cv_error.comp.txt
sed -e 's/[)(]//g' cv_error.comp.txt > cv_err.txt
cat cv_err.txt | tr -d ":" > Full_ADMIXTURE.cv_error.txt
rm cv_err*.txt

###Plot cv error plot in R####
cv = read.table("Full_ADMIXTURE.cv_error.txt", header = FALSE)
head(cv)
cv$K=cv$V3
cv$K=gsub("K=","",cv$K)
head(cv)
cv$K=as.numeric(cv$K)

no_runs=1  #how many times was ADMIXTURE repeated (J value)
no_K=24    #how many Ks were there? the value of K in ADMIXTURE minus 1 because K   starts from 2
runs=c()
runs2=c()
for (i in 1:no_runs){runs[i]=paste("Run",i,sep = "")}
for (i in 1:no_runs){x=(rep(runs[i],no_K));runs2=c(runs2,x)}
runs2  
cv$V1=runs2
head(cv)
cv$V3=as.factor(cv$V3)

color = grDevices::colors()[grep('gr(a|e)y', grDevices::colors(), invert = T)]
kolor = sample(color,no_K)

pdf("M1T.Ref_allHim_CVerrorplot.pdf", width=6, height=3)
library(ggplot2)
ggplot(cv, aes(x=K,y=V4)) + geom_point(alpha=0.5,color="black") +  
  xlab("") + ylab("CV_error") + geom_line(color="gray27", linetype = "dashed", alpha=0.5) +
  scale_x_continuous(breaks = c(1:25)) +
  scale_color_manual(values=kolor) + theme_bw() + 
  theme(legend.position = "none")
dev.off()


###### for barplotsof admixture####3



library(pophelper)
library(tidyr)
library(dplyr)

#check version
packageDescription("pophelper",fields="Version")

#read q files  
afiles=c("SIG1.2.Q",
         "SIG1.3.Q",
         "SIG1.4.Q",
         "SIG1.5.Q",
         "SIG1.6.Q",
         "SIG1.7.Q",
         "SIG1.8.Q",
         "SIG1.9.Q",
         "SIG1.10.Q",
         "SIG1.11.Q",
         "SIG1.12.Q",
         "SIG1.13.Q",
         "SIG1.14.Q",
         "SIG1.15.Q",
         "SIG1.16.Q",
         "SIG1.17.Q",
         "SIG1.18.Q",
         "SIG1.19.Q",
         "SIG1.20.Q",
         "SIG1.21.Q",
         "SIG1.22.Q",
         "SIG1.23.Q",
         "SIG1.24.Q",
         "SIG1.25.Q")


#generate admixture metadata
famfile=read.table("SIG1.fam",header=FALSE)
head(famfile);dim(famfile)
famfile=famfile[,1:2]
famfile$IndOrder=c(1:dim(famfile)[1])
famfile$SampleID=paste(famfile$V1,famfile$V2,sep="_")
head(famfile)
popinfo=readxl::read_xlsx("POPINFO_152023.xlsx")


popinfo2=merge(popinfo,famfile,by.x = "FID", by.y = "V2") ##FORCHANGED##

###### popinfo$present <- famfile$V2 %in% popinfo$IID
#sample_list <- subset(popinfo, present == "FALSE")
#sample_list <- sample_list[, c(1:10)]

popinfo2=popinfo2[order(popinfo2$IndOrder),]

head(popinfo2);names(popinfo2)

#creating a dataframe object for plot labels
names(popinfo2)

labs <- popinfo2[,c("POPULATION","COUNTRY_ORIGIN","LANGUAGE","CONTINENT","IndOrder")]###ABI K LIYE
labs$POPULATION=as.character(labs$POPULATION)
labs$CONTINENT=factor(labs$CONTINENT)
labs$CONTINENT=as.character(labs$CONTINENT)
labs$POPORD=paste(labs$CONTINENT,labs$SUPERPOPULATION,labs$COUNTRY_ORIGIN,labs$POPULATION,sep="_")
names(labs)

labs=data.frame(labs[,c(1,6)]) 
str(labs)

#Plotting all Ks
alist <- readQ(files=afiles[1:24])
plotQ(alist,imgoutput="sep",showindlab=F,sharedindlab=F,
      grplab=labs,grplabangle=90,linesize=0.1,pointsize=0.5,divtype=1,divsize=0.1,
      grplabjust="right",grpmean=F,showgrplab=T,
      selgrp="POPORD",divgrp="POPULATION",
      ordergrp=T,showlegend=F,showtitle=F,showsubtitle=F,
      titlelab=NA,subtitlelab=NA,
      clustercol=NA, 
      indlabsize=0,indlabheight=0,indlabspacer=-1,
      grplabsize=4.0,grplabpos=0.7,grplabheight=30,barsize=1,
      panelratio=c(1,1),width=150, 
      returnplot=T,exportplot=T,showyaxis=T,showsp=F,sortind="Cluster1",
      barbordercolour="white",barbordersize=0,outputfilename=NA,imgtype="png",exportpath=getwd())

