rm(list=ls())
install.packages("PopGenome",repos = "http://cran.us.r-project.org")
library("PopGenome")
gene_list<-read.table("gene_list.txt",header=TRUE)
genes1<-gene_list$gene
genes<-genes1[c(-68,-90)] #drop AR and CHRDL1 from chrX
CHB<-readData("AADAC/CHB",format ="VCF" ) 
CEU<-readData("AADAC/CEU",format="VCF")  
YRI<-readData("AADAC/YRI",format="VCF")  
CHB<-c(get.individuals(CHB)[[1]]) #pop1
CEU<-c(get.individuals(CEU)[[1]]) #pop2
YRI<-c(get.individuals(YRI)[[1]]) #pop3

results<-data.frame(gene=character(),Fst=numeric(),p_Fst=numeric(),Fst1=numeric(),p_Fst1=numeric(),Fst2=numeric(),p_Fst2=numeric(),Fst3=numeric(),p_Fst3=numeric(),PBS1=numeric(),p_PBS1=numeric(),PBS2=numeric(),p_PBS2=numeric(),PBS3=numeric(),p_PBS3=numeric())
for (gene in genes){
  setwd(gene)
  ALL<-readData("ALL",format="VCF")  
  ALL<-set.populations(ALL,list(CHB,CEU,YRI))
  ALL<-F_ST.stats(ALL, mode="nucleotide")
  Fst<-as.numeric(ALL@nucleotide.F_ST)
  stat1<-as.data.frame(ALL@nuc.F_ST.pairwise)
  Fst1<-stat1[1,1] #pop1/pop2 (CHB/CEU)
  Fst2<-stat1[2,1] #pop1/pop3 (CHB/YRI)
  Fst3<-stat1[3,1] #pop2/pop3 (CEU/YRI)
  T_1A2<- -log(1-Fst1) 
  T_1A3<- -log(1-Fst2)
  T_2A3<- -log(1-Fst3)
  PBS1<- (T_1A2+T_1A3-T_2A3)/2 #CHB and last common ancestor
  PBS2<- (T_1A2+T_2A3-T_1A3)/2 #CEU and last common ancestor
  PBS3<- (T_1A3+T_2A3-T_1A2)/2 #YRI and last common ancestor
  
  ms<-readMS(paste(gene,".ms.gz",sep=""))
  ms<-set.populations(ms,list(c(get.individuals(ms)[[1]][415:620]),c(get.individuals(ms)[[1]][217:414]),c(get.individuals(ms)[[1]][1:216])))
  ms<-F_ST.stats(ms,mode="nucleotide")
  Fst_1<-as.data.frame(ms@nucleotide.F_ST)
  p_Fst<-sum(Fst_1 >= Fst,na.rm = TRUE)/sum(!is.na(Fst_1))
  stat1_1<-as.data.frame(ms@nuc.F_ST.pairwise)
  Fst1_1<-stat1_1[1,] #pop1/pop2 (CHB/CEU)
  Fst2_1<-stat1_1[2,] #pop1/pop3 (CHB/YRI)
  Fst3_1<-stat1_1[3,] #pop2/pop3 (CEU/YRI)
  T_1A2_1<- -log(1-Fst1_1) 
  T_1A3_1<- -log(1-Fst2_1)
  T_2A3_1<- -log(1-Fst3_1)
  PBS1_1<- (T_1A2_1+T_1A3_1-T_2A3_1)/2 #CHB and last common ancestor
  PBS2_1<- (T_1A2_1+T_2A3_1-T_1A3_1)/2 #CEU and last common ancestor
  PBS3_1<- (T_1A3_1+T_2A3_1-T_1A2_1)/2 #YRI and last common ancestor
  p_Fst1<-sum(Fst1_1 >= Fst1,na.rm = TRUE)/sum(!is.na(Fst1_1))
  p_Fst2<-sum(Fst2_1 >= Fst2,na.rm = TRUE)/sum(!is.na(Fst2_1))
  p_Fst3<-sum(Fst3_1 >= Fst3,na.rm = TRUE)/sum(!is.na(Fst3_1))
  p_PBS1<-sum(PBS1_1 >= PBS1,na.rm = TRUE)/sum(!is.na(PBS1_1))
  p_PBS2<-sum(PBS2_1 >= PBS2,na.rm = TRUE)/sum(!is.na(PBS2_1))
  p_PBS3<-sum(PBS3_1 >= PBS3,na.rm = TRUE)/sum(!is.na(PBS3_1))
  
  temp<-data.frame(gene=gene,Fst=Fst,p_Fst=p_Fst,Fst1=Fst1,p_Fst1=p_Fst1,Fst2=Fst2,p_Fst2=p_Fst2,Fst3=Fst3,p_Fst3=p_Fst3,PBS1=PBS1,p_PBS1=p_PBS1,PBS2=PBS2,p_PBS2=p_PBS2,PBS3=PBS3,p_PBS3=p_PBS3)
  results<-rbind(results,temp)
  
  setwd("..")
  }

write.table(results,"Fst_results.txt",sep = "\t",row.names=FALSE)
 
