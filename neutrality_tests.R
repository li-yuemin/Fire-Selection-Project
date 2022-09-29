rm(list=ls())
install.packages("PopGenome",repos = "http://cran.us.r-project.org")
library("PopGenome")
gene_list<-read.table("gene_list.txt",header=TRUE)
genes1<-gene_list$gene
genes<-genes1[c(-68,-90)] #drop AR and CHRDL1 from chrX
pops<-c("CHB","CEU","YRI")

results<-data.frame(gene=character(),pop=character(),Tajima.D=numeric(),p_D=numeric(),Fu.Li.F=numeric(),p_F=numeric(),Fu.Li.D=numeric(),p_D1=numeric())
for (gene in genes){
  setwd(gene)
  ms<-readMS(paste(gene,".ms.gz",sep=""))
  ms<-set.populations(ms,list(c(get.individuals(ms)[[1]][415:620]),c(get.individuals(ms)[[1]][217:414]),c(get.individuals(ms)[[1]][1:216])))
  ms<-neutrality.stats(ms)
  for(pop in pops){
    VCF<-readData(pop,format ="VCF" )
    VCF<-neutrality.stats(VCF)
    stat<-as.data.frame(get.neutrality(VCF)[[1]])
    if(pop=="CHB"){
      stat1<-as.data.frame(get.neutrality(ms)[[1]])
    }
    if(pop=="CEU"){
      stat1<-as.data.frame(get.neutrality(ms)[[2]])
    }
    if(pop=="YRI"){
      stat1<-as.data.frame(get.neutrality(ms)[[3]])
    }
    p_D<-sum(stat1$Tajima.D >= stat$Tajima.D,na.rm = TRUE)/sum(!is.na(stat1$Tajima.D))
    p_F<-sum(stat1$Fu.Li.F >= stat$Fu.Li.F,na.rm = TRUE)/sum(!is.na(stat1$Fu.Li.F))
    p_D1<-sum(stat1$Fu.Li.D >= stat$Fu.Li.D,na.rm = TRUE)/sum(!is.na(stat1$Fu.Li.D))
    temp<-data.frame(gene=gene,pop=pop,Tajima.D=stat$Tajima.D,p_D=p_D,Fu.Li.F=stat$Fu.Li.F,p_F=p_F,Fu.Li.D=stat$Fu.Li.D,p_D1=p_D1)
    results<-rbind(results,temp)
  }
  setwd("..")
}

write.table(results,"neutrality_results.txt",sep = "\t",row.names=FALSE)




