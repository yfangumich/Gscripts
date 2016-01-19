library(foreign)
library(plyr)
#### phenotypes ####
## Phenotypes
pheno<-read.spss("Z:/Data Analysis/Yu Fang/HRS/Phenotype/YuSrijanVars.sav",to.data.frame = TRUE)
pheno$LOCAL_ID<-as.integer(as.character(pheno$LOCAL_ID))
## Principal Components
PCs<-read.csv("Z:/Data Analysis/Yu Fang/HRS/GenoInfo/Principal_components.csv")
PCcols<-c("EV1","EV2","EV3","EV4","EV5","EV6","EV7","EV8","EV9","EV10")
## geno fam
genfam<-read.table("Z:./././././Data Analysis/Yu Fang/HRS/GeneImputed/HRS_chr22_awk.fam")
## annotation
Anno<-read.csv("Z:./././././Data Analysis/Yu Fang/HRS/GenoInfo/Sample_annotation.csv")
Annosub<-Anno[c("subjectID","local.id")]
AnnoGenfam<-merge(genfam,Annosub,by.x="V2",by.y="subjectID")
AnnoGenfam<-merge(AnnoGenfam,PCs,by.x="V2",by.y="subjectID")
## phenotype file preparation
PhenoPrep=function(targetpheno){
testvar=c(targetpheno)
pheno_use=pheno[c("hhidpn","LOCAL_ID","Age","Gender","Hispanic","Race","Education",testvar)]
# merge 
AnnoGenfamPheno<-merge(AnnoGenfam,pheno_use,by.x="local.id",by.y="LOCAL_ID")
# pick EAonly
AnnoGenfamPheno_EA<-AnnoGenfamPheno[which(AnnoGenfamPheno$Race==unique(AnnoGenfamPheno$Race)[1]),]
# pick pheno
Pheno_picked<-AnnoGenfamPheno_EA[which(!is.na(AnnoGenfamPheno_EA[,paste(targetpheno)])),]
# output
Pheno_neuroticism_out<-Pheno_picked[c("V2",targetpheno,PCcols)]
names(Pheno_neuroticism_out)[1]<-"ID"
write.table(Pheno_neuroticism_out,paste("Z:./Data Analysis/Yu Fang/HRS/Phenotype/EA_",targetpheno,".pheno",sep=""),quote=F,row.names=F,col.names = T)
}
PhenoPrep("neuroticism")

#### base ####
base_mdd=read.table("Z:././././././Data Analysis/Yu Fang/HRS/PGC/pgc.mdd.2012-04/pgc.mdd.clump.2012-04.txt",header=T)
base_mdd$seq=seq(from=1,to=nrow(base_mdd))
ped_mdd=data.frame(Chromosome=base_mdd$hg18chr,chromstart=base_mdd$bp-1,chromend=base_mdd$bp,snpid=base_mdd$snpid,seq=base_mdd$seq)
ped_mdd$Chromosome=paste("chr",ped_mdd$Chromosome,sep="")
ped_mdd$chromstart <- as.integer(ped_mdd$chromstart)
ped_mdd$chromend<- as.integer(ped_mdd$chromend)
write.table(ped_mdd,"Z:././././././Data Analysis/Yu Fang/HRS/PGC/pgc.mdd.2012-04/mdd.ped",col.names=F,row.names=F,quote=F)
write.table(ped_mdd$snpid,"Z:/Data Analysis/Yu Fang/HRS/PGC/pgc.mdd.2012-04/mdd.rs",col.names = F,row.names=F,quote = F)
# geno map
genmap<-read.table("Z:/Data Analysis/Yu Fang/HRS/GeneImputed/HRS_auto.map")
# liftover ran in linux #
lift_mdd=read.table("Z:/Data Analysis/Yu Fang/HRS/PGC/pgc.mdd.2012-04/mdd_hg19.ped")
new_mdd=merge(lift_mdd,base_mdd,by.x="V4",by.y="seq")
new_mdd$snpid=paste(as.character(new_mdd$hg18chr),":",as.character(new_mdd$V3),sep="")
new_mdd$bp=new_mdd$V3
new_mdd=new_mdd[c("snpid","hg18chr","bp","a1","a2","or","se","pval")]
new_mdd=rename(new_mdd,c("snpid"="SNP","hg18chr"="CHR","bp"="BP","a1"="A1","a2"="A2","or"="OR","se"="SE","pval"="P"))
write.table(new_mdd,"Z:/Data Analysis/Yu Fang/HRS/PGC/pgc.mdd.2012-04/base_mdd.assoc",quote=FALSE,row.names=FALSE,col.names=TRUE)
