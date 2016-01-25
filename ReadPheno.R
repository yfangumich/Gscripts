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
genfam<-read.table("Z:/Data Analysis/Yu Fang/HRS/GeneImputed/HRS_chr22_awk.fam")
## annotation
Anno<-read.csv("Z:/Data Analysis/Yu Fang/HRS/GenoInfo/Sample_annotation.csv")
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
Pheno_picked_out<-Pheno_picked[c("V2",targetpheno,PCcols)]
names(Pheno_picked_out)[1]<-"ID"
write.table(Pheno_neuroticism_out,paste("Z:/Data Analysis/Yu Fang/HRS/Phenotype/EA_",targetpheno,".pheno",sep=""),quote=F,row.names=F,col.names = T)
Pheno_picked_out_covary<-Pheno_picked[c("V2",PCcols)]
names(Pheno_picked_out_covary)[1]="IID"
write.table(Pheno_picked_out_covary,"Z:/Data Analysis/Yu Fang/HRS/Phenotype/Top10PC.covary",quote=F,row.names=F,col.names = T)
}
PhenoPrep("neuroticism")

#### base ####
# geno map
genmap<-read.table("Z:/Data Analysis/Yu Fang/HRS/GeneImputed/HRS_auto.map")
# read base
base_ds<-read.table("Z:/Data Analysis/Yu Fang/HRS/Public/pgc.mdd.2012-04/pgc.mdd.clump.2012-04.txt",header=T)
base_col<-c("hg18chr","bp","snpid")
BasePrep=function(base_ds,base_col){
  base_ds$seq=seq(from=1,to=nrow(base_ds))
  ped_ds=data.frame(Chromosome=base_ds[,paste(base_col[1])],chromstart=base_ds[,paste(base_col[2])]-1,chromend=base_ds[,paste(base_col[2])],snpid=base_ds[,paste(base_col[3])],seq=base_ds$seq)
  ped_ds$Chromosome=paste("chr",ped_ds$Chromosome,sep="")
  ped_ds$chromstart <- as.integer(ped_ds$chromstart)
  ped_ds$chromend<- as.integer(ped_ds$chromend)
  write.table(ped_ds,"Z:././././././Data Analysis/Yu Fang/HRS/Public/pgc.mdd.2012-04/mdd.ped",col.names=F,row.names=F,quote=F)
  write.table(ped_ds$snpid,"Z:/Data Analysis/Yu Fang/HRS/Public/pgc.mdd.2012-04/mdd.rs",col.names = F,row.names=F,quote = F)
}
BasePrep(base_ds,base_col)

# pause: liftover ran in linux #

# after lift or lift is not needed

# lift_ds<-read.table("Z:/Data Analysis/Yu Fang/HRS/Public/pgc.mdd.2012-04/mdd_hg19.ped")
# lift_cols<-c("snpid","hg18chr","bp","a1","a2","or","se","pval")
# lifted<-T
# new_cols<-c("SNP","CHR","BP","A1","A2","oR","SE","P")
# outpath_ds<-"Z:/Data Analysis/Yu Fang/HRS/Public/pgc.mdd.2012-04/base_mdd.assoc"
lift_ds<-read.table("Z:/Data Analysis/Yu Fang/HRS/Public/GPC-2.NEUROTICISM/GPC-2.NEUROTICISM.full.txt")
lift_ds<-lift_ds[which(lift_ds$V10>0.02 & lift_ds$V10<0.98),]
lift_cols<-c("V1","V2","V3","V4","V5","V6","V7","V8")
lifted<-F
new_cols<-c("SNP","CHR","BP","A1","A2","BETA","SE","P")
outpath_ds<-"Z:/Data Analysis/Yu Fang/HRS/Public/GPC-2.NEUROTICISM/GPC_neuroticism.assoc"


LiftedPrep=function(lift_ds){
 if (lifted){
  new_ds=merge(lift_ds,base_ds,by.x="V4",by.y="seq")
  new_ds$bp=new_ds$V3 
 } else{  new_ds=lift_ds }
 # commonsnp=intersect(new_ds$lift_cols[1],genmap$V2)
 #  print(nrow(new_ds))
 #  print(length(commonsnp))
  new_ds=new_ds[lift_cols]
  names(new_ds)=new_cols
  write.table(new_ds,outpath_ds,quote=FALSE,row.names=FALSE,col.names=TRUE)
}
LiftedPrep(lift_ds)






