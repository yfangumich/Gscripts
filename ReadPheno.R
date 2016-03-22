library(foreign)
library(plyr)
library("Hmisc")
#### phenotypes ####
## Phenotypes
pheno<-read.spss("Z:/Data Analysis/Yu Fang/HRS/Phenotype/YuSrijanVars.sav",to.data.frame = TRUE)
pheno$LOCAL_ID<-as.integer(as.character(pheno$LOCAL_ID))
#### mental health ####
pheno=pheno[c("hhidpn","LOCAL_ID","Age","Gender","Hispanic","Race","Education",
              "Depression","anxiety","anger","extraversion","conscientiousness","agreeableness","openness")]
#### CIDI ####
CIDI=read.spss("Z:././././././Data Analysis/Yu Fang/HRS/Phenotype/MasterCIDI_3subscales.sav",to.data.frame = TRUE)
CIDI=CIDI[c("hhidpn","CIDIDepression","dysphoriaEVER","anhedoniaEVER","DysphoriaSeverity","AnhedoniaSeverity","CIDISeverity")]
#### Height ####
Height=read.spss("Z:./././././Data Analysis/Yu Fang/HRS/Phenotype/Height.sav",to.data.frame=TRUE)
#### merge CIDI to pheno ####
pheno<-merge(pheno,CIDI,by="hhidpn",all.x=TRUE,all.y=TRUE)
#### merge Height to pheno ####
pheno<-merge(pheno,Height,by="hhidpn",all.y=TRUE)
democols<-c("Age","Gender")
## Principal Components
PCs<-read.csv("Z:/Data Analysis/Yu Fang/HRS/GenoInfo/Principal_components.csv")
PCcols<-c("EV1","EV2","EV3","EV4","EV5","EV6","EV7","EV8","EV9","EV10")
## geno fam
genfam<-read.table("Z:/Data Analysis/Yu Fang/HRS/GeneImputed/HRS_imputed.fam")
## annotation
Anno<-read.csv("Z:/Data Analysis/Yu Fang/HRS/GenoInfo/Sample_annotation.csv")
Annosub<-Anno[c("subjectID","local.id")]
AnnoGenfam<-merge(genfam,Annosub,by.x="V2",by.y="subjectID")
AnnoGenfam<-merge(AnnoGenfam,PCs,by.x="V2",by.y="subjectID")

## phenotype file preparation
PhenoPrep=function(targetpheno,filename){
testvar=c(targetpheno)
pheno_use=pheno[c("hhidpn","LOCAL_ID","Age","Gender","Hispanic","Race","Education",testvar)]
# merge 
AnnoGenfamPheno<-merge(AnnoGenfam,pheno_use,by.x="local.id",by.y="LOCAL_ID")
# pick by race
AnnoGenfamPheno_EA<-AnnoGenfamPheno[which(AnnoGenfamPheno$Race==unique(AnnoGenfamPheno$Race)[1]),]
AnnoGenfamPheno_AA<-AnnoGenfamPheno[which(AnnoGenfamPheno$Race==unique(AnnoGenfamPheno$Race)[3]),]
# pick pheno
Pheno_EA_picked<-na.omit(AnnoGenfamPheno_EA[,c("V2",PCcols,democols,targetpheno)])
Pheno_AA_picked<-na.omit(AnnoGenfamPheno_AA[,c("V2",PCcols,democols,targetpheno)])
# output
Pheno_EA_picked_out<-Pheno_EA_picked[c("V2",targetpheno)]
Pheno_AA_picked_out<-Pheno_AA_picked[c("V2",targetpheno)]
names(Pheno_EA_picked_out)[1]<-"ID"
names(Pheno_AA_picked_out)[1]<-"ID"
write.table(Pheno_EA_picked_out,paste("Z:/Data Analysis/Yu Fang/HRS/Phenotype/EA_",filename,".pheno",sep=""),quote=F,row.names=F,col.names = T)
write.table(Pheno_AA_picked_out,paste("Z:/Data Analysis/Yu Fang/HRS/Phenotype/AA_",filename,".pheno",sep=""),quote=F,row.names=F,col.names = T)
Pheno_EA_picked_out_covary<-Pheno_EA_picked[c("V2",PCcols,democols)]
Pheno_AA_picked_out_covary<-Pheno_AA_picked[c("V2",PCcols,democols)]
names(Pheno_EA_picked_out_covary)[1]="IID"
names(Pheno_AA_picked_out_covary)[1]="IID"
write.table(Pheno_EA_picked_out_covary,paste("Z:/Data Analysis/Yu Fang/HRS/Phenotype/EA_",filename,".covary",sep=""),quote=F,row.names=F,col.names = T)
write.table(Pheno_AA_picked_out_covary,paste("Z:/Data Analysis/Yu Fang/HRS/Phenotype/AA_",filename,".covary",sep=""),quote=F,row.names=F,col.names = T)
}
#PhenoPrep("neuroticism","neuroticism")
#PhenoPrep("CIDIDepression","CIDIDepression")
#PhenoPrep(c("AnhedoniaSeverity","DysphoriaSeverity"),"CIDIADSeverity")
#PhenoPrep("Height","Height")
PhenoPrep(c("Depression","anxiety","anger","extraversion","conscientiousness","agreeableness","openness"),"mental")

######## pheno correlation ########
Pheno_mental=read.table("Z:/Data Analysis/Yu Fang/HRS/Phenotype/EA_mental.pheno",header=T)
Pheno_neuroticism=read.table("Z:/Data Analysis/Yu Fang/HRS/Phenotype/EA_neuroticism.pheno",header=T)
Pheno_all=merge(Pheno_mental,Pheno_neuroticism,by="ID")
Pheno_cor=Pheno_all[c(-1)]
a=rcorr(as.matrix(Pheno_cor))
write.csv(a$r,"work/Fitbit/FBscripts/tmp.csv")
df=merge(prof,Pheno_all,by.x="IID",by.y="ID")

covariates <- "EV1,EV2,EV3,EV4,EV5,EV6,EV7,EV8,EV9,EV10,Age,Gender"
covariates <- strsplit(covariates, split=",")[[1]]
outcomes=c("Depression","anxiety","anger","extraversion","conscientiousness","agreeableness","openness")
pheno.results<-data.frame(name=character(length(outcomes)),r2ratio=numeric(length(outcomes)),
                          r2score=numeric(length(outcomes)),r2pheno=numeric(length(outcomes)),
                          pscore=numeric(length(outcomes)),ppheno=numeric(length(outcomes)),
                          estimatescore=numeric(length(outcomes)),estimatepheno=numeric(length(outcomes)))
pheno.results$name<-as.character(pheno.results$name)
for (i in 1:length(outcomes)){
  outcome<-outcomes[i]
  pheno.results$name[i]<-outcome
  frm<-as.formula(paste(outcome,"~ ."))
  model.full  <- glm(formula=frm, family="gaussian", data= df[,c(outcome,"SCORE","neuroticism", covariates)])
  model.logit <- glm(formula=frm, family="gaussian", data = df[,c(outcome, "SCORE", covariates)])
  model.null  <-  glm(formula=frm, family="gaussian", data = df[,c(outcome, covariates)])
  r2full  <- 1-model.full$deviance/model.full$null.deviance
  r2logit <- 1-model.logit$deviance/model.logit$null.deviance
  r2null  <- 1-model.null$deviance/model.null$null.deviance
  pheno.results$r2pheno[i]<-r2full - r2logit
  pheno.results$r2score[i]<-r2logit - r2null
  pheno.results$r2ratio[i]<-pheno.results$r2score[i] / (r2full - r2null)  
  pheno.results$pscore[i]<-summary(model.full)$coefficients[2,4]
  pheno.results$estimatescore[i]<-summary(model.full)$coefficients[2,1]
  pheno.results$ppheno[i]<-summary(model.full)$coefficients[3,4]
  pheno.results$estimatepheno[i]<-summary(model.full)$coefficients[3,1]
}
attach(pheno.results)
write.table(pheno.results[order(r2ratio),],"Z:/Data Analysis/Yu Fang/HRS/PRSoutput/tmp.csv",col.names=T,row.names=F,quote=F,sep=",")
score.results<-read.csv("Z:/Data Analysis/Yu Fang/HRS/PRSoutput/2016-03-10_neuroticism_mental/PRSice_TOP_SCORES_ACROSS_PHENOTYPES.txt",header=T)
pheno.results$ratio<-score.results$top.thresh.r2/pheno.results$r2
m1<-glm(formula=frm, family="gaussian", data = df[,c(outcome,"SCORE","neuroticism", covariates)])
summary(m1)

######## base ########
# geno map
genmap<-read.table("Z:/Data Analysis/Yu Fang/HRS/GeneImputed/HRS_auto.map")
# read base - mdd
base_ds1<-read.table("Z:/Data Analysis/Yu Fang/HRS/Public/pgc.mdd.2012-04/pgc.mdd.clump.2012-04.txt",header=T)
base_ds1_full<-read.table("Z:/Data Analysis/Yu Fang/HRS/Public/pgc.mdd.2012-04/pgc.mdd.full.2012-04.txt",header=T)
base_col<-c("hg18chr","bp","snpid")
base_ds2<-read.table("Z:/Data Analysis/Yu Fang/HRS/Public/pgc.cross.mdd/pgc.cross.MDD9.2013-05.txt",header=T)
base_ds3<-read.table("Z:./././././Data Analysis/Yu Fang/HRS/Public/HanWomenMDD/HanWomenMDDSep2015.txt",header=T)
# read base - height
LangoAllen2010<-read.table("Z:/Data Analysis/Yu Fang/HRS/Public/GIANT/GIANT_HEIGHT_LangoAllen2010_publicrelease_HapMapCeuFreq.txt",header=T)
Wood2014<-read.table("Z:/Data Analysis/Yu Fang/HRS/Public/GIANT/GIANT_HEIGHT_Wood_et_al_2014_publicrelease_HapMapCeuFreq.txt",header=T)
Yang2012<-read.table("Z:/Data Analysis/Yu Fang/HRS/Public/GIANT/GIANT_Yang2012Nature_publicrelease_HapMapCeuFreq_Height.txt",header=T)
Randall2013MEN<-read.table("Z:/Data Analysis/Yu Fang/HRS/Public/GIANT/GIANT_Randall2013PlosGenet_stage1_publicrelease_HapMapCeuFreq_HEIGHT_MEN_N.txt",header=T)
Randall2013WOMEN<-read.table("Z:/Data Analysis/Yu Fang/HRS/Public/GIANT/GIANT_Randall2013PlosGenet_stage1_publicrelease_HapMapCeuFreq_HEIGHT_WOMEN_N.txt",header=T)
Berndt2013EXTREME<-read.table("Z:/Data Analysis/Yu Fang/HRS/Public/GIANT/GIANT_EXTREME_HEIGHT_Stage1_Berndt2013_publicrelease_HapMapCeuFreq.txt",header=T)
nrow(LangoAllen2010)
nrow(Wood2014)
nrow(Yang2012)
nrow(Randall2013MEN)
nrow(Randall2013WOMEN)
nrow(Berndt2013EXTREME)

#base_ds<-Randall2013MEN
#base_ds<-Randall2013WOMEN
#base_ds<-Wood2014
base_ds<-base_ds2
#base_col<-c("MarkerName")
base_col<-c("hg18chr","bp","snpid")
#base_rspath<-"Z:/Data Analysis/Yu Fang/HRS/Public/GIANT/Randall2013MEN.rs"
#base_rspath<-"Z:/Data Analysis/Yu Fang/HRS/Public/GIANT/Randall2013WOMEN.rs"
#base_rspath<-"Z:/Data Analysis/Yu Fang/HRS/Public/GIANT/Wood2014.rs"
base_rspath<-"Z:/Data Analysis/Yu Fang/HRS/Public/pgc.cross.mdd/mdd.ped"

## Prepare for lift, if necessary ##
BasePrep=function(base_ds,base_col){
  base_ds$seq=seq(from=1,to=nrow(base_ds))
  #ped_ds=data.frame(snpid=base_ds[,paste(base_col[1])],seq=base_ds$seq)
  ped_ds=data.frame(Chromosome=base_ds[,paste(base_col[1])],chromstart=base_ds[,paste(base_col[2])]-1,chromend=base_ds[,paste(base_col[2])],snpid=base_ds[,paste(base_col[3])],seq=base_ds$seq)
  ped_ds$Chromosome=paste("chr",ped_ds$Chromosome,sep="")
  ped_ds$chromstart <- as.integer(ped_ds$chromstart)
  ped_ds$chromend<- as.integer(ped_ds$chromend)
  write.table(ped_ds,base_rspath,col.names = F,row.names=F,quote = F)
  #write.table(ped_ds$snpid,base_rspath,col.names = F,row.names=F,quote = F)
  return(ped_ds)
}
ped_ds=BasePrep(base_ds,base_col)
## pause: liftover ran in linux ##
## after lift or lift is not needed ##
 #lift_ds<-read.table("Z:/Data Analysis/Yu Fang/HRS/Public/pgc.mdd.2012-04/mdd_hg19.ped")
 lift_ds<-read.table("Z:/Data Analysis/Yu Fang/HRS/Public/pgc.cross.mdd/mdd_hg19.ped")
 lift_cols<-c("snpid","hg18chr","bp","a1","a2","or","se","pval")
 new_cols<-c("SNP","CHR","BP","A1","A2","oR","SE","P")

## save base file ##
#lift_ds<-read.table("Z:/Data Analysis/Yu Fang/HRS/Public/GPC-2.NEUROTICISM/GPC-2.NEUROTICISM.full.txt")
#lift_ds<-Wood2014
#lift_ds<-read.table("Z:/Data Analysis/Yu Fang/HRS/Public/GIANT/Randall2013MEN_lifted.rs")
#lift_ds<-Randall2013MEN
#lift_ds<-Randall2013WOMEN
lift_ds<-base_ds3
#lift_ds<-lift_ds[which(lift_ds$V10>0.02 & lift_ds$V10<0.98),]
#lift_cols<-c("V1","V2","V3","V4","V5","V6","V7","V8")
#lift_cols<-c("MarkerName","Allele1","Allele2","b","SE","p")
#lift_cols<-c("MarkerName","A1","A2","BETA","SE.2gc","P.2gc")
#lift_cols<-c("RSID","ALT_ALLELE","REF_ALLELE","OR.logistic","SE","P.lmm")
lift_cols<-c("V4","hg18chr","bp","a1","a2","or","se","pval")
#new_cols<-c("SNP","A1","A2","BETA","SE","P")
new_cols<-c("SNP","CHR","BP","A1","A2","OR","SE","P")  # PRSice treats A1 as effect(alternative) allele
#outpath_ds<-"Z:/Data Analysis/Yu Fang/HRS/Public/GIANT/Height_Wood2014.assoc"
#outpath_ds<-"Z:/Data Analysis/Yu Fang/HRS/Public/GIANT/Height_Randall2013MEN.assoc"
#outpath_ds<-"Z:/Data Analysis/Yu Fang/HRS/Public/GIANT/Height_Randall2013WOMEN.assoc"
#outpath_ds<-"Z:/Data Analysis/Yu Fang/HRS/Public/HanWomenMDD/HanWomenMDDSep2015.assoc"
#outpath_ds<-"Z:/Data Analysis/Yu Fang/HRS/Public/pgc.mdd.2012-04/base_mdd.assoc"
outpath_ds<-"Z:/Data Analysis/Yu Fang/HRS/Public/pgc.cross.mdd/base_mdd.assoc"
lifted<-T
needtoupper<-F

LiftedPrep=function(lift_ds){
 if (lifted){
  new_ds=merge(lift_ds,base_ds,by.x=c("V4","V5"),by.y=c("snpid","seq"))
  new_ds$bp=new_ds$V3 
 } else{  new_ds=lift_ds }
 # commonsnp=intersect(new_ds$lift_cols[1],genmap$V2)
 #  print(nrow(new_ds))
 #  print(length(commonsnp))
  new_ds=new_ds[lift_cols]
  names(new_ds)=new_cols
 if(needtoupper){
  new_ds$A1<-toupper(lift_ds$A1)
  new_ds$A2<-toupper(lift_ds$A2)
 }
  write.table(new_ds,outpath_ds,quote=FALSE,row.names=FALSE,col.names=TRUE)
 return(new_ds)
}
new_ds=LiftedPrep(lift_ds)






