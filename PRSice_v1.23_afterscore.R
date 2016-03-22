library(batch)
start.time <- proc.time()[3]
options(echo = FALSE)
options(warn=-1)
cat(" ################################# \n # \n # \n # \n # \n # PRSice: Polygenic Risk Score software \n # \n # Jack Euesden, Cathryn M. Lewis, Paul F. O'Reilly 2014 \n # \n # \n # If you use PRSice in published work, please cite: \n # \n # \"PRSice: Polygenic Risk Score software\" \n # Euesden, Lewis, O'Reilly, Bioinformatics (2015) 31 (9):1466-1468 \n # \n # \n # \n #  \n ################################# \n")
#################################
#
#  Default options
#
#################################

## essential
target <-  NA
base <-   NA

## preferable
plink <-  NA
order.cols <- "SNP,CHR,BP,A1,A2,OR,SE,P"
supplied.order <- F

# phenotype options
pheno.file <-   NA
binary.target <-  T

# covariate options
covary <- T
covariates <- "C1,C2"
user.covariate.file <-  NA
ancestry.dim <- "MDS"

# graphical parameters
ggfig <-  T
barchart.levels <- "0.001, 0.05, 0.1, 0.2, 0.3, 0.4, 0.5"
barpalatte <- "YlOrRd"
best.thresh.on.bar <- F
scatter.R2 <- F						
figname <- "PRSice"
bar.col.is.pval <- T
bar.col.is.pval.lowcol <- "dodgerblue"
bar.col.is.pval.highcol <- "firebrick"

# clumping							
clump.snps <- T
clump.p1 <- 1
clump.p2 <- 1
clump.r2 <- 0.1
clump.kb <- 250

# pruning
prune.snps <- F
prune.kb.wind <- 50
prune.kb.step <- 2
prune.kb.r2 <- 0.8

							
# high density scoring							
slower <- 0.0001
supper <- 0.5
sinc <-  0.00005
fastscore <- F
							
# dosage
dosage <- F
dosage.format <- "gen"
dos.skip0 <- 1
dos.skip1 <- 1
dos.coding <- 1
dos.format <- 3
dos.sep.fam <- NA
dos.fam.is.samp <- F
dos.impute2 <- F
dos.path.to.lists <- NA
dos.list.file <- NA

# alternate genotype formats
geno.is.ped <- F
geno.as.list <- F

# miscellaneous options
wd <-  "./"							
print.time <-  T
cleanup <- T
plinkpath <-  "./"
remove.mhc <- F
for.meta <- T    # FY: want to get the coefficients
plink.silent <- T
allow.no.sex <- F

# A Better Coefficient of Determination
calculate.abc.r2 <- F
pi0 <- 0
n.ca.base <- NA
n.co.base <- NA
n.ca.targ <- NA
n.co.targ <- NA
prev.base <- NA
prev.targ <- NA

# quantiles
quantiles <- F
num.quantiles <- 5
quant.ref <- NA

mend.score <- F
mend.score.len <- 100
score.at.1 <- F
report.individual.scores <- T
report.best.score.only <- T
no.regression <- F
debug.mode <- F

# sumsum
sumsum <- F
clump.ref <- NA
size.targ <- NA

# multiple phenotype options
multiple.target.phenotypes <- F
target.phenotypes <-  NA #  "V2,V3"
target.phenotypes.binary <- NA # "T,T"  ## NB:: 'QT' or 'BIN' 
multiple.base.phenotypes <- F
base.phenotypes.names <- NA  ## sub in for PHEN.NAME
heat.r2 <- F

if(!fastscore){
  best.thresh.on.bar <- T
}

if(dosage){
  quantiles <- F
}



#cat(" ################################# \n # \n #  Begin Script \n # \n ################################# \n")


if(Sys.info()[1] == "Darwin"){
	os <- "mac"
}
if(Sys.info()[1] == "Linux"){
	os <- "linux"
}
if(Sys.info()[1] == "Windows"){
	os <- "windows"
}

#if(os == "windows"){
#	print("ERROR: Windows not supported")
#	quit()
#}
													

if(dosage){
	sinc <- 0.001
	slower <- 0.001
	supper <- 0.5
}


cat(" ################################# \n # \n #  Read in Command Line Arguments & interpret \n # \n ################################# \n")

##### manually input the customized options ####

# frequently changed
#base <- "Z:/Data Analysis/Yu Fang/HRS/Public/pgc.mdd.2012-04/base_mdd.assoc"
#base <- "Z:/Data Analysis/Yu Fang/HRS/Public/GPC-2.NEUROTICISM/GPC_neuroticism.assoc"
#base <- "Z:/Data Analysis/Yu Fang/HRS/Publig/GIANT/Height_Wood2014.assoc"
base <- "Z:/Data Analysis/Yu Fang/HRS/Public/HanWomenMDD/HanWomenMDDSep2015_use.assoc"

#pheno.file <- "Z:/Data Analysis/Yu Fang/HRS/Phenotype/EA_neuroticism.pheno"
#pheno.file <- "Z:/Data Analysis/Yu Fang/HRS/Phenotype/EA_CIDISeverity.pheno"
#pheno.file <- "Z:/Data Analysis/Yu Fang/HRS/Phenotype/EA_CIDIDepression.pheno"
#pheno.file <- "Z:/Data Analysis/Yu Fang/HRS/Phenotype/EA_Height.pheno"
#pheno.file <- "Z:/Data Analysis/Yu Fang/HRS/Phenotype/EA_mental.pheno"

#target.phenotypes <- "neuroticism"
#target.phenotypes <- "CIDISeverity"
#target.phenotypes <- "CIDIDepression"
#target.phenotypes <- "Height"
#target.phenotypes <- c("Depression","anxiety","anger","extraversion","conscientiousness","agreeableness","openness")
#target.phenotypes <- "Depression"

#target.phenotypes.binary <- c(F,F,F,F,F,F,F)
target.phenotypes.binary <- T

#user.covariate.file <- "Z:/Data Analysis/Yu Fang/HRS/Phenotype/EA_Height.covary"
#user.covariate.file <- "Z:/Data Analysis/Yu Fang/HRS/Phenotype/EA_neuroticism.covary"
user.covariate.file <- "Z:/Data Analysis/Yu Fang/HRS/Phenotype/EA_CIDIDepression.covary"
#user.covariate.file <- "Z:/Data Analysis/Yu Fang/HRS/Phenotype/EA_mental.covary"

#wd <- "Z:/Data Analysis/Yu Fang/HRS/PRSoutput/2016-02-29_height_height_EA_r2"
#wd <- "Z:/Data Analysis/Yu Fang/HRS/PRSoutput/2016-03-07_neuroticism_neuroticism"
#wd <- "Z:/Data Analysis/Yu Fang/HRS/PRSoutput/2016-03-08_neuroticism_CIDIDepression"
wd <- "Z:/Data Analysis/Yu Fang/HRS/PRSoutput/2016-03-22_mdd_mdd"

# non-frequently changed but need to pay attention
covariates <- "EV1,EV2,EV3,EV4,EV5,EV6,EV7,EV8,EV9,EV10,Age,Gender"
barchart.levels <- "0.00001, 0.0001, 0.001, 0.01, 0.1, 0.2, 0.3, 0.4, 0.5"

# non-frequently changed
plink <- "C:/Users/yfang/Documents/work/Gene/Software/plink-1.07-x86_64/plink"
dosage <- T
target <- "/home/srijan_lab/yfang/HRS/GeneImputed/phg000264.v1.CIDR_HRS.genotype-imputed-data.c1/HRS_auto.gprobs"
dos.sep.fam <- "Z:/Data Analysis/Yu Fang/HRS/GeneImputed/HRS_imputed.fam"
dos.sep.map <- "Z:/Data Analysis/Yu Fang/HRS/GeneImputed/HRS_auto.map"
dos.impute2 <- T
slower <- 0.00001
supper <- 0.1
sinc <- 10
debug.mode <- T
clump.snps <- F
prune.snps <- F
multiple.target.phenotypes <- T
remove.mhc <- T
score.at.1 <- F
quantiles <- T
print.time <- T
covary <- T
cleanup <- F
allow.no.sex <- T


# manually input end

#####
#supper <- supper+sinc
#slower <- slower - sinc
covariates <- strsplit(covariates, split=",")[[1]]
order.cols <- strsplit(order.cols, split=",")[[1]]
barchart.levels <- as.numeric(strsplit(barchart.levels, split=",")[[1]])
ggfig <- as.logical(ggfig)
covary <- as.logical(covary)
fastscore <- as.logical(fastscore)
print.time <- as.logical(print.time)
supplied.order <- as.logical(supplied.order)
binary.target <- as.logical(binary.target)
best.thresh.on.bar <- as.logical(best.thresh.on.bar)
cleanup <- as.logical(cleanup)
clump.snps <- as.logical(clump.snps)
prune.snps <- as.logical(prune.snps)
scatter.R2 <- as.logical(scatter.R2)
dosage <- as.logical(dosage)
dos.fam.is.samp <- as.logical(dos.fam.is.samp)
dos.impute2 <- as.logical(dos.impute2)
remove.mhc <- as.logical(remove.mhc)
bar.col.is.pval <- as.logical(bar.col.is.pval)
for.meta <- as.logical(for.meta)
geno.is.ped <- as.logical(geno.is.ped)
geno.as.list <- as.logical(geno.as.list)
mend.score <- as.logical(mend.score)
score.at.1 <- as.logical(score.at.1)
report.individual.scores <- as.logical(report.individual.scores)
report.best.score.only <- as.logical(report.best.score.only)
plink.silent <- as.logical(plink.silent)
no.regression <- as.logical(no.regression)
multiple.target.phenotypes <- as.logical(multiple.target.phenotypes)
multiple.base.phenotypes <- as.logical(multiple.base.phenotypes)
debug.mode <- as.logical(debug.mode)
quantiles <- as.logical(quantiles)
sumsum <- as.logical(sumsum)
calculate.abc.r2 <- as.logical(calculate.abc.r2)
heat.r2 <- as.logical(heat.r2)

barchart.levels.old <- barchart.levels

if(multiple.base.phenotypes){
  base.phenotypes.names <- strsplit(gsub(" ", "", base.phenotypes.names), split=",")[[1]]
  if(is.na(base.phenotypes.names)){
    cat("ERROR: Please select base phenotypes to use, using base.phenotypes.names \n or set multiple.base.phenotypes F. \n Quitting")
    quit()
  }
}

if(multiple.target.phenotypes){
  #target.phenotypes <- strsplit(gsub(" ", "", target.phenotypes), split=",")[[1]]   # FY: this doesn't seem right...
  #target.phenotypes <- strsplit(gsub(" ", "", target.phenotypes), split=",")
  #target.phenotypes.binary  <- as.logical(strsplit(gsub(" ", "", target.phenotypes.binary), split=",")[[1]]) # FY: this doesn't seem right...
  #target.phenotypes.binary  <- as.logical(strsplit(gsub(" ", "", target.phenotypes.binary), split=",")) 
  if(is.na(target.phenotypes)){
    cat("ERROR: Please select target phenotypes to use, using target.phenotypes \n or set multiple.target.phenotypes F. \n Quitting")
    quit()
  }
  if(is.na(target.phenotypes.binary)){
    cat("ERROR: Please specify whether target phenotypes are binary or QT \n or set multiple.target.phenotypes F. \n Quitting")
    quit()
  }
  if(length(target.phenotypes) != length(target.phenotypes.binary)){
  	cat("ERROR: Different number of target phenotype names and target phenotype types specified \n Check these lists are the same length \n Quitting")
  	quit()
  } 
}

if(!sumsum & quantiles){
 if(is.na(quant.ref)){
 	quant.ref <- ceiling(num.quantiles/2)
 }
 if(!is.na(quant.ref)){
 	if(quant.ref > num.quantiles){
   	  quant.ref <- ceiling(num.quantiles/2)
 	  cat(paste("WARNING: reference quantile", quant.ref, "is greater than number of quantiles", num.quantiles, "\n Using middle quantile by default"))
    }
  }
}

if(calculate.abc.r2 & binary.target){
  n1 <- n.ca.base + n.co.base
  n2 <- n.ca.targ + n.co.targ 
  sampling1 <- n.ca.base / n1
  sampling2 <- n.ca.targ / n2
  if(is.na(n.ca.base) | is.na(n.co.base) | is.na(n.ca.targ) | is.na(n.co.targ) ){
  	cat("ERROR: Need to supply information on \n Base and Target sizes \n to use R2 on Liability Scale \n Defaulting to Nagelkerke's Pseudo R2 \n")
    calculate.abc.r2 <- F
  }
  prevalence1 <- prev.base
  prevalence2 <- prev.targ
  if(is.na(prev.targ) | is.na(prev.base) ){
  	cat("ERROR: Need to supply information on \n Base and Target phenotype prevalences \n to use R2 on Liability Scale \n Defaulting to Nagelkerke's Pseudo R2 \n")
    calculate.abc.r2 <- F
  }
}

if(debug.mode){
  options(warn=0)
  plink.silent <- F
}

if(!is.na(pheno.file)){
  ext.phen <- T
}
if(is.na(pheno.file)){
  ext.phen <- F
}
if(!is.na(user.covariate.file)){
  covary <- T
}

if(!sumsum){
if(prune.snps & clump.snps){
  cat("ERROR:: Please select either clumping or pruning (or neither) \n Defaulting to Clumping \n")
  prune.snps <- F
}
}

if(multiple.target.phenotypes){
  ext.phen <- T
  cat("Using data on multiple phenotypes for target data \n")
}



######################################
#
#
#  Compile Functions from polygenescore.R (Dudbridge 2013): End
#
#
#######################################
if(!sumsum){
for(basePhen in 1:length(base.phenotypes.names)){
	output <- as.vector(1)
	cat(" ################################# \n # \n #  Check options match \n # \n ################################# \n")
	use.beta <- F
	if(is.na(plink)){
	  cat("ERROR: Please supply a path to PLINK. \n For Download Links, see www.PRSice.info \n NB: PLINK-1.07 is required for dosage data, PLINK2 is required for genotype data \n Quitting \n"); quit()
	}
	if(wd == "./"){
	  print("Using current directory as working directory")
	}
	setwd(wd)
	if(ggfig){
	  library(ggplot2)
	  library(plyr)
	}
	if(binary.target){
	  library(fmsb)
	}
	options("scipen"=100,"digits"=4)
  
	# FY: add list file here
	lists <- read.table("rangelist.txt", head = F)
	
	if(covary & is.na(user.covariate.file)){
	  cat(" ################################# \n # \n #   Covary by generated dimensions \n # \n ################################# \n")
	}
	
	if(covary & is.na(user.covariate.file)){
	  if(!dosage){
	  	if(geno.as.list){
	  		cat("ERROR: Please merge genome to calculate ancenstry-informative dimensions \n");quit()
	  	}
	  	if(!(geno.as.list)){
	      system("echo 6 26000000 33000000 mhc > mhc.txt")
	      system(paste(plink, " --noweb --bfile ",target, "         --exclude range  mhc.txt   --make-bed  --out target_no_mhc "))
	      system(paste(plink, " --noweb --bfile   target_no_mhc     --indep-pairwise 100 25 0.2     --out prune_target"))
	      system(paste(plink, " --noweb --bfile ",target, "         --extract prune_target.prune.in    --make-bed   --out prune_target_out"))
	      system(paste(plink, " --noweb --bfile   prune_target_out  --genome     --out prune_target_out"))
	      if(ancestry.dim == "MDS"){
	        system(paste(plink, " --noweb --bfile   prune_target_out  --read-genome prune_target_out.genome   --cluster --mds-plot", length(covariates), " --out 	  ANCESTRY_INFORMATIVE_DIMENSIONS"))
	        pc.file <- "ANCESTRY_INFORMATIVE_DIMENSIONS.mds"
	      }
	      if(ancestry.dim == "PCA"){
	        system(paste(plink, " --noweb --bfile   prune_target_out  --read-genome prune_target_out.genome  --cluster --pca ", length(covariates), " header --out 		ANCESTRY_INFORMATIVE_DIMENSIONS"))
	        system("sed -e '1s/P//g' ANCESTRY_INFORMATIVE_DIMENSIONS.eigenvec > ANCESTRY_INFORMATIVE_DIMENSIONS.readable")
	        pc.file <- "ANCESTRY_INFORMATIVE_DIMENSIONS.readable"
	      }
	    } 
	  }
	}
	
	if(covary & !is.na(user.covariate.file)){
	  pc.file<- user.covariate.file
	}
	
	if(ext.phen & !multiple.target.phenotypes){
	  pheno.data <- read.table(pheno.file, head = F)
	}
	if(ext.phen & multiple.target.phenotypes){
	  pheno.data <- read.table(pheno.file, head = T)
	}
	
	####################################################
	## Extract phenotype data and save it somewhere::
    ####################################################
    if(quantiles & dosage & is.na(ext.phen)){
      cat("WARNING: Cannot produce quantiles plot for dosage data \n unless external phenotype file is supplied \n")
      quantiles <- F
    }
	if(quantiles){
	  if(ext.phen){
	    phen.file.internal <- pheno.data 
	    # debug
	    print(names(phen.file.internal))
	  }
	  if(!ext.phen){
	    if(!dosage){
	      if(!geno.as.list){
            phen.file.internal <- read.table(paste(target, "fam", sep  = "."), head = F)[,c("V2", "V6")] 
	      }
	      if(geno.as.list){
            phen.file.internal <- read.table("flipped_target_1.fam", head = F)[,c("V2", "V6")]     	
	      }
	    }
      }
      names(phen.file.internal) <- c("ID", "PHEN")
      if(levels(as.factor(phen.file.internal$PHEN))[1]=="1"&levels(as.factor(phen.file.internal$PHEN))[2]=="2"){ 
        phen.file.internal$PHEN <- phen.file.internal$PHEN - 1
	  }  
    phen.file.internal <- phen.file.internal[phen.file.internal$PHEN != -9 , ]
    }
	
	if(covary){	
	  pcs <- read.table(pc.file, head = T)
	}
	profile.list <- read.table("profile_list", head = F)
	lists.full <- lists	
		
	## SOMETIMES: an empty list will not print. this section is important - it takes the lists which did print and retains only these from the list of names to read and test model fit etc
	if(!dosage){
	  if(geno.as.list){
	    for(i in 1:22){
	      if(i == 1){
	        lists.temp <- lists[paste(i,profile.list$V1,sep="") %in% list.files(pattern='*PROFILES.*\\.profile'),]
	        lists.out <- lists.temp
	      }
	      if(i > 1){
	        lists.temp <- lists[paste(i,profile.list$V1,sep="") %in% list.files(pattern='*PROFILES.*\\.profile'),]
	        lists.out <- rbind(lists.out, lists.temp)
	      }
	    }
	    lists <- lists[!duplicated(lists$V1),]  
	    profile.list.temp <- data.frame(profile.list, lists.full)
	    names(profile.list.temp) <- c("V1", "V2", "V3","V4")
	    reduced.list <- data.frame(profile.list.temp$V1[profile.list.temp$V2 %in% lists$V1])
	    names(reduced.list) <- c("V1")
	  }
	  if(!geno.as.list){
	    lists <- lists[profile.list$V1 %in% list.files(pattern='PROFILES.*\\.profile'),]
	    reduced.list <- profile.list$V1[profile.list$V1 %in% list.files(pattern='PROFILES.*\\.profile')]
	  }
	}
	if(dosage){
		lists <- lists[profile.list$V1 %in% list.files(pattern='PROFILES.*\\.profile'),]
		reduced.list <- profile.list$V1[profile.list$V1 %in% list.files(pattern='PROFILES.*\\.profile')]
	}
	if(!no.regression){
	  cat(" ################################# \n # \n #   Regression Models \n # \n ################################# \n")
	}

	ranges.list <- read.table("rangelist_ranges",head=F)

	if(!multiple.target.phenotypes){
	  p.out <- as.vector(1)
	  r2.out <- as.vector(1)
    r2full<-as.vector(1);r2null<-as.vector(1)
	  nsnps <- as.vector(1)
	  coefficient <- as.vector(1)
	  s.err <- as.vector(1)
	  prof.temp <- as.vector(1)
	  prev.files <- F
	  abc.r2 <- as.vector(1)
	  ## logistic regression on phenotype
	  if(binary.target){
	    for(i in 1:length(lists$V1)){
	      if(dosage & (!is.na(dos.path.to.lists) | !is.na(dos.list.file))){
	        if(!is.na(dos.list.file)){
	          max.files <- dim(read.table(dos.list.file))[1]
	        }
	        if(!is.na(dos.path.to.lists)){
	          max.files <- 22
	        }
	        for(j in 1:max.files){
	          if(paste(j, profile.list[i,1], sep = "") %in% list.files()){
	            prof.temp <- read.table(paste(j, profile.list[i,1], sep = ""), head = T)[,c("IID","SCORE","PHENO")]
	            names(prof.temp) <- c("temp.IID", "temp.SCORE","temp.PHENO")
	            if(prev.files){
	              prof.temp <- merge(x = prof.temp, y = prof.out, by.x = "temp.IID", by.y = "IID")
	              prof.temp$SCORE <- prof.temp$SCORE + prof.temp$temp.SCORE    
	              prof.temp <- prof.temp[,c("IID", "SCORE","PHENO")]
	            }
	            prof.out <- prof.temp
	            names(prof.out) <- c("IID","SCORE","PHENO")
	            prof <- prof.out
	            prev.files <- T
	          }
	        }
	      }
	      if(geno.as.list){
	         for(j in 1:22){
	          if(paste(j, profile.list[i,1], sep = "") %in% list.files()){
	            prof.temp <- read.table(paste(j, profile.list[i,1], sep = ""), head = T)[,c("IID","SCORE","PHENO")]
	            names(prof.temp) <- c("temp.IID", "temp.SCORE","temp.PHENO")
	            if(prev.files){
	              prof.temp <- merge(x = prof.temp, y = prof.out, by.x = "temp.IID", by.y = "IID")
	              prof.temp$SCORE <- prof.temp$SCORE + prof.temp$temp.SCORE    
	              prof.temp <- prof.temp[,c("IID", "SCORE","PHENO")]
	            }
	            prof.out <- prof.temp
	            names(prof.out) <- c("IID","SCORE","PHENO")
	            prof <- prof.out
	            prev.files <- T         
	          }
	        prev.files <- F
	        }
	      }
	      if((!dosage & !geno.as.list)  | (dosage & is.na(dos.path.to.lists) & is.na(dos.list.file))){			
		    ## Default - no fastscore
		    prof <- read.table(paste(reduced.list[i], sep = ""), head  =T)
	        prof <- prof[,c("IID", "PHENO", "SCORE")]
	      }
	      if(ext.phen){
	        prof <- subset(prof, select=-c(PHENO))
		    prof <- merge(x = prof, by.x = "IID", y = pheno.data, by.y = "V1")
	      }
		  if(!ext.phen) {
	        prof <- prof[prof$PHENO != -9,] 
	        prof <- prof[!is.na(prof$PHENO),]
	        if(levels(as.factor(prof$PHENO))[1]=="1"&levels(as.factor(prof$PHENO))[2]=="2"){ 
		      prof$PHENO <- prof$PHENO - 1
	        }
	        if(!(levels(as.factor(prof$PHENO))[1]=="0"&levels(as.factor(prof$PHENO))[2]=="1")){ 
		      print("ERROR: Unrecognised values for binary phenotype.");quit()
	        }
	      }
	      if(length(levels(as.factor(prof$PHENO)) )< 2 & !ext.phen){
	        cat("ERROR: Phenotype does not have more than one level. Will not perform regression \n")
	        no.regression <- T
	      }
	      if(!no.regression){
	        if(covary & ext.phen){
	          prof <- merge(x = prof, by.x = "IID", y = pcs, by.y = "IID")
	          model.logit <- glm(V2 ~., family="binomial", data = prof[,c("V2","SCORE", covariates)])
	          model.null <-  glm(V2 ~., family="binomial", data = prof[,c("V2", covariates)])
	          p.out[i]  <- summary(model.logit)$coefficients[2,4]
	          if(for.meta){
	            coefficient[i]  <- summary(model.logit)$coefficients[2,1]
	            s.err[i] <- summary(model.logit)$coefficients[2,2]
	          }
		      r2.out[i] <- NagelkerkeR2(model.logit)$R2 - NagelkerkeR2(model.null)$R2  
	        }
	        if(covary & !ext.phen){
	          prof <- merge(x = prof, by.x = "IID", y = pcs, by.y = "IID")
		      model.logit <- glm(PHENO ~ ., family="binomial", data = prof[,c("PHENO", "SCORE", covariates)])
		      model.null <-  glm(PHENO  ~ ., family="binomial", data = prof[,c("PHENO", covariates)])
		      p.out[i]  <- summary(model.logit)$coefficients[2,4]
	          if(for.meta){
	            coefficient[i]  <- summary(model.logit)$coefficients[2,1]
	            s.err[i] <- summary(model.logit)$coefficients[2,2]
	          }
		      r2.out[i] <- NagelkerkeR2(model.logit)$R2 - NagelkerkeR2(model.null)$R2
	        }
	        if(!covary & ext.phen){
	          p.out[i]  <- summary(with(prof, glm(V2  ~ SCORE, family="binomial")))$coefficients[2,4]
	          if(for.meta){
	            coefficient[i]  <- summary(with(prof, glm(V2  ~ SCORE, family="binomial")))$coefficients[2,1]
	            s.err[i] <- summary(with(prof, glm(V2  ~ SCORE, family="binomial")))$coefficients[2,2]
	          }
	          r2.out[i] <- NagelkerkeR2(with(prof, glm(V2 ~ SCORE, family="binomial")))$R2
	        }
	        if(!covary & !ext.phen){
	          p.out[i]  <- summary(with(prof, glm(PHENO  ~ SCORE, family="binomial")))$coefficients[2,4]
	           if(for.meta){
	            coefficient[i]  <- summary(with(prof, glm(PHENO  ~ SCORE, family="binomial")))$coefficients[2,1]
	            s.err[i] <- summary(with(prof, glm(PHENO  ~ SCORE, family="binomial")))$coefficients[2,2]
	          }
	          r2.out[i] <- NagelkerkeR2(with(prof, glm(PHENO ~ SCORE, family="binomial")))$R2
	        }
	      }
	      prof.red <- prof[,c("IID","SCORE")]
	      names(prof.red) <- c("IID", paste("pT", lists$V1[i], sep = "_"))
	      if(report.individual.scores){
	        if(i == 1){
	          prof.all.scores <- prof.red
	        }
	        if(i > 1){
	           prof.all.scores <- merge(x = prof.all.scores, by.x = "IID", y = prof.red, by.y = "IID")    
	        }
	      }
	      nsnps[i] <- dim(ranges.list[ranges.list$V2 < lists$V1[i],])[1]
	      if(i %in% (seq(1, length(lists$V1), 1)[seq(1, length(lists$V1), 1) %in% ceiling(length(lists$V1)*seq(0, 1,0.1))] )){
	        cat(paste("Regression Models: " ,round(100*(i/length(lists$V1)),digits=-1), "% Complete \n",sep=""))
	      }
          if(calculate.abc.r2){
          	abc.r2[i] <- estimateVg2FromP(prevalence1=prevalence1,
          	  prevalence2=prevalence2,
          	  n1=n1, 
          	  sampling1= sampling1, 
          	  n2= n2, 
          	  sampling2= sampling2, 
          	  nsnp= nsnps[i], 
          	  plower = 0, 
          	  pupper = lists$V1[i], 
          	  binary=binary.target, 
          	  p=p.out[i],
          	  nullfraction=pi0
          	)$vg
          }
	    }
	  }	
	  ## linear regression on phenotype
	  if(!binary.target){
	    for(i in 1:length(lists$V1)){
	      if(dosage & (!is.na(dos.path.to.lists) | !is.na(dos.list.file))){
	        for(j in 1:22){
	          prof.temp <- read.table(paste(j, profile.list[j,1], sep = ""), head = T)[,c("IID","SCORE","PHENO")]
	          names(prof.temp) <- c("temp.IID", "temp.SCORE")
	          if(j > 1){
	            prof.temp <- merge(x = prof.temp, y = prof.out, by.x = "temp.IID", by.y = "IID")
	            prof.temp$SCORE <- prof.temp$SCORE + prof.temp$temp.SCORE    
	            prof.temp <- prof.temp[,c("IID", "SCORE","PHENO")]
	          }
	          prof.out <- prof.temp
	          names(prof.out) <- c("IID","SCORE","PHENO")
	        }
	        prof <- prof.out
	      }
	      if(!dosage | (is.na(dos.path.to.lists) & is.na(dos.list.file))){			
	        ## Default - no fastscore
	        prof <- read.table(paste(reduced.list[i], sep = ""), head  =T)
	        prof <- prof[,c("IID", "PHENO", "SCORE")]
	      }
	      if(ext.phen){
	        prof <- subset(prof, select=-c(PHENO))
	        prof <- merge(x = prof, by.x = "IID", y = pheno.data, by.y = "V1")
	      } 
	      if(length(levels(as.factor(prof$PHENO)) )< 2 & !ext.phen){
	        cat("ERROR: Phenotype does not have more than one level. Will not perform regression \n")
	        no.regression <- T
	      }
	      if(!no.regression){
	        if(covary & ext.phen){
	          prof <- merge(x = prof, by.x = "IID", y = pcs, by.y = "IID")
		      model.logit <- glm(V2  ~ ., family="gaussian", data = prof[,c("V2", "SCORE", covariates)])
		      model.null <-  glm(V2  ~ ., family="gaussian", data = prof[,c("V2", covariates)])
		      p.out[i]  <- summary(model.logit)$coefficients[2,4]
	          if(for.meta){
	            coefficient[i]  <- summary(model.logit)$coefficients[2,1]
	            s.err[i] <- summary(model.logit)$coefficients[2,2]
	          }
		      #r2.out[i] <- (1 - (sum( (prof$V2-predict(model.logit))^2 ) / sum( (prof$V2-mean(prof$V2))^2 ) ) )  - (1 - ( sum( (prof$V2-predict(model.null))^2 ) / sum( (prof$V2-mean(prof$V2))^2 ) ) )  
		      r2full[i] <- 1-model.logit$deviance/model.logit$null.deviance
          r2null[i] <- 1-model.null$deviance/model.null$null.deviance
		      r2.out[i] <- r2full[i] - r2null[i]
	     }
		    if(covary & !ext.phen){
	          prof <- merge(x = prof, by.x = "IID", y = pcs, by.y = "IID")
	          model.logit <- glm(PHENO ~ ., family="gaussian", data = prof[,c("PHENO", "SCORE", covariates)])
	          model.null <-  glm(PHENO ~ ., family="gaussian", data = prof[,c("PHENO", covariates)])
	          p.out[i]  <- summary(model.logit)$coefficients[2,4]
	          if(for.meta){
	            coefficient[i]  <- summary(model.logit)$coefficients[2,1]
	            s.err[i] <- summary(model.logit)$coefficients[2,2]
	          }
	         #r2.out[i] <- (1 - ( sum( (prof$PHENO-predict(model.logit))^2 ) / sum( (prof$PHENO-mean(prof$PHENO))^2 ) ) )  - (1 - ( sum( (prof$PHENO-predict(model.null))^2 ) / sum( (prof$PHENO-mean(prof$PHENO))^2 ) ) ) 
	         r2full[i] <- 1-model.logit$deviance/model.logit$null.deviance
	         r2null[i] <- 1-model.null$deviance/model.null$null.deviance
	         r2.out[i] <- r2full[i] - r2null[i]
		    }
	        if(!covary & ext.phen){
	          p.out[i]  <- summary(with(prof, glm(V2  ~ SCORE, family="gaussian")))$coefficients[2,4]
	          if(for.meta){
	            coefficient[i]  <- summary(with(prof, glm(V2  ~ SCORE, family="gaussian")))$coefficients[2,1]
	            s.err[i] <- summary(with(prof, glm(V2  ~ SCORE, family="gaussian")))$coefficients[2,2]
	          }
	          r2.out[i] <-  1 - (sum( (prof$V2-predict(with(prof, glm(V2 ~ SCORE, family="gaussian"))))^2 )/sum( (prof$V2-mean(prof$V2))^2))   			
	        }
	        if(!covary & !ext.phen){
	          p.out[i]  <- summary(with(prof, glm(PHENO  ~ SCORE, family="gaussian")))$coefficients[2,4]
	          if(for.meta){
	            coefficient[i]  <- summary(with(prof, glm(PHENO  ~ SCORE, family="gaussian")))$coefficients[2,1]
	            s.err[i] <- summary(with(prof, glm(PHENO  ~ SCORE, family="gaussian")))$coefficients[2,2]
	          }
	          r2.out[i] <- 1 - (sum( (prof$PHENO-predict(with(prof, glm(PHENO ~ SCORE, family="gaussian"))))^2 )/sum( (prof$PHENO-mean(prof$PHENO))^2 )) 
	        }
	      }
	      prof.red <- prof[,c("IID","SCORE")]
	      names(prof.red) <- c("IID", paste("pT", lists$V1[i], sep = "_"))
	      if(report.individual.scores){
	        if(i == 1){
	          prof.all.scores <- prof.red
	        }
	        if(i > 1){
	          prof.all.scores <- merge(x = prof.all.scores, by.x = "IID", y = prof.red, by.y = "IID")    
	        }
	      }
	      nsnps[i] <- dim(ranges.list[ranges.list$V2 < lists$V1[i],])[1]  
	      if(i %in% (seq(1, length(lists$V1), 1)[seq(1, length(lists$V1), 1) %in% ceiling(length(lists$V1)*seq(0, 1,0.1))] )){
	        cat(paste("Regression Models: " ,round(100*(i/length(lists$V1)),digits=-1), "% Complete \n",sep=""))
	      }
          if(calculate.abc.r2){
          	abc.r2[i] <- estimateVg2FromP(prevalence1=prevalence1,
          	  prevalence2=prevalence2,
          	  n1=n1, 
          	  sampling1= sampling1, 
          	  n2= n2, 
          	  sampling2= sampling2, 
          	  nsnp= nsnps[i], 
          	  plower = 0, 
          	  pupper = lists$V1[i], 
          	  binary=binary.target, 
          	  p=p.out[i],
          	  nullfraction=pi0
          	)$vg
          }
	      
	    }
	  }
	}

	## Alternately, regression with multiple phenotypes::
	if(multiple.target.phenotypes){
	  top.thresh <- as.vector(1)
	  top.thresh.pval <- as.vector(1)
      top.thresh.r2 <- as.vector(1)
	  for(k in 1:length(target.phenotypes)){
	  	cat(paste("Polygenic Scores Regressed on ", target.phenotypes[k] ," Status \n"))
	    binary.target <- target.phenotypes.binary[k]
	    p.out <- as.vector(1)
	    r2.out <- as.vector(1)
	  	r2full<-as.vector(1);r2null<-as.vector(1)
	    nsnps <- as.vector(1)
	    coefficient <- as.vector(1)
	    s.err <- as.vector(1)
	    prof.temp <- as.vector(1)
	    prev.files <- F
	    ## logistic regression on phenotype
	    if(binary.target){
	      for(i in 1:length(lists$V1)){
	        if(dosage & (!is.na(dos.path.to.lists) | !is.na(dos.list.file))){
	          if(!is.na(dos.list.file)){
	            max.files <- dim(read.table(dos.list.file))[1]
	          }
	          if(!is.na(dos.path.to.lists)){
	            max.files <- 22
	          }
	          for(j in 1:max.files){
	            if(paste(j, profile.list[i,1], sep = "") %in% list.files()){
	              prof.temp <- read.table(paste(j, profile.list[i,1], sep = ""), head = T)[,c("IID","SCORE","PHENO")]
	              names(prof.temp) <- c("temp.IID", "temp.SCORE","temp.PHENO")
	              if(prev.files){
	                prof.temp <- merge(x = prof.temp, y = prof.out, by.x = "temp.IID", by.y = "IID")
	                prof.temp$SCORE <- prof.temp$SCORE + prof.temp$temp.SCORE    
	                prof.temp <- prof.temp[,c("IID", "SCORE","PHENO")]
	              }
	              prof.out <- prof.temp
	              names(prof.out) <- c("IID","SCORE","PHENO")
	              prof <- prof.out
	              prev.files <- T
	            }
	          }
	        }
	        if(geno.as.list){
	           for(j in 1:22){
	            if(paste(j, profile.list[i,1], sep = "") %in% list.files()){
	              prof.temp <- read.table(paste(j, profile.list[i,1], sep = ""), head = T)[,c("IID","SCORE","PHENO")]
	              names(prof.temp) <- c("temp.IID", "temp.SCORE","temp.PHENO")
	              if(prev.files){
	                prof.temp <- merge(x = prof.temp, y = prof.out, by.x = "temp.IID", by.y = "IID")
	                prof.temp$SCORE <- prof.temp$SCORE + prof.temp$temp.SCORE    
	                prof.temp <- prof.temp[,c("IID", "SCORE","PHENO")]
	              }
	              prof.out <- prof.temp
	              names(prof.out) <- c("IID","SCORE","PHENO")
	              prof <- prof.out
	              prev.files <- T         
	            }
	          prev.files <- F
	          }
	        }
	        if((!dosage & !geno.as.list)  | (dosage & is.na(dos.path.to.lists) & is.na(dos.list.file))){			
		      ## Default - no fastscore
		      prof <- read.table(paste(reduced.list[i], sep = ""), head  =T)
	          prof <- prof[,c("IID", "PHENO", "SCORE")]
	        }
	        prof <- subset(prof, select=-c(PHENO))
	        temp.pheno.data <- pheno.data[,c("ID", target.phenotypes[k])]
	        names(temp.pheno.data) <- c("V1", "V2")
	        prof <- merge(x = prof, by.x = "IID", y = temp.pheno.data, by.y = "V1")
	        if(!no.regression){
	          if(covary){
	            prof <- merge(x = prof, by.x = "IID", y = pcs, by.y = "IID")
                prof <- prof[prof$V2 != "" ,]
                prof <- prof[prof$V2 != -9 ,]
                prof <- prof[!is.na(prof$V2) ,]
	            model.logit <- glm(V2 ~., family="binomial", data = prof[,c("V2","SCORE", covariates)])
	            model.null <-  glm(V2 ~., family="binomial", data = prof[,c("V2", covariates)])
              
	            # patch: sometimes want to correct for gender and/or age
	           #   ds=prof[,c("V2", "SCORE", covariates)]
	           # 	model.pheno<-glm(V2 ~ Gender+Age, family="binomial",data=ds)
	           # 	model.res<-residuals(model.pheno)
	           # 	cov.new=c("EV1","EV2","EV3","EV4","EV5","EV6","EV7","EV8","EV9","EV10")
	           # 	ds.new=prof[,c("SCORE", cov.new)]
	           # 	ds.new$V2=model.res
	           # 	model.logit <- glm(V2  ~ ., family="gaussian", data = ds.new)
	           # 	model.null <-  glm(V2  ~ ., family="gaussian", data = ds.new[,c("V2", cov.new)])
              
	            p.out[i]  <- summary(model.logit)$coefficients[2,4]
	            if(for.meta){
	              coefficient[i]  <- summary(model.logit)$coefficients[2,1]
	              s.err[i] <- summary(model.logit)$coefficients[2,2]
	            }
	            r2full[i] <- NagelkerkeR2(model.logit)$R2
	            r2null[i] <- NagelkerkeR2(model.null)$R2
	          
		        r2.out[i] <- r2full[i] - r2null[i] 
	          }
	          if(!covary){
                prof <- prof[prof$V2 != "" ,]
                prof <- prof[prof$V2 != -9 ,]
                prof <- prof[!is.na(prof$V2) ,]
	            p.out[i]  <- summary(with(prof, glm(V2  ~ SCORE, family="binomial")))$coefficients[2,4]
	            if(for.meta){
	              coefficient[i]  <- summary(with(prof, glm(V2  ~ SCORE, family="binomial")))$coefficients[2,1]
	              s.err[i] <- summary(with(prof, glm(V2  ~ SCORE, family="binomial")))$coefficients[2,2]
	            }
	            r2.out[i] <- NagelkerkeR2(with(prof, glm(V2 ~ SCORE, family="binomial")))$R2
	          }
	        }
	        prof.red <- prof[,c("IID","SCORE")]
	        names(prof.red) <- c("IID", paste("pT", lists$V1[i], sep = "_"))
	        if(report.individual.scores){
	          if(i == 1){
	            prof.all.scores <- prof.red
	          }
	          if(i > 1){
	             prof.all.scores <- merge(x = prof.all.scores, by.x = "IID", y = prof.red, by.y = "IID")    
	          }
	        }
	        nsnps[i] <- dim(ranges.list[ranges.list$V2 < lists$V1[i],])[1]
	        if(i %in% (seq(1, length(lists$V1), 1)[seq(1, length(lists$V1), 1) %in% ceiling(length(lists$V1)*seq(0, 1,0.1))] )){
	          cat(paste("Regression Models: " ,round(100*(i/length(lists$V1)),digits=-1), "% Complete \n",sep=""))
	        }
          if(calculate.abc.r2){
          	abc.r2[i] <- estimateVg2FromP(prevalence1=prevalence1,
          	  prevalence2=prevalence2,
          	  n1=n1, 
          	  sampling1= sampling1, 
          	  n2= n2, 
          	  sampling2= sampling2, 
          	  nsnp= nsnps[i], 
          	  plower = 0, 
          	  pupper = lists$V1[i], 
          	  binary=binary.target, 
          	  p=p.out[i],
          	  nullfraction=pi0
          	)$vg
          }
	      }
	    }	
	    ## linear regression on phenotype
	    if(!binary.target){
	      for(i in 1:length(lists$V1)){
	        if(dosage & (!is.na(dos.path.to.lists) | !is.na(dos.list.file))){
	          for(j in 1:22){
	            prof.temp <- read.table(paste(j, profile.list[j,1], sep = ""), head = T)[,c("IID","SCORE","PHENO")]
	            names(prof.temp) <- c("temp.IID", "temp.SCORE")
	            if(j > 1){
	              prof.temp <- merge(x = prof.temp, y = prof.out, by.x = "temp.IID", by.y = "IID")
	              prof.temp$SCORE <- prof.temp$SCORE + prof.temp$temp.SCORE    
	              prof.temp <- prof.temp[,c("IID", "SCORE","PHENO")]
	            }
	            prof.out <- prof.temp
	            names(prof.out) <- c("IID","SCORE","PHENO")
	          }
	          prof <- prof.out
	        }
	        if(!dosage | (is.na(dos.path.to.lists) & is.na(dos.list.file))){			
	          ## Default - no fastscore
	          prof <- read.table(paste(reduced.list[i], sep = ""), head  =T)
	          prof <- prof[,c("IID", "PHENO", "SCORE")]
	        }
	        prof <- subset(prof, select=-c(PHENO))
	        temp.pheno.data <- pheno.data[,c("ID", target.phenotypes[k])]
	        names(temp.pheno.data) <- c("V1", "V2")
	        prof <- merge(x = prof, by.x = "IID", y = temp.pheno.data, by.y = "V1")
	        if(!no.regression){
	          if(covary){
	            prof <- merge(x = prof, by.x = "IID", y = pcs, by.y = "IID")
                prof <- prof[prof$V2 != "" ,]
                prof <- prof[!is.na(prof$V2) ,]
		        model.logit <- glm(V2  ~ ., family="gaussian", data = prof[,c("V2", "SCORE", covariates)])
		        model.null <-  glm(V2  ~ ., family="gaussian", data = prof[,c("V2", covariates)])
 
  # patch: sometimes want to correct for gender and/or age
#	ds=prof[,c("V2", "SCORE", covariates)]
#	model.pheno<-glm(V2 ~ Gender+Age, family="gaussian",data=ds)
#	model.res<-residuals(model.pheno)
#	cov.new=c("EV1","EV2","EV3","EV4","EV5","EV6","EV7","EV8","EV9","EV10")
#	ds.new=prof[,c("SCORE", cov.new)]
#	ds.new$V2=model.res
#	model.logit <- glm(V2  ~ ., family="gaussian", data = ds.new)
#	model.null <-  glm(V2  ~ ., family="gaussian", data = ds.new[,c("V2", cov.new)])

  	        p.out[i]  <- summary(model.logit)$coefficients[2,4]
	            if(for.meta){
	              coefficient[i]  <- summary(model.logit)$coefficients[2,1]
	              s.err[i] <- summary(model.logit)$coefficients[2,2]
	            }
		        #r2.out[i] <- (1 - (sum( (prof$V2-predict(model.logit))^2 ) / sum( (prof$V2-mean(prof$V2))^2 ) ) )  - (1 - ( sum( (prof$V2-predict(model.null))^2 ) / sum( (prof$V2-mean(prof$V2))^2 ) ) )  
		        r2full[i] <- 1-model.logit$deviance/model.logit$null.deviance
		        r2null[i] <- 1-model.null$deviance/model.null$null.deviance
		        r2.out[i] <- r2full[i] - r2null[i]
            # above two are the same
	       }
	          if(!covary){
                prof <- prof[prof$V2 != "" ,]
                prof <- prof[!is.na(prof$V2) ,]
	            p.out[i]  <- summary(with(prof, glm(V2  ~ SCORE, family="gaussian")))$coefficients[2,4]
	            if(for.meta){
	              coefficient[i]  <- summary(with(prof, glm(V2  ~ SCORE, family="gaussian")))$coefficients[2,1]
	              s.err[i] <- summary(with(prof, glm(V2  ~ SCORE, family="gaussian")))$coefficients[2,2]
	            }
	            #r2.out[i] <-  1 - (sum( (prof$V2-predict(with(prof, glm(V2 ~ SCORE, family="gaussian"))))^2 )/sum( (prof$V2-mean(prof$V2))^2)) 
	            r2full[i] <- 1-model.logit$deviance/model.logit$null.deviance
	            r2null[i] <- 1-model.null$deviance/model.null$null.deviance
	            r2.out[i] <- r2full[i] - r2null[i]
	          }
	        }
	        prof.red <- prof[,c("IID","SCORE")]
	        names(prof.red) <- c("IID", paste("pT", lists$V1[i], sep = "_"))
	        if(report.individual.scores){
	          if(i == 1){
	            prof.all.scores <- prof.red
	          }
	          if(i > 1){
	            prof.all.scores <- merge(x = prof.all.scores, by.x = "IID", y = prof.red, by.y = "IID")    
	          }
	        }
	        nsnps[i] <- dim(ranges.list[ranges.list$V2 < lists$V1[i],])[1]  
	        if(i %in% (seq(1, length(lists$V1), 1)[seq(1, length(lists$V1), 1) %in% ceiling(length(lists$V1)*seq(0, 1,0.1))] )){
	          cat(paste("Regression Models: " ,round(100*(i/length(lists$V1)),digits=-1), "% Complete \n",sep=""))
	        }
          if(calculate.abc.r2){
          	abc.r2[i] <- estimateVg2FromP(prevalence1=prevalence1,
          	  prevalence2=prevalence2,
          	  n1=n1, 
          	  sampling1= sampling1, 
          	  n2= n2, 
          	  sampling2= sampling2, 
          	  nsnp= nsnps[i], 
          	  plower = 0, 
          	  pupper = lists$V1[i], 
          	  binary=binary.target, 
          	  p=p.out[i],
          	  nullfraction=pi0
          	)$vg
          }
	      }
	    }
	    top.thresh[k] <- lists$V1[which.min(p.out)]
  	    top.thresh.pval[k] <- p.out[which.min(p.out)]
  	    top.thresh.r2[k] <- r2.out[which.min(p.out)]
	    thresh <- lists$V1 
        if(!calculate.abc.r2){
	      if(!for.meta){
	        output <- data.frame(thresh, p.out, r2.out, nsnps)
	      }
	      if(for.meta){
	        output <- data.frame(thresh, p.out, r2.out, nsnps, coefficient, s.err)
	      }
	    }
        if(calculate.abc.r2){
	      if(!for.meta){
	        output <- data.frame(thresh, p.out, abc.r2, nsnps)
	      }
	      if(for.meta){
	        output <- data.frame(thresh, p.out,abc.r2, nsnps, coefficient, s.err)
	      }
	    }
	    if(!multiple.base.phenotypes){
  	      write.table(output, paste(figname, target.phenotypes[k], "RAW_RESULTS_DATA.txt", sep = "_"), col.names = T, row.names = F, quote  =F)
	    }
	    if(multiple.base.phenotypes){
  	      write.table(output, paste(figname, base.phenotypes.names[basePhen], "PREDICTING", target.phenotypes[k], "RAW_RESULTS_DATA.txt", sep = "_"), col.names = T, row.names = F, quote  =F)
	    }
	    if(report.individual.scores){
	      if(!report.best.score.only & !multiple.base.phenotypes){
	        write.table(prof.all.scores, paste(figname, target.phenotypes[k], "SCORES_AT_ALL_THRESHOLDS.txt", sep = "_"), col.names = T, row.names = F, quote = F)
	      }
	      if(!report.best.score.only & multiple.base.phenotypes){
	        write.table(prof.all.scores, paste(figname, base.phenotypes.names[basePhen], "PREDICTING", target.phenotypes[k],   "SCORES_AT_ALL_THRESHOLDS.txt", sep = "_"), col.names = T, row.names = F, quote = F)
	      }
	      if(report.best.score.only & !multiple.base.phenotypes){
	        write.table(prof.all.scores[,c("IID", paste("pT", output$thresh[which.min(output$p.out)],sep="_"))], paste(figname, target.phenotypes[k], "SCORES_AT_BEST-FIT-PRS.txt", sep = "_"), col.names = T, row.names = F, quote = F)
	      }
	      if(report.best.score.only & multiple.base.phenotypes){
	        write.table(prof.all.scores[,c("IID", paste("pT", output$thresh[which.min(output$p.out)],sep="_"))], paste(figname,base.phenotypes.names[basePhen], "PREDICTING", target.phenotypes[k], "SCORES_AT_BEST-FIT-PRS.txt", sep = "_"), col.names = T, row.names = F, quote = F)
	      }
	    }
	  }
	  ## end of loop through target phenotypes
	if(!multiple.base.phenotypes){  
    	write.table(data.frame(target.phenotypes, target.phenotypes.binary, top.thresh, top.thresh.pval, top.thresh.r2), paste(figname, "TOP_SCORES_ACROSS_PHENOTYPES.txt",sep="_"),
	  col.names=T, row.names=F, quote=F)
    }	
	if(multiple.base.phenotypes){  
    	write.table(data.frame(target.phenotypes, target.phenotypes.binary, top.thresh, top.thresh.pval, top.thresh.r2), paste(figname, base.phenotypes.names[basePhen], "PREDICTING_TOP_SCORES_ACROSS_TARGET_PHENOTYPES.txt",sep="_"),
	  col.names=T, row.names=F, quote=F)
    }	

  }
  

	if(no.regression & !multiple.target.phenotypes){
	  if(report.individual.scores & !multiple.base.phenotypes){
	      write.table(prof.all.scores, paste(figname, "SCORES_AT_ALL_THRESHOLDS.txt", sep = "_"), col.names = T, row.names = F, quote = F)
	  }
	  if(report.individual.scores & multiple.base.phenotypes){
	      write.table(prof.all.scores, paste(figname,base.phenotypes.names[basePhen], "SCORES_AT_ALL_THRESHOLDS.txt", sep = "_"), col.names = T, row.names = F, quote = F)
	  }
	  if(cleanup){
	    cat(" ################################# \n # \n #   Cleanup \n # \n ################################# \n")
	    system("rm PROFILES.0.5*.profile", ignore.stdout=T,ignore.stderr=T)	
	    system("rm PROFILES.0.4*.profile", ignore.stdout=T,ignore.stderr=T)
	    system("rm PROFILES.0.3*.profile", ignore.stdout=T,ignore.stderr=T)
	    system("rm PROFILES.0.2*.profile", ignore.stdout=T,ignore.stderr=T)
	    system("rm PROFILES.0.1*.profile", ignore.stdout=T,ignore.stderr=T)
	    system("rm PROFILES.0.0*.profile", ignore.stdout=T,ignore.stderr=T)
	    system("rm PROFILES.*.profile", ignore.stdout=T,ignore.stderr=T)
	    system("rm PROFILES.S1*.profile", ignore.stdout=T,ignore.stderr=T)	
	    system("rm PROFILES.S2*.profile", ignore.stdout=T,ignore.stderr=T)
	    system("rm PROFILES.S3*.profile", ignore.stdout=T,ignore.stderr=T)
	    system("rm PROFILES.S4*.profile", ignore.stdout=T,ignore.stderr=T)
	    system("rm PROFILES.S5*.profile", ignore.stdout=T,ignore.stderr=T)
	    system("rm PROFILES.S6*.profile", ignore.stdout=T,ignore.stderr=T)
	    system("rm PROFILES.S7*.profile", ignore.stdout=T,ignore.stderr=T)
	    system("rm PROFILES.S8*.profile", ignore.stdout=T,ignore.stderr=T)
	    system("rm PROFILES.S9*.profile", ignore.stdout=T,ignore.stderr=T)
	    system("rm *PROFILES.S*.profile", ignore.stdout=T,ignore.stderr=T)
	    system("rm *PROFILES.0.5*.profile", ignore.stdout=T,ignore.stderr=T)	
	    system("rm *PROFILES.0.4*.profile", ignore.stdout=T,ignore.stderr=T)
	    system("rm *PROFILES.0.3*.profile", ignore.stdout=T,ignore.stderr=T)
	    system("rm *PROFILES.0.2*.profile", ignore.stdout=T,ignore.stderr=T)
	    system("rm *PROFILES.0.1*.profile", ignore.stdout=T,ignore.stderr=T)
	    system("rm *PROFILES.0.0*.profile", ignore.stdout=T,ignore.stderr=T)
	    system("rm *PROFILES.*.profile", ignore.stdout=T,ignore.stderr=T)
	    system("rm *PROFILES.S1*.profile", ignore.stdout=T,ignore.stderr=T)	
	    system("rm *PROFILES.S2*.profile", ignore.stdout=T,ignore.stderr=T)
	    system("rm *PROFILES.S3*.profile", ignore.stdout=T,ignore.stderr=T)
	    system("rm *PROFILES.S4*.profile", ignore.stdout=T,ignore.stderr=T)
	    system("rm *PROFILES.S5*.profile", ignore.stdout=T,ignore.stderr=T)
	    system("rm *PROFILES.S6*.profile", ignore.stdout=T,ignore.stderr=T) 
	    system("rm *PROFILES.S7*.profile", ignore.stdout=T,ignore.stderr=T)
	    system("rm *PROFILES.S8*.profile", ignore.stdout=T,ignore.stderr=T)
	    system("rm *PROFILES.S9*.profile", ignore.stdout=T,ignore.stderr=T)
	    system("rm *PROFILES.S*.profile", ignore.stdout=T,ignore.stderr=T) 
	    system("rm flip*", ignore.stdout=T,ignore.stderr=T)
	    system("rm clean*", ignore.stdout=T,ignore.stderr=T)
	    system("rm profile_list", ignore.stdout=T,ignore.stderr=T)
	    system("rm rangelist.txt", ignore.stdout=T,ignore.stderr=T)
	    system("rm rangelist_ranges", ignore.stdout=T,ignore.stderr=T)
	    system("rm rawfile.raw", ignore.stdout=T,ignore.stderr=T)
	    system("rm reordered_base", ignore.stdout=T,ignore.stderr=T)
	    system("rm target_no_mhc*", ignore.stdout=T,ignore.stderr=T)
	    system("rm ANCESTRY_INFORMATIVE_DIMENSIONS.cluster*", ignore.stdout=T,ignore.stderr=T)
	    system("rm reordered_base", ignore.stdout=T,ignore.stderr=T)
	    system("rm ANCESTRY_INFORMATIVE_DIMENSIONS.log", ignore.stdout=T,ignore.stderr=T)
	    system("rm head_disc", ignore.stdout=T,ignore.stderr=T)
	    system("rm LE_SNPs", ignore.stdout=T,ignore.stderr=T)
	    system("rm prune_target*", ignore.stdout=T,ignore.stderr=T)
	    system("rm plink.log", ignore.stdout=T,ignore.stderr=T)
	    system("rm HEADER", ignore.stdout=T,ignore.stderr=T)
	    system("rm PROFILES.hh", ignore.stdout=T,ignore.stderr=T)
	    system("rm PROFILES.log", ignore.stdout=T,ignore.stderr=T)
	    system("rm synonymous_snps*", ignore.stdout=T,ignore.stderr=T)
	    system("rm synonymous_snps", ignore.stdout=T,ignore.stderr=T) 
	    system("rm non_synonymous_snps_only*", ignore.stdout=T,ignore.stderr=T)
	    system("rm mhc.txt", ignore.stdout=T,ignore.stderr=T)
	    system("rm Rplots.pdf", ignore.stdout=T,ignore.stderr=T)
	    system("rm temp.raw", ignore.stdout=T,ignore.stderr=T)
	    system("rm TARGET_SNPs", ignore.stdout=T,ignore.stderr=T)
	    system("rm base_SNPS", ignore.stdout=T,ignore.stderr=T)  
	  }
	  if(print.time){
	    cat(" ################################# \n # \n #   Print time \n # \n ################################# \n")
	  }
	  running.time <- proc.time()[3] - start.time 
	  if(running.time < 60){
	    out.run.time <- round(running.time, digits = 2)
	    out.time.units <- "seconds"
	  }
	  if(running.time > 60 & running.time < 3600){
	    out.run.time <- round((running.time/60), digits = 2)
	    out.time.units <- "minutes"
	  }
	  if(running.time > 3600){
	    out.run.time <- round((running.time/3600), digits = 2)
	    out.time.units <- "hours"
	  }
	  if(print.time){
	    print(paste("RUNNING TIME: ", 	out.run.time, out.time.units, sep = " "))
	  }
	  quit()
	}
	
	if(!multiple.target.phenotypes){
	  thresh <- lists$V1 
        if(!calculate.abc.r2){
	      if(!for.meta){
	        output <- data.frame(thresh, p.out, r2.out, nsnps)
	      }
	      if(for.meta){
	        output <- data.frame(thresh, p.out, r2.out, nsnps, coefficient, s.err)
	      }
	    }
        if(calculate.abc.r2){
	      if(!for.meta){
	        output <- data.frame(thresh, p.out, abc.r2, nsnps)
	      }
	      if(for.meta){
	        output <- data.frame(thresh, p.out,abc.r2, nsnps, coefficient, s.err)
	      }
	    }
	  if(!multiple.base.phenotypes){
  	    write.table(output, paste(figname, "RAW_RESULTS_DATA.txt", sep = "_"), col.names = T, row.names = F, quote  =F)
  	  }
	  if(multiple.base.phenotypes){
  	    write.table(output, paste(figname,base.phenotypes.names[basePhen], "RAW_RESULTS_DATA.txt", sep = "_"), col.names = T, row.names = F, quote  =F)
  	  }
  	  
	  if(report.individual.scores){
	    if(!report.best.score.only & !multiple.base.phenotypes){
	      write.table(prof.all.scores, paste(figname, "SCORES_AT_ALL_THRESHOLDS.txt", sep = "_"), col.names = T, row.names = F, quote = F)
	    }
	    if(!report.best.score.only & multiple.base.phenotypes){
	      write.table(prof.all.scores, paste(figname,base.phenotypes.names[basePhen],"SCORES_AT_ALL_THRESHOLDS.txt", sep = "_"), col.names = T, row.names = F, quote = F)
	    }
	    if(report.best.score.only & !multiple.base.phenotypes){
	      write.table(prof.all.scores[,c("IID", paste("pT", output$thresh[which.min(output$p.out)],sep="_"))], paste(figname, "SCORES_AT_BEST-FIT-PRS.txt", sep = "_"), col.names = T, row.names = F, quote = F)
	    }
	    if(report.best.score.only & multiple.base.phenotypes){
	      write.table(prof.all.scores[,c("IID", paste("pT", output$thresh[which.min(output$p.out)],sep="_"))], paste(figname,base.phenotypes.names[basePhen], "SCORES_AT_BEST-FIT-PRS.txt", sep = "_"), col.names = T, row.names = F, quote = F)
	    }
	  }
	}
    
	if(!multiple.base.phenotypes){
	  if(!multiple.target.phenotypes){
            if(quantiles & ggfig){
    cat(" ################################# \n # \n #   Quantiles Plots \n # \n ################################# \n")
      if(binary.target){
        if(report.individual.scores & !report.best.score.only){
          scores.internal <- prof.all.scores[,c("IID", paste("pT", output$thresh[which.min(output$p.out)],sep="_"))]
        }
        if(report.individual.scores & report.best.score.only){
          scores.internal <- read.table(paste(figname, "SCORES_AT_BEST-FIT-PRS.txt", sep = "_"), head = T)
        }
        names(scores.internal) <- c("ID", "SCORE")
        for.quantiles <- merge(x = phen.file.internal, by.x = "ID", y = scores.internal, by.y = "ID")
        for.quantiles$quantiles <- as.numeric(cut(for.quantiles$SCORE, breaks = quantile(for.quantiles$SCORE, probs = seq(0, 1, 1/num.quantiles)), include.lowest=T))
        for.quantiles$quantiles <- factor(for.quantiles$quantiles, levels = c(quant.ref, seq(1, num.quantiles, 1)[seq(1, num.quantiles, 1) != quant.ref]))

        if(covary == F){
          or.quantiles <- exp(summary(glm(PHEN ~ quantiles, family="binomial", data = for.quantiles))$coefficients[1:num.quantiles,1])
          ci.quantiles.u <- exp(summary(glm(PHEN ~ quantiles, family="binomial", data = for.quantiles))$coefficients[1:num.quantiles,1] + (1.96*summary(glm(PHEN ~ quantiles, family="binomial", data = for.quantiles))$coefficients[1:num.quantiles,2]))
          ci.quantiles.l <- exp(summary(glm(PHEN ~ quantiles, family="binomial", data = for.quantiles))$coefficients[1:num.quantiles,1] - (1.96*summary(glm(PHEN ~ quantiles, family="binomial", data = for.quantiles))$coefficients[1:num.quantiles,2]))
        }
        if(covary == T){
          for.quantiles <- merge(x = for.quantiles, by.x = "ID", y = pcs, by.y = "IID")

          or.quantiles <- exp(summary(glm(PHEN ~ ., family="binomial", data = for.quantiles[,c("PHEN","quantiles",covariates)]))$coefficients[1:num.quantiles,1])
          ci.quantiles.u <- exp(summary(glm(PHEN ~ ., family="binomial", data = for.quantiles[,c("PHEN","quantiles",covariates)]))$coefficients[1:num.quantiles,1] + (1.96*summary(glm(PHEN ~ quantiles, family="binomial", data = for.quantiles))$coefficients[1:num.quantiles,2]))
          ci.quantiles.l <- exp(summary(glm(PHEN ~ ., family="binomial", data = for.quantiles[,c("PHEN","quantiles",covariates)]))$coefficients[1:num.quantiles,1] - (1.96*summary(glm(PHEN ~ quantiles, family="binomial", data = for.quantiles))$coefficients[1:num.quantiles,2]))
        }
        or.quantiles[1] <- 1 
        ci.quantiles.u[1] <- 1
        ci.quantiles.l[1] <- 1
        quant.list <- seq(1, num.quantiles, 1)
        quant.list <- quant.list[quant.list != quant.ref]     
        quantiles.for.table <- c(quant.ref, quant.list)
        quantiles.df <- data.frame(or.quantiles, ci.quantiles.u, ci.quantiles.l, quantiles.for.table)
        names(quantiles.df) <- c("OR", "CI.U", "CI.L", "DEC")
        quantiles.df <- quantiles.df[order(quantiles.df$DEC),]
        quantiles.plot <- ggplot(quantiles.df) + 
          geom_point(aes(x = DEC, y = OR), colour = "royalblue2", size=4) + 
          theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(), panel.background = element_blank(), 
          axis.line = element_line(colour = "black",size=0.5)) + 
          geom_pointrange(aes(ymin = CI.L,ymax = CI.U, y = OR, x = DEC), colour = "royalblue2", size = 0.9) + 
          ylab("Odds Ratio for Score on Phenotype") + 
          xlab("Quantiles for Polygenic Score") + 
          scale_x_continuous(breaks=seq(0, num.quantiles, 1))   
        }  
      if(!binary.target){
        if(report.individual.scores & !report.best.score.only){
          scores.internal <- prof.all.scores[,c("IID", paste("pT", output$thresh[which.min(output$p.out)],sep="_"))]
        }
        if(report.individual.scores & report.best.score.only){
          scores.internal <- read.table(paste(figname, "SCORES_AT_BEST-FIT-PRS.txt", sep = "_"), head = T)
        }
        names(scores.internal) <- c("ID", "SCORE")
        for.quantiles <- merge(x = phen.file.internal, by.x = "ID", y = scores.internal, by.y = "ID")
        for.quantiles$quantiles <- as.numeric(cut(for.quantiles$SCORE, breaks = quantile(for.quantiles$SCORE, probs = seq(0, 1, 1/num.quantiles)), include.lowest=T))
        for.quantiles$quantiles <- factor(for.quantiles$quantiles, levels = c(quant.ref, seq(1, num.quantiles, 1)[seq(1, num.quantiles, 1) != quant.ref]))

        if(covary == F){
          coef.quantiles <- summary(glm(PHEN ~ quantiles, family="gaussian", data = for.quantiles))$coefficients[1:num.quantiles,1]
          ci.quantiles.u <- summary(glm(PHEN ~ quantiles, family="gaussian", data = for.quantiles))$coefficients[1:num.quantiles,1] + (1.96*summary(glm(PHEN ~ quantiles, family="gaussian", data = for.quantiles))$coefficients[1:num.quantiles,2])
          ci.quantiles.l <- summary(glm(PHEN ~ quantiles, family="gaussian", data = for.quantiles))$coefficients[1:num.quantiles,1] - (1.96*summary(glm(PHEN ~ quantiles, family="gaussian", data = for.quantiles))$coefficients[1:num.quantiles,2])
        }
        if(covary == T){
          for.quantiles <- merge(x = for.quantiles, by.x = "ID", y = pcs, by.y = "IID")
          coef.quantiles <- summary(glm(PHEN ~ ., family="gaussian", data = for.quantiles[,c("PHEN","quantiles",covariates)]))$coefficients[1:num.quantiles,1]
          ci.quantiles.u <- summary(glm(PHEN ~ ., family="gaussian", data = for.quantiles[,c("PHEN","quantiles",covariates)]))$coefficients[1:num.quantiles,1] + (1.96*summary(glm(PHEN ~ ., family="gaussian", data = for.quantiles[,c("PHEN","quantiles",covariates)]))$coefficients[1:num.quantiles,2])
          ci.quantiles.l <- summary(glm(PHEN ~ ., family="gaussian", data = for.quantiles[,c("PHEN","quantiles",covariates)]))$coefficients[1:num.quantiles,1] - (1.96*summary(glm(PHEN ~ ., family="gaussian", data = for.quantiles[,c("PHEN","quantiles",covariates)]))$coefficients[1:num.quantiles,2])
        }
        coef.quantiles[1] <- 0 
        ci.quantiles.u[1] <- 0
        ci.quantiles.l[1] <- 0
        quant.list <- seq(1, num.quantiles, 1)
        quant.list <- quant.list[quant.list != quant.ref]     
        quantiles.for.table <- c(quant.ref, quant.list)
        quantiles.df <- data.frame(coef.quantiles, ci.quantiles.u, ci.quantiles.l, quantiles.for.table)
        names(quantiles.df) <- c("Coef", "CI.U", "CI.L", "DEC")
        quantiles.df <- quantiles.df[order(quantiles.df$DEC),]
        quantiles.plot <- ggplot(quantiles.df) + 
          geom_point(aes(x = DEC, y = Coef), colour = "royalblue2", size=4) + 
          theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(), panel.background = element_blank(), 
          axis.line = element_line(colour = "black",size=0.5)) + 
          geom_pointrange(aes(ymin = CI.L,ymax = CI.U, y = Coef, x = DEC), colour = "royalblue2", size = 0.9) + 
          ylab("Change in Phenotype given score in quantiles") + 
          xlab("Quantiles for Polygenic Score") + 
          scale_x_continuous(breaks=seq(0, num.quantiles, 1))   
        }  
      quantiles.plot
      ggsave(filename=paste(figname, "QUANTILES_PLOT.png", sep = "_"),width=7,height=7)  
    }	

	    if(!fastscore){
	      cat(" ################################# \n # \n #   High Density Plots \n # \n ################################# \n")
	    }
#	    barchart.levels.old <- barchart.levels
	    if(best.thresh.on.bar){
	      #add max thresh
	      barchart.levels <- c(barchart.levels.old, output$thresh[which.min(output$p.out)])
	      barchart.levels <- barchart.levels[!duplicated(barchart.levels)]	
	      barchart.levels <- sort(barchart.levels, decreasing = F)
	    }
	    if(!fastscore & !ggfig){
	      if(!scatter.R2){
	        png(paste(figname,"_HIGH-RES_PLOT_", Sys.Date(), ".png", sep = ""))
	          with(output, plot(	x = output$thresh,
	          y = -log10(output$p.out),
	          xlab = expression(italic(P)-value~threshold),
	          ylab = bquote(-log[10]~model~italic(P)-value),
	          pch = 19,
	          type = "b"
	        ))
	        with(output[output$thresh %in% barchart.levels.old,], points(thresh, -log10(p.out), col = "green", pch = 19, type= "b"))
	      }
	      if(scatter.R2){
	        png(paste(figname,"_HIGH-RES_PLOT_", Sys.Date(), ".png", sep = ""))
	        with(output, plot(	x = output$thresh,
	          y = r2.out,
	          xlab = expression(italic(P)-value~threshold),
	          ylab = bquote(-log[10]~model~italic(P)-value),
	          pch = 19,
	          type = "b"
	        ))
	        with(output[output$thresh %in% barchart.levels.old,], points(thresh, r2.out, col = "green", pch = 19, type= "b"))
	      }	
	      dev.off()
	    }
	    if(!fastscore & ggfig){
	      if(!scatter.R2){
	        ggfig.points <- ggplot(data = output) + geom_point(aes(x = thresh, y = -log10(p.out))) +
	          geom_line(aes(x = thresh, y = -log10(p.out))) +
              xlab(expression(italic(P)-value~threshold~(italic(P)[T])))  +
              ylab(bquote(PRS~model~fit:~italic(P)-value~(-log[10]))) +
	          geom_line(aes(thresh,  -log10(p.out)), colour = "green",
	          subset = .(thresh %in% barchart.levels.old))  +
	          geom_point(aes(thresh,  -log10(p.out)), colour = "green",
	          subset = .(thresh %in% barchart.levels.old))  + 
	          theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(), panel.background = element_blank(), 
          axis.line = element_line(colour = "black",size=0.5))
	        ggfig.points
	        ggsave(filename=paste(figname,"_HIGH-RES_PLOT_", Sys.Date(), ".png", sep = ""),width=7,height=7)
	      }
	      if(scatter.R2){
	        ggfig.points <- ggplot(data = output) + geom_point(aes(x = thresh, y = r2.out)) +
	          geom_line(aes(x = thresh, y = r2.out)) +
              xlab(expression(italic(P)-value~threshold~(italic(P)[T])))  +
	          geom_line(aes(thresh,  r2.out), colour = "green",
	          subset = .(thresh %in% barchart.levels.old))  +
	          geom_point(aes(thresh,  r2.out), colour = "green",
	          subset = .(thresh %in% barchart.levels.old))  + 
	          theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(), panel.background = element_blank(), 
          axis.line = element_line(colour = "black",size=0.5))
   if(binary.target){
   	  if(!calculate.abc.r2){
        ggfig.points <- ggfig.points + ylab(expression(paste("PRS model fit:  ", R^2, " (Nagelkerke's)", sep = " ")))
      }
   	  if(calculate.abc.r2){
        ggfig.points <- ggfig.points + ylab(expression(paste("PRS model fit:  ", R^2, " (Liability Scale)", sep = " ")))
      }
    }
    if(!binary.target){
      ggfig.points <- ggfig.points + ylab(paste("PRS model fit:  ",expression(paste(R^2))))
    }
	        ggfig.points
	        ggsave(filename=paste(figname,"_HIGH-RES_PLOT_", Sys.Date(), ".png", sep = ""),width=7,height=7)
	      }	 
	    }
	    options(scipen=0,digits=7)
	    output <- read.table(paste(figname, "RAW_RESULTS_DATA.txt", sep = "_"),head=T)
	    output <- output[,c(1:3)]
	    names(output) <- c("thresh", "p.out", "r2.out")
	    output <- output[output$thresh %in% barchart.levels,]
	    output$print.p[round(output$p.out, digits = 3) != 0] <- round(output$p.out[round(output$p.out, digits = 3) != 0], digits = 3)
	    output$print.p[round(output$p.out, digits = 3) == 0] <- format(output$p.out[round(output$p.out, digits = 3) == 0], digits=2)
	    output$print.p <- sub("e", "*x*10^", output$print.p)
	    cat(" ################################# \n # \n #   Barplots \n # Bars for inclusion can be changed using the barchart.levels option \n ################################# \n")
	    if(ggfig){
	      if(!bar.col.is.pval){
	        ggfig.plot <-   ggplot(data = output)  + geom_bar(aes(x = factor(thresh), y = r2.out, fill = factor(thresh)), stat="identity") +     scale_fill_brewer(palette=barpalatte, name = expression(italic(P)-value~threshold)) + xlab(expression(italic(P)-value~threshold~(italic(P)[T])))   
	            }
	      if(bar.col.is.pval){
	        ggfig.plot <-   ggplot(data = output)  + geom_bar(aes(x = factor(thresh), y = r2.out, fill = -log10(p.out)), stat="identity") +     scale_fill_gradient(low= bar.col.is.pval.lowcol, high= bar.col.is.pval.highcol, name =bquote(atop(-log[10]~model,italic(P)-value),)) +  xlab(expression(italic(P)-value~threshold~(italic(P)[T]))) 	      }
	      ggfig.plot <- ggfig.plot + geom_text(aes(x = factor(thresh), y = r2.out, label = paste(print.p)), vjust = -1.5, hjust = 0, angle = 45, cex = 2.8, parse=T) +
	        scale_y_continuous(limits = c(0, max(output$r2.out)*1.25)) +
	         theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(), panel.background = element_blank(), 
          axis.line = element_line(colour = "black",size=0.5))
   if(binary.target){
   	  if(!calculate.abc.r2){
        ggfig.plot <- ggfig.plot + ylab(expression(paste("PRS model fit:  ", R^2, " (Nagelkerke's)", sep = " ")))
      }
   	  if(calculate.abc.r2){
        ggfig.plot <- ggfig.plot + ylab(expression(paste("PRS model fit:  ", R^2, " (Liability Scale)", sep = " ")))
      }
    }
            if(!binary.target){
              ggfig.plot <- ggfig.plot + ylab(expression(paste("PRS model fit:  ",R^2)))
            }
	        ggfig.plot
	        ggsave(filename=paste(figname, "_BARPLOT_", Sys.Date(), ".png", sep = ""),width=7,height=7)		
	    }  
	    if(!ggfig) {
	      png(paste(figname, "_BARPLOT_", Sys.Date(), ".png", sep = ""))
	      plot.fig  <- with(output, barplot(r2.out, 
	        names = thresh,
	        main = "", 
	        col = "red",
	        xlab = expression(italic(P)[T]),
	        ylab = expression(R^2),
	        ylim = c(0,  max(output$r2.out)*1.25)  ))
	        text( parse(text=paste(
	        output$print.p)), 
	        x = plot.fig+0.1, 
	        y =  output$r2.out+ (max(output$r2.out)*1.125-max(output$r2.out)), 
	        srt = 45)
	      dev.off()
	    }
	  }
	  if(multiple.target.phenotypes){
	    for(k in 1:length(target.phenotypes)){
	      output <- read.table(paste(figname, target.phenotypes[k], "RAW_RESULTS_DATA.txt", sep = "_"), head=T)
	      
	          if(quantiles & ggfig){
    cat(" ################################# \n # \n #   Quantiles Plots \n # \n ################################# \n")
      if(binary.target){
        if(report.individual.scores & !report.best.score.only){
          scores.internal <- prof.all.scores[,c("IID", paste("pT", output$thresh[which.min(output$p.out)],sep="_"))]
        }
        if(report.individual.scores & report.best.score.only){
          scores.internal <- read.table(paste(figname,target.phenotypes[k], "SCORES_AT_BEST-FIT-PRS.txt", sep = "_"), head = T)
        }
        names(scores.internal) <- c("ID", "SCORE")
                for.quantiles <- merge(x = phen.file.internal, by.x = "ID", y = scores.internal, by.y = "ID")
        for.quantiles$quantiles <- as.numeric(cut(for.quantiles$SCORE, breaks = quantile(for.quantiles$SCORE, probs = seq(0, 1, 1/num.quantiles)), include.lowest=T))
        for.quantiles$quantiles <- factor(for.quantiles$quantiles, levels = c(quant.ref, seq(1, num.quantiles, 1)[seq(1, num.quantiles, 1) != quant.ref]))

        if(covary == F){
          or.quantiles <- exp(summary(glm(PHEN ~ quantiles, family="binomial", data = for.quantiles))$coefficients[1:num.quantiles,1])
          ci.quantiles.u <- exp(summary(glm(PHEN ~ quantiles, family="binomial", data = for.quantiles))$coefficients[1:num.quantiles,1] + (1.96*summary(glm(PHEN ~ quantiles, family="binomial", data = for.quantiles))$coefficients[1:num.quantiles,2]))
          ci.quantiles.l <- exp(summary(glm(PHEN ~ quantiles, family="binomial", data = for.quantiles))$coefficients[1:num.quantiles,1] - (1.96*summary(glm(PHEN ~ quantiles, family="binomial", data = for.quantiles))$coefficients[1:num.quantiles,2]))
        }
        if(covary == T){
          for.quantiles <- merge(x = for.quantiles, by.x = "ID", y = pcs, by.y = "IID")

          or.quantiles <- exp(summary(glm(PHEN ~ ., family="binomial", data = for.quantiles[,c("PHEN","quantiles",covariates)]))$coefficients[1:num.quantiles,1])
          ci.quantiles.u <- exp(summary(glm(PHEN ~ ., family="binomial", data = for.quantiles[,c("PHEN","quantiles",covariates)]))$coefficients[1:num.quantiles,1] + (1.96*summary(glm(PHEN ~ ., family="binomial", data = for.quantiles[,c("PHEN","quantiles",covariates)]))$coefficients[1:num.quantiles,2]))
          ci.quantiles.l <- exp(summary(glm(PHEN ~ ., family="binomial", data = for.quantiles[,c("PHEN","quantiles",covariates)]))$coefficients[1:num.quantiles,1] - (1.96*summary(glm(PHEN ~ ., family="binomial", data = for.quantiles[,c("PHEN","quantiles",covariates)]))$coefficients[1:num.quantiles,2]))
        }
        or.quantiles[1] <- 1 
        ci.quantiles.u[1] <- 1
        ci.quantiles.l[1] <- 1
        quant.list <- seq(1, num.quantiles, 1)
        quant.list <- quant.list[quant.list != quant.ref]     
        quantiles.for.table <- c(quant.ref, quant.list)
        quantiles.df <- data.frame(or.quantiles, ci.quantiles.u, ci.quantiles.l, quantiles.for.table)
        names(quantiles.df) <- c("OR", "CI.U", "CI.L", "DEC")
        quantiles.df <- quantiles.df[order(quantiles.df$DEC),]
        quantiles.plot <- ggplot(quantiles.df) + 
          geom_point(aes(x = DEC, y = OR), colour = "royalblue2", size=4) + 
          theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(), panel.background = element_blank(), 
          axis.line = element_line(colour = "black",size=0.5)) + 
          geom_pointrange(aes(ymin = CI.L,ymax = CI.U, y = OR, x = DEC), colour = "royalblue2", size = 0.9) + 
          ylab("Odds Ratio for Score on Phenotype") + 
          xlab("quantiles for Polygenic Score") + 
          scale_x_continuous(breaks=seq(0, num.quantiles, 1))   
        }  
      if(!binary.target){
        if(report.individual.scores & !report.best.score.only){
          scores.internal <- prof.all.scores[,c("IID", paste("pT", output$thresh[which.min(output$p.out)],sep="_"))]
        }
        if(report.individual.scores & report.best.score.only){
          scores.internal <- read.table(paste(figname,target.phenotypes[k], "SCORES_AT_BEST-FIT-PRS.txt", sep = "_"), head = T)
        }
        names(scores.internal) <- c("ID", "SCORE")
        for.quantiles <- merge(x = phen.file.internal, by.x = "ID", y = scores.internal, by.y = "ID")
        for.quantiles$quantiles <- as.numeric(cut(for.quantiles$SCORE, breaks = quantile(for.quantiles$SCORE, probs = seq(0, 1, 1/num.quantiles)), include.lowest=T))
        for.quantiles$quantiles <- factor(for.quantiles$quantiles, levels = c(quant.ref, seq(1, num.quantiles, 1)[seq(1, num.quantiles, 1) != quant.ref]))

        if(covary == F){
          coef.quantiles <- summary(glm(PHEN ~ quantiles, family="gaussian", data = for.quantiles))$coefficients[1:num.quantiles,1]
          ci.quantiles.u <- summary(glm(PHEN ~ quantiles, family="gaussian", data = for.quantiles))$coefficients[1:num.quantiles,1] + (1.96*summary(glm(PHEN ~ quantiles, family="gaussian", data = for.quantiles))$coefficients[1:num.quantiles,2])
          ci.quantiles.l <- summary(glm(PHEN ~ quantiles, family="gaussian", data = for.quantiles))$coefficients[1:num.quantiles,1] - (1.96*summary(glm(PHEN ~ quantiles, family="gaussian", data = for.quantiles))$coefficients[1:num.quantiles,2])
        }
        if(covary == T){
          for.quantiles <- merge(x = for.quantiles, by.x = "ID", y = pcs, by.y = "IID")
          coef.quantiles <- summary(glm(PHEN ~ ., family="gaussian", data = for.quantiles[,c("PHEN","quantiles",covariates)]))$coefficients[1:num.quantiles,1]
          ci.quantiles.u <- summary(glm(PHEN ~ ., family="gaussian", data = for.quantiles[,c("PHEN","quantiles",covariates)]))$coefficients[1:num.quantiles,1] + (1.96*summary(glm(PHEN ~ ., family="gaussian", data = for.quantiles[,c("PHEN","quantiles",covariates)]))$coefficients[1:num.quantiles,2])
          ci.quantiles.l <- summary(glm(PHEN ~ ., family="gaussian", data = for.quantiles[,c("PHEN","quantiles",covariates)]))$coefficients[1:num.quantiles,1] - (1.96*summary(glm(PHEN ~ ., family="gaussian", data = for.quantiles[,c("PHEN","quantiles",covariates)]))$coefficients[1:num.quantiles,2])
        }
        coef.quantiles[1] <- 0 
        ci.quantiles.u[1] <- 0
        ci.quantiles.l[1] <- 0
        quant.list <- seq(1, num.quantiles, 1)
        quant.list <- quant.list[quant.list != quant.ref]     
        quantiles.for.table <- c(quant.ref, quant.list)
        quantiles.df <- data.frame(coef.quantiles, ci.quantiles.u, ci.quantiles.l, quantiles.for.table)
        names(quantiles.df) <- c("Coef", "CI.U", "CI.L", "DEC")
        quantiles.df <- quantiles.df[order(quantiles.df$DEC),]
        quantiles.plot <- ggplot(quantiles.df) + 
          geom_point(aes(x = DEC, y = Coef), colour = "royalblue2", size=4) + 
          theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(), panel.background = element_blank(), 
          axis.line = element_line(colour = "black",size=0.5)) + 
          geom_pointrange(aes(ymin = CI.L,ymax = CI.U, y = Coef, x = DEC), colour = "royalblue2", size = 0.9) + 
          ylab("Change in Phenotype given score in quantiles") + 
          xlab("Quantiles for Polygenic Score") + 
          scale_x_continuous(breaks=seq(0, num.quantiles, 1))   
        }  
      quantiles.plot
      ggsave(filename=paste(figname,target.phenotypes[k], "QUANTILES_PLOT.png", sep = "_"),width=7,height=7)  
    }	

	      
	      if(!fastscore){
	        cat(" ################################# \n # \n #   High Density Plots \n # \n ################################# \n")
	      }
	      if(best.thresh.on.bar){
	        #add max thresh
	        barchart.levels <- c(barchart.levels.old, output$thresh[which.min(output$p.out)])
	        barchart.levels <- barchart.levels[!duplicated(barchart.levels)]	
	        barchart.levels <- sort(barchart.levels, decreasing = F)
	      }
	      if(!fastscore & !ggfig){
	        if(!scatter.R2){
	        png(paste(figname,"_",target.phenotypes[k],"_HIGH-RES_PLOT_", Sys.Date(), ".png", sep = ""))
	          with(output, plot(	x = output$thresh,
	            y = -log10(output$p.out),
	            xlab = expression(italic(P)-value~threshold),
	            ylab = bquote(-log[10]~model~italic(P)-value),
	            pch = 19,
	            type = "b"
	          ))
	          with(output[output$thresh %in% barchart.levels.old,], points(thresh, -log10(p.out), col = "green", pch = 19, type= "b"))
	        }
	        if(scatter.R2){
	          png(paste(figname,"_",target.phenotypes[k],"_HIGH-RES_PLOT_", Sys.Date(), ".png", sep = ""))
	          with(output, plot(	x = output$thresh,
	            y = r2.out,
	            xlab = expression(italic(P)-value~threshold),
	            ylab = bquote(-log[10]~model~italic(P)-value),
	            pch = 19,
	            type = "b"
	          ))
	          with(output[output$thresh %in% barchart.levels.old,], points(thresh, r2.out, col = "green", pch = 19, type= "b"))
	        }	
	        dev.off()
	       }  
	      if(!fastscore & ggfig){
	        if(!scatter.R2){
	          ggfig.points <- ggplot(data = output) + geom_point(aes(x = thresh, y = -log10(p.out))) +
	            geom_line(aes(x = thresh, y = -log10(p.out))) +
                xlab(expression(italic(P)-value~threshold~(italic(P)[T])))  +
                ylab(bquote(PRS~model~fit:~italic(P)-value~(-log[10]))) +
	            geom_line(aes(thresh,  -log10(p.out)), colour = "green",
	            subset = .(thresh %in% barchart.levels.old))  +
	            geom_point(aes(thresh,  -log10(p.out)), colour = "green",
	            subset = .(thresh %in% barchart.levels.old))  + 
	          theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(), panel.background = element_blank(), 
          axis.line = element_line(colour = "black",size=0.5))
	          ggfig.points
	          ggsave(filename=paste(figname,"_",target.phenotypes[k],"_HIGH-RES_PLOT_", Sys.Date(), ".png", sep = ""),width=7,height=7)
	        }
	        if(scatter.R2){
	          ggfig.points <- ggplot(data = output) + geom_point(aes(x = thresh, y = r2.out)) +
	            geom_line(aes(x = thresh, y = r2.out)) +
                xlab(expression(italic(P)-value~threshold~(italic(P)[T])))  +
	            geom_line(aes(thresh,  r2.out), colour = "green",
	            subset = .(thresh %in% barchart.levels.old))  +
	            geom_point(aes(thresh,  r2.out), colour = "green",
	            subset = .(thresh %in% barchart.levels.old))  + 
	          theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(), panel.background = element_blank(), 
          axis.line = element_line(colour = "black",size=0.5))
   if(binary.target){
   	  if(!calculate.abc.r2){
        ggfig.points <- ggfig.points + ylab(expression(paste("PRS model fit:  ", R^2, " (Nagelkerke's)", sep = " ")))
      }
   	  if(calculate.abc.r2){
        ggfig.points <- ggfig.points + ylab(expression(paste("PRS model fit:  ", R^2, " (Liability Scale)", sep = " ")))
      }
    }
              if(!binary.target){
                ggfig.points <- ggfig.points + ylab(paste("PRS model fit:  ",expression(paste(R^2))))
              }
	          ggfig.points
	          ggsave(filename=paste(figname,"_",target.phenotypes[k],"_HIGH-RES_PLOT_", Sys.Date(), ".png", sep = ""),width=7,height=7)
	        }	 
	      }
	      options(scipen=0,digits=7)
	      #output <- read.table(paste(figname, "RAW_RESULTS_DATA.txt", sep = "_"),head=T)
	    output <- output[,c(1:3)]
	    names(output) <- c("thresh", "p.out", "r2.out")
	    output <- output[output$thresh %in% barchart.levels,]
	    output$print.p[round(output$p.out, digits = 3) != 0] <- round(output$p.out[round(output$p.out, digits = 3) != 0], digits = 3)
	    output$print.p[round(output$p.out, digits = 3) == 0] <- format(output$p.out[round(output$p.out, digits = 3) == 0], digits=2)
	    output$print.p <- sub("e", "*x*10^", output$print.p)
	      cat(" ################################# \n # \n #   Barplots \n # Bars for inclusion can be changed using the barchart.levels option \n ################################# \n")
	      if(ggfig){
	        if(!bar.col.is.pval){
	          ggfig.plot <-   ggplot(data = output)  + geom_bar(aes(x = factor(thresh), y = r2.out, fill = factor(thresh)), stat="identity") +     scale_fill_brewer(palette=barpalatte, name = expression(italic(P)-value~threshold)) + theme(axis.text.x=element_blank(), axis.title.x=element_blank())  +   xlab(expression(italic(P)-value~threshold~(italic(P)[T])))  
	        }
	        if(bar.col.is.pval){
	          ggfig.plot <-   ggplot(data = output)  + geom_bar(aes(x = factor(thresh), y = r2.out, fill = -log10(p.out)), stat="identity") +     scale_fill_gradient(low= bar.col.is.pval.lowcol, high= bar.col.is.pval.highcol, name =bquote(atop(-log[10]~model,italic(P)-value),)) +  xlab(expression(italic(P)-value~threshold~(italic(P)[T])))  
	        }
	        ggfig.plot <- ggfig.plot + geom_text(aes(x = factor(thresh), y = r2.out, label = paste(print.p)), vjust = -1.5, hjust = 0, angle = 45, cex = 2.8, parse=T) +
	          scale_y_continuous(limits = c(0, max(output$r2.out)*1.25))  + 
	          theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(), panel.background = element_blank(), 
          axis.line = element_line(colour = "black",size=0.5))
   if(binary.target){
   	  if(!calculate.abc.r2){
        ggfig.plot <- ggfig.plot + ylab(expression(paste("PRS model fit:  ", R^2, " (Nagelkerke's)", sep = " ")))
      }
   	  if(calculate.abc.r2){
        ggfig.plot <- ggfig.plot + ylab(expression(paste("PRS model fit:  ", R^2, " (Liability Scale)", sep = " ")))
      }
    }
    if(!binary.target){
      ggfig.plot <- ggfig.plot + ylab(expression(paste("PRS model fit:  ",R^2)))
    }
	          ggfig.plot
	          ggsave(filename=paste(figname,"_",target.phenotypes[k],"_BARPLOT_", Sys.Date(), ".png", sep = ""),width=7,height=7)		
	      }  
	      if(!ggfig) {
	        png(paste(figname,"_",target.phenotypes[k],"_BARPLOT_", Sys.Date(), ".png", sep = ""))
	        plot.fig  <- with(output, barplot(r2.out, 
	          names = thresh,
	          main = "", 
	          col = "red",
	          xlab = expression(italic(P)[T]),
	          ylab = expression(R^2),
	          ylim = c(0,  max(output$r2.out)*1.25)  ))
	          text( parse(text=paste(
	          output$print.p)), 
	          x = plot.fig+0.1, 
	          y =  output$r2.out+ (max(output$r2.out)*1.125-max(output$r2.out)), 
	          srt = 45)
	        dev.off()
	      }
	    }
	  }
	}
	if(multiple.base.phenotypes){
	  if(!multiple.target.phenotypes){
	    output <- read.table(paste(figname,base.phenotypes.names[basePhen], "RAW_RESULTS_DATA.txt", sep = "_"), head=T)
            if(quantiles & ggfig){
    cat(" ################################# \n # \n #   Quantiles Plots \n # \n ################################# \n")
      if(binary.target){
        if(report.individual.scores & !report.best.score.only){
          scores.internal <- prof.all.scores[,c("IID", paste("pT", output$thresh[which.min(output$p.out)],sep="_"))]
        }
        if(report.individual.scores & report.best.score.only){
          scores.internal <- read.table(paste(figname,base.phenotypes.names[basePhen], "SCORES_AT_BEST-FIT-PRS.txt", sep = "_"), head = T)
        }
        names(scores.internal) <- c("ID", "SCORE")
                for.quantiles <- merge(x = phen.file.internal, by.x = "ID", y = scores.internal, by.y = "ID")
        for.quantiles$quantiles <- as.numeric(cut(for.quantiles$SCORE, breaks = quantile(for.quantiles$SCORE, probs = seq(0, 1, 1/num.quantiles)), include.lowest=T))
        for.quantiles$quantiles <- factor(for.quantiles$quantiles, levels = c(quant.ref, seq(1, num.quantiles, 1)[seq(1, num.quantiles, 1) != quant.ref]))

        if(covary == F){
          or.quantiles <- exp(summary(glm(PHEN ~ quantiles, family="binomial", data = for.quantiles))$coefficients[1:num.quantiles,1])
          ci.quantiles.u <- exp(summary(glm(PHEN ~ quantiles, family="binomial", data = for.quantiles))$coefficients[1:num.quantiles,1] + (1.96*summary(glm(PHEN ~ quantiles, family="binomial", data = for.quantiles))$coefficients[1:num.quantiles,2]))
          ci.quantiles.l <- exp(summary(glm(PHEN ~ quantiles, family="binomial", data = for.quantiles))$coefficients[1:num.quantiles,1] - (1.96*summary(glm(PHEN ~ quantiles, family="binomial", data = for.quantiles))$coefficients[1:num.quantiles,2]))
        }
        if(covary == T){
          for.quantiles <- merge(x = for.quantiles, by.x = "ID", y = pcs, by.y = "IID")

          or.quantiles <- exp(summary(glm(PHEN ~ ., family="binomial", data = for.quantiles[,c("PHEN","quantiles",covariates)]))$coefficients[1:num.quantiles,1])
          ci.quantiles.u <- exp(summary(glm(PHEN ~ ., family="binomial", data = for.quantiles[,c("PHEN","quantiles",covariates)]))$coefficients[1:num.quantiles,1] + (1.96*summary(glm(PHEN ~ ., family="binomial", data = for.quantiles[,c("PHEN","quantiles",covariates)]))$coefficients[1:num.quantiles,2]))
          ci.quantiles.l <- exp(summary(glm(PHEN ~ ., family="binomial", data = for.quantiles[,c("PHEN","quantiles",covariates)]))$coefficients[1:num.quantiles,1] - (1.96*summary(glm(PHEN ~ ., family="binomial", data = for.quantiles[,c("PHEN","quantiles",covariates)]))$coefficients[1:num.quantiles,2]))
        }
        or.quantiles[1] <- 1 
        ci.quantiles.u[1] <- 1
        ci.quantiles.l[1] <- 1
        quant.list <- seq(1, num.quantiles, 1)
        quant.list <- quant.list[quant.list != quant.ref]     
        quantiles.for.table <- c(quant.ref, quant.list)
        quantiles.df <- data.frame(or.quantiles, ci.quantiles.u, ci.quantiles.l, quantiles.for.table)
        names(quantiles.df) <- c("OR", "CI.U", "CI.L", "DEC")
        quantiles.df <- quantiles.df[order(quantiles.df$DEC),]
        quantiles.plot <- ggplot(quantiles.df) + 
          geom_point(aes(x = DEC, y = OR), colour = "royalblue2", size=4) + 
          theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(), panel.background = element_blank(), 
          axis.line = element_line(colour = "black",size=0.5)) + 
          geom_pointrange(aes(ymin = CI.L,ymax = CI.U, y = OR, x = DEC), colour = "royalblue2", size = 0.9) + 
          ylab("Odds Ratio for Score on Phenotype") + 
          xlab("quantiles for Polygenic Score") + 
          scale_x_continuous(breaks=seq(0, num.quantiles, 1))   
        }  
      if(!binary.target){
        if(report.individual.scores & !report.best.score.only){
          scores.internal <- prof.all.scores[,c("IID", paste("pT", output$thresh[which.min(output$p.out)],sep="_"))]
        }
        if(report.individual.scores & report.best.score.only){
          scores.internal <- read.table(paste(figname, base.phenotypes.names[basePhen],"SCORES_AT_BEST-FIT-PRS.txt", sep = "_"), head = T)
        }
        names(scores.internal) <- c("ID", "SCORE")
        for.quantiles <- merge(x = phen.file.internal, by.x = "ID", y = scores.internal, by.y = "ID")
        for.quantiles$quantiles <- as.numeric(cut(for.quantiles$SCORE, breaks = quantile(for.quantiles$SCORE, probs = seq(0, 1, 1/num.quantiles)), include.lowest=T))
        for.quantiles$quantiles <- factor(for.quantiles$quantiles, levels = c(quant.ref, seq(1, num.quantiles, 1)[seq(1, num.quantiles, 1) != quant.ref]))

        if(covary == F){
          coef.quantiles <- summary(glm(PHEN ~ quantiles, family="gaussian", data = for.quantiles))$coefficients[1:num.quantiles,1]
          ci.quantiles.u <- summary(glm(PHEN ~ quantiles, family="gaussian", data = for.quantiles))$coefficients[1:num.quantiles,1] + (1.96*summary(glm(PHEN ~ quantiles, family="gaussian", data = for.quantiles))$coefficients[1:num.quantiles,2])
          ci.quantiles.l <- summary(glm(PHEN ~ quantiles, family="gaussian", data = for.quantiles))$coefficients[1:num.quantiles,1] - (1.96*summary(glm(PHEN ~ quantiles, family="gaussian", data = for.quantiles))$coefficients[1:num.quantiles,2])
        }
        if(covary == T){
          for.quantiles <- merge(x = for.quantiles, by.x = "ID", y = pcs, by.y = "IID")
          coef.quantiles <- summary(glm(PHEN ~ ., family="gaussian", data = for.quantiles[,c("PHEN","quantiles",covariates)]))$coefficients[1:num.quantiles,1]
          ci.quantiles.u <- summary(glm(PHEN ~ ., family="gaussian", data = for.quantiles[,c("PHEN","quantiles",covariates)]))$coefficients[1:num.quantiles,1] + (1.96*summary(glm(PHEN ~ ., family="gaussian", data = for.quantiles[,c("PHEN","quantiles",covariates)]))$coefficients[1:num.quantiles,2])
          ci.quantiles.l <- summary(glm(PHEN ~ ., family="gaussian", data = for.quantiles[,c("PHEN","quantiles",covariates)]))$coefficients[1:num.quantiles,1] - (1.96*summary(glm(PHEN ~ ., family="gaussian", data = for.quantiles[,c("PHEN","quantiles",covariates)]))$coefficients[1:num.quantiles,2])
        }
        coef.quantiles[1] <- 0 
        ci.quantiles.u[1] <- 0
        ci.quantiles.l[1] <- 0
        quant.list <- seq(1, num.quantiles, 1)
        quant.list <- quant.list[quant.list != quant.ref]     
        quantiles.for.table <- c(quant.ref, quant.list)
        quantiles.df <- data.frame(coef.quantiles, ci.quantiles.u, ci.quantiles.l, quantiles.for.table)
        names(quantiles.df) <- c("Coef", "CI.U", "CI.L", "DEC")
        quantiles.df <- quantiles.df[order(quantiles.df$DEC),]
        quantiles.plot <- ggplot(quantiles.df) + 
          geom_point(aes(x = DEC, y = Coef), colour = "royalblue2", size=4) + 
          theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(), panel.background = element_blank(), 
          axis.line = element_line(colour = "black",size=0.5)) + 
          geom_pointrange(aes(ymin = CI.L,ymax = CI.U, y = Coef, x = DEC), colour = "royalblue2", size = 0.9) + 
          ylab("Change in Phenotype given score in quantiles") + 
          xlab("quantiles for Polygenic Score") + 
          scale_x_continuous(breaks=seq(0, num.quantiles, 1))   
        }  
      quantiles.plot
      ggsave(filename=paste(figname, base.phenotypes.names[basePhen],"QUANTILES_PLOT.png", sep = "_"),width=7,height=7)  
    }	

	    if(!fastscore){
	      cat(" ################################# \n # \n #   High Density Plots \n # \n ################################# \n")
	    }
#	    barchart.levels.old <- barchart.levels
	    if(best.thresh.on.bar){
	      #add max thresh
	      barchart.levels <- c(barchart.levels.old, output$thresh[which.min(output$p.out)])
	      barchart.levels <- barchart.levels[!duplicated(barchart.levels)]	
	      barchart.levels <- sort(barchart.levels, decreasing = F)
	    }
	    if(!fastscore & !ggfig){
	      if(!scatter.R2){
	        png(paste(figname,"_",base.phenotypes.names[basePhen],"_HIGH-RES_PLOT_", Sys.Date(), ".png", sep = ""))
	          with(output, plot(	x = output$thresh,
	          y = -log10(output$p.out),
	          xlab = expression(italic(P)-value~threshold),
	          ylab = bquote(-log[10]~model~italic(P)-value),
	          pch = 19,
	          type = "b"
	        ))
	        with(output[output$thresh %in% barchart.levels.old,], points(thresh, -log10(p.out), col = "green", pch = 19, type= "b"))
	      }
	      if(scatter.R2){
	        png(paste(figname,"_",base.phenotypes.names[basePhen],"_HIGH-RES_PLOT_", Sys.Date(), ".png", sep = ""))
	        with(output, plot(	x = output$thresh,
	          y = r2.out,
	          xlab = expression(italic(P)-value~threshold),
	          ylab = bquote(-log[10]~model~italic(P)-value),
	          pch = 19,
	          type = "b"
	        ))
	        with(output[output$thresh %in% barchart.levels.old,], points(thresh, r2.out, col = "green", pch = 19, type= "b"))
	      }	
	      dev.off()
	    }
	    if(!fastscore & ggfig){
	      if(!scatter.R2){
	        ggfig.points <- ggplot(data = output) + geom_point(aes(x = thresh, y = -log10(p.out))) +
	          geom_line(aes(x = thresh, y = -log10(p.out))) +
              xlab(expression(italic(P)-value~threshold~(italic(P)[T])))  +
              ylab(bquote(PRS~model~fit:~italic(P)-value~(-log[10]))) +
	          geom_line(aes(thresh,  -log10(p.out)), colour = "green",
	          ,subset = .(thresh %in% barchart.levels.old))  +
	          geom_point(aes(thresh,  -log10(p.out)), colour = "green",
	          ,subset = .(thresh %in% barchart.levels.old))  + 
	          theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(), panel.background = element_blank(), 
          axis.line = element_line(colour = "black",size=0.5))
	        ggfig.points
	        ggsave(filename=paste(figname,"_",base.phenotypes.names[basePhen],"_HIGH-RES_PLOT_", Sys.Date(), ".png", sep = ""),width=7,height=7)
	      }
	      if(scatter.R2){
	        ggfig.points <- ggplot(data = output) + geom_point(aes(x = thresh, y = r2.out)) +
	          geom_line(aes(x = thresh, y = r2.out)) +
              xlab(expression(italic(P)-value~threshold~(italic(P)[T])))  +
	          geom_line(aes(thresh,  r2.out), colour = "green",
	          ,subset = .(thresh %in% barchart.levels.old))  +
	          geom_point(aes(thresh,  r2.out), colour = "green",
	          ,subset = .(thresh %in% barchart.levels.old))  + 
	          theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(), panel.background = element_blank(), 
          axis.line = element_line(colour = "black",size=0.5))
   if(binary.target){
   	  if(!calculate.abc.r2){
        ggfig.points <- ggfig.points + ylab(expression(paste("PRS model fit:  ", R^2, " (Nagelkerke's)", sep = " ")))
      }
   	  if(calculate.abc.r2){
        ggfig.points <- ggfig.points + ylab(expression(paste("PRS model fit:  ", R^2, " (Liability Scale)", sep = " ")))
      }
    }
            if(!binary.target){
              ggfig.points <- ggfig.points + ylab(paste("PRS model fit:  ",expression(paste(R^2))))
            }
	        ggfig.points
	        ggsave(filename=paste(figname,"_",base.phenotypes.names[basePhen],"_HIGH-RES_PLOT_", Sys.Date(), ".png", sep = ""),width=7,height=7)
	      }	 
	    }
	    options(scipen=0,digits=7)
	    output <- read.table(paste(figname,base.phenotypes.names[basePhen], "RAW_RESULTS_DATA.txt", sep = "_"),head=T)
	    output <- output[,c(1:3)]
	    names(output) <- c("thresh", "p.out", "r2.out")
	    output <- output[output$thresh %in% barchart.levels,]
	    output$print.p[round(output$p.out, digits = 3) != 0] <- round(output$p.out[round(output$p.out, digits = 3) != 0], digits = 3)
	    output$print.p[round(output$p.out, digits = 3) == 0] <- format(output$p.out[round(output$p.out, digits = 3) == 0], digits=2)
	    output$print.p <- sub("e", "*x*10^", output$print.p)
	    cat(" ################################# \n # \n #   Barplots \n # Bars for inclusion can be changed using the barchart.levels option \n ################################# \n")
	    if(ggfig){
	      if(!bar.col.is.pval){
	        ggfig.plot <-   ggplot(data = output)  + geom_bar(aes(x = factor(thresh), y = r2.out, fill = factor(thresh)), stat="identity") +     scale_fill_brewer(palette=barpalatte, name = expression(italic(P)-value~threshold)) + theme(axis.text.x=element_blank(), axis.title.x=element_blank())  +   xlab(expression(italic(P)-value~threshold~(italic(P)[T])))  
	      }
	      if(bar.col.is.pval){
	        ggfig.plot <-   ggplot(data = output)  + geom_bar(aes(x = factor(thresh), y = r2.out, fill = -log10(p.out)), stat="identity") +     scale_fill_gradient(low= bar.col.is.pval.lowcol, high= bar.col.is.pval.highcol, name =bquote(atop(-log[10]~model,italic(P)-value),)) +  xlab(expression(italic(P)-value~threshold~(italic(P)[T])))  	      }
	      ggfig.plot <- ggfig.plot + geom_text(aes(x = factor(thresh), y = r2.out, label = paste(print.p)), vjust = -1.5, hjust = 0, angle = 45, cex = 2.8, parse=T) +
	        scale_y_continuous(limits = c(0, max(output$r2.out)*1.25))  + 
	          theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(), panel.background = element_blank(), 
          axis.line = element_line(colour = "black",size=0.5))
   if(binary.target){
   	  if(!calculate.abc.r2){
        ggfig.plot <- ggfig.plot + ylab(expression(paste("PRS model fit:  ", R^2, " (Nagelkerke's)", sep = " ")))
      }
   	  if(calculate.abc.r2){
        ggfig.plot <- ggfig.plot + ylab(expression(paste("PRS model fit:  ", R^2, " (Liability Scale)", sep = " ")))
      }
    }
            if(!binary.target){
              ggfig.plot <- ggfig.plot + ylab(expression(paste("PRS model fit:  ",R^2)))
            }
	        ggfig.plot
	        ggsave(filename=paste(figname,"_",base.phenotypes.names[basePhen],"_BARPLOT_", Sys.Date(), ".png", sep = ""),width=7,height=7)		
	    }  
	    if(!ggfig) {
	      png(paste(figname,"_",base.phenotypes.names[basePhen],"_BARPLOT_", Sys.Date(), ".png", sep = ""))
	      plot.fig  <- with(output, barplot(r2.out, 
	        names = thresh,
	        main = "", 
	        col = "red",
	        xlab = expression(italic(P)[T]),
	        ylab = expression(R^2),
	        ylim = c(0,  max(output$r2.out)*1.25)  ))
	        text( parse(text=paste(
	        output$print.p)), 
	        x = plot.fig+0.1, 
	        y =  output$r2.out+ (max(output$r2.out)*1.125-max(output$r2.out)), 
	        srt = 45)
	      dev.off()
	    }
	  }
	  if(multiple.target.phenotypes){
	    for(k in 1:length(target.phenotypes)){
	      output <- read.table(paste(figname, base.phenotypes.names[basePhen], "PREDICTING", target.phenotypes[k], "RAW_RESULTS_DATA.txt", sep = "_"), head=T)
            if(quantiles & ggfig){
    cat(" ################################# \n # \n #   Quantiles Plots \n # \n ################################# \n")
      if(binary.target){
        if(report.individual.scores & !report.best.score.only){
          scores.internal <- prof.all.scores[,c("IID", paste("pT", output$thresh[which.min(output$p.out)],sep="_"))]
        }
        if(report.individual.scores & report.best.score.only){
          scores.internal <- read.table(paste(figname,base.phenotypes.names[basePhen], "PREDICTING", target.phenotypes[k], "SCORES_AT_BEST-FIT-PRS.txt", sep = "_"), head = T)
        }
        names(scores.internal) <- c("ID", "SCORE")
        for.quantiles <- merge(x = phen.file.internal, by.x = "ID", y = scores.internal, by.y = "ID")
        for.quantiles$quantiles <- as.numeric(cut(for.quantiles$SCORE, breaks = quantile(for.quantiles$SCORE, probs = seq(0, 1, 1/num.quantiles)), include.lowest=T))
        for.quantiles$quantiles <- factor(for.quantiles$quantiles, levels = c(quant.ref, seq(1, num.quantiles, 1)[seq(1, num.quantiles, 1) != quant.ref]))

        if(covary == F){
          or.quantiles <- exp(summary(glm(PHEN ~ quantiles, family="binomial", data = for.quantiles))$coefficients[1:num.quantiles,1])
          ci.quantiles.u <- exp(summary(glm(PHEN ~ quantiles, family="binomial", data = for.quantiles))$coefficients[1:num.quantiles,1] + (1.96*summary(glm(PHEN ~ quantiles, family="binomial", data = for.quantiles))$coefficients[1:num.quantiles,2]))
          ci.quantiles.l <- exp(summary(glm(PHEN ~ quantiles, family="binomial", data = for.quantiles))$coefficients[1:num.quantiles,1] - (1.96*summary(glm(PHEN ~ quantiles, family="binomial", data = for.quantiles))$coefficients[1:num.quantiles,2]))
        }
        if(covary == T){
          for.quantiles <- merge(x = for.quantiles, by.x = "ID", y = pcs, by.y = "IID")

          or.quantiles <- exp(summary(glm(PHEN ~ ., family="binomial", data = for.quantiles[,c("PHEN","quantiles",covariates)]))$coefficients[1:num.quantiles,1])
          ci.quantiles.u <- exp(summary(glm(PHEN ~ ., family="binomial", data = for.quantiles[,c("PHEN","quantiles",covariates)]))$coefficients[1:num.quantiles,1] + (1.96*summary(glm(PHEN ~ ., family="binomial", data = for.quantiles[,c("PHEN","quantiles",covariates)]))$coefficients[1:num.quantiles,2]))
          ci.quantiles.l <- exp(summary(glm(PHEN ~ ., family="binomial", data = for.quantiles[,c("PHEN","quantiles",covariates)]))$coefficients[1:num.quantiles,1] - (1.96*summary(glm(PHEN ~ ., family="binomial", data = for.quantiles[,c("PHEN","quantiles",covariates)]))$coefficients[1:num.quantiles,2]))
        }
        or.quantiles[1] <- 1 
        ci.quantiles.u[1] <- 1
        ci.quantiles.l[1] <- 1
        quant.list <- seq(1, num.quantiles, 1)
        quant.list <- quant.list[quant.list != quant.ref]     
        quantiles.for.table <- c(quant.ref, quant.list)
        quantiles.df <- data.frame(or.quantiles, ci.quantiles.u, ci.quantiles.l, quantiles.for.table)
        names(quantiles.df) <- c("OR", "CI.U", "CI.L", "DEC")
        quantiles.df <- quantiles.df[order(quantiles.df$DEC),]
        quantiles.plot <- ggplot(quantiles.df) + 
          geom_point(aes(x = DEC, y = OR), colour = "royalblue2", size=4) + 
          theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(), panel.background = element_blank(), 
          axis.line = element_line(colour = "black",size=0.5)) + 
          geom_pointrange(aes(ymin = CI.L,ymax = CI.U, y = OR, x = DEC), colour = "royalblue2", size = 0.9) + 
          ylab("Odds Ratio for Score on Phenotype") + 
          xlab("quantiles for Polygenic Score") + 
          scale_x_continuous(breaks=seq(0, num.quantiles, 1))   
        }  
      if(!binary.target){
        if(report.individual.scores & !report.best.score.only){
          scores.internal <- prof.all.scores[,c("IID", paste("pT", output$thresh[which.min(output$p.out)],sep="_"))]
        }
        if(report.individual.scores & report.best.score.only){
          scores.internal <- read.table(paste(figname, base.phenotypes.names[basePhen], "PREDICTING", target.phenotypes[k], "SCORES_AT_BEST-FIT-PRS.txt", sep = "_"), head = T)
        }
        names(scores.internal) <- c("ID", "SCORE")
        for.quantiles <- merge(x = phen.file.internal, by.x = "ID", y = scores.internal, by.y = "ID")
        for.quantiles$quantiles <- as.numeric(cut(for.quantiles$SCORE, breaks = quantile(for.quantiles$SCORE, probs = seq(0, 1, 1/num.quantiles)), include.lowest=T))
        for.quantiles$quantiles <- factor(for.quantiles$quantiles, levels = c(quant.ref, seq(1, num.quantiles, 1)[seq(1, num.quantiles, 1) != quant.ref]))

        if(covary == F){
          coef.quantiles <- summary(glm(PHEN ~ quantiles, family="gaussian", data = for.quantiles))$coefficients[1:num.quantiles,1]
          ci.quantiles.u <- summary(glm(PHEN ~ quantiles, family="gaussian", data = for.quantiles))$coefficients[1:num.quantiles,1] + (1.96*summary(glm(PHEN ~ quantiles, family="gaussian", data = for.quantiles))$coefficients[1:num.quantiles,2])
          ci.quantiles.l <- summary(glm(PHEN ~ quantiles, family="gaussian", data = for.quantiles))$coefficients[1:num.quantiles,1] - (1.96*summary(glm(PHEN ~ quantiles, family="gaussian", data = for.quantiles))$coefficients[1:num.quantiles,2])
        }
        if(covary == T){
          for.quantiles <- merge(x = for.quantiles, by.x = "ID", y = pcs, by.y = "IID")
          coef.quantiles <- summary(glm(PHEN ~ ., family="gaussian", data = for.quantiles[,c("PHEN","quantiles",covariates)]))$coefficients[1:num.quantiles,1]
          ci.quantiles.u <- summary(glm(PHEN ~ ., family="gaussian", data = for.quantiles[,c("PHEN","quantiles",covariates)]))$coefficients[1:num.quantiles,1] + (1.96*summary(glm(PHEN ~ ., family="gaussian", data = for.quantiles[,c("PHEN","quantiles",covariates)]))$coefficients[1:num.quantiles,2])
          ci.quantiles.l <- summary(glm(PHEN ~ ., family="gaussian", data = for.quantiles[,c("PHEN","quantiles",covariates)]))$coefficients[1:num.quantiles,1] - (1.96*summary(glm(PHEN ~ ., family="gaussian", data = for.quantiles[,c("PHEN","quantiles",covariates)]))$coefficients[1:num.quantiles,2])
        }
        coef.quantiles[1] <- 0 
        ci.quantiles.u[1] <- 0
        ci.quantiles.l[1] <- 0
        quant.list <- seq(1, num.quantiles, 1)
        quant.list <- quant.list[quant.list != quant.ref]     
        quantiles.for.table <- c(quant.ref, quant.list)
        quantiles.df <- data.frame(coef.quantiles, ci.quantiles.u, ci.quantiles.l, quantiles.for.table)
        names(quantiles.df) <- c("Coef", "CI.U", "CI.L", "DEC")
        quantiles.df <- quantiles.df[order(quantiles.df$DEC),]
        quantiles.plot <- ggplot(quantiles.df) + 
          geom_point(aes(x = DEC, y = Coef), colour = "royalblue2", size=4) + 
          theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(), panel.background = element_blank(), 
          axis.line = element_line(colour = "black",size=0.5)) + 
          geom_pointrange(aes(ymin = CI.L,ymax = CI.U, y = Coef, x = DEC), colour = "royalblue2", size = 0.9) + 
          ylab("Change in Phenotype given score in quantiles") + 
          xlab("quantiles for Polygenic Score") + 
          scale_x_continuous(breaks=seq(0, num.quantiles, 1))   
        }  
      quantiles.plot
      ggsave(filename=paste(figname,base.phenotypes.names[basePhen], "PREDICTING", target.phenotypes[k], "QUANTILES_PLOT.png", sep = "_"),,width=7,height=7)  
    }	

	      if(!fastscore){
	        cat(" ################################# \n # \n #   High Density Plots \n # \n ################################# \n")
	      }
#	      barchart.levels.old <- barchart.levels
	      if(best.thresh.on.bar){
	        #add max thresh
	        barchart.levels <- c(barchart.levels.old, output$thresh[which.min(output$p.out)])
	        barchart.levels <- barchart.levels[!duplicated(barchart.levels)]	
	        barchart.levels <- sort(barchart.levels, decreasing = F)
	      }
	      if(!fastscore & !ggfig){
	        if(!scatter.R2){
	        png(paste(figname,"_",base.phenotypes.names[basePhen], "PREDICTING", target.phenotypes[k],"_HIGH-RES_PLOT_", Sys.Date(), ".png", sep = ""))
	          with(output, plot(	x = output$thresh,
	            y = -log10(output$p.out),
	            xlab = expression(italic(P)-value~threshold),
	            ylab = bquote(-log[10]~model~italic(P)-value),
	            pch = 19,
	            type = "b"
	          ))
	          with(output[output$thresh %in% barchart.levels.old,], points(thresh, -log10(p.out), col = "green", pch = 19, type= "b"))
	        }
	        if(scatter.R2){
	          png(paste(figname,"_",base.phenotypes.names[basePhen], "PREDICTING", target.phenotypes[k],"_HIGH-RES_PLOT_", Sys.Date(), ".png", sep = ""))
	          with(output, plot(	x = output$thresh,
	            y = r2.out,
	            xlab = expression(italic(P)-value~threshold),
	            ylab = bquote(Variance~Explained:~Pseudo~R^2),
	            pch = 19,
	            type = "b"
	          ))
	          with(output[output$thresh %in% barchart.levels.old,], points(thresh, r2.out, col = "green", pch = 19, type= "b"))
	        }	
	        dev.off()
	       }  
	      if(!fastscore & ggfig){
	        if(!scatter.R2){
	          ggfig.points <- ggplot(data = output) + geom_point(aes(x = thresh, y = -log10(p.out))) +
	            geom_line(aes(x = thresh, y = -log10(p.out))) +
                xlab(expression(italic(P)-value~threshold~(italic(P)[T])))  +
                ylab(bquote(PRS~model~fit:~italic(P)-value~(-log[10]))) +
	            geom_line(aes(thresh,  -log10(p.out)), colour = "green",
	            subset = .(thresh %in% barchart.levels.old))  +
	            geom_point(aes(thresh,  -log10(p.out)), colour = "green",
	            subset = .(thresh %in% barchart.levels.old))  + 
	          theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(), panel.background = element_blank(), 
          axis.line = element_line(colour = "black",size=0.5))
	          ggfig.points
	          ggsave(filename=paste(figname,"_",base.phenotypes.names[basePhen], "PREDICTING", target.phenotypes[k],"_HIGH-RES_PLOT_", Sys.Date(), ".png", sep = ""),width=7,height=7)
	        }
	        if(scatter.R2){
	          ggfig.points <- ggplot(data = output) + geom_point(aes(x = thresh, y = r2.out)) +
	            geom_line(aes(x = thresh, y = r2.out)) +
                xlab(expression(italic(P)-value~threshold~(italic(P)[T])))  +
	            geom_line(aes(thresh,  r2.out), colour = "green",
	            subset = .(thresh %in% barchart.levels.old))  +
	            geom_point(aes(thresh,  r2.out), colour = "green",
	            subset = .(thresh %in% barchart.levels.old))  + 
	          theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(), panel.background = element_blank(), 
          axis.line = element_line(colour = "black",size=0.5))
   if(binary.target){
   	  if(!calculate.abc.r2){
        ggfig.points <- ggfig.points + ylab(expression(paste("PRS model fit:  ", R^2, " (Nagelkerke's)", sep = " ")))
      }
   	  if(calculate.abc.r2){
        ggfig.points <- ggfig.points + ylab(expression(paste("PRS model fit:  ", R^2, " (Liability Scale)", sep = " ")))
      }
    }
              if(!binary.target){
                ggfig.points <- ggfig.points + ylab(paste("PRS model fit:  ",expression(paste(R^2))))
              }
	          ggfig.points
	          ggsave(filename=paste(figname,"_",base.phenotypes.names[basePhen], "PREDICTING", target.phenotypes[k],"_HIGH-RES_PLOT_", Sys.Date(), ".png", sep = ""),width=7,height=7)
	        }	 
	      }
	      options(scipen=0,digits=7)
	    output <- output[,c(1:3)]
	    names(output) <- c("thresh", "p.out", "r2.out")
	      output <- output[output$thresh %in% barchart.levels,]
	      output$print.p[round(output$p.out, digits = 3) != 0] <- round(output$p.out[round(output$p.out, digits = 3) != 0], digits = 3)
	      output$print.p[round(output$p.out, digits = 3) == 0] <- format(output$p.out[round(output$p.out, digits = 3) == 0], digits=2)
	      output$print.p <- sub("e", "*x*10^", output$print.p)
	      cat(" ################################# \n # \n #   Barplots \n # Bars for inclusion can be changed using the barchart.levels option \n ################################# \n")
	      if(ggfig){
	        if(!bar.col.is.pval){
	          ggfig.plot <-   ggplot(data = output)  + geom_bar(aes(x = factor(thresh), y = r2.out, fill = factor(thresh)), stat="identity") +     scale_fill_brewer(palette=barpalatte, name = expression(italic(P)-value~threshold)) + theme(axis.text.x=element_blank(), axis.title.x=element_blank())  +   xlab(expression(italic(P)-value~threshold~(italic(P)[T])))  
	        }
	        if(bar.col.is.pval){
	          ggfig.plot <-   ggplot(data = output)  + geom_bar(aes(x = factor(thresh), y = r2.out, fill = -log10(p.out)), stat="identity") +     scale_fill_gradient(low= bar.col.is.pval.lowcol, high= bar.col.is.pval.highcol, name =bquote(atop(-log[10]~model,italic(P)-value),)) +   xlab(expression(italic(P)-value~threshold~(italic(P)[T])))  
	        }
	        ggfig.plot <- ggfig.plot + geom_text(aes(x = factor(thresh), y = r2.out, label = paste(print.p)), vjust = -1.5, hjust = 0, angle = 45, cex = 2.8, parse=T) +
	          scale_y_continuous(limits = c(0, max(output$r2.out)*1.25))  + 
	          theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(), panel.background = element_blank(), 
          axis.line = element_line(colour = "black",size=0.5))
   if(binary.target){
   	  if(!calculate.abc.r2){
        ggfig.plot <- ggfig.plot + ylab(expression(paste("PRS model fit:  ", R^2, " (Nagelkerke's)", sep = " ")))
      }
   	  if(calculate.abc.r2){
        ggfig.plot <- ggfig.plot + ylab(expression(paste("PRS model fit:  ", R^2, " (Liability Scale)", sep = " ")))
      }
    }
              if(!binary.target){
                ggfig.plot <- ggfig.plot + ylab(expression(paste("PRS model fit:  ",R^2)))
              }
	          ggfig.plot
	          ggsave(filename=paste(figname,"_",base.phenotypes.names[basePhen], "PREDICTING", target.phenotypes[k],"_BARPLOT_", Sys.Date(), ".png", sep = ""),width=7,height=7)		
	      }  
	      if(!ggfig) {
	        png(paste(figname,"_",base.phenotypes.names[basePhen], "PREDICTING", target.phenotypes[k],"_BARPLOT_", Sys.Date(), ".png", sep = ""))
	        plot.fig  <- with(output, barplot(r2.out, 
	          names = thresh,
	          main = "", 
	          col = "red",
	          xlab = expression(italic(P)[T]),
	          ylab = expression(R^2),
	          ylim = c(0,  max(output$r2.out)*1.25)  ))
	          text( parse(text=paste(
	          output$print.p)), 
	          x = plot.fig+0.1, 
	          y =  output$r2.out+ (max(output$r2.out)*1.125-max(output$r2.out)), 
	          srt = 45)
	        dev.off()
	      }
	    }
	  }
	}
  if(multiple.base.phenotypes & (basePhen == 1) & !heat.r2){
    out.multibase <- data.frame(top.thresh.pval)
  }
  if(multiple.base.phenotypes & (basePhen > 1) & !heat.r2){
    out.multibase <- data.frame(out.multibase, top.thresh.pval)
  }
  if(multiple.base.phenotypes & (basePhen == 1) & heat.r2){
    out.multibase <- data.frame(top.thresh.r2)
  }
  if(multiple.base.phenotypes & (basePhen > 1) & heat.r2){
    out.multibase <- data.frame(out.multibase, top.thresh.r2)
  }
  if(multiple.base.phenotypes){
	    cat(" ################################# \n # \n #   Remove temp files \n # \n ################################# \n")
	    system("rm PROFILES.0.5*.profile", ignore.stdout=T,ignore.stderr=T)	
	    system("rm PROFILES.0.4*.profile", ignore.stdout=T,ignore.stderr=T)
	    system("rm PROFILES.0.3*.profile", ignore.stdout=T,ignore.stderr=T)
	    system("rm PROFILES.0.2*.profile", ignore.stdout=T,ignore.stderr=T)
	    system("rm PROFILES.0.1*.profile", ignore.stdout=T,ignore.stderr=T)
	    system("rm PROFILES.0.0*.profile", ignore.stdout=T,ignore.stderr=T)
	    system("rm PROFILES.*.profile", ignore.stdout=T,ignore.stderr=T)
	    system("rm PROFILES.S1*.profile", ignore.stdout=T,ignore.stderr=T)	
	    system("rm PROFILES.S2*.profile", ignore.stdout=T,ignore.stderr=T)
	    system("rm PROFILES.S3*.profile", ignore.stdout=T,ignore.stderr=T)
	    system("rm PROFILES.S4*.profile", ignore.stdout=T,ignore.stderr=T)
	    system("rm PROFILES.S5*.profile", ignore.stdout=T,ignore.stderr=T)
	    system("rm PROFILES.S6*.profile", ignore.stdout=T,ignore.stderr=T)
	    system("rm PROFILES.S7*.profile", ignore.stdout=T,ignore.stderr=T)
	    system("rm PROFILES.S8*.profile", ignore.stdout=T,ignore.stderr=T)
	    system("rm PROFILES.S9*.profile", ignore.stdout=T,ignore.stderr=T)
	    system("rm *PROFILES.S*.profile", ignore.stdout=T,ignore.stderr=T)
	    system("rm *PROFILES.0.5*.profile", ignore.stdout=T,ignore.stderr=T)	
	    system("rm *PROFILES.0.4*.profile", ignore.stdout=T,ignore.stderr=T)
	    system("rm *PROFILES.0.3*.profile", ignore.stdout=T,ignore.stderr=T)
	    system("rm *PROFILES.0.2*.profile", ignore.stdout=T,ignore.stderr=T)
	    system("rm *PROFILES.0.1*.profile", ignore.stdout=T,ignore.stderr=T)
	    system("rm *PROFILES.0.0*.profile", ignore.stdout=T,ignore.stderr=T)
	    system("rm *PROFILES.*.profile", ignore.stdout=T,ignore.stderr=T)
	    system("rm *PROFILES.S1*.profile", ignore.stdout=T,ignore.stderr=T)	
	    system("rm *PROFILES.S2*.profile", ignore.stdout=T,ignore.stderr=T)
	    system("rm *PROFILES.S3*.profile", ignore.stdout=T,ignore.stderr=T)
	    system("rm *PROFILES.S4*.profile", ignore.stdout=T,ignore.stderr=T)
	    system("rm *PROFILES.S5*.profile", ignore.stdout=T,ignore.stderr=T)
	    system("rm *PROFILES.S6*.profile", ignore.stdout=T,ignore.stderr=T) 
	    system("rm *PROFILES.S7*.profile", ignore.stdout=T,ignore.stderr=T)
	    system("rm *PROFILES.S8*.profile", ignore.stdout=T,ignore.stderr=T)
	    system("rm *PROFILES.S9*.profile", ignore.stdout=T,ignore.stderr=T)
	    system("rm *PROFILES.S*.profile", ignore.stdout=T,ignore.stderr=T) 
	    system("rm flip*", ignore.stdout=T,ignore.stderr=T)
	    system("rm clean*", ignore.stdout=T,ignore.stderr=T)
	    system("rm profile_list", ignore.stdout=T,ignore.stderr=T)
	    system("rm rangelist.txt", ignore.stdout=T,ignore.stderr=T)
	    system("rm rangelist_ranges", ignore.stdout=T,ignore.stderr=T)
	    system("rm rawfile.raw", ignore.stdout=T,ignore.stderr=T)
	    system("rm reordered_base", ignore.stdout=T,ignore.stderr=T)
	    system("rm target_no_mhc*", ignore.stdout=T,ignore.stderr=T)
	    system("rm ANCESTRY_INFORMATIVE_DIMENSIONS.cluster*", ignore.stdout=T,ignore.stderr=T)
	    system("rm reordered_base", ignore.stdout=T,ignore.stderr=T)
	    system("rm ANCESTRY_INFORMATIVE_DIMENSIONS.log", ignore.stdout=T,ignore.stderr=T)
	    system("rm head_disc", ignore.stdout=T,ignore.stderr=T)
	    system("rm LE_SNPs", ignore.stdout=T,ignore.stderr=T)
	    system("rm prune_target*", ignore.stdout=T,ignore.stderr=T)
	    system("rm plink.log", ignore.stdout=T,ignore.stderr=T)
	    system("rm HEADER", ignore.stdout=T,ignore.stderr=T)
	    system("rm PROFILES.hh", ignore.stdout=T,ignore.stderr=T)
	    system("rm PROFILES.log", ignore.stdout=T,ignore.stderr=T)
	    system("rm synonymous_snps*", ignore.stdout=T,ignore.stderr=T)
	    system("rm synonymous_snps", ignore.stdout=T,ignore.stderr=T) 
	    system("rm non_synonymous_snps_only*", ignore.stdout=T,ignore.stderr=T)
	    system("rm mhc.txt", ignore.stdout=T,ignore.stderr=T)
	    system("rm Rplots.pdf", ignore.stdout=T,ignore.stderr=T)
	    system("rm temp.raw", ignore.stdout=T,ignore.stderr=T)
	    system("rm TARGET_SNPs", ignore.stdout=T,ignore.stderr=T)
	    system("rm base_SNPS", ignore.stdout=T,ignore.stderr=T)  
    }
}


if(multiple.base.phenotypes & multiple.target.phenotypes){
	angle.heat <- 45
	if(length(base.phenotypes.names) > 20){
		angle.heat <- 90
	}
  names(out.multibase) <- base.phenotypes.names
  row.names(out.multibase) <- target.phenotypes
  write.table(out.multibase, paste(figname, "ALL_BEST_THRESHOLDS_BASE_AND_TARGET.txt",sep="_"), col.names=T, row.names=T, quote= F)  
  if(ggfig){
    if(!heat.r2){
      trans.multibase <- data.frame(Base.Phenotype = rep(colnames(out.multibase), each = nrow(out.multibase)), Target.Phenotype = row.names(out.multibase), PRS.P.value = unlist(-log10(out.multibase)))
      }
    if(heat.r2){
      trans.multibase <- data.frame(Base.Phenotype = rep(colnames(out.multibase), each = nrow(out.multibase)), Target.Phenotype = row.names(out.multibase), PRS.P.value = unlist(out.multibase))
      }
    tile.plot <- ggplot(trans.multibase, aes(x = Base.Phenotype, y = Target.Phenotype, fill = PRS.P.value)) + geom_tile() +
      xlab("Base Phenotype") +
      ylab("Target Phenotype") + 
      theme(axis.text.x=element_text(angle=angle.heat, hjust=1))
    if(!heat.r2){
    	tile.plot <- tile.plot + 
    	  scale_fill_gradient(bquote(atop(Best~-log[10]~model,italic(P)-value),), low= bar.col.is.pval.lowcol, high= bar.col.is.pval.highcol)
    }
    if(heat.r2){
    	tile.plot <- tile.plot + 
    	  scale_fill_gradient(bquote(Best~R^2), low= bar.col.is.pval.lowcol, high= bar.col.is.pval.highcol)
    }

    tile.plot
    ggsave(filename=paste(figname, "_HEATMAP.png", sep=""),,width=7,height=7)
  }
  if(!ggfig){
  	png(paste(figname, "_HEATMAP.png", sep=""))
  	heatmap(t(data.matrix(-log10(out.multibase))), Rowv=NA, Colv=NA, margins=c(10,10), cexCol=0.8, cexRow=0.8)
  	dev.off()
  }
}



if(cleanup){
cat(" ################################# \n # \n #   Cleanup \n # \n ################################# \n")
  system("rm PROFILES.0.5*.profile", ignore.stdout=T,ignore.stderr=T)	
  system("rm PROFILES.0.4*.profile", ignore.stdout=T,ignore.stderr=T)
  system("rm PROFILES.0.3*.profile", ignore.stdout=T,ignore.stderr=T)
  system("rm PROFILES.0.2*.profile", ignore.stdout=T,ignore.stderr=T)
  system("rm PROFILES.0.1*.profile", ignore.stdout=T,ignore.stderr=T)
  system("rm PROFILES.0.0*.profile", ignore.stdout=T,ignore.stderr=T)
  system("rm PROFILES.*.profile", ignore.stdout=T,ignore.stderr=T)
  system("rm PROFILES.S1*.profile", ignore.stdout=T,ignore.stderr=T)	
  system("rm PROFILES.S2*.profile", ignore.stdout=T,ignore.stderr=T)
  system("rm PROFILES.S3*.profile", ignore.stdout=T,ignore.stderr=T)
  system("rm PROFILES.S4*.profile", ignore.stdout=T,ignore.stderr=T)
  system("rm PROFILES.S5*.profile", ignore.stdout=T,ignore.stderr=T)
  system("rm PROFILES.S6*.profile", ignore.stdout=T,ignore.stderr=T)
  system("rm PROFILES.S7*.profile", ignore.stdout=T,ignore.stderr=T)
  system("rm PROFILES.S8*.profile", ignore.stdout=T,ignore.stderr=T)
  system("rm PROFILES.S9*.profile", ignore.stdout=T,ignore.stderr=T)
  system("rm *PROFILES.S*.profile", ignore.stdout=T,ignore.stderr=T)
  system("rm *PROFILES.0.5*.profile", ignore.stdout=T,ignore.stderr=T)	
  system("rm *PROFILES.0.4*.profile", ignore.stdout=T,ignore.stderr=T)
  system("rm *PROFILES.0.3*.profile", ignore.stdout=T,ignore.stderr=T)
  system("rm *PROFILES.0.2*.profile", ignore.stdout=T,ignore.stderr=T)
  system("rm *PROFILES.0.1*.profile", ignore.stdout=T,ignore.stderr=T)
  system("rm *PROFILES.0.0*.profile", ignore.stdout=T,ignore.stderr=T)
  system("rm *PROFILES.*.profile", ignore.stdout=T,ignore.stderr=T)
  system("rm *PROFILES.S1*.profile", ignore.stdout=T,ignore.stderr=T)	
  system("rm *PROFILES.S2*.profile", ignore.stdout=T,ignore.stderr=T)
  system("rm *PROFILES.S3*.profile", ignore.stdout=T,ignore.stderr=T)
  system("rm *PROFILES.S4*.profile", ignore.stdout=T,ignore.stderr=T)
  system("rm *PROFILES.S5*.profile", ignore.stdout=T,ignore.stderr=T)
  system("rm *PROFILES.S6*.profile", ignore.stdout=T,ignore.stderr=T)
  system("rm *PROFILES.S7*.profile", ignore.stdout=T,ignore.stderr=T)
  system("rm *PROFILES.S8*.profile", ignore.stdout=T,ignore.stderr=T)
  system("rm *PROFILES.S9*.profile", ignore.stdout=T,ignore.stderr=T)
  system("rm *PROFILES.S*.profile", ignore.stdout=T,ignore.stderr=T) 
  system("rm flip*", ignore.stdout=T,ignore.stderr=T)
  system("rm clean*", ignore.stdout=T,ignore.stderr=T)
  system("rm profile_list", ignore.stdout=T,ignore.stderr=T)
  system("rm rangelist.txt", ignore.stdout=T,ignore.stderr=T)
  system("rm rangelist_ranges", ignore.stdout=T,ignore.stderr=T)
  system("rm rawfile.raw", ignore.stdout=T,ignore.stderr=T)
  system("rm reordered_base", ignore.stdout=T,ignore.stderr=T)
  system("rm target_no_mhc*", ignore.stdout=T,ignore.stderr=T)
  system("rm ANCESTRY_INFORMATIVE_DIMENSIONS.cluster*", ignore.stdout=T,ignore.stderr=T)
  system("rm reordered_base", ignore.stdout=T,ignore.stderr=T)
  system("rm ANCESTRY_INFORMATIVE_DIMENSIONS.log", ignore.stdout=T,ignore.stderr=T)
  system("rm head_disc", ignore.stdout=T,ignore.stderr=T)
  system("rm LE_SNPs", ignore.stdout=T,ignore.stderr=T)
  system("rm prune_target*", ignore.stdout=T,ignore.stderr=T)
  system("rm plink.log", ignore.stdout=T,ignore.stderr=T)
  system("rm HEADER", ignore.stdout=T,ignore.stderr=T)
  system("rm PROFILES.hh", ignore.stdout=T,ignore.stderr=T)
  system("rm PROFILES.log", ignore.stdout=T,ignore.stderr=T)
  system("rm synonymous_snps*", ignore.stdout=T,ignore.stderr=T)
  system("rm synonymous_snps", ignore.stdout=T,ignore.stderr=T) 
  system("rm non_synonymous_snps_only*", ignore.stdout=T,ignore.stderr=T)
  system("rm mhc.txt", ignore.stdout=T,ignore.stderr=T)
  system("rm Rplots.pdf", ignore.stdout=T,ignore.stderr=T)
  system("rm temp.raw", ignore.stdout=T,ignore.stderr=T)
  system("rm TARGET_SNPs", ignore.stdout=T,ignore.stderr=T)
  system("rm base_SNPS", ignore.stdout=T,ignore.stderr=T)  
}


}

# FY: won't use it for now
if(sumsum){
  cat(" ################################# \n # \n #   Using summary statistic - summary statistic method of Johnson et al \n # \n ################################# \n")

if(plink.silent){
  plink <- paste(plink, " --silent ")
}

#supper <- supper-sinc
#slower <- slower+sinc


library(gtx)
if(ggfig){
  library(ggplot2)  
  library(plyr)
}
options(stringsAsFactors=F)

targ.gwas <- read.table(target,head=T)

if(is.na(size.targ)){
  cat("ERROR: Please supply sample size for Target GWAS \n Quitting \n")
  quit()
}

if(clump.snps){
	if(is.na(clump.ref)){
		cat("ERROR: No reference panel supplied for clumping \n Quitting \n")
		quit()
	}
system(paste("awk '{print $2}' ", clump.ref, ".bim | sort -k1,1 > clump_panel_snps.txt", sep = ""))
system(paste("head -n +1 ", base, " > head.base"))
system(paste("head -n +1 ", target, " > head.targ"))

base.names <- read.table("head.base",head=T)
targ.names <- read.table("head.targ",head=T)

system(paste("tail -n +2 ", base, " | awk '{print $", which(colnames(base.names) == "SNP"), "}' | sort -k1,1 > base_snps",sep="")) 
system(paste("tail -n +2 ", target, " | awk '{print $", which(colnames(targ.names) == "SNP"), "}' | sort -k1,1 > targ_snps",sep="")) 
system(paste("join -1 1 -2 1 clump_panel_snps.txt base_snps | join -1 1 -2 1 ``-'' targ_snps > clean_snps", sep = ""))


if("BETA" %in% colnames(base.names)){
  system(paste("tail -n +2 ", base, " | awk '{print $", which(colnames(base.names) == "SNP"),",$", which(colnames(base.names) == "A1"),",$", which(colnames(base.names) == "A2"),",$", which(colnames(base.names) == "BETA"),",$", which(colnames(base.names) == "P")," }' > reordered_base", sep=""))
  system("echo SNP A1 A2 BETA P > HEADER")
}
if(!("BETA" %in% colnames(base.names))){
  system(paste("tail -n +2 ", base, " | awk '{print $", which(colnames(base.names) == "SNP"),",$", which(colnames(base.names) == "A1"),",$", which(colnames(base.names) == "A2"),",$", which(colnames(base.names) == "OR"),",$", which(colnames(base.names) == "P")," }' > reordered_base", sep=""))
  system("echo SNP A1 A2 OR P > HEADER")
}


system(paste("sort -k1,1 reordered_base | join -1 1 -2 1 ``-'' clean_snps | cat HEADER ``-'' > cleaned_base", sep = ""))


system(paste(plink," --noweb --bfile ", clump.ref, "  --clump cleaned_base --clump-p1 ",clump.p1," --clump-p2 ",clump.p2," --clump-r2 ",clump.r2," --clump-kb ",clump.kb," --out clumped_base", sep=" "))
			system("tail -n +2 clumped_base.clumped | awk '{print $3}'  | sort -k1,1 | awk '($1 != \"\"){print}'  > LE_SNPs")


le.snps <- read.table("LE_SNPs",head=F)
targ.gwas <- targ.gwas[targ.gwas$SNP %in% le.snps$V1 , ]
}
if(!clump.snps){
  system(paste("cp ", base, " cleaned_base"))
}
base.gwas <- read.table("cleaned_base",head=T)
targ.gwas <- targ.gwas[targ.gwas$SNP %in% base.gwas$SNP , ]
base.gwas <- base.gwas[base.gwas$SNP %in% targ.gwas$SNP , ]

if("OR" %in% names(targ.gwas)){
  targ.gwas$BETA <- log(targ.gwas$OR)
}
if("OR" %in% names(base.gwas)){
  base.gwas$BETA <- log(base.gwas$OR)
}

cat(" ################################# \n # \n #   Check input format \n # \n ################################# \n")

if(!("SNP" %in% names(base.gwas))){
  cat("ERROR: No SNP Name Column in Base GWAS \n Quitting \n")
  quit()
}
if(!("P" %in% names(base.gwas))){
  cat("ERROR: No P-value Column in Base GWAS \n Quitting \n")
  quit()
}
if(!("A1" %in% names(base.gwas))){
  cat("ERROR: No A1 Column in Base GWAS \n Quitting \n")
  quit()
}
if(!("A2" %in% names(base.gwas))){
  cat("ERROR: No A2 Column in Base GWAS \n Quitting \n")
  quit()
}
if(!("SNP" %in% names(targ.gwas))){
  cat("ERROR: No SNP Name Column in Target GWAS \n Quitting \n")
  quit()
}
if(!("SE" %in% names(targ.gwas))){
  cat("ERROR: No SE Column in Target GWAS \n Quitting \n")
  quit()
}
if(!("A1" %in% names(targ.gwas))){
  cat("ERROR: No A2 Column in Target GWAS \n Quitting \n")
  quit()
}
if(!("A2" %in% names(targ.gwas))){
  cat("ERROR: No A2 Column in Target GWAS \n Quitting \n")
  quit()
}

if(!("OR" %in% names(base.gwas) | "BETA" %in% names(base.gwas))){
  cat("ERROR: No Effect Size Column in Base GWAS \n ie OR or BETA \n Quitting \n")
  quit()
}
if(!("OR" %in% names(targ.gwas) | "BETA" %in% names(targ.gwas))){
  cat("ERROR: No Effect Size Column in Target GWAS \n ie OR or BETA \n Quitting \n")
  quit()
}



base.gwas <- base.gwas[order(base.gwas$SNP) , ]
targ.gwas <- targ.gwas[order(targ.gwas$SNP) , ]

base.gwas$BETA[base.gwas$A1 == targ.gwas$A2] <- -1*base.gwas$BETA[base.gwas$A1 == targ.gwas$A2]
base.gwas <- base.gwas[base.gwas$A1 == targ.gwas$A1 | base.gwas$A1 == targ.gwas$A2  , ]


#high.res <- seq(slower, supper, sinc)
high.res <- slower*sinc^(0:log(supper/slower,sinc))
pval.out <- as.vector(1)
r2.out <- as.vector(1)
r2full<-as.vector(1);
r2null<-as.vector(1)
nsnps <- as.vector(1)

high.res <- c(high.res, barchart.levels)
high.res <- sort(high.res)
high.res <- high.res[!duplicated(high.res)]

for(i in 1:length(high.res)){
  base.gwas.temp <- base.gwas[base.gwas$P < high.res[i] , ]
  targ.gwas.temp <- targ.gwas[targ.gwas$SNP %in% base.gwas.temp$SNP , ]
#  pval.out[i] <- grs.summary(w = base.gwas.temp$BETA, b = targ.gwas.temp$BETA, s = targ.gwas.temp$SE, n = size.targ)$pval
  pval.out[i] <- pnorm(abs(grs.summary(w = base.gwas.temp$BETA, b = targ.gwas.temp$BETA, s = targ.gwas.temp$SE, n = size.targ)$ahat/ grs.summary(w = base.gwas.temp$BETA, b = targ.gwas.temp$BETA, s = targ.gwas.temp$SE, n = size.targ)$aSE), lower.tail = F)
  r2.out[i] <- grs.summary(w = base.gwas.temp$BETA, b = targ.gwas.temp$BETA, s = targ.gwas.temp$SE, n = size.targ)$R2rs
  r2full[i] <- grs.summary(w = base.gwas.temp$BETA, b = targ.gwas.temp$BETA, s = targ.gwas.temp$SE, n = size.targ)$R2rs
  r2null[i] <- grs.summary(w = base.gwas.temp$BETA, b = targ.gwas.temp$BETA, s = targ.gwas.temp$SE, n = size.targ)$R2rs
  nsnps[i] <- dim(base.gwas.temp)[1]
}

output <- data.frame(high.res, pval.out, r2.out, nsnps)
output$pval.out[output$pval.out == 0 ] <- 2.2e-16


temp.output <- output
names(temp.output) <- c("thresh", "pval", "r2", "nsnps")
temp.output <- temp.output[!duplicated(temp.output$thresh),]
temp.output <- temp.output[!is.na(temp.output$pval),]
write.table(temp.output, paste(figname, "RAW_RESULTS_DATA.txt", sep = "_"), col.names=T, row.names=F,quote=F)

cat(" ################################# \n # \n #   High Density Plots \n # \n ################################# \n")

	      if(best.thresh.on.bar){
	        #add max thresh
	        barchart.levels <- c(barchart.levels.old, output$thresh[which.min(output$p.out)])
	        barchart.levels <- barchart.levels[!duplicated(barchart.levels)]	
	        barchart.levels <- sort(barchart.levels, decreasing = F)
	      }

if(ggfig){
  ggfig.points <- ggplot(data = output) + geom_point(aes(x = high.res, y = -log10(pval.out))) +
	            geom_line(aes(x = high.res, y = -log10(pval.out))) +
                xlab(expression(italic(P)-value~threshold~(italic(P)[T])))  +
                ylab(bquote(Evidence~For~Shared~Genetic~Architecture:~italic(P)-value~(-log[10]))) +
	            geom_line(aes(high.res,  -log10(pval.out)), colour = "green",
	            subset = .(high.res %in% barchart.levels))  +
	            geom_point(aes(high.res,  -log10(pval.out)), colour = "green",
	            subset = .(high.res %in% barchart.levels))  + 
	          theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(), panel.background = element_blank(), 
          axis.line = element_line(colour = "black",size=0.5))
	          ggfig.points
	          ggsave(filename=paste(figname,"_HIGH-RES_PLOT_", Sys.Date(), ".png", sep = ""),width=7,height=7)
}
if(!ggfig){
        png(paste(figname,"_HIGH-RES_PLOT_", Sys.Date(), ".png", sep = ""))
	          with(temp.output, plot(	x = temp.output$thresh,
	          y = -log10(temp.output$pval),
	          xlab = expression(italic(P)-value~threshold),
              ylab = bquote(Evidence~For~Shared~Genetic~Architecture:~italic(P)-value~(-log[10])),
	          pch = 19,
	          type = "b"
	        ))
	        with(temp.output[temp.output$thresh %in% barchart.levels,], points(thresh, -log10(pval), col = "green", pch = 19, type= "b"))
}

cat(" ################################# \n # \n #   Barplots \n # Bars for inclusion can be changed using the barchart.levels option \n ################################# \n")

barchart.levels <- c(barchart.levels.old, output$high.res[which.min(output$pval.out)])
barchart.levels <- sort(barchart.levels)

output <- output[!duplicated(output$high.res),]
output <- output[!is.na(output$pval.out),]


options(scipen=0,digits=7)
output <- output[output$high.res %in% barchart.levels,]
output$print.p[round(output$pval.out, digits = 3) != 0] <- round(output$pval.out[round(output$pval.out, digits = 3) != 0], digits = 3)
output$print.p[round(output$pval.out, digits = 3) == 0] <- format(output$pval.out[round(output$pval.out, digits = 3) == 0], digits=2)
output$print.p <- sub("e", "*x*10^", output$print.p)
if(ggfig){
  ggfig.plot <-   ggplot(data = output)  + geom_bar(aes(x = factor(high.res), y = r2.out, fill = -log10(pval.out)), stat="identity") +     
  scale_fill_gradient(low= bar.col.is.pval.lowcol, high= bar.col.is.pval.highcol, name =bquote(atop(-log[10]~model,italic(P)-value),)) +  
  xlab(expression(italic(P)-value~threshold~(italic(P)[T])))  
  ggfig.plot <- ggfig.plot + geom_text(aes(x = factor(high.res), y = r2.out, label = paste(print.p)), vjust = -1.5, hjust = 0, angle = 45, cex = 2.8, parse=T) +
  scale_y_continuous(limits = c(0, max(output$r2.out)*1.25))  + 
	          theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(), panel.background = element_blank(), 
          axis.line = element_line(colour = "black",size=0.5))
  if(binary.target){
	ggfig.plot <- ggfig.plot + ylab(bquote(Variance~Explained:~Pseudo~R^2)) 
  } 
  if(binary.target){
	ggfig.plot <- ggfig.plot + ylab(bquote(Variance~Explained:~R^2)) 
  } 
  ggfig.plot
  ggsave(filename=paste(figname,"_BARPLOT_", Sys.Date(), ".png", sep = ""),width=7,height=7)		
}
if(!ggfig){
        png(paste(figname,"_BARPLOT_", Sys.Date(), ".png", sep = ""))
       if(binary.target){
       plot.fig  <- with(output, barplot(r2.out, 
         names = high.res,
         main = "", 
         col = "red",
         xlab = expression(italic(P)[T]),
         ylab = bquote(Pseudo~R^2),
         ylim = c(0,  max(output$r2.out)*1.25)  ))
          text( parse(text=paste(
         output$print.p)), 
         x = plot.fig+0.1, 
        y =  output$r2.out+ (max(output$r2.out)*1.125-max(output$r2.out)), 
         srt = 45)
       }
       if(!binary.target){
       plot.fig  <- with(output, barplot(r2.out, 
         names = high.res,
         main = "", 
         col = "red",
         xlab = expression(italic(P)[T]),
         ylab = expression(R^2),
         ylim = c(0,  max(output$r2.out)*1.25)  ))
          text( parse(text=paste(
         output$print.p)), 
         x = plot.fig+0.1, 
        y =  output$r2.out+ (max(output$r2.out)*1.125-max(output$r2.out)), 
         srt = 45)
       }
       dev.off()
}

if(cleanup){
	system("rm cleaned_base_temp", ignore.stdout=T,ignore.stderr=T)
	system("rm clump_panel_snps.txt", ignore.stdout=T,ignore.stderr=T)
	system("rm clumped_base.*", ignore.stdout=T,ignore.stderr=T)
	system("rm head.base", ignore.stdout=T,ignore.stderr=T)
	system("rm head.targ", ignore.stdout=T,ignore.stderr=T)
	system("rm targ_snps", ignore.stdout=T,ignore.stderr=T)
	system("rm HEADER", ignore.stdout=T,ignore.stderr=T)
	system("rm cleaned_base", ignore.stdout=T,ignore.stderr=T)
	system("rm clean_snps", ignore.stdout=T,ignore.stderr=T)
	system("rm base_snps", ignore.stdout=T,ignore.stderr=T)
	system("rm LE_SNPs", ignore.stdout=T,ignore.stderr=T)
}
	

}


if(print.time){
  cat(" ################################# \n # \n #   Print time \n # \n ################################# \n")
}

running.time <- proc.time()[3] - start.time 
if(running.time < 60){
	out.run.time <- round(running.time, digits = 2)
	out.time.units <- "seconds"
}
if(running.time > 60 & running.time < 3600){
	out.run.time <- round((running.time/60), digits = 2)
	out.time.units <- "minutes"
}
if(running.time > 3600){
	out.run.time <- round((running.time/3600), digits = 2)
	out.time.units <- "hours"
}


if(print.time){
  cat(paste("RUNNING TIME: ", 	out.run.time, out.time.units, "\n", sep = " "))
}



