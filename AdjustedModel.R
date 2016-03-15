ds=prof[,c("V2", "SCORE", covariates)]
model.gender<-glm(V2 ~ Gender, family="gaussian",data=ds)
height.res<-residuals(model.gender)
cov.new=c("EV1","EV2","EV3","EV4","EV5","EV6","EV7","EV8","EV9","EV10","Gender")
ds.new=prof[,c("SCORE", cov.new)]
ds.new$V2=height.res
model.logit <- glm(V2  ~ ., family="gaussian", data = ds.new)
model.null <-  glm(V2  ~ ., family="gaussian", data = ds.new[,c("V2", cov.new)])
p.out  <- summary(model.logit)$coefficients[2,4]
r2full <- 1-model.logit$deviance/model.logit$null.deviance
r2null <- 1-model.null$deviance/model.null$null.deviance
r2.out <- r2full - r2null
