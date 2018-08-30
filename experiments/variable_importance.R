# clinicalnona

matchday <- match(clinicalnona$samplecode,allsamples$samplecode)
days <- allsamples$days[matchday]
y <- factor(days >= median(days))
X <- clinicalnona[,-1]
df <- data.frame(y,X)
names(df) <- c("y", "pathology_T","pathology_N", "pathology_stage", "smoking_history", "age", "gender", "radiation_therapy", "molecular_therapy")
fit = randomForest(y~., data=df)
isdead <- clinicalnona$samplecode %in% deadsamples$samplecode
survobj <- Surv(days, isdead)
# library('caret')
# varImp(fit)

varImpPlot(fit,type=2)
# 
# importanceOrder=order(-fit$importance)
# names=rownames(fit$importance)[importanceOrder][1:15]
# par(mfrow=c(4, 2), xpd=NA)
# for (name in names)
#    partialPlot(fit, df, eval(name), main=name, xlab=name,ylim=c(-.2,.9))

library(rpart)
fit=rpart(factor(y)~., df)
plot(fit)
text(fit)

# marsModel <- earth(y ~ ., data=df) # build model
# ev <- evimp (marsModel) # estimate variable importance


# dummy variables
names(X) <-  c("pathology_T","pathology_N", "pathology_stage", "smoking_history", "age", "gender", "radiation_therapy", "molecular_therapy")
temp <-  as.data.frame(sapply(X, function(x) factor(x)))
clibinaries <- dummy.data.frame(temp)

matchday <- match(clinicalnona$samplecode,allsamples$samplecode)
days <- allsamples$days[matchday]
y <- factor(days >= median(days))
df <- data.frame(y,clibinaries)
fit = randomForest(y~., data=df, importance=TRUE)
isdead <- clinicalnona$samplecode %in% deadsamples$samplecode
survobj <- Surv(days, isdead)
# library('caret')
# varImp(fit)

varImpPlot(fit,type=2,pch=19, col=1, cex=.8, main="")
# 
# importanceOrder=order(-fit$importance)
# names=rownames(fit$importance)[importanceOrder][1:15]
# par(mfrow=c(4, 2), xpd=NA)
# for (name in names)
#    partialPlot(fit, df, eval(name), main=name, xlab=name,ylim=c(-.2,.9))

# agenumbers <-clinical$age_at_initial_pathologic_diagnosis[match(clinicalnona$samplecode, clinical$samplecode)] 
X2 <- X
X2$age <- agenumbers
coxpvalues <- sapply(as.data.frame(X2), function(x) coef(summary(coxph(survobj ~ x)))[5])
aucvalues <- sapply(as.data.frame(X2),  function(x)  roc(y,x)$auc)
spearmanvalues <-  sapply(as.data.frame(X2),  function(x)  cor(as.numeric(y),x,method="spearman"))
clinicalvars <- t(rbind(spearmanvalues,aucvalues,coxpvalues))

library(rpart)
fit=rpart(survobj~., X)
plot(fit, compress = TRUE, uniform=TRUE)
text(fit, use.n = TRUE)


