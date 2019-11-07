## ---- load_data
diabetes.data <- read.csv(
  file="pima-indians-diabetes-database.csv",
  header=TRUE,
  col.names=c("Pregnancies",
              "Glucose",
              "BloodPressure",
              "SkinThickness",
              "Insulin",
              "BMI",
              "DiabetesPedigreeFunction",
              "Age",
              "Outcome"))

for (v in c("Glucose", "BloodPressure", "SkinThickness", "Insulin", "BMI")) {
  diabetes.data[v][diabetes.data[v] == 0] = NA
}

pMiss <- function(x){sum(is.na(x))/length(x)*100}

library(dplyr)
max.1.missing.diabetes.data <- dplyr::filter(diabetes.data, rowSums(is.na(diabetes.data)) <= 1)

mean.imputation.diabetes.data <- max.1.missing.diabetes.data
for (v in c("Glucose", "BloodPressure", "SkinThickness", "Insulin", "BMI")) {
  mean.imputation.diabetes.data[v][is.na(max.1.missing.diabetes.data[v])] = colMeans(max.1.missing.diabetes.data[v], na.rm = TRUE)
}

clean.diabetes.data <- mean.imputation.diabetes.data
clean.diabetes.data$Insulin <- NULL
clean.diabetes.data$Outcome <- as.factor(clean.diabetes.data$Outcome)
for (v in c("Pregnancies", "Glucose", "BloodPressure", "DiabetesPedigreeFunction", "SkinThickness", "BMI", "Age")) {
  clean.diabetes.data[v] <- as.numeric(clean.diabetes.data[[v]])
}

scaled.diabetes.data <- data.frame(
  scale(clean.diabetes.data[c("Pregnancies",
                              "Glucose",
                              "BloodPressure",
                              "DiabetesPedigreeFunction",
                              "SkinThickness",
                              "BMI",
                              "Age")]),
  Outcome=clean.diabetes.data$Outcome)


## ---- raw_data_histograms

diabetes.raw.data <- read.csv(
  file="pima-indians-diabetes-database.csv",
  header=TRUE,
  col.names=c("Pregnancies",
              "Glucose",
              "BloodPressure",
              "SkinThickness",
              "Insulin",
              "BMI",
              "DiabetesPedigreeFunction",
              "Age",
              "Outcome"))

vars = c("Pregnancies",
         "Glucose",
         "BloodPressure",
         "SkinThickness",
         "Insulin",
         "BMI",
         "DiabetesPedigreeFunction",
         "Age")
par(mfrow=c(2,4), mar=c(2,2,5,2))
for (i in seq_along(vars)) {
  mains = c("Pregnancies",
            "Glucose",
            "Blood Pressure",
            "Skin Thickness",
            "Insulin",
            "BMI",
            "Diabetes Pedigree\n\nFunction",
            "Age")
  if (vars[i] == "Pregnancies") {
    breaks <- 15
  } else {
    breaks <- 30
  }
  hist(diabetes.raw.data[[i]],
       main=strwrap(mains[i], width=20),
       xlab="",
       breaks=breaks,
       cex.axis=2,
       cex.main=2)
}

## ---- clean_data_histograms

vars = c("Pregnancies",
         "Glucose",
         "BloodPressure",
         "SkinThickness",
         "BMI",
         "DiabetesPedigreeFunction",
         "Age")

# 1 2 3 4
#  5 6 7 
layout(matrix(c(1, 0, 1, 5, 2, 5, 2, 6, 3, 6, 3, 7, 4, 7, 4, 0), nrow=2, ncol=8))
for (i in seq_along(vars)) {
  mains = c("Pregnancies",
            "Glucose",
            "Blood Pressure",
            "Skin Thickness",
            "BMI",
            "Diabetes Pedigree\n\nFunction",
            "Age")
  if (vars[i] == "Pregnancies") {
    breaks <- 15
  } else {
    breaks <- 30
  }
  hist(clean.diabetes.data[[i]],
       main=strwrap(mains[i], width=20),
       xlab="",
       breaks=breaks,
       cex.axis=2,
       cex.main=2)
}


## ---- scatter_matrix

pairs(scaled.diabetes.data,
      labels=c("Pregnancies",
               "Glucose",
               "Blood\nPressure",
               "DPF",
               "Skin\nThickness",
               "BMI",
               "Age",
               "Outcome"),
      xaxt="n", yaxt="n",
      cex.labels=1.5,
      row1attop=FALSE)

## ---- correlation_plot

library(corrplot)
cor.mtest <- function(mat, ...) {
  # Get the p-values of the correlation matrix
  mat <- as.matrix(mat)
  n <- ncol(mat)
  p.mat<- matrix(NA, n, n)
  diag(p.mat) <- 0
  for (i in 1:(n - 1)) {
    for (j in (i + 1):n) {
      tmp <- cor.test(mat[, i], mat[, j], ...)
      p.mat[i, j] <- p.mat[j, i] <- tmp$p.value
    }
  }
  colnames(p.mat) <- rownames(p.mat) <- colnames(mat)
  p.mat
}

correlation.matrix <- cor(scaled.diabetes.data[c("Pregnancies",
                                                 "Glucose",
                                                 "BloodPressure",
                                                 "DiabetesPedigreeFunction",
                                                 "SkinThickness",
                                                 "BMI",
                                                 "Age")])

rownames(correlation.matrix) <- colnames(correlation.matrix) <- c("Pregnancies",
                                                                  "Glucose",
                                                                  "Blood Pressure",
                                                                  "DPF",
                                                                  "Skin Thickness",
                                                                  "BMI",
                                                                  "Age")

corrplot(correlation.matrix,
         type="lower",
         order="hclust",
         method="color",
         tl.col="black",
         tl.srt=45,
         tl.cex=1.5,
         cl.cex=1.5,
         cl.lim=c(0,1),
         cl.length=6)

## ---- cross_validation

shuffled.data <- dplyr::slice(scaled.diabetes.data, sample(1:n()))

library(caret)
diabetes.data.with.folds <- scaled.diabetes.data
diabetes.data.with.folds$Fold <- createFolds(scaled.diabetes.data$Outcome, k=10, list=FALSE)

cross_val_model <- function(train_method, ...) {
  models = list()
  acc = list()
  val.acc = list()
  for (f in unique(diabetes.data.with.folds$Fold)) {
    test.data = dplyr::filter(diabetes.data.with.folds, diabetes.data.with.folds$Fold == f)
    train.data = dplyr::filter(diabetes.data.with.folds, diabetes.data.with.folds$Fold != f)
    
    model = train_method(train.data, ...)
    
    train.prob <- predict(model, newdata=train.data, type="response")
    if (is.factor(train.prob)) {
      train.prob <- as.numeric(levels(train.prob))[train.prob]
    }
    train.pred <- round(train.prob)
    train.target <- as.numeric(levels(train.data$Outcome))[train.data$Outcome]
    train.acc <- mean(train.pred == train.target)
    
    test.prob <- predict(model, newdata=test.data, type="response")
    if (is.factor(test.prob)) {
      test.prob <- as.numeric(levels(test.prob))[test.prob]
    }
    test.pred <- round(test.prob)
    test.target <- as.numeric(levels(test.data$Outcome))[test.data$Outcome]
    test.acc <- mean(test.pred == test.target)
    
    models <- c(models, list(model))
    acc <- c(acc, list(train.acc))
    val.acc <- c(val.acc, list(test.acc))
  }
  
  return(list("models" = models,
              "acc" = as.vector(acc, "numeric"),
              "val.acc" = as.vector(val.acc, "numeric")))
}


## ---- logistic_regression

library(MASS)
library(pscl)

train.logistic.regression <- function(train.data, formula) {
  lr <- glm(formula,
            family=binomial(),
            data=train.data)
  
  return(lr)
}

logistic.regression = list()
logistic.regression$baseline <- glm(Outcome ~ Pregnancies + Glucose + BloodPressure + SkinThickness + BMI + DiabetesPedigreeFunction + Age,
                                    family=binomial(),
                                    data=scaled.diabetes.data)

logistic.regression$step <- stepAIC(logistic.regression$baseline,
                                    scope= ~ (Pregnancies + Glucose + BloodPressure + SkinThickness + BMI + DiabetesPedigreeFunction + Age)**2,
                                    trace=0)


logistic.regression <- c(cross_val_model(train.logistic.regression, logistic.regression$step$formula),
                         logistic.regression)

# Extract coefficients and p-values of cross-validated models
logistic.regression$coefs = data.frame()
logistic.regression$p.values = data.frame()
for (model in logistic.regression$models) {
  logistic.regression$coefs <- rbind(logistic.regression$coefs, as.data.frame(matrix(coef(model), nrow=1)))
  logistic.regression$p.values <- rbind(logistic.regression$p.values, as.data.frame(matrix(summary(model)$coefficients[,4], nrow=1)))
}
colnames(logistic.regression$coefs) <- names(coef(logistic.regression$step))
colnames(logistic.regression$p.values) <- names(coef(logistic.regression$step))

get_estimate <- function(param) {
  return(summary(logistic.regression$step)$coefficients[param, 1])
}

get_std_err <- function(param) {
  return(summary(logistic.regression$step)$coefficients[param, 2])
}

get_p_value <- function(param) {
  return(summary(logistic.regression$step)$coefficients[param, 4])
}

step.confint <- confint(logistic.regression$step)
get_conf_int <- function(param, updown) {
  step.confint[param, updown]
}

## ---- logistic_regression_forestplot

library(forestplot)
forestplot(labeltext=rbind(c("Coefficient", "Estimate"),
                           cbind(names(exp(coef(logistic.regression$step))),
                                 prettyNum(exp(coef(logistic.regression$step)),
                                           format="fg",
                                           digits=3))),
           align=c("r", "c"),
           mean=c(NA,exp(coef(logistic.regression$step))),
           lower=c(NA,exp(step.confint[,1])),
           upper=c(NA,exp(step.confint[,2])),
           boxsize=0.2,
           lwd.ci=3,
           zero=1,
           ci.vertices=TRUE,
           xlab="Odds ratios",
           txt_gp=fpTxtGp(label=gpar(cex=1),
                          ticks=gpar(cex=1),
                          xlab=gpar(cex=1)),
           col=fpColors(box="gray50", lines="black"))

## ---- random_forest

forest = list()

forest$formula <- Outcome ~ Pregnancies + Glucose + BloodPressure + SkinThickness + BMI + DiabetesPedigreeFunction + Age
# forest$formula <- logistic.regression$step$formula

forest$baseline <- randomForest(forest$formula,
                                data=scaled.diabetes.data,
                                importance=TRUE,
                                proximity=TRUE)

train_rf <- function(train.data, formula, ntree=10, maxnodes=NULL) {
  return(randomForest(formula,
                      data=scaled.diabetes.data,
                      ntree=ntree,
                      maxnodes=maxnodes,
                      importance=TRUE))
}

forest <- c(forest,
            cross_val_model(train_rf,
                            forest$formula,
                            ntree=10,
                            maxnodes=50))

## ---- rf_partial_dependence_plots

# Plots $\tilde{f}(x) = \frac{1}{n} \sum_{i=1}^n f(x, x_{iC})$, with $f(x) = \log p_k(x) - \frac{1}{K} \sum_{j=1}^K \log p_j(x)$. $K$ is the total number of classes (in this case $K=2$), $k$ is the class of interest (i.e. the for which we are measuring the dependence), $n$ is the number of data points, $x$ is the variable whoses dependence we are measuring, and $x_{iC}$ are the other variables.

layout(matrix(c(1, 0, 1, 5, 2, 5, 2, 6, 3, 6, 3, 7, 4, 7, 4, 0), nrow=2, ncol=8))
for (i in seq_along(vars)) {
  if (forest$formula == logistic.regression$baseline$formula) {
    mains = c("Pregnancies",
              "Glucose",
              "Blood Pressure",
              "Skin Thickness",
              "BMI",
              "Diabetes Pedigree\n\nFunction",
              "Age")
  } else {
    mains = c("Pregnancies",
              "Glucose",
              "BMI",
              "DPF",
              "Age",
              "Glucose-DPF",
              "Pregnancies-Age")
  }
  par(mar=c(2,3,4,2))
  mains = strsplit(as.character(forest$formula[3]), split=" ")[[1]][seq(1,13,2)]
  partialPlot(forest$baseline,
              pred.data=scaled.diabetes.data,
              x.var=mains[i],
              which.class=1,
              xlab="",
              main=mains[i])
}

## ---- logistic_regression_roc_curve

library(ROCR)
prob <- predict(logistic.regression$finalModel, newdata=scaled.diabetes.data, type="response")
pred <- prediction(prob, scaled.diabetes.data$Outcome)

auc <- performance(pred, measure = "auc")
auc <- auc@y.values[[1]]

auc_perf <- performance(pred, measure="tpr", x.measure="fpr")
plot(auc_perf, xlab="1 - Specificity", ylab="Sensitivity")
