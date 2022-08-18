## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  collapse = FALSE,
  comment = "#>"
)

options(rmarkdown.html_vignette.check_title = FALSE)

## ---- eval = FALSE------------------------------------------------------------
#  install.packages("multid")

## ---- eval = FALSE------------------------------------------------------------
#  install.packages("devtools")
#  devtools::install_github("vjilmari/multid")

## ----setup,warning = FALSE, message = FALSE-----------------------------------
library(rio)
library(multid)
library(dplyr)
library(overlapping)

## ----importdata---------------------------------------------------------------
# create a temporary directory
td <- tempdir()

# create a temporary file
tf <- tempfile(tmpdir=td, fileext=".zip")

# download file from internet into temporary location
download.file("http://openpsychometrics.org/_rawdata/BIG5.zip", tf)

# list zip archive
file_names <- unzip(tf, list=TRUE)

# extract files from zip file
unzip(tf, exdir=td, overwrite=TRUE)

# use when zip file has only one file
dat <- import(file.path(td, "BIG5/data.csv"))

# delete the files and directories
unlink(td)

## -----------------------------------------------------------------------------
# save names of personality items to a vector
per.items<-names(dat)[which(names(dat)=="E1"):
                        which(names(dat)=="O10")]

# code item response 0 (zero) to not available (NA) on all these items
dat[,per.items]<-
  sapply(dat[,per.items],function(x){ifelse(x==0,NA,x)})

# check that the numerical range makes sense
range(dat[,per.items],na.rm=T)

# calculate sum scores for Big Five (reverse coded items subtracted from 6)
# see codebook from https://openpsychometrics.org for item content and keys

dat$N<-rowMeans(
  cbind(dat[,c("N1","N3","N5","N6","N7","N8","N9","N10")],
        6-dat[,c("N2","N4")]),na.rm=T)

dat$E<-rowMeans(
  cbind(dat[,c("E1","E3","E5","E7","E9")],
        6-dat[,c("E2","E4","E6","E8","E10")]),na.rm=T)

dat$O<-rowMeans(
  cbind(dat[,c("O1","O3","O5","O7","O8","O9","O10")],
        6-dat[,c("O2","O4","O6")]),na.rm=T)

dat$A<-rowMeans(
  cbind(dat[,c("A2","A4","A6","A8","A9","A10")],
        6-dat[,c("A1","A3","A5","A7")]),na.rm=T)

dat$C<-rowMeans(
  cbind(dat[,c("C1","C3","C5","C7","C9","C10")],
        6-dat[,c("C2","C4","C6","C8")]),na.rm=T)


## -----------------------------------------------------------------------------
US.dat<-
  na.omit(dat[dat$country=="US" & (dat$gender==1 | dat$gender==2),
              c("gender",per.items,"N","E","O","A","C")])

# code categorical gender variables
US.dat$gender.cat<-
  ifelse(US.dat$gender==1,"Male","Female")


## -----------------------------------------------------------------------------
d.N<-d_pooled_sd(data = US.dat,var = "N",
            group.var = "gender.cat",group.values = c("Male","Female"),
            infer=T)

d.E<-d_pooled_sd(data = US.dat,var = "E",
            group.var = "gender.cat",group.values = c("Male","Female"),
            infer=T)

d.O<-d_pooled_sd(data = US.dat,var = "O",
            group.var = "gender.cat",group.values = c("Male","Female"),
            infer=T)

d.A<-d_pooled_sd(data = US.dat,var = "A",
            group.var = "gender.cat",group.values = c("Male","Female"),
            infer=T)

d.C<-d_pooled_sd(data = US.dat,var = "C",
            group.var = "gender.cat",group.values = c("Male","Female"),
            infer=T)

# compile to a table

uni.d<-rbind(d.N,d.E,d.O,d.A,d.C)
rownames(uni.d)<-c("N","E","O","A","C")

round(uni.d,2)


## -----------------------------------------------------------------------------

# let's use the previously obtained univariate d-values
d_vector<-uni.d[,"D"]
# calculate intercorrelations
cor_mat<-cor(US.dat[,c("N","E","O","A","C")])
# test that the ordering of these variables matches
names(d_vector)==colnames(cor_mat)

# print the correlation matrix
round(cor_mat,2)

# Calculate D
D_maha<-sqrt(t(d_vector) %*% solve(cor_mat) %*% d_vector)
D_maha

## -----------------------------------------------------------------------------
# set seed number for reproducibility
set.seed(37127)

# decide the size of the regularization (training data fold)
# here, half of the smaller of the male and female groups is used
size=round(min(table(US.dat$gender.cat))/2)

# run D_regularized
D_out<-
  D_regularized(data=US.dat,
                mv.vars=c("N","E","O","A","C"),
                group.var="gender.cat",
                group.values = c("Male","Female"),
                out = T,
                size=size,
                pcc = T,
                auc = T,
                pred.prob = T)

# print the output
round(D_out$D,2)

## -----------------------------------------------------------------------------
round(100*D_out$P.table,1)

## -----------------------------------------------------------------------------
coefficients(D_out$cv.mod,
             s = "lambda.min")

coefficients(D_out$cv.mod,
             s = "lambda.1se")

plot(D_out$cv.mod)

## ---- fig.height=4, fig.width=8-----------------------------------------------

# Obtain the D-estimate
D<-unname(D_out$D[,"D"])

# Overlap relative to a single distribution (OVL)
OVL<-2*pnorm((-D/2))

# Proportion of overlap relative to the joint distribution
OVL2<-OVL/(2-OVL)

# non-parametric overlap with overlapping package

np.overlap<-
  overlap(x = list(D_out$pred.dat[
  D_out$pred.dat$group=="Male","pred"],
  D_out$pred.dat[
  D_out$pred.dat$group=="Female","pred"]),
  plot=T)

# this corresponds to Proportion of overlap relative to the joint distribution (OVL2)
np.OVL2<-unname(np.overlap$OV)

# from which Proportion of overlap relative to a single distribution (OVL) is approximated at
np.OVL<-(2*np.OVL2)/(1+np.OVL2)

# comparison of different overlap-estimates
round(cbind(OVL,np.OVL,OVL2,np.OVL2),2)


## -----------------------------------------------------------------------------
# number of redraws
reps=10

# placeholder for D-estimates
D_out_boot<-rep(NA,reps)

# repeat the procedure to obtain a distribution

t1<-Sys.time()

for (i in 1:reps){
  # draw replaced samples with the same size as the original data
  boot.dat<-sample_n(US.dat,size=nrow(US.dat),replace=T)
  
  # run D_regularized for each sample
  D_out_boot[i]<-
    D_regularized(data=boot.dat,
                mv.vars=c("N","E","O","A","C"),
              group.var="gender.cat",
              group.values = c("Male","Female"),
              out = T,
              size=size)$D[,"D"]
  
}
t2<-Sys.time()
t2-t1

# print 95% CI percentile interval

round(quantile(D_out_boot,c(.025,.975)),2)


## -----------------------------------------------------------------------------
# place-holder vector for univariate Cohen d's for items
d_vector_items<-rep(NA,length(per.items))

# loop through the items and obtain d for each item


for (i in 1:length(per.items)){
  d_vector_items[i]<-
    d_pooled_sd(data=US.dat,
              var=per.items[i],
            group.var="gender.cat",
            group.values=c("Male","Female"))[,"D"]
  
}

cor_mat_items<-cor(US.dat[,per.items])

# Calculate D
D_maha_items<-
  sqrt(t(d_vector_items) %*% solve(cor_mat_items) %*% d_vector_items)
D_maha_items

## -----------------------------------------------------------------------------

D_items_out<-
  D_regularized(data=US.dat,
                mv.vars=per.items,
              group.var="gender.cat",
              group.values = c("Male","Female"),
              out = T,
              size=size,pcc = T,auc = T,pred.prob = T)

# Item-based D
round(D_items_out$D,2)

# Big Five domain-based D
round(D_out$D,2)

## -----------------------------------------------------------------------------
round(100*D_items_out$P.table,1)
round(100*D_out$P.table,1)

## ---- fig.height=4, fig.width=8-----------------------------------------------

# Obtain the D-estimate
D_items<-unname(D_items_out$D[,"D"])

# Overlap relative to a single distribution (OVL)
OVL_items<-2*pnorm((-D_items/2))

# Proportion of overlap relative to the joint distribution
OVL2_items<-OVL_items/(2-OVL_items)

# non-parametric overlap

np.overlap_items<-
  overlap(x = list(D_items_out$pred.dat[
  D_items_out$pred.dat$group=="Male","pred"],
  D_items_out$pred.dat[
  D_items_out$pred.dat$group=="Female","pred"]),
  plot=T)

# this corresponds to Proportion of overlap relative to the joint distribution (OVL2)
np.OVL2_items<-unname(np.overlap_items$OV)

# from which Proportion of overlap relative to a single distribution (OVL) is approximated at
np.OVL_items<-(2*np.OVL2_items)/(1+np.OVL2_items)

# compare the different overlap-estimates

round(cbind(OVL_items,np.OVL_items,
            OVL2_items,np.OVL2_items),2)

# print Big Five domain based overlaps for comparison

round(cbind(OVL,np.OVL,OVL2,np.OVL2),2)

## -----------------------------------------------------------------------------
sessionInfo()

