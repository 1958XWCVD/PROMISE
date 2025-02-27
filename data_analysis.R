

#############################################################
# output
#############################################################
blue="#00A6C7";red="#C30B21";purple="#886FB8";green="#69B843"


paths = c("C:/Users/Administrator/Desktop/RICA/")
path = paths[dir.exists(paths)]
code.path = paste0(path, "codes/")
setwd(code.path)

library(dplyr)
library(survival)
library(lubridate)
library(rms)
library(MASS)
library(xlsx)
source("data_analysis_utilities.R")

data.path = paste0(path, "data/")
out.path = paste0(path, "results/")


############################################
# read data
############################################
treat.code = 2
control.code = ifelse(treat.code==1, 2, 1)

# read data
setwd(data.path)
data.ori = read.xlsx("data.xlsx", sheetIndex=1)
data = data.ori
data$FBG = as.numeric(data$FBG)
data$ABCD2_score = as.numeric(data$ABCD2_score)
data$mRS_score = as.numeric(data$mRS_score)

# prepare data 
cols_name = colnames(data)
event_names=c("CIS")
dt.new = data_used(data, cols_name, event_names)

# if the event of a subject occurs after 1 year, make it censoring
idx.1y = which(dt.new$time>1 & dt.new$event==1)
dt.new$event[idx.1y] = 0
dt.new$time[dt.new$time>1] = 1
dt.new$time = dt.new$time*12
dt.new = reformat_data(dt.new)

# split data
dt.modeling.overall.train = dt.new[which(dt.new$trainset==1),]
dt.modeling.overall.valid = dt.new[which(dt.new$trainset==0),]
dt.modeling.treat.train = dt.modeling.overall.train[which(dt.modeling.overall.train$group_id==treat.code),]
dt.modeling.treat.valid = dt.modeling.overall.valid[which(dt.modeling.overall.valid$group_id==treat.code),]
dt.modeling.control.train = dt.modeling.overall.train[which(dt.modeling.overall.train$group_id==control.code),]
dt.modeling.control.valid = dt.modeling.overall.valid[which(dt.modeling.overall.valid$group_id==control.code),]

setwd(out.path)
if( file.exists("results.xlsx") ) file.remove("results.xlsx")

############################################
# baseline
############################################

# overall table
overall.train.tab = tab_baseline(dt.modeling.overall.train, dig=1, dig.P=3)
overall.valid.tab = tab_baseline(dt.modeling.overall.valid, dig=1, dig.P=3)
overall.diff.overall.test = diff_test_overall_all(dt.modeling.overall.train, dt.modeling.overall.valid, dig.P=3)
baseline.overall.tab = cbind(overall.train.tab, overall.valid.tab, overall.diff.overall.test)

# treatment table
treat.train.tab = tab_baseline(dt.modeling.treat.train, dig=1, dig.P=3)
treat.valid.tab = tab_baseline(dt.modeling.treat.valid, dig=1, dig.P=3)
treat.diff.overall.test = diff_test_overall_all(dt.modeling.treat.train, dt.modeling.treat.valid, dig.P=3)
baseline.treat.tab = cbind(treat.train.tab, treat.valid.tab, treat.diff.overall.test)

# control table
control.train.tab = tab_baseline(dt.modeling.control.train, dig=1, dig.P=3)
control.valid.tab = tab_baseline(dt.modeling.control.valid, dig=1, dig.P=3)
control.diff.overall.test = diff_test_overall_all(dt.modeling.control.train, dt.modeling.control.valid, dig.P=3)
baseline.control.tab = cbind(control.train.tab, control.valid.tab, control.diff.overall.test)

# output
empty.line = rep("----------", ncol(baseline.overall.tab))
baseline.stat = rbind(baseline.overall.tab, empty.line, baseline.treat.tab, empty.line, baseline.control.tab)

write.xlsx(baseline.overall.tab, file="results.xlsx", sheetName="baseline_overall")
write.xlsx(baseline.treat.tab, file="results.xlsx", sheetName="baseline_treat", append=TRUE)
write.xlsx(baseline.control.tab, file="results.xlsx", sheetName="baseline_control", append=TRUE)


############################################
# modeling
############################################

# overall HR
model.overall.tab = modeling(dt.modeling.overall.train, backward=TRUE)

# treatment HR
model.treat.tab = modeling(dt.modeling.treat.train, backward=TRUE)

# control HR
model.control.tab = modeling(dt.modeling.control.train, backward=TRUE)

# output
write.xlsx(model.overall.tab, file="results.xlsx", sheetName="HR_overall", append=TRUE)
write.xlsx(model.treat.tab, file="results.xlsx", sheetName="HR_treat", append=TRUE)
write.xlsx(model.control.tab, file="results.xlsx", sheetName="HR_control", append=TRUE)


############################################
# C-statistics
############################################

library(boot)
library(compareC)
dig = 2
dig.P = 4
B = 1000
vars.all = c("age", "BMI", "hypertension", "LDLC", "diabetes", "smoke", "zeren", "xiazhai")
C.index.overall.train = C_index_CI(dt.modeling.overall.train, vars.all, B=B, dig=dig, dig.P=dig.P)
C.index.treat.train = C_index_CI(dt.modeling.treat.train, vars.all, B=B, dig=dig, dig.P=dig.P)
C.index.control.train = C_index_CI(dt.modeling.control.train, vars.all, B=B, dig=dig, dig.P=dig.P)
C.index.overall.valid = C_index_CI(dt.modeling.overall.valid, vars.all, B=B, dig=dig, dig.P=dig.P)
C.index.treat.valid = C_index_CI(dt.modeling.treat.valid, vars.all, B=B, dig=dig, dig.P=dig.P)
C.index.control.valid = C_index_CI(dt.modeling.control.valid, vars.all, B=B, dig=dig, dig.P=dig.P)
C.overall.tabs = cbind(C.index.overall.train, C.index.treat.train, C.index.control.train,
			C.index.overall.valid, C.index.treat.valid, C.index.control.valid)
write.xlsx(C.overall.tabs, file="results.xlsx", sheetName="C_overall", append=TRUE)


############################################
# An extension of Hosmer - Lemeshow test
############################################

library(survMisc)
var.names = c("age", "BMI", "hypertension", "LDLC", "diabetes", "smoke", "zeren", "xiazhai")
formula = as.formula(paste("Surv(time, event)~", paste(var.names, collapse="+")))

# overall
dat = dt.modeling.overall.train
cox = coxph(formula, data=dat, x=TRUE, y=TRUE)
pval.overallgof.overall = survMisc::gof(cox)$lrTest[2,4]

# treat
dat = dt.modeling.treat.train
cox = coxph(formula, data=dat, x=TRUE, y=TRUE)
pval.overallgof.treat = survMisc::gof(cox)$lrTest[2,4]

# control
dat = dt.modeling.control.train
cox = coxph(formula, data=dat, x=TRUE, y=TRUE)
pval.overallgof.control = survMisc::gof(cox)$lrTest[2,4]

HL.tab = c(pval.overallgof.overall, pval.overallgof.treat, pval.overallgof.control)
names(HL.tab) = c("Overall-Model-overall", "Overall-Model-RIC", "Overall-Model-control")
write.xlsx(HL.tab, file="results.xlsx", sheetName="HL-test", append=TRUE)


############################################
# Calibration
############################################
col.names = c("group_id", "time", "event", "age", "BMI", "hypertension", "LDLC", "diabetes", 
		"smoke", "zeren", "xiazhai")
cols = which(colnames(dt.modeling.control.train) %in% col.names)
vars = c("age", "BMI", "hypertension", "LDLC", "diabetes", "smoke", "zeren", "xiazhai")

# set parameters
# u is the month
u=12;B=10;H=5;W=6.5;L=300;res=100;pt=40;lwd=5;lwd.pdf=3
file.type = "tif"
if( u==3 ){
	num.label="A"
}else if( u==6 ){
	num.label="B"
}else if( u==12 ){
	num.label="C"
}else{
	num.label=""
}

# overall
file.name = paste0("calibration_", num.label, "_总模型_总队列_month-", u, ".", file.type)
if( file.type=="pdf" ){
	pdf(file.name, width=W, height=H)
	lwd=lwd.pdf
}else{
	tiff(file.name, width=W*L, height=H*L, pointsize=pt, res=res, compression="lzw")
}
dd=datadist(dt.modeling.overall.train);option=options(datadist="dd")
plot_cal(dt.modeling.overall.train, dt.modeling.overall.valid, vars, u=u, lwd=lwd, col.train=blue, col.val=red, col.ideal=green, num.label=num.label)
dev.off()

# treat
file.name = paste0("calibration_", num.label, "_总模型_RIC队列_month-", u, ".", file.type)
if( file.type=="pdf" ){
	pdf(file.name, width=W, height=H)
	lwd=lwd.pdf
}else{
	tiff(file.name, width=W*L, height=H*L, pointsize=pt, res=res, compression="lzw")
}
dd=datadist(dt.modeling.treat.train);option=options(datadist="dd")
plot_cal(dt.modeling.treat.train, dt.modeling.treat.valid, vars, u=u, lwd=lwd, col.train=blue, col.val=red, col.ideal=green, num.label=num.label)
dev.off()

# control
file.name = paste0("calibration_", num.label, "_总模型_Control队列_month-", u, ".", file.type)
if( file.type=="pdf" ){
	pdf(file.name, width=W, height=H)
	lwd=lwd.pdf
}else{
	tiff(file.name, width=W*L, height=H*L, pointsize=pt, res=res, compression="lzw")
}
dd=datadist(dt.modeling.control.train);option=options(datadist="dd")
plot_cal(dt.modeling.control.train, dt.modeling.control.valid, vars, u=u, lwd=lwd, col.train=blue, col.val=red, col.ideal=green, num.label=num.label)
dev.off()


## merged plot
file.name = paste0("calibration_总模型_总队列.", file.type)
if( file.type=="pdf" ){
	pdf(file.name, width=W, height=H)
	lwd=lwd.pdf
}else{
	tiff(file.name, width=W*L*3, height=H*L, pointsize=pt*1.5, res=res, compression="lzw")
}
par(mfrow=c(1, 3))
g=-0.04-0.02;cex.m=2/3
dd=datadist(dt.modeling.overall.train);option=options(datadist="dd")
plot_cal(dt.modeling.overall.train, dt.modeling.overall.valid, vars, u=3, lwd=lwd, col.train=blue, col.val=red, col.ideal=green, num.label="", no.y.labels=FALSE, merged=TRUE)
mtext("A", side=2, line=0.5+g, at=1.1, adj=1, cex=cex.m, las=2)
plot_cal(dt.modeling.overall.train, dt.modeling.overall.valid, vars, u=6, lwd=lwd, col.train=blue, col.val=red, col.ideal=green, num.label="", no.y.labels=TRUE, merged=TRUE)
mtext("B", side=2, line=0.5+g, at=1.1, adj=1, cex=cex.m, las=2)
plot_cal(dt.modeling.overall.train, dt.modeling.overall.valid, vars, u=12, lwd=lwd, col.train=blue, col.val=red, col.ideal=green, num.label="", no.y.labels=TRUE, merged=TRUE)
mtext("C", side=2, line=0.5+g, at=1.1, adj=1, cex=cex.m, las=2)
dev.off()

file.name = paste0("calibration_总模型_RIC队列.", file.type)
if( file.type=="pdf" ){
	pdf(file.name, width=W, height=H)
	lwd=lwd.pdf
}else{
	tiff(file.name, width=W*L*3, height=H*L, pointsize=pt*1.5, res=res, compression="lzw")
}
par(mfrow=c(1, 3))
g=-0.04-0.02;cex.m=2/3
dd=datadist(dt.modeling.treat.train);option=options(datadist="dd")
plot_cal(dt.modeling.treat.train, dt.modeling.treat.valid, vars, u=3, lwd=lwd, col.train=blue, col.val=red, col.ideal=green, num.label="", no.y.labels=FALSE, merged=TRUE)
mtext("A", side=2, line=0.5+g, at=1.1, adj=1, cex=cex.m, las=2)
plot_cal(dt.modeling.treat.train, dt.modeling.treat.valid, vars, u=6, lwd=lwd, col.train=blue, col.val=red, col.ideal=green, num.label="", no.y.labels=TRUE, merged=TRUE)
mtext("B", side=2, line=0.5+g, at=1.1, adj=1, cex=cex.m, las=2)
plot_cal(dt.modeling.treat.train, dt.modeling.treat.valid, vars, u=12, lwd=lwd, col.train=blue, col.val=red, col.ideal=green, num.label="", no.y.labels=TRUE, merged=TRUE)
mtext("C", side=2, line=0.5+g, at=1.1, adj=1, cex=cex.m, las=2)
dev.off()

file.name = paste0("calibration_总模型_Control队列.", file.type)
if( file.type=="pdf" ){
	pdf(file.name, width=W, height=H)
	lwd=lwd.pdf
}else{
	tiff(file.name, width=W*L*3, height=H*L, pointsize=pt*1.5, res=res, compression="lzw")
}
par(mfrow=c(1, 3))
g=-0.04-0.02;cex.m=2/3
dd=datadist(dt.modeling.control.train);option=options(datadist="dd")
plot_cal(dt.modeling.control.train, dt.modeling.control.valid, vars, u=3, lwd=lwd, col.train=blue, col.val=red, col.ideal=green, num.label="", no.y.labels=FALSE, merged=TRUE)
mtext("A", side=2, line=0.5+g, at=1.1, adj=1, cex=cex.m, las=2)
plot_cal(dt.modeling.control.train, dt.modeling.control.valid, vars, u=6, lwd=lwd, col.train=blue, col.val=red, col.ideal=green, num.label="", no.y.labels=TRUE, merged=TRUE)
mtext("B", side=2, line=0.5+g, at=1.1, adj=1, cex=cex.m, las=2)
plot_cal(dt.modeling.control.train, dt.modeling.control.valid, vars, u=12, lwd=lwd, col.train=blue, col.val=red, col.ideal=green, num.label="", no.y.labels=TRUE, merged=TRUE)
mtext("C", side=2, line=0.5+g, at=1.1, adj=1, cex=cex.m, las=2)
dev.off()


############################################
# DCA
############################################
library(dcurves)
library(ggplot2)

# set parameters
u=12;W=6.5;H=5;L=300;res=100;pt=40;lwd=5;lwd.pdf=3
file.type = "tif"
var.names = c("age", "BMI", "hypertension", "LDLC", "diabetes", "smoke", "zeren", "xiazhai")

# overall train
file.name = paste0("DCA_总模型_总队列-", label_name(u)[1], ".", file.type)
if( file.type=="pdf" ){
	pdf(file.name, width=W, height=H)
}else{
	tiff(file.name, width=W*L, height=H*L, pointsize=pt, res=res, compression="lzw")
}
dat = dt.modeling.overall.train
formula = as.formula(paste("Surv(time, event)~", paste(var.names, collapse="+")))
cox = coxph(formula, data=dat, x=TRUE, y=TRUE) 
plot_DCA(dat, cox, u=u, thresholds=seq(0, 1, 0.01), lwd=lwd, col.all=blue, col.none=red, col.model=green, num.label=label_name(u)[1])
dev.off()

# overall valid
file.name = paste0("DCA_总模型_总队列-", label_name(u)[2], ".", file.type)
if( file.type=="pdf" ){
	pdf(file.name, width=W, height=H)
}else{
	tiff(file.name, width=W*L, height=H*L, pointsize=pt, res=res, compression="lzw")
}
dat = dt.modeling.overall.valid
formula = as.formula(paste("Surv(time, event)~", paste(var.names, collapse="+")))
cox = coxph(formula, data=dat, x=TRUE, y=TRUE)
plot_DCA(dat, cox, u=u, thresholds=seq(0, 1, 0.01), lwd=lwd, col.all=blue, col.none=red, col.model=green, num.label=label_name(u)[2])
dev.off()

# treat train
file.name = paste0("DCA_总模型_RIC队列-", label_name(u)[1], ".", file.type)
if( file.type=="pdf" ){
	pdf(file.name, width=W, height=H)
}else{
	tiff(file.name, width=W*L, height=H*L, pointsize=pt, res=res, compression="lzw")
}
dat = dt.modeling.treat.train
formula = as.formula(paste("Surv(time, event)~", paste(var.names, collapse="+")))
cox = coxph(formula, data=dat, x=TRUE, y=TRUE) 
plot_DCA(dat, cox, u=u, thresholds=seq(0, 1, 0.01), lwd=lwd, col.all=blue, col.none=red, col.model=green, num.label=label_name(u)[1])
dev.off()

# treat valid
file.name = paste0("DCA_总模型_RIC队列-", label_name(u)[2], ".", file.type)
if( file.type=="pdf" ){
	pdf(file.name, width=W, height=H)
}else{
	tiff(file.name, width=W*L, height=H*L, pointsize=pt, res=res, compression="lzw")
}
dat = dt.modeling.treat.valid
formula = as.formula(paste("Surv(time, event)~", paste(var.names, collapse="+")))
cox = coxph(formula, data=dat, x=TRUE, y=TRUE)
plot_DCA(dat, cox, u=u, thresholds=seq(0, 1, 0.01), lwd=lwd, col.all=blue, col.none=red, col.model=green, num.label=label_name(u)[2])
dev.off()


# control train
file.name = paste0("DCA_总模型_Control队列-", label_name(u)[1], ".", file.type)
if( file.type=="pdf" ){
	pdf(file.name, width=W, height=H)
}else{
	tiff(file.name, width=W*L, height=H*L, pointsize=pt, res=res, compression="lzw")
}
dat = dt.modeling.control.train
formula = as.formula(paste("Surv(time, event)~", paste(var.names, collapse="+")))
cox = coxph(formula, data=dat, x=TRUE, y=TRUE) 
plot_DCA(dat, cox, u=u, thresholds=seq(0, 1, 0.01), lwd=lwd, col.all=blue, col.none=red, col.model=green, num.label=label_name(u)[1])
dev.off()

# control valid
file.name = paste0("DCA_总模型_Control队列-", label_name(u)[2], ".", file.type)
if( file.type=="pdf" ){
	pdf(file.name, width=W, height=H)
}else{
	tiff(file.name, width=W*L, height=H*L, pointsize=pt, res=res, compression="lzw")
}
dat = dt.modeling.control.valid
formula = as.formula(paste("Surv(time, event)~", paste(var.names, collapse="+")))
cox = coxph(formula, data=dat, x=TRUE, y=TRUE)
#dev.new(width=7, height=5)
plot_DCA(dat, cox, u=u, thresholds=seq(0, 1, 0.01), lwd=lwd, col.all=blue, col.none=red, col.model=green, num.label=label_name(u)[2])
dev.off()


########################## 
## 3 columns
########################## 

############# merged plot, overall
file.name = paste0("DCA_总模型_总队列.", file.type)
if( file.type=="pdf" ){
	pdf(file.name, width=W, height=H)
}else{
	tiff(file.name, width=W*L*3, height=H*L*2, pointsize=pt*1.5, res=res, compression="lzw")
}

par(mfrow=c(2, 3))
g=-0.09;cex.m=2/3
#### train
dat = dt.modeling.overall.train
formula = as.formula(paste("Surv(time, event)~", paste(var.names, collapse="+")))
cox = coxph(formula, data=dat, x=TRUE, y=TRUE) 
plot_DCA(dat, cox, u=3, thresholds=seq(0, 1, 0.01), lwd=lwd, col.all=blue, col.none=red, col.model=green, num.label="", merged=TRUE, no.ylab=TRUE)
mtext(label_name(3)[1], side=2, at=0.05+0.01, line=1.3+g, cex=cex.m, las=2, adj=0)
mtext("Net benefit", side=2, line=2.5, cex=cex.m)
plot_DCA(dat, cox, u=6, thresholds=seq(0, 1, 0.01), lwd=lwd, col.all=blue, col.none=red, col.model=green, num.label="", no.y.labels=TRUE, merged=TRUE, no.ylab=TRUE)
mtext(label_name(6)[1], side=2, at=0.1+0.015, line=1.3+g, cex=cex.m, las=2, adj=0)
plot_DCA(dat, cox, u=12, thresholds=seq(0, 1, 0.01), lwd=lwd, col.all=blue, col.none=red, col.model=green, num.label="", no.y.labels=TRUE, merged=TRUE, no.ylab=TRUE)
mtext(label_name(12)[1], side=2, at=0.15+0.02, line=1.3+g, cex=cex.m, las=2, adj=0)
#### valid
dat = dt.modeling.overall.valid
formula = as.formula(paste("Surv(time, event)~", paste(var.names, collapse="+")))
cox = coxph(formula, data=dat, x=TRUE, y=TRUE)
plot_DCA(dat, cox, u=3, thresholds=seq(0, 1, 0.01), lwd=lwd, col.all=blue, col.none=red, col.model=green, num.label="", no.ylab=TRUE)
mtext(label_name(3)[2], side=2, at=0.05+0.01, line=1.3+g, cex=cex.m, las=2, adj=0)
mtext("Net benefit", side=2, line=2.5, cex=cex.m)
plot_DCA(dat, cox, u=6, thresholds=seq(0, 1, 0.01), lwd=lwd, col.all=blue, col.none=red, col.model=green, num.label="", no.y.labels=TRUE, no.ylab=TRUE)
mtext(label_name(6)[2], side=2, at=0.1+0.015, line=1.3+g, cex=cex.m, las=2, adj=0)
plot_DCA(dat, cox, u=12, thresholds=seq(0, 1, 0.01), lwd=lwd, col.all=blue, col.none=red, col.model=green, num.label="", no.y.labels=TRUE, no.ylab=TRUE)
mtext(label_name(12)[2], side=2, at=0.15+0.02, line=1.3+g, cex=cex.m, las=2, adj=0)
dev.off()

############# merged plot, RIC
file.name = paste0("DCA_总模型_RIC队列.", file.type)
if( file.type=="pdf" ){
	pdf(file.name, width=W, height=H)
}else{
	tiff(file.name, width=W*L*3, height=H*L*2, pointsize=pt*1.5, res=res, compression="lzw")
}

par(mfrow=c(2, 3))
g=-0.09;cex.m=2/3
####
dat = dt.modeling.treat.train
formula = as.formula(paste("Surv(time, event)~", paste(var.names, collapse="+")))
cox = coxph(formula, data=dat, x=TRUE, y=TRUE) 
plot_DCA(dat, cox, u=3, thresholds=seq(0, 1, 0.01), lwd=lwd, col.all=blue, col.none=red, col.model=green, num.label="", merged=TRUE, no.ylab=TRUE)
mtext(label_name(3)[1], side=2, at=0.05+0.01, line=1.3+g, cex=cex.m, las=2, adj=0)
mtext("Net benefit", side=2, line=2.5, cex=cex.m)
plot_DCA(dat, cox, u=6, thresholds=seq(0, 1, 0.01), lwd=lwd, col.all=blue, col.none=red, col.model=green, num.label="", no.y.labels=TRUE, merged=TRUE, no.ylab=TRUE)
mtext(label_name(6)[1], side=2, at=0.1+0.015, line=1.3+g, cex=cex.m, las=2, adj=0)
plot_DCA(dat, cox, u=12, thresholds=seq(0, 1, 0.01), lwd=lwd, col.all=blue, col.none=red, col.model=green, num.label="", no.y.labels=TRUE, merged=TRUE, no.ylab=TRUE)
mtext(label_name(12)[1], side=2, at=0.15+0.02, line=1.3+g, cex=cex.m, las=2, adj=0)
####
dat = dt.modeling.treat.valid
formula = as.formula(paste("Surv(time, event)~", paste(var.names, collapse="+")))
cox = coxph(formula, data=dat, x=TRUE, y=TRUE)
plot_DCA(dat, cox, u=3, thresholds=seq(0, 1, 0.01), lwd=lwd, col.all=blue, col.none=red, col.model=green, num.label="", no.ylab=TRUE)
mtext(label_name(3)[2], side=2, at=0.05+0.01, line=1.3+g, cex=cex.m, las=2, adj=0)
mtext("Net benefit", side=2, line=2.5, cex=cex.m)
plot_DCA(dat, cox, u=6, thresholds=seq(0, 1, 0.01), lwd=lwd, col.all=blue, col.none=red, col.model=green, num.label="", no.y.labels=TRUE, no.ylab=TRUE)
mtext(label_name(6)[2], side=2, at=0.1+0.015, line=1.3+g, cex=cex.m, las=2, adj=0)
plot_DCA(dat, cox, u=12, thresholds=seq(0, 1, 0.01), lwd=lwd, col.all=blue, col.none=red, col.model=green, num.label="", no.y.labels=TRUE, no.ylab=TRUE)
mtext(label_name(12)[2], side=2, at=0.15+0.02, line=1.3+g, cex=cex.m, las=2, adj=0)
dev.off()

############# merged plot, Control
file.name = paste0("DCA_总模型_Control队列.", file.type)
if( file.type=="pdf" ){
	pdf(file.name, width=W, height=H)
}else{
	tiff(file.name, width=W*L*3, height=H*L*2, pointsize=pt*1.5, res=res, compression="lzw")
}

par(mfrow=c(2, 3))
g=-0.09;cex.m=2/3
#### train
dat = dt.modeling.control.train
formula = as.formula(paste("Surv(time, event)~", paste(var.names, collapse="+")))
cox = coxph(formula, data=dat, x=TRUE, y=TRUE) 
plot_DCA(dat, cox, u=3, thresholds=seq(0, 1, 0.01), lwd=lwd, col.all=blue, col.none=red, col.model=green, num.label="", merged=TRUE, no.ylab=TRUE)
mtext(label_name(3)[1], side=2, at=0.05+0.01, line=1.3+g, cex=cex.m, las=2, adj=0)
mtext("Net benefit", side=2, line=2.5, cex=cex.m)
plot_DCA(dat, cox, u=6, thresholds=seq(0, 1, 0.01), lwd=lwd, col.all=blue, col.none=red, col.model=green, num.label="", no.y.labels=TRUE, merged=TRUE, no.ylab=TRUE)
mtext(label_name(6)[1], side=2, at=0.1+0.015, line=1.3+g, cex=cex.m, las=2, adj=0)
plot_DCA(dat, cox, u=12, thresholds=seq(0, 1, 0.01), lwd=lwd, col.all=blue, col.none=red, col.model=green, num.label="", no.y.labels=TRUE, merged=TRUE, no.ylab=TRUE)
mtext(label_name(12)[1], side=2, at=0.15+0.02, line=1.3+g, cex=cex.m, las=2, adj=0)
#### valid
dat = dt.modeling.control.valid
formula = as.formula(paste("Surv(time, event)~", paste(var.names, collapse="+")))
cox = coxph(formula, data=dat, x=TRUE, y=TRUE)
plot_DCA(dat, cox, u=3, thresholds=seq(0, 1, 0.01), lwd=lwd, col.all=blue, col.none=red, col.model=green, num.label="", no.ylab=TRUE)
mtext(label_name(3)[2], side=2, at=0.05+0.01, line=1.3+g, cex=cex.m, las=2, adj=0)
mtext("Net benefit", side=2, line=2.5, cex=cex.m)
plot_DCA(dat, cox, u=6, thresholds=seq(0, 1, 0.01), lwd=lwd, col.all=blue, col.none=red, col.model=green, num.label="", no.y.labels=TRUE, no.ylab=TRUE)
mtext(label_name(6)[2], side=2, at=0.1+0.015, line=1.3+g, cex=cex.m, las=2, adj=0)
plot_DCA(dat, cox, u=12, thresholds=seq(0, 1, 0.01), lwd=lwd, col.all=blue, col.none=red, col.model=green, num.label="", no.y.labels=TRUE, no.ylab=TRUE)
mtext(label_name(12)[2], side=2, at=0.15+0.02, line=1.3+g, cex=cex.m, las=2, adj=0)
dev.off()


########################## 
## 2 columns
########################## 

############# merged plot, overall
file.name = paste0("DCA_总模型_总队列.", file.type)
if( file.type=="pdf" ){
	pdf(file.name, width=W, height=H)
}else{
	tiff(file.name, width=W*L*2, height=H*L*3, pointsize=pt*1.5, res=res, compression="lzw")
}

par(mfrow=c(3, 2))
g=-0.00;cex.m=2/3
#### train
dat = dt.modeling.overall.train
formula = as.formula(paste("Surv(time, event)~", paste(var.names, collapse="+")))
cox.col.1 = coxph(formula, data=dat, x=TRUE, y=TRUE) 
dat.col.1 = dat
#### valid
dat = dt.modeling.overall.valid
formula = as.formula(paste("Surv(time, event)~", paste(var.names, collapse="+")))
cox.col.2 = coxph(formula, data=dat, x=TRUE, y=TRUE)
dat.col.2 = dat

plot_DCA(dat.col.1, cox.col.1, u=3, thresholds=seq(0, 1, 0.01), lwd=lwd, col.all=blue, col.none=red, col.model=green, num.label="", merged=TRUE, no.ylab=TRUE)
mtext(label_name(3)[1], side=2, at=0.05+0.01, line=1.3+g, cex=cex.m, las=2, adj=0)
mtext("Net benefit", side=2, line=2.5, cex=cex.m)
plot_DCA(dat.col.2, cox.col.2, u=3, thresholds=seq(0, 1, 0.01), lwd=lwd, col.all=blue, col.none=red, col.model=green, num.label="", no.y.labels=TRUE, merged=TRUE, no.ylab=TRUE)
mtext(label_name(3)[2], side=2, at=0.05+0.01, line=1.3+g, cex=cex.m, las=2, adj=0)
plot_DCA(dat.col.1, cox.col.1, u=6, thresholds=seq(0, 1, 0.01), lwd=lwd, col.all=blue, col.none=red, col.model=green, num.label="", merged=TRUE, no.ylab=TRUE)
mtext(label_name(6)[1], side=2, at=0.1+0.015, line=1.3+g, cex=cex.m, las=2, adj=0)
mtext("Net benefit", side=2, line=2.5, cex=cex.m)
plot_DCA(dat.col.2, cox.col.2, u=6, thresholds=seq(0, 1, 0.01), lwd=lwd, col.all=blue, col.none=red, col.model=green, num.label="", no.y.labels=TRUE, merged=TRUE, no.ylab=TRUE)
mtext(label_name(6)[2], side=2, at=0.1+0.015, line=1.3+g, cex=cex.m, las=2, adj=0)
plot_DCA(dat.col.1, cox.col.1, u=12, thresholds=seq(0, 1, 0.01), lwd=lwd, col.all=blue, col.none=red, col.model=green, num.label="", no.ylab=TRUE)
mtext("Net benefit", side=2, line=2.5, cex=cex.m)
mtext(label_name(12)[1], side=2, at=0.15+0.02, line=1.3+g, cex=cex.m, las=2, adj=0)
plot_DCA(dat.col.2, cox.col.2, u=12, thresholds=seq(0, 1, 0.01), lwd=lwd, col.all=blue, col.none=red, col.model=green, num.label="", no.y.labels=TRUE, no.ylab=TRUE)
mtext(label_name(12)[2], side=2, at=0.15+0.02, line=1.3+g, cex=cex.m, las=2, adj=0)
dev.off()


############# merged plot, RIC
file.name = paste0("DCA_总模型_RIC队列.", file.type)
if( file.type=="pdf" ){
	pdf(file.name, width=W, height=H)
}else{
	tiff(file.name, width=W*L*2, height=H*L*3, pointsize=pt*1.5, res=res, compression="lzw")
}

par(mfrow=c(3, 2))
g=-0.00;cex.m=2/3
#### train
dat = dt.modeling.treat.train
formula = as.formula(paste("Surv(time, event)~", paste(var.names, collapse="+")))
cox.col.1 = coxph(formula, data=dat, x=TRUE, y=TRUE)
dat.col.1 = dat
#### valid
dat = dt.modeling.treat.valid
formula = as.formula(paste("Surv(time, event)~", paste(var.names, collapse="+")))
cox.col.2 = coxph(formula, data=dat, x=TRUE, y=TRUE)
dat.col.2 = dat

plot_DCA(dat.col.1, cox.col.1, u=3, thresholds=seq(0, 1, 0.01), lwd=lwd, col.all=blue, col.none=red, col.model=green, num.label="", merged=TRUE, no.ylab=TRUE)
mtext(label_name(3)[1], side=2, at=0.05+0.01, line=1.3+g, cex=cex.m, las=2, adj=0)
mtext("Net benefit", side=2, line=2.5, cex=cex.m)
plot_DCA(dat.col.2, cox.col.2, u=3, thresholds=seq(0, 1, 0.01), lwd=lwd, col.all=blue, col.none=red, col.model=green, num.label="", no.y.labels=TRUE, merged=TRUE, no.ylab=TRUE)
mtext(label_name(3)[2], side=2, at=0.05+0.01, line=1.3+g, cex=cex.m, las=2, adj=0)
plot_DCA(dat.col.1, cox.col.1, u=6, thresholds=seq(0, 1, 0.01), lwd=lwd, col.all=blue, col.none=red, col.model=green, num.label="", merged=TRUE, no.ylab=TRUE)
mtext(label_name(6)[1], side=2, at=0.1+0.015, line=1.3+g, cex=cex.m, las=2, adj=0)
mtext("Net benefit", side=2, line=2.5, cex=cex.m)
plot_DCA(dat.col.2, cox.col.2, u=6, thresholds=seq(0, 1, 0.01), lwd=lwd, col.all=blue, col.none=red, col.model=green, num.label="", no.y.labels=TRUE, merged=TRUE, no.ylab=TRUE)
mtext(label_name(6)[2], side=2, at=0.1+0.015, line=1.3+g, cex=cex.m, las=2, adj=0)
plot_DCA(dat.col.1, cox.col.1, u=12, thresholds=seq(0, 1, 0.01), lwd=lwd, col.all=blue, col.none=red, col.model=green, num.label="", no.ylab=TRUE)
mtext("Net benefit", side=2, line=2.5, cex=cex.m)
mtext(label_name(12)[1], side=2, at=0.15+0.02, line=1.3+g, cex=cex.m, las=2, adj=0)
plot_DCA(dat.col.2, cox.col.2, u=12, thresholds=seq(0, 1, 0.01), lwd=lwd, col.all=blue, col.none=red, col.model=green, num.label="", no.y.labels=TRUE, no.ylab=TRUE)
mtext(label_name(12)[2], side=2, at=0.15+0.02, line=1.3+g, cex=cex.m, las=2, adj=0)
dev.off()


############# merged plot, Control
file.name = paste0("DCA_总模型_Control队列.", file.type)
if( file.type=="pdf" ){
	pdf(file.name, width=W, height=H)
}else{
	tiff(file.name, width=W*L*2, height=H*L*3, pointsize=pt*1.5, res=res, compression="lzw")
}

par(mfrow=c(3, 2))
g=-0.00;cex.m=2/3
#### train
dat = dt.modeling.control.train
formula = as.formula(paste("Surv(time, event)~", paste(var.names, collapse="+")))
cox.col.1 = coxph(formula, data=dat, x=TRUE, y=TRUE)
dat.col.1 = dat
#### valid
dat = dt.modeling.control.valid
formula = as.formula(paste("Surv(time, event)~", paste(var.names, collapse="+")))
cox.col.2 = coxph(formula, data=dat, x=TRUE, y=TRUE)
dat.col.2 = dat

plot_DCA(dat.col.1, cox.col.1, u=3, thresholds=seq(0, 1, 0.01), lwd=lwd, col.all=blue, col.none=red, col.model=green, num.label="", merged=TRUE, no.ylab=TRUE)
mtext(label_name(3)[1], side=2, at=0.05+0.01, line=1.3+g, cex=cex.m, las=2, adj=0)
mtext("Net benefit", side=2, line=2.5, cex=cex.m)
plot_DCA(dat.col.2, cox.col.2, u=3, thresholds=seq(0, 1, 0.01), lwd=lwd, col.all=blue, col.none=red, col.model=green, num.label="", no.y.labels=TRUE, merged=TRUE, no.ylab=TRUE)
mtext(label_name(3)[2], side=2, at=0.05+0.01, line=1.3+g, cex=cex.m, las=2, adj=0)
plot_DCA(dat.col.1, cox.col.1, u=6, thresholds=seq(0, 1, 0.01), lwd=lwd, col.all=blue, col.none=red, col.model=green, num.label="", merged=TRUE, no.ylab=TRUE)
mtext(label_name(6)[1], side=2, at=0.1+0.015, line=1.3+g, cex=cex.m, las=2, adj=0)
mtext("Net benefit", side=2, line=2.5, cex=cex.m)
plot_DCA(dat.col.2, cox.col.2, u=6, thresholds=seq(0, 1, 0.01), lwd=lwd, col.all=blue, col.none=red, col.model=green, num.label="", no.y.labels=TRUE, merged=TRUE, no.ylab=TRUE)
mtext(label_name(6)[2], side=2, at=0.1+0.015, line=1.3+g, cex=cex.m, las=2, adj=0)
plot_DCA(dat.col.1, cox.col.1, u=12, thresholds=seq(0, 1, 0.01), lwd=lwd, col.all=blue, col.none=red, col.model=green, num.label="", no.ylab=TRUE)
mtext("Net benefit", side=2, line=2.5, cex=cex.m)
mtext(label_name(12)[1], side=2, at=0.15+0.02, line=1.3+g, cex=cex.m, las=2, adj=0)
plot_DCA(dat.col.2, cox.col.2, u=12, thresholds=seq(0, 1, 0.01), lwd=lwd, col.all=blue, col.none=red, col.model=green, num.label="", no.y.labels=TRUE, no.ylab=TRUE)
mtext(label_name(12)[2], side=2, at=0.15+0.02, line=1.3+g, cex=cex.m, las=2, adj=0)
dev.off()


############################################
# risk groups
############################################

library(survival)
library(survminer)

# set parameters
W=6.5;H=5+1;L=300;res=100;pt=40;lwd=5
file.type = "tif"
dig = 4

HR.tab = NULL
tab.names = NULL
var.names = c("age", "BMI", "hypertension", "LDLC", "diabetes", "smoke", "zeren", "xiazhai")

###### overall
#### train
file.name = paste0("KM-总模型_总队列_训练集.", file.type)
if( file.type=="pdf" ){
	pdf(file.name, width=W*1.5, height=H*1.5, onefile=FALSE)
}else{
	tiff(file.name, width=W*L, height=H*L, pointsize=pt, res=res, compression="lzw")
}
dat = dt.modeling.overall.train
dat = group_dat(dat, var.names, type="event")
fit = survfit(Surv(time, event)~group, data=dat)
surv_diff = survdiff(Surv(time, event)~group, data=dat)
pval = 1 - pchisq(surv_diff$chisq, length(surv_diff$n) - 1)
pval = ifelse(pval<0.0001, "P<0.0001", paste0("P=", pval))
cox = coxph(Surv(time, event)~as.factor(group), data=dat)
HR = cbind(t(sapply(cox_out(cox, dig=2)[,2], function(x)unlist(strsplit(x, split=" \\(|\\)")))), cox_out(cox, dig=3)[,1])
plot_KM(fit, num.label="A", lwd=lwd, col.low=green, col.mid=blue, col.high=red, pval=pval, HR=HR)
dev.off()
#### valid
file.name = paste0("KM-总模型_总队列_验证集.", file.type)
if( file.type=="pdf" ){
	pdf(file.name, width=W*1.5, height=H*1.5, onefile=FALSE)
}else{
	tiff(file.name, width=W*L, height=H*L, pointsize=pt, res=res, compression="lzw")
}
dat = dt.modeling.overall.valid
dat = group_dat(dat, var.names, type="event")
fit = survfit(Surv(time, event)~group, data=dat)
surv_diff = survdiff(Surv(time, event)~group, data=dat)
pval = 1 - pchisq(surv_diff$chisq, length(surv_diff$n) - 1)
pval = ifelse(pval<0.0001, "P<0.0001", paste0("P=", pval))
cox = coxph(Surv(time, event)~as.factor(group), data=dat)
HR = cbind(t(sapply(cox_out(cox, dig=2)[,2], function(x)unlist(strsplit(x, split=" \\(|\\)")))), cox_out(cox, dig=3)[,1])
plot_KM(fit, num.label="B", lwd=lwd, col.low=green, col.mid=blue, col.high=red, pval=pval, HR=HR)
dev.off()


###### treat
#### train
file.name = paste0("KM-总模型_RIC队列_训练集.", file.type)
if( file.type=="pdf" ){
	pdf(file.name, width=W*1.5, height=H*1.5, onefile=FALSE)
}else{
	tiff(file.name, width=W*L, height=H*L, pointsize=pt, res=res, compression="lzw")
}
dat = dt.modeling.treat.train
dat = group_dat(dat, var.names)
fit = survfit(Surv(time, event)~group, data=dat)
surv_diff = survdiff(Surv(time, event)~group, data=dat)
pval = 1 - pchisq(surv_diff$chisq, length(surv_diff$n) - 1)
pval = ifelse(pval<0.0001, "P<0.0001", paste0("P=", pval))
cox = coxph(Surv(time, event)~as.factor(group), data=dat)
HR = cbind(t(sapply(cox_out(cox, dig=2)[,2], function(x)unlist(strsplit(x, split=" \\(|\\)")))), cox_out(cox, dig=3)[,1])
plot_KM(fit, num.label="A", lwd=lwd, col.low=green, col.mid=blue, col.high=red, pval=pval, HR=HR)
dev.off()
#### valid
file.name = paste0("KM-总模型_RIC队列_验证集.", file.type)
if( file.type=="pdf" ){
	pdf(file.name, width=W*1.5, height=H*1.5, onefile=FALSE)
}else{
	tiff(file.name, width=W*L, height=H*L, pointsize=pt, res=res, compression="lzw")
}
dat = dt.modeling.treat.valid
dat = group_dat(dat, var.names)
fit = survfit(Surv(time, event)~group, data=dat)
surv_diff = survdiff(Surv(time, event)~group, data=dat)
pval = 1 - pchisq(surv_diff$chisq, length(surv_diff$n) - 1)
pval = ifelse(pval<0.0001, "P<0.0001", paste0("P=", pval))
cox = coxph(Surv(time, event)~as.factor(group), data=dat)
HR = cbind(t(sapply(cox_out(cox, dig=2)[,2], function(x)unlist(strsplit(x, split=" \\(|\\)")))), cox_out(cox, dig=3)[,1])
plot_KM(fit, num.label="B", lwd=lwd, col.low=green, col.mid=blue, col.high=red, pval=pval, HR=HR)
dev.off()


###### control
#### train
file.name = paste0("KM-总模型_Control队列_训练集.", file.type)
if( file.type=="pdf" ){
	pdf(file.name, width=W*1.5, height=H*1.5, onefile=FALSE)
}else{
	tiff(file.name, width=W*L, height=H*L, pointsize=pt, res=res, compression="lzw")
}
dat = dt.modeling.control.train
dat = group_dat(dat, var.names)
fit = survfit(Surv(time, event)~group, data=dat)
surv_diff = survdiff(Surv(time, event)~group, data=dat)
pval = 1 - pchisq(surv_diff$chisq, length(surv_diff$n) - 1)
pval = ifelse(pval<0.0001, "P<0.0001", paste0("P=", pval))
cox = coxph(Surv(time, event)~as.factor(group), data=dat)
HR = cbind(t(sapply(cox_out(cox, dig=2)[,2], function(x)unlist(strsplit(x, split=" \\(|\\)")))), cox_out(cox, dig=3)[,1])
plot_KM(fit, num.label="A", lwd=lwd, col.low=green, col.mid=blue, col.high=red, pval=pval, HR=HR)
dev.off()
#### valid
file.name = paste0("KM-总模型_Control队列_验证集.", file.type)
if( file.type=="pdf" ){
	pdf(file.name, width=W*1.5, height=H*1.5, onefile=FALSE)
}else{
	tiff(file.name, width=W*L, height=H*L, pointsize=pt, res=res, compression="lzw")
}
dat = dt.modeling.control.valid
dat = group_dat(dat, var.names)
fit = survfit(Surv(time, event)~group, data=dat)
surv_diff = survdiff(Surv(time, event)~group, data=dat)
pval = 1 - pchisq(surv_diff$chisq, length(surv_diff$n) - 1)
pval = ifelse(pval<0.0001, "P<0.0001", paste0("P=", pval))
cox = coxph(Surv(time, event)~as.factor(group), data=dat)
HR = cbind(t(sapply(cox_out(cox, dig=2)[,2], function(x)unlist(strsplit(x, split=" \\(|\\)")))), cox_out(cox, dig=3)[,1])
plot_KM(fit, num.label="B", lwd=lwd, col.low=green, col.mid=blue, col.high=red, pval=pval, HR=HR)
dev.off()


########################## 
## 1 column
########################## 
############# merged plot, overll
file.name = paste0("KM_总模型_总队列_1c.", file.type)
if( file.type=="pdf" ){
	pdf(file.name, width=W*1.5, height=H*1.5, onefile=FALSE)
}else{
	tiff(file.name, width=W*L, height=H*L*2, pointsize=pt, res=res, compression="lzw")
}
par(mfrow=c(2, 1))
#### train
dat = dt.modeling.overall.train
dat = group_dat(dat, var.names)
fit = survfit(Surv(time, event)~group, data=dat)
surv_diff = survdiff(Surv(time, event)~group, data=dat)
pval = 1 - pchisq(surv_diff$chisq, length(surv_diff$n) - 1)
pval = ifelse(pval<0.0001, "P<0.0001", paste0("P=", pval))
cox = coxph(Surv(time, event)~as.factor(group), data=dat)
HR = cbind(t(sapply(cox_out(cox, dig=2)[,2], function(x)unlist(strsplit(x, split=" \\(|\\)")))), cox_out(cox, dig=3)[,1])
plot_KM(fit, num.label="A", lwd=lwd, col.low=green, col.mid=blue, col.high=red, pval=pval, HR=HR, no.xlab=TRUE)
#### valid
dat = dt.modeling.overall.valid
dat = group_dat(dat, var.names)
fit = survfit(Surv(time, event)~group, data=dat)
surv_diff = survdiff(Surv(time, event)~group, data=dat)
pval = 1 - pchisq(surv_diff$chisq, length(surv_diff$n) - 1)
pval = ifelse(pval<0.0001, "P<0.0001", paste0("P=", pval))
cox = coxph(Surv(time, event)~as.factor(group), data=dat)
HR = cbind(t(sapply(cox_out(cox, dig=2)[,2], function(x)unlist(strsplit(x, split=" \\(|\\)")))), cox_out(cox, dig=3)[,1])
plot_KM(fit, num.label="B", lwd=lwd, col.low=green, col.mid=blue, col.high=red, pval=pval, HR=HR)
dev.off()


############# merged plot, RIC
file.name = paste0("KM_总模型_RIC队列_1c.", file.type)
if( file.type=="pdf" ){
	pdf(file.name, width=W*1.5, height=H*1.5, onefile=FALSE)
}else{
	tiff(file.name, width=W*L, height=H*L*2, pointsize=pt, res=res, compression="lzw")
}
par(mfrow=c(2, 1))
#### train
dat = dt.modeling.treat.train
dat = group_dat(dat, var.names)
fit = survfit(Surv(time, event)~group, data=dat)
surv_diff = survdiff(Surv(time, event)~group, data=dat)
pval = 1 - pchisq(surv_diff$chisq, length(surv_diff$n) - 1)
pval = ifelse(pval<0.0001, "P<0.0001", paste0("P=", pval))
cox = coxph(Surv(time, event)~as.factor(group), data=dat)
HR = cbind(t(sapply(cox_out(cox, dig=2)[,2], function(x)unlist(strsplit(x, split=" \\(|\\)")))), cox_out(cox, dig=3)[,1])
plot_KM(fit, num.label="A", lwd=lwd, col.low=green, col.mid=blue, col.high=red, pval=pval, HR=HR, no.xlab=TRUE)
#### valid
dat = dt.modeling.treat.valid
dat = group_dat(dat, var.names)
fit = survfit(Surv(time, event)~group, data=dat)
surv_diff = survdiff(Surv(time, event)~group, data=dat)
pval = 1 - pchisq(surv_diff$chisq, length(surv_diff$n) - 1)
pval = ifelse(pval<0.0001, "P<0.0001", paste0("P=", pval))
cox = coxph(Surv(time, event)~as.factor(group), data=dat)
HR = cbind(t(sapply(cox_out(cox, dig=2)[,2], function(x)unlist(strsplit(x, split=" \\(|\\)")))), cox_out(cox, dig=3)[,1])
plot_KM(fit, num.label="B", lwd=lwd, col.low=green, col.mid=blue, col.high=red, pval=pval, HR=HR)
dev.off()


############# merged plot, Control
file.name = paste0("KM_总模型_Control队列_1c.", file.type)
if( file.type=="pdf" ){
	pdf(file.name, width=W*1.5, height=H*1.5, onefile=FALSE)
}else{
	tiff(file.name, width=W*L, height=H*L*2, pointsize=pt, res=res, compression="lzw")
}
par(mfrow=c(2, 1))
#### train
dat = dt.modeling.control.train
dat = group_dat(dat, var.names)
fit = survfit(Surv(time, event)~group, data=dat)
surv_diff = survdiff(Surv(time, event)~group, data=dat)
pval = 1 - pchisq(surv_diff$chisq, length(surv_diff$n) - 1)
pval = ifelse(pval<0.0001, "P<0.0001", paste0("P=", pval))
cox = coxph(Surv(time, event)~as.factor(group), data=dat)
HR = cbind(t(sapply(cox_out(cox, dig=2)[,2], function(x)unlist(strsplit(x, split=" \\(|\\)")))), cox_out(cox, dig=3)[,1])
plot_KM(fit, num.label="A", lwd=lwd, col.low=green, col.mid=blue, col.high=red, pval=pval, HR=HR, no.xlab=TRUE)
#### valid
dat = dt.modeling.control.valid
dat = group_dat(dat, var.names)
fit = survfit(Surv(time, event)~group, data=dat)
surv_diff = survdiff(Surv(time, event)~group, data=dat)
pval = 1 - pchisq(surv_diff$chisq, length(surv_diff$n) - 1)
pval = ifelse(pval<0.0001, "P<0.0001", paste0("P=", pval))
cox = coxph(Surv(time, event)~as.factor(group), data=dat)
HR = cbind(t(sapply(cox_out(cox, dig=2)[,2], function(x)unlist(strsplit(x, split=" \\(|\\)")))), cox_out(cox, dig=3)[,1])
plot_KM(fit, num.label="B", lwd=lwd, col.low=green, col.mid=blue, col.high=red, pval=pval, HR=HR)
dev.off()


########################## 
## 2 column
########################## 
############# merged plot, overll
file.name = paste0("KM_总模型_总队列_2c.", file.type)
if( file.type=="pdf" ){
	pdf(file.name, width=W*1.5, height=H*1.5, onefile=FALSE)
}else{
	tiff(file.name, width=W*L*2, height=H*L, pointsize=pt, res=res, compression="lzw")
}
par(mfrow=c(1, 2))
#### train
dat = dt.modeling.overall.train
dat = group_dat(dat, var.names)
fit = survfit(Surv(time, event)~group, data=dat)
surv_diff = survdiff(Surv(time, event)~group, data=dat)
pval = 1 - pchisq(surv_diff$chisq, length(surv_diff$n) - 1)
pval = ifelse(pval<0.0001, "P<0.0001", paste0("P=", pval))
cox = coxph(Surv(time, event)~as.factor(group), data=dat)
HR = cbind(t(sapply(cox_out(cox, dig=2)[,2], function(x)unlist(strsplit(x, split=" \\(|\\)")))), cox_out(cox, dig=3)[,1])
plot_KM(fit, num.label="A", lwd=lwd, col.low=green, col.mid=blue, col.high=red, pval=pval, HR=HR)
#### valid
dat = dt.modeling.overall.valid
dat = group_dat(dat, var.names)
fit = survfit(Surv(time, event)~group, data=dat)
surv_diff = survdiff(Surv(time, event)~group, data=dat)
pval = 1 - pchisq(surv_diff$chisq, length(surv_diff$n) - 1)
pval = ifelse(pval<0.0001, "P<0.0001", paste0("P=", pval))
cox = coxph(Surv(time, event)~as.factor(group), data=dat)
HR = cbind(t(sapply(cox_out(cox, dig=2)[,2], function(x)unlist(strsplit(x, split=" \\(|\\)")))), cox_out(cox, dig=3)[,1])
plot_KM(fit, num.label="B", lwd=lwd, col.low=green, col.mid=blue, col.high=red, pval=pval, HR=HR, no.ylab=TRUE)
dev.off()


############# merged plot, RIC
file.name = paste0("KM_总模型_RIC队列_2c.", file.type)
if( file.type=="pdf" ){
	pdf(file.name, width=W*1.5, height=H*1.5, onefile=FALSE)
}else{
	tiff(file.name, width=W*L*2, height=H*L, pointsize=pt, res=res, compression="lzw")
}
par(mfrow=c(1, 2))
#### train
dat = dt.modeling.treat.train
dat = group_dat(dat, var.names)
fit = survfit(Surv(time, event)~group, data=dat)
surv_diff = survdiff(Surv(time, event)~group, data=dat)
pval = 1 - pchisq(surv_diff$chisq, length(surv_diff$n) - 1)
pval = ifelse(pval<0.0001, "P<0.0001", paste0("P=", pval))
cox = coxph(Surv(time, event)~as.factor(group), data=dat)
HR = cbind(t(sapply(cox_out(cox, dig=2)[,2], function(x)unlist(strsplit(x, split=" \\(|\\)")))), cox_out(cox, dig=3)[,1])
plot_KM(fit, num.label="A", lwd=lwd, col.low=green, col.mid=blue, col.high=red, pval=pval, HR=HR)
#### valid
dat = dt.modeling.treat.valid
dat = group_dat(dat, var.names)
fit = survfit(Surv(time, event)~group, data=dat)
surv_diff = survdiff(Surv(time, event)~group, data=dat)
pval = 1 - pchisq(surv_diff$chisq, length(surv_diff$n) - 1)
pval = ifelse(pval<0.0001, "P<0.0001", paste0("P=", pval))
cox = coxph(Surv(time, event)~as.factor(group), data=dat)
HR = cbind(t(sapply(cox_out(cox, dig=2)[,2], function(x)unlist(strsplit(x, split=" \\(|\\)")))), cox_out(cox, dig=3)[,1])
plot_KM(fit, num.label="B", lwd=lwd, col.low=green, col.mid=blue, col.high=red, pval=pval, HR=HR, no.ylab=TRUE)
dev.off()


############# merged plot, Control
file.name = paste0("KM_总模型_Control队列_2c.", file.type)
if( file.type=="pdf" ){
	pdf(file.name, width=W*1.5, height=H*1.5, onefile=FALSE)
}else{
	tiff(file.name, width=W*L*2, height=H*L, pointsize=pt, res=res, compression="lzw")
}
par(mfrow=c(1, 2))
#### train
dat = dt.modeling.control.train
dat = group_dat(dat, var.names)
fit = survfit(Surv(time, event)~group, data=dat)
surv_diff = survdiff(Surv(time, event)~group, data=dat)
pval = 1 - pchisq(surv_diff$chisq, length(surv_diff$n) - 1)
pval = ifelse(pval<0.0001, "P<0.0001", paste0("P=", pval))
cox = coxph(Surv(time, event)~as.factor(group), data=dat)
HR = cbind(t(sapply(cox_out(cox, dig=2)[,2], function(x)unlist(strsplit(x, split=" \\(|\\)")))), cox_out(cox, dig=3)[,1])
plot_KM(fit, num.label="A", lwd=lwd, col.low=green, col.mid=blue, col.high=red, pval=pval, HR=HR)
#### valid
dat = dt.modeling.control.valid
dat = group_dat(dat, var.names)
fit = survfit(Surv(time, event)~group, data=dat)
surv_diff = survdiff(Surv(time, event)~group, data=dat)
pval = 1 - pchisq(surv_diff$chisq, length(surv_diff$n) - 1)
pval = ifelse(pval<0.0001, "P<0.0001", paste0("P=", pval))
cox = coxph(Surv(time, event)~as.factor(group), data=dat)
HR = cbind(t(sapply(cox_out(cox, dig=2)[,2], function(x)unlist(strsplit(x, split=" \\(|\\)")))), cox_out(cox, dig=3)[,1])
plot_KM(fit, num.label="B", lwd=lwd, col.low=green, col.mid=blue, col.high=red, pval=pval, HR=HR, no.ylab=TRUE)
dev.off()
