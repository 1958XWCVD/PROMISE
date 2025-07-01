

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
# VIF analysis
############################################
dat = dt.modeling.overall.train
var.names = c("age", "BMI", "hypertension", "LDLC", "diabetes", "smoke", "zeren", "xiazhai", 
			"gaozhi", "SP", "TIA", "AM")
formula = as.formula(paste("Surv(time, event)~", paste(var.names, collapse="+")))
cox = coxph(formula, data=dat)

library(car)
vif_values <- vif(cox)
print(vif_values)
write.xlsx(vif_values, file="results.xlsx", sheetName="VIF", append=TRUE)


############################################
# Sensitivity Analysis on the CTA subgroup
############################################

idx = which(dt.modeling.overall.train$CTA==1)
dat.CTA.train = dt.modeling.overall.train[idx,]

idx = which(dt.modeling.overall.valid$CTA==1)
dat.CTA.valid = dt.modeling.overall.valid[idx,]

# C-index for the CTA subgroup
B = 1000; dig = 2; dig.P = 5
vars.all = c("age", "BMI", "hypertension", "LDLC", "diabetes", "smoke", "zeren", "xiazhai")
set.seed(111)	#132-B100
C.index.CTA.train = C_index_CI(dat.CTA.train, vars.all, B=B, dig=dig, dig.P=dig.P)
C.index.CTA.valid = C_index_CI(dat.CTA.valid, vars.all, B=B, dig=dig, dig.P=dig.P)
C.CTA.tabs = cbind(C.index.CTA.train, C.index.CTA.valid)
rownames(C.CTA.tabs) = c("Overall", vars.all)
colnames(C.CTA.tabs) = c("C-index-train", "P-val-train", "C-index-valid", "P-val-valid")
write.xlsx(C.CTA.tabs, file="results.xlsx", sheetName="C_CTA", append=TRUE)

# The modified H-L test for the CTA subgroup
var.names = c("age", "BMI", "hypertension", "LDLC", "diabetes", "smoke", "zeren", "xiazhai")
formula = as.formula(paste("Surv(time, event)~", paste(var.names, collapse="+")))
cox = coxph(formula, data=dat.CTA.train, x=TRUE, y=TRUE)
survMisc::gof(cox)$lrTest


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
u=12;B=10;H=5;W=6.5;L=300;res=100;pt=40;lwd=5
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
if( file.type=="tif" ){
	tiff(file.name, width=W*L, height=H*L, pointsize=pt, res=res, compression="lzw")
}else if( file.type=="eps" ){
	pt=20;lwd=4
	postscript(file.name, width=W*1.5, height=H*1.5, pointsize=pt, horizontal=FALSE, onefile=FALSE, paper="special")
}
dd=datadist(dt.modeling.overall.train);option=options(datadist="dd")
plot_cal(dt.modeling.overall.train, dt.modeling.overall.valid, vars, u=u, lwd=lwd, col.train=blue, col.val=red, col.ideal=green, num.label=num.label)
dev.off()

# treat
file.name = paste0("calibration_", num.label, "_总模型_RIC队列_month-", u, ".", file.type)
if( file.type=="tif" ){
	tiff(file.name, width=W*L, height=H*L, pointsize=pt, res=res, compression="lzw")
}else if( file.type=="eps" ){
	pt=20;lwd=4
	postscript(file.name, width=W*1.5, height=H*1.5, pointsize=pt, horizontal=FALSE, onefile=FALSE, paper="special")
}
dd=datadist(dt.modeling.treat.train);option=options(datadist="dd")
plot_cal(dt.modeling.treat.train, dt.modeling.treat.valid, vars, u=u, lwd=lwd, col.train=blue, col.val=red, col.ideal=green, num.label=num.label)
dev.off()

# control
file.name = paste0("calibration_", num.label, "_总模型_Control队列_month-", u, ".", file.type)
if( file.type=="tif" ){
	tiff(file.name, width=W*L, height=H*L, pointsize=pt, res=res, compression="lzw")
}else if( file.type=="eps" ){
	pt=20;lwd=4
	postscript(file.name, width=W*1.5, height=H*1.5, pointsize=pt, horizontal=FALSE, onefile=FALSE, paper="special")
}
dd=datadist(dt.modeling.control.train);option=options(datadist="dd")
plot_cal(dt.modeling.control.train, dt.modeling.control.valid, vars, u=u, lwd=lwd, col.train=blue, col.val=red, col.ideal=green, num.label=num.label)
dev.off()


## merged plot
file.name = paste0("calibration_总模型_总队列.", file.type)
if( file.type=="tif" ){
	tiff(file.name, width=W*L*3, height=H*L, pointsize=pt, res=res, compression="lzw")
}else if( file.type=="eps" ){
	pt=20;lwd=4
	postscript(file.name, width=W*1.5*3, height=H*1.5, pointsize=pt, horizontal=FALSE, onefile=FALSE, paper="special")
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
if( file.type=="tif" ){
	tiff(file.name, width=W*L*3, height=H*L, pointsize=pt, res=res, compression="lzw")
}else if( file.type=="eps" ){
	pt=20;lwd=4
	postscript(file.name, width=W*1.5*3, height=H*1.5, pointsize=pt, horizontal=FALSE, onefile=FALSE, paper="special")
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
if( file.type=="tif" ){
	tiff(file.name, width=W*L*3, height=H*L, pointsize=pt, res=res, compression="lzw")
}else if( file.type=="eps" ){
	pt=20;lwd=4
	postscript(file.name, width=W*1.5*3, height=H*1.5, pointsize=pt, horizontal=FALSE, onefile=FALSE, paper="special")
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
if( file.type=="tif" ){
	tiff(file.name, width=W*L, height=H*L, pointsize=pt, res=res, compression="lzw")
}else if( file.type=="eps" ){
	pt=20;lwd=4
	postscript(file.name, width=W*1.5, height=H*1.5, pointsize=pt, horizontal=FALSE, onefile=FALSE, paper="special")
}
dat = dt.modeling.overall.train
formula = as.formula(paste("Surv(time, event)~", paste(var.names, collapse="+")))
cox = coxph(formula, data=dat, x=TRUE, y=TRUE) 
plot_DCA(dat, cox, u=u, thresholds=seq(0, 1, 0.01), lwd=lwd, col.all=blue, col.none=red, col.model=green, num.label=label_name(u)[1])
dev.off()

# overall valid
file.name = paste0("DCA_总模型_总队列-", label_name(u)[2], ".", file.type)
if( file.type=="tif" ){
	tiff(file.name, width=W*L, height=H*L, pointsize=pt, res=res, compression="lzw")
}else if( file.type=="eps" ){
	pt=20;lwd=4
	postscript(file.name, width=W*1.5, height=H*1.5, pointsize=pt, horizontal=FALSE, onefile=FALSE, paper="special")
}
dat = dt.modeling.overall.valid
formula = as.formula(paste("Surv(time, event)~", paste(var.names, collapse="+")))
cox = coxph(formula, data=dat, x=TRUE, y=TRUE)
plot_DCA(dat, cox, u=u, thresholds=seq(0, 1, 0.01), lwd=lwd, col.all=blue, col.none=red, col.model=green, num.label=label_name(u)[2])
dev.off()

# treat train
file.name = paste0("DCA_总模型_RIC队列-", label_name(u)[1], ".", file.type)
if( file.type=="tif" ){
	tiff(file.name, width=W*L, height=H*L, pointsize=pt, res=res, compression="lzw")
}else if( file.type=="eps" ){
	pt=20;lwd=4
	postscript(file.name, width=W*1.5, height=H*1.5, pointsize=pt, horizontal=FALSE, onefile=FALSE, paper="special")
}
dat = dt.modeling.treat.train
formula = as.formula(paste("Surv(time, event)~", paste(var.names, collapse="+")))
cox = coxph(formula, data=dat, x=TRUE, y=TRUE) 
plot_DCA(dat, cox, u=u, thresholds=seq(0, 1, 0.01), lwd=lwd, col.all=blue, col.none=red, col.model=green, num.label=label_name(u)[1])
dev.off()

# treat valid
file.name = paste0("DCA_总模型_RIC队列-", label_name(u)[2], ".", file.type)
if( file.type=="tif" ){
	tiff(file.name, width=W*L, height=H*L, pointsize=pt, res=res, compression="lzw")
}else if( file.type=="eps" ){
	pt=20;lwd=4
	postscript(file.name, width=W*1.5, height=H*1.5, pointsize=pt, horizontal=FALSE, onefile=FALSE, paper="special")
}
dat = dt.modeling.treat.valid
formula = as.formula(paste("Surv(time, event)~", paste(var.names, collapse="+")))
cox = coxph(formula, data=dat, x=TRUE, y=TRUE)
plot_DCA(dat, cox, u=u, thresholds=seq(0, 1, 0.01), lwd=lwd, col.all=blue, col.none=red, col.model=green, num.label=label_name(u)[2])
dev.off()


# control train
file.name = paste0("DCA_总模型_Control队列-", label_name(u)[1], ".", file.type)
if( file.type=="tif" ){
	tiff(file.name, width=W*L, height=H*L, pointsize=pt, res=res, compression="lzw")
}else if( file.type=="eps" ){
	pt=20;lwd=4
	postscript(file.name, width=W*1.5, height=H*1.5, pointsize=pt, horizontal=FALSE, onefile=FALSE, paper="special")
}
dat = dt.modeling.control.train
formula = as.formula(paste("Surv(time, event)~", paste(var.names, collapse="+")))
cox = coxph(formula, data=dat, x=TRUE, y=TRUE) 
plot_DCA(dat, cox, u=u, thresholds=seq(0, 1, 0.01), lwd=lwd, col.all=blue, col.none=red, col.model=green, num.label=label_name(u)[1])
dev.off()

# control valid
file.name = paste0("DCA_总模型_Control队列-", label_name(u)[2], ".", file.type)
if( file.type=="tif" ){
	tiff(file.name, width=W*L, height=H*L, pointsize=pt, res=res, compression="lzw")
}else if( file.type=="eps" ){
	pt=20;lwd=4
	postscript(file.name, width=W*1.5, height=H*1.5, pointsize=pt, horizontal=FALSE, onefile=FALSE, paper="special")
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
if( file.type=="tif" ){
	tiff(file.name, width=W*L*3, height=H*L*2, pointsize=pt, res=res, compression="lzw")
}else if( file.type=="eps" ){
	pt=20;lwd=4
	postscript(file.name, width=W*1.5*3, height=H*1.5*2, pointsize=pt, horizontal=FALSE, onefile=FALSE, paper="special")
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
if( file.type=="tif" ){
	tiff(file.name, width=W*L*3, height=H*L*2, pointsize=pt, res=res, compression="lzw")
}else if( file.type=="eps" ){
	pt=20;lwd=4
	postscript(file.name, width=W*1.5*3, height=H*1.5*2, pointsize=pt, horizontal=FALSE, onefile=FALSE, paper="special")
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
if( file.type=="tif" ){
	tiff(file.name, width=W*L*3, height=H*L*2, pointsize=pt, res=res, compression="lzw")
}else if( file.type=="eps" ){
	pt=20;lwd=4
	postscript(file.name, width=W*1.5*3, height=H*1.5*2, pointsize=pt, horizontal=FALSE, onefile=FALSE, paper="special")
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
if( file.type=="tif" ){
	tiff(file.name, width=W*L, height=H*L, pointsize=pt, res=res, compression="lzw")
}else if( file.type=="eps" ){
	pt=20;lwd=4
	postscript(file.name, width=W*1.5, height=H*1.5, pointsize=pt, horizontal=FALSE, onefile=FALSE, paper="special")
}
dat = dt.modeling.overall.train
dat = group_dat(dat, var.names, type="event")
table(dat$group)
cutoff.1 = (dat$risk_score[tail(which(dat$group==1), 1)] + dat$risk_score[which(dat$group==2)[1]]) / 2
cutoff.2 = (dat$risk_score[tail(which(dat$group==2), 1)] + dat$risk_score[which(dat$group==3)[1]]) / 2
sum(dat$event[dat$group==1]);round(sum(dat$event[dat$group==1]) / sum(dat$group==1) * 100, 2)
sum(dat$event[dat$group==2]);round(sum(dat$event[dat$group==2]) / sum(dat$group==2) * 100, 2)
sum(dat$event[dat$group==3]);round(sum(dat$event[dat$group==3]) / sum(dat$group==3) * 100, 2)
fit = survfit(Surv(time, event)~group, data=dat)
surv_diff = survdiff(Surv(time, event)~group, data=dat)
pval = 1 - pchisq(surv_diff$chisq, length(surv_diff$n) - 1)
pval = ifelse(pval<0.0001, "P<0.0001", paste0("P=", pval))
cox = coxph(Surv(time, event)~as.factor(group), data=dat)
HR = cbind(t(sapply(cox_out(cox, dig=2)[,2], function(x)unlist(strsplit(x, split=" \\(|\\)")))), cox_out(cox, dig=3)[,1])
HR
plot_KM(fit, num.label="A", lwd=lwd, col.low=green, col.mid=blue, col.high=red, pval=pval, HR=HR)
dev.off()
#### valid
file.name = paste0("KM-总模型_总队列_验证集.", file.type)
if( file.type=="tif" ){
	tiff(file.name, width=W*L, height=H*L, pointsize=pt, res=res, compression="lzw")
}else if( file.type=="eps" ){
	pt=20;lwd=4
	postscript(file.name, width=W*1.5, height=H*1.5, pointsize=pt, horizontal=FALSE, onefile=FALSE, paper="special")
}
dat.train = dt.modeling.overall.train
formula = as.formula(paste("Surv(time, event)~", paste(var.names, collapse="+")))
cox = coxph(formula, data=dat.train, x=TRUE, y=TRUE) 
dat = dt.modeling.overall.valid
risk.scores = predict(cox, newdata=dat, type="risk")
group.label = cut(risk.scores, breaks=c(-Inf, cutoff.1, cutoff.2, Inf), labels=1:3)
dat$group = as.numeric(group.label)
sum(dat$event[group.label==1]);round(sum(dat$event[group.label==1]) / sum(group.label==1) * 100, 2)
sum(dat$event[group.label==2]);round(sum(dat$event[group.label==2]) / sum(group.label==2) * 100, 2)
sum(dat$event[group.label==3]);round(sum(dat$event[group.label==3]) / sum(group.label==3) * 100, 2)
fit = survfit(Surv(time, event)~group, data=dat)
surv_diff = survdiff(Surv(time, event)~group, data=dat)
pval = 1 - pchisq(surv_diff$chisq, length(surv_diff$n) - 1)
pval = ifelse(pval<0.0001, "P<0.0001", paste0("P=", pval))
cox = coxph(Surv(time, event)~as.factor(group), data=dat)
HR = cbind(t(sapply(cox_out(cox, dig=2)[,2], function(x)unlist(strsplit(x, split=" \\(|\\)")))), cox_out(cox, dig=3)[,1])
HR
plot_KM(fit, num.label="B", lwd=lwd, col.low=green, col.mid=blue, col.high=red, pval=pval, HR=HR)
dev.off()


###### treat
#### train
file.name = paste0("KM-总模型_RIC队列_训练集.", file.type)
if( file.type=="tif" ){
	tiff(file.name, width=W*L, height=H*L, pointsize=pt, res=res, compression="lzw")
}else if( file.type=="eps" ){
	pt=20;lwd=4
	postscript(file.name, width=W*1.5, height=H*1.5, pointsize=pt, horizontal=FALSE, onefile=FALSE, paper="special")
}
dat.train = dt.modeling.overall.train
formula = as.formula(paste("Surv(time, event)~", paste(var.names, collapse="+")))
cox = coxph(formula, data=dat.train, x=TRUE, y=TRUE) 
dat = dt.modeling.treat.train
risk.scores = predict(cox, newdata=dat, type="risk")
group.label = cut(risk.scores, breaks=c(-Inf, cutoff.1, cutoff.2, Inf), labels=1:3)
dat$group = as.numeric(group.label)
sum(dat$event[group.label==1]);round(sum(dat$event[group.label==1]) / sum(group.label==1) * 100, 2)
sum(dat$event[group.label==2]);round(sum(dat$event[group.label==2]) / sum(group.label==2) * 100, 2)
sum(dat$event[group.label==3]);round(sum(dat$event[group.label==3]) / sum(group.label==3) * 100, 2)
fit = survfit(Surv(time, event)~group, data=dat)
surv_diff = survdiff(Surv(time, event)~group, data=dat)
pval = 1 - pchisq(surv_diff$chisq, length(surv_diff$n) - 1)
pval = ifelse(pval<0.0001, "P<0.0001", paste0("P=", pval))
cox = coxph(Surv(time, event)~as.factor(group), data=dat)
HR = cbind(t(sapply(cox_out(cox, dig=2)[,2], function(x)unlist(strsplit(x, split=" \\(|\\)")))), cox_out(cox, dig=3)[,1])
HR
plot_KM(fit, num.label="A", lwd=lwd, col.low=green, col.mid=blue, col.high=red, pval=pval, HR=HR)
dev.off()
#### valid
file.name = paste0("KM-总模型_RIC队列_验证集.", file.type)
if( file.type=="tif" ){
	tiff(file.name, width=W*L, height=H*L, pointsize=pt, res=res, compression="lzw")
}else if( file.type=="eps" ){
	pt=20;lwd=4
	postscript(file.name, width=W*1.5, height=H*1.5, pointsize=pt, horizontal=FALSE, onefile=FALSE, paper="special")
}
dat.train = dt.modeling.overall.train
formula = as.formula(paste("Surv(time, event)~", paste(var.names, collapse="+")))
cox = coxph(formula, data=dat.train, x=TRUE, y=TRUE) 
dat = dt.modeling.treat.valid
risk.scores = predict(cox, newdata=dat, type="risk")
group.label = cut(risk.scores, breaks=c(-Inf, cutoff.1, cutoff.2, Inf), labels=1:3)
dat$group = as.numeric(group.label)
sum(dat$event[group.label==1]);round(sum(dat$event[group.label==1]) / sum(group.label==1) * 100, 2)
sum(dat$event[group.label==2]);round(sum(dat$event[group.label==2]) / sum(group.label==2) * 100, 2)
sum(dat$event[group.label==3]);round(sum(dat$event[group.label==3]) / sum(group.label==3) * 100, 2)
fit = survfit(Surv(time, event)~group, data=dat)
surv_diff = survdiff(Surv(time, event)~group, data=dat)
pval = 1 - pchisq(surv_diff$chisq, length(surv_diff$n) - 1)
pval = ifelse(pval<0.0001, "P<0.0001", paste0("P=", pval))
cox = coxph(Surv(time, event)~as.factor(group), data=dat)
HR = cbind(t(sapply(cox_out(cox, dig=2)[,2], function(x)unlist(strsplit(x, split=" \\(|\\)")))), cox_out(cox, dig=3)[,1])
HR
plot_KM(fit, num.label="B", lwd=lwd, col.low=green, col.mid=blue, col.high=red, pval=pval, HR=HR)
dev.off()

###### control
#### train
file.name = paste0("KM-总模型_Control队列_训练集.", file.type)
if( file.type=="tif" ){
	tiff(file.name, width=W*L, height=H*L, pointsize=pt, res=res, compression="lzw")
}else if( file.type=="eps" ){
	pt=20;lwd=4
	postscript(file.name, width=W*1.5, height=H*1.5, pointsize=pt, horizontal=FALSE, onefile=FALSE, paper="special")
}
dat.train = dt.modeling.overall.train
formula = as.formula(paste("Surv(time, event)~", paste(var.names, collapse="+")))
cox = coxph(formula, data=dat.train, x=TRUE, y=TRUE) 
dat = dt.modeling.control.train
risk.scores = predict(cox, newdata=dat, type="risk")
group.label = cut(risk.scores, breaks=c(-Inf, cutoff.1, cutoff.2, Inf), labels=1:3)
dat$group = as.numeric(group.label)
sum(dat$event[group.label==1]);round(sum(dat$event[group.label==1]) / sum(group.label==1) * 100, 2)
sum(dat$event[group.label==2]);round(sum(dat$event[group.label==2]) / sum(group.label==2) * 100, 2)
sum(dat$event[group.label==3]);round(sum(dat$event[group.label==3]) / sum(group.label==3) * 100, 2)
fit = survfit(Surv(time, event)~group, data=dat)
surv_diff = survdiff(Surv(time, event)~group, data=dat)
pval = 1 - pchisq(surv_diff$chisq, length(surv_diff$n) - 1)
pval = ifelse(pval<0.0001, "P<0.0001", paste0("P=", pval))
cox = coxph(Surv(time, event)~as.factor(group), data=dat)
HR = cbind(t(sapply(cox_out(cox, dig=2)[,2], function(x)unlist(strsplit(x, split=" \\(|\\)")))), cox_out(cox, dig=3)[,1])
HR
plot_KM(fit, num.label="A", lwd=lwd, col.low=green, col.mid=blue, col.high=red, pval=pval, HR=HR)
dev.off()
#### valid
file.name = paste0("KM-总模型_Control队列_验证集.", file.type)
if( file.type=="tif" ){
	tiff(file.name, width=W*L, height=H*L, pointsize=pt, res=res, compression="lzw")
}else if( file.type=="eps" ){
	pt=20;lwd=4
	postscript(file.name, width=W*1.5, height=H*1.5, pointsize=pt, horizontal=FALSE, onefile=FALSE, paper="special")
}
dat.train = dt.modeling.overall.train
formula = as.formula(paste("Surv(time, event)~", paste(var.names, collapse="+")))
cox = coxph(formula, data=dat.train, x=TRUE, y=TRUE) 
dat = dt.modeling.control.valid
risk.scores = predict(cox, newdata=dat, type="risk")
group.label = cut(risk.scores, breaks=c(-Inf, cutoff.1, cutoff.2, Inf), labels=1:3)
dat$group = as.numeric(group.label)
sum(dat$event[group.label==1]);round(sum(dat$event[group.label==1]) / sum(group.label==1) * 100, 2)
sum(dat$event[group.label==2]);round(sum(dat$event[group.label==2]) / sum(group.label==2) * 100, 2)
sum(dat$event[group.label==3]);round(sum(dat$event[group.label==3]) / sum(group.label==3) * 100, 2)
fit = survfit(Surv(time, event)~group, data=dat)
surv_diff = survdiff(Surv(time, event)~group, data=dat)
pval = 1 - pchisq(surv_diff$chisq, length(surv_diff$n) - 1)
pval = ifelse(pval<0.0001, "P<0.0001", paste0("P=", pval))
cox = coxph(Surv(time, event)~as.factor(group), data=dat)
HR = cbind(t(sapply(cox_out(cox, dig=2)[,2], function(x)unlist(strsplit(x, split=" \\(|\\)")))), cox_out(cox, dig=3)[,1])
HR
plot_KM(fit, num.label="B", lwd=lwd, col.low=green, col.mid=blue, col.high=red, pval=pval, HR=HR)
dev.off()


########################## 
## 1 column
########################## 
############# merged plot, overll
file.name = paste0("KM_总模型_总队列_1c.", file.type)
if( file.type=="tif" ){
	tiff(file.name, width=W*L, height=H*L*2, pointsize=pt, res=res, compression="lzw")
}else if( file.type=="eps" ){
	pt=20;lwd=4
	postscript(file.name, width=W*1.5, height=H*1.5*2, pointsize=pt, horizontal=FALSE, onefile=FALSE, paper="special")
}
par(mfrow=c(2, 1))
#### train
dat = dt.modeling.overall.train
dat = group_dat(dat, var.names, type="event")
cutoff.1 = (dat$risk_score[tail(which(dat$group==1), 1)] + dat$risk_score[which(dat$group==2)[1]]) / 2
cutoff.2 = (dat$risk_score[tail(which(dat$group==2), 1)] + dat$risk_score[which(dat$group==3)[1]]) / 2
fit = survfit(Surv(time, event)~group, data=dat)
surv_diff = survdiff(Surv(time, event)~group, data=dat)
pval = 1 - pchisq(surv_diff$chisq, length(surv_diff$n) - 1)
pval = ifelse(pval<0.0001, "P<0.0001", paste0("P=", pval))
cox = coxph(Surv(time, event)~as.factor(group), data=dat)
HR = cbind(t(sapply(cox_out(cox, dig=2)[,2], function(x)unlist(strsplit(x, split=" \\(|\\)")))), cox_out(cox, dig=3)[,1])
plot_KM(fit, num.label="A", lwd=lwd, col.low=green, col.mid=blue, col.high=red, pval=pval, HR=HR)
#### valid
dat.train = dt.modeling.overall.train
formula = as.formula(paste("Surv(time, event)~", paste(var.names, collapse="+")))
cox = coxph(formula, data=dat.train, x=TRUE, y=TRUE) 
dat = dt.modeling.overall.valid
risk.scores = predict(cox, newdata=dat, type="risk")
group.label = cut(risk.scores, breaks=c(-Inf, cutoff.1, cutoff.2, Inf), labels=1:3)
dat$group = as.numeric(group.label)
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
############# merged plot, RIC
file.name = paste0("KM_总模型_RIC队列_2c.", file.type)
if( file.type=="tif" ){
	tiff(file.name, width=W*L*2, height=H*L, pointsize=pt, res=res, compression="lzw")
}else if( file.type=="eps" ){
	pt=20;lwd=4
	postscript(file.name, width=W*1.5*2, height=H*1.5, pointsize=pt, horizontal=FALSE, onefile=FALSE, paper="special")
}
par(mfrow=c(1, 2))
#### train
dat.train = dt.modeling.overall.train
formula = as.formula(paste("Surv(time, event)~", paste(var.names, collapse="+")))
cox = coxph(formula, data=dat.train, x=TRUE, y=TRUE) 
dat = dt.modeling.treat.train
risk.scores = predict(cox, newdata=dat, type="risk")
group.label = cut(risk.scores, breaks=c(-Inf, cutoff.1, cutoff.2, Inf), labels=1:3)
dat$group = as.numeric(group.label)
fit = survfit(Surv(time, event)~group, data=dat)
surv_diff = survdiff(Surv(time, event)~group, data=dat)
pval = 1 - pchisq(surv_diff$chisq, length(surv_diff$n) - 1)
pval = ifelse(pval<0.0001, "P<0.0001", paste0("P=", pval))
cox = coxph(Surv(time, event)~as.factor(group), data=dat)
HR = cbind(t(sapply(cox_out(cox, dig=2)[,2], function(x)unlist(strsplit(x, split=" \\(|\\)")))), cox_out(cox, dig=3)[,1])
plot_KM(fit, num.label="A", lwd=lwd, col.low=green, col.mid=blue, col.high=red, pval=pval, HR=HR)
#### valid
dat.train = dt.modeling.overall.train
formula = as.formula(paste("Surv(time, event)~", paste(var.names, collapse="+")))
cox = coxph(formula, data=dat.train, x=TRUE, y=TRUE) 
dat = dt.modeling.treat.valid
risk.scores = predict(cox, newdata=dat, type="risk")
group.label = cut(risk.scores, breaks=c(-Inf, cutoff.1, cutoff.2, Inf), labels=1:3)
dat$group = as.numeric(group.label)
fit = survfit(Surv(time, event)~group, data=dat)
surv_diff = survdiff(Surv(time, event)~group, data=dat)
pval = 1 - pchisq(surv_diff$chisq, length(surv_diff$n) - 1)
pval = ifelse(pval<0.0001, "P<0.0001", paste0("P=", pval))
cox = coxph(Surv(time, event)~as.factor(group), data=dat)
HR = cbind(t(sapply(cox_out(cox, dig=2)[,2], function(x)unlist(strsplit(x, split=" \\(|\\)")))), cox_out(cox, dig=3)[,1])
plot_KM(fit, num.label="B", lwd=lwd, col.low=green, col.mid=blue, col.high=red, pval=pval, HR=HR)
dev.off()


############# merged plot, Control
file.name = paste0("KM_总模型_Control队列_2c.", file.type)
if( file.type=="tif" ){
	tiff(file.name, width=W*L*2, height=H*L, pointsize=pt, res=res, compression="lzw")
}else if( file.type=="eps" ){
	pt=20;lwd=4
	postscript(file.name, width=W*1.5*2, height=H*1.5, pointsize=pt, horizontal=FALSE, onefile=FALSE, paper="special")
}
par(mfrow=c(1, 2))
#### train
dat.train = dt.modeling.overall.train
formula = as.formula(paste("Surv(time, event)~", paste(var.names, collapse="+")))
cox = coxph(formula, data=dat.train, x=TRUE, y=TRUE) 
dat = dt.modeling.control.train
risk.scores = predict(cox, newdata=dat, type="risk")
group.label = cut(risk.scores, breaks=c(-Inf, cutoff.1, cutoff.2, Inf), labels=1:3)
dat$group = as.numeric(group.label)
fit = survfit(Surv(time, event)~group, data=dat)
surv_diff = survdiff(Surv(time, event)~group, data=dat)
pval = 1 - pchisq(surv_diff$chisq, length(surv_diff$n) - 1)
pval = ifelse(pval<0.0001, "P<0.0001", paste0("P=", pval))
cox = coxph(Surv(time, event)~as.factor(group), data=dat)
HR = cbind(t(sapply(cox_out(cox, dig=2)[,2], function(x)unlist(strsplit(x, split=" \\(|\\)")))), cox_out(cox, dig=3)[,1])
plot_KM(fit, num.label="A", lwd=lwd, col.low=green, col.mid=blue, col.high=red, pval=pval, HR=HR)
#### valid
dat.train = dt.modeling.overall.train
formula = as.formula(paste("Surv(time, event)~", paste(var.names, collapse="+")))
cox = coxph(formula, data=dat.train, x=TRUE, y=TRUE) 
dat = dt.modeling.control.valid
risk.scores = predict(cox, newdata=dat, type="risk")
group.label = cut(risk.scores, breaks=c(-Inf, cutoff.1, cutoff.2, Inf), labels=1:3)
dat$group = as.numeric(group.label)
fit = survfit(Surv(time, event)~group, data=dat)
surv_diff = survdiff(Surv(time, event)~group, data=dat)
pval = 1 - pchisq(surv_diff$chisq, length(surv_diff$n) - 1)
pval = ifelse(pval<0.0001, "P<0.0001", paste0("P=", pval))
cox = coxph(Surv(time, event)~as.factor(group), data=dat)
HR = cbind(t(sapply(cox_out(cox, dig=2)[,2], function(x)unlist(strsplit(x, split=" \\(|\\)")))), cox_out(cox, dig=3)[,1])
plot_KM(fit, num.label="B", lwd=lwd, col.low=green, col.mid=blue, col.high=red, pval=pval, HR=HR)
dev.off()



###########################################################################################
## methods comparisons
###########################################################################################


library(survcomp)
library(survival)
library(pec)


cox.Cs <- cb.Cs <- rsf.Cs <- gbm.Cs <- svm.Cs <- NULL
cox.rmses <- cb.rmses <- rsf.rmses <- gbm.rmses <- svm.rmses <- NULL
cox.ibss <- cb.ibss <- rsf.ibss <- gbm.ibss <- svm.ibss <- NULL

M = 10
for(i in 1:M)
{
print(i)

n = nrow(dt.modeling.overall.train)
set.seed(i)
train.idx = sample(1:n, size=round(0.7*n), replace=FALSE)

var.selected = c("age", "BMI", "hypertension", "LDLC", "diabetes", "smoke", "zeren", "xiazhai")
formula = as.formula(paste("Surv(time, event)~", paste(var.selected, collapse="+")))
dt.cur = dt.modeling.overall.train[train.idx,]
dt.cur.valid = dt.modeling.overall.train[-train.idx,]

print("Cox")
dat = dt.cur
dat.valid = dt.cur.valid
cox = coxph(formula, data=dat, x=TRUE, y=TRUE)
lp.cox = predict(cox, newdata=dat.valid, type="lp")
C.index.cox = concordance.index(x=lp.cox, surv.time=dat.valid$time, surv.event=dat.valid$event, method="noether")
cox.C.cur = C.index.cox$c.index
cox.Cs = c(cox.Cs, cox.C.cur)

time_grid <- seq(quantile(dat.valid$time, 0.05), quantile(dat.valid$time, 0.95), length.out=100)

# IBS
set.seed(100*i)
ibs_result <- pec(object=list("CoxPH"=cox), formula=Surv(time, event)~1, data=dat.valid, cens.model="marginal", splitMethod = "none",
		times=time_grid)
cox.ibs = ibs(ibs_result)[2,]
cox.ibss = c(cox.ibss, cox.ibs)

#RMSE
set.seed(200*i)
sf = survfit(cox, newdata=dat.valid)
pred_surv_mean <- rowMeans(sf$surv)
km_fit <- survfit(Surv(time, event) ~ 1, data=dat.valid)
km_surv_interp <- approx(km_fit$time, km_fit$surv, xout=time_grid, rule = 2)$y
pred_surv_interp <- approx(sf$time, pred_surv_mean, xout = time_grid, rule = 2)$y
cox.rmse <- sqrt(mean((pred_surv_interp - km_surv_interp)^2))
cox.rmses = c(cox.rmses, cox.rmse)


##################
# CoxBoost
##################
library(CoxBoost)

print("CoxBoost")
dat <- as.data.frame(lapply(dt.cur, function(x) if (is.character(x)) as.factor(x) else x))
dat.valid <- as.data.frame(lapply(dt.cur.valid, function(x) if (is.character(x)) as.factor(x) else x))
x = data.matrix(dat[,which(colnames(dat) %in% var.selected)])
x.valid = data.matrix(dat.valid[,which(colnames(dat.valid) %in% var.selected)])
fit.cb = CoxBoost(time=dat$time, status=dat$event, x=x, stepno=100)
risk.scores.cb = as.vector(predict(fit.cb, newdata=x.valid, type = "lp"))
C.index.cb = concordance.index(x=risk.scores.cb, surv.time=dat.valid$time, surv.event=dat.valid$event, method="noether")
cb.C.cur = C.index.cb$c.index
cb.Cs = c(cb.Cs, cb.C.cur)

# IBS
time_grid <- seq(quantile(dat.valid$time, 0.05), quantile(dat.valid$time, 0.95), length.out=100)
risk_scores = as.vector(predict(fit.cb, newdata=x.valid, type = "lp"))
cox_baseline <- coxph(Surv(dat.valid$time, dat.valid$event)~offset(risk_scores))
base_surv <- survfit(cox_baseline)
base_surv_interp <- approx(base_surv$time, base_surv$surv, xout = time_grid, rule = 2)$y
surv_probs <- outer(base_surv_interp, exp(risk_scores), `^`)  # [time x subject]
surv_probs <- t(surv_probs)  # [subject x time]
km_censor <- survfit(Surv(time, 1 - event) ~ 1, data=dat.valid)
G_t <- approx(km_censor$time, km_censor$surv, xout = time_grid, rule = 2)$y
brier_scores <- numeric(length(time_grid))
for (k in seq_along(time_grid)) {
  t <- time_grid[k]
  s_hat <- surv_probs[, k]
  
  # Observed outcome at time t
  y_t <- ifelse(dat.valid$time <= t & dat.valid$event == 1, 0,
                ifelse(dat.valid$time > t, 1, NA))
  
  # IPCW weights
  G_val <- approx(km_censor$time, km_censor$surv, xout = pmin(dat.valid$time, t), rule = 2)$y
  w <- 1 / G_val
  w[is.na(y_t)] <- 0  # subjects censored before t → don't contribute
  y_t[is.na(y_t)] <- 0  # placeholder; their weight is 0
  
  brier_scores[k] <- mean(w * (s_hat - y_t)^2, na.rm = TRUE)
}
cb.ibs <- sum(diff(time_grid) * (head(brier_scores, -1) + tail(brier_scores, -1)) / 2) / (max(time_grid) - min(time_grid))
cb.ibss = c(cb.ibss, cb.ibs)


# RMSE
risk_scores = as.vector(predict(fit.cb, newdata=x.valid, type = "lp"))
cox_baseline <- coxph(Surv(dat.valid$time, dat.valid$event)~offset(risk_scores))
base_surv <- survfit(cox_baseline)
base_surv_probs <- base_surv$surv

# Compute predicted survival for each subject
pred_surv_mat <- sapply(risk_scores, function(eta) {
  base_surv_probs ^ exp(eta)
})
pred_surv_mean <- rowMeans(pred_surv_mat)
km_fit <- survfit(Surv(time, event) ~ 1, data=dat.valid)
time_grid <- seq(quantile(dat.valid$time, 0.05), quantile(dat.valid$time, 0.95), length.out=100)
km_interp <- approx(km_fit$time, km_fit$surv, xout = time_grid, rule = 2)$y
pred_interp <- approx(base_surv$time, pred_surv_mean, xout = time_grid, rule = 2)$y
cb.rmse <- sqrt(mean((pred_interp - km_interp)^2))
cb.rmses = c(cb.rmses, cb.rmse)





##################
## Random Survival Forests
##################
library(randomForestSRC)
library(survival)
library(pec)
library(prodlim)

print("Random Survival Forests")
dat <- as.data.frame(lapply(dt.cur, function(x) if (is.character(x)) as.factor(x) else x))
dat.valid <- as.data.frame(lapply(dt.cur.valid, function(x) if (is.character(x)) as.factor(x) else x))
set.seed(11)
rsf.fit <- rfsrc(formula, data=dat, ntree=50)
risk.scores.rsf = predict(rsf.fit, newdata=dat.valid)$predicted
C.index.rsf = concordance.index(x=risk.scores.rsf, surv.time=dat.valid$time, surv.event=dat.valid$event, method="noether")
rsf.C.cur = C.index.rsf$c.index
rsf.Cs = c(rsf.Cs, rsf.C.cur)

# rsf_pred$survival is a matrix [n x T], each row is an individual survival curve
surv_probs <- predict(rsf.fit, newdata=dat.valid)$survival   # subject x time

# IBS
km_cens <- survfit(Surv(time, 1 - event) ~ 1, data = dat.valid)
brier_scores <- numeric(length(time_grid))
for (k in seq_along(time_grid)) {
  t <- time_grid[k]
  s_hat <- surv_probs[, k]
  
  # Observed outcome at time t
  y_t <- ifelse(dat.valid$time <= t & dat.valid$event == 1, 0,
                ifelse(dat.valid$time > t, 1, NA))
  
  # IPCW weights
  G_val <- approx(km_censor$time, km_censor$surv, xout = pmin(dat.valid$time, t), rule = 2)$y
  w <- 1 / G_val
  w[is.na(y_t)] <- 0  # subjects censored before t → don't contribute
  y_t[is.na(y_t)] <- 0  # placeholder; their weight is 0
  
  brier_scores[k] <- mean(w * (s_hat - y_t)^2, na.rm = TRUE)
}
rsf.ibs <- sum(diff(time_grid) * (head(brier_scores, -1) + tail(brier_scores, -1)) / 2) / 
       (max(time_grid) - min(time_grid))
rsf.ibss = c(rsf.ibss, rsf.ibs)

# RMSE
pred_surv_mean <- colMeans(surv_probs)  # average survival at each time
km_fit <- survfit(Surv(time, event)~1, data = dat.valid)
km_interp <- approx(km_fit$time, km_fit$surv, xout = time_grid, rule = 2)$y
rsf.rmse <- sqrt(mean((pred_surv_mean - km_interp)^2))
rsf.rmses = c(rsf.rmses, rsf.rmse)



##################
# GBM with Cox Loss
##################
library(gbm)

print("GBM with Cox Loss")
dat <- as.data.frame(lapply(dt.cur, function(x) if (is.character(x)) as.factor(x) else x))
dat.valid <- as.data.frame(lapply(dt.cur.valid, function(x) if (is.character(x)) as.factor(x) else x))
set.seed(123)
gbm_model <- gbm(formula, data=dat, distribution="coxph", n.trees=50, verbose=FALSE)
risk.scores.gbm <- predict(gbm_model, newdata=dat.valid, n.trees=gbm_model$n.trees, type="link")
C.index.gbm = concordance.index(x=risk.scores.gbm, surv.time=dat.valid$time, surv.event=dat.valid$event, method="noether")
gbm.C.cur = C.index.gbm$c.index
gbm.Cs = c(gbm.Cs, gbm.C.cur)


time_grid <- seq(quantile(dat.valid$time, 0.05), quantile(dat.valid$time, 0.95), length.out=100)
risk_scores = predict(gbm_model, newdata=dat.valid, n.trees=gbm_model$n.trees, type="link")
cox_baseline <- coxph(Surv(time, event) ~ offset(risk_scores), data = dat.valid)
base_surv <- survfit(cox_baseline)
base_surv_interp <- approx(base_surv$time, base_surv$surv, xout = time_grid, rule = 2)$y
surv_probs <- outer(base_surv_interp, exp(risk_scores), `^`)  # [time x subject]
surv_probs <- t(surv_probs)  # [subject x time]

# IBS
km_cens <- survfit(Surv(time, 1 - event) ~ 1, data = dat.valid)
brier_scores <- numeric(length(time_grid))
for (k in seq_along(time_grid)) {
  t <- time_grid[k]
  s_hat <- surv_probs[, k]
  
  # Observed outcome at time t
  y_t <- ifelse(dat.valid$time <= t & dat.valid$event == 1, 0,
                ifelse(dat.valid$time > t, 1, NA))
  
  # IPCW weights
  G_val <- approx(km_censor$time, km_censor$surv, xout = pmin(dat.valid$time, t), rule = 2)$y
  w <- 1 / G_val
  w[is.na(y_t)] <- 0  # subjects censored before t → don't contribute
  y_t[is.na(y_t)] <- 0  # placeholder; their weight is 0
  
  brier_scores[k] <- mean(w * (s_hat - y_t)^2, na.rm = TRUE)
}
gbm.ibs <- sum(diff(time_grid) * (head(brier_scores, -1) + tail(brier_scores, -1)) / 2) / 
       (max(time_grid) - min(time_grid))
gbm.ibss = c(gbm.ibss, gbm.ibs)

# RMSE
pred_surv_mean <- colMeans(surv_probs)
km_fit <- survfit(Surv(time, event) ~ 1, data = dat.valid)
km_interp <- approx(km_fit$time, km_fit$surv, xout = time_grid, rule = 2)$y
gbm.rmse <- sqrt(mean((pred_surv_mean - km_interp)^2))
gbm.rmses = c(gbm.rmses, gbm.rmse)




##################
# Survival-SVM
##################
library(survivalsvm)
library(quadprog)
library(survcomp)

print("Survival-SVM")
survsvm.reg <- survivalsvm(formula=formula, data=dat, gamma.mu=0.1, eig.tol=1e-03, conv.tol=1e-04, posd.tol=1e-05)
pred.survsvm.reg <- predict(object=survsvm.reg, newdata=dat.valid)
risk.scores.svm = as.vector(pred.survsvm.reg$predicted)
C.index.svm = concordance.index(x=risk.scores.svm, surv.time=dat.valid$time, surv.event=dat.valid$event, method="noether")
svm.C.cur = C.index.svm$c.index
svm.Cs = c(svm.Cs, svm.C.cur)

# IBS
time_grid <- seq(quantile(dat.valid$time, 0.05), quantile(dat.valid$time, 0.95), length.out=100)
risk_scores = risk.scores.svm
Q = 5
dat.valid$risk_group <- cut(risk_scores, breaks = quantile(risk_scores, probs = seq(0, 1, by=1/Q)),
			include.lowest = TRUE, labels = paste0("Q", 1:Q))

km_Q1 <- survfit(Surv(time, event) ~ 1, data = dat.valid[dat.valid$risk_group == "Q1", ])
km_Q2 <- survfit(Surv(time, event) ~ 1, data = dat.valid[dat.valid$risk_group == "Q2", ])
km_Q3 <- survfit(Surv(time, event) ~ 1, data = dat.valid[dat.valid$risk_group == "Q3", ])
km_Q4 <- survfit(Surv(time, event) ~ 1, data = dat.valid[dat.valid$risk_group == "Q4", ])
km_Q5 <- survfit(Surv(time, event) ~ 1, data = dat.valid[dat.valid$risk_group == "Q5", ])


# Build subject-by-time predicted survival matrix
assign_surv_probs <- function(km_fit, times) {
  approx(km_fit$time, km_fit$surv, xout = times, rule = 2)$y
}

surv_Q1 <- assign_surv_probs(km_Q1, time_grid)
surv_Q2 <- assign_surv_probs(km_Q2, time_grid)
surv_Q3 <- assign_surv_probs(km_Q3, time_grid)
surv_Q4 <- assign_surv_probs(km_Q4, time_grid)
surv_Q5 <- assign_surv_probs(km_Q5, time_grid)

surv_mat <- matrix(NA, nrow = nrow(dat.valid), ncol = length(time_grid))
for (i in seq_len(nrow(dat.valid))) {
  if (dat.valid$risk_group[i] == "Q1") {
    surv_mat[i, ] <- surv_Q1
  } else if (dat.valid$risk_group[i] == "Q2") {
    surv_mat[i, ] <- surv_Q2
  } else if (dat.valid$risk_group[i] == "Q3") {
    surv_mat[i, ] <- surv_Q3
  } else if (dat.valid$risk_group[i] == "Q4") {
    surv_mat[i, ] <- surv_Q4
  } else {
    surv_mat[i, ] <- surv_Q5
  }
}


brier_scores <- numeric(length(time_grid))
km_cens <- survfit(Surv(time, 1 - event) ~ 1, data = dat.valid)
for (k in seq_along(time_grid)) {
  t <- time_grid[k]
  s_hat <- surv_mat[, k]
  
  # Binary outcome
  y_t <- ifelse(dat.valid$time <= t & dat.valid$event == 1, 0,
                ifelse(dat.valid$time > t, 1, NA))
  
  # IPCW weights
  G_val <- approx(km_cens$time, km_cens$surv, xout = pmin(dat.valid$time, t), rule = 2)$y
  w <- 1 / G_val
  w[is.na(y_t)] <- 0
  y_t[is.na(y_t)] <- 0
  
  brier_scores[k] <- mean(w * (s_hat - y_t)^2, na.rm = TRUE)
}
svm.ibs <- sum(diff(time_grid) * (head(brier_scores, -1) + tail(brier_scores, -1)) / 2) / (max(time_grid) - min(time_grid))
svm.ibss = c(svm.ibss, svm.ibs)

# RMSE
km_all <- survfit(Surv(time, event) ~ 1, data = dat.valid)
km_all_interp <- approx(km_all$time, km_all$surv, xout = time_grid, rule = 2)$y
pred_surv_mean <- colMeans(surv_mat)
svm.rmse <- sqrt(mean((pred_surv_mean - km_all_interp)^2))
svm.rmses = c(svm.rmses, svm.rmse)
}


stat_mean_sd = function(x, dig=3)
{
	res = sprintf(paste0("%.", dig, "f (%.", dig, "f)"), mean(x), sd(x))
	return(res)
}

dig = 3
Cs = c(stat_mean_sd(cox.Cs, dig=dig), stat_mean_sd(cb.Cs, dig=dig), stat_mean_sd(rsf.Cs, dig=dig), stat_mean_sd(gbm.Cs, dig=dig), stat_mean_sd(svm.Cs, dig=dig))
ibss = c(stat_mean_sd(cox.ibss, dig=dig), stat_mean_sd(cb.ibss, dig=dig), stat_mean_sd(rsf.ibss, dig=dig), stat_mean_sd(gbm.ibss, dig=dig), stat_mean_sd(svm.ibss, dig=dig))
rmses = c(stat_mean_sd(cox.rmses, dig=dig), stat_mean_sd(cb.rmses, dig=dig), stat_mean_sd(rsf.rmses, dig=dig), stat_mean_sd(gbm.rmses, dig=dig), stat_mean_sd(svm.rmses, dig=dig))
res = cbind(Cs, ibss, rmses)
colnames(res) = c("C-index", "IBS", "RMSE")
rownames(res) = c("Cox", "CoxBoot", "Random Survival Forest", "GBM", "Survival-SVM")
write.csv(res, file="comparison.csv")

