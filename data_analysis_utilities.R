time_zero = function(time)
{
	set = unlist(strsplit(as.character(time), split="-"))
	if( "00" %in% set | "0000" %in% set )  return(TRUE)
	return(FALSE)
}


date_missed_single = function(date)
{
	date = as.character(date[[1]])
	return( date == "" | time_zero(date) )
}


date_missed = function(line, targets_col)
{
	return( any( sapply(targets_col, function(x){date_missed_single(line[x])}) ) )
}


data_used = function(data, cols_name, event_names, truncate.year=999)
{
	dt = data[, as.numeric(sapply(cols_name, function(x){which(names(data)==x)}))]

	events = NULL
	dates = NULL
	missing_dates = NULL
	times = NULL
	for(i in 1:nrow(dt))
	{
		line = dt[i,]
		missing_date = 0
		event_cols = sapply(event_names, function(x)which(names(dt) == paste0(x, "_yes")))

		if( any(as.logical(line[event_cols])) ){
			event = 1
			events_happened = event_names[which(line[event_cols]==1)]
			events_happened_date_cols = sapply(events_happened, function(x)which(names(dt) == paste0(x, "_date")))
			if( date_missed(line, events_happened_date_cols) ){
				date = ""
				missing_date = 1
			}else{
				X = line[events_happened_date_cols]
				date = as.Date(as.character(X[[1]]))
				for(x in X) date = min(date, as.Date(as.character(x)))
				date = as.character(date)
			}
		}else{
			event = 0
			if( line$death_yes ){
				date = as.character(line$death_date)
			}else if( line$chexiao ){
				date = as.character(line$chexiao.date)
			}else if( line$missings ){
				date = as.character(line$missings.date)
			}else{
				date = as.character(line$end_date)
			}
			if( date == "" | time_zero(date) ) missing_date = 1
		}

		if( missing_date ){
			time = 0
		}else{
			time = lubridate::time_length(difftime(as.Date(date), as.Date(line$enroll_date), units="days"), "years")
			if( time>truncate.year ){
				time = truncate.year
				event = 0
				date = ymd(as.Date(line$enroll_date)) + years(truncate.year)
			}
		}

		events = c(events, event)
		dates = c(dates, date)
		missing_dates = c(missing_dates, missing_date)
		times = c(times, time)
	}

	dt$event = events
	dt$missing_date = missing_dates
	dt$time = times
	return(dt)
}


reformat_data = function(dt)
{
	dt$age = round(dt$age)
	dt$BMI = round(dt$BMI,1)

	zeren = as.vector(sapply(dt$zeren, function(x)unlist(strsplit(x, split="（"))[1]))
	zeren[which(zeren %in% c("颈内动脉", "大脑中动脉"))] = 1		#前循环
	zeren[which(zeren %in% c("基底动脉", "椎动脉"))] = 2		#后循环
	dt$zeren = zeren

	# xiazhai="无" is >= 70%
	dt$xiazhai[which(dt$xiazhai=="有")] = 1
	dt$xiazhai[which(dt$xiazhai=="无")] = 2

	dt$hypertension[which(dt$hypertension=="无")] = 1
	dt$hypertension[which(dt$hypertension=="有")] = 2

	dt$diabetes[which(dt$diabetes=="无")] = 1
	dt$diabetes[which(dt$diabetes=="有")] = 2

	dt$smoke[which(dt$smoke=="无")] = 1
	dt$smoke[which(dt$smoke=="有")] = 2

	dt$TIA[which(dt$TIA=="无")] = 1
	dt$TIA[which(dt$TIA=="有")] = 2

	dt$CIS[which(dt$CIS=="无")] = 1
	dt$CIS[which(dt$CIS=="有")] = 2

	dt$AMI[which(dt$AMI=="无")] = 1
	dt$AMI[which(dt$AMI=="有")] = 2

	dt$gaozhi[which(dt$gaozhi=="无")] = 1
	dt$gaozhi[which(dt$gaozhi=="有")] = 2

	dt$gender[which(dt$gender=="男")] = 1
	dt$gender[which(dt$gender=="女")] = 2

	dt$LDLC = round(dt$LDLC)
	dt$HDLC = round(dt$HDLC)
	dt$TCHO = round(dt$TCHO)

	return(dt)
}


tab_baseline = function(dat, dig=1, dig.P=3)
{
	event.idx = which(dat$event==1)
	baseline = baseline_stat(dat, dig=dig)
	event.baseline = baseline_stat(dat[event.idx,], dig=dig)
	nonevent.baseline = baseline_stat(dat[-event.idx,], dig=dig)
	out.tab = cbind(baseline, event.baseline, nonevent.baseline, diff_test_all(dat, dig.P=dig.P))

	return(out.tab)
}


baseline_stat = function(dat, dig=1, na=FALSE)
{
	n = nrow(dat)

	age.stat = stat_conti(dat$age, dig, na, topic="Age")	# age
	gender.stat = stat_count(dat$gender, order=c("2", "1"), n=n, dig=dig, topic="Gender")	# gender
	BMI.stat = stat_conti(dat$BMI, dig, na, topic="BMI")	# BMI
	enroll.stat = stat_count(dat$enroll_event, order=c(1, 2), n=n, dig=dig, topic="Enroll.Event")	# enroll-event
	TIA.gap.time.stat = stat_conti(dat$gap_time[which(dat$enroll_event==1)], dig, na, topic="Time-to-Rand(TIA)")
	CIS.gap.time.stat = stat_conti(dat$gap_time[which(dat$enroll_event==2)], dig, na, topic="Time-to-Rand(CIS)")
	ABCD2.stat = stat_conti(dat[which(dat$enroll_event==1),]$ABCD2_score, dig, na, topic="ABCD2")	# ABCD2
	mRS.stat = stat_conti(dat[which(dat$enroll_event==2),]$mRS_score, dig, na, topic="mRS")	# mRS
	stroke.stat = stat_count(dat$CIS, order="2", n=n, dig=dig, topic="Stroke")	## stroke
	TIA.stat = stat_count(dat$TIA, order="2", n=n, dig=dig, topic="TIA")	## TIA
	AMI.stat = stat_count(dat$AMI, order="2", n=n, dig=dig, topic="AMI")	## AMI
	hypertension.stat = stat_count(dat$hypertension, order="2", n=n, dig=dig, topic="Hypertension")	## Ischemic stroke
	gaozhi.stat = stat_count(dat$gaozhi, order="2", n=n, dig=dig, topic="Gaozhi")	## Gaozhixue Zheng
	diabetes.stat = stat_count(dat$diabetes, order="2", n=n, dig=dig, topic="Diabetes")	## diabetes
	smoke.stat = stat_count(dat$smoke, order="2", n=n, dig=dig, topic="Smoke")	# smoke,　是否有吸烟史
	dongmai.stat = stat_count(dat$zeren, order="2", n=n, dig=dig, topic="Xueguan")	# 后循环
	xiazhai.stat = stat_count(dat$xiazhai, order="2", n=n, dig=dig, topic="Xiazhai")	## Xiazhai Chengdu
	SP.stat = stat_conti(dat$SP, dig, na, topic="SP")	# Systolic pressure
	DP.stat = stat_conti(dat$DP, dig, na, topic="DP")	# Diastolic pressure
	FBG.stat = stat_conti(dat$FBG, dig, na, topic="FBG")	# FBG
	LDLC.stat = stat_conti(dat$LDLC, dig, na, topic="LDLC")	# LDLC
	HDLC.stat = stat_conti(dat$HDLC, dig, na, topic="HDLC")	# HDLC
	TCHO.stat = stat_conti(dat$TCHO, dig, na, topic="TCHO")	# TCHO

	stats = c(n, "", age.stat[1], "", gender.stat, BMI.stat[1], "", enroll.stat, "", TIA.gap.time.stat[2], 
			CIS.gap.time.stat[2], "", ABCD2.stat[2], mRS.stat[2], "", stroke.stat, TIA.stat,
			AMI.stat, hypertension.stat, gaozhi.stat, diabetes.stat, smoke.stat, dongmai.stat, xiazhai.stat, FBG.stat[1], "", SP.stat[1], DP.stat[1], 
			"", LDLC.stat[1], HDLC.stat[1], TCHO.stat[1])
	return(stats)
}


stat_count = function(seq, order, n, dig=1, topic=NULL)
{
	if( length(order)>1 ){
		seq.no = table(seq)
		ranks = sapply(names(seq.no), function(x){which(x==order)})
		seq.no = seq.no[ranks]	
	}else{
		seq.no = sum(seq==order)
	}
	seq.rate = round(seq.no / n * 100, dig=dig)

	res = NULL
	for(i in 1:length(seq.no)) res = c(res, paste0(seq.no[i], " (", seq.rate[i], "%)"))
	names(res) = paste0(topic, "-", order, "-No(%)")

	return(res)
}


stat_conti = function(seq, dig=1, na=TRUE, topic=NULL)
{
	seq.mid = round(median(seq, na.rm=TRUE), dig=dig)
	seq.quats = round(quantile(seq, probs=c(0.25, 0.75), na.rm=TRUE), dig=dig)
	seq.mean = round(mean(seq, na.rm=TRUE), dig=dig)
	seq.sd = round(sd(seq, na.rm=TRUE), dig=dig)
	if( na ){
		seq.na = sum(is.na(seq))
		res = c(seq.mid, seq.quats, seq.mean, seq.sd, seq.na)
		names(res) = c(paste0(topic, "-Median"), paste0(topic, "-25%"), paste0(topic, "-75%"), 
						paste0(topic, "-Mean"), paste0(topic, "-SD"), paste0(topic, "-#NA"))
	}else{
		res = c(paste0(seq.mean, " (", seq.sd, ")"), paste0(seq.mid, " (", seq.quats[1], "-", seq.quats[2], ")"))
		names(res) = c(paste0(topic, "-Mean(SD)"), paste0(topic, "-Median(IQR)"))
	}

	return(res)
}


diff_test_all = function(dat, dig.P=3)
{
	event.idx = which(dat$event==1)

	age.test = diff_test(dat$age[event.idx], dat$age[-event.idx], type="mean", dig.P)
	gender.test = diff_test(dat$gender[event.idx], dat$gender[-event.idx], type="proportion", dig.P)
	BMI.test = diff_test(dat$BMI[event.idx], dat$BMI[-event.idx], type="mean", dig.P)
	enroll.test = diff_test(dat$enroll_event[event.idx], dat$enroll_event[-event.idx], type="proportion", dig.P)
	TIA.gap.time.test = diff_test(dat$gap_time[which(dat$enroll_event==1&dat$event==1)], dat$gap_time[which(dat$enroll_event==1&dat$event==0)], type="median", dig.P)
	CIS.gap.time.test = diff_test(dat$gap_time[which(dat$enroll_event==2&dat$event==1)], dat$gap_time[which(dat$enroll_event==2&dat$event==0)], type="median", dig.P)
	ABCD2.test = diff_test(dat[which(dat$enroll_event==1&dat$event==1),]$ABCD2_score, dat[which(dat$enroll_event==1&dat$event==0),]$ABCD2_score, type="median", dig.P)
	mRS.test = diff_test(dat[which(dat$enroll_event==2&dat$event==1),]$mRS_score, dat[which(dat$enroll_event==2&dat$event==0),]$mRS_score, type="median", dig.P)
	stroke.test = diff_test(dat$CIS[event.idx], dat$CIS[-event.idx], type="proportion", dig.P)
	TIA.test = diff_test(dat$TIA[event.idx], dat$TIA[-event.idx], type="proportion", dig.P)
	AMI.test = diff_test(dat$AMI[event.idx], dat$AMI[-event.idx], type="proportion", dig.P)
	hypertension.test = diff_test(dat$hypertension[event.idx], dat$hypertension[-event.idx], type="proportion", dig.P)
	gaozhi.test = diff_test(dat$gaozhi[event.idx], dat$gaozhi[-event.idx], type="proportion", dig.P)
	diabetes.test = diff_test(dat$diabetes[event.idx], dat$diabetes[-event.idx], type="proportion", dig.P)
	smoke.test = diff_test(dat$smoke[event.idx], dat$smoke[-event.idx], type="proportion", dig.P)
	xiazhai.test = diff_test(dat$xiazhai[event.idx], dat$xiazhai[-event.idx], type="proportion", dig.P)
	SP.test = diff_test(dat$SP[event.idx], dat$SP[-event.idx], type="mean", dig.P)
	DP.test = diff_test(dat$DP[event.idx], dat$DP[-event.idx], type="mean", dig.P)
	FBG.test = diff_test(dat$FBG[event.idx], dat$FBG[-event.idx], type="mean", dig.P)
	LDLC.test = diff_test(dat$LDLC[event.idx], dat$LDLC[-event.idx], type="mean", dig.P)
	HDLC.test = diff_test(dat$HDLC[event.idx], dat$HDLC[-event.idx], type="mean", dig.P)
	TCHO.test = diff_test(dat$TCHO[event.idx], dat$TCHO[-event.idx], type="mean", dig.P)
	dongmai.test = diff_test(dat$zeren[event.idx], dat$zeren[-event.idx], type="proportion", dig.P)

	stats = c("", "", age.test, "", "", gender.test, BMI.test, "", "", enroll.test, 
			"", TIA.gap.time.test, CIS.gap.time.test, "", ABCD2.test, mRS.test, 
			"", stroke.test, TIA.test, AMI.test, hypertension.test, gaozhi.test, diabetes.test, smoke.test, 
			dongmai.test, xiazhai.test, FBG.test, 
			"", SP.test, DP.test, 
			"", LDLC.test, HDLC.test, TCHO.test)
	return(stats)
}


diff_test = function(group1, group2, type="median", dig=3)
{
	if( type=="median" ){
		test.res = wilcox.test(group1, group2)$p.value
	}else if( type=="mean" ){
		test.res = t.test(group1, group2)$p.value
	}else if( type=="proportion" ){
		contingency.table = rbind(table(group1), table(group2))
		test.res = chisq.test(contingency.table)$p.value
	}
	res = ifelse(round(test.res, dig)==0, paste0(c("<0.", rep(0, dig), "1"), collapse=''), sprintf(paste0("%.", dig, "f"), test.res))

	return(res)
}


diff_test_overall_all = function(dat.train, dat.valid, dig.P=3)
{
	age.test = diff_test(dat.train$age, dat.valid$age, type="mean", dig.P)
	gender.test = diff_test(dat.train$gender, dat.valid$gender, type="proportion", dig.P)
	BMI.test = diff_test(dat.train$BMI, dat.valid$BMI, type="mean", dig.P)
	enroll.test = diff_test(dat.train$enroll_event, dat.valid$enroll_event, type="proportion", dig.P)
	TIA.gap.time.test = diff_test(dat.train$gap_time[which(dat.train$enroll_event==1)], dat.valid$gap_time[which(dat.valid$enroll_event==1)], type="median", dig.P)
	CIS.gap.time.test = diff_test(dat.train$gap_time[which(dat.train$enroll_event==2)], dat.valid$gap_time[which(dat.valid$enroll_event==2)], type="median", dig.P)
	ABCD2.test = diff_test(dat.train[which(dat.train$enroll_event==1),]$ABCD2_score, dat.valid[which(dat.valid$enroll_event==1),]$ABCD2_score, type="median", dig.P)
	mRS.test = diff_test(dat.train[which(dat.train$enroll_event==2),]$mRS_score, dat.valid[which(dat.valid$enroll_event==2),]$mRS_score, type="median", dig.P)
	stroke.test = diff_test(dat.train$CIS, dat.valid$CIS, type="proportion", dig.P)
	TIA.test = diff_test(dat.train$TIA, dat.valid$TIA, type="proportion", dig.P)
	AMI.test = diff_test(dat.train$AMI, dat.valid$AMI, type="proportion", dig.P)
	hypertension.test = diff_test(dat.train$hypertension, dat.valid$hypertension, type="proportion", dig.P)
	gaozhi.test = diff_test(dat.train$gaozhi, dat.valid$gaozhi, type="proportion", dig.P)
	diabetes.test = diff_test(dat.train$diabetes, dat.valid$diabetes, type="proportion", dig.P)
	smoke.test = diff_test(dat.train$smoke, dat.valid$smoke, type="proportion", dig.P)
	xiazhai.test = diff_test(dat.train$xiazhai, dat.valid$xiazhai, type="proportion", dig.P)
	SP.test = diff_test(dat.train$SP, dat.valid$SP, type="mean", dig.P)
	DP.test = diff_test(dat.train$DP, dat.valid$DP, type="mean", dig.P)
	FBG.test = diff_test(dat.train$FBG, dat.valid$FBG, type="mean", dig.P)
	LDLC.test = diff_test(dat.train$LDLC, dat.valid$LDLC, type="mean", dig.P)
	HDLC.test = diff_test(dat.train$HDLC, dat.valid$HDLC, type="mean", dig.P)
	TCHO.test = diff_test(dat.train$TCHO, dat.valid$TCHO, type="mean", dig.P)
	dongmai.test = diff_test(dat.train$zeren, dat.valid$zeren, type="proportion", dig.P)

	stats = c("", "", age.test, "", "", gender.test, BMI.test, "", "", enroll.test, 
			"", TIA.gap.time.test, CIS.gap.time.test, "", ABCD2.test, mRS.test, 
			"", stroke.test, TIA.test, AMI.test, hypertension.test, gaozhi.test, diabetes.test, smoke.test, 
			dongmai.test, xiazhai.test, FBG.test, 
			"", SP.test, DP.test, 
			"", LDLC.test, HDLC.test, TCHO.test)
	return(stats)
}


modeling = function(dt, backward=TRUE, dig=3)
{
	var.names = c("age", "BMI", "hypertension", "LDLC", "diabetes", "smoke", "zeren", "xiazhai", 
			"gaozhi", "SP", "TIA", "TCHO", 
			"DP", "gender", "enroll_event", "CIS", "AMI", "FBG", "HDLC")
	# univariate regression
	single.res = NULL
	for( var in var.names )
	{
		formula = as.formula(paste("Surv(time, event)~", var))
		cox = coxph(formula, data=dt)
		single.res = rbind(single.res, cox_out(cox, dig))
	}

	# multiple regression
	var.names = c("age", "BMI", "hypertension", "LDLC", "diabetes", "smoke", "zeren", "xiazhai", 
			"gaozhi", "SP", "TIA", "TCHO")
	formula = as.formula(paste("Surv(time, event)~", paste(var.names, collapse="+")))
	cox = coxph(formula, data=dt)
	multi.res = cox_out(cox, dig)

	res = cbind(single.res, 
		rbind(cbind(rownames(multi.res), multi.res), matrix("", ncol=3, nrow=nrow(single.res)-nrow(multi.res)))
	)
	colnames(res) = c("Single-Pvalue", "Single-HR", 
		"", "Multi-Pvalue", "Multi-HR")


	if( backward ){
		# backwards regression with AIC
		stepwise_model <- stepAIC(cox, direction="backward", k=2, trace=0)
		stepres.AIC = cox_out(stepwise_model, dig)

		res = cbind(single.res, 
			rbind(cbind(rownames(multi.res), multi.res), matrix("", ncol=3, nrow=nrow(single.res)-nrow(multi.res))),
			rbind(cbind(rownames(stepres.AIC), stepres.AIC), matrix("", ncol=3, nrow=nrow(single.res)-nrow(stepres.AIC)))
		)

		colnames(res) = c("Single-Pvalue", "Single-HR", 
			"", "Multi-Pvalue", "Multi-HR", 
			"", "AIC-Pvalue", "AIC-HR")
	}
	
	return(res)
}


cox_out = function(cox, dig=3)
{
	cox_summary = summary(cox)

	names = rownames(cox_summary$coef)
	pos = sapply(round(cox_summary$coef[,5], dig), function(x)identical(x,0))
	pvals = sprintf(paste0("%.", dig, "f"), cox_summary$coef[,5])
	pvals[pos] = paste0("<0.", strrep("0", dig), "1")
	HR = sprintf(paste0("%.", dig, "f"), cox_summary$coef[,2])
	HR_ci_lower = sprintf(paste0("%.", dig, "f"), cox_summary$conf.int[,"lower .95"])
	HR_ci_upper = sprintf(paste0("%.", dig, "f"), cox_summary$conf.int[,"upper .95"])
	HR = paste0(HR, " (", HR_ci_lower, "-", HR_ci_upper, ")")
	res = cbind(pvals, HR)
	colnames(res) = c("P-value", "HR")
	rownames(res) = names

	return(res)
}


C_index_CI = function(dat, vars.all, B=100, dig=3, dig.P=4, only.overall=FALSE)
{
	boot.res = boot(data=dat, statistic=C_index_fun, var.names=vars.all, only.overall=only.overall, R=B)
	C.index = boot.res$t0
	C.CI = apply(boot.res$t, 2, quantile, probs=c(0.025, 0.975))
	res = sapply(1:length(C.index), function(i)sprintf(paste0("%.", dig, "f (%.", dig, "f-%.", dig, "f)"), C.index[i], C.CI[1,i], C.CI[2,i]))
	var.selected = vars.all[1:8]

	if( !only.overall ){
		# p-values
		formula = as.formula(paste("Surv(time, event)~", paste(var.selected, collapse="+")))
		cox = coxph(formula, data=dat, x=TRUE, y=TRUE)
		lp1 = predict(cox, type="lp")

		pvals = NULL
		for(var in var.selected){
			#print(var)
			formula = as.formula(paste("Surv(time, event)~", paste(var, collapse="+")))
			cox.tmp = coxph(formula, data=dat, x=TRUE, y=TRUE)
			lp2 = predict(cox.tmp, type="lp")
			test.tmp = compareC(dat$time, dat$event, lp1, lp2)$pval
			pval.tmp = ifelse(identical(round(test.tmp, dig.P), 0), paste0("<0.", paste0(rep(0, 3), collapse=""), "1"), sprintf(paste0("%.", dig, "f"), test.tmp))
			pvals = c(pvals, pval.tmp)
		}
		res = cbind(res, c("", pvals))
	}

	return(res)
}


C_index_fun = function(dat, indices, var.names, only.overall=FALSE)
{
	dat = dat[indices,]
	formula = as.formula(paste("Surv(time, event)~", paste(var.names, collapse="+")))
	cox = coxph(formula, data=dat, x=TRUE, y=TRUE)
	C.index = cox$concordance[names(cox$concordance)=="concordance"]

	if( !only.overall ){
		for(var in var.names[1:8]){
			formula = as.formula(paste("Surv(time, event)~", paste(var, collapse="+")))
			cox.tmp = coxph(formula, data=dat, x=TRUE, y=TRUE)
			C.tmp = cox.tmp$concordance[names(cox.tmp$concordance)=="concordance"]
			C.index = c(C.index, C.tmp)
		}
	}

	return(C.index)
}


plot_cal = function(dat, dat.val, vars, u=12, xlim=c(0,1), ylim=c(0,1), lwd=2, 
			col.train="#5B88C3", col.val="#C7543C", col.ideal="black", 
			num.label="", no.y.labels=FALSE, merged=FALSE)
{
	formula = as.formula(paste("Surv(time, event)~", paste(vars, collapse="+")))
	cox.rms = cph(formula, data=dat, x=TRUE, y=TRUE, surv=TRUE, time.inc=u)
	cal = calibrate(cox.rms, method="boot", B=B, u=u, pr=FALSE)

	par(mgp=c(2, 1, 0))
	par(mar=c(4, 4, 2, 2))
	xlab=paste0("Predicted ", u, "-month occurrence-free survival probability")
	ylab="Actual probability"
	if( no.y.labels ) ylab=""
	plot(cal[,"pred"], cal[,"calibrated"], xlim=xlim, ylim=ylim, xlab=xlab, ylab=ylab, type="n", axes=F)
	segments(x0=0, y0=0, x1=1, y1=1, lty="22", lwd=lwd, col=col.ideal)
	seq.whole = seq(0, 1, by=0.2)
	axis(side=1, at=seq.whole, labels=c("0", seq.whole[-c(1,6)], "1.0"), pos=0, lwd=lwd, lwd.ticks=lwd)	
	if( no.y.labels ){
		axis(side=2, at=seq.whole, labels=c("0", seq.whole[-c(1,6)], "1.0"), pos=0, lwd=lwd, lwd.ticks=lwd, las=2)
	}else{
		axis(side=2, at=seq.whole, labels=c("0", seq.whole[-c(1,6)], "1.0"), pos=0, lwd=lwd, lwd.ticks=lwd, las=2)
	}
	mtext(num.label, side=2, line=0.5-0.07, at=1.1, adj=1, cex=1, las=2)

	x.cal = cal[,"pred"]
	y.cal = cal[,"calibrated"]
	lines(x.cal, y.cal, col=col.train, lwd=lwd)

	S = Surv(dat.val$time, dat.val$event)
	cal.val = val.surv(cox.rms, newdata=dat.val, S=S, u=u)
	x.cal.val = cal.val$pseq
	y.cal.val = cal.val$actualseq
	lines(x.cal.val, y.cal.val, col=col.val, lwd=lwd)

	legend("bottomright", c("Training set", "External validation set", "Ideal"), lwd=lwd, 
	      lty=c(NA,NA,"22"), col=c(col.train, col.val, col.ideal), bty="n", 
		inset=c(0.05, 0.05), cex=0.7)
	legend("bottomright", c("Training set", "External validation set", "Ideal"), lwd=lwd, 
	      lty=c(1,1,NA), col=c(col.train, col.val, col.ideal), bty="n", 
		inset=c(0.05, 0.05), cex=0.7)
	#abline(v=0.6)
}


plot_DCA = function(dt, cox, u=12, thresholds=seq(0, 1, 0.01), span=0.3, lwd=2, 
			col.all="#5B88C3", col.none="#C7543C", col.model="black", 
			num.label="", no.y.labels=FALSE, merged=FALSE, no.ylab=FALSE)
{
	dt$fitted = c(1 - (summary(survfit(cox,newdata=dt), times=u)$surv))
	p.0 = dca(Surv(time, event)~fitted,
	    data=dt,
	    time=u,
	    thresholds=thresholds,
	    label=list(fitted="Prediction Model")
	) 

	rm.idx = which(p.0$dca$label=="Prediction Model" & (sapply(p.0$dca$net_benefit, identical, y=0) | is.na(p.0$dca$net_benefit)))
	if( length(rm.idx)>0 ) p.0$dca = p.0$dca[-rm.idx,]
	p.0 =  plot(p.0, smooth=TRUE, span=span)

	g = ggplot_build(p.0)$data[[1]]
	all.y = g$y[g$group==1];all.x = g$x[g$group==1]
	all.x = all.x[all.y>-0.05];all.y = all.y[all.y>-0.05]
	none.y = g$y[g$group==2];none.x = g$x[g$group==2]
	model.y = g$y[g$group==3];model.x = g$x[g$group==3]

	par(mgp=c(2, 1, 0))
	gap.vec = rep(0,4)
	if( merged ) gap.vec[1] = 2
	if( no.y.labels ) gap.vec[2] = 2
	par(mar=c(4, 4, 2, 2)-gap.vec)
	tmp = seq(0, 0.5, by=0.05)
	y.up = (tmp[tmp>max(g$y)])[1]
	xlim = c(0,1); ylim=c(-0.05, y.up)
	xlab = ifelse(merged, "", "Threshold probability"); ylab=ifelse(no.y.labels, "", "Net benefit")
	plot(none.x, none.y, xlim=xlim, ylim=ylim, xlab=xlab, ylab="", type="n", axes=FALSE)
	if(!no.ylab) mtext(ylab, side=2, line=2.5)
	lines(none.x, none.y, col=col.none, lwd=lwd)
	lines(all.x, all.y, col=col.all, lwd=lwd)
	lines(model.x, model.y, col=col.model, lwd=lwd)

	seq.whole = seq(0, 1, by=0.2)
	if(u==3){
		F = 0.01
	}else if(u==6){
		F = 0.015
	}else if(u==12){
		F = 0.02
	}

	axis(side=1, at=seq.whole, labels=c("0", seq.whole[-c(1,6)], "1.0"), pos=-0.050, lwd=lwd, lwd.ticks=lwd)	
	if( no.y.labels ){
		axis(side=2, at=seq(-0.05, y.up, by=0.05), pos=0, lwd=lwd, lwd.ticks=lwd, las=2)
	}else{
		axis(side=2, at=seq(-0.05, y.up, by=0.05), pos=0, lwd=lwd, lwd.ticks=lwd, las=2)
	}
	mtext(num.label, side=2, at=y.up+F, line=1.3, cex=1, las=2, adj=0)
	
	legend("topright", c("Treat all", "Treat none", "Prediction model"), lwd=lwd, 
	      lty=1, col=c(col.all, col.none, col.model), bty="n", cex=1)
}


label_name = function(u)
{
	if( u==3 ){
		num.label="A Training set (3 month)"
		num.label.val="D External validation set (3 month)"
	}else if( u==6 ){
		num.label="B Training set (6 month)"
		num.label.val="E External validation set (6 month)"
	}else if( u==12 ){
		num.label="C Training set (12 month)"
		num.label.val="F External validation set (12 month)"
	}else{
		num.label=""
		num.label.val=""
	}

	res = c(num.label, num.label.val)
	return(res)
}


plot_KM = function(f, num.label="A", lwd=2, col.low="green", col.mid="blue", col.high="red", pval="<0.0001", 
			HR=NULL, no.ylab=FALSE, no.xlab=FALSE)
{
	x = f$time;y = f$surv

	starts = sapply(f$n, function(x)which(f$n.risk==x))
	idx.low = starts[1]:(starts[2]-1)
	idx.mid = starts[2]:(starts[3]-1)
	idx.high = starts[3]:length(f$n.risk)
	x.low=x[idx.low];x.mid=x[idx.mid];x.high=x[idx.high]
	y.low=y[idx.low];y.mid=y[idx.mid];y.high=y[idx.high]

	n.risks = f$n.risk
	n.risks.low = n.risks[idx.low]
	n.risks.mid = n.risks[idx.mid]
	n.risks.high = n.risks[idx.high]

	n.censored = f$n.censor
	n.censored.low = n.censored[idx.low]
	n.censored.mid = n.censored[idx.mid]
	n.censored.high = n.censored[idx.high]

	tt = seq(0, 12, by=3)
	n.risks.show = NULL
	n.censored.show = NULL
	for( i in 1:length(tt) )
	{
		n.risk.high = n.risks.high[which(x.high-tt[i]>=0)[1]]
		n.risk.mid = n.risks.mid[which(x.mid-tt[i]>=0)[1]]
		n.risk.low = n.risks.low[which(x.low-tt[i]>=0)[1]]
		n.risks.show = cbind(n.risks.show, c(n.risk.high, n.risk.mid, n.risk.low))

		n.censor.high = sum(n.censored.high[which(x.high-tt[i]<0)])
		n.censor.mid = sum(n.censored.mid[which(x.mid-tt[i]<0)])
		n.censor.low = sum(n.censored.low[which(x.low-tt[i]<0)])
		n.censored.show = cbind(n.censored.show, c(n.censor.high, n.censor.mid, n.censor.low))
	}
	rownames(n.risks.show) = c("High risk", "Medium risk", "Low risk")
	colnames(n.risks.show) = tt
	rownames(n.censored.show) = c("High risk", "Medium risk", "Low risk")
	colnames(n.censored.show) = tt
	for(i in 1:3) n.risks.show[i,] = sapply(1:ncol(n.risks.show), function(j)paste0(n.risks.show[i,j], "(", n.censored.show[i,j], ")"))

	par(mar=c(10.5, 9.5-1.5, 3, 2))
	par(mgp=c(2.5, 1, 0))
	if( no.ylab ){
		y.lab = ""
	}else{
		y.lab = "Occurrence-free survival probability"
	}
	if( no.xlab ){
		x.lab = ""
	}else{
		x.lab = "Follow-up time (months)"
	}
	plot(c(0, x.low), c(1, y.low), type="s", lwd=lwd, col=col.low, xaxt="n", yaxt="n", xaxs="i", yaxs="i", axes=FALSE, 
		xlab=x.lab, ylab="", xlim=c(0,12), ylim=c(0, 1))
	mtext(y.lab, side=2, line=3)
	lines(c(0, x.mid), c(1, y.mid), lwd=lwd, col=col.mid, type='s')
	lines(c(0, x.high), c(1, y.high), lwd=lwd, col=col.high, type='s')
	axis(1, at=tt, lwd=lwd, lwd.ticks=lwd)
	axis(2, at=seq(0, 1, by=0.2), labels=c("0", seq(0.2,0.8,by=0.2), "1.0"), las=2, lwd=lwd, lwd.ticks=lwd)
	LEG = legend("bottomleft", legend=c("Low risk", "Medium risk", "High risk"), 
			col=c(col.low, col.mid, col.high), lty=1, lwd=lwd, bty='n', plot=FALSE)
	legend(x=LEG$rect$left, y=LEG$rect$top+0.225, legend=c("Low risk", "Medium risk", "High risk"), 
		col=c(col.low, col.mid, col.high), lty=1, lwd=lwd, bty='n')

	HR.text.1 = paste0("HR med vs low ", HR[1,1], " (95%CI ", HR[1,2], "), p", ifelse(HR[1,3]=="<0.0001", HR[1,3], paste0("=", HR[1,3])))
	HR.text.2 = paste0("HR high vs low ", HR[2,1], " (95%CI ", HR[2,2], "), p", ifelse(HR[2,3]=="<0.0001", HR[2,3], paste0("=", HR[2,3])))
	mtext(HR.text.1, side=1, line=-1.5*2, at=0.25, adj=0, font=1)
	mtext(HR.text.2, side=1, line=-1.5, at=0.25, adj=0, font=1)

	mtext(c("Number at risk (number censored)"), side=1, line=4.5, at=-5+0.35, adj=0, font=1)
	mtext(c("Low risk", "Medium risk", "High risk"), side=1, line=c(6, 7.5, 9), at=-1.25, adj=1)
	mtext(n.risks.show[3,], side=1, line=6, at=tt)
	mtext(n.risks.show[2,], side=1, line=7.5, at=tt)
	mtext(n.risks.show[1,], side=1, line=9, at=tt)
	mtext(num.label, side=2, line=1.3, at=1.1, adj=1, cex=1, las=2)
}


group_dat = function(dat, var.names, type="event")
{
	formula = as.formula(paste("Surv(time, event)~", paste(var.names, collapse="+")))
	cox = coxph(formula, data=dat, x=TRUE, y=TRUE) 
	risk.scores = predict(cox, type = "risk")
	dat$risk_score = risk.scores

	if( type=="event" ){
		dat <- dat[order(dat$risk_score), ]
		dat$cum_events <- cumsum(dat$event)
		total_events <- sum(dat$event)
		
		num_groups <- 3  # Define the number of groups
		events_per_group <- total_events / num_groups
		group <- 1
		event_count <- 0
		dat$group <- 1
		for (i in 1:nrow(dat)) {
		event_count <- event_count + dat$event[i]
		dat$group[i] <- group
		if (event_count >= events_per_group && group < num_groups) {
			group <- group + 1
			event_count <- 0
		}
		}
	}else{
		cutpoints <- surv_cutpoint(dat, time = "time", event = "event", variables = "risk_score")
		risk_group <- cut(dat$risk_score, 
					breaks = c(-Inf, cutpoints$cutpoint[1], cutpoints$cutpoint[2], Inf), 
					labels = c(1:3))
		dat$group = risk_group
	}

	return(dat)
}