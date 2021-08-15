library(dplyr)
library(stringr)
library(reshape2)
library(purrr)
library(zoo)
#US GDP Deflator and Capital User Cost (Optional these deflators are not industry specific)
macro <- read.csv('/Users/justindoty/Documents/Research/Dissertation/Nonlinear_Production_Function_QR/Data/macro_vars.csv')[,-1]
#Data for other deflators from NBER
#NAICS
nbernaics <- read.csv('/Users/justindoty/Documents/Research/Dissertation/Nonlinear_Production_Function_QR/Data/naicsdef.csv', header=TRUE)
#SIC
nbersic <- read.csv('/Users/justindoty/Documents/Research/Dissertation/Nonlinear_Production_Function_QR/Data/sicdef.csv', header=TRUE)
#RnD Deflators from the BEA 
rnd <- read.csv('/Users/justindoty/Documents/Research/Dissertation/Nonlinear_Production_Function_QR/Data/rndbea.csv', header=TRUE)
names(rnd) <- c("year", "rprice")
rnd <- rnd %>% transmute(year=as.numeric(substr(year, 1, 4)), rprice=rprice/100)
#Average Sales per NAICS code (used for an SIC to NAICS conversion)
wsales <- read.csv('/Users/justindoty/Documents/Research/Dissertation/Nonlinear_Production_Function_QR/Data/compustat.csv') %>% select(gvkey, naics, sic, sale) %>% na.omit() %>%
	group_by(naics) %>% summarise(wsales=mean(sale))
wsales <- data.frame(wsales)
#SIC to NAICS mapping
SICto1997 <- read.csv('/Users/justindoty/Documents/Research/Dissertation/Nonlinear_Production_Function_QR/Data/1987_SIC_to_1997_NAICS.csv') %>% na.omit() 
names(SICto1997) <- c("naics", "SIC")
SICto1997 <- SICto1997 %>% left_join(wsales, by="naics") %>% select(SIC, naics, wsales) %>% arrange(SIC) %>% group_by(SIC) %>%
	summarise(sictonaics=naics[which.max(wsales)]) %>% rename(sic=SIC)
SICto1997 <- data.frame(SICto1997)
# Clean the NBER data a bit
#NAICS VERSION
naics3def <- nbernaics %>% select(naics, year, emp, pay, cap, invest, piinv, piship, pimat, vship, matcost, invest) %>% na.omit() %>%  mutate(Y=vship/piship, M=matcost/pimat, I=invest/piinv) %>% 
group_by(naics) %>% mutate(dep=ifelse(year==first(year), 0, (cap-(lag(invest/piinv)))/lag(cap)), avgpay=(pay)/(emp), naics3=str_extract(as.character(naics), "^.{3}")) %>% ungroup() %>%
group_by(naics3, year) %>% summarise(lprice=mean(avgpay), drate=mean(dep), yprice=sum(vship)/sum(Y), mprice=sum(matcost)/sum(M), iprice=sum(invest)/sum(I))
naics3def <- data.frame(naics3def)
names(naics3def) <- c("naics3", "year", "lprice3", "drate3", "yprice3", "mprice3", "iprice3")
#At the 2 digit level
naics2def <- nbernaics %>% select(naics, year, emp, pay, cap, invest, piinv, piship, pimat, vship, matcost, invest) %>% na.omit() %>%  mutate(Y=vship/piship, M=matcost/pimat, I=invest/piinv) %>% 
group_by(naics) %>% mutate(dep=ifelse(year==first(year), 0, (cap-(lag(invest/piinv)))/lag(cap)), avgpay=(pay)/(emp), naics2=str_extract(as.character(naics), "^.{2}")) %>% ungroup() %>%
group_by(naics2, year) %>% summarise(lprice=mean(avgpay), drate=mean(dep), yprice=sum(vship)/sum(Y), mprice=sum(matcost)/sum(M), iprice=sum(invest)/sum(I))
naics2def <- data.frame(naics2def)
names(naics2def) <- c("naics2", "year", "lprice2", "drate2", "yprice2", "mprice2", "iprice2")
#For capital
#I create a different database for capital delfators
kdef3 <- naics3def %>% select(naics3, year, iprice3) %>% rename(kyear=year, kprice3=iprice3)
kdef2 <- naics2def %>% select(naics2, year, iprice2) %>% rename(kyear=year, kprice2=iprice2)
# #SIC VERSION
# nberdef <- nbersic %>% select(sic, year, emp, pay, cap, invest, piinv, piship, pimat) %>% rename(sic4=sic) %>% mutate(sic4=as.character(sic4)) %>% na.omit() %>% group_by(sic4) %>% 
# mutate(dep=ifelse(year==first(year), 0, (cap-(lag(invest/piinv)))/lag(cap)), avgpay=pay/emp) %>% ungroup() %>%
# select(sic4, year, pay, dep, piship, pimat, piinv)
# sic4def <- data.frame(nberdef)
# names(sic4def) <- c("sic4", "year", "lprice4", "drate4", "yprice4", "mprice4", "iprice4")
# #At the 3 digit level
# sic3def <- sic4def %>% mutate(sic3=str_extract(sic4, "^.{3}")) %>% group_by(sic3, year) %>% summarise(lprice=mean(lprice4), drate=mean(drate4), yprice=mean(yprice4), mprice=mean(mprice4), iprice=mean(iprice4))
# sic3def <- data.frame(sic3def)
# names(sic3def) <- c("sic3", "year", "lprice3", "drate3", "yprice3", "mprice3", "iprice3")
# #At the 2 digit level
# sic2def <- sic3def %>% mutate(sic2=str_extract(sic3, "^.{2}")) %>% group_by(sic2, year) %>% summarise(lprice=mean(lprice3), drate=mean(drate3), yprice=mean(yprice3), mprice=mean(mprice3), iprice=mean(iprice3))
# sic2def <- data.frame(sic2def)
# names(sic2def) <- c("sic2", "year", "lprice2", "drate2", "yprice2", "mprice2", "iprice2")
# #At the 1 digit level
# sic1def <- sic2def %>% mutate(sic1=str_extract(sic2, "^.{1}")) %>% group_by(sic1, year) %>% summarise(lprice=mean(lprice2), drate=mean(drate2), yprice=mean(yprice2), mprice=mean(mprice2), iprice=mean(iprice2))
# sic1def <- data.frame(sic1def)
# names(sic1def) <- c("sic1", "year", "lprice1", "drate1", "yprice1", "mprice1", "iprice1")
# #For capital
# #I create a different database for capital delfators
# kdef4 <- sic4def %>% select(sic4, year, iprice4) %>% rename(kyear=year, kprice4=iprice4)
# kdef3 <- sic3def %>% select(sic3, year, iprice3) %>% rename(kyear=year, kprice3=iprice3)
# kdef2 <- sic2def %>% select(sic2, year, iprice2) %>% rename(kyear=year, kprice2=iprice2)
# kdef1 <- sic1def %>% select(sic1, year, iprice1) %>% rename(kyear=year, kprice1=iprice1)
#Cleaning and merging computstat data with price deflator data
# NAICS VERSION
compstat <- read.csv('/Users/justindoty/Documents/Research/Dissertation/Nonlinear_Production_Function_QR/Data/compustat.csv', header=TRUE) %>% rename(year=fyear, employ=emp) %>% 
	#Standard screening condition
	filter(indfmt=="INDL", consol=="C", datafmt=="STD", fic=="USA", curcd=="USD") %>%
	# # Use Historical NAICS when available, otherwise use most current NAICS code
 # 	mutate(naicshfill=ifelse(is.na(naicsh)&!is.na(naics), naics, naicsh)) %>% 
 # 	# For firms without current or historical NAICS code, covert SIC to NAICS
 # 	left_join(SICto1997, by="sic") %>%
 # 	mutate(sictonaics=ifelse(is.na(naicshfill), sictonaics, NA)) %>% 
 # 	# Create new NAICS codes with current/historical and converted SIC codes
 # 	mutate(naics=case_when(is.na(naicshfill)&!is.na(sictonaics)~sictonaics, !is.na(naicshfill)&is.na(sictonaics)~naicshfill)) %>%
	# # # Select Variables and Construct Age
	select(gvkey, year, sale, employ, ppegt, ppent, cogs, xsga, capx, naics, sic, fic, oibdp, dp, dpact) %>% group_by(gvkey) %>% mutate(age=year-first(year)+1) %>% ungroup() %>%
	#Manufacutring Industries
	filter(str_detect(naics, "^3")) %>% rename(id=gvkey) %>%
	#4-3-2 digit NAICS codes 
	mutate(naics3=ifelse(nchar(as.character(naics))>=3, substr(as.character(naics), 1, 3), NA), naics2=ifelse(nchar(as.character(naics))==2, substr(as.character(naics), 1, 2), NA)) %>%
	#Merge the deflator data
	left_join(naics3def, by=c("year", "naics3")) %>% left_join(naics2def, by=c("year", "naics2")) %>%
	#Coalesce the deflator data
	mutate(yprice=coalesce(yprice3, yprice2), lprice=coalesce(lprice3, lprice2), mprice=coalesce(mprice3, mprice2), iprice=coalesce(iprice3, iprice2)) %>%
	select(id, year, sale, employ, ppegt, ppent, cogs, xsga, capx, oibdp, dp, dpact, age, yprice, lprice, mprice, iprice, naics, sic, naics3, naics2) %>% filter(!is.na(dpact), !is.na(dp)) %>%
	#Calculate Average Age of Capital
	mutate(kage=dpact/dp) %>% group_by(id) %>% mutate(ksmooth=if(n()>=3) rollmean(kage, 3, fill=first(kage), align="right") else kage) %>% ungroup() %>% mutate(kyear=year-round(ksmooth)) %>%
	#Merge in capital deflators by average age of capital
	left_join(kdef3, by=c("kyear", "naics3")) %>% left_join(kdef2, by=c("kyear", "naics2")) %>%
	#Coalesce
	mutate(kprice=coalesce(kprice3, kprice2)) %>% inner_join(rnd, by="year") %>% inner_join(macro, by="year") %>% mutate(naics=coalesce(naics3, naics2)) %>%
	select(id, year, sale, employ, ppegt, ppent, cogs, xsga, capx, oibdp, age, yprice, lprice, mprice, iprice, kprice, rprice, USGDP, naics, sic) %>%
	na.omit() %>% mutate(Y=sale/yprice, K1=ppegt/kprice, K2=ppent/kprice, L=employ, I=capx/iprice, lexp=employ*lprice) %>% 
	filter(cogs>0, xsga>0) %>%
	mutate(mexp=(sale-oibdp)-lexp) %>% mutate(M=mexp/mprice) %>% mutate(VA=(oibdp+lexp)/yprice) %>% select(id, year, Y, VA, K1, K2, L, M, I, age, naics, sic, lexp, mexp) %>%
	filter(Y>0, K1>0, K2>0, L>0, M>0, I>0, VA>0) 
#SIC VERSION
# compstat <- read.csv('/Users/justindoty/Documents/Research/Dissertation/Nonlinear_Production_Function_QR/Data/compustat.csv', header=TRUE) %>% rename(year=fyear, employ=emp) %>% 
# 	#Standard screening condition
# 	filter(indfmt=="INDL", consol=="C", datafmt=="STD", fic=="USA", curcd=="USD") %>%
# 	# Select Variables and Construct Age
# 	select(gvkey, year, sale, employ, ppegt, ppent, cogs, xsga, capx, sic, fic, oibdp, dp, dpact, xrd) %>% group_by(gvkey) %>% mutate(age=year-first(year)+1) %>% ungroup() %>% rename(id=gvkey, rd=xrd) %>%
# 	#4-3-2-1 digit SIC codes 
# 	mutate(sic1=ifelse(sic%%1000==0, substr(as.character(sic), 1, 1), NA), sic2=ifelse(sic%%1000!=0&sic%%100==0, substr(as.character(sic), 1, 2), NA), sic3=ifelse(sic%%1000!=0&sic%%100!=0&sic%%10==0, substr(as.character(sic), 1, 3), NA), sic4=ifelse(sic%%10!=0, substr(as.character(sic), 1, 4), NA)) %>%
# 	#Merge the deflator data
# 	left_join(sic4def, by=c("year", "sic4")) %>% left_join(sic3def, by=c("year", "sic3")) %>% left_join(sic2def, by=c("year", "sic2")) %>% left_join(sic1def, by=c("year", "sic1")) %>%
# 	#Coalesce the deflator data
# 	mutate(yprice=coalesce(yprice4, yprice3, yprice2, yprice1), lprice=coalesce(lprice4, lprice3, lprice2, lprice1), mprice=coalesce(mprice4, mprice3, mprice2, mprice1), iprice=coalesce(iprice4, iprice3, iprice2, iprice1)) %>%
# 	select(id, year, sale, employ, ppegt, ppent, cogs, xsga, capx, oibdp, dp, dpact, rd, age, yprice, lprice, mprice, iprice, sic, sic4, sic3, sic2, sic1) %>%
# 	filter(!is.na(dpact), !is.na(dp)) %>%
# 	#Calculate Average Age of Capital
# 	mutate(kage=dpact/dp) %>% group_by(id) %>% mutate(ksmooth=if(n()>=3) rollmean(kage, 3, fill=kage) else kage) %>% ungroup() %>% mutate(kyear=year-round(ksmooth)) %>%
# 	#Merge in capital deflators by average age of capital
# 	left_join(kdef4, by=c("kyear", "sic4")) %>% left_join(kdef3, by=c("kyear", "sic3")) %>% left_join(kdef2, by=c("kyear", "sic2")) %>% left_join(kdef1, by=c("kyear", "sic1")) %>%
# 	#Coalesce
# 	mutate(kprice=coalesce(kprice4, kprice3, kprice2, kprice1)) %>% inner_join(rnd, by="year") %>% inner_join(macro, by="year") %>% mutate(USGDP=USGDP/100) %>%
# 	select(id, year, sale, employ, ppegt, ppent, cogs, xsga, capx, oibdp, dp, rd, age, yprice, USGDP, lprice, mprice, iprice, kprice, rprice, sic) %>%
# 	na.omit() %>% mutate(Y=sale*1e6/USGDP, K1=ppegt*1e6/kprice, K2=ppent*1e6/kprice, L=employ*1e3, I=capx*1e6/iprice, R=rd*1e6/rprice) %>% mutate(lexp=L*lprice) %>% 
# 	mutate(mexp=(sale-oibdp)*1e6-lexp) %>% mutate(M=mexp/USGDP) %>% mutate(VA=(oibdp+lexp)/yprice) %>% select(id, year, Y, VA, K1, K2, L, M, I, R, age, sic, lexp, mexp) %>%
# 	filter(Y>0, K1>0, K2>0, L>0, M>0, I>0, R>0) 
####################################################################################	
US <- compstat %>% mutate(naics2=str_extract(naics, "^.{2}")) %>% group_by(id) %>% filter(n()>=3)
US <- data.frame(US)
USest <- US %>% transmute(id=id, year=year, naics=naics, naics2=naics2, lny=log(Y), lnva=log(VA), lnk1=log(K1), lnk2=log(K2), lnl=log(L), lnm=log(M), lni=log(I), age=age)
#Rnd and Advertising
RAD <- read.csv('/Users/justindoty/Documents/Research/Dissertation/Nonlinear_Production_Function_QR/Data/compustat.csv', header=TRUE) %>% 
	#Standard screening condition
	filter(indfmt=="INDL", consol=="C", datafmt=="STD", fic=="USA", curcd=="USD") %>%
	select(gvkey, fyear, sale, xrd, xad) %>% rename(id=gvkey, year=fyear) %>% 
	inner_join(rnd, by="year") %>% inner_join(macro, by="year") %>%
	mutate(R=xrd/sale, RB=ifelse(xrd==0, 0, 1), ADV=xad/sale, ADVB=ifelse(xad==0, 0, 1)) %>% select(id, year, R, RB, ADV, ADVB)
USest <- USest %>% left_join(RAD, by=c("id", "year"))
####################################################################################################
#Summary statistics for the cleaned data set
print(summary(USest))
#Average number of firms per year
avg_firms <- round(as.numeric(colMeans(group_by(USest, year) %>% summarise(firms=n()))[2]))
print(avg_firms)
#Total number of firms in the dataset
unique_firms <- length(unique(USest$id))
print(unique_firms)
#Panel time length
panelT <- max(USest$year)-min(USest$year)+1
print(panelT)
#Average number of years a firm in the data
avg_time <- summary(group_by(USest, id) %>% summarise(time=n()))
print(avg_time)
#Observations per year
nsize <- as.data.frame(table(USest$year))$Freq
print(nsize)
#Total observations
print(nrow(USest))
#Save clean data
path_out <- '/Users/justindoty/Documents/Research/Dissertation/Nonlinear_Production_Function_QR/Data/'
# fileName <- paste(path_out, 'USdata.csv',sep = '')
# write.csv(USest, fileName, row.names=FALSE)
# USind <- USest 
USind <- USest %>% filter(year>=1997) %>% group_by(id) %>% filter(n()>=3)
#Prodest Test
require(prodest)
GO <- prodestLP(Y=USind$lny, fX=cbind(USind$lnl, USind$lnm), sX=USind$lnk1, pX=USind$lnm, idvar=USind$id, timevar=USind$year)
VA <- prodestLP(Y=USind$lnva, fX=USind$lnl, sX=USind$lnk1, pX=USind$lnm, idvar=USind$id, timevar=USind$year)
print(GO)
print(VA)

GOTFP <- USind$lny-cbind(USind$lnl, USind$lnm, USind$lnk1)%*%as.numeric(GO@Estimates$pars)-GO@Data$FSresiduals
GOTFP <- (GOTFP-mean(GOTFP))
