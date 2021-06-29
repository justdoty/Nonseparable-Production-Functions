library(dplyr)
library(stringr)
library(reshape2)
library(purrr)
#US GDP Deflator and Capital User Cost
macro <- read.csv('/Users/justindoty/Documents/Research/Dissertation/Nonlinear_Production_Function_QR/Data/macro_vars.csv')[,-1]
#Data for other deflators from NBER
nber <- read.csv('/Users/justindoty/Documents/Research/Dissertation/Nonlinear_Production_Function_QR/Data/naicsdef.csv', header=TRUE)
#Clean the NBER data a bit
#I keep emp and pay to construct industry level average wage bill since XLR (staff expense) is often missing in compustat
#EMP is measured in thousdands
#PAY is measured in millions
#I keep cap, invest, and piinv to construct industry level average depreciation rates on capital
nberdef <- nber %>% select(naics, year, emp, pay, cap, invest, piinv, piship, pimat) %>% na.omit() %>% group_by(naics) %>% 
mutate(dep=ifelse(year==first(year), 0, (cap-(lag(invest/piinv)))/lag(cap)), avgpay=(pay)/(emp), naics3=str_extract(as.character(naics), "^.{3}")) %>%
group_by(naics3, year) %>% summarise(lprice=mean(avgpay), drate=mean(dep), yprice=mean(piship), mprice=mean(pimat), iprice=mean(piinv))
#I create a different database for capital delfators
kdef <- nber %>% select(naics, year, piinv) %>% rename(kyear=year) %>% na.omit() %>% group_by(naics) %>% 
mutate(naics3=str_extract(as.character(naics), "^.{3}")) %>%
group_by(naics3, kyear) %>% summarise(kprice=mean(piinv))
#Cleaning and merging computstat data with price deflator data
#This is the main dataset, it includes all variables needed for production function estimation using OP, LP, ACF, or GNR
compstat <- read.csv('/Users/justindoty/Documents/Research/Dissertation/Nonlinear_Production_Function_QR/Data/compustat.csv', header=TRUE) %>% rename(year=fyear, employ=emp) %>% 
	select(gvkey, year, sale, employ, ppegt, ppent, cogs, xsga, capx, naics, fic, oibdp, dp, dpact, ipodate) %>% group_by(gvkey) %>% mutate(age=year-first(year)+1) %>% ungroup() %>%
	filter(str_detect(naics, "^31|^32|^33"), fic=="USA") %>%
	transmute(id=gvkey, year=year, sale=sale*1e3, oibdp=oibdp*1e3, cogs=cogs*1e3, xsga=xsga*1e3, ppegt=ppegt*1e3, ppent=ppent*1e3, employ=employ*1e3, capx=capx*1e3, dpact=dpact*1e3, dp=dp*1e3, naics3=str_extract(as.character(naics), "^.{3}"), age=age, ipoyear=as.numeric(str_extract(as.character(ipodate), "^.{4}"))) %>%
	#Merge with US GDP deflator and NBER Data
	inner_join(macro, "year") %>% inner_join(nberdef, c("naics3", "year")) %>% mutate(lexp=employ*lprice, kage=dpact/dp) %>% mutate(mexp=cogs+xsga-dp-lexp, kyear=year-round(kage)) %>%
	mutate(id=id, year=year, Y=(sale/USGDP)*100, K=(ppegt/USGDP)*100, K2=(ppent/USGDP)*100, M=(mexp/USGDP)*100, I=(capx/USGDP)*100, L=employ) %>% 
	mutate(VA=((sale-mexp)/USGDP)*100, S=mexp/sale) %>% 
	group_by(id) %>% na.omit() %>% filter(Y>0, K>0, K2>0, M>0, I>0, L>0, VA>0)
compstat <- compstat %>% select(id, year, Y, K, K2, M, I, L, VA, S, ipoyear, age, naics3, drate) %>% filter(year>=1990, year<=2016) %>% group_by(id) %>% filter(n()>=4)
####################################################################################################
#Summary statistics for the cleaned data set
print(summary(compstat))
#Average number of firms per year
avg_firms <- round(as.numeric(colMeans(group_by(compstat, year) %>% summarise(firms=n()))[2]))
print(avg_firms)
#Total number of firms in the dataset
unique_firms <- length(unique(compstat$id))
print(unique_firms)
#Panel time length
panelT <- max(compstat$year)-min(compstat$year)+1
print(panelT)
#Average number of years a firm in the data
avg_time <- summary(group_by(compstat, id) %>% summarise(time=n()))
print(avg_time)
#Total observations
print(nrow(compstat))
#Save clean data
path_out <- '/Users/justindoty/Documents/Research/Dissertation/Nonlinear_Production_Function_QR/Data/'
fileName <- paste(path_out, 'USdata.csv',sep = '')
write.csv(compstat,fileName, row.names=FALSE)

