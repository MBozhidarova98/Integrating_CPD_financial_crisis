##################### Do for the 5 highest dependent variables #######################
negative_log_returns <- read.csv("negative_log_returns.csv")
library(evd)
library(extRemes)
library(ReIns)
library(tea)
library(evir)

#Collect all data before 2008:
start_date='2002-01-17' #index 3808
end_date='2007-06-30' #index 3208

#Collect all data before 2020:
start_date='2002-01-17' #index 619
end_date='2020-02-29' #index 257

negative_log_returns=negative_log_returns[negative_log_returns$Date >= start_date & negative_log_returns$Date <= end_date, ]

#remove the date column:
negative_log_returns=negative_log_returns[, !(colnames(negative_log_returns) %in% c('Date','X'))]


is=c()
js=c()
chis=c()
comb=combn(colnames(negative_log_returns), 2, simplify = FALSE)
for(i in 1:length(comb)){
  x=negative_log_returns[,comb[[i]][1]]
  y=negative_log_returns[,comb[[i]][2]]
  df=data.frame(x,y)
  df=na.omit(df)
  if(dim(df)[1]>500){
    t=taildep(x,y,0.95,type=c('chi'),na.rm=TRUE)
    is=c(is,comb[[i]][1])
    js=c(js,comb[[i]][2])
    chis=c(chis,t[[1]])
  }
}  
df=data.frame(is,js,chis)
df=df[order(df$chis),]
df<-df[dim(df)[1]:1,]
df<-na.omit(df)
write.csv(x=df, file='chis_95.csv')


