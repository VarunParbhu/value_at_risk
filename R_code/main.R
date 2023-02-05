library(vctrs)
library(dplyr)
library(tidyquant)
library(quantmod)
library(PerformanceAnalytics)
library(PortfolioAnalytics)
library(ROI)
library(car)
library(psych)
library(geometry)
library(ggplot2)
library(reshape2)
library(xts)
library(tseries)
library(zoo)
library(FinTS)
library(rugarch)

options("getSymbols.warning4.0"=FALSE)
options("getSymbols.yahoo.warning"=FALSE)

tickers = c("SAP", "ORCL", "WDAY", "CDAY", "CACI", "ADP", "PCTY", "PAYC")

#weights = c(0.125, 0.125 ,0.125 , 0.125 , 0.125 , 0.125 , 0.125 ,0.125)
# SAP  Sap SuccessFactors
# ORCL Oracle HCM
# WDAY Workday HCM
# CDAY Ceridian Dayforce HCM
# ADP  Automatic Data Processing
# CACI CACI International, Inc.
# PCTY Paylocity Holding Corp
# PAYC Paycom

#Getting Data
getSymbols(tickers)

#NA Present
# Market Analysis
###Figure 1 : Stock Price of 8 HCM Company [January 2007 - June 2021]
Port.prices.NA = merge(Ad(SAP),Ad(ORCL),Ad(WDAY),Ad(CDAY),Ad(CACI),Ad(ADP),Ad(PCTY),Ad(PAYC))
plot(Port.prices.NA, lwd =1, main = "Stock Price ($)",major.ticks= "years",
     minor.ticks = 'years')
addLegend("topleft", on=1, 
          legend.names = c("SAP", "ORCL", "WDAY", "CDAY", "CACI", "ADP", "PCTY", "PAYC"), col=c(),lty=c(1, 1),cex=0.75)

###Figure 2 : Stock Prices of 8 HCM Companies [January 2017 - June 2021]
Port.prices.NA = merge(Ad(SAP),Ad(ORCL),Ad(WDAY),Ad(CDAY),Ad(CACI),Ad(ADP),Ad(PCTY),Ad(PAYC))
plot(Port.prices.NA["20170101/20210601"], lwd =1, main = "Stock Price ($)",major.ticks= "years",
        minor.ticks = 'years')
addLegend("topleft", on=1, 
          legend.names = c("SAP", "ORCL", "WDAY", "CDAY", "CACI", "ADP", "PCTY", "PAYC"), col=c(),lty=c(1, 1),cex=0.75)

#NA Deleted and Log Returns
Port.prices = na.omit(merge(Ad(SAP),Ad(ORCL),Ad(WDAY),Ad(CDAY),Ad(CACI),Ad(ADP),Ad(PCTY),Ad(PAYC)))
port.logreturns = diff(log(Port.prices))[-1]
colnames(port.logreturns) <- tickers

chart.Correlation(port.logreturns) #Correlation among stock value Chart

Port.returns = port.logreturns #Using LogReturns

#Simple Portfolio Optimization to get optimal weights
portf <- portfolio.spec(colnames(Port.returns))

##CONSTRAINTS
### No Shorts is allowed (#100% of Initial value is invested)
portf <- add.constraint(portf,type="weight_sum", min_sum = 1, max_sum = 1) 
### 0% can be invested in a stock; not more than 50% should be invested in a single stock
portf <- add.constraint(portf, type='box',min=0,max = 0.5)
### Maximize the Return
portf <- add.objective(portf, type="return", name="mean")
### Minimize the Variance
portf <- add.objective(portf, type="risk", name="StdDev")
### Optimizing Portfolio using ROI (simple optimization)
optPort <- optimize.portfolio(Port.returns, portf, optimize_method = "ROI")

#Optimal Weights
#SAP           ORCL          WDAY          CDAY           CACI           ADP            PCTY          PAYC 
#7.920902e-17  3.498359e-01  2.815557e-17  2.25029e-17    1.519507e-01   7.297861e-01   4.252348e-01  2.293632e-18
#0%            35%           0%            0%             15.2%          7.3%           42.5%         0%  


#Histogram of Historical Portfolio Returns using the Optimal Weights
opt.rtrns = Port.returns%*%optPort$weights
hist(opt.rtrns, breaks = 100, main = "Historical Portfolio Log Returns",xlab="Log Returns")
curve(dnorm(x, mean=mean(opt.rtrns), sd=sd(opt.rtrns)), add=TRUE, col='red')

### CONFIDENCE INTERVAL p
p=0.99

## Value at Risk##
#Individual VaR at 99%
  VAR.hist = VaR(Port.returns,p=p,weights = NULL, portfolio_method = "single",method = "historical")
  VAR.gaus = VaR(Port.returns,p=p,weights = NULL, portfolio_method = "single",method = "gaussian")
  VAR.modi = VaR(Port.returns,p=p,weights = NULL, portfolio_method = "single",method = "modified")
  All.VAR = data.frame(rbind(VAR.hist,VAR.gaus,VAR.modi))
  rownames(All.VAR) <- c("Hist","Gauss","Modified")
  
  
#Adding PortFolio VaR at 99% in data frame
  PortVAR.hist=VaR(Port.returns,p=p,weights = optPort$weights, portfolio_method = "component",method = "historical")$hVaR
  PortVAR.gaus=VaR(Port.returns,p=p,weights = optPort$weights, portfolio_method = "component",method = "gaussian")$VaR
  PortVAR.modi=VaR(Port.returns,p=p,weights = optPort$weights, portfolio_method = "component",method = "modified")$MVaR
  All.VAR$Portfolio <- c(PortVAR.hist,PortVAR.gaus,PortVAR.modi)
  All.VAR = -1*abs(All.VAR)
                         
All.VAR = (exp(All.VAR)-1)*10000

## Conditional Value at Risk##
#Individual CVaR at 99%
  CVAR.hist = CVaR(Port.returns,p=p,weights = NULL, portfolio_method = "single",method = "historical")
  CVAR.gaus = CVaR(Port.returns,p=p,weights = NULL, portfolio_method = "single",method = "gaussian")
  CVAR.modi = CVaR(Port.returns,p=p,weights = NULL, portfolio_method = "single",method = "modified")
  All.CVAR <- data.frame(rbind(CVAR.hist,CVAR.gaus,CVAR.modi))
  rownames(All.CVAR) <- c("Hist","Gauss","Modified")

  
#Adding PortFolio CVaR at 99% in data frame
  CPortVAR.hist= CVaR(Port.returns,p=p,weights = optPort$weights, portfolio_method = "component",method = "historical")[1]
  CPortVAR.gaus= CVaR(Port.returns,p=p,weights = optPort$weights, portfolio_method = "component",method = "gaussian")[1]
  CPortVAR.modi= CVaR(Port.returns,p=p,weights = optPort$weights, portfolio_method = "component",method = "modified")[1]
  All.CVAR$Portfolio <- c(CPortVAR.hist,CPortVAR.gaus,CPortVAR.modi)
  All.CVAR$Portfolio = -1*(as.numeric(All.CVAR$Portfolio))
  
All.CVAR = (exp(All.CVAR)-1)*10000  

#  Monte Carlo Method
## Multivariate Normal Distribution

OptimalWeights = matrix(optPort$weights)
  
COVM.Port = cov(Port.returns) # Covariance Matrix of Port Folio
MEAN.PortAssets = matrix(colMeans(Port.returns)) # Mean of Individual Stocks
STD.PortAssets =sqrt(diag(COVM.Port)) #STD of each stock
PortCorrelation = cor(Port.returns) #PortFolio Correlation

STD.Weighted = OptimalWeights * STD.PortAssets
Portfolio.std = sqrt(t(OptimalWeights) %*% (COVM.Port %*% OptimalWeights)) #Portfolio Standard Deviation

### Cholesky Decomposition
ColCOVM.port = matrix(chol(COVM.Port),nrow=8,ncol=8)
ColTCOVM.port = t(matrix(chol(COVM.Port),nrow=8,ncol=8))

### Setting up MC for 10000 Simulation for 30 Trading Days
mc_sims = 10000
days = 30

Mean.Assets = matrix( MEAN.PortAssets , length(MEAN.PortAssets) , 30 )
Varianceby2 = matrix(0.5*diag(COVM.Port) , length(0.5*diag(COVM.Port)) , 30 )

simulation = data.frame(Days = 0:30)
InitialPortfolio = 10000

seed = 1
for (i in 1:mc_sims){
                      DailyLogReturns = ColTCOVM.port %*% matrix(rnorm(8*30, mean=0, sd=1),8,30) + Mean.Assets - Varianceby2
                      
                      sim = t(OptimalWeights) %*% DailyLogReturns 
                      sim2  = rbind(rep(0,1), matrix(cumsum(sim)))
                      simulation[i+1] = InitialPortfolio*exp(sim2)
}

#Graphing of first 250 Simulation

YY = simulation [2,]
YYR = t(tail(YY,n=1)[-1]) - 10000
YYRS = sort(YYR, decreasing = TRUE)
MVAR = YYRS[9900]
MCVAR = mean(YYRS[9500:length(PortfolioResults_sorted)])


simulation.m <- as.data.frame.array(melt(simulation[1:250], id.var="Days"))
ggplot(simulation.m, aes(x =Days, y=value,colour = variable)) +
  geom_line() +theme(legend.position = "none") + xlab("Days") +
  ylab("Portfolio Value") + 
  ggtitle("First 250 Simulations") + 
  scale_color_hue(h = c(0,200), l = 70, c = 50) 


# 1 Day VaR MVN
DAY_1 = simulation [2,]
DAY_1_Return = t(tail(DAY_1,n=1)[-1]) - 10000
DAY_1_Return_Sort = sort(DAY_1_Return, decreasing = TRUE)
MVAR_MVT = DAY_1_Return_Sort[9900]
MCVAR_MVT = mean(YYRS[9900:length(PortfolioResults_sorted)])




#VaR and CVaR Monte Carlo Simulation
PortfolioResults = t(tail(simulation,n=1)[-1]) - 10000
PortfolioResults_sorted = sort(PortfolioResults, decreasing = TRUE)

## 99% Confidence Position = 9900 at 30th day
MVAR = PortfolioResults_sorted[9900]
MCVAR = mean(PortfolioResults_sorted[9900:length(PortfolioResults_sorted)])

#VAR of Portfolio from Historical Data
## 30 Day Portfolio VAR

var30 = Portfolio.std*InitialPortfolio*sqrt(30)*2.33 #High Correlation; high standard Deviation
hist(PortfolioResults/10000, breaks=20)
##### -> Calculate future VAR from the existing one


####GARCH Simulation
Portfolio_series.return = Port.returns%*%optPort$weights
date = index(Port.returns)
Port.timeseries = xts(x = Portfolio_series.return, order.by = date)

par(mfrow=c(2,2),xpd=FALSE)
Acf(Portfolio_series.return,main="ACF of Portfolio Return",lag=100)
pacf(Portfolio_series.return, main="PACF of Portfolio Return",lag=100)
Acf(((Portfolio_series.return - mean(Portfolio_series.return))^2),main = "ACF of Squared Portfolio Error",lag=100)
pacf((Portfolio_series.return - mean(Portfolio_series.return))^2,main ="PACF of Squared Portfolio Error",lag=100)

#Fitting the GARCH(1,1) Model using rugarch package
##Model Specification
x = ugarchspec(variance.model = list(garchOrder=c(1,1)),mean.model = list(armaOrder=c(0,0)))
##Fitting Model Specification
x_fit = ugarchfit(x,data=Port.timeseries)

#Setting up Monte Carlo Simulation for Fitted GARCH(1,1) for 30 Trading Days
xfinal = x
setfixed(xfinal)<-as.list(coef(x_fit))

sim_G = ugarchpath(spec = xfinal,m.sim=10000,n.sim=1*30,rseed=1)

p_G = data.frame(Days = 0:30)
p_G = cbind(p_G,rbind(T=10000,InitialPortfolio*apply(fitted(sim_G),2,'cumsum')+ InitialPortfolio))


simulation.G <- as.data.frame.array(melt(p_G[1:250], id.var="Days"))
ggplot(simulation.G, aes(x =Days, y=value,colour = variable)) +
  geom_line() +theme(legend.position = "none") + xlab("Days") +
  ylab("Portfolio Value") + 
  ggtitle("First 250 Simulations") + 
  scale_color_hue(h = c(0,200), l = 70, c = 50) 


#1 Day VaR GARCH
DAY_1_G = p_G[2,]
DAY_1_Return_G = t(tail(DAY_1_G,n=1)[-1]) - 10000
DAY_1_Return_Sort_G = sort(DAY_1_Return_G, decreasing = TRUE)
MVAR_G1 = DAY_1_Return_Sort_G[9900]
MCVAR_G1 = mean(DAY_1_Return_Sort_G[9900:length(PortfolioResults_sorted)])


#VaR and CVaR Monte Carlo Simulation
PortfolioResults_G = t(tail(p,n=1)[-1]) - 10000
PortfolioResults_sorted_G = sort(PortfolioResults_G, decreasing = TRUE)

## 95% Confidence Position = 9500 at 30th day
MVAR_G = PortfolioResults_sorted_G[9900]
MCVAR_G = mean(PortfolioResults_sorted_G[9900:length(PortfolioResults_sorted)])
