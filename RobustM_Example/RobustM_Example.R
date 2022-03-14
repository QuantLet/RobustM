
##### SETTINGS #####

# Clean the environment 
graphics.off()
rm(list = ls(all = TRUE))

# Load Functions and other Files
source('./PackagesRobustMV.R')
source('./FunctionsRobustMV.R')

#Choose dataset to analyse
load("SP100.RData")

#### Efficient frontiers: seinsitivity for input parpameters ####
tickers = colnames(SP100$LogRets)
tickers_short = c('AAPL', 'IBM', 'NVDA', 'AMZN', 'MSFT')
start_date = "2015-01-01"
returns.data = SP100$LogRets[,match(tickers_short,tickers)]
end_date   = tail(index(returns.data), n = 1)
returns.data = returns.data[paste(start_date, end_date, sep="/"), ]

#Adjusted mean
returns.data1 = returns.data
returns.data1[ ,1] = returns.data1[ ,1]-0.0002

asset.names = tickers_short
er = colMeans(returns.data)
names(er) = asset.names
covmat = cov(returns.data)

er1 = colMeans(returns.data1)
names(er1) = asset.names
covmat1 = cov(returns.data1)

r.free = 0
dimnames(covmat) = list(asset.names, asset.names)
dimnames(covmat1) = list(asset.names, asset.names)

# tangency portfolio
tan.port = tangency.portfolio(er, covmat, r.free)
tan.port1 = tangency.portfolio(er1, covmat1, r.free)
# compute global minimum variance portfolio
gmin.port = globalMin.portfolio(er, covmat)
gmin.port1 = globalMin.portfolio(er1, covmat1)

# compute portfolio frontier
ef  = efficient.frontier(er, covmat, alpha.min=-2,
                         alpha.max=1.5, nport=20)
ef1 = efficient.frontier(er1, covmat1, alpha.min=-2,
                         alpha.max=1.5, nport=20)
#Ploting of the efficient frontiers
pngname =  paste0("figures/Efficient_frontiers.png")
png(file = pngname)
plot(ef)
plot(ef$sd, ef$er, col = "blue", pch = 16, xlab = "Standard Deviation", ylab="Mean")
points(ef1$sd,  ef$er, col = "red", pch = 16)
points(tan.port1$sd, tan.port1$er, col = "magenta", pch = 16, cex = 2)
points(tan.port$sd, tan.port$er, col = "green", pch = 16, cex = 2)
w_string  = toString(round(tan.port$weights,2))
w1_string = toString(round(tan.port1$weights, 2))
text(tan.port1$sd, tan.port1$er, labels= w_string, pos = 2)
text(tan.port$sd, tan.port$er, labels=w1_string, pos = 2)
dev.off()
