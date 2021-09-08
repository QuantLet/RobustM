##### SETTINGS #####
# Clean the environment 
graphics.off()
rm(list = ls(all = TRUE))

# Load Functions and other Files
source('./RobustM_Packages.R')
source('./RobustM_Functions.R')

# Working directory
print('Please use the code from the project in order to maintain a correct functionality of the code')

#Choose dataset to analyse
df = "Russell3000"# "SP100"
#load("SP100.RData")
load("Russell3000.RData")

#Choose the starting point of analysis  
start_point = "2010-01-01" #for the Russell3000 dataset
myData = switch(df,
                "SP100" = SP100,
                "Russell3000" = Russell3000)
insample = F
random   = T #if TRUE then it picks 600 random stocks for Russell3000 dataset

#### Python environment ####

use_virtualenv("r-reticulate")
path_to_python = ""#define the path to Pythonpython"
Sys.which("python")
Sys.setenv(RETICULATE_PYTHON = path_to_python)
use_python(path_to_python,  required = TRUE)
reticulate::py_config()


#Uncomment if set up the Python-R connection for the first time
# create a new environment if you run the code for the first time
#conda_create("r-reticulate")
#install Packages to the environment
# conda_install("r-reticulate", "scipy")
# conda_install("r-reticulate", "pandas")
# conda_install("r-reticulate", "matplotlib")
# conda_install("r-reticulate", "scikit-learn", pip = T)
# conda_install("r-reticulate", "statsmodels", pip = T)
# conda_install("r-reticulate", "seaborn")
# conda_install("r-reticulate", "cvxopt", pip = T)#from cvxopt import solvers, matrix, spmatrix

np     = import("numpy")
pd     = import("pandas")
mpl    = import("matplotlib")
plt    = import("matplotlib.pyplot")
cvxopt = import("cvxopt")
mrkw   = import("markowitz")

##Create Vectors of prices, returns 
if (df == "Russell3000" & random == T) { 
  random_stocks = sample(1:dim(myData$Data)[2], 600)
  rets_log    = myData$LogRets[, random_stocks]
  rets_log    = rets_log[which(index(rets_log) >= start_point),]
  prices      = as.matrix(myData$Data[which(index(rets_log) >= start_point), random_stocks])
  commonDateR = index(myData$Data)[which(index(myData$Data) >= start_point)]
  commonDate  = as.numeric(commonDateR)
}else {
  rets        = myData$Rets[, ]
  rets_log    = myData$LogRets[, ]
  commonDateR = index(myData$Data)
  prices      = as.matrix(myData$Data)
  commonDate  = as.numeric(commonDateR)
}

####Parameters for Portfolio Strategies####
ret                          = rets_log
freq                         = "months"
transi                       = 0.005 # transactional costs ratio
end                          = dim(ret)[1]
lambda                       = 2
hel                          = ret;
hel[]                        = 0
#gammas = 1.0*np$exp(-np$linspace(0, 4, as.integer(2000)))

####Portfolios constraction####
if (df == "Russell3000") {
     w_length = 2 # number of years of observtions to estimate input parameters
} else {
  w_length = 1
}
if (insample) {
  {first_signal = endpoints(commonDateR,on = "years")[1]}
} else {first_signal = endpoints(commonDateR,on = "years")[w_length + 1]}
ep = endpoints(commonDateR,on = freq) # rebalancing dates
ep = ep[(ep >= first_signal) & (ep < (end - 1))]
ret_data   = list()
ret_random = list()

if (df == "Russell3000") {
  for (t in commonDateR[ep]) {#building data for estimation of parameters
    random_stocks = sample(1:dim(myData$Data)[2], 600)
    t = as.POSIXlt(as.Date(t)) 
    start_date = t
    start_date$year = t$year - w_length # date of begining of the estimation window
    end_date = t
    end_date$year = t$year + 1
    ret_data[[as.character(as.Date(t))]] = ret[which(index(ret) > as.Date(start_date) & index(ret) <= as.Date(t)), ]#random_stocks
    ret_random[[as.character(as.Date(t))]] = ret[which(index(ret) > as.Date(t) & index(ret) <= as.Date(end_date)), ]#random_stocks
  }}else {
    for (t in commonDateR[ep]) {#building data for estimation of parameters
    t = as.POSIXlt(as.Date(t)) 
    start_date = t
    start_date$year = t$year - w_length # date of begining of the estimation window
    ret_data[[as.character(as.Date(t))]] = ret[which(index(ret) > as.Date(start_date) & index(ret) <= as.Date(t)),]
  }}
if (df == "Russell3000") {
  strats = c("EW",  "GMV_robust", "GMV_long", 
            "GMV_lin", "GMV_nlin", "GMV_mcd", "GMV_mve", "GMV_ogk", "GMV_nnve") 
} else {
  strats = c("EW",  "GMV_robust", "GMV_long", "GMV",
             "GMV_lin", "GMV_nlin", "GMV_mcd", "GMV_mve", "GMV_ogk", "GMV_nnve")
}

#Settings for parallel computing
rescoll     = list()
weightscoll = list()
UseCores  = detectCores() - 1
cl = parallel::makeCluster(UseCores, setup_strategy = "sequential")
registerDoParallel(cl)
ptm = proc.time()

results = foreach::foreach(stratloop = strats, .packages=libraries, .combine = cbind) %dopar% {
  weightcoll = matrix()
  retttscoll = list()
  ind =  which(stratloop %in% strats)

  end = dim(ret)[1]
  weight = hel;
  period = NA
  period = c(first_signal:(end - 1))
  if (insample) period = c(1:(end))
  for (t in period) {
    if (!is.na(match(t,ep))) { 
      retT = ret_data[[as.character(as.Date(commonDateR[t]))]]
      n = dim(retT)[2]
      rettt = retT
      w1 = NA
      mu = meanEstimation(na.omit(rettt))
      # Covariance estimation
      Sigma = covEstimation(na.omit(rettt))
      #for robust Markowitz
      block_n = as.integer(5)
      n = length(mu)
      if (qr(Sigma)$rank < n) Sigma = covEstimation(na.omit(rettt), control = list(type = 'cor'))
      
      if (stratloop == "EW")  w1 = rep(1/dim(rettt)[2], dim(rettt)[2])
      if (stratloop == "GMV")  w1 = optimalPortfolio(Sigma = Sigma,  
                                                     control = list(type = 'minvol', constraint = 'none'))
      if (stratloop == "GMV_long")  w1 = optimalPortfolio(Sigma = Sigma,  
                                                          control = list(type = 'minvol', constraint = 'lo'))
      if (stratloop == "GMV_robust") {
          mrkw   = import("markowitz") 
          gammas = mrkw$recommended_gamma(Sigma) 
          w1 = mrkw$robust_qp(rettt, block_n, gamma = gammas, lmbd = 0)}
      if (stratloop == "GMV_lin") {          
          Sigma_lin = linshrink_cov(na.omit(rettt), k = 0)
          w1 = optimalPortfolio(Sigma = Sigma_lin, 
                                control = list(type = 'minvol'))
        }
      if (stratloop == "GMV_nlin"){
          Sigma_nlin = nlshrink_cov(na.omit(rettt), k = 0)
          w1 = optimalPortfolio(Sigma = Sigma_nlin, 
                                control = list(type = 'minvol'))
        } 

      if (stratloop == "GMV_mcd") {Sigma_mcd = assetsMeanCov(na.omit(rettt), method = "mcd") 
          w1 = optimalPortfolio(Sigma = Sigma_mcd, mu = mu, control = list(type = 'minvol'))}
      if (stratloop == "GMV_mve")  {Sigma_mve = assetsMeanCov(na.omit(rettt), method = "mve")
          w1 = optimalPortfolio(Sigma = Sigma_mve, mu = mu, control = list(type = 'minvol'))}
      if (stratloop == "GMV_ogk")  {Sigma_ogk = assetsMeanCov(na.omit(rettt), method = "ogk")
          w1 = optimalPortfolio(Sigma = Sigma_ogk, mu = mu, control = list(type = 'minvol'))}
      if (stratloop == "GMV_nnve")  {Sigma_nnve = assetsMeanCov(na.omit(rettt), method = "nnve")
          w1 = optimalPortfolio(Sigma = Sigma_nnve, mu = mu, control = list(type = 'minvol'))}
          weight[t,] = w1 
    } else{    
      weight[t,] = weight[(t - 1),]
    }
  }
  res                                   = as.vector(row_sums(weight*ret)) 
  rescoll[[stratloop]]                  = res;
  weightscoll[[stratloop]]              = weight; 
  print(list(returns_str = res, weights = weight)) 
}
calculation_time = proc.time() - ptm
weightscoll                       = results[seq(2,length(results), 2)]
names(weightscoll)                = strats
rescoll                           = as.data.frame(results[seq(1,length(results), 2)])
rownames(rescoll)                 = index(weightscoll[[1]])
colnames(rescoll)                 = strats 
#### Compute performance####
#long period
collNumbers                            = NA
collNumbers_tc                         = NA
cp                                     = NA
weightsNumbers                         = NA
start                                  = first_signal
end                                    = dim(ret)[1] 
cp                                     = calculatePerformanceMeasures(start,end);
collNumbers                            = cp[[1]];
collres                                = cp[[2]]; 
weightsNumbers                         = cp[[3]];
collNumbers_tc                         = cp[[4]]; 
colnames(collNumbers)                  = strats; 
rownames(collNumbers)                  = c("CumWealth", "SR", "TTO", "TO")
collNumbers                            = t(collNumbers)
colnames(collNumbers_tc)               = strats; 
rownames(collNumbers_tc)               = c("CumWealth-TC", 
                                           "Sharpe.Ann-TC",  "TTO", "TO")
collNumbers_tc                         = t(collNumbers_tc)
colnames(collres)                      = strats 
n                                      = length(strats)
cols                                   = rainbow(length(strats), s = 1, v = 1, start = 0, end = max(1, n - 1)/n, alpha = 1)#c("red","blue","pink","darkgreen","lightgreen","lightblue","black")
weightsNumbers = t(weightsNumbers)
colnames(weightsNumbers) = c("min", "max", "sd", "mad-ew", "max-min");
rownames(weightsNumbers) = strats#[b][-4];

#### Visualization of performace####
if (df == "Russell3000") {
    ind = match(as.vector(c("EW","GMV_robust",  "GMV_long",
                           "GMV_lin", "GMV_nlin")), strats)
}else{
  ind = match(as.vector(c("EW",  "GMV_robust", "GMV_long", 
                          "GMV_lin", "GMV_nlin", "GMV")), strats) 
}
 
my_colors = RColorBrewer::brewer.pal(7,"Set1")[c(1,3:5,7)]
# EW - "#E41A1C" RColorBrewer::brewer.pal(7,"Set1")[1]
# GMV - "#377EB8" RColorBrewer::brewer.pal(7,"Set1")[2]
# GMV_lin - "#4DAF4A" RColorBrewer::brewer.pal(7,"Set1")[3]
# GMV_long - "#984EA3" RColorBrewer::brewer.pal(7,"Set1")[4]
# GMV_nlin - "#FF7F00" RColorBrewer::brewer.pal(7,"Set1")[5]
# GMV_robust - "#A65628" RColorBrewer::brewer.pal(7,"Set1")[7]
cumWealth = 1 + cumsum(collres[,ind])

plot13 = tidy(cumWealth) %>% ggplot(aes(x = index, y = value, color = series)) + 
  geom_line() + theme_bw() + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        legend.position = "none") +
  ylab("Cumulative Wealth") + xlab("Time") + 
  scale_color_manual(values = my_colors) 

pdfname =  paste0("CumWealth", df, ".pdf")
pngname =  paste0("CumWealth", df, ".png")
pdf(file = pdfname)
plot13
dev.off()
png(file = pngname)
plot13
dev.off()

#### Weights vizualization ####
GMV_lin = as.data.frame(weightscoll$GMV_lin[first_signal:end])
GMV_lin$Date = as.Date(index(weightscoll$GMV_lin[first_signal:end]))
weights_GMV_lin = gather(GMV_lin, Ticker, P.weight, 1:dim(rets)[2])
plot6 = ggplot(weights_GMV_lin, aes(x = Date, y = P.weight, fill = Ticker)) +
  geom_area() + theme_bw() + 
  ggtitle("GMV_lin") + theme(legend.position = "none", 
                             plot.title = element_text(color = RColorBrewer::brewer.pal(7,"Set1")[3])) +
  labs(x = "Year", y = "Weights") + ylim(-4, 4)

GMV_nlin = as.data.frame(weightscoll$GMV_nlin[first_signal:end])
GMV_nlin$Date = as.Date(index(weightscoll$GMV_nlin[first_signal:end]))
weights_GMV_nlin = gather(GMV_nlin, Ticker, P.weight, 1:dim(rets)[2])
plot7 = ggplot(weights_GMV_nlin, aes(x = Date, y = P.weight, fill = Ticker)) +
  geom_area() +
  theme_bw() + 
  ggtitle("GMV_nlin") + theme(legend.position = "none", 
                              plot.title = element_text(color = RColorBrewer::brewer.pal(7,"Set1")[5])) +
  labs(x = "Year", y = "Weights") + ylim(-4, 4)

GMV_long = as.data.frame(weightscoll$GMV_long[first_signal:end])
GMV_long$Date = as.Date(index(weightscoll$GMV_long[first_signal:end]))
weights_GMV_long = gather(GMV_long, Ticker, P.weight, 1:dim(rets)[2])
plot9 = ggplot(weights_GMV_long, aes(x = Date, y = P.weight, fill = Ticker)) +
  geom_area() +
  theme_bw() + ggtitle("GMV_long") +
  theme(legend.position = "none", plot.title = element_text(color = RColorBrewer::brewer.pal(7,"Set1")[4])) +
  labs(x = "Year", y = "Weights") + ylim(-4, 4)

GMV_robust = as.data.frame(weightscoll$GMV_robust[first_signal:end])
GMV_robust$Date = as.Date(index(weightscoll$GMV_robust[first_signal:end]))
weights_GMV_robust = gather(GMV_robust, Ticker, P.weight, 1:dim(rets)[2])
plot10 = ggplot(weights_GMV_robust, aes(x = Date, y = P.weight, fill = Ticker)) +
  geom_area() +
  theme_bw() + ggtitle("GMV_robust") + 
  theme(legend.position = "none", plot.title = element_text(color = RColorBrewer::brewer.pal(7,"Set1")[7])) +
  labs(x = "Year", y = "Weights") + ylim(-4, 4)

EW = as.data.frame(weightscoll$EW[first_signal:end])
EW$Date = as.Date(index(weightscoll$EW[first_signal:end]))
weights_EW = gather(EW, Ticker, P.weight, 1:dim(rets)[2])
plot11 = ggplot(weights_EW, aes(x = Date, y = P.weight, fill = Ticker)) +
  geom_area() +   ggtitle("EW") +
  theme_bw() + theme(legend.position = "none", plot.title = element_text(color = RColorBrewer::brewer.pal(7,"Set1")[1])) +
  labs(x = "Year", y = "Weights") + ylim(-4, 4)

GMV = as.data.frame(weightscoll$GMV[first_signal:end] )
GMV$Date = as.Date(index(weightscoll$GMV[first_signal:end] ))
weights_GMV = gather(GMV, Ticker, P.weight, 1:dim(rets)[2])
plot12 = ggplot(weights_GMV, aes(x = Date, y = P.weight, fill = Ticker)) +
  geom_area() +  ggtitle("GMV") + theme_bw() +
  theme(legend.position = "none", plot.title = element_text(color = RColorBrewer::brewer.pal(7,"Set1")[2])) +
  labs(x = "Year", y = "Weights") + ylim(-4, 4)

pdfname5 = paste0("Weights_GMV_all",df,".pdf")
pngname5 = paste0("Weights_GMV_all",df,".png")
if (df == "SP100") {
  g = arrangeGrob(plot11, plot10, plot9, plot12, plot6, plot7,
                  ncol = 2, nrow = 3) 
  ggsave(file = pdfname5, g)
  ggsave(file = pngname5, g)
}else{
  g = arrangeGrob(plot10, plot9, plot6, plot7,
                  ncol = 2, nrow = 2)
  ggsave(file = pdfname5, g)
  ggsave(file = pngname5, g)
}
dev.off()

#### Create latex table ####
performance.table.latex = xtable(collNumbers[ind,], digits = 4)
print(x = performance.table.latex, file = paste0("Performance_",  df, ".tex") ,
      include.rownames = T, booktabs = T, floating = F)

performance_tc.table.latex = xtable(collNumbers_tc[ind,], digits = 4)
print(x = performance_tc.table.latex, file = paste0("Performance_tc_", df, ".tex") ,
      include.rownames = T, booktabs = T, floating = F)

weights.table.latex = xtable(weightsNumbers[ind,], digits = 4)
print(x = weights.table.latex, file = paste0("Weights_", df, ".tex") ,
      include.rownames = T, booktabs = T, floating = F)

####Save  Results ####

save.image(file = paste0(Sys.Date(),"Results_RobustM",  df,  ".RData"))  
#Tracking of packages used versions#
writeLines(capture.output(sessionInfo()), paste("sessionInfo.txt"))



