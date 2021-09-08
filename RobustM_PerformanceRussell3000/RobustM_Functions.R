
pturnoverDN  = function(weights, rets, freq){
  results = Return.portfolio(R = rets,   
                             weights = weights, 
                             rebalance_on = freq, verbose = T)
  bop = results$BOP.Weight #beginning of period weights
  bop
  eop = results$EOP.Weight #end of period weights
  eop
  f = abs(bop - eop)
  out = sum(f)*(1/(length(ep) - 1)) #  
  return(out)
}

mad_ew = function(x){
  a = abs(x - 1/ncol(ret))
  return(a)
}
trans_cost = function(weights, rets, freq, c){
  results = Return.portfolio(R = rets,   
                             weights = weights, 
                             rebalance_on = freq, verbose = T)
  
  bop = results$BOP.Weight #beginning of period weights
  bop
  eop = results$EOP.Weight #end of period weights
  eop
  out = c*row_sums(abs(eop - bop))
  return(out)
} 

calculatePerformanceMeasures = function(start,end){
  collNumbers = vector()
  collNumbers_tc = vector()
  collres = xts()
  weightsNumbers = vector()
  for (stratloop in 1:length(strats)){
    Rb = as.xts(rescoll[,  which(strats %in% c("EqualWeight","EW"))], 
                order.by = as.Date(rownames(rescoll)))
    portfolioret_net = na.omit(rescoll[,stratloop])
    strat_weights = weightscoll[[stratloop]] 
    strat_weights[is.nan(strat_weights)] = 0.0
    portfolioret_net_xts = as.xts(as.matrix(na.omit(rescoll[,stratloop])), 
                                  order.by = as.Date(na.omit(rownames(rescoll))))
    portfolioEquity_net = 1 + cumsum(portfolioret_net)
    cumWealth = tail(portfolioEquity_net, 1)
    firstsignal = start
    rettt = portfolioret_net[firstsignal:end]
    rettt_xts = portfolioret_net_xts[firstsignal:end]
    ret_data = rets_log
    stock_rets = ret_data[firstsignal:end]
    tc = trans_cost(strat_weights[firstsignal:end,], stock_rets, freq, c = transi)
    portfolioret_net_tc = portfolioret_net[firstsignal:(end - 1)] - tc
    portfolioEquity_net_tc = 1 + cumsum(portfolioret_net_tc)
    cumWealth_tc = tail(portfolioEquity_net_tc, 1)
    T = (commonDate[end] - commonDate[firstsignal])/365
    Return.ann = (portfolioEquity_net[end]/portfolioEquity_net[firstsignal - 1])^(1/T) - 1
    Return.ann_tc = (tail(portfolioEquity_net_tc, 1)/portfolioEquity_net_tc[firstsignal - 1])^(1/T) - 1
    Vola.ann = sd(rettt)*sqrt(252);
    Vola.ann_tc = sd(portfolioret_net_tc)*sqrt(252);
    Sharpe.ann = Return.ann/Vola.ann
    Sharpe.ann_tc = Return.ann_tc/Vola.ann_tc
    target_turnover = vector();
    for (i in 2:dim(strat_weights)[1]) {
      target_turnover[i] = sum(abs(matrix(strat_weights[i, ]) - matrix(strat_weights[i - 1,])))/dim(strat_weights)[2] 
    }
    Turnover = mean(na.omit(target_turnover))
    value = portfolioEquity_net 
    ma = unlist(lapply(c(2:length(value)),function(x) max(value[1:x])))
    dddisc = value[-1]/ma - 1
    datums = commonDateR[firstsignal:end]
    num = as.numeric(tail(datums,1)-datums[1])
    PR = as.vector(PainRatio(rettt_xts))
    TurnoverDM  = pturnoverDN(strat_weights[firstsignal:end,], stock_rets, freq)
    Return_annual = as.vector(Return.annualized(rettt_xts, geometric = F))
    AverageDrawdown = as.numeric(AverageDrawdown(rettt_xts))
    Sharpe =  as.numeric(SharpeRatio(rettt_xts))[1]
    StdDev.annualized = as.numeric(StdDev.annualized(rettt_xts))
    collNumbers = cbind(collNumbers, as.vector(c(cumWealth, 100*Sharpe.ann,  
                                                 Turnover, TurnoverDM)))#, 
    collNumbers_tc = cbind(collNumbers_tc, as.vector(c(cumWealth_tc,  
                                                       100*Sharpe.ann_tc, Turnover, TurnoverDM)))
    weightcoll = as.data.frame(strat_weights[first_signal:end])
    weightsNumbers = cbind(weightsNumbers, as.vector(c((mean(apply(weightcoll, 1, min))), mean(apply(weightcoll, 1, max)), 
      mean(apply(weightcoll, 1, sd)), mean(apply(weightcoll, 1, mad_ew)), mean(diff(apply(weightcoll, 1, range))))))
    collNumbers = round(collNumbers,4)
    weightsNumbers = round(weightsNumbers,4)
    res = as.xts(portfolioret_net[first_signal:end],order.by=index(ret)[first_signal:end])
    collres = cbind(collres,res)
    first(index(res))
    last(index(res))
    
  }
  return(list(collNumbers,collres,weightsNumbers, collNumbers_tc))
}





