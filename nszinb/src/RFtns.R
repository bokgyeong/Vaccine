
# Randomized quantile residuals (RQRs)
stdRQR_ZINB = function(y, lambda, theta, pii){
  N = length(y)
  u = runif(N)
  pzinb = dzinb = numeric(N)
  indexzero = which(y == 0)
  indexpos = which(y > 0)
  pzinb[indexpos] = (1 - pii[indexpos]) + pii[indexpos] * pnbinom(y[indexpos]-1, mu = lambda[indexpos], size = theta[indexpos])
  
  dzinb[indexzero] = (1 - pii[indexzero]) + pii[indexzero] * dnbinom(0, mu = lambda[indexzero], size = theta[indexzero])
  dzinb[indexpos] = pii[indexpos] * dnbinom(y[indexpos], mu = lambda[indexpos], size = theta[indexpos])
  RQR = qnorm(pzinb + u * dzinb)
  result = (RQR - mean(RQR[!(is.infinite(RQR) | is.nan(RQR))])) / sd(RQR[!(is.infinite(RQR) | is.nan(RQR))])
  return(result)
}


stdRQR_ZIP = function(y, lambda, pii){
  N = length(y)
  u = runif(N)
  pzip = dzip = numeric(N)
  indexzero = which(y == 0)
  indexpos = which(y > 0)
  pzip[indexpos] = (1 - pii[indexpos]) + pii[indexpos] * ppois(y[indexpos]-1, lambda = lambda[indexpos])
  
  dzip[indexzero] = (1 - pii[indexzero]) + pii[indexzero] * dpois(0, lambda = lambda[indexzero])
  dzip[indexpos] = pii[indexpos] * dpois(y[indexpos], lambda = lambda[indexpos])
  RQR = qnorm(pzip + u * dzip)
  result = (RQR - mean(RQR[!(is.infinite(RQR) | is.nan(RQR))])) / sd(RQR[!(is.infinite(RQR) | is.nan(RQR))])
  return(result)
}

# from the following blog:
# https://beckmw.wordpress.com/2013/02/05/collinearity-and-stepwise-vif-selection/
# stepwise VIF function used below
vif_func<-function(in_frame,thresh=10,trace=T,...){
  require(fmsb)
  if(class(in_frame) != 'data.frame') in_frame<-data.frame(in_frame)
  #get initial vif value for all comparisons of variables
  vif_init<-NULL
  for(val in names(in_frame)){
    form_in<-formula(paste(val,' ~ .'))
    vif_init<-rbind(vif_init,c(val,VIF(lm(form_in,data=in_frame,...))))
  }
  vif_max<-max(as.numeric(vif_init[,2]))
  if(vif_max < thresh){
    if(trace==T){ #print output of each iteration
      prmatrix(vif_init,collab=c('var','vif'),rowlab=rep('',nrow(vif_init)),
               quote=F)
      cat('\n')
      cat(paste('All variables have VIF < ', thresh,', max VIF ',
                 round(vif_max,2), sep=''),'\n\n')
    }
    return(names(in_frame))
  }
  else{
    in_dat<-in_frame
    #backwards selection of explanatory variables
    #stops when all VIF values are below ’thresh’
    while(vif_max >= thresh){
      vif_vals<-NULL
      for(val in names(in_dat)){
        form_in<-formula(paste(val,' ~ .'))
        vif_add<-VIF(lm(form_in,data=in_dat,...))
        vif_vals<-rbind(vif_vals,c(val,vif_add))
      }
      max_row<-which(vif_vals[,2] == max(as.numeric(vif_vals[,2])))[1]
      vif_max<-as.numeric(vif_vals[max_row,2])
      if(vif_max<thresh) break
      if(trace==T){ #print output of each iteration
        prmatrix(vif_vals,collab=c('var','vif'),rowlab=rep('',nrow(vif_vals)),
                 quote=F)
        cat('\n')
        cat('removed: ',vif_vals[max_row,1],vif_max,'\n\n')
        flush.console()
      }
      in_dat<-in_dat[,!names(in_dat) %in% vif_vals[max_row,1]]
    }
    return(names(in_dat))
  }
}