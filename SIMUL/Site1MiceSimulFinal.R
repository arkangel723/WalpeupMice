#SIMULACION 
simula8<-function(No,Mobj,years,Nt,V1,V2){
  
  # librerias
  require(nlstools)
  require(MASS)
  
  # Extract parameters from the modelo object
  param <- c(Mobj$m$getAllPars())
  set.seed(6)
  pr <- mvrnorm(5000, param, vcov(Mobj))
  
  
  # Create object for output matrices and make the parametric sampling
  p1   <- No
  N2   <- matrix(NA, nrow = 22, ncol = 5000)
  N3   <- matrix(NA, nrow = 22, ncol = 5000)
  qqN3 <- matrix(NA, nrow = 22, ncol = 3)
  qqN2 <- matrix(NA, nrow = 22, ncol = 3)
  
  # Nested loop for CI and predictions based in the 5000 sampled parameters OJO la t es i
  for(t in 2:length(V1)){
    N2[1,]<-p1
    N1b<-Nt[t-1]
    for(j in 1:length(pr[,1])){
      
      # model prediction with all generated parraneters N distributed 
      N2[t,j] <- N2[t-1,j]*exp(2.5*(1-(N2[t-1,j]/pr[j,1])^pr[j,2]))+pr[j,3]*V1[t-1]+pr[j,4]*V2[t-1]
                         
      # model prediction with all generated parraneters N distributed
      N3[t,j] <- N1b*exp(2.5*(1-(N1b/pr[j,1])^pr[j,2]))+pr[j,3]*V1[t-1]+pr[j,4]*V2[t-1]
      
      # CI estimation trough gettin the .25, .5 and .975
      # quantiles of the normal distribution for each prediction (year)
      qqN3[t,] <- quantile((N3[t,]),c(0.025,0.5,0.975), na.rm = T)
      qqN2[t,] <- quantile((N2[t,]),c(0.025,0.5,0.975), na.rm = T)
p1<-qqN2[t,2]
    }
  }
  na.pad <- function(x,len){
    x[1:len]
  }
  makePaddedDataFrame <- function(l,...){
    maxlen <- max(sapply(l,length))
    data.frame(lapply(l,na.pad,len=maxlen),...)
  }
  
  # create data frame for de predictions and CI bands
  out.est <- makePaddedDataFrame(list(years = years, 
                                      Nt = Nt,
                                      Sim1 = qqN3[,2],
                                      Sim2 = qqN2[,2], 
                                      Sim1.CI.L = qqN3[,1],
                                      Sim1.CI.U = qqN3[,3],
                                      Sim2.CI.L = qqN2[,1],
                                      Sim2.CI.U = qqN2[,3]))
print(out.est)
  
  # Plotting the Nt data
  par(mar=c(4.5,4,1,2.5),lwd=2,bty="l",cex=0.6)
  plot(out.est[,2]~out.est[,1], xlab="Time", ylab="Abundance",
       cex.lab=0.9, cex.axis=0.9, tck=-0.03,
       type="p",
       pch=16,
       ylim=c(min(out.est[,2],na.rm=TRUE)-1,1.1*max(out.est[,2],na.rm=TRUE)),
       xlim=c(min(out.est[,1],na.rm=TRUE),max(out.est[,1],na.rm=TRUE)))
  
 
  # Ploting prediction and CI curves
  points(out.est[,2]~out.est[,1],type="p", pch=16)
  # N2 prediction
  #points(out.est[,4]~out.est[,1],col="red",type="l")
  # N3 prediction
  points(out.est[,3]~out.est[,1],col="blue",type="l")
  
  # Poligon for the 95% CI region
  
  # CI Simulation I
  polygon(c(out.est[,1],rev(out.est[,1])),
           c(out.est$Sim1.CI.U, rev(out.est$Sim1.CI.L)),
           col=rgb(.224, .224, .224,0.3), border=NA)
  
  # CI Simulation II
  #polygon(c(out.est[,1],rev(out.est[,1])),
          #c(out.est$Sim2.CI.U, rev(out.est$Sim2.CI.L)),
          #col=rgb(.224, .224, .224,0.3), border=NA)
  return(out.est)
}

simm1<-simula8(8.60,Walpeup8model8,seq(1983,2004),WalpeupNew$N,WalpeupNew$LLAG,WalpeupNew$ELAG)
param <- c(Walpeup8model8$m$getAllPars())
