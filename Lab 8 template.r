#####################################################################
#Lab 8: Bayesian model of yelloweye rockfish
#Template
#Trevor A. Branch, tbranch@uw.edu
#1 April 2015
#####################################################################

#============getNLL=====================================
#The negative log-likelihood of the model given the 
#data. Assumed exponential model with base CV=0.02
#plus additional variance. 
#=======================================================
getNLL <- function(catches, CPUE1, CPUE2, r, K, q1, q2) {
   start.yr <- 1970
   end.yr <- 2011
   years <- start.yr:end.yr
   nyr <- length(years)
   B <- vector(length=nyr)
   B[1] <- K*0.8  #stock depleted to 80% of K in 1970
   NLLprior <- 0
   for (i in 1:(nyr-1)) {
      #logistic model, the max() part ensures it does not go extinct
      B[i+1] = max(0.01, B[i] + r*B[i]*(1-B[i]/K) - catches[i,2]) 
   }
   #complicated bit where the model figures out if there is a 
   #CPUE data point in the particular year of interest or not
   for (i in 1:nyr) {
      current.yr <- years[i]  #find the current year of the model
      #if there is a CPUE1 in that year then add NLL
      if (current.yr %in% CPUE1[,1]) {
         obsCPUE <- CPUE1[CPUE1[,1]==current.yr,2]
         obsCV <- CPUE1[CPUE1[,1]==current.yr,3]
         predCPUE <- B[i]*q1
         #print(paste(obsCPUE, predCPUE, obsCV))
         NLLprior = NLLprior + (log(obsCPUE/predCPUE))^2/(2*obsCV^2)
      }
      #if there is a CPUE2 in that year then add NLL
      if (current.yr %in% CPUE2[,1]) {
         obsCPUE <- CPUE2[CPUE2[,1]==current.yr,2]
         obsCV <- CPUE2[CPUE2[,1]==current.yr,3]
         predCPUE <- B[i]*q2
         NLLprior = NLLprior + (log(obsCPUE/predCPUE))^2/(2*obsCV^2)
      }
   }
   #add -ln(prior) for r, K, q1, q2
   #r~U[0,0.2], K~U[9766,100000], q1~U[0.00001,0.01], q2~U[1e-7, 1e-3]
   NLLprior <- NLLprior -log(1/0.2) +  #prior on r 
      -log(1/(100000-9766)) +   #prior on K
      -log(1/(0.01-0.00001)) +  #prior on q1
      -log(1/(1e-3 - 1e-7))     #prior on q2
   
   return(c(NLLprior,B[nyr]))  #return the NLL and final B
}

#=================================================
#Run the MCMC chain
#=================================================
runMCMC <- function(Cfile, CPUE1file, CPUE2file, 
                    ndraws,
                    rinit, Kinit, q1init, q2init) {
   #read in the data (only want to do this once!)
   catches <- read.csv(file=Cfile, header=T, colClasses="numeric")
   CPUE1 <-   read.csv(file=CPUE1file, header=T, colClasses="numeric")
   CPUE2 <-   read.csv(file=CPUE2file, header=T, colClasses="numeric")
   
   #create matrix to store the sampled draws
   posterior <- matrix(nrow=ndraws, ncol=6, 
                       dimnames=list(paste("Draw",1:ndraws),
                                     c("r","K","q1","q2", "NLLprior", "B2011")))
   
   #calculate the NLLprior for the initial sample
   NLLinit <- getNLL(catches=catches, CPUE1=CPUE1, CPUE2=CPUE2,
                     r=rinit, K=Kinit, q1=q1init, q2=q2init)
   posterior[1,] <- c(rinit, Kinit, q1init, q2init, 
                      NLLinit[1], NLLinit[2])
   naccepted <- 1
   noutbounds <- 1
   
   #r~U[0,0.2], K~U[9766,100000], q1~U[0.00001,0.01], q2~U[1e-7, 1e-3]
   for (i in 2:ndraws) {
      if (i %% 1000 == 0) {  #%% is function returning remainder of i/1000
         print(paste("draw",i,"of",ndraws))
      }
      #the jump function
      rstar <- posterior[i-1,1] + runif(n=1, min=-0.0005, max=0.0005)
      Kstar <- posterior[i-1,2] + runif(n=1, min=-30,max=30)
      q1star <- posterior[i-1,3] + runif(n=1, min=-0.0001,max=+0.0001)
      q2star <- posterior[i-1,4] + runif(n=1, min=-5e-6, max=5e-6)
      
      #check to see if draw outside parameter bounds, if so, stay with previous draw
      if (rstar<=0 || rstar>=0.2 || Kstar<=9766 || Kstar>=100000 || 
             q1star<=0.00001 || q1star>=0.01 || q2star<=1e-7 || q2star>=1e-3) {
         noutbounds <- noutbounds+1
         posterior[i,] <- posterior[i-1,]  #outside bounds, reject draw
      }
      else {
         NLLstar <- getNLL(catches=catches, CPUE1=CPUE1, CPUE2=CPUE2,
                           r=rstar, K=Kstar, q1=q1star, q2=q2star)
         
         #THE KEY MCMC step: whether to accept or reject the new values
         ratio <- exp(posterior[i-1,5] - NLLstar[1])
         
         if (ratio >= 1) {  #better place, accept draw automatically
            posterior[i,] <- c(rstar, Kstar, q1star, q2star, NLLstar[1], NLLstar[2])
            naccepted <- naccepted + 1
         } else {
            if (runif(n=1, min=0, max=1) < ratio) { #accept with probability
               posterior[i,] <- c(rstar, Kstar, q1star, q2star, NLLstar[1], NLLstar[2])
               naccepted <- naccepted + 1
            }
            else {
               posterior[i,] <- posterior[i-1,]  #reject and keep last draw
            }
         }
      }
   }
   #plot the posterior
   par(mfrow=c(3,2), oma=c(1,1,1,1), mar=c(2,4,1,1))
   plot(posterior[,1], type="l", ylab="r")
   plot(posterior[,2], type="l", ylab="K")
   plot(posterior[,3], type="l", ylab="q1")
   plot(posterior[,4], type="l", ylab="q2")
   plot(posterior[,5], type="l", ylab="NLL+NLprior")
   plot(posterior[,6]/posterior[,2], type="l", ylab="B2011 / K")  #finalB/K = how depleted stock is
   
   print("accepted divided by draws")
   print(naccepted/ndraws)
   print("draws outside bounds")
   print(noutbounds/ndraws)
   
   #return the results
   return(posterior)
}
x <- runMCMC(Cfile="8 catches.csv", CPUE1file="8 early CPUE.csv", 
             CPUE2file="8 late CPUE.csv", ndraws=1000,
             rinit=0.05, Kinit=12000, q1init=0.002, q2init=1e-4)
acf(x[,1])[1]
acf(x[,2])[1]
acf(x[,3])[1]
acf(x[,4])[1]
Y <- x[seq(20100,to=100000,by=500),]  #burn in and thinning
acf(Y[,1])[1]
acf(Y[,2])[1]
acf(Y[,3])[1]
acf(Y[,4])[1]

x <- runMCMC(Cfile="8 catches.csv", CPUE1file="8 early CPUE.csv", 
             CPUE2file="8 late CPUE.csv", ndraws=100000,
             rinit=0.05, Kinit=12000, q1init=0.002, q2init=1e-4)
