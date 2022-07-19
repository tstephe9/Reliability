library(dplyr)
library(in2extRemes)

#Time Code
ptm <- proc.time()

#Load annual maximum flood series
peakQ <- read.table("E:/UGA/Files/RAS/LittleSugar/Flows/PeakQ.txt", header = TRUE, 
                    sep = "\t",colClasses = "character")
DA <- 42.6 #Drainage area in sq miles
DA <- 42.6*2.59 #Drainage area sq km

##Set start year
start <- 1946
peakQ$ams <- as.numeric(peakQ$Qcfs)
peakQ$ams <- (peakQ$ams/35.3147) #cms/sq km
peakQ$I <- as.numeric(peakQ$I)


#Fit a stationary set of model parameters
#sta <- fevd(ams, data = peakQ, type = "Gumbel", method = "MLE", time.units = "years", period.basis = "year")
# sta_ci <- ci(sta, alpha = 0.05, type = c("return.level"),
#    return.period = RI, which.par, R = 1000, method = "boot")
# 
# write.table(sta_ci, file = "C:/Users/tas48127/Desktop/Files/RAS/LittleSugar/STA_CI.txt",
#             row.names=FALSE, col.names = FALSE, na="", sep="\t")


#For Little Sugar Creek, everything prior to 1945 is stationary
npeakQ <- peakQ[23:96,]
npeakQ$I <- seq(0,length(npeakQ$I)-1,1)

#Fit a nonstationary set of model parameters
nsta <- fevd(ams, data = npeakQ,location.fun = ~I,scale.fun = ~I, type = "Gumbel", time.units = "years", period.basis = "year")
z <- summary(nsta, silent = TRUE)


tobs <- npeakQ$I #years of observation
#set the planning horizon
ph <- 30
yearone <- 2014 #2014 water year
yearind <- 67 #index for year 2014
nsend <- 2024 #The end of nonstationarity, stationary assumption applies afterward.
yearend <- yearone + ph - 1
nsind <- yearind + nsend - yearone #index for end of nonstationarity
nslen <- nsend - yearone #length of nonstationary period
ustalen <- ph - nslen #length of updated stationary period
realind <- yearind + ph - 1

#Create matrix of years to quantify reliability
tfut <- c(seq(yearind, nsind), rep(nsind, ustalen - 1))
check <- ph - (nslen + ustalen) + ph - length(tfut) #make sure all indexing correct
#record mut for each year
nsfit <- z$par[1] + (z$par[2]*tfut)


# #Plot the fitted model
# #Actual data
# plot(npeakQ$I + start, npeakQ$ams, ylab = "Discharge (cfs)", xlab = "year", xlim = c(min(tobs)+start, max(tfut)+start))
# #Fitted model for observations
# lines(tobs + start, nsfit, col = "red")
# #Model projected to t
# lines(tfut + start, nsfuture, col = "blue")


#x = mean resistance
#Rl = reliability to solve for
#nsta = nonstationary model from in2Extremes
#covr = coeff. variation in resistance
#ph = planning horizon
 #x <- 21500
 #covr <- 0.1
 #Rl <- 0.74

fr <- function(x,covr, Rl, nsta, t){
  #Get resistance distribution parameters
  sdy <- sqrt(log(1+(covr^2)))
  muy <- log((x/((1+(covr^2))^0.5)))
  z <- summary(nsta, silent = TRUE)
  #Check to make sure not solving for present (ph = 0)
  
  #years reliability evaluated over
  t <- t
  
  #variable to store annual reliability
  Rel_store <- rep(NA, length(t))
  
  for (i in 1:length(t)){
    #Calculate parameters for loading distribution at year ti
    mutt <- z$par[1] + (z$par[2]*t[i])
    sdtt <- z$par[3] + (z$par[4]*t[i])
    diff <- 1
    ii <- 0
    n <- 5 #set the incremental discharge
    while (diff > 0.0001){
      ii <- ii + 1
      if (ii == 1){
        r <- seq(n,90000,n)
        fr <- (1/(r*sdy*sqrt(2*pi)))*exp(-0.5*(((log(r)-muy)/sdy)^2))
        Fl <- exp(-exp(-((r-mutt)/sdtt)))
        Flfr <- fr*Fl
        a <- Flfr[1:(length(Flfr)-1)]
        b <- Flfr[2:length(Flfr)]
        h <- ((b+a)/2)*n
        approx <- sum(h)
      } else if (ii == 2) {
        r <- seq(n,95000,n)
        Flfr <- (1/(sqrt(2*pi)*sdy*r)) * exp((-0.5*(((log(r)-muy)/sdy)^2))) * exp(-exp(-((r-mutt)/sdtt)))
        a <- Flfr[1:(length(Flfr)-1)]
        b <- Flfr[2:length(Flfr)]
        h <- ((b+a)/2)*n
        diff <- sum(h) - approx
        approx <- sum(h)
      } else {
        r <- seq(n,9*(10^(ii+3)),n)
        Flfr <- (1/(sqrt(2*pi)*sdy*r)) * exp((-0.5*(((log(r)-muy)/sdy)^2))) * exp(-exp(-((r-mutt)/sdtt)))
        a <- Flfr[1:(length(Flfr)-1)]
        b <- Flfr[2:length(Flfr)]
        h <- ((b+a)/2)*n
        diff <- sum(h) - approx
        approx <- sum(h)
      }
    }
    Rel_store[i] <- approx
  }
  Rli <- prod(Rel_store[1:i])
  #Rli <- prod(1-(1-exp(-exp(-(x[1] - (z$par[1] + z$par[2]*t))/(z$par[3] + z$par[4]*t))))) #calculate reliability according Salas and Obeysekera, 2014
  ((Rli*10) - (Rl*10))^2 #The difference between the calculated (Rli) and desired (Rl) reliability
  #The above equation is optimized to a minimum by changing input x (i.e. Q) - the Q to providie
  #desired reliability over ph years into the future based on the nonstationary model
}
##############################################################################################
#bootstrap method
nyears <- npeakQ$I #year of block maxima
nQs <- npeakQ$ams #value of block maxima


mut <- z$par[1] + z$par[2]*tfut[1] #fitted mean flow at time t
sdt <- z$par[3] + z$par[4]*tfut[1] #fitted std at time t
#shape <- z$par[5]

###  Reliability to Solve For  ######
Rset <- 0.9
###################################

#Set up
nb <- 1000 #number of bootstraps
boot <- as.data.frame(matrix(NA, nrow = nb, ncol = length(Rset))) #dataframe for the solved Q's
colnames(boot) <- Rset #Set the column names to write the results to
conv <- boot #dataframe to write the convergence codes to
meansample <- as.data.frame(matrix(data = NA, nrow = nb, ncol = 2)) #record location parameters
colnames(meansample) <- c("mu0","trend")
sdsample <- as.data.frame(matrix(data = NA, nrow = nb, ncol = 2)) #record scale parameters
colnames(sdsample) <- c("sd0", "trend")

#Fitted nonstationary model
allmu <- z$par[1] + z$par[2]*tobs #fitted mean flow at time t
allsd <- z$par[3] + z$par[4]*tobs #fitted std at time t
normresid <- (nQs - allmu)/allsd  #Normalized residuals for bootstrapping
#normresid <- (1/shape)*log(1+(shape*((nQs - mut)/sdt))) #GEV normalized residuals based on fitted model at time t
#normresid <- (nQs - mut)/sdt #Gumbel normalized residuals based on fitted model at time t

#Plot the actual data
plot(npeakQ$I, npeakQ$ams, ylab = "Discharge (cms)", xlab = "year", xlim = c(min(tobs), realind))

#Conduct the bootstrap
for (zz in 1:nb){
  sampled <- sample(normresid, size = length(normresid), replace = TRUE) #sample with replacement
  transformed <- ((sampled*allsd)+allmu) #back-transform normalized residuals
  #transformed <- (((exp(sampled*shape)-1)/shape)*sdt)+mut #"back-transformed" residuals
  bpeakQ <- data.frame(transformed, nyears) #create dataframe of sampled (i.e. back-transformed) flood series
  nsta <- fevd(transformed, data = bpeakQ,location.fun = ~nyears,scale.fun = ~nyears, type = "Gumbel", time.units = "years", period.basis = "year")
  sampledloading <- summary(nsta, silent = TRUE) #Summary of nonstationary model fitted to bootstrapped sample
  meansample[zz,1] <- nsta$results$par[1] #mu0
  meansample[zz,2] <- nsta$results$par[2] #mu slope
  sdsample[zz,1] <- nsta$results$par[3] #sd0
  sdsample[zz,2] <- nsta$results$par[4]#sdslope
  obsts <- sampledloading$par[1] + (sampledloading$par[2]*tobs)
  relts <- sampledloading$par[1] + (sampledloading$par[2]*tfut)
  lines(tobs, obsts, col = "lightgray") #plot the observed values
  lines(seq(yearind,realind,1), relts, col = "lightgray") #plot the time analyzed for reliability
  #########################################
  #Solve for design Q of Bootstrapped sample
  setQ <- 500 #the initial guess
  covr <- 0.1 #COV of resistance
  #fr <- function(x,covr, Rl, nsta, ph)
  #test <- optim(setQ, fr, Rl = Rset, nsta = nsta, covr = covr, ph = length(tfut)) #optimize the function
  test <- optimize(fr,c(1, 2000), Rl = Rset, nsta = nsta, covr = covr, t = tfut, tol = 0.001, maximum = FALSE) #optimize the function
  boot[zz,1] <- test$minimum #Store solution (design Q) for optimize
  conv[zz,1] <- test$objective #Store convergence for optimize
  # boot[zz,1] <- test$par[1] #Store solution (design Q) for optim
  # conv[zz,1] <- test$convergence[1] #Store convergence code for optim
  ########################################

}
# p <- tobs[length(tobs)]
# sdp <- sdsample[,1] + (sdsample[,2]*p)
# mup <- meansample[,1] + (meansample[,2]*p)
# Q100 <- mup-(sdp*log(-log(1-0.01)))
# 
# #Fitted model for observations
# points(npeakQ$I, npeakQ$ams)
# lines(tobs, nsfit, col = "red")
# #Model projected to t
# lines(tfut, nsfuture, col = "blue")
proc.time() - ptm

cdf <- quantile(boot[,1], probs = c(0.1, 0.5, 0.9))
cdf
Rset
covr

out <- cbind(boot, meansample, sdsample)
colnames(out)<- c("R0.9", "meana", "meanb", "sda", "sdb")

write.table(out, file = "E:/UGA/Files/RAS/LittleSugar/DesignQ_R0.90_30yrph_NS2014to2024_cov0.1.txt",
            row.names=FALSE, col.names = TRUE, na="", sep="\t")



ci_results <- as.data.frame(matrix(NA, nrow = length(Rset), ncol = 3)) #Create dataframe to write out results
colnames(ci_results) <- c("pt", "Qestimate", "sd")
ci_results$pt <- Rset 
ci_results$sd <- sapply(boot[,], FUN = sd) #Calculate the sd of Q for each AEP of the bootstrap (n = nb above, rows = number of AEPs)
mut <- z$par[1] + z$par[2]*max(tfut) #fitted mean flow at time t
sdt <- z$par[3] + z$par[4]*max(tfut) #fitted std at time t
ci_results$Qestimate <- mut - (sdt*log(-log(Rset))) #estimate of mean annaual flow at future time t for each AEP







