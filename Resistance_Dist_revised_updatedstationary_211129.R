library(dplyr)
library(in2extRemes)

#Time Code
ptm <- proc.time()

#Set the River station to analyze
RiverSTA <- 34056
#Set the elevation to evaluate (e.g. a base flood elevation)
BFE <- 567.73/3.28
newdata <- data.frame("lny" = log(BFE)) #ln of BFE for predicting Q
#Load monte-carlo results
results <- read.table("E:/UGA/Files/RAS/LittleSugar/MonteCarlo/LS_MC_Rel_1000NS0.txt", header = TRUE, 
                    sep = "\t",colClasses = "character")

# #Remove structure locations
results <- results[results$nodeType == "", ]
results<- sapply(results[,], as.numeric) #Convert values to numeric
results <- as.data.frame(results) #Convert back to df

#create a dataframe of the results for the specified station
a <- results[results$STA == RiverSTA,]
a$Q <- a$Q/35.3147
a$WSE <- a$WSE/3.28
a$lnQ <- log(a$Q) #ln Q
a$lny <- log(a$WSE) #ln WSE

#plot all observation of stage discharge
plot(a$WSE, a$Q, log = "y", xlab = "Water Surface Elevation (ft)", ylab = "Discharge (cfs)",
     col = "grey") 
Q <- rep(NA, max(a$simulation))
plot(1, type="n", xlim=c(560/3.28,580/3.28), ylim=c(min(a$Q), (max(a$Q))+50),
     ylab = "Discharge (cms)", xlab = "Water Surface Elevation (m)")

for (i in 1:max(a$simulation)){
  sim <- a[a$simulation == i,]
  fit<- lm(lnQ~poly(lny,2), sim)
  z <- summary(fit)
  R2 <- z$adj.r.squared
  
  if (!is.nan(R2) & R2 < 0.95) {
    fit<- lm(lnQ~poly(lny,3), sim)
    z <- summary(fit)
    R2 <- z$adj.r.squared
    if (!is.nan(R2) & R2 < 0.95) {
      next
    }
  }
  pp <- exp(fitted(fit))
  pp <- pp#/35.3147
  lines(sim$WSE,pp, col = "gray80")
  Q[i] <- predict(fit, newdata)
}



quants <- seq(from = 0.01, to = 0.99, by = 0.01)
ecdfQ <- quantile(Q, probs = quants)
plot(exp(ecdfQ), quants, type = "l", lwd = 2, xlab = "Discharge (cms)", ylab = "Quantile", yaxt = 'n')
grid(ny = 9, col = "lightgray", lty = "solid")
axis(side = 2, at = seq(0,1,0.2), labels = seq(0,1,0.2))
lines(exp(ecdfQ), quants, type = "l", lwd = 2)

qq <- qqnorm(Q, pch = 1, frame = FALSE)
plot(qq$qnorm, qq$data, pch = 1, xlab = "Standard Normal Quantiles", ylab = "Sample Quantiles")
qqline(Q, col = "steelblue", lwd = 2)

mur <- mean(exp(Q))
sdr <- sd(exp(Q))
COVr <- sdr/mur

#######################################################################################################
########################################################################################################

#fit a nonstationary model to peak flow data
#Load annual maximum flood series
peakQ <- read.table("E:/UGA/Files/RAS/LittleSugar/Flows/PeakQ.txt", header = TRUE, 
                    sep = "\t",colClasses = "character")
DA <- 42.6 #Drainage area in sq miles
peakQ$ams <- as.numeric(peakQ$Qcfs)
#Adjust to DA of STA
peakQ$ams <- (peakQ$ams/DA)*45.9 #cfs/sq mi
peakQ$ams <- peakQ$ams/35.3147
peakQ$I <- as.numeric(peakQ$I)
#For Little Sugar Creek, everything prior to 1945 is stationary
npeakQ <- peakQ[23:96,]
npeakQ$I <- seq(0,length(npeakQ$I)-1,1)
#Fit a nonstationary set of model parameters
nsta <- fevd(ams, data = npeakQ,location.fun = ~I,scale.fun = ~I, type = "Gumbel", time.units = "years", period.basis = "year")
z <- summary(nsta, silent = TRUE)

#################################################################################################################################
#################################################################################################################################
# #set the number of years in the future
#tfut <- seq(from = max(npeakQ$I)+1,to = max(npeakQ$I)+100, by = 1) #Little Sugar Record starts at 74
tfut <- seq(from = 67,to = 67+50 - 1, by = 1) #Little Sugar Record starts at 74
tobs <- npeakQ$I
#Set year 1 for reliability calculation
yearone <- 2014 #1970 water year
yearind <- 67 #index for year 2014
yearend <- 2024 #The end of nonstationarity, stationary assumption applies afterward.
nsind <- yearind + (yearend - yearone) + 1 #index for end of nonstationarity


#Perform a bootsrap to determine uncertainity
#bootstrap method
nyears <- npeakQ$I #year of block maxima
nQs <- npeakQ$ams #value of block maxima
nb <- 1000 #number of bootstraps
boot <- as.data.frame(matrix(NA, nrow = nb, ncol = length(tfut))) #dataframe for model fitting
boot_stats <- as.data.frame(matrix(NA, nrow = nb, ncol = 4)) #dataframe for model fitting
nstaboot <- boot #dataframe for nsta results excluding channel capacity uncertainty
colnames(boot) <- c(tfut)
colnames(boot_stats) <- c("mu", "muslope", "sd", "sdslope") 
#dataframe  rows = number of bootstraps
#dataframe columns = number of year

sdy <- sqrt(log(1+(COVr^2))) #lognormal resistance sd
muy <- log((mur/(sqrt((1+(COVr^2)))))) #lognormal resistance mean

#Trapezoidal rule integration fucntion
#r is resistance
#sdy is lognormal sd
#muy lognormal mean
#mut is Gumbel mean of annual flood peaks
#sdt is Gumbel sd of annual flood peaks
#n is the delta r (dr) for integration
integrand <- function(r, sdy, muy, mut, sdt, n) {
  Flfr <- (1/(sqrt(2*pi)*sdy*r)) * exp((-0.5*(((log(r)-muy)/sdy)^2))) * exp(-exp(-((r-mut)/sdt)))
  a <- Flfr[1:(length(Flfr)-1)]
  b <- Flfr[2:length(Flfr)]
  h <- ((b+a)/2)*n
  approx <- sum(h)
  return(approx)
}
# parameterize dr
n <- 1
#Create resistances (r) to integrate over, ensure upper bound is large enough
Qs <- seq(from = 5, to = 50000, by = n)

#Fitted nonstationary model
allmu <- z$par[1] + z$par[2]*tobs #fitted mean flow at time t
nsmu <- z$par[1] + z$par[2]*seq(tobs[length(tobs)], nsind)
nsfuture <- z$par[1] + z$par[2]*rep(nsind, length(tfut)- (yearend - yearone)+1)
allsd <- z$par[3] + z$par[4]*tobs #fitted std at time t
normresid <- (npeakQ$ams - allmu)/allsd  #Normalized residuals for bootstrapping

# mut <- z$par[1] + z$par[2]*max(tfut[1]) #fitted mean flow at time t
# sdt <- z$par[3] + z$par[4]*max(tfut[1]) #fitted std at time t

#Plot the actual data
plot((npeakQ$I+1946), npeakQ$ams, ylab = "Discharge (cfs)", xlab = "year", xlim = c(min(tobs+1946), max(tfut+1946)))
lines(tobs+1946, allmu, col = "red")
#Model projected to t
lines(seq(tobs[length(tobs)]+1946, yearend), nsmu, col = "blue")
lines(seq(yearend, length(tfut) + yearone), nsfuture, col = "blue")



#Conduct the bootstrap
for (i in 1:nb){
  sampled <- sample(normresid, size = length(normresid), replace = TRUE) #sample with replacement
  transformed <- ((sampled*allsd)+allmu) #back-transform normalized residuals
  bpeakQ <- data.frame(transformed, nyears) #create dataframe of sampled (i.e. back-transformed) flood series
  nsta <- fevd(transformed, data = bpeakQ,location.fun = ~nyears, scale.fun =  ~nyears, type = "Gumbel", time.units = "years", period.basis = "year")
  zz <- summary(nsta, silent = TRUE) #store parameters of fitted bootstrap Gumbel model
  boot_stats[i,1] <- zz$par[1]
  boot_stats[i,2] <- zz$par[2]
  boot_stats[i,3] <- zz$par[3]
  boot_stats[i,4] <- zz$par[4]
  Rk <- rep(NA, length(tfut)) #vector to write each year's non-exceedance probability to
  Rks <- rep(NA, length(tfut)) #vector to write each year's non-exceedance probability to not considering uncertainty
  for (k in 1:length(tfut)){ #loop through each year
    
    if(yearone - 1 + k <= yearend){
      #use nonstationary
      mutt <- zz$par[1] + zz$par[2]*(yearind - 1 + k) #fitted mean flow at last observed
      sdtt <- zz$par[3] + zz$par[4]*(yearind - 1 + k) #fitted std at time last observed
      Rel <- integrand(Qs, sdy, muy, mutt, sdtt, n) #Integrate reliability function
      #Rel <- 1
      Rk[k] <- Rel #store that year's reliability
      Rks[k] <-1-(1-(exp(-exp(-((mur-mutt)/sdtt)))))#reliability no channel capacity uncertainty
    }else{
      #use updated stationary
      mutt <- zz$par[1] + zz$par[2]*(nsind) #fitted mean flow at last observed
      sdtt <- zz$par[3] + zz$par[4]*(nsind) #fitted std at time last observed
      Rel <- integrand(Qs, sdy, muy, mutt, sdtt, n) #Integrate reliability function
      #Rel <- 1
      Rk[k] <- Rel #store that year's reliability
      Rks[k] <- 1-(1-(exp(-exp(-((mur-mutt)/sdtt)))))#reliability no channel capacity uncertainty
    }
    if (k == 1){
      boot[i,k] <- Rel
      nstaboot[i,k] <- Rks[k]
    } else {
      boot[i,k] <- Rel*prod(Rk[1:k-1]) #based on Reliability eqn of Salas and Obeysekera (2014)
      nstaboot[i,k] <- Rks[k]*prod(Rks[1:k-1])
    }
    
  }
  allts <- c(tobs,seq(tobs[length(tobs)]+1, nsind),rep(nsind+1, length(tfut)- (nsind - yearind))) #variable of all times
  plotsample <- zz$par[1] + (zz$par[2]*allts) #vector to plot time varying mean (location parameter) of bootstrapped sample
  lines(seq(1946, 1946+length(allts)-1), plotsample, col = "lightgray") #plot time varying mean (location parameter) of bootstrapped sample
}

lines(tobs+1946, allmu, col = "red")
#Model projected to t
lines(seq(tobs[length(tobs)]+1946, yearend), nsmu, col = "blue")
lines(seq(yearend, length(tfut) + yearone), nsfuture, col = "blue")
proc.time() - ptm
# #Remove faulty integration
# out <- boot[boot[,40] != max(boot[,40]),]
out <- boot
nstaout <- nstaboot

#Plot the results
plot(seq(1,length(tfut)),out[1,], typ = "l", col = "gray",xlab = "Planning Horizon", ylab = "Reliability", ylim = c(0,1))
for (i in 2:length(boot[,1])){
  lines(seq(1,length(tfut)), out[i,], col = "gray")
}



#Get the quantiles of reliability
ci_results <- apply(out,2, quantile, probs=c(0.025,0.5,.975), na.rm=TRUE)
ci90 <- apply(out,2, quantile, probs=c(0.1), na.rm=TRUE)
#
nstaci_results <- apply(nstaout,2, quantile, probs=c(0.025,0.5,.975), na.rm=TRUE)

#Fit a stationary distribution and calculate reliability
sta <- fevd(ams, data = npeakQ, type = "Gumbel", method = "MLE", time.units = "years", period.basis = "year") #fit the model
sz <- summary(sta) #store the results
stap <-1-(exp(-exp(-((mur-sz$par[1])/sz$par[2])))) #calculate non-exceedance probability (1-pt)
stapt <- seq(from = 1, to = length(tfut), by = 1) #Create a vector to calculate stationary reliability for each year in a planning horizon
staR <- (1-stap)^stapt #Calculate stationary reliability

#Parameters to create a polygon for 95% CI of reliability
min_a <- pmin(ci_results[1,], ci_results[3,])
max_a <- pmax(ci_results[1,], ci_results[3,])

#Plot the results
plot(seq(1,length(tfut)),ci_results[2,], type = "l", col = "black", xlab = "Planning Horizon", ylab = "Reliability", ylim = c(0,1)) #avg reliability w/ uncertainty
polygon(c(seq(1,length(tfut)), rev(seq(1,length(tfut)))), c(ci_results[1,],rev(ci_results[3,])), col=adjustcolor("gray89", alpha=0.5), border = "gray89") #95% CI of reliability
lines(seq(1,length(tfut)),ci_results[2,], col = "black", lwd = 2) #redraw avg reliability w/ uncertainty
lines(seq(1,length(tfut)),nstaci_results[2,], col = "black", lwd = 2, lty=2) #redraw avg nonstationary reliability w/ out uncertainty
lines(seq(1,length(tfut)),ci90, col = "blue", lwd = 2) #90% confidence in reliability (i.e. 10% quantile)
lines(seq(1,length(tfut)),staR, col= "red", lwd = 2) #Avg. Stationary reliability without uncertainty
legend(55, 1, legend=c("90% Conf.", "50% Conf.", "Avg. Stationary"),
       col=c("blue", "black", "red"), lty=c(1,1,1), cex=0.8, lwd = 2)

#Plot the pdf of reliability
d<- quantile(out[,30], probs = seq(0.01, 0.99, 0.01))
plot(d,seq(0.01, 0.99, 0.01), col = "black", main = "30-year Planning Horizon", 
     xlab = "Reliability", ylab = "CDF",
     type = "l", lwd = 2, xlim = c(0,0.2))

#Save any results
write.table(boot, file = "E:/UGA/Files/RAS/LittleSugar/Rel_Dist_STA34056_BFE_revised220709.txt",
            row.names=FALSE, col.names = FALSE, na="", sep="\t")
write.table(nstaboot, file = "E:/UGA/Files/RAS/LittleSugar/Rel_Dist_STA34056_BFE_revised_noUncertainty_220709.txt",
            row.names=FALSE, col.names = FALSE, na="", sep="\t")
write.table(boot_stats, file = "E:/UGA/Files/RAS/LittleSugar/PeakQ_Moments_Distribution_220709.txt",
            row.names=FALSE, col.names = TRUE, na="", sep="\t")

