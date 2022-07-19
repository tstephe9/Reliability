#Gumbel distribution of Annual Flood peaks

#Time Code
ptm <- proc.time()

#Define number of years
yr <- seq(1,500,1)

#Define p0
p0 <- 0.01

#Determine trend magnitude
#decadal magnification factor
m <- c(1, 1.05, 1.2, 1.5, 2)

#Assign COV of annual flood peak distribution
lcov <- c(0.5, 1, 1.5)
#Assign COV of channel capacity
rcov <- c(0.01, 0.1, 0.4)

#Calculate mean and sd of flood peaks
mul <- 50 #mean annual flood

#Calculate stationary hydrologic reliability
Rs <- (1-p0)^yr

sims <- length(lcov)*length(rcov)*length(m)*length(yr)
results <- as.data.frame(matrix(data = NA, nrow = sims, ncol = 10))
colnames(results) <- c("Rel","yr", "m", "rcov", "lcov", "mut", "sdt", "mur", "sdr", "p0")
Rel_store <- rep(NA,length(yr))
STARel <- as.data.frame(matrix(ncol = 3, nrow = max(yr), data = NA))
colnames(STARel)<-c("yr", "Rel", "p0")
STARel$yr <- seq(1,500,1)
STARel$p0 <- p0
STARel$Rel <- (1-p0)^yr
# m <- c(1, 1.25)
# rcov <- 0.01
# lcov <- 1

v <- 1
#Loop through the COV of loading
for (i in 1:length(lcov)){
  sdl <- mul*lcov[i] #sd of mean annual floods
  mur <- mul-(sdl*log(-log(1-p0))) #mean resistance
  
  #Loop throug the cov of resistance
  for (k in 1:length(rcov)){
  sdy <- sqrt(log(1+(rcov[k]^2)))
  sdr <- sqrt(((mur^2)*((exp(sdy^2))-1)))
  muy <- log((mur/((1+(rcov[k]^2))^0.5)))
    #Loop through the magnification factors
    for (j in 1:length(m)){
      beta <- ((m[j]*mul)-mul)/10 #trend slope based on decadal magnification factor
      #Loop through each year
      for (z in 1:length(yr)){
        n <- 1
        mutt <- (beta*yr[z])+mul
        sdtt <- sdl
        Rel_store[z] <- integrand(sdy, muy, mutt, sdtt, n) #Integrate reliability function
        if (z >1){
          Rel <- prod(Rel_store[1:z])
        } else {
          Rel <- Rel_store[z]
        }
        results[v,1] <- Rel
        results[v,2] <- yr[z]
        results[v,3] <- m[j]
        results[v,4] <- rcov[k]
        results[v,5] <- lcov[i]
        results[v,6] <- mutt
        results[v,7] <- sdtt
        results[v,8] <- mur
        results[v,9] <- sdr
        results[v,10] <- p0
        v <- v + 1
      }
    }
  }
}

proc.time() - ptm


write.table(results, file = "C:/Users/tas48127/Desktop/ReliabilityAnalysis_newCOV_25yr.txt",
            row.names=FALSE, col.names = TRUE, na="", sep="\t")

results <- read.table("E:/UGA/Files/RAS/LittleSugar/ReliabilityAnalysis_newCOV_p01.txt",
                  header = TRUE, sep = "\t",colClasses = "character")
results <- as.data.frame(sapply(results[,], as.numeric))
urcov <- unique(results$rcov)
ulcov <- unique(results$lcov)
results$M <- as.factor(results$m)
mm <- unique(results$m)
ilab <- c("(a)", "(b)", "(c)", "(d)", "(e)", "(f)", "(g)", "(h)", "(i)")

pdf("E:/UGA/Files/RAS/LittleSugar/Plots/Gen_Plot_newCOV_p01_220715.pdf", width = 7, height = 7)
#png("C:/Users/tas48127/Desktop/Gen_Plot.png", width = 550, height = 550, res = 83)
par(mfrow=c(3,3),
            oma = c(3,3,0,0) + 0.1,
            mar = c(3,1.5,1,1) + 0.1)
ccc <- c("blue", "orange", "brown", "red")
ltt <- c(2,1,4,5,1)

kk <- 1
for (i in 1:length(ulcov)){
  for (k in 1:length(urcov)){
    plot1 <- results[results$rcov == urcov[k] & results$lcov == ulcov[i] & results$m == 1,]
    text <- bquote(.(ilab[kk])~~COV[l]==.(ulcov[i])~~COV[r]==.(urcov[k]))
    plot(plot1$yr, plot1$Rel, xlab = "", ylab = "",
         xlim = c(0,200), ylim = c(0,1),type = "l", lty = 2, lwd = 1.3, main = text,
         col = "black")
    lines(STARel$yr, STARel$Rel, col = "gray",lty = 1, lwd = 1.3)
    lines(plot1$yr, plot1$Rel, type = "l", lty = 2, lwd = 1.3, col = "black")
    for (z in 2:length(mm)){
      plot1 <- results[results$rcov == urcov[k] & results$lcov == ulcov[i] & results$m == mm[z],]
      lines(plot1$yr, plot1$Rel, lty = ltt[z], lwd = 1.3, col = ccc[z-1])
    }
    legend("topright", legend= c("1, NEP", mm),lwd = 1.3,col = c("gray","black","blue", "orange", "brown", "red"), title = "M", lty = c(1,ltt))
    kk <- kk + 1
  }
}
mtext(text = "Reliability", side = 2, line = 1, outer = TRUE)
mtext(text = "Planning Horizon (years)", side = 1, line = 0, outer =TRUE)
dev.off()










p <- list()
z <- 1
for (i in 1:length(ulcov)){
  for (k in 1:length(urcov)){
    plot1 <- results[results$rcov == urcov[k] & results$lcov == ulcov[i],]
    text <- bquote(COV[l]==.(ulcov[i])~~COV[r]==.(urcov[k]))
    p[[z]] <- ggplot() +
      geom_line(data = plot1, aes(x = yr, group = m, y = Rel, linetype = M), size = 1) +
      geom_line(data = STARel, aes(x = yr, y = Rel))+
      xlim(0,500)+
      labs(y = "Reliability", x = "Planning Horizon (yr)") +
      ggtitle(text)+
      theme_bw() +
      theme(plot.background = element_blank())+
      theme(legend.position = c(0.8,0.6))+
      theme(plot.title = element_text(size = 11, hjust = 0.5))
      z <- z + 1
  }
}

pp <- do.call(grid.arrange, p)
ggsave("C:/Users/tas48127/Desktop/Gen_Plot.pdf", pp, width = 10, height = 10)


plot1 <- results[!is.na(results$yr),]
ggplot(data = plot1, aes(x = yr, y = Rel, group = m)) +
  geom_line(aes(col = m), size = 1) 
  # labs(y = "Top Width Error (m)", x = "Station (km)", color = "Quantile")  +
  # scale_y_continuous(limits = c(-250,250), breaks = seq(-250,250, 50)) +
  # theme_bw() +
  # theme(plot.background = element_blank(),
  #       panel.grid.major = element_blank(),
  #       panel.grid.minor = element_blank())#,





#Evaluate integral for nonstationary case
##############################################################################################################
#Trapezoidal rule integration fucntion to solve reliability integral
##############################################################################################################

#r is resistance
#sdy is lognormal sd
#muy lognormal mean
#mut is Gumbel mean of annual flood peaks
#sdt is Gumbel sd of annual flood peaks
#n is the delta r (dr) for integration
integrand <- function(sdy, muy, mut, sdt,n) {
  diff <- 1
  ii <- 0
  while (diff > 0.00001){
    ii <- ii + 1
    if (ii == 1){
      r <- seq(1,90000,n)
      fr <- (1/(r*sdy*sqrt(2*pi)))*exp(-0.5*(((log(r)-muy)/sdy)^2))
      Fl <- exp(-exp(-((r-mut)/sdt)))
      Flfr <- fr*Fl
       
      #Flfr <- (1/(sqrt(2*pi)*sdy*r)) * exp((-0.5*(((log(r)-muy)/sdy)^2))) * exp(-exp(-((r-mut)/sdt)))
      a <- Flfr[1:(length(Flfr)-1)]
      b <- Flfr[2:length(Flfr)]
      h <- ((b+a)/2)*n
      approx <- sum(h)
    } else {
      r <- seq(10,9*(10^(ii+3)),n)
      Flfr <- (1/(sqrt(2*pi)*sdy*r)) * exp((-0.5*(((log(r)-muy)/sdy)^2))) * exp(-exp(-((r-mut)/sdt)))
      a <- Flfr[1:(length(Flfr)-1)]
      b <- Flfr[2:length(Flfr)]
      h <- ((b+a)/2)*n
      diff <- sum(h) - approx
      approx <- sum(h)
    }
  }
  return(approx)
}


