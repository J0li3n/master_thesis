#---------------------------------------------------------------
#Jolien Ketelaar
#Master thesis
#---------------------------------------------------------------

#Activate required packages
library(MASS)
library(lavaan)

#Set working directory
#Specify here the working directory where the Research Archive is saved
setwd("U:/My Documents/Master thesis")
setwd("C:/Users/Jolien/Documents/Universiteit Utrecht")
#Load environment from Research Archive when preferred
load("Results.RData")

#All code below is annotated with explanations of what is done
#Code to produce the tables and figures in the thesis are annotated as such
#For example: (Figure 1 in thesis)

#SIMULATION STUDY
#---------------------------------------------------------------

#REQUIRED FUNCTIONS
#---------------------------------------------------------------

#Sphericity measure
#------------------

#Function to calculate the epsilon of a covariance matrix
#Note: This function is not in general form yet and only works when p = 3
epsilonf <- function(S){
  #Transformation matrix (p-1) * p: lower part of Helmert matrix
  C <- matrix(c(1/(sqrt(2)), -1/(sqrt(2)), 0,
                1/(sqrt(6)), 1/(sqrt(6)), -2/(sqrt(6))), nrow = 2, byrow=T)
  #Calculate transformed covariance matrix
  S_star <- C%*%S%*%t(C)
  #Retrieve the eigenvalues of this matrix
  eigval <- eigen(S_star)$values
  
  #Measure of sphericity (Box's epsilon, see Cornell, Young, & Bratcher, 1990)
  epsilon <- (sum(eigval))^2 /
    ((p-1) * sum(eigval^2))
  list(epsilon = epsilon) #eigval = eigval)
}

#R transformation matrix
#-----------------------

#Function to make an R matrix given p for the MANOVA
Rmatrix <- function(p){
  Rg <- matrix(0, nrow = (p-1), ncol = p)
  
  for(i in 1:(p-1)){
    Rg[i, i] <- 1
    Rg[i, i+1] <- -1
  }
  Rg
}

#ANOVA function for within-subjects and between-subjects factor (Gives value for sphericity assumed, no corrections)
ANOVA <- function(X, p, n){
  #Convert to a data frame
  data <- as.data.frame(X)
  #Reshape the data into a long format
  #Specify that 1:3 are the different measurements, name of the variable and the time variable
  dataL <- reshape(data = data, varying = 2:(p+1), v.names = "outcome", timevar = "time", direction = "long")
  
  #Performing the repeated measures
  test <- (aov(dataL$outcome ~ factor(dataL$V1) + factor(dataL$time) + factor(dataL$V1)*factor(dataL$time) + Error(factor(dataL$id))))
  #Retrieve the p-value from the function
  summary(test)$ "Error: Within" [[1]] $"Pr" [[2]]
}

#MANOVA function for within-subjects and between-subjects (two-groups)
#Greenhouse, Geisser (1959) Profile data, p. 107
MANOVA <- function(X, p, n){
  #Calculate the covariance matrix of the data
  S1 <- (n[[1]]-1)*cov(X[X[,1]==1,2:(p+1)])/n[[1]]
  S2 <- (n[[2]]-1)*cov(X[X[,1]==2,2:(p+1)])/n[[2]]
  #Calculate the pooled sample covariance matrix
  Spool <- ((n[[1]])*S1 + (n[[2]])*S2) / (n[[1]] + n[[2]] -2)
  #Calculate x bar of the data
  xb1 <- matrix(apply(X[X[,1]==1,2:(p+1)],2,mean), p, 1)
  #For group 2
  xb2 <- matrix(apply(X[X[,1]==2,2:(p+1)],2,mean), p, 1)
  #Combine in a list
  xb_j <- list(xb1, xb2)
  #To transfer the matrices into difference scores
  #Function to make an R matrix given p
  R <- Rmatrix(p)
  #Calculate the T2 with the R and pooled covariance matrix
  T2 <- ((n[[1]]*n[[2]])/(n[[1]]+n[[2]]))*t(xb_j[[1]] - xb_j[[2]])%*%t(R)%*%solve(R%*%Spool%*%t(R))%*%R%*%(xb_j[[1]] - xb_j[[2]])
  #Calculate the F value
  Fvalue <- (n[[1]] + n[[2]] - p) / ((n[[1]]+n[[2]]-2)*(p-1)) * T2
  #Calculate the p-value
  1-pf(Fvalue,p-1,(n[[1]]+n[[2]]-p))
}


#GLR within-subjects and between-subjects
#Using the lavaan package
GLR <- function(X, p, n){
  D <- as.data.frame(X)
  names(D)<-c("g","X1","X2", "X3")
  
  # het gesatureerde model
  #eerst groep, dan tijdstip
  model1<-"X1~c(m11,m21)*1
  X2~c(m12,m22)*1
  X3~c(m13,m23)*1
  X1~~NA*X2
  X1~~NA*X3
  X2~~NA*X3
  "
  fit1<-sem(model1,data=D,group="g")
  
  # het gerestricteerde model (geen interactie)
  model2<-"X1~c(m11,m21)*1
  X2~c(m12,m22)*1
  X3~c(m13,m23)*1
  X1~~NA*X2
  X1~~NA*X3
  X2~~NA*X3
  m22 == m12 + m21 - m11
  m23 == m13 + m21 - m11
  "
  
  fit2<-sem(model2,data=D,group="g")
  
  X2 <- -2*(logLik(fit2)[[1]] - logLik(fit1)[[1]])
  1-pchisq(X2,2)
}


#The function for the simulation
#-------------------------------

#Note: This simulation is now based on an alpha level of .05
#This simulation only works with two groups
#Entrance for p should be scalar, n list, mu list, sigma list, sim scalar
simulation <- function(p, n, mu, Sigma, sim) {
  #Empty vectors for the pvalues of the GLR, ANOVA and the MANOVA
  GLRres <- c()
  MANres <- c()
  ANOres <- c()
  
  #Empty vectors that indicate later whether there were significant results
  fpGLR <- c()
  fpMAN <- c()
  fpANO <- c()
  
  
  #Simulation to calculate the robustness and power for the GLR, ANOVA and MANOVA
  for(i in 1:sim){
    #Sample data from a multinormal distribution with specified mu and Sigma
    X <- matrix(NA, (n[[1]]+n[[2]]), (p+1))
    X[1:n[[1]],2:(p+1)] <- mvrnorm(n[[1]],mu[[1]],Sigma[[1]])
    X[(n[[1]]+1):(n[[1]]+n[[2]]),2:(p+1)] <- mvrnorm(n[[2]],mu[[2]],Sigma[[2]])
    
    X[,1] <- c(rep(1,n[[1]]), rep(2, n[[2]]))

    #Perform the different tests
    GLRres[i] <- GLR(X, p, n)
    MANres[i] <- MANOVA(X, p, n)
    ANOres[i] <- ANOVA(X, p, n)
    #Count the number of significant and insignificant results for the different tests
    if(GLRres[i] <= 0.05){
      fpGLR[i] <- 1
    }
    if(GLRres[i] > 0.05){
      fpGLR[i] <- 0
    }
    if(ANOres[i] <= 0.05){
      fpANO[i] <- 1
    }
    if(ANOres[i] > 0.05){
      fpANO[i] <- 0
    }
    if(MANres[i] <= 0.05){
      fpMAN[i] <- 1
    }
    if(MANres[i] > 0.05){
      fpMAN[i] <- 0
    }
  }
  
  
  #Calculate the robustness or power
  list(GLR = mean(fpGLR), ANO =  mean(fpANO), MAN = mean(fpMAN))
}


#PERFORMING THE SIMULATION
#---------------------------------------------------------------

#Specifying the different conditions
#-----------------------------------
sim <- 1000
#Number of measurements
p <- 3
#Number of groups
k <- 2
#Determine the n's of the groups
n1 <- 25
n2 <- 50
n3 <- 100
n4 <- 200

n25 <- 25
n50 <- 50
n75 <- 75
n100 <- 100
n125 <- 125
n150 <- 150
n175 <- 175
n200 <- 200

#Population values
m <- 0
a <- c(-1,0,1)
b <- c(-1, 1)
gam1 <- c(-0.1, 0, 0.1)
gam2 <- c(0.1, 0, -0.1)

mu1<-rep(m,p)+a+rep(b[1],p)
mu2<-rep(m,p)+a+rep(b[2],p)
mu1i <- rep(m,p)+a+rep(b[1],p) + gam1 #met interactie effect
mu2i <- rep(m,p)+a+rep(b[2],p) + gam2 #met interactie effect

#Setting the different conditions for the covariance matrix 
Sigma1 <- matrix(c(1,0.5,0.5,
                   0.5,1,0.5,
                   0.5,0.5,1),nrow=3)
Sigma2 <- matrix(c(1,0.355,0.5,
                   0.355,1,0.645,
                   0.5,0.645,1),nrow=3)
Sigma3 <- matrix(c(1,0.285,0.5,
                   0.285,1,0.715,
                   0.5,0.715,1),nrow=3)
Sigma4 <- matrix(c(1,0.215,0.5,
                   0.215,1,0.785,
                   0.5,0.785,1),nrow=3)
Sigma5 <- matrix(c(1,0.15,0.5,
                   0.15,1,0.85,
                   0.5,0.85,1),nrow=3)

#Check Box's epsilon for the specified Sigma's
epsilonf(Sigma1) #no violation of sphericity, 1
epsilonf(Sigma2) # 0.9
epsilonf(Sigma3) #moderate violation of sphericity, 0.8
epsilonf(Sigma4) #0.7
epsilonf(Sigma5) #0.6, (0.5396, the minimum where it is still positive (semi-)definite)

#Making sure they are positive (semi-)definite.
eigen(Sigma1)
eigen(Sigma2)
eigen(Sigma3)
eigen(Sigma4)
eigen(Sigma5)

#Plot functions
#--------------

#Function to make plots for robustness as a function of the sphericity violations
plotfunction <- function(sim1, sim2, sim3, sim4, sim5) {
  plot1 <- matrix(c(sim1, sim2, sim3, sim4, sim5), nrow = 5, ncol = 3, byrow=T)
  plot2 <- cbind(c(epsilonf(Sigma1)[[1]],
                   epsilonf(Sigma2)[[1]],
                   epsilonf(Sigma3)[[1]],
                   epsilonf(Sigma4)[[1]],
                   epsilonf(Sigma5)[[1]]), plot1)
  
  plot(plot2[,1], plot2[,2], type = "l", col="black", lwd = 1.5, ylim=c(0.02,0.080),
       xlab = "epsilon", ylab = "rejection rate", axes = FALSE)
  axis(side=1, at=seq(0.6, 1, by=0.10))
  axis(side=2, at=seq(0.02, 0.08, by=0.01))
  lines(plot2[,1], plot2[,3], col = "gray30", lty = 2, lwd = 1.5)
  lines(plot2[,1], plot2[,4], col = "gray60", lty = 3, lwd = 1.5)
  legend(0.9, 0.038, legend=c("GLR", "ANOVA", "MANOVA"),
         col=c("black", "gray30", "gray60"), lty=c(1, 2, 3), lwd = 1.5, cex=0.8)
}

#Function to make plots for the power as a function of the sphericity violations
plotPowerfunction <- function(sim1, sim2, sim3, sim4, sim5) {
  plot1 <- matrix(c(sim1, sim2, sim3, sim4, sim5), nrow = 5, ncol = 3, byrow=T)
  plot2 <- cbind(c(epsilonf(Sigma1)[[1]],
                   epsilonf(Sigma2)[[1]],
                   epsilonf(Sigma3)[[1]],
                   epsilonf(Sigma4)[[1]],
                   epsilonf(Sigma5)[[1]]), plot1)
  
  plot(plot2[,1], plot2[,2], type = "l", col="black", lwd = 1.5, ylim=c(0.10,1.00), xlab = "epsilon", ylab = "detection rate", axes = FALSE)
  axis(side=1, at=seq(0.6, 1, by=0.10))
  axis(side=2, at=seq(0.00, 1.00, by=0.2))
  lines(plot2[,1], plot2[,3], col = "gray30", lty = 2, lwd = 1.5)
  lines(plot2[,1], plot2[,4], col = "gray60", lty = 3, lwd = 1.5)
  legend(0.9, 0.35, legend=c("GLR", "ANOVA", "MANOVA"),
         col=c("black", "gray30", "gray60"), lty=c(1, 2, 3), lwd = 1.5, cex=0.8)
}


#Plot function for robustness for sample size differences, N should be a vector of the sample sizes
plotSampleSize200 <- function(sim1, sim2, sim3, sim4, sim5, sim6, sim7, sim8, N){
  plot1AA <- matrix(c(sim1, sim2, sim3, sim4, sim5, sim6, sim7, sim8), nrow = 8, ncol = 3, byrow=T)
  plot2AA <- cbind(N, plot1AA)
  
  plot(plot2AA[,1], plot2AA[,2], type = "l", col="black", lwd = 1.5, ylim=c(0.02,0.08), xlab = "sample size group 1", ylab = "rejection rate", axes = FALSE)
  axis(side=1, at=seq(25, 200, by=25))
  axis(side=2, at=seq(0.02, 0.08, by=0.01))
  lines(plot2AA[,1], plot2AA[,3], col = "gray30", lty = 2, lwd = 1.5)
  lines(plot2AA[,1], plot2AA[,4], col = "gray60", lty = 3, lwd = 1.5)
  legend(155, 0.038, legend=c("GLR", "ANOVA", "MANOVA"),
         col=c("black", "gray30", "gray60"), lty=c(1, 2, 3), lwd = 1.5, cex=0.8)
}

#Plot function for power for sample size differences of the power, N should be a vector of the sample sizes
plotPowerSampleSize200 <- function(sim1, sim2, sim3, sim4, sim5, sim6, sim7, sim8, N){
  plot1AA <- matrix(c(sim1, sim2, sim3, sim4, sim5, sim6, sim7, sim8), nrow = 8, ncol = 3, byrow=T)
  plot2AA <- cbind(N, plot1AA)
  plot(plot2AA[,1], plot2AA[,2], type = "l", col="black", lwd = 1.5, ylim=c(0.10,1.00), xlab = "sample size group 1", ylab = "detection rate", axes = FALSE)
  axis(side=1, at=seq(25, 200, by=25))
  axis(side=2, at=seq(0.00, 1.00, by=0.2))
  lines(plot2AA[,1], plot2AA[,3], col = "gray30", lty = 2, lwd = 1.5)
  lines(plot2AA[,1], plot2AA[,4], col = "gray60", lty = 3, lwd = 1.5)
  legend(155, 0.35, legend=c("GLR", "ANOVA", "MANOVA"),
         col=c("black", "gray30", "gray60"), lty=c(1, 2, 3), lwd = 1.5, cex=0.8)
}

#Round function for the matrix of the results
#And to delete the list inside
RoundMatrix <- function(M, digits = 2){
  MR <- matrix(NA, nrow(M), ncol(M))
  
  for(r in 1:nrow(M)){
    for(c in 1:ncol(M)){
      MR[r, c] <- round(M[r, c][[1]], digits = digits)
    }
  }
  MR
}


#Make lists of the combinations
#------------------------------
#Ns
n.1.1 <- list(n1, n1)
n.2.2 <- list(n2, n2)
n.3.3 <- list(n3, n3)
n.4.4 <- list(n4, n4)

#Unequal group sizes
n.1.4 <- list(n1, n4)
n.2.4 <- list(n2, n4)
n.3.4 <- list(n3, n4)

#Means
mu.1.1 <- list(mu1, mu1)
mu.1.2 <- list(mu1, mu2)

#Equal covariance matrices per group
Sigma.1.1 <- list(Sigma1, Sigma1)
Sigma.2.2 <- list(Sigma2, Sigma2)
Sigma.3.3 <- list(Sigma3, Sigma3)
Sigma.4.4 <- list(Sigma4, Sigma4)
Sigma.5.5 <- list(Sigma5, Sigma5)

#Running the simulations
#-----------------------

#Robustness
#----------

#With equal sample sizes and no violation
set.seed(221)
sim221 <- simulation(p, list(n25, n25), mu.1.2, Sigma.1.1, 1000)
set.seed(222)
sim222 <- simulation(p, list(n50, n50), mu.1.2, Sigma.1.1, 1000)
set.seed(223)
sim223 <- simulation(p, list(n75, n75), mu.1.2, Sigma.1.1, 1000)
set.seed(224)
sim224 <- simulation(p, list(n100, n100), mu.1.2, Sigma.1.1, 1000)
set.seed(225)
sim225 <- simulation(p, list(n125, n125), mu.1.2, Sigma.1.1, 1000)
set.seed(226)
sim226 <- simulation(p, list(n150, n150), mu.1.2, Sigma.1.1, 1000)
set.seed(227)
sim227 <- simulation(p, list(n175, n175), mu.1.2, Sigma.1.1, 1000)
set.seed(228)
sim228 <- simulation(p, list(n200, n200), mu.1.2, Sigma.1.1, 1000)

#Plot for the robustness with no violations (Figure 1 in thesis)
plot1AA <- matrix(c(sim221, sim222, sim223, sim224, sim225, sim226, sim227, sim228), nrow = 8, ncol = 3, byrow=T)
plot2AA <- cbind(c(25, 50, 75, 100, 125, 150, 175, 200), plot1AA)

plot(plot2AA[,1], plot2AA[,2], type = "l", col="black", lwd = 1.5, ylim=c(0.02,0.08), xlab = "sample size in both groups", ylab = "rejection rate", axes = FALSE)
axis(side=1, at=seq(25, 200, by=25))
axis(side=2, at=seq(0.02, 0.08, by=0.01))
lines(plot2AA[,1], plot2AA[,3], col = "gray30", lty = 2, lwd = 1.5)
lines(plot2AA[,1], plot2AA[,4], col = "gray60", lty = 3, lwd = 1.5)
legend(155, 0.038, legend=c("GLR", "ANOVA", "MANOVA"),
       col=c("black", "gray30", "gray60"), lty=c(1, 2, 3), lwd = 1.5, cex=0.8)

#n = 200 both groups, sphericity violation the same in both groups, mu's different
set.seed(61)
sim61 <- simulation(p, n.4.4, mu.1.2, Sigma.1.1, 1000) 
set.seed(62)
sim62 <- simulation(p, n.4.4, mu.1.2, Sigma.2.2, 1000)
set.seed(63)
sim63 <- simulation(p, n.4.4, mu.1.2, Sigma.3.3, 1000) 
set.seed(64)
sim64 <- simulation(p, n.4.4, mu.1.2, Sigma.4.4, 1000) 
set.seed(65)
sim65 <- simulation(p, n.4.4, mu.1.2, Sigma.5.5, 1000)

plotfunction(sim61, sim62, sim63, sim64, sim65) #(Figure 2 in thesis)

#n = 100 both groups, sphericity violation the same in both groups, mu's different
set.seed(71)
sim71 <- simulation(p, n.3.3, mu.1.2, Sigma.1.1, 1000) 
set.seed(72)
sim72 <- simulation(p, n.3.3, mu.1.2, Sigma.2.2, 1000)
set.seed(73)
sim73 <- simulation(p, n.3.3, mu.1.2, Sigma.3.3, 1000) 
set.seed(74)
sim74 <- simulation(p, n.3.3, mu.1.2, Sigma.4.4, 1000) 
set.seed(75)
sim75 <- simulation(p, n.3.3, mu.1.2, Sigma.5.5, 1000)

plotfunction(sim71, sim72, sim73, sim74, sim75)


#n = 50 both groups, sphericity violation the same in both groups, mu's different
set.seed(81)
sim81 <- simulation(p, n.2.2, mu.1.2, Sigma.1.1, 1000) 
set.seed(82)
sim82 <- simulation(p, n.2.2, mu.1.2, Sigma.2.2, 1000)
set.seed(83)
sim83 <- simulation(p, n.2.2, mu.1.2, Sigma.3.3, 1000) 
set.seed(84)
sim84 <- simulation(p, n.2.2, mu.1.2, Sigma.4.4, 1000) 
set.seed(85)
sim85 <- simulation(p, n.2.2, mu.1.2, Sigma.5.5, 1000)

plotfunction(sim81, sim82, sim83, sim84, sim85)


#n = 25 both groups, sphericity violation the same in both groups, mu's different
set.seed(91)
sim91 <- simulation(p, n.1.1, mu.1.2, Sigma.1.1, 1000) 
set.seed(92)
sim92 <- simulation(p, n.1.1, mu.1.2, Sigma.2.2, 1000)
set.seed(93)
sim93 <- simulation(p, n.1.1, mu.1.2, Sigma.3.3, 1000) 
set.seed(94)
sim94 <- simulation(p, n.1.1, mu.1.2, Sigma.4.4, 1000) 
set.seed(95)
sim95 <- simulation(p, n.1.1, mu.1.2, Sigma.5.5, 1000)

plotfunction(sim91, sim92, sim93, sim94, sim95)

#Making a table for sphericity violated

SpViol200 <- matrix(c(sim61, sim62, sim63, sim64, sim65), ncol = 3, byrow=T)
SpViol100 <- matrix(c(sim71, sim72, sim73, sim74, sim75), ncol = 3, byrow=T)
SpViol50 <- matrix(c(sim81, sim82, sim83, sim84, sim85), ncol = 3, byrow=T)
SpViol25 <- matrix(c(sim91, sim92, sim93, sim94, sim95), ncol = 3, byrow=T)

SpViol <- cbind(SpViol25, SpViol50, SpViol100, SpViol200)

SpViolR <- RoundMatrix(SpViol, 2) #(Table 6 in thesis)

#Homogeneity violated

#0.5Sigma
set.seed(0411)
sim0411 <- simulation(p, list(n25, n200), mu.1.2, list((1/2)*Sigma1, Sigma1), 1000)
set.seed(0412)
sim0412 <- simulation(p, list(n50, n200), mu.1.2, list((1/2)*Sigma1, Sigma1), 1000)
set.seed(0413)
sim0413 <- simulation(p, list(n75, n200), mu.1.2, list((1/2)*Sigma1, Sigma1), 1000)
set.seed(0414)
sim0414 <- simulation(p, list(n100, n200), mu.1.2, list((1/2)*Sigma1, Sigma1), 1000)
set.seed(0415)
sim0415 <- simulation(p, list(n125, n200), mu.1.2, list((1/2)*Sigma1, Sigma1), 1000)
set.seed(0416)
sim0416 <- simulation(p, list(n150, n200), mu.1.2, list((1/2)*Sigma1, Sigma1), 1000)
set.seed(0417)
sim0417 <- simulation(p, list(n175, n200), mu.1.2, list((1/2)*Sigma1, Sigma1), 1000)
set.seed(0418)
sim0418 <- simulation(p, list(n200, n200), mu.1.2, list((1/2)*Sigma1, Sigma1), 1000)

plotSampleSize200(sim0411, sim0412, sim0413, sim0414, sim0415, sim0416, sim0417, sim0418, c(25, 50, 75, 100, 125, 150, 175, 200))
Homogen0.5 <- rbind(sim0411, sim0412, sim0413, sim0414, sim0415, sim0416, sim0417, sim0418)

#0.75Sigma
set.seed(0511)
sim0511 <- simulation(p, list(n25, n200), mu.1.2, list((3/4)*Sigma1, Sigma1), 1000)
set.seed(0512)
sim0512 <- simulation(p, list(n50, n200), mu.1.2, list((3/4)*Sigma1, Sigma1), 1000)
set.seed(0513)
sim0513 <- simulation(p, list(n75, n200), mu.1.2, list((3/4)*Sigma1, Sigma1), 1000)
set.seed(0514)
sim0514 <- simulation(p, list(n100, n200), mu.1.2, list((3/4)*Sigma1, Sigma1), 1000)
set.seed(0515)
sim0515 <- simulation(p, list(n125, n200), mu.1.2, list((3/4)*Sigma1, Sigma1), 1000)
set.seed(0516)
sim0516 <- simulation(p, list(n150, n200), mu.1.2, list((3/4)*Sigma1, Sigma1), 1000)
set.seed(0517)
sim0517 <- simulation(p, list(n175, n200), mu.1.2, list((3/4)*Sigma1, Sigma1), 1000)
set.seed(0518)
sim0518 <- simulation(p, list(n200, n200), mu.1.2, list((3/4)*Sigma1, Sigma1), 1000)

plotSampleSize200(sim0511, sim0512, sim0513, sim0514, sim0515, sim0516, sim0517, sim0518, c(25, 50, 75, 100, 125, 150, 175, 200))
Homogen0.75<- rbind(sim0511, sim0512, sim0513, sim0514, sim0515, sim0516, sim0517, sim0518)

#1.25Sigma
set.seed(0521)
sim0521 <- simulation(p, list(n25, n200), mu.1.2, list((1.25)*Sigma1, Sigma1), 1000)
set.seed(0522)
sim0522 <- simulation(p, list(n50, n200), mu.1.2, list((1.25)*Sigma1, Sigma1), 1000)
set.seed(0523)
sim0523 <- simulation(p, list(n75, n200), mu.1.2, list((1.25)*Sigma1, Sigma1), 1000)
set.seed(0524)
sim0524 <- simulation(p, list(n100, n200), mu.1.2, list((1.25)*Sigma1, Sigma1), 1000)
set.seed(0525)
sim0525 <- simulation(p, list(n125, n200), mu.1.2, list((1.25)*Sigma1, Sigma1), 1000)
set.seed(0526)
sim0526 <- simulation(p, list(n150, n200), mu.1.2, list((1.25)*Sigma1, Sigma1), 1000)
set.seed(0527)
sim0527 <- simulation(p, list(n175, n200), mu.1.2, list((1.25)*Sigma1, Sigma1), 1000)
set.seed(0528)
sim0528 <- simulation(p, list(n200, n200), mu.1.2, list((1.25)*Sigma1, Sigma1), 1000)

plotSampleSize200(sim0521, sim0522, sim0523, sim0524, sim0525, sim0526, sim0527, sim0528, c(25, 50, 75, 100, 125, 150, 175, 200))
Homogen1.25<- rbind(sim0521, sim0522, sim0523, sim0524, sim0525, sim0526, sim0527, sim0528)

#1.50Sigma
set.seed(0611)
sim0611 <- simulation(p, list(n25, n200), mu.1.2, list((1.5)*Sigma1, Sigma1), 1000)
set.seed(0612)
sim0612 <- simulation(p, list(n50, n200), mu.1.2, list((1.5)*Sigma1, Sigma1), 1000)
set.seed(0613)
sim0613 <- simulation(p, list(n75, n200), mu.1.2, list((1.5)*Sigma1, Sigma1), 1000)
set.seed(0614)
sim0614 <- simulation(p, list(n100, n200), mu.1.2, list((1.5)*Sigma1, Sigma1), 1000)
set.seed(0615)
sim0615 <- simulation(p, list(n125, n200), mu.1.2, list((1.5)*Sigma1, Sigma1), 1000)
set.seed(0616)
sim0616 <- simulation(p, list(n150, n200), mu.1.2, list((1.5)*Sigma1, Sigma1), 1000)
set.seed(0617)
sim0617 <- simulation(p, list(n175, n200), mu.1.2, list((1.5)*Sigma1, Sigma1), 1000)
set.seed(0618)
sim0618 <- simulation(p, list(n200, n200), mu.1.2, list((1.5)*Sigma1, Sigma1), 1000)

plotSampleSize200(sim0611, sim0612, sim0613, sim0614, sim0615, sim0616, sim0617, sim0618, c(25, 50, 75, 100, 125, 150, 175, 200))
Homogen1.5 <- rbind(sim0611, sim0612, sim0613, sim0614, sim0615, sim0616, sim0617, sim0618)

#1.75Sigma
set.seed(0621)
sim0621 <- simulation(p, list(n25, n200), mu.1.2, list((1.75)*Sigma1, Sigma1), 1000)
set.seed(0629) #Different seed because of covergence problems
sim0622 <- simulation(p, list(n50, n200), mu.1.2, list((1.75)*Sigma1, Sigma1), 1000)
set.seed(0623)
sim0623 <- simulation(p, list(n75, n200), mu.1.2, list((1.75)*Sigma1, Sigma1), 1000)
set.seed(0624)
sim0624 <- simulation(p, list(n100, n200), mu.1.2, list((1.75)*Sigma1, Sigma1), 1000)
set.seed(0625)
sim0625 <- simulation(p, list(n125, n200), mu.1.2, list((1.75)*Sigma1, Sigma1), 1000)
set.seed(0626)
sim0626 <- simulation(p, list(n150, n200), mu.1.2, list((1.75)*Sigma1, Sigma1), 1000)
set.seed(0627)
sim0627 <- simulation(p, list(n175, n200), mu.1.2, list((1.75)*Sigma1, Sigma1), 1000)
set.seed(0628)
sim0628 <- simulation(p, list(n200, n200), mu.1.2, list((1.75)*Sigma1, Sigma1), 1000)

plotSampleSize200(sim0621, sim0622, sim0623, sim0624, sim0625, sim0626, sim0627, sim0628, c(25, 50, 75, 100, 125, 150, 175, 200))
Homogen1.75 <- rbind(sim0621, sim0622, sim0623, sim0624, sim0625, sim0626, sim0627, sim0628)

#2.00Sigma
set.seed(0711)
sim0711 <- simulation(p, list(n25, n200), mu.1.2, list((2)*Sigma1, Sigma1), 1000)
set.seed(0712)
sim0712 <- simulation(p, list(n50, n200), mu.1.2, list((2)*Sigma1, Sigma1), 1000)
set.seed(0713)
sim0713 <- simulation(p, list(n75, n200), mu.1.2, list((2)*Sigma1, Sigma1), 1000)
set.seed(0714)
sim0714 <- simulation(p, list(n100, n200), mu.1.2, list((2)*Sigma1, Sigma1), 1000)
set.seed(0715)
sim0715 <- simulation(p, list(n125, n200), mu.1.2, list((2)*Sigma1, Sigma1), 1000)
set.seed(0716)
sim0716 <- simulation(p, list(n150, n200), mu.1.2, list((2)*Sigma1, Sigma1), 1000)
set.seed(0717)
sim0717 <- simulation(p, list(n175, n200), mu.1.2, list((2)*Sigma1, Sigma1), 1000)
set.seed(0718)
sim0718 <- simulation(p, list(n200, n200), mu.1.2, list((2)*Sigma1, Sigma1), 1000)

plotSampleSize200(sim0711, sim0712, sim0713, sim0714, sim0715, sim0716, sim0717, sim0718, c(25, 50, 75, 100, 125, 150, 175, 200))
Homogen2 <- rbind(sim0711, sim0712, sim0713, sim0714, sim0715, sim0716, sim0717, sim0718)

#Make a matrix for all the homogeneity results
Homogen <- cbind(Homogen0.5, Homogen0.75, Homogen1.25, Homogen1.5, Homogen1.75, Homogen2)
#Make separate matrices for smaller and larger covariance matrices
HomogenA <- cbind(Homogen0.5, Homogen0.75)
HomogenB <- cbind(Homogen1.25, Homogen1.5, Homogen1.75, Homogen2)

#Round the matrices
HomogenR <- RoundMatrix(Homogen)
HomogenAR <- RoundMatrix(HomogenA) #(Table 7 in thesis)
HomogenBR <- RoundMatrix(HomogenB) #(Table 8 in thesis)

#Both violated

#0.75 Sigma, 0.8 epsilon
set.seed(0811)
sim0811 <- simulation(p, list(n25, n200), mu.1.2, list((3/4)*Sigma3, Sigma3), 1000)
set.seed(0812)
sim0812 <- simulation(p, list(n50, n200), mu.1.2, list((3/4)*Sigma3, Sigma3), 1000)
set.seed(0813)
sim0813 <- simulation(p, list(n75, n200), mu.1.2, list((3/4)*Sigma3, Sigma3), 1000)
set.seed(0814)
sim0814 <- simulation(p, list(n100, n200), mu.1.2, list((3/4)*Sigma3, Sigma3), 1000)
set.seed(0815)
sim0815 <- simulation(p, list(n125, n200), mu.1.2, list((3/4)*Sigma3, Sigma3), 1000)
set.seed(0816)
sim0816 <- simulation(p, list(n150, n200), mu.1.2, list((3/4)*Sigma3, Sigma3), 1000)
set.seed(0817)
sim0817 <- simulation(p, list(n175, n200), mu.1.2, list((3/4)*Sigma3, Sigma3), 1000)
set.seed(0818)
sim0818 <- simulation(p, list(n200, n200), mu.1.2, list((3/4)*Sigma3, Sigma3), 1000)

plotSampleSize200(sim0811, sim0812, sim0813, sim0814, sim0815, sim0816, sim0817, sim0818, c(25, 50, 75, 100, 125, 150, 175, 200))
Homogen0.75Sph0.8 <- rbind(sim0811, sim0812, sim0813, sim0814, sim0815, sim0816, sim0817, sim0818)

#1.25 Sigma, 0.8 epsilon
set.seed(0821)
sim0821 <- simulation(p, list(n25, n200), mu.1.2, list((1.25)*Sigma3, Sigma3), 1000)
set.seed(0822)
sim0822 <- simulation(p, list(n50, n200), mu.1.2, list((1.25)*Sigma3, Sigma3), 1000)
set.seed(0823)
sim0823 <- simulation(p, list(n75, n200), mu.1.2, list((1.25)*Sigma3, Sigma3), 1000)
set.seed(0824)
sim0824 <- simulation(p, list(n100, n200), mu.1.2, list((1.25)*Sigma3, Sigma3), 1000)
set.seed(0825)
sim0825 <- simulation(p, list(n125, n200), mu.1.2, list((1.25)*Sigma3, Sigma3), 1000)
set.seed(0826)
sim0826 <- simulation(p, list(n150, n200), mu.1.2, list((1.25)*Sigma3, Sigma3), 1000)
set.seed(0827)
sim0827 <- simulation(p, list(n175, n200), mu.1.2, list((1.25)*Sigma3, Sigma3), 1000)
set.seed(0828)
sim0828 <- simulation(p, list(n200, n200), mu.1.2, list((1.25)*Sigma3, Sigma3), 1000)

plotSampleSize200(sim0821, sim0822, sim0823, sim0824, sim0825, sim0826, sim0827, sim0828, c(25, 50, 75, 100, 125, 150, 175, 200))
Homogen1.25Sph0.8 <- rbind(sim0821, sim0822, sim0823, sim0824, sim0825, sim0826, sim0827, sim0828)

#1.50 Sigma, 0.8 epsilon
set.seed(0851)
sim0851 <- simulation(p, list(n25, n200), mu.1.2, list((1.5)*Sigma3, Sigma3), 1000)
set.seed(0852)
sim0852 <- simulation(p, list(n50, n200), mu.1.2, list((1.5)*Sigma3, Sigma3), 1000)
set.seed(0853)
sim0853 <- simulation(p, list(n75, n200), mu.1.2, list((1.5)*Sigma3, Sigma3), 1000)
set.seed(0854)
sim0854 <- simulation(p, list(n100, n200), mu.1.2, list((1.5)*Sigma3, Sigma3), 1000)
set.seed(0855)
sim0855 <- simulation(p, list(n125, n200), mu.1.2, list((1.5)*Sigma3, Sigma3), 1000)
set.seed(0856)
sim0856 <- simulation(p, list(n150, n200), mu.1.2, list((1.5)*Sigma3, Sigma3), 1000)
set.seed(0857)
sim0857 <- simulation(p, list(n175, n200), mu.1.2, list((1.5)*Sigma3, Sigma3), 1000)
set.seed(0858)
sim0858 <- simulation(p, list(n200, n200), mu.1.2, list((1.5)*Sigma3, Sigma3), 1000)

plotSampleSize200(sim0851, sim0852, sim0853, sim0854, sim0855, sim0856, sim0857, sim0858, c(25, 50, 75, 100, 125, 150, 175, 200))
Homogen1.5Sph0.8 <- rbind(sim0851, sim0852, sim0853, sim0854, sim0855, sim0856, sim0857, sim0858)

#0.75 Sigma, 0.7 epsilon
set.seed(0831)
sim0831 <- simulation(p, list(n25, n200), mu.1.2, list((3/4)*Sigma4, Sigma4), 1000)
set.seed(0832)
sim0832 <- simulation(p, list(n50, n200), mu.1.2, list((3/4)*Sigma4, Sigma4), 1000)
set.seed(0833)
sim0833 <- simulation(p, list(n75, n200), mu.1.2, list((3/4)*Sigma4, Sigma4), 1000)
set.seed(0834)
sim0834 <- simulation(p, list(n100, n200), mu.1.2, list((3/4)*Sigma4, Sigma4), 1000)
set.seed(0835)
sim0835 <- simulation(p, list(n125, n200), mu.1.2, list((3/4)*Sigma4, Sigma4), 1000)
set.seed(0836)
sim0836 <- simulation(p, list(n150, n200), mu.1.2, list((3/4)*Sigma4, Sigma4), 1000)
set.seed(0837)
sim0837 <- simulation(p, list(n175, n200), mu.1.2, list((3/4)*Sigma4, Sigma4), 1000)
set.seed(0838)
sim0838 <- simulation(p, list(n200, n200), mu.1.2, list((3/4)*Sigma4, Sigma4), 1000)

plotSampleSize200(sim0831, sim0832, sim0833, sim0834, sim0835, sim0836, sim0837, sim0838, c(25, 50, 75, 100, 125, 150, 175, 200))
Homogen0.75Sph0.7 <- rbind(sim0831, sim0832, sim0833, sim0834, sim0835, sim0836, sim0837, sim0838)

#1.25 Sigma, 0.7 epsilon
set.seed(0841)
sim0841 <- simulation(p, list(n25, n200), mu.1.2, list((1.25)*Sigma4, Sigma4), 1000)
set.seed(0842)
sim0842 <- simulation(p, list(n50, n200), mu.1.2, list((1.25)*Sigma4, Sigma4), 1000)
set.seed(0843)
sim0843 <- simulation(p, list(n75, n200), mu.1.2, list((1.25)*Sigma4, Sigma4), 1000)
set.seed(0844)
sim0844 <- simulation(p, list(n100, n200), mu.1.2, list((1.25)*Sigma4, Sigma4), 1000)
set.seed(0845)
sim0845 <- simulation(p, list(n125, n200), mu.1.2, list((1.25)*Sigma4, Sigma4), 1000)
set.seed(0846)
sim0846 <- simulation(p, list(n150, n200), mu.1.2, list((1.25)*Sigma4, Sigma4), 1000)
set.seed(0847)
sim0847 <- simulation(p, list(n175, n200), mu.1.2, list((1.25)*Sigma4, Sigma4), 1000)
set.seed(0848)
sim0848 <- simulation(p, list(n200, n200), mu.1.2, list((1.25)*Sigma4, Sigma4), 1000)

plotSampleSize200(sim0841, sim0842, sim0843, sim0844, sim0845, sim0846, sim0847, sim0848, c(25, 50, 75, 100, 125, 150, 175, 200)) #(Figure 3 in thesis)
Homogen1.25Sph0.7 <- rbind(sim0841, sim0842, sim0843, sim0844, sim0845, sim0846, sim0847, sim0848)

#1.50 Sigma, 0.7 epsilon
set.seed(0861)
sim0861 <- simulation(p, list(n25, n200), mu.1.2, list((1.5)*Sigma4, Sigma4), 1000)
set.seed(0862)
sim0862 <- simulation(p, list(n50, n200), mu.1.2, list((1.5)*Sigma4, Sigma4), 1000)
set.seed(0863)
sim0863 <- simulation(p, list(n75, n200), mu.1.2, list((1.5)*Sigma4, Sigma4), 1000)
set.seed(0864)
sim0864 <- simulation(p, list(n100, n200), mu.1.2, list((1.5)*Sigma4, Sigma4), 1000)
set.seed(0865)
sim0865 <- simulation(p, list(n125, n200), mu.1.2, list((1.5)*Sigma4, Sigma4), 1000)
set.seed(0866)
sim0866 <- simulation(p, list(n150, n200), mu.1.2, list((1.5)*Sigma4, Sigma4), 1000)
set.seed(0867)
sim0867 <- simulation(p, list(n175, n200), mu.1.2, list((1.5)*Sigma4, Sigma4), 1000)
set.seed(0868)
sim0868 <- simulation(p, list(n200, n200), mu.1.2, list((1.5)*Sigma4, Sigma4), 1000)

plotSampleSize200(sim0861, sim0862, sim0863, sim0864, sim0865, sim0866, sim0867, sim0868, c(25, 50, 75, 100, 125, 150, 175, 200))
Homogen1.5Sph0.7 <- rbind(sim0861, sim0862, sim0863, sim0864, sim0865, sim0866, sim0867, sim0868)

#Combine all the results of both assumptions violated
HomoSphV<- cbind(Homogen0.75Sph0.7, Homogen1.25Sph0.7, Homogen1.5Sph0.7,
                   Homogen0.75Sph0.8, Homogen1.25Sph0.8, Homogen1.5Sph0.8)

#Round the matrix
HomoSphVR <- RoundMatrix(HomoSphV) #(Table 9 in thesis)


#Power
#-----
#With equal sample sizes and no violation
set.seed(1001)
sim1001 <- simulation(p, list(n25, n25), list(mu1i, mu2i), Sigma.1.1, 1000)
set.seed(1002)
sim1002 <- simulation(p, list(n50, n50), list(mu1i, mu2i), Sigma.1.1, 1000)
set.seed(1003)
sim1003 <- simulation(p, list(n75, n75), list(mu1i, mu2i), Sigma.1.1, 1000)
set.seed(1004)
sim1004 <- simulation(p, list(n100, n100), list(mu1i, mu2i), Sigma.1.1, 1000)
set.seed(1005)
sim1005 <- simulation(p, list(n125, n125), list(mu1i, mu2i), Sigma.1.1, 1000)
set.seed(1006)
sim1006 <- simulation(p, list(n150, n150), list(mu1i, mu2i), Sigma.1.1, 1000)
set.seed(1007)
sim1007 <- simulation(p, list(n175, n175), list(mu1i, mu2i), Sigma.1.1, 1000)
set.seed(1008)
sim1008 <- simulation(p, list(n200, n200), list(mu1i, mu2i), Sigma.1.1, 1000)

#Plot for power without violations #(Figure 4 in thesis)
plot1AA <- matrix(c(sim1001, sim1002, sim1003, sim1004, sim1005, sim1006, sim1007, sim1008), nrow = 8, ncol = 3, byrow=T)
plot2AA <- cbind(c(25, 50, 75, 100, 125, 150, 175, 200), plot1AA)

plot(plot2AA[,1], plot2AA[,2], type = "l", col="black", lwd = 1.5, ylim=c(0.10,1.00), xlab = "sample size in both groups", ylab = "detection rate", axes = FALSE)
axis(side=1, at=seq(25, 200, by=25))
axis(side=2, at=seq(0.00, 1.00, by=0.2))
lines(plot2AA[,1], plot2AA[,3], col = "gray30", lty = 2, lwd = 1.5)
lines(plot2AA[,1], plot2AA[,4], col = "gray60", lty = 3, lwd = 1.5)
legend(155, 0.35, legend=c("GLR", "ANOVA", "MANOVA"),
       col=c("black", "gray30", "gray60"), lty=c(1, 2, 3), lwd = 1.5, cex=0.8)

#n = 200 both groups, sphericity violation the same in both groups, mu's different
set.seed(1011)
sim1011 <- simulation(p, n.4.4, list(mu1i, mu2i), Sigma.1.1, 1000) 
set.seed(1012)
sim1012 <- simulation(p, n.4.4, list(mu1i, mu2i), Sigma.2.2, 1000)
set.seed(1013)
sim1013 <- simulation(p, n.4.4, list(mu1i, mu2i), Sigma.3.3, 1000) 
set.seed(1014)
sim1014 <- simulation(p, n.4.4, list(mu1i, mu2i), Sigma.4.4, 1000) 
set.seed(1015)
sim1015 <- simulation(p, n.4.4, list(mu1i, mu2i), Sigma.5.5, 1000)

plotPowerfunction(sim1011, sim1012, sim1013, sim1014, sim1015)

#n = 100 both groups, sphericity violation the same in both groups, mu's different
set.seed(1021)
sim1021 <- simulation(p, n.3.3, list(mu1i, mu2i), Sigma.1.1, 1000) 
set.seed(1022)
sim1022 <- simulation(p, n.3.3, list(mu1i, mu2i), Sigma.2.2, 1000)
set.seed(1023)
sim1023 <- simulation(p, n.3.3, list(mu1i, mu2i), Sigma.3.3, 1000) 
set.seed(1024)
sim1024 <- simulation(p, n.3.3, list(mu1i, mu2i), Sigma.4.4, 1000) 
set.seed(1025)
sim1025 <- simulation(p, n.3.3, list(mu1i, mu2i), Sigma.5.5, 1000)

plotPowerfunction(sim1021, sim1022, sim1023, sim1024, sim1025)

#n = 50 both groups, sphericity violation the same in both groups, mu's different
set.seed(1031)
sim1031 <- simulation(p, n.2.2, list(mu1i, mu2i), Sigma.1.1, 1000) 
set.seed(1032)
sim1032 <- simulation(p, n.2.2, list(mu1i, mu2i), Sigma.2.2, 1000)
set.seed(1033)
sim1033 <- simulation(p, n.2.2, list(mu1i, mu2i), Sigma.3.3, 1000) 
set.seed(1034)
sim1034 <- simulation(p, n.2.2, list(mu1i, mu2i), Sigma.4.4, 1000) 
set.seed(1035)
sim1035 <- simulation(p, n.2.2, list(mu1i, mu2i), Sigma.5.5, 1000)

plotPowerfunction(sim1031, sim1032, sim1033, sim1034, sim1035)

#n = 25 both groups, sphericity violation the same in both groups, mu's different
set.seed(1041)
sim1041 <- simulation(p, n.1.1, list(mu1i, mu2i), Sigma.1.1, 1000) 
set.seed(1042)
sim1042 <- simulation(p, n.1.1, list(mu1i, mu2i), Sigma.2.2, 1000)
set.seed(1043)
sim1043 <- simulation(p, n.1.1, list(mu1i, mu2i), Sigma.3.3, 1000) 
set.seed(1044)
sim1044 <- simulation(p, n.1.1, list(mu1i, mu2i), Sigma.4.4, 1000) 
set.seed(1045)
sim1045 <- simulation(p, n.1.1, list(mu1i, mu2i), Sigma.5.5, 1000)

plotPowerfunction(sim1041, sim1042, sim1043, sim1044, sim1045)

#Make matrices of the results
SpViol200P <- matrix(c(sim1011, sim1012, sim1013, sim1014, sim1015), ncol = 3, byrow=T)
SpViol100P <- matrix(c(sim1021, sim1022, sim1023, sim1024, sim1025), ncol = 3, byrow=T)
SpViol50P <- matrix(c(sim1031, sim1032, sim1033, sim1034, sim1035), ncol = 3, byrow=T)
SpViol25P <- matrix(c(sim1041, sim1042, sim1043, sim1044, sim1045), ncol = 3, byrow=T)

#Make a combined matrix
SpViolP <- cbind(SpViol25P, SpViol50P, SpViol100P, SpViol200P)
plotfunction(sim61, sim62, sim63, sim64, sim65)

SpViolPR <- RoundMatrix(SpViolP, 2) #(Table 10 in thesis)

#Homogeneity, power
#0.5Sigma
set.seed(1411)
sim1411 <- simulation(p, list(n25, n200), list(mu1i, mu2i), list((1/2)*Sigma1, Sigma1), 1000)
set.seed(1412)
sim1412 <- simulation(p, list(n50, n200), list(mu1i, mu2i), list((1/2)*Sigma1, Sigma1), 1000)
set.seed(1413)
sim1413 <- simulation(p, list(n75, n200), list(mu1i, mu2i), list((1/2)*Sigma1, Sigma1), 1000)
set.seed(1414)
sim1414 <- simulation(p, list(n100, n200), list(mu1i, mu2i), list((1/2)*Sigma1, Sigma1), 1000)
set.seed(1415)
sim1415 <- simulation(p, list(n125, n200), list(mu1i, mu2i), list((1/2)*Sigma1, Sigma1), 1000)
set.seed(1416)
sim1416 <- simulation(p, list(n150, n200), list(mu1i, mu2i), list((1/2)*Sigma1, Sigma1), 1000)
set.seed(1417)
sim1417 <- simulation(p, list(n175, n200), list(mu1i, mu2i), list((1/2)*Sigma1, Sigma1), 1000)
set.seed(1418)
sim1418 <- simulation(p, list(n200, n200), list(mu1i, mu2i), list((1/2)*Sigma1, Sigma1), 1000)

plotPowerSampleSize200(sim1411, sim1412, sim1413, sim1414, sim1415, sim1416, sim1417, sim1418, c(25, 50, 75, 100, 125, 150, 175, 200))
Homogen0.5P <- rbind(sim1411, sim1412, sim1413, sim1414, sim1415, sim1416, sim1417, sim1418)

#0.75Sigma
set.seed(1511)
sim1511 <- simulation(p, list(n25, n200), list(mu1i, mu2i), list((3/4)*Sigma1, Sigma1), 1000)
set.seed(1512)
sim1512 <- simulation(p, list(n50, n200), list(mu1i, mu2i), list((3/4)*Sigma1, Sigma1), 1000)
set.seed(1513)
sim1513 <- simulation(p, list(n75, n200), list(mu1i, mu2i), list((3/4)*Sigma1, Sigma1), 1000)
set.seed(1514)
sim1514 <- simulation(p, list(n100, n200), list(mu1i, mu2i), list((3/4)*Sigma1, Sigma1), 1000)
set.seed(1515)
sim1515 <- simulation(p, list(n125, n200), list(mu1i, mu2i), list((3/4)*Sigma1, Sigma1), 1000)
set.seed(1516)
sim1516 <- simulation(p, list(n150, n200), list(mu1i, mu2i), list((3/4)*Sigma1, Sigma1), 1000)
set.seed(1517)
sim1517 <- simulation(p, list(n175, n200), list(mu1i, mu2i), list((3/4)*Sigma1, Sigma1), 1000)
set.seed(1518)
sim1518 <- simulation(p, list(n200, n200), list(mu1i, mu2i), list((3/4)*Sigma1, Sigma1), 1000)

plotPowerSampleSize200(sim1511, sim1512, sim1513, sim1514, sim1515, sim1516, sim1517, sim1518, c(25, 50, 75, 100, 125, 150, 175, 200))
Homogen0.75P <- rbind(sim1511, sim1512, sim1513, sim1514, sim1515, sim1516, sim1517, sim1518)

#1.25Sigma
set.seed(1521)
sim1521 <- simulation(p, list(n25, n200), list(mu1i, mu2i), list((1.25)*Sigma1, Sigma1), 1000)
set.seed(1522)
sim1522 <- simulation(p, list(n50, n200), list(mu1i, mu2i), list((1.25)*Sigma1, Sigma1), 1000)
set.seed(1523)
sim1523 <- simulation(p, list(n75, n200), list(mu1i, mu2i), list((1.25)*Sigma1, Sigma1), 1000)
set.seed(1524)
sim1524 <- simulation(p, list(n100, n200), list(mu1i, mu2i), list((1.25)*Sigma1, Sigma1), 1000)
set.seed(1525)
sim1525 <- simulation(p, list(n125, n200), list(mu1i, mu2i), list((1.25)*Sigma1, Sigma1), 1000)
set.seed(1526)
sim1526 <- simulation(p, list(n150, n200), list(mu1i, mu2i), list((1.25)*Sigma1, Sigma1), 1000)
set.seed(1527)
sim1527 <- simulation(p, list(n175, n200), list(mu1i, mu2i), list((1.25)*Sigma1, Sigma1), 1000)
set.seed(1528)
sim1528 <- simulation(p, list(n200, n200), list(mu1i, mu2i), list((1.25)*Sigma1, Sigma1), 1000)

plotPowerSampleSize200(sim1521, sim1522, sim1523, sim1524, sim1525, sim1526, sim1527, sim1528, c(25, 50, 75, 100, 125, 150, 175, 200))
Homogen1.25P <- rbind(sim1521, sim1522, sim1523, sim1524, sim1525, sim1526, sim1527, sim1528)

#1.5Sigma
set.seed(1611)
sim1611 <- simulation(p, list(n25, n200), list(mu1i, mu2i), list((1.5)*Sigma1, Sigma1), 1000)
set.seed(1612)
sim1612 <- simulation(p, list(n50, n200), list(mu1i, mu2i), list((1.5)*Sigma1, Sigma1), 1000)
set.seed(1613)
sim1613 <- simulation(p, list(n75, n200), list(mu1i, mu2i), list((1.5)*Sigma1, Sigma1), 1000)
set.seed(1614)
sim1614 <- simulation(p, list(n100, n200), list(mu1i, mu2i), list((1.5)*Sigma1, Sigma1), 1000)
set.seed(1615)
sim1615 <- simulation(p, list(n125, n200), list(mu1i, mu2i), list((1.5)*Sigma1, Sigma1), 1000)
set.seed(1616)
sim1616 <- simulation(p, list(n150, n200), list(mu1i, mu2i), list((1.5)*Sigma1, Sigma1), 1000)
set.seed(1617)
sim1617 <- simulation(p, list(n175, n200), list(mu1i, mu2i), list((1.5)*Sigma1, Sigma1), 1000)
set.seed(1618)
sim1618 <- simulation(p, list(n200, n200), list(mu1i, mu2i), list((1.5)*Sigma1, Sigma1), 1000)

plotPowerSampleSize200(sim1611, sim1612, sim1613, sim1614, sim1615, sim1616, sim1617, sim1618, c(25, 50, 75, 100, 125, 150, 175, 200))
Homogen1.5P <- rbind(sim1611, sim1612, sim1613, sim1614, sim1615, sim1616, sim1617, sim1618)

#1.75Sigma
set.seed(1621)
sim1621 <- simulation(p, list(n25, n200), list(mu1i, mu2i), list((1.75)*Sigma1, Sigma1), 1000)
set.seed(1622)
sim1622 <- simulation(p, list(n50, n200), list(mu1i, mu2i), list((1.75)*Sigma1, Sigma1), 1000)
set.seed(1623)
sim1623 <- simulation(p, list(n75, n200), list(mu1i, mu2i), list((1.75)*Sigma1, Sigma1), 1000)
set.seed(1624)
sim1624 <- simulation(p, list(n100, n200), list(mu1i, mu2i), list((1.75)*Sigma1, Sigma1), 1000)
set.seed(1625)
sim1625 <- simulation(p, list(n125, n200), list(mu1i, mu2i), list((1.75)*Sigma1, Sigma1), 1000)
set.seed(1626)
sim1626 <- simulation(p, list(n150, n200), list(mu1i, mu2i), list((1.75)*Sigma1, Sigma1), 1000)
set.seed(1627)
sim1627 <- simulation(p, list(n175, n200), list(mu1i, mu2i), list((1.75)*Sigma1, Sigma1), 1000)
set.seed(1628)
sim1628 <- simulation(p, list(n200, n200), list(mu1i, mu2i), list((1.75)*Sigma1, Sigma1), 1000)

plotPowerSampleSize200(sim1621, sim1622, sim1623, sim1624, sim1625, sim1626, sim1627, sim1628, c(25, 50, 75, 100, 125, 150, 175, 200))
Homogen1.75P <- rbind(sim1621, sim1622, sim1623, sim1624, sim1625, sim1626, sim1627, sim1628)

#2Sigma
set.seed(1711)
sim1711 <- simulation(p, list(n25, n200), list(mu1i, mu2i), list((2)*Sigma1, Sigma1), 1000)
set.seed(1712)
sim1712 <- simulation(p, list(n50, n200), list(mu1i, mu2i), list((2)*Sigma1, Sigma1), 1000)
set.seed(1713)
sim1713 <- simulation(p, list(n75, n200), list(mu1i, mu2i), list((2)*Sigma1, Sigma1), 1000)
set.seed(1714)
sim1714 <- simulation(p, list(n100, n200), list(mu1i, mu2i), list((2)*Sigma1, Sigma1), 1000)
set.seed(1715)
sim1715 <- simulation(p, list(n125, n200), list(mu1i, mu2i), list((2)*Sigma1, Sigma1), 1000)
set.seed(1716)
sim1716 <- simulation(p, list(n150, n200), list(mu1i, mu2i), list((2)*Sigma1, Sigma1), 1000)
set.seed(1717)
sim1717 <- simulation(p, list(n175, n200), list(mu1i, mu2i), list((2)*Sigma1, Sigma1), 1000)
set.seed(1718)
sim1718 <- simulation(p, list(n200, n200), list(mu1i, mu2i), list((2)*Sigma1, Sigma1), 1000)

plotPowerSampleSize200(sim1711, sim1712, sim1713, sim1714, sim1715, sim1716, sim1717, sim1718, c(25, 50, 75, 100, 125, 150, 175, 200))
Homogen2P <- rbind(sim1711, sim1712, sim1713, sim1714, sim1715, sim1716, sim1717, sim1718)

#Combine the results in matrices
HomogenP <- cbind(Homogen0.5P, Homogen0.75P, Homogen1.25P, Homogen1.5P, Homogen1.75P, Homogen2P)
HomogenPA <- cbind(Homogen0.5P, Homogen0.75P)
HomogenPB <- cbind(Homogen1.25P, Homogen1.5P, Homogen1.75P, Homogen2P)

#Round the matrices
HomogenPR <- RoundMatrix(HomogenP)
HomogenPAR <- RoundMatrix(HomogenPA) #(Table 11 in thesis)
HomogenPBR <- RoundMatrix(HomogenPB) #(Table 12 in thesis)


#Both violated, power
#0.75 Sigma, 0.8 epsilon
set.seed(1811)
sim1811 <- simulation(p, list(n25, n200), list(mu1i, mu2i), list((3/4)*Sigma3, Sigma3), 1000)
set.seed(1812)
sim1812 <- simulation(p, list(n50, n200), list(mu1i, mu2i), list((3/4)*Sigma3, Sigma3), 1000)
set.seed(1813)
sim1813 <- simulation(p, list(n75, n200), list(mu1i, mu2i), list((3/4)*Sigma3, Sigma3), 1000)
set.seed(1814)
sim1814 <- simulation(p, list(n100, n200), list(mu1i, mu2i), list((3/4)*Sigma3, Sigma3), 1000)
set.seed(1815)
sim1815 <- simulation(p, list(n125, n200), list(mu1i, mu2i), list((3/4)*Sigma3, Sigma3), 1000)
set.seed(1816)
sim1816 <- simulation(p, list(n150, n200), list(mu1i, mu2i), list((3/4)*Sigma3, Sigma3), 1000)
set.seed(1817)
sim1817 <- simulation(p, list(n175, n200), list(mu1i, mu2i), list((3/4)*Sigma3, Sigma3), 1000)
set.seed(1818)
sim1818 <- simulation(p, list(n200, n200), list(mu1i, mu2i), list((3/4)*Sigma3, Sigma3), 1000)

plotPowerSampleSize200(sim1811, sim1812, sim1813, sim1814, sim1815, sim1816, sim1817, sim1818, c(25, 50, 75, 100, 125, 150, 175, 200))
Homogen0.75Sph0.8P <- rbind(sim1811, sim1812, sim1813, sim1814, sim1815, sim1816, sim1817, sim1818)

#1.25 Sigma, 0.8 epsilon
set.seed(1821)
sim1821 <- simulation(p, list(n25, n200), list(mu1i, mu2i), list((1.25)*Sigma3, Sigma3), 1000)
set.seed(1822)
sim1822 <- simulation(p, list(n50, n200), list(mu1i, mu2i), list((1.25)*Sigma3, Sigma3), 1000)
set.seed(1823)
sim1823 <- simulation(p, list(n75, n200), list(mu1i, mu2i), list((1.25)*Sigma3, Sigma3), 1000)
set.seed(1824)
sim1824 <- simulation(p, list(n100, n200), list(mu1i, mu2i), list((1.25)*Sigma3, Sigma3), 1000)
set.seed(1825)
sim1825 <- simulation(p, list(n125, n200), list(mu1i, mu2i), list((1.25)*Sigma3, Sigma3), 1000)
set.seed(1826)
sim1826 <- simulation(p, list(n150, n200), list(mu1i, mu2i), list((1.25)*Sigma3, Sigma3), 1000)
set.seed(1827)
sim1827 <- simulation(p, list(n175, n200), list(mu1i, mu2i), list((1.25)*Sigma3, Sigma3), 1000)
set.seed(1828)
sim1828 <- simulation(p, list(n200, n200), list(mu1i, mu2i), list((1.25)*Sigma3, Sigma3), 1000)

plotPowerSampleSize200(sim1821, sim1822, sim1823, sim1824, sim1825, sim1826, sim1827, sim1828, c(25, 50, 75, 100, 125, 150, 175, 200))
Homogen1.25Sph0.8P <- rbind(sim1821, sim1822, sim1823, sim1824, sim1825, sim1826, sim1827, sim1828)

#1.5 Sigma, 0.8 epsilon
set.seed(1851)
sim1851 <- simulation(p, list(n25, n200), list(mu1i, mu2i), list((1.5)*Sigma3, Sigma3), 1000)
set.seed(1852)
sim1852 <- simulation(p, list(n50, n200), list(mu1i, mu2i), list((1.5)*Sigma3, Sigma3), 1000)
set.seed(1853)
sim1853 <- simulation(p, list(n75, n200), list(mu1i, mu2i), list((1.5)*Sigma3, Sigma3), 1000)
set.seed(1854)
sim1854 <- simulation(p, list(n100, n200), list(mu1i, mu2i), list((1.5)*Sigma3, Sigma3), 1000)
set.seed(1855)
sim1855 <- simulation(p, list(n125, n200), list(mu1i, mu2i), list((1.5)*Sigma3, Sigma3), 1000)
set.seed(1856)
sim1856 <- simulation(p, list(n150, n200), list(mu1i, mu2i), list((1.5)*Sigma3, Sigma3), 1000)
set.seed(1857)
sim1857 <- simulation(p, list(n175, n200), list(mu1i, mu2i), list((1.5)*Sigma3, Sigma3), 1000)
set.seed(1858)
sim1858 <- simulation(p, list(n200, n200), list(mu1i, mu2i), list((1.5)*Sigma3, Sigma3), 1000)

plotPowerSampleSize200(sim1851, sim1852, sim1853, sim1854, sim1855, sim1856, sim1857, sim1858, c(25, 50, 75, 100, 125, 150, 175, 200))
Homogen1.5Sph0.8P <- rbind(sim1851, sim1852, sim1853, sim1854, sim1855, sim1856, sim1857, sim1858)

#0.75 Sigma, 0.7 epsilon
set.seed(1831)
sim1831 <- simulation(p, list(n25, n200), list(mu1i, mu2i), list((3/4)*Sigma4, Sigma4), 1000)
set.seed(1832)
sim1832 <- simulation(p, list(n50, n200), list(mu1i, mu2i), list((3/4)*Sigma4, Sigma4), 1000)
set.seed(1833)
sim1833 <- simulation(p, list(n75, n200), list(mu1i, mu2i), list((3/4)*Sigma4, Sigma4), 1000)
set.seed(1834)
sim1834 <- simulation(p, list(n100, n200), list(mu1i, mu2i), list((3/4)*Sigma4, Sigma4), 1000)
set.seed(1835)
sim1835 <- simulation(p, list(n125, n200), list(mu1i, mu2i), list((3/4)*Sigma4, Sigma4), 1000)
set.seed(1836)
sim1836 <- simulation(p, list(n150, n200), list(mu1i, mu2i), list((3/4)*Sigma4, Sigma4), 1000)
set.seed(1837)
sim1837 <- simulation(p, list(n175, n200), list(mu1i, mu2i), list((3/4)*Sigma4, Sigma4), 1000)
set.seed(1838)
sim1838 <- simulation(p, list(n200, n200), list(mu1i, mu2i), list((3/4)*Sigma4, Sigma4), 1000)

plotPowerSampleSize200(sim1831, sim1832, sim1833, sim1834, sim1835, sim1836, sim1837, sim1838, c(25, 50, 75, 100, 125, 150, 175, 200))
Homogen0.75Sph0.7P <- rbind(sim1831, sim1832, sim1833, sim1834, sim1835, sim1836, sim1837, sim1838)

#1.25 Sigma, 0.7 epsilon
set.seed(1841)
sim1841 <- simulation(p, list(n25, n200), list(mu1i, mu2i), list((1.25)*Sigma4, Sigma4), 1000)
set.seed(1842)
sim1842 <- simulation(p, list(n50, n200), list(mu1i, mu2i), list((1.25)*Sigma4, Sigma4), 1000)
set.seed(1843)
sim1843 <- simulation(p, list(n75, n200), list(mu1i, mu2i), list((1.25)*Sigma4, Sigma4), 1000)
set.seed(1844)
sim1844 <- simulation(p, list(n100, n200), list(mu1i, mu2i), list((1.25)*Sigma4, Sigma4), 1000)
set.seed(1845)
sim1845 <- simulation(p, list(n125, n200), list(mu1i, mu2i), list((1.25)*Sigma4, Sigma4), 1000)
set.seed(1846)
sim1846 <- simulation(p, list(n150, n200), list(mu1i, mu2i), list((1.25)*Sigma4, Sigma4), 1000)
set.seed(1847)
sim1847 <- simulation(p, list(n175, n200), list(mu1i, mu2i), list((1.25)*Sigma4, Sigma4), 1000)
set.seed(1848)
sim1848 <- simulation(p, list(n200, n200), list(mu1i, mu2i), list((1.25)*Sigma4, Sigma4), 1000)

plotPowerSampleSize200(sim1841, sim1842, sim1843, sim1844, sim1845, sim1846, sim1847, sim1848, c(25, 50, 75, 100, 125, 150, 175, 200))
Homogen1.25Sph0.7P <- rbind(sim1841, sim1842, sim1843, sim1844, sim1845, sim1846, sim1847, sim1848)

#1.50 Sigma, 0.8 epsilon
set.seed(1861)
sim1861 <- simulation(p, list(n25, n200), list(mu1i, mu2i), list((1.5)*Sigma4, Sigma4), 1000)
set.seed(1862)
sim1862 <- simulation(p, list(n50, n200), list(mu1i, mu2i), list((1.5)*Sigma4, Sigma4), 1000)
set.seed(1863)
sim1863 <- simulation(p, list(n75, n200), list(mu1i, mu2i), list((1.5)*Sigma4, Sigma4), 1000)
set.seed(1864)
sim1864 <- simulation(p, list(n100, n200), list(mu1i, mu2i), list((1.5)*Sigma4, Sigma4), 1000)
set.seed(1865)
sim1865 <- simulation(p, list(n125, n200), list(mu1i, mu2i), list((1.5)*Sigma4, Sigma4), 1000)
set.seed(1866)
sim1866 <- simulation(p, list(n150, n200), list(mu1i, mu2i), list((1.5)*Sigma4, Sigma4), 1000)
set.seed(1867)
sim1867 <- simulation(p, list(n175, n200), list(mu1i, mu2i), list((1.5)*Sigma4, Sigma4), 1000)
set.seed(1868)
sim1868 <- simulation(p, list(n200, n200), list(mu1i, mu2i), list((1.5)*Sigma4, Sigma4), 1000)

plotPowerSampleSize200(sim1861, sim1862, sim1863, sim1864, sim1865, sim1866, sim1867, sim1868, c(25, 50, 75, 100, 125, 150, 175, 200))
Homogen1.5Sph0.7P <- rbind(sim1861, sim1862, sim1863, sim1864, sim1865, sim1866, sim1867, sim1868)

#Combine all three times two results

HomoSphVP <- cbind(Homogen0.75Sph0.7P, Homogen1.25Sph0.7P, Homogen1.5Sph0.7P,
                  Homogen0.75Sph0.8P, Homogen1.25Sph0.8P, Homogen1.5Sph0.8P)

HomoSphVPR <- RoundMatrix(HomoSphVP) #(Table 13 in thesis)


#APPLICATION EXISTING DATA
#---------------------------------------------------------------
#Activate required packages
library(haven)
library(biotools)

#Retrieve data
dataset <- read_sav("Master thesis/DatVil24erv123JKNewIDMergedCleaned.sav")
View(dataset)

#Create the dataset with the required variables
X2 <- dataset[,c(13, 6, 7, 8, 9, 10, 11)]

#Value epsilon
epsilonf(cov(X2[,2:4]))
epsilonf(cov(X2[,5:7]))

#Calculate Box's M
boxM(X2[,2:4], X2[,1])

#Define the n's for every group
nData <- list(30, 35)

#Replace the zeros and ones with ones and twos
X2[,1][X2[,1]==1] <- 2
X2[,1][X2[,1]==0] <- 1

#Adapt the GLR function so it also returns Chi-square value
GLRdata <- function(X, p, n){
  D <- as.data.frame(X)
  names(D)<-c("g","X1","X2", "X3")
  
  # het gesatureerde model
  #eerst groep, dan tijdstip
  model1<-"X1~c(m11,m21)*1
  X2~c(m12,m22)*1
  X3~c(m13,m23)*1
  X1~~NA*X2
  X1~~NA*X3
  X2~~NA*X3
  "
  fit1<-sem(model1,data=D,group="g")
  
  # het gerestricteerde model (geen interactie)
  model2<-"X1~c(m11,m21)*1
  X2~c(m12,m22)*1
  X3~c(m13,m23)*1
  X1~~NA*X2
  X1~~NA*X3
  X2~~NA*X3
  m22 == m12 + m21 - m11
  m23 == m13 + m21 - m11
  "
  
  fit2<-sem(model2,data=D,group="g")
  
  X2 <- -2*(logLik(fit2)[[1]] - logLik(fit1)[[1]])
  list(X2, 1-pchisq(X2,2))
}

GLRdata(X2, p, nData) #Also gives Chi-square

#Descriptives
#Means and standard deviations of agency #(Table 14 in thesis)
X2means1 <- apply(X2[which(X2[,1]==1),2:4], 2, mean)
X2means2 <- apply(X2[which(X2[,1]==2),2:4], 2, mean)

X2sd1 <- apply(X2[which(X2[,1]==1),2:4], 2, sd)
X2sd2 <- apply(X2[which(X2[,1]==2),2:4], 2, sd)

#Check the covariances
cov(X2[which(X2[,1]==1),2:4])
cov(X2[which(X2[,1]==2),2:4])

cor(X2[which(X2[,1]==1),2:4])
cor(X2[which(X2[,1]==2),2:4])

#Overall covariance
cor(X2[,2:4])

#Means and standard deviations of communion #(Table 14 in thesis)
X2means1C <- apply(X2[which(X2[,1]==1),5:7], 2, mean)
X2means2C <- apply(X2[which(X2[,1]==2),5:7], 2, mean)

X2sd1C <- apply(X2[which(X2[,1]==1),5:7], 2, sd)
X2sd2C <- apply(X2[which(X2[,1]==2),5:7], 2, sd)

#Plot agency #(Figure 5 in thesis)
plot(c(1, 2, 3), X2means1, type = "l", col="grey30", ylim=c(0.00,0.20), xlab = "year", ylab = "agency", lty=1, lwd=1.5, axes=FALSE)
axis(side=1, at=c(0:3))
axis(side=2, at=seq(0.00, 0.20, by=0.05))
lines(c(1, 2, 3), X2means2, col = "black", lty = 2, lwd = 1.5)
legend(2.2, 0.05, legend=c("Low on communion", "High on communion"),
       col=c("grey30", "black"), lty=c(1, 2), lwd = 1.5, cex=0.8)

#Plot communion
plot(c(1, 2, 3), X2means1C, type = "l", col="grey30", ylim=c(0.00,0.40), xlab = "year", ylab = "communion", lty=1, lwd=1.5, axes=FALSE)
axis(side=1, at=c(0:3))
axis(side=2, at=seq(0.00, 0.40, by=0.05))
lines(c(1, 2, 3), X2means2C, col = "black", lty = 2, lwd = 1.5)
legend(2.2, 0.10, legend=c("Low on communion", "High on communion"),
       col=c("grey30", "black"), lty=c(1, 2), lwd = 1.5, cex=0.8)

#Remove the data for the application from the environment
rm(dataset)
rm(X2)

#Save the environment
save.image(file = "Results.RData")

#-------------------------------------------------------------------------------------------------------
#-------------------------------------------------------------------------------------------------------
#-------------------------------------------------------------------------------------------------------