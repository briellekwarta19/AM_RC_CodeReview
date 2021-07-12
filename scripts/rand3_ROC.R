#libraries
library(nimble)
library(here)
library(gtools)
library(MCMCvis)
library(abind)
library(gtools)
library(tidyverse)
library(Rlab)


#Problems with this code:
#Sometimes I need to restart R in order to get the code to run correctly. 
#I want to have occupancy data at 3 random sites and NAs else where…
  #BUT! But I fail to integrate data in NIMBLE when I have missing data… 
  #So right now in the script I’m just collecting occupancy and count data everywhere. I need help!! 
  
#### Script Description ####
#The script can be divided into a few steps: 

## A. Model Set up ##
# STEP 1: Initial simulation of true abundance= where we simulate true system dynamics (R)
#We simulate population abundance and data collection (data collection feeds into step 2)

# STEP 2: First learning process= where we generate our belief of the system (NIMBLE)

#TO DO: I need to write the code for steps 3-6. but this is similar to the abund3_R script#
# STEP 3: First decision process= where we select our management and data collection rule
#This is based on our belief of the system
#Management alternative= place 2 traps (n.trap) at 3 random segments 
#Data collection rule= collect removal data (with the management alternative specified above)
  #and collect occupancy and count data (we want to collect this at 3 random segments)
  #But... my code isn't working when I have NA data at the locations where I do not collect information. 

## B. General Procedure: ##
# STEP 4: Simulate Population = simulating true system dynamics
# STEP 5: Learning process = where we generate our belief of the system using data we collected
# STEP 6: Update decision process = based on our belief, we update management and data collection rules

#For me NIMBLE only runs if I do this... if you downloaded nimble correctly you probably shouldnt do this!
path <- Sys.getenv('PATH')
newPath <- paste("C:\\Rtools\\bin;C:\\Rtools\\mingw_64\\bin;",
                 path, sep = "")
Sys.setenv(PATH = newPath)

#### Nimble Model ####
#index i (or h) = segment
#index j = primary sampling period
#index k = secondary sampling period

modcode <- nimbleCode({
  ##----------Likelihood----------##
  ## Data Collection Model: ##
  for (i in 1:I) {
    #Initial abundance: 
    N.belief[i, 1, 1] ~ dbin(psi,10000) #tries to get initial conditions near 5000 ish
    for (j in 2:J) {
      #Abundance for 1st removal of the primary period (aka. 1st secondary removal)
      N.belief[i, j, 1] <- D.after[i,j-1] 
    }
    for (j in 1:J) {
      for (k in 2:K) {
        # Abundance after 1st secondary removal
        N.index[i,j,k] <- step(N.belief[i,j,k-1] - Y[i,j,k-1])#if negative population N.index = 0, else 1
        N.belief[i,j,k] <- N.index[i,j,k]*(N.belief[i,j,k-1] - Y[i,j,k-1]) 
      }
    }
  }
  
  for(j in 1:J){ #index j = primary sampling period
    for(i in 1:I){ #index i = segment
      for(k in 1: K){
        # Removal data:
        #Calculating max individuals that could be trapped: 
        #traps can capture between 0 to 50 individuals
        #But if population < 50, our max capture per trap is N.belief/n.trap
        max.trap[i,j,k] <- min(50,(N.belief[i,j,k]/n.trap)) #n.trap = number of traps per segment where trapping occurs
        
        #Y.trap = number of individuals trapped and removed per trap
        Y.trap[i,j,k] ~ dunif(0, max.trap[i,j,k])
        
        #Y = number of individuals trapped and removed (whole number)
        Y[i,j,k] <- ceiling(Y.trap[i,j,k]*n.trap*trap.index[i,j]) 
        
        #Detection non-detection data: 
        O[i,j,k] ~ dbern(o.site[i,j]) #detection/non-detection
        
        #Count data: 
        C[i,j,k] ~ dbin(c.site, N.belief[i,j,k]) #count data
      }
      # Abundance at the end of primary removal period j:
      R.index[i,j] <- step(N.belief[i,j,K] - Y[i,j,K]) #if negative population R.index = 0, else 1
      R[i,j] <- R.index[i,j]*(N.belief[i,j,K] - Y[i,j,K])
      
      ##--Population growth--##
      #Surviving individuals
      S[i,j] ~ dbin(phi, R[i,j]) 
      
      #Recruited individuals
      G[i,j] ~ dpois(lambda[j]) 
      D[i,j] <- S[i,j] + R[i,j] #population available for dispersal
      
      #--Movement process--#
      # Movemement probabilities:
      pi[i,1:I,j] ~ ddirch(alpha[1:I]) #movement probabilities
      #Calculating M[i,h,j] = number of individuals moving from i to h at j
      
      # num of individs that move to seg 1
      M[i,1,j] ~ dbin(pi[i,1,j],D[i,j])
      for(h in 2:(I-1)){
        # Calculating num of individs that move to all other segments (not final seg)
        # Probability of going to segment h given not in any previous segment
        pi.prob[i,h,j] <- pi[i,h,j]/(1-sum(pi[i,1:(h-1),j])) # movement prob.
        ind.move[i,h,j] <- (D[i,j] - sum(M[i,1:(h-1),j])) # num of individs moving
        M[i,h,j] ~ dbin(pi.prob[i,h,j],ind.move[i,h,j])
      }
      # Calculating the num of individs that move to the final segment  
      ind.move.last[i,j] <- sum(M[i,1:(I-1),j])
      M[i,I,j] <- (D[i,j]  - ind.move.last[i,j])
      
      D.after[i,j] <- sum(M[1:I,i,j]) # number of individuals after movement (feeds into removal model)
      
    } #Ends segment loop I
  } #Ends primary loop J
  
  
  ##----------Priors----------##
  #Parameter for initial population
  psi ~ dunif(0.49,0.51) 
  
  #survival rate (constant through time)
  phi ~ dunif(0.7, 0.95)
  
  #fecundity
  for(j in 1:J){
    lambda[j] ~dunif(0,10)
  }
})

#############################################################################
#### A. Model set up ####
#### STEP 1: Initial simulation of true system dynamics ####
#Simulate 2 years of data
J <- 2 # Number of primary periods
K <- 2 # Number of secondary periods (abundance measurements per site)
I <- 5 # Number of total spatial units
n.trap <- 2 #Number of traps

## Initial arrays ## 
Y.trap <- array(0, dim = c(I,J,K)) #Y.trap = individuals removed per trap
Y <- array(0, dim = c(I,J,K)) #Y = removal data (integer)

o.base <- runif(1,0.5, 1) # baseline detection from occupancy data
o.site <- array(0, dim = c(I,J)) #detection from occupancy data at each site
O <- array(0, dim = c(I,J,K)) #O = detection/nondetection data

c.site <- runif(1, 0.2, 0.5) #detection rate from count data
C <- array(0, dim = c(I,J,K)) #C = count data

N.truth <- array(0, dim = c(I,J,K)) #N.truth = true population abundance
N.truth[,1,1] <- 5000 #Let initial abundance be 1000
R <- array(0, dim = c(I,J)) #R = population remaining after removals at j

#For the first primary select 3 random segments because all initial abundances are the same
max.trap <- array(0,dim = c(I,J,K))
site.trap1 <- sample(1:I, 3) #Select 3 random trapping sites for j = 1 
site.trap2 <- sample(1:I, 3) #Select 3 random trapping sites for j = 2 
trap.index<- array(0, dim = c(I, J)) #Trapping index: either trap or not

for(i in 1:I){
  if(i %in% site.trap1){trap.index[i,1] <- 1} else{
    trap.index[i,1] <- 0
  }
  if(i %in% site.trap2){trap.index[i,2] <- 1} else{
    trap.index[i,2] <- 0
  }
}

site.occ1 <- sample(1:I, 3) #Select 3 random trapping sites for j = 1 
site.occ2 <- sample(1:I, 3) #Select 3 random trapping sites for j = 2 

site.count1 <- sample(1:I, 3) #Select 3 random trapping sites for j = 1 
site.count2 <- sample(1:I, 3) #Select 3 random trapping sites for j = 2 

phi <- runif(1,0.7, 0.95) #Survival rate for each primary period
S <- array(0, dim = c(I,J)) #Individuals that survive

lambda <- array(0, dim = c(J)) #Fecundity
G <- array(0, dim = c(I,J)) #G = Population recruited
D <- array(0, dim = c(I,J)) #D = individuals available for dispersal 

pi <- array(0, dim = c(I,I,J)) #pi = movement probability to each segment at j
alpha <- rep(1,I) #goes into dirichlet prior
M <- array(0, dim = c(I,I,J)) #M = number of individuals that move 

pi.prob <- array(0, dim = c(I,I,J)) # Probability of going to segment h given not in any previous segment
ind.move<- array(0, dim = c(I,I,J)) # Number of individuals moving
ind.move.last <- array(0, dim = c(I,J)) # Number of individuals that move to the last segment

D.after <- array(0, dim = c(I,J)) #D.after = individuals that move to i 

#filling in the initial simulation data
reps <- 50
for(rep in 1:reps) {
  for(j in 1:J){
    for(i in 1:I){
      ## Data Collection Model: ##
      for(k in 1: K){
        #Removal data
        max.trap[i,j,k] <- min(50,(N.truth[i,j,k]/n.trap))
        Y.trap[i,j,k] <- runif(1,0, max.trap[i,j,k])
        Y[i,j,k] <- ceiling(Y.trap[i,j,k]*n.trap*trap.index[i,j])
      
        #Detection/non-detection data
        o.site[i,j] <- (1-((1-o.base)^N.truth[i,j,1])) #calculating site detection
        
        #I want it to be something like this (the uncommented script) but it doesn't work!
        #if(i %in% site.occ1){
          O[i,1,k] <- rbern(1,o.site[i,1])
        #} else{
         # O[i,1,k] <- NA
        #}
        
        #if(i %in% site.occ2){
          O[i,2,k] <- rbern(1,o.site[i,2])
        #} else{
         # O[i,2,k] <- NA
        #}
      
        #Count data
        # if (i %in% site.count1){
          C[i,1,k] <- rbinom(1, N.truth[i,1,k], c.site)
        # } else {
        #   C[i,1,k] <- NA
        # }
        
        #if (i %in% site.count2){
          C[i,2,k] <- rbinom(1, N.truth[i,2,k], c.site)
        #} else {
         # C[i,2,k] <- NA
        #}
      
      }
      
      #Abundance after 1st secondary removal
      for(k in 2:K){ 
        if((N.truth[i,j,k-1] - Y[i,j,k-1]) < 0){
          N.truth[i,j,k] <- 0
        } else{
          N.truth[i,j,k] <- (N.truth[i,j,k-1] - Y[i,j,k-1])
        }
      }
      
      #Abundance at the end of primary removal period j:
      if((N.truth[i,j,K] - Y[i,j,K]) < 0){
        R[i,j] <- 0
      } else{
        R[i,j] <- (N.truth[i,j,K] - Y[i,j,K]) 
      }
      
      #Abundance for 1st removal of the primary period (aka. 1st secondary removal)
      #This abundance arises from movement into the segment at the end of the previous primary period
      N.truth[i,2:J, 1] <- D.after[i,1:(J-1)]
      
      ## Population Change Model: ##  
      lambda[j] <- runif(1,0,10) #Population gains
      S[i,j] <- rbinom(1,R[i,j], phi) #surviving individuals
      G[i,j] <- rpois(1,lambda[j]) #population available for dispersal 
      
      D[i,j] <- G[i,j] + S[i,j]
      
      ## Movement Model: ##
      #Movement probabilities:
      alpha <- rep(1:I)
      
      #Calculating movement probabilities: 
      pi[i,1:I,j] <- rdirichlet(1, alpha[1:I]) #good congugate prior for multinomial
      
      # Calculating the number of individuals that move to segment 1
      M[i,1,j] <- rbinom(1, D[i,j], pi[i,1,j])
      
      # Calculating the number of individuals that move to all other segments (not final segment)
      # prob of going to segment h given not in any previous segment
      for(h in 2:(I-1)){
        pi.prob[i,h,j] <- pi[i,h,j]/(1-sum(pi[i,1:(h-1),j])) #calculating movement probability
        ind.move[i,h,j] <- (D[i,j] - sum(M[i,(1:(h-1)),j])) #counting individuals that are moving
        M[i,h,j] <- rbinom(1,ind.move[i,h,j], pi.prob[i,h,j]) 
      }
      
      # Calculating the number of individuals that move to the final segment
      ind.move.last[i,j] <- sum(M[i,(1:(I-1)),j])
      M[i,I,j] <- (D[i,j]  - ind.move.last[i,j])
      
      D.after[i,j] <- sum(M[,i,j]) #number of individuals after movement (feeds into removal model)
      
    } #Ends segment loop I
  } #Ends primary loop J
}

#### STEP 2: First learning process ####
#Here we run our NIMBLE model to estimate our belief of the population
#Constants
my.constants <- list(I = I, #Number of segments
                     J = J, #Number of primary periods
                     K = K, #Number of secondary periods
                     n.trap = n.trap, #Number of traps
                     trap.index = trap.index, #which segments have traps
                     o.site= o.site, #detection probability from occupancy data for each site
                     c.site = c.site, #detection probability from count data for each site
                     alpha = rep(1, I)) #Data for dirichlet process
#Data
my.data <- list(Y.trap = Y.trap, O = O, C = C) #removal data

# Initial values
N.start <- array(c(rep(5000,I),rep(5000, I + 2*(J-1)*(I))), dim = c(I,J,K))
initial.values <- function(){list(N.belief = N.start)}

# Parameters monitored
parameters.to.save <- c("N.belief", "lambda", "pi","D.after") #should be D.after

# MCMC settings
n.iter <- 5000
n.burnin <- 1000
n.chains <- 3

mcmc <- nimbleMCMC(code = modcode, 
                   constants = my.constants,
                   data = my.data,              
                   inits = initial.values,
                   monitors = parameters.to.save,
                   niter = n.iter,
                   nburnin = n.burnin, 
                   nchains = n.chains, 
                   summary = TRUE)

############################################################################

# Inspect results
#Check to make sure N[i,j,k] = D.after[i,j-1]
results <- mcmc$summary$all.chains
results <- as.data.frame(results)
results <- cbind(rownames(results), results)

colnames(results)[1] <- "param"


results %>% 
  filter(str_detect(param, "^D|^N")) #extracting D.after and N samples