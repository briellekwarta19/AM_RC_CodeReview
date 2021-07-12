#libraries
library(nimble)
library(here)
library(gtools)
library(MCMCvis)
library(abind)
library(tidyverse)

#Problems with this code: 
#I need tricks to make my model run faster!
#I need tricks to make my model work with low number of individuals.
  #Specifically, I get: logProb is -Inf issues

#### Script Description ####
#The script can be divided into a few steps: 

## A. Model Set up ##
# STEP 1: Initial simulation of true abundance= where we simulate true system dynamics (R)
#We simulate population abundance and data collection (data collection feeds into step 2)

# STEP 2: First learning process= where we generate our belief of the system (NIMBLE)

# STEP 3: First decision process= where we select our management and data collection rule
#This is based on our belief of the system
#Management alternative= place 2 traps (n.trap) at the 3 segments where we believe the population is the largest       
#Data collection rule= only collect removal data (with the management alternative specified above)

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
      
      # Removal data:
      for(k in 1: K){
        #Calculating max individuals that could be trapped: 
        #traps can capture between 0 to 50 individuals
        #But if population < 50, our max capture per trap is N.belief/n.trap
        max.trap[i,j,k] <- min(50,(N.belief[i,j,k]/n.trap)) #n.trap = number of traps per segment where trapping occurs
        
        #Y.trap = number of individuals trapped and removed per trap
        Y.trap[i,j,k] ~ dunif(0, max.trap[i,j,k]) #subtract 1 here? -from mark
        
        #Y = number of individuals trapped and removed (whole number)
        Y[i,j,k] <- ceiling(Y.trap[i,j,k]*n.trap*trap.index[i,j]) 
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
N.truth <- array(0, dim = c(I,J,K)) #N.truth = true population abundance
N.truth[,1,1] <- 5000 #Let initial abundance be 1000
R <- array(0, dim = c(I,J)) #R = population remaining after removals at j

#For the first primary select 3 random segments because all initial abundances are the same
max.trap <- array(0,dim = c(I,J,K))
site.trap1 <- sample(1:I, 3) #Select 3 trapping sites for j = 1 
trap.index<- array(0, dim = c(I, J)) #Trapping index: either trap or not

for(i in 1:I){
  if(i %in% site.trap1){trap.index[i,1] <- 1} else{
    trap.index[i,1] <- 0
  }
}

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
      #Sorting the index of abundance array, first element is index (segment location) of most abundance segment
      ind <- sort((N.truth[,j,1]), index.return=TRUE, decreasing=TRUE)$ix #extracting index
      site.trap <- ind[1:3] #where trapping occurs
      
      if(i %in% site.trap){trap.index[i,2] <- 1} else{
        trap.index[i,2] <- 0
      }
      
      #Removal data
      for(k in 1: K){
        max.trap[i,j,k] <- min(50,(N.truth[i,j,k]/n.trap))
        Y.trap[i,j,k] <- runif(1,0, max.trap[i,j,k])
        Y[i,j,k] <- ceiling(Y.trap[i,j,k]*n.trap*trap.index[i,j])
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
                     n.trap = n.trap,#Number of traps
                     trap.index = trap.index, #which segments have traps
                     alpha = rep(1, I)) #Data for dirichlet process
#Data
my.data <- list(Y.trap = Y.trap) #removal data

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

#### STEP 3: First decision process ####
#DECISION RULE: place traps at the 3 segments where population is the largest

#We manage based on what we believe the population is doing
#So, we need to extract data from NIMBLE to help us decide management
#Specifically we need to know population abundance prior to our next removal at j = 3
#We need D.after[i,j] where j = 2

#Extracting information for D.after:
results <- mcmc$summary$all.chains
results <- as.data.frame(results)
results <- cbind(rownames(results), results)
colnames(results)[1] <- "param"

Nimble.results <- results %>%
  filter(str_detect(param, "^D")) #extracting D.after samples

#extracting the last D.after variables (where j = 2, D.after[i,2]
Nimble.results <- c(tail(Nimble.results$Mean, I)) 

#adding a new column to trap.index (named track.add) to represent j = 3
trap.add <- abind(trap.index, array(NA, replace(dim(trap.index), 2, 1)), along = 2)

#Sorting segments from largest to smallest population abundance
ind <- sort((Nimble.results), index.return=TRUE, decreasing=TRUE)$ix 

#Extracting the sites that have the top 3 largest abundance
site.trap <- ind[1:3] 

#Indexing the sites that trapping will occur
for(i in 1:I){
  if(i %in% site.trap){trap.add[i,J+1] <- 1} else{
    trap.add[i,J+1] <- 0
  }
}

trap.index <- trap.add

####################################################################################
#### B. General Procedure: ####
#repeats steps 4-6
final.run <- 12
# Parameter values
for (run in 3:final.run) {
  #Indices: 
  J <- run #Number of primary periods
  
  #### STEP 4: Simulate population: ####
  # Generate TRUE abundance at J
  #first add a column to each simulated data
  N.truth <- abind(N.truth, array(0, replace(dim(N.truth), 2, 1)), along = 2)
  Y.trap <- abind(Y.trap, array(0, replace(dim(Y), 2, 1)), along = 2)
  max.trap <- abind(max.trap, array(0, replace(dim(max.trap), 2, 1)), along = 2)
  Y <- abind(Y, array(0, replace(dim(Y), 2, 1)), along = 2)
  lambda <- abind(lambda, array(0, replace(dim(lambda), 2, 1)), along = 2)#assume growth is the same across space
  G <- abind(G, array(0, replace(dim(G), 2, 1)), along = 2)
  S <- abind(S, array(0, replace(dim(S), 2, 1)), along = 2)
  D <- abind(D, array(0, replace(dim(D), 2, 1)), along = 2)
  R <- abind(R, array(0, replace(dim(R), 2, 1)), along = 2)
  
  pi <- abind(pi, matrix(0, nrow = I, ncol = I), along = 3) #pi_ihj
  
  matrix(7,nrow=2,ncol=3)
  
  M <- abind(M, matrix(0, nrow = I, ncol = I), along = 3)#M_ihj
  
  pi.prob <- abind(pi.prob, matrix(0, nrow = I, ncol = I), along = 3) #pi.prob_ihj
  ind.move <- abind(ind.move, matrix(0, nrow = I, ncol = I), along = 3) #ind.move_ihj
  ind.move.last <- abind(ind.move.last, array(0, replace(dim(G), 2, 1)), along = 2) #ind.move.last_ij
  
  D.after <- abind(D.after, array(0, replace(dim(D.after), 2, 1)), along = 2)#D_hj
  
  #Simulating the truth:
  reps <- 50
  for(rep in 1:reps) {
    for(i in 1:I){
      
      ## Data Collection Model: ##
      #Removal data
      for(k in 1: K){
        max.trap[i,J,k] <- min(50,(N.truth[i,J,k]/n.trap))
        Y.trap[i,J,k] <- runif(1,0, max.trap[i,J,k])
        Y[i,J,k] <- ceiling(Y.trap[i,J,k]*n.trap* trap.index[i,J])
      }
      
      #Abundance after 1st secondary removal
      for(k in 2:K){ 
        N.truth[i,J,k] <- (N.truth[i,J,k-1] - Y[i,J,k-1])
      }
      
      #Abundance at the end of primary removal period J:
      R[i,J] <- (N.truth[i,J,K] - Y[i,J,K]) 
      
      #Abundance for 1st removal of the primary period (aka. 1st secondary removal)
      #This abundance arises from movement into the segment at the end of the previous primary period
      N.truth[i,J, 1] <- D.after[i,(J-1)]
      
      ## Population Change Model: ##  
      lambda[J] <- runif(1,0,10) #Population gains
      G[i,J] <- rpois(1,lambda[J]) #population available for dispersal 
      S[i,J] <- rbinom(1,R[i,J], phi) 
      D[i,J] <- G[i,J] + S[i,J] #population available for dispersal 
      
      ## Movement Model: ##
      #Movement probabilities:
      alpha <- rep(1:I)
      
      #Calculating movement probabilities: 
      pi[i,1:I,J] <- rdirichlet(1, alpha[1:I]) #good congugate prior for multinomial
      
      # Calculating the number of individuals that move to segment 1
      M[i,1,J] <- rbinom(1, D[i,J], pi[i,1,J])
      
      # Calculating the number of individuals that move to all other segments (not final segment)
      # prob of going to segment h given not in any previous segment
      for(h in 2:(I-1)){
        pi.prob[i,h,J] <- pi[i,h,J]/(1-sum(pi[i,1:(h-1),J])) #calculating movement probability
        ind.move[i,h,J] <- (D[i,J] - sum(M[i,(1:(h-1)),J])) #counting individuals that are moving
        M[i,h,J] <- rbinom(1,ind.move[i,h,J], pi.prob[i,h,J]) 
      }
      
      # Calculating the number of individuals that move to the final segment
      ind.move.last[i,J] <- sum(M[i,(1:(I-1)),J])
      M[i,I,J] <- (D[i,J]  - ind.move.last[i,J])
      
      D.after[i,J] <- sum(M[,i,J]) #number of individuals after movement (feeds into removal model)
      
      
    } #Ends segment loop
  }
  
  #### STEP 5: Learning process ####
  #Running Nimble to generate our belief of the system
  my.constants <- list(I = I, 
                       J = J, 
                       K = K, 
                       n.trap = n.trap,
                       trap.index = trap.index, 
                       alpha = rep(1, I))
  
  my.data <- list(Y.trap = Y.trap)
  
  
  # Initial values
  N.start <- array(c(rep(5000,I),rep(5000, I + 2*(J-1)*(I))), dim = c(I,J,K))
  initial.values <- function(){list(N.belief = N.start)}
  
  # Parameters monitored
  parameters.to.save <- c("N.belief", "lambda", "pi","D.after") 
  
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
  
  
  #### STEP 6: Update decision process ####
  #DECISION RULE: place traps at the 3 segments where population is the largest
  
  #We manage based on what we believe the population is doing
  #So, we need to extract data from NIMBLE to help us decide management
  #Specifically we need to know population abundance prior to our next removal at j = J+1
  #We need D.after[i,j] where j = J
  
  #Extracting information for D.after:
  results <- mcmc$summary$all.chains
  results <- as.data.frame(results)
  results <- cbind(rownames(results), results)
  colnames(results)[1] <- "param"
  
  Nimble.results <- results %>%
    filter(str_detect(param, "^D")) #extracting D.after samples
  
  #extracting the last D.after variables (where j = J, D.after[i,J]
  Nimble.results <- c(tail(Nimble.results$Mean, I))
  
  #adding a new column to trap.index (named track.add) to represent j = J + 1
  trap.add <- abind(trap.index, array(NA, replace(dim(trap.index), 2, 1)), along = 2)
  
  #Sorting segments from largest to smallest population abundance
  ind <- sort((Nimble.results), index.return=TRUE, decreasing=TRUE)$ix #extracting index of most abundance segment
  
  #Extracting the sites that have the top 3 largest abundance
  site.trap <- ind[1:3]
  
  #Indexing the sites that trapping will occur
  for(i in 1:I){
    if(i %in% site.trap){trap.add[i,J+1] <- 1} else{
      trap.add[i,J+1] <- 0
    }
  }
  
  trap.index <- trap.add
  
}

############################################################################

# Inspect results
#Check to make sure N[i,j,k] = D.after[i,j-1]
results <- mcmc$summary$all.chains
results <- as.data.frame(results)
results <- cbind(rownames(results), results)

colnames(results)[1] <- "param"


results %>% 
  filter(str_detect(param, "^D|^N")) #extracting D.after and N samples
