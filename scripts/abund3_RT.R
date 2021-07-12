library(tidyverse)
library(here)
library(LaplacesDemon)
library(nimble)
library(MCMCvis)
library(abind)

# Problems with this code: logProb infinite values..

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
#and collect radio telemetry data

## B. General Procedure: ##
# STEP 4: Simulate Population = simulating true system dynamics
# STEP 5: Learning process = where we generate our belief of the system using data we collected
# STEP 6: Update decision process = based on our belief, we update management and data collection rules


#For me NIMBLE only runs if I do this... if you downloaded nimble correctly you probably shouldnt do this!
path <- Sys.getenv('PATH')
newPath <- paste("C:\\Rtools\\bin;C:\\Rtools\\mingw_64\\bin;",
                 path, sep = "")
Sys.setenv(PATH = newPath)

#### NIMBLE model ####
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
        Y.trap[i,j,k] ~ dunif(0, max.trap[i,j,k])
        
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
  
    ##--Telemetry data--##
  
    for(j in 1:J){
      # trobabilities of state z(j+1) given z(j)
      gamma[1,1,j] <- phi * pi[1,1,j] #survived and stayed at 1
      gamma[1,2,j] <- phi * pi[1,2,j] #survived and moved to 2
      gamma[1,3,j] <- phi * pi[1,3,j] #survived and moved to 3
      gamma[1,4,j] <- phi * pi[1,4,j] #survived and moved to 4
      gamma[1,5,j] <- phi * pi[1,5,j] #survived and moved to 5
      gamma[1,6,j] <- 1 - phi      # died
    
      gamma[2,1,j] <- phi * pi[2,1,j] #survived and moved to 1
      gamma[2,2,j] <- phi * pi[2,2,j] #survived and stayed at 2
      gamma[2,3,j] <- phi * pi[2,3,j] #survived and moved to 3
      gamma[2,4,j] <- phi * pi[2,4,j] #survived and moved to 4
      gamma[2,5,j] <- phi * pi[2,5,j] #survived and moved to 5
      gamma[2,6,j] <- 1 - phi          #died
    
      gamma[3,1,j] <- phi * pi[3,1,j] #survived and moved to 1
      gamma[3,2,j] <- phi * pi[3,2,j] #survived and moved to 2
      gamma[3,3,j] <- phi * pi[3,3,j] #survived and stayed at 3
      gamma[3,4,j] <- phi * pi[3,4,j] #survived and moved to 4
      gamma[3,5,j] <- phi * pi[3,5,j] #survived and moved to 5
      gamma[3,6,j] <- 1 - phi         #died
    
      gamma[4,1,j] <- phi * pi[4,1,j] #survived and moved to 1
      gamma[4,2,j] <- phi * pi[4,2,j] #survived and moved to 2
      gamma[4,3,j] <- phi * pi[4,3,j] #survived and moved to 3
      gamma[4,4,j] <- phi * pi[4,4,j] #survived and stayed at 4
      gamma[4,5,j] <- phi * pi[4,5,j] #survived and moved to 5
      gamma[4,6,j] <- 1 - phi         #died
    
      gamma[5,1,j] <- phi * pi[5,1,j] #survived and moved to 1
      gamma[5,2,j] <- phi * pi[5,2,j] #survived and moved to 2
      gamma[5,3,j] <- phi * pi[5,3,j] #survived and moved to 3
      gamma[5,4,j] <- phi * pi[5,4,j] #survived and moved to 4
      gamma[5,5,j] <- phi * pi[5,5,j] #survived and stayed at 5
      gamma[5,6,j] <- 1 - phi         #died
    
      gamma[6,1,j] <- 0            #died and moved to 1 - impossible
      gamma[6,2,j] <- 0            #died and moved to 2 - impossible
      gamma[6,3,j] <- 0            #died and moved to 3 - impossible
      gamma[6,4,j] <- 0            #died and moved to 4 - impossible
      gamma[6,5,j] <- 0            #died and moved to 5 - impossible
      gamma[6,6,j] <- 1            #died and stayed dead - the only possible transition 
    
      # probabilities of y(j) given z(j)
      omega[1,1,j] <- 1 - t[1,j] # Pr(alive 1 j -> non-detected j)
      omega[1,2,j] <- t[1,j]     # Pr(alive 1 j -> detected 1 j)
      omega[1,3,j] <- 0          # Pr(alive 1 j -> detected 2 j)
      omega[1,4,j] <- 0          # Pr(alive 1 j -> detected 3 j)
      omega[1,5,j] <- 0          # Pr(alive 1 j -> detected 4 j)
      omega[1,6,j] <- 0          # Pr(alive 1 j -> detected 5 j)
    
      omega[2,1,j] <- 1 - t[2,j] # Pr(alive 2 j -> non-detected j)
      omega[2,2,j] <- 0          # Pr(alive 2 j -> detected 1 j)
      omega[2,3,j] <- t[2,j]     # Pr(alive 2 j -> detected 2 j)
      omega[2,4,j] <- 0          # Pr(alive 2 j -> detected 3 j)
      omega[2,5,j] <- 0          # Pr(alive 2 j -> detected 4 j)
      omega[2,6,j] <- 0          # Pr(alive 2 j -> detected 5 j)
    
      omega[3,1,j] <- 1 - t[3,j] # Pr(alive 3 j -> non-detected j)
      omega[3,2,j] <- 0          # Pr(alive 3 j -> detected 1 j)
      omega[3,3,j] <- 0          # Pr(alive 3 j -> detected 2 j)
      omega[3,4,j] <- t[3,j]     # Pr(alive 3 j -> detected 3 j)
      omega[3,5,j] <- 0          # Pr(alive 3 j -> detected 4 j)
      omega[3,6,j] <- 0          # Pr(alive 3 j -> detected 5 j)
    
      omega[4,1,j] <- 1 - t[4,j] # Pr(alive 4 j -> non-detected j)
      omega[4,2,j] <- 0          # Pr(alive 4 j -> detected 1 j)
      omega[4,3,j] <- 0          # Pr(alive 4 j -> detected 2 j)
      omega[4,4,j] <- 0          # Pr(alive 4 j -> detected 3 j)
      omega[4,5,j] <- t[4,j]     # Pr(alive 4 j -> detected 4 j)
      omega[4,6,j] <- 0          # Pr(alive 4 j -> detected 5 j)
    
      omega[5,1,j] <- 1 - t[5,j] # Pr(alive 5 j -> non-detected j)
      omega[5,2,j] <- 0          # Pr(alive 5 j -> detected 1 j)
      omega[5,3,j] <- 0          # Pr(alive 5 j -> detected 2 j)
      omega[5,4,j] <- 0          # Pr(alive 5 j -> detected 3 j)
      omega[5,5,j] <- 0          # Pr(alive 5 j -> detected 4 j)
      omega[5,6,j] <- t[5,j]     # Pr(alive 5 j -> detected 5 j)
    
      omega[6,1,j] <- 1          # Pr(dead j -> non-detected j)
      omega[6,2,j] <- 0          # Pr(dead j -> detected 1 j)
      omega[6,3,j] <- 0          # Pr(dead j -> detected 2 j)
      omega[6,4,j] <- 0          # Pr(dead j -> detected 3 j)
      omega[6,5,j] <- 0          # Pr(dead j -> detected 4 j)
      omega[6,6,j] <- 0          # Pr(dead j -> detected 5 j)
    }
  
    # likelihood 
    for (n in 1:N){
      # latent state at first capture
      z[n,1] <- y[n,1] - 1
    
      for (j in 2:J){
        z[n,j] ~ dcat(gamma[z[n,j-1],1:6,j]) #z[n,j,k]
        # y(j) given z(j)
        y[n,j] ~ dcat(omega[z[n,j],1:6,j])   #y[n,j,k]
      }
    }
  
  ##----------Priors----------##
  #Parameter for initial population
  psi ~ dunif(0.49,0.51) 
  #survival rate (constant through time)
  phi ~ dunif(0.7, 0.95)
  #fecundity
  for(j in 1:J){
    lambda[j] ~dunif(0,10)
  }
  
}) #end model

##########################################################################
#### A. Model set up ####
#### STEP 1a: Initial simulation of true system dynamics ####
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

#### STEP 1b: Simulate telemetry data ####
#parameters
t <- array(NA, dim = c(I,J)) #detection probability from telemetry data
alpha <- rep(1,I)

for(i in 1:I){
  for(j in 1:J){
    t[i,j] <- runif(1,0.3,0.4) #detection at each site
  }
}

#rewrite 6-10 (June -Oct) -Don't run this in this script -run later on when time is longer
# for(j in 6:10){
#   t[1:I,j] <- runif(1,0.6, 0.7) #detection at each site
# }


# probabilities of state z(j+1) given z(j)
gamma <- array(NA, dim = c(I+1,I+1,J))
# probabilities of y(j) given z(j)
omega <- array(NA, dim = c(I+1,I+1,J))

for(j in 1:J){
  # probabilities of state z(j+1) given z(j)
  gamma[1,1,j] <- phi * pi[1,1,j] #survived and stayed at 1
  gamma[1,2,j] <- phi * pi[1,2,j] #survived and moved to 2
  gamma[1,3,j] <- phi * pi[1,3,j] #survived and moved to 3
  gamma[1,4,j] <- phi * pi[1,4,j] #survived and moved to 4
  gamma[1,5,j] <- phi * pi[1,5,j] #survived and moved to 5
  gamma[1,6,j] <- 1 - phi      # died
  
  gamma[2,1,j] <- phi * pi[2,1,j] #survived and moved to 1
  gamma[2,2,j] <- phi * pi[2,2,j] #survived and stayed at 2
  gamma[2,3,j] <- phi * pi[2,3,j] #survived and moved to 3
  gamma[2,4,j] <- phi * pi[2,4,j] #survived and moved to 4
  gamma[2,5,j] <- phi * pi[2,5,j] #survived and moved to 5
  gamma[2,6,j] <- 1 - phi          #died
  
  gamma[3,1,j] <- phi * pi[3,1,j] #survived and moved to 1
  gamma[3,2,j] <- phi * pi[3,2,j] #survived and moved to 2
  gamma[3,3,j] <- phi * pi[3,3,j] #survived and stayed at 3
  gamma[3,4,j] <- phi * pi[3,4,j] #survived and moved to 4
  gamma[3,5,j] <- phi * pi[3,5,j] #survived and moved to 5
  gamma[3,6,j] <- 1 - phi         #died
  
  gamma[4,1,j] <- phi * pi[4,1,j] #survived and moved to 1
  gamma[4,2,j] <- phi * pi[4,2,j] #survived and moved to 2
  gamma[4,3,j] <- phi * pi[4,3,j] #survived and moved to 3
  gamma[4,4,j] <- phi * pi[4,4,j] #survived and stayed at 4
  gamma[4,5,j] <- phi * pi[4,5,j] #survived and moved to 5
  gamma[4,6,j] <- 1 - phi         #died
  
  gamma[5,1,j] <- phi * pi[5,1,j] #survived and moved to 1
  gamma[5,2,j] <- phi * pi[5,2,j] #survived and moved to 2
  gamma[5,3,j] <- phi * pi[5,3,j] #survived and moved to 3
  gamma[5,4,j] <- phi * pi[5,4,j] #survived and moved to 4
  gamma[5,5,j] <- phi * pi[5,5,j] #survived and stayed at 5
  gamma[5,6,j] <- 1 - phi         #died
  
  gamma[6,1,j] <- 0            #died and moved to 1 - impossible
  gamma[6,2,j] <- 0            #died and moved to 2 - impossible
  gamma[6,3,j] <- 0            #died and moved to 3 - impossible
  gamma[6,4,j] <- 0            #died and moved to 4 - impossible
  gamma[6,5,j] <- 0            #died and moved to 5 - impossible
  gamma[6,6,j] <- 1            #died and stayed dead - the only possible transition 
  
  # probabilities of y(j) given z(j)
  omega[1,1,j] <- 1 - t[1,j] # Pr(alive 1 j -> non-detected j)
  omega[1,2,j] <- t[1,j]     # Pr(alive 1 j -> detected 1 j)
  omega[1,3,j] <- 0          # Pr(alive 1 j -> detected 2 j)
  omega[1,4,j] <- 0          # Pr(alive 1 j -> detected 3 j)
  omega[1,5,j] <- 0          # Pr(alive 1 j -> detected 4 j)
  omega[1,6,j] <- 0          # Pr(alive 1 j -> detected 5 j)
  
  omega[2,1,j] <- 1 - t[2,j] # Pr(alive 2 j -> non-detected j)
  omega[2,2,j] <- 0          # Pr(alive 2 j -> detected 1 j)
  omega[2,3,j] <- t[2,j]     # Pr(alive 2 j -> detected 2 j)
  omega[2,4,j] <- 0          # Pr(alive 2 j -> detected 3 j)
  omega[2,5,j] <- 0          # Pr(alive 2 j -> detected 4 j)
  omega[2,6,j] <- 0          # Pr(alive 2 j -> detected 5 j)
  
  omega[3,1,j] <- 1 - t[3,j] # Pr(alive 3 j -> non-detected j)
  omega[3,2,j] <- 0          # Pr(alive 3 j -> detected 1 j)
  omega[3,3,j] <- 0          # Pr(alive 3 j -> detected 2 j)
  omega[3,4,j] <- t[3,j]     # Pr(alive 3 j -> detected 3 j)
  omega[3,5,j] <- 0          # Pr(alive 3 j -> detected 4 j)
  omega[3,6,j] <- 0          # Pr(alive 3 j -> detected 5 j)
  
  omega[4,1,j] <- 1 - t[4,j] # Pr(alive 4 j -> non-detected j)
  omega[4,2,j] <- 0          # Pr(alive 4 j -> detected 1 j)
  omega[4,3,j] <- 0          # Pr(alive 4 j -> detected 2 j)
  omega[4,4,j] <- 0          # Pr(alive 4 j -> detected 3 j)
  omega[4,5,j] <- t[4,j]     # Pr(alive 4 j -> detected 4 j)
  omega[4,6,j] <- 0          # Pr(alive 4 j -> detected 5 j)
  
  omega[5,1,j] <- 1 - t[5,j] # Pr(alive 5 j -> non-detected j)
  omega[5,2,j] <- 0          # Pr(alive 5 j -> detected 1 j)
  omega[5,3,j] <- 0          # Pr(alive 5 j -> detected 2 j)
  omega[5,4,j] <- 0          # Pr(alive 5 j -> detected 3 j)
  omega[5,5,j] <- 0          # Pr(alive 5 j -> detected 4 j)
  omega[5,6,j] <- t[5,j]     # Pr(alive 5 j -> detected 5 j)
  
  omega[6,1,j] <- 1          # Pr(dead j -> non-detected j)
  omega[6,2,j] <- 0          # Pr(dead j -> detected 1 j)
  omega[6,3,j] <- 0          # Pr(dead j -> detected 2 j)
  omega[6,4,j] <- 0          # Pr(dead j -> detected 3 j)
  omega[6,5,j] <- 0          # Pr(dead j -> detected 4 j)
  omega[6,6,j] <- 0          # Pr(dead j -> detected 5 j)
}


# likelihood 
N <- 20 #n individuals
#K <- 2 #n occasions
z <- array(1,dim = c(N,J))
y <- array(1,dim = c(N,J))
first <- rep(1, N) # single cohort

for (n in 1:N){
  # latent state at first capture
  z[n,first[n]] <- y[n,first[n]] - 1
  for (j in 1:J){
    z[n,j] <- rcat(1,gamma[z[n,j-1],1:6,j]) #z[n,j,k]
    # y(j) given z(j)
    y[n,j] <- rcat(1,omega[z[n,j],1:6,j])   #y[n,j,k]
  }
}

#############################################################################
#### STEP2: First Learning ####
# Data in a list. 
my.data <- list(y = y + 1, Y = Y)

my.constants <- list(I = I, #Number of segments
                     J = dim(y)[2], #Number of primary periods
                     K = K, #Number of secondary periods
                     N = nrow(y),
                     t = t,#Number of i ndividuals tracked in telemetry data
                     n.trap = n.trap,#Number of traps
                     trap.index = trap.index, #which segments have traps
                     alpha = c(1, 1, 1, 1, 1))

#Initial values. 
z.start <- y
z.start[z.start==0] <- sample(c(1,2,3,4,5), sum(z.start==0), replace = TRUE)

t.start <- array(NA, dim = c(I,J))
for(i in 1:I){
  for(j in 1:J){
    t.start[i,j] <- runif(1,0.3, 0.4)
  }
}

#rewrite 6-10 (June -Oct) #ignore this line for now... don't have enough timesteps yet
# for(j in 6:10){
#   t.start[1:I,j] <- runif(1,0.6, 0.7) #detection at each site
# }

N.start <- array(c(rep(5000,I),rep(5000, I + 2*(J-1)*(I))), dim = c(I,J,K))

initial.values <- function(){list(phi = runif(1, 0, 1), 
                                  pi = array(1/I, dim = c(I,I,J)),
                                  z = z.start, 
                                  N.belief = N.start)}

# Parameters monitored
parameters.to.save <- c("phi", "pi", "N.belief", "lambda", "pi","D.after") #should be D.after

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
