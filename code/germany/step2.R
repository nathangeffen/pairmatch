
if (Sys.info()['sysname'] == 'Linux'){path <- paste(getwd(), '/', sep='')}
if (Sys.info()['sysname'] == 'Windows'){path <- paste(as.vector(strsplit(getwd(), '/')[[1]][1:length(strsplit(getwd(), '/')[[1]])-1]), '/', sep="", collapse="")}

##### READ IN DATA #####
### CORRECTED VALUES OF SINGLES AND AGENTS
data <- read.csv2(paste(path, "data.csv", sep=""), header=T)
singles <- read.csv(paste(path, "singles.csv", sep=""), header=T)
if (as.integer(commandArgs()[7]) == 0) {
  partners <- read.csv(paste(path, "partners.csv", sep=""), header=T)
} else {
  partners <- read.csv(paste(path, "partnersACP.csv", sep=""), header=T)
}

### MATRICES OF AGE DIFFERENCES
mat.mm <- read.csv2(paste(path, "mat.msm.csv", sep=""), header=F)
mat.ww <- read.csv2(paste(path, "mat.wsw.csv", sep=""), header=F)
mat.mw <- read.csv2(paste(path, "mat.msw.csv", sep=""), header=F)
mat.wm <- read.csv2(paste(path, "mat.wsm.csv", sep=""), header=F)

### DETERMINE GENERAL PARAMETERS
N <- dim(partners)[1] # Number of Age Groups in partner dataset
X <- as.integer(commandArgs()[6]) # Number of agents to be simulated
S <- round(X*sum(data$femratio*(1-data$relw_share)*data$ageshare + (1-data$femratio)*(1-data$relm_share)*data$ageshare)) # Calculate number of single agents
set.seed(5)

### CREATE EMPTY DATAFRAME FOR AGENTS
pop <- matrix(ncol = 10, nrow = X, dimnames=list(1:X, c("id", "age", "sex", "sexor", "rel", "pid", "page", "psex", "psexor", "page1")))


# Function for creating agents, with the following arguments:
#   pop       : (empty) matrix to store the values of agents created by the function
#   fromAgent : Index of agents to start producing single agent value sets (e.g. if singles are created first, = 1)
#   toAgent   : Index of agents to stop producing single agent value sets (e.g. total number of singles + fromAgent)
#   ageRange  : Range of ages agents can be created
#   ageShare  : Vector of probabilities for being in age a (same length as ageRange)
#   femRatio  : Vector of probabilities for being female given age a (same length as ageRange)
#   wswRate   : Vector of probabilities for female being homosexual at age a (same length as ageRange)
#   msmRate   : Vector of probabilities for male being homosexual at age a (same length as ageRange)

createSingles <- function(pop, fromAgent, toAgent,  ageRange, ageShare, femRatio, wswRate, msmRate){
  ### TAKE TIME
  tic <- Sys.time()
  i <- fromAgent
  ### CREATE SINGLE AGENTS
  repeat{
    # assign id:
    pop[i,1] <- i
    # sample age:
    pop[i,2] <- sample(ageRange, 1, replace = T, prob = ageShare)
    # sample sex dependent on age (0:=female, 1:=male):
    pop[i,3] <- sample(c(0,1), 1, replace = T, prob = c(femRatio[pop[i,2]+1],(1 - femRatio[pop[i,2]+1])))
    # sample sexor dependent on sex
    pop[i,4] <- if (pop[i,3] == 0) {sample(c(1,0), 1, replace = T, prob = c(wswRate[pop[i,2]+1], (1 - wswRate[pop[i,2]+1])))}
                else {sample(c(1,0), 1, replace=T, prob=c(msmRate[pop[i,2]+1], (1 - msmRate[pop[i,2]+1])))}
    # relationship status
    pop[i,5] <- 0
    # values of non-existent partner
    pop[i,6] <- 0
    pop[i,7] <- 0
    pop[i,8] <- 0
    pop[i,9] <- 0
    pop[i,10] <- NA
    i <- i + 1
    if (i > toAgent) break
  }
  cat(paste("The creation of", toAgent - fromAgent, "single agents took", round(Sys.time()-tic, 2), "seconds, \n"))
  return(pop)
}

# Function for creating agents, with the following arguments:
#   pop       : (empty) matrix to store the values of agents created by the function
#   fromAgent : Index of agents to start producing single agent value sets (e.g. if singles are created first, = 1)
#   toAgent   : Index of agents to stop producing single agent value sets (e.g. total number of singles + fromAgent)
#   ageRange  : Range of ages agents can be created
#   ageShare  : Vector of probabilities for being in age a
#   femRatio  : Vector of probabilities for being female given age a
#   wswRate   : Vector of probabilities for female being homosexual at age a
#   msmRate   : Vector of probabilities for male being homosexual at age a
#   matWW     : Matrix of probabilities for partner of homosexual, female agent of age a, being 0 to 100 years
#   matMW     : Matrix of probabilities for partner of heterosexual, male agent of age a, being 0 to 100 years
#   matWM     : Matrix of probabilities for partner of heterosexual, female agent of age a, being 0 to 100 years
#   matMM     : Matrix of probabilities for partner of homosexual, male agent of age a, being 0 to 100 years

# This function creates agents in partnerships, who do not know who their partner is
createPartnersMatch <- function(pop, fromAgent, toAgent, ageRange,  ageShare, femRatio, wswRate, msmRate, matWW, matMW, matWM, matMM){
  ### TAKE TIME
  tic <- Sys.time()
  i <- fromAgent
  ### CREATE PARTNER AGENTS
  repeat{
    # assign id:
    pop[i,1]<-i
    # sample age:
    pop[i,2] <- sample(ageRange, 1, replace = T, prob = ageShare)
    # sample sex dependent on age (0:=female, 1:=male):
    pop[i,3] <- sample(c(0,1), 1, replace = T, prob = c(femRatio[pop[i,2]-11], (1 - femRatio[pop[i,2]-11])))
    # sample sexor dependent on sex (0:=heterosexual, 1:=homosexual)
    if (pop[i,3]==0) {pop[i,4] <- sample(c(1,0), 1, replace = T, prob = c(wswRate[pop[i,2]-11], (1 - wswRate[pop[i,2]-11])))}
    if (pop[i,3]==1) {pop[i,4] <- sample(c(1,0), 1, replace = T, prob = c(msmRate[pop[i,2]-11], (1 - msmRate[pop[i,2]-11])))}
    # relationship status
    pop[i,5] <- 1
    # Set no partner id
    pop[i, 6] <- NA
    # preferences in partnership
    # preference in age corresponding to the dependent age distribution dependent on sex and sexual orientation:
    if (pop[i,3] == 0 & pop[i,4] == 1){pop[i,7] <- sample(c(12:100), 1, replace = T, matWW[,pop[i,2]-11])}
    if (pop[i,3] == 0 & pop[i,4] == 0){pop[i,7] <- sample(c(12:100), 1, replace = T, matMW[,pop[i,2]-11])}
    if (pop[i,3] == 1 & pop[i,4] == 0){pop[i,7] <- sample(c(12:100), 1, replace = T, matWM[,pop[i,2]-11])}
    if (pop[i,3] == 1 & pop[i,4] == 1){pop[i,7] <- sample(c(12:100), 1, replace = T, matMM[,pop[i,2]-11])}
    # preference in sex corresponding to the agents sexual orientation
    if (pop[i,4] == 0) {pop[i,8] <- abs(1 - pop[i,3])}
    if (pop[i,4] == 1) {pop[i,8] <- pop[i,3]}
    # set partner sexor corresponding to the agents sexual orientation
    pop[i,9] <- pop[i,4]
    pop[i,10] <- pop[i,7]
    i <- i + 1
    if (i > toAgent) break
  }
  cat(paste("The creation of", toAgent - fromAgent, "agents in partnerships took", round(Sys.time()-tic, 2), "seconds, \n"))
  return(pop)
}

# This function creates agents in partnerships, who already know their partner
createPartnersACP <- function(pop, fromAgent, toAgent, ageRange,  ageShare, femRatio, wswRate, msmRate, matWW, matMW, matWM, matMM){
  ### TAKE TIME
  tic <- Sys.time()
  i <- fromAgent
  ### CREATE PARTNER AGENTS
  repeat{
    if (i + 1 > toAgent) break
    # assign id:
    pop[i,1]<-i
    # sample age:
    pop[i,2] <- sample(ageRange, 1, replace = T, prob = ageShare)
    # sample sex dependent on age (0:=female, 1:=male):
    pop[i,3] <- sample(c(0,1), 1, replace = T, prob = c(femRatio[pop[i,2]-11], (1 - femRatio[pop[i,2]-11])))
    # sample sexor dependent on sex
    if (pop[i,3]==0) {pop[i,4] <- sample(c(1,0), 1, replace = T, prob = c(wswRate[pop[i,2]-11], (1 - wswRate[pop[i,2]-11])))}
    if (pop[i,3]==1) {pop[i,4] <- sample(c(1,0), 1, replace = T, prob = c(msmRate[pop[i,2]-11], (1 - msmRate[pop[i,2]-11])))}
    pop[i,5] <- 1
    # assign partner id of the partner to be created in the next steps
    pop[i,6] <- i+1
    # create partner
    # assign id:
    pop[i+1,1] <- i+1
    # sexural orientation is the same as the partners' sexual orientation
    pop[i+1,4] <- pop[i,4]
    # partners sex corresponding to the creators sexual orientation
    if (pop[i,4] == 1) {pop[i+1,3] <- pop[i,3]}
    if (pop[i,4] == 0) {pop[i+1,3] <- abs(1 - pop[i,3])}
    # partners age corresponding to the dependent age distribution dependent on sex and sexual orientation:
    if (pop[i+1,3] == 0 & pop[i+1,4] == 1){pop[i+1,2] <- sample(c(12:100), 1, replace = T, matWW[,pop[i,2]-11])}
    if (pop[i+1,3] == 0 & pop[i+1,4] == 0){pop[i+1,2] <- sample(c(12:100), 1, replace = T, matWM[,pop[i,2]-11])}
    if (pop[i+1,3] == 1 & pop[i+1,4] == 0){pop[i+1,2] <- sample(c(12:100), 1, replace = T, matMW[,pop[i,2]-11])}
    if (pop[i+1,3] == 1 & pop[i+1,4] == 1){pop[i+1,2] <- sample(c(12:100), 1, replace = T, matMM[,pop[i,2]-11])}
    # partner is always in a relationship
    pop[i+1,5] <- 1
    # partner information:
    pop[i+1,6] <- pop[i,1]
    # partner age
    pop[i+1,7] <- pop[i,2]
    pop[i,7] <- pop[i+1,2]
    # partner sex
    pop[i+1,8] <- pop[i,3]
    pop[i,8] <- pop[i+1,3]
    # partner sexor
    pop[i+1,9] <- pop[i,4]
    pop[i,9] <- pop[i+1,4]
    # Unneccessary variable for this version, set to NA
    pop[i,10] <- NA
    pop[i+1,10] <- NA
    i <- i+2
  }
  cat(paste("The creation of", toAgent - fromAgent, "agents in partnerships took", round(Sys.time()-tic, 2), "seconds, \n"))
  return(pop)
}

test <- createSingles(pop, 1, S, seq(0, 100, 1), singles$ageshare, singles$femratio, singles$wsw_rate, singles$msm_rate)

if (as.integer(commandArgs()[7]) == 0) {
  test <- createPartnersMatch(test, S, X, seq(12, 100, 1), partners$ageshare, partners$femratio, partners$wsw_rate, partners$msm_rate, mat.ww, mat.mw, mat.wm, mat.mm)
} else {
  test <- createPartnersACP(test, S, X, seq(12, 100, 1), partners$ageshare, partners$femratio, partners$wsw_rate, partners$msm_rate, mat.ww, mat.mw, mat.wm, mat.mm)
}

pop <- as.data.frame(test)
write.csv(pop, 'agent_list.csv', row.names=F)
