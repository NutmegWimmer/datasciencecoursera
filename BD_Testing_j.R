rm(list = ls()) # clear R Environment

# load required libraries
library("purrr")
library("dplyr") 

# runSEIQR takes the SEIR model parameters and initial condition 
# and returns a time series for each group
runSEIQR <- function(beta, # transmission rate
                    sigma, # 1/incubation time 
                    gamma, # recovery rate
                    theta, # test turn around time (e.g., 2 days)
                    omega, # 1/time in quarantine 
                    tau.s, # proportion of S tested
                    tau.e, # proportion of E tested
                    tau.i, # proportion of I tested
                    tau.r, # proportion of R tested
                    phi.p, # specificity (true negative rate)
                    phi.e, # sensitivity (true positive rate)
                    phi.a, # specificity w.r.t. culture (?)
                    q,     # reduction in transmissibility for infected awaiting test results
                    initial.state, 
                    step.size = 1,
                    freq.dependent=TRUE,
                    final.only=FALSE) {


  #If the model is frequency dependent we modify beta based on the total population size
  beta.divisor <- ifelse(freq.dependent,
                         initial.state["S"]+initial.state["E"]+initial.state["I"]+initial.state["Q"]+initial.state["R"]+
                         initial.state["Ts"]+initial.state["Te"]+initial.state["Ti"]+initial.state["Tr"],
                         1)

  #create the parameter vector.
  param <- c(beta=beta/beta.divisor,
             sigma=sigma, gamma=gamma, theta=theta, omega=omega, tau.s=tau.s, tau.e=tau.e, 
             tau.i=tau.i, tau.r=tau.r, theta=theta, phi.p=phi.p, phi.e=phi.e, phi.a=phi.a, q=q)
  
  #create the true positive reporting vector
  
  #Since we are not using a fancy solver we will need to run this on our own.
  #note that the simulation ends once there are 0 people in groups E or I
  t <- 0
  y <- initial.state
  
  seir.output <- matrix(ncol=13, nrow=1)
  colnames(seir.output) <- c("time", "S", "E", "I","Q","R","Ts", "Te", "Ti", "Tr","incident","missed","total.tests")
  seir.output[1,] <- c(t,y)
  
  ## add a thing for retesting # flagged; didn't check this
  if(test.type=="retest negatives"){
    retest<-data.frame(time=0:1000,s=0,e=0,i=0,r=0)
  }
  
  while (y["E"]>0 | y["I"]>0 | y["Te"]>0 | y["Ti"]>0) {
  #while (t < 360){  # if you comment the line above, and uncomment this line, it will run simulation for fixed time
    
    t <- t+step.size 
    print(t)
   
    #calculate the probability of expose, infection, and recovery in this time step
    # flagged: in this model, is it ok that someone in T_E always either goes to Q or E, and never to I? Presumably they could have become infectious while waiting for test results?
    ### SEIR disease course ### 
    p.expose <- 1-exp(-step.size*(param["beta"]*y["I"] + param["beta"]*param["q"]*y["Ti"])) # flagged: I changed this line 
    p.infect <- 1-exp(-step.size*param["sigma"])
    p.recover <- 1-exp(-step.size*param["gamma"])
    # added symptomatic
    p.symptomatic<- 0.6 #1-exp(-step.size*.6) ## probability of being symptomatic # flagged: I changed this to just an actual probability
    
    ## Testing; probability that a person is tested in each state
    p.test.s <- 1-exp(-step.size*param["tau.s"]) 
    p.test.e <- 1-exp(-step.size*param["tau.e"])
    p.test.i <- 1-exp(-step.size*param["tau.i"])
    p.test.r <- 1-exp(-step.size*param["tau.r"])

    ## Probability of people getting test results # flagged; I changed all of these from what was here before 
    # see note below about how this section simultaneously considers two random processes:
    # testing positive/negative and also the number of people getting test results back
    p.positive.s <- (1-param["phi.p"])*(1-exp(-step.size*param["theta"])) ## probability of a person testing + but who is susceptible (false positives)
    p.negative.s <- param["phi.p"]*(1-exp(-step.size*param["theta"])) ## probability of a person testing - who is susceptible (true negatives)
    p.positive.e <- (1-param["phi.p"])*(1-exp(-step.size*param["theta"])) ## true positives
    p.negative.e <- param["phi.p"]*(1-exp(-step.size*param["theta"])) ## false negatives
    p.positive.i <- param["phi.e"]*(1-exp(-step.size*param["theta"])) ## true positives
    p.negative.i <- (1-param["phi.e"])*(1-exp(-step.size*param["theta"])) ## false negatives
    p.positive.r <- (1-param["phi.p"])*(1-exp(-step.size*param["theta"])) ## false positives
    p.negative.r <- param["phi.p"]*(1-exp(-step.size*param["theta"])) ## true negatives
   
    ## probability of leaving/ending quarantine
    p.end.isolate <- 1-exp(-step.size*param["omega"]) ## probability of quarantined person moving to recovered
    
    ## Tested disease course (same as p.expose above, except with an additional q)
    p.test.expose <- 1-exp(-step.size*(param["q"]*param["beta"]*y["I"] + param["q"]*param["beta"]*param["q"]*y["Ti"])) # flagged: I changed this line to be consistent with beta definition at top
    p.test.infect <- 1-exp(-step.size*param["sigma"])
    p.test.recover <- 1-exp(-step.size*param["gamma"])
    
    p.test.symptomatic<-1-exp(-step.size*.6) ## probability of being symptomatic; # flagged; didn't check this
    
    #draw random variable from binomial distribution for new number of exposed, incident, and recovered cases
    ### Not tested disease course
    exposed.cases <- rbinom(1, y["S"],p.expose) # number of exposed people in this time step
    incident.cases <- rbinom(1, y["E"],p.infect)
    
    # added symptomatic classification for reflexing; # flagged; didn't check this
    symptomatic.i.cases<-rbinom(1, y["E"],p.infect*p.symptomatic)
    asymptomatic.i.cases<-rbinom(1, y["E"],p.infect*(1-p.symptomatic))
    
    recovered.cases <- rbinom(1, y["I"], p.recover)

    ## Testing
    tested.susceptibles <- rbinom(1, y["S"]-exposed.cases, p.test.s) ### number of susc. people who are tested
    tested.exposed <- rbinom(1, y["E"]-incident.cases, p.test.e)
    tested.infected <- rbinom(1, y["I"]-recovered.cases, p.test.i)
    tested.recovered <- rbinom(1, y["R"], p.test.r)
    
    ## but overwrite testing if we are retesting negatives; # flagged; didn't check this
    if(test.type=="retest negatives"){
      if(t>2){ #this line to avoid error when pulling from negative day
      ## add in cases to retest from each compartment from two days ago
      tested.susceptibles <- rbinom(1, y["S"]-exposed.cases, p.test.s)+retest$s[t-2] ### number of susc. people who are tested
      tested.exposed <- rbinom(1, y["E"]-incident.cases, p.test.e)+retest$e[t-2]
      tested.infected <- rbinom(1, y["I"]-recovered.cases, p.test.i)+retest$i[t-2]
      tested.recovered <- rbinom(1, y["R"], p.test.r)+retest$r[t-2]
      }
    }
    
    ## Tested disease course
    exposed.tested.cases <- rbinom(1, y["Ts"],p.test.expose)
    incident.tested.cases <- rbinom(1, y["Te"],p.test.infect)
    recovered.tested.cases <- rbinom(1, y["Ti"], p.test.recover)
    
    ## added symptomatic classification for reflexing; # flagged; didn't check this
    symptomatic.i.tested.cases<-rbinom(1, y["Te"],p.test.infect*p.test.symptomatic)
    asymptomatic.i.tested.cases<-rbinom(1, y["Te"],p.test.infect*(1-p.test.symptomatic))
    
    
    ## Test results - n individuals in each time step # flagged: note, in the code below and related code above, 
    # the two stochastic processes of getting test results back and then whether those results were positive or negative were treated as
    # one random process; this should probably be changed to consider them seperately 
    positive.susceptibles <- rbinom(1, y["Ts"]-exposed.tested.cases, p.positive.s) ## false positives # flagged: I think there was a typo in this line before, should be good now
    positive.exposed <- rbinom(1, y["Te"]-incident.tested.cases, p.positive.e) ## true positives
    positive.infected <- rbinom(1, y["Ti"]-recovered.tested.cases, p.positive.i) ## true positives
    positive.recovered <- rbinom(1, y["Tr"], p.positive.r) ## false positives... but kinda not? This is people who were infected, made it to recovered, then were tested and tested positive.

    # flagged: I changed the following four lines, removing the "-positive.susceptibles" terms as I think they skewed the results in a weird way and those two processes should happen simultaneously
    negative.susceptibles <- rbinom(1, y["Ts"]-exposed.tested.cases, p.negative.s) ## true negatives
    negative.exposed <- rbinom(1, y["Te"]-incident.tested.cases, p.negative.e) ## false negatives
    negative.infected <- rbinom(1, y["Ti"]-recovered.tested.cases, p.negative.i) ## false negatives
    negative.recovered <- rbinom(1, y["Tr"], p.negative.r) ## true negatives.... but actually not. 
    ## This is people who were infected, made it to recovered, then were tested and tested negative. 
    ## So actually this could be considered a false negative... but it really depends on what your basis for truth is. 
    ## A recovered case is not infectious, so if you care about finding infectious cases, then you have missed nothing here. 
    ## But according to the guidelines set by Megan and Lindsay yesterday, this is one part of the "missed cases"
    
    ##### EDIT \/ \/ \/ \/ 
    missed.cases<-negative.recovered+recovered.cases ## flagged: THIS WON'T WORK CORRECTLY! see note on next line:
    # since people in class R can be tested theoretically many times, every time they reenter class R, they get counted again as a missed case
    # so your missed.cases numbers were WAY BIGGER than they should have been, especially late in the simulation
    # when there are lots of people in the R bin. 

    ## when retesting negatives, fill row in retest data frame with numbers of seir who need retesting after delay
    # flagged; didn't check this
    if(test.type=="retest negatives"){
    retest$s[t]<-negative.susceptibles
    retest$e[t]<-negative.exposed
    retest$i[t]<-negative.infected
    retest$r[t]<-negative.recovered
    }

    ## End quarantine
    end.isolated.cases <- rbinom(1, y["Q"], p.end.isolate) ## people who have lived out quarantine
    
    #Find the deltas for each compartment
    dS <- -exposed.cases - tested.susceptibles + negative.susceptibles
    dE <- exposed.cases - incident.cases - tested.exposed + negative.exposed
    dI <- incident.cases - recovered.cases - tested.infected + negative.infected
    dR <- recovered.cases + end.isolated.cases - tested.recovered + negative.recovered ### 
    dQ <- positive.susceptibles +  positive.exposed + positive.infected + positive.recovered - end.isolated.cases 
    
    dTs <- -exposed.tested.cases + tested.susceptibles - positive.susceptibles - negative.susceptibles
    dTe <- exposed.tested.cases  - incident.tested.cases + tested.exposed - positive.exposed - negative.exposed
    dTi <- incident.tested.cases - recovered.tested.cases + tested.infected - positive.infected - negative.infected
    dTr <- recovered.tested.cases + tested.recovered - positive.recovered - negative.recovered
      
    ##### Lindsay had a line to calculate change in incident, But it seems not useful for our  current purpose (we just want to report how many of each)
    # dc <- -y["incident"] ## subtract previous incident cases (so that they don't just add up)??
    # dm <- -y["missed"]
    # dtest<- -y["total.tests"]
    
    delta <- c(dS, dE, dI, dQ, dR, dTs, dTe, dTi, dTr, NA, NA, NA) # calculate step sizes ### NAs will be replaced in immediately following lines
    
    y <- y+delta
    y["incident"] <- incident.cases ## pull out day's new cases into vector; flagged; not sure if this is correctly defined; seems a bit weird to me but not sure what you're trying to do with it
    y["missed"] <- missed.cases ## 
    ## pull out number of tests
    y["total.tests"]<-tested.susceptibles+tested.exposed+tested.infected+tested.recovered

    
    print(y)
    # 
#trick to speed up the code 
    if (!final.only) {
      seir.output<-rbind(seir.output, c(t,y))
    }
    
    
  }############## end loop of days
  
  if(final.only) {
    seir.output[1,]<-c(t,y)
  }
  
  return(seir.output)
  
}################ end function

## Run this model with parameters for the following situations:
## 1) nursing home (86 residents + .18 staff per resident = 102)
#(86*1.18)

## 2) dorm ### size to come
##
## and track things with cost: 
## a) infections
## b) true and false quaratines
## c) 



# create empty storage objects
output.table <- list()
used.param <- list()
results<-list()
rep <- 1 # number of stochastic simulations you want to run

for(i in 1:rep){
  site<-"not NH"
  test.type<-"PCR" #"retest negatives" 
  pop.size<-ifelse(site=="NH",round(86*1.18),10000) # total population size ## change 1000 to actual pop in dorm
  n.seed.events <- 20 # number of infected people at time t=0
  initial.state <- c(S=pop.size-n.seed.events, E=0, I=n.seed.events, Q=0, R=0, 
                     Ts= 0, Te=0, Ti=0, Tr=0,
                     incident=0, missed=0,total.tests=0)
  
  # define SEIR parameters
  #note: if you want them to be pulled from a distribution, 
  #replace values with the rnorm code included as a comment on each line
  
  sigma <- 1/5.2  ## incubation (1/days to symptom onset) #rnorm(1,mean=1/5.2,sd=0.01) 
  gamma <- 1/runif(1,10,21)  ## recovery rate (1/infectious period) #rnorm(1,mean=1/6,sd=0.01) #flagged
  R0 <- runif(1,1.2,1.5)   #rnorm(1,mean=2,sd=0.01) 
  beta <- gamma*R0
  omega <- 1/14 ## quarantine factor (1/time in quarantine) ## also simulate with 10 days
  theta <- 1/ifelse(test.type=="PCR",3,0) ## test turn around (1/processing time for test) ## this will change by test type
  # flagged: if test.type isn't "PCR", then theta=inf; double check that is what you want
  tau.s <- -log(1-0.01) ## -log(proportion of state compartment tested) \/
  tau.e <- -log(1-0.01) ## "
  tau.i <- -log(1-0.01) ## "
  tau.r <- -log(1-0.01) ## "
  phi.p <- 0.99 ## specificity ## will change by test type: ifelse(test.type=="PCR",something,somethingelse)
  phi.e <- 0.99 ## sensitivity ## will change by test type
  phi.a <- 0.80 ## specificity with respect to culture ## will change by test type
  q  <- 0.5 ## reduction in beta due to testing ## try with .25 and .75 as well 
  

  # run the SEIR model function defined above
  output.table[[i]] <- as.data.frame(runSEIQR(
                                              beta = beta, 
                                              sigma = sigma,
                                              gamma = gamma,
                                              theta = theta, 
                                              omega = omega, 
                                              tau.s = tau.s, 
                                              tau.e = tau.e, 
                                              tau.i = tau.i,
                                              tau.r = tau.r, 
                                              phi.p = phi.p, 
                                              phi.e = phi.e, 
                                              phi.a = phi.a, 
                                              q = q, 
                                              initial.state = initial.state,
                                              step.size = 1,
                                              freq.dependent=TRUE,
                                              final.only=FALSE
                                              ))
  
  used.param[[i]] <- c(beta,sigma, gamma, theta, omega, tau.s, tau.e, tau.i, tau.r, theta, phi.p, 
                       phi.e, phi.a, q)
   
  output.table[[i]]<- output.table[[i]]


results<-output.table 

} ## end for loop of simulations


 
# make a new plot for each simulation
for(i in 1:rep){
  results<-output.table
  plot(x=results[[i]]$time, y=results[[i]]$S, type="n",
       xlab="time (days)", ylab="Number of People",ylim=c(0,pop.size),
       main="COVID-19 SEIR \n \n \n")
  lines(results[[i]][,2], col=1)
  lines(results[[i]][,3], col=2)
  lines(results[[i]][,4], col=3)
  lines(results[[i]][,5], col=4)
  lines(results[[i]][,6], col=5)
  lines(results[[i]]$incident,col="purple" )
  legend("right", legend=c("S", "E", "I","Q", "R"), col=1:5, lty=1)
  title(main=paste("beta=",round(used.param[[i]][1],2)," sigma=",round(used.param[[i]][2],2)," \n gamma=",round(used.param[[i]][3],2),sep=""), line=0.75)
}


## take a look here::
results1<-results[[1]]
plot((incident-missed)~time,results1,xlim=c(0,100))
plot(missed~time,results1,xlim=c(0,100))
#View(results1)
total.infections<-sum(results1$incident)
total.infections
days.quarantine.necessary<-total.infections*(1/omega)
days.quarantine.actual<-sum(results1$Q)
days.quarantine.unnecessary<-days.quarantine.actual-days.quarantine.necessary
days.quarantine.unnecessary
total.hospitalizations<-total.infections*0.039 ## hospitalization rate
total.icu<-total.hospitalizations*0.38 ## ICU entrance rate (of hospitalized cases)
total.deaths<-total.infections*0.0039 ## death rate
total.deaths
total.missed.infections<-sum(results1$missed)
total.missed.infections
total.tests<-sum(results1$total.tests)
total.tests
plot(total.tests~time,results1)
plot(missed~time,results1,xlim=c(0,100))
plot(incident~time,results1,col="red",pch=16)




# 
# ####### mess around with some stuff::
# sigma <- 1/5.2   #rnorm(1,mean=1/5.2,sd=0.01) 
# gamma <- 1/14   #rnorm(1,mean=1/6,sd=0.01) 
# R0 <- 1   #rnorm(1,mean=2,sd=0.01) 
# beta <- gamma*R0
# theta <- 1/3
# omega <- 1/14
# tau.s <- 0.01
# tau.e <- 0.01
# tau.i <- 0.01
# tau.r <- 0.01
# phi.p <- 0.99
# phi.e <- 0.99
# phi.a <- 0.80
# q  <- 0.25
# n.seed.events<-10 # number of infected people at time t=0
# pop.size<-10000
# initial.state<-c(S=pop.size-n.seed.events, E=0, I=n.seed.events, Q=0, R=0, Ts= 0, Te=0, Ti=0, Tr=0)
# runSEIQR(
#   beta = beta, 
#   sigma = sigma,
#   gamma = gamma,
#   theta = theta, 
#   omega = omega, 
#   tau.s = tau.s, 
#   tau.e = tau.e, 
#   tau.i = tau.i,
#   tau.r = tau.r, 
#   phi.p = phi.p, 
#   phi.e = phi.e, 
#   phi.a = phi.a, 
#   q = q, 
#   initial.state = initial.state,
#   step.size = 1,
#   freq.dependent=TRUE,
#   final.only=FALSE
# )
# 
# 
# ## Assumptions for campus::
# 
# # - Hospitalization rate of 3.9%
# # - ICU Rate of 38% of Hospitalizations
# # - Vent Rate of 60% of ICUs
# # - Death Rate of 0.39% of infections (10x hospitalized die)
# # - Transmission rate for online classes based on July transmssion rate for Utah: Uniform(1.04,1.24)
# # - Transmission rate for in person classes based on mask wearing: Uniform(0.375,1.66)
# 
# ################# finding out how many people are in ICU/on ventilator/dead at end
# #.0399 Derived from U School Data and UDOH hospitalization Rate
# # df_final$icu <- rbinom(n = nrow(df_final),size = df_final$hosp,prob = .38)
# # df_final$vent <- rbinom(n = nrow(df_final),size = df_final$icu,prob = .60)
# # df_final$death <- rbinom(n = nrow(df_final),size = df_final$Incident,prob = .0039) #Half a percent from Lindsay's source...
# 






