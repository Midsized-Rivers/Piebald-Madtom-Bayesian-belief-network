################################ Notes #########################################
## This script recreates the Bayesian belief network used in Dunn et al. for Piebald Madtom

#### !! set working director !!! 
workingdirectory <-  " " # set filepath to working director between quotations
setwd(workingdirectory)


### Install and load packages
packages <- c("data.table", "openxlsx")

installed_packages <- packages %in% rownames(installed.packages())
if(any(installed_packages == FALSE)) {
  install.packages(packages[!installed_packages])
}
lapply(packages, require, character.only = T)

### Functions
enumerate_states <- function(x){1:nrow(x)}


#### indexing for nested loop
### ecological settings
e_settings    <- 4                                         
es_lab <- c("Obion", "Hatchie_Wolf", "Yazoo", "BigBlack")  # labels for ecological settings

### populations/HUC10s 
populations <- c(4, 10, 3, 4) # Populations within ecological settings. Hierarchy needs to align with Excel 

### table that aligns loop over excel columns with ecological settings
pop_indexing <- data.frame("e_settings" = 1:e_settings, # "shift" aligns loop with Excel columns
                          "populations" = populations, 
                          "shift" = c(0, cumsum(populations)[1:(length(populations)-1)])) 

#### storage structures
### each vector corresponds to an ecological setting: Need to automate this
pop_resiliency_l<-lapply(populations,function(x){vector("list", length = x)})

es_resiliency_l <- vector("list", length = e_settings)  # resiliency of ecological settings
es_redundancy_l <- vector("list", length = e_settings)  # redundancy of ecological settings
es_extinction_l <- vector("list", length = e_settings)  # extirpation risk of ecological settings


############################# Population resiliency ############################
### priors of parent nodes
parent_priors  <- read.xlsx("bbn_inputs_pmt_cr.xlsx", sheet = "parent_priors")

### Observed states of parent nodes
parent_nodes   <- read.xlsx("bbn_inputs_pmt_cr.xlsx", sheet = "observed_resiliency")

### cpts of child nodes
threats_cpt     <- read.xlsx("bbn_inputs_pmt_cr.xlsx", sheet = "cpt_other_threats")
p_structure_cpt <- read.xlsx("bbn_inputs_pmt_cr.xlsx", sheet = "cpt_population_structure")
distribution_cpt<- read.xlsx("bbn_inputs_pmt_cr.xlsx", sheet = "cpt_local_distribution")
resiliency_cpt  <- read.xlsx("bbn_inputs_pmt_cr.xlsx", sheet = "cpt_pop_resiliency")

# i <- 1 # Used for debugging
# j <- 3 # Used for debugging
for(j in 1:e_settings){                   # outer loop over ecologcial settings
for(i in 1:populations[j]){               # inner loop over HUC10s

##### Child node: Other Threats
focal_cpt <- threats_cpt

### labels from nodes excel tab
parent_labs <- c("Nonnative Species (count)", "Hybridization") # these must match cpt order

### grab possible states from parent nodes: will need automation if expanding threats
parents <- list(                                             # 1, 2 correspond with parent nodes
  parent_priors[parent_priors$node == parent_labs[1], 1:2],  
  parent_priors[parent_priors$node == parent_labs[2], 1:2])
names(parents) <- parent_labs

### Index of observed state combinations: will need automation if expanding threats
s1 <- parent_nodes[parent_nodes$node == parent_labs[1] & parent_nodes[,2 + pop_indexing[j, 3] + i] == 1, 2] # indexed to excel columns
s2 <- parent_nodes[parent_nodes$node == parent_labs[2] & parent_nodes[,2 + pop_indexing[j, 3] + i] == 1, 2]

### indexing of observed states: will need automation if making threats flexible
observed_states <- ifelse(focal_cpt[, 1] == s1 &
                          focal_cpt[, 2] == s2, 1, 0) 

### create marginal probability from parent nodes:
tmp <- lapply(parents, enumerate_states)  # numerated states of parent nodes
tmp <- data.frame(setorder(expand.grid(tmp)))  # sort all combinations

### replace state with its prior probability sorted by tmp
parent_probs <- tmp          # copy dimensions of tmp
parent_probs[,] <- NA        # for quality control
for(k in 1:length(parent_probs)){                     
  parent_probs[, k] <- parents[[k]][tmp[, k],]$prob
}

## checks
#sum(apply(parent_probs,1,prod))                       # sums to 1, probability of each combination of states occurring
#length(apply(parent_probs,1,prod)) == nrow(focal_cpt) # should be true

### marginal probability child = t(conditional probability parent * observed state) * marginal probability of parent node
child_mp <- t((focal_cpt[,-c(1:length(parents))])*observed_states)%*%apply(parent_probs,1,prod)
child_mp <- child_mp/sum(child_mp)                                                       # put on probability scale
child_mp <- data.frame(state = names(focal_cpt[,-c(1:length(parents))]),prob = child_mp) # clean up to mesh with child nodes
threats_mp <- child_mp  # marginal probability table


##### Child node: Population structure
focal_cpt <- p_structure_cpt

parent_labs <- c("Naive Occupancy (Current)", "Qualitative Abundance (count)",
                 "Time Since Last Encounter (years)", "Naive Occupancy Trend (%)")

occupancy_p   <- parent_priors[parent_priors$node == parent_labs[1], 1:2]
abundance_p   <- parent_priors[parent_priors$node == parent_labs[2], 1:2]
lastcaught_p  <- parent_priors[parent_priors$node == parent_labs[3], 1:2]
trend_p       <- parent_priors[parent_priors$node == parent_labs[4], 1:2]

parents <- list(occupancy_p, abundance_p, lastcaught_p, trend_p) 

### indexing of observed state combinations
parent_labs <- c("Naive Occupancy (Current)", "Qualitative Abundance (count)",
                 "Time Since Last Encounter (years)", "Naive Occupancy Trend (%)")
s1 <- parent_nodes[parent_nodes$node == parent_labs[1] & parent_nodes[, 2 + pop_indexing[j, 3] + i] == 1, 2]
s2 <- parent_nodes[parent_nodes$node == parent_labs[2] & parent_nodes[, 2 + pop_indexing[j, 3] + i] == 1, 2]
s3 <- parent_nodes[parent_nodes$node == parent_labs[3] & parent_nodes[, 2 + pop_indexing[j, 3] + i] == 1, 2]
s4 <- parent_nodes[parent_nodes$node == parent_labs[4] & parent_nodes[, 2 + pop_indexing[j, 3] + i] == 1, 2]

### indexing for observed states
observed_states <- ifelse(focal_cpt[, 1] %in% s1 & 
                          focal_cpt[, 2] %in% s2 & 
                          focal_cpt[, 3] %in% s3 &
                          focal_cpt[, 4] %in% s4, 1, 0)  

### create marginal probability from two parent nodes
tmp <- lapply(parents, enumerate_states) 
tmp <- data.frame(setorder(expand.grid(tmp)))

parent_probs <- tmp 
parent_probs[,] <- NA
for(k in 1:length(parent_probs)){                     
  parent_probs[, k] <- parents[[k]][tmp[,k],]$prob
}

### Checks
#sum(apply(parent_probs,1,prod))                         # sums to 1, probability of each combination of states occurring
#length(apply(parent_probs,1,prod)) == nrow(focal_cpt)   # Should be true

### marginal probability
child_mp <- t((focal_cpt[,-c(1:length(parents))])*observed_states)%*%apply(parent_probs,1,prod)
child_mp <- child_mp/sum(child_mp)             # Put on probability scale
child_mp <- data.frame(state = names(focal_cpt[,-c(1:length(parents))]),prob = child_mp) # clean up to mesh with child nodes
p_structure_mp <- child_mp


##### child node: Local distribution
focal_cpt <- distribution_cpt

parent_labs <- c("Occupied Stream Length (km)", "Occupied Segments", "Network Complexity")

parents <- list(
  parent_priors[parent_priors$node == parent_labs[1], 1:2],
  parent_priors[parent_priors$node == parent_labs[2], 1:2],
  parent_priors[parent_priors$node == parent_labs[3], 1:2])
names(parents) <- parent_labs

### Index of observed state combinations
s1 <- parent_nodes[parent_nodes$node == parent_labs[1] & parent_nodes[, 2 + pop_indexing[j, 3] + i] == 1, 2]
s2 <- parent_nodes[parent_nodes$node == parent_labs[2] & parent_nodes[, 2 + pop_indexing[j, 3] + i] == 1, 2]
s3 <- parent_nodes[parent_nodes$node == parent_labs[3] & parent_nodes[, 2 + pop_indexing[j, 3] + i] == 1, 2]

### Indexing of observed states
observed_states <- ifelse(focal_cpt[, 1] == s1 & 
                          focal_cpt[, 2] == s2 & 
                          focal_cpt[, 3] == s3, 1, 0)  

### create marginal probability from two parent nodes
tmp <- lapply(parents, enumerate_states) 
tmp <- data.frame(setorder(expand.grid(tmp)))

parent_probs <- tmp 
parent_probs[,] <- NA
for(k in 1:length(parent_probs)){              
  parent_probs[, k] <- parents[[k]][tmp[,k],]$prob
}

## checks
# sum(apply(parent_probs,1,prod))                       # sums to 1, probability of each combination of states occurring
# length(apply(parent_probs,1,prod)) == nrow(focal_cpt) # Should be true

### marginal probability
child_mp <- t((focal_cpt[,-c(1:length(parents))])*observed_states)%*%apply(parent_probs,1,prod)
child_mp <- child_mp/sum(child_mp)             # Put on probability scale
child_mp <- data.frame(state = names(focal_cpt[,-c(1:length(parents))]),prob = child_mp) # clean up to mesh with child nodes
distribution_mp <- child_mp


##### child node: Population resiliency
focal_cpt <- resiliency_cpt

parents <- list(distribution_mp, p_structure_mp, threats_mp)

### create marginal probability from two parent nodes
tmp <- lapply(parents, enumerate_states) 
tmp <- data.frame(setorder(expand.grid(tmp)))

parent_probs <- tmp 
parent_probs[,] <- NA
for(k in 1:length(parent_probs)){                     
  parent_probs[, k] <- parents[[k]][tmp[,k],]$prob
}

## checks
sum(apply(parent_probs,1,prod))                       # sums to 1, probability of each combination of states occurring
length(apply(parent_probs,1,prod)) == nrow(focal_cpt) # Should be true

### marginal probability
child_mp <- t(focal_cpt[,-c(1:length(parents))])%*%apply(parent_probs,1,prod)
child_mp <- child_mp/sum(child_mp)             # Put on probability scale
child_mp <- data.frame(state = names(focal_cpt[,-c(1:length(parents))]),prob = child_mp) # clean up to mesh with child nodes
pop_resiliency_mp <- child_mp
pop_resiliency_l[[j]][[i]]          <- pop_resiliency_mp
names(pop_resiliency_l[[j]])[i] <- names(parent_nodes)[2 + pop_indexing[j, 3] + i]
names(pop_resiliency_l)[j] <- es_lab[j]
}    # end inner loop over HUC10s
# } # end outer loop if need be over ecological settings


######################### Mean resiliency: ecological setting ##################
#### create conditional probability table for ES resiliency on the fly
es_resiliency_states <- c("Secure", "At_Risk")                                         
states_pop_resilience <- data.frame(matrix(nrow = 2, ncol = populations[j], 
                            rep(es_resiliency_states, times = populations[j])))
colnames(states_pop_resilience) <- sapply(1:populations[j], function(x)paste("Resiliency", x, sep = ""))

es_resiliency_cpt <- setDT(expand.grid(states_pop_resilience))  # state combinations for each population
es_resiliency_cpt <- es_resiliency_cpt[, lapply(.SD, factor, levels = es_resiliency_states)]
es_resiliency_cpt <- data.frame(setorder(es_resiliency_cpt))
es_resiliency_cpt$Secure <- round((apply(es_resiliency_cpt, 1,                         # Could replicate this function for each state
                                         function(x)sum((x == "Secure")))/populations[j]), digits = 2)
es_resiliency_cpt$At_Risk <- 1 - es_resiliency_cpt$Secure

### focal cpt
focal_cpt <- es_resiliency_cpt

### indexing of observed states
observed_states <- rep(1, length = nrow(focal_cpt))

#### create marginal probability from two parent nodes
#tmp <- lapply(pop_resiliency_l[[j]],function(x){1:nrow(x)})
tmp <- lapply(pop_resiliency_l[[j]], enumerate_states) 
tmp <- data.frame(setorder(expand.grid(tmp)))

parent_probs <- tmp 
parent_probs[,] <- NA
parents_populations <- pop_resiliency_l[[j]]      
for(k in 1:length(parent_probs)){                                      
  parent_probs[, k]  <- parents_populations[[k]][tmp[,k],]$prob
}

## checks
#sum(apply(parent_probs,1,prod))                       # sums to 1, probability of each combination of states occurring
#length(apply(parent_probs,1,prod)) == nrow(focal_cpt) # Should be true


### marginal probability
child_mp <- t((focal_cpt[,-c(1:length(pop_resiliency_l[[j]]))]))%*%apply(parent_probs,1,prod)
child_mp <- child_mp/sum(child_mp)             # Put on probability scale
child_mp <- data.frame(state = names(focal_cpt[,-c(1:length(pop_resiliency_l[[j]]))]),prob = child_mp) # clean up to mesh with child nodes
es_resiliency_mp <- child_mp
es_resiliency_l[[j]] <- es_resiliency_mp
names(es_resiliency_l)[j] <- es_lab[j]


############################# Population Redundancy ############################
##### cpts of child nodes
pop_connectivity_cpt <- read.xlsx("bbn_inputs_pmt_cr.xlsx", sheet = "cpt_pop_connectivity")
es_redundancy_cpt    <- read.xlsx("bbn_inputs_pmt_cr.xlsx", sheet = "cpt_ES_redundancy")

##### observed states
node_redundancy <- read.xlsx("bbn_inputs_pmt_cr.xlsx", sheet = "observed_redundancy")

##### population connectivity
focal_cpt <- pop_connectivity_cpt

parent_labs <- c("Network Connectivity (km)", "Population Isolation (km)", "Ranging Movements (km)")

parents <- list(
  parent_priors[parent_priors$node == parent_labs[1], 1:2],
  parent_priors[parent_priors$node == parent_labs[2], 1:2],
  parent_priors[parent_priors$node == parent_labs[3], 1:2])
names(parents) <- parent_labs

### index of observed state combinations: index this
s1 <- node_redundancy[node_redundancy$node == parent_labs[1] & node_redundancy[, 2 + j] == 1, 2]
s2 <- node_redundancy[node_redundancy$node == parent_labs[2] & node_redundancy[, 2 + j] == 1, 2]
s3 <- node_redundancy[node_redundancy$node == parent_labs[3] & node_redundancy[, 2 + j] == 1, 2]

### indexing of observed states
observed_states <- ifelse(focal_cpt[, 1] %in% s1 & 
                          focal_cpt[, 2] %in% s2 & 
                          focal_cpt[, 3] %in% s3, 1, 0)  

### create marginal probability from two parent nodes
tmp <- lapply(parents, enumerate_states) 
tmp <- data.frame(setorder(expand.grid(tmp)))

parent_probs <- tmp 
parent_probs[,] <- NA
for(k in 1:length(parent_probs)){
  parent_probs[, k] <- parents[[k]][tmp[,k],]$prob
}

## checks
#sum(apply(parent_probs,1,prod))                       # sums to 1, probability of each combination of states occurring
#length(apply(parent_probs,1,prod)) == nrow(focal_cpt) # Should be true

### marginal probability
child_mp <- t((focal_cpt[,-c(1:length(parents))])*observed_states)%*%apply(parent_probs,1,prod)
child_mp <- child_mp/sum(child_mp)             # Put on probability scale
child_mp <- data.frame(state = names(focal_cpt[,-c(1:length(parents))]),prob = child_mp) # clean up to mesh with child nodes
pop_connectivity_mp <- child_mp


##### child node: Redundancy
focal_cpt <- es_redundancy_cpt

parent_labs <- c("Proportion Extant (count)", "Extant Populations", "pop_connectivity_mp")

parents <- list(parent_priors[parent_priors$node == parent_labs[1], 1:2],
                parent_priors[parent_priors$node == parent_labs[2], 1:2],
                pop_connectivity_mp)

### indicator of observed state combinations
s1 <- node_redundancy[node_redundancy$node == parent_labs[1] & node_redundancy[, 2 + j] == 1, 2]
s2 <- node_redundancy[node_redundancy$node == parent_labs[2] & node_redundancy[, 2 + j] == 1, 2]
s3 <- pop_connectivity_mp[, 1]  # Need to propogate marginal probabilities of states

## indexing of observed states
observed_states <- ifelse(focal_cpt[, 1] %in% s1 & 
                          focal_cpt[, 2] %in% s2 & 
                          focal_cpt[, 3] %in% s3, 1, 0)  

### create marginal probability from two parent nodes
tmp <- lapply(parents, enumerate_states)
tmp <- data.frame(setorder(expand.grid(tmp)))

parent_probs <- tmp 
parent_probs[,] <- NA
for(k in 1:length(parent_probs)){
  parent_probs[, k] <- parents[[k]][tmp[,k],]$prob
}

## checks
#sum(apply(parent_probs,1,prod))                       # sums to 1, probability of each combination of states occurring
#length(apply(parent_probs,1,prod)) == nrow(focal_cpt) # Should be true

### marginal probability
child_mp <- t((focal_cpt[,-c(1:length(parents))])*observed_states)%*%apply(parent_probs,1,prod)
child_mp <- child_mp/sum(child_mp)             # Put on probability scale
child_mp <- data.frame(state = names(focal_cpt[,-c(1:length(parents))]),prob = child_mp) # clean up to mesh with child nodes
es_redundancy_mp <- child_mp
es_redundancy_l[[j]] <- es_redundancy_mp
names(es_redundancy_l)[j] <- es_lab[j]


############################# Taxon Vulnerability ############################
##### child node: Other Threats
##### cpts of child nodes
specialization_cpt <- read.xlsx("bbn_inputs_pmt_cr.xlsx", sheet = "cpt_specialization")

##### observed states
node_vulnerability <- read.xlsx("bbn_inputs_pmt_cr.xlsx", sheet = "observed_vulnerability")   # houses observed states

##### Specialization
focal_cpt <- specialization_cpt

parent_labs <- c("Adult Feeding Guild", "Benthic Dependency",                 # Needs to match parent_priors tab 
                 "Drift Dependency", "Lotic Dependency")

parents <- list(
  parent_priors[parent_priors$node == parent_labs[1], 1:2],
  parent_priors[parent_priors$node == parent_labs[2], 1:2],
  parent_priors[parent_priors$node == parent_labs[3], 1:2],
  parent_priors[parent_priors$node == parent_labs[4], 1:2])
names(parents) <- parent_labs

### index of observed state combinations
s1 <- node_vulnerability[node_vulnerability$node == parent_labs[1] & node_vulnerability[,3] == 1, 2]
s2 <- node_vulnerability[node_vulnerability$node == parent_labs[2] & node_vulnerability[,3] == 1, 2]
s3 <- node_vulnerability[node_vulnerability$node == parent_labs[3] & node_vulnerability[,3] == 1, 2]
s4 <- node_vulnerability[node_vulnerability$node == parent_labs[4] & node_vulnerability[,3] == 1, 2]

### indexing of observed states
observed_states <- ifelse(focal_cpt[, 1] %in% s1 & 
                            focal_cpt[, 2] %in% s2 & 
                            focal_cpt[, 3] %in% s3 &
                            focal_cpt[, 4] %in% s4, 1, 0)  

#### create marginal probability from two parent nodes
tmp <- lapply(parents, enumerate_states)
tmp <- data.frame(setorder(expand.grid(tmp)))

parent_probs <- tmp 
parent_probs[,] <- NA
for(k in 1:length(parent_probs)){                    
  parent_probs[, k] <- parents[[k]][tmp[,k],]$prob
}

## checks
# sum(apply(parent_probs,1,prod))                       # sums to 1, probability of each combination of states occurring
# length(apply(parent_probs,1,prod)) == nrow(focal_cpt) # Should be true

### marginal probability
child_mp <- t((focal_cpt[,-c(1:length(parents))])*observed_states)%*%apply(parent_probs,1,prod)
child_mp <- child_mp/sum(child_mp)             # Put on probability scale
child_mp <- data.frame(state = names(focal_cpt[,-c(1:length(parents))]),prob = child_mp) # clean up to mesh with child nodes
specialization_mp <- child_mp


##### child node: vulnerability
vulnerability_cpt <- read.xlsx("bbn_inputs_pmt_cr.xlsx", sheet = "cpt_taxon_vulnerability")

### observed states
node_vulnerability <- read.xlsx("bbn_inputs_pmt_cr.xlsx", sheet = "observed_vulnerability")   # houses observed states

focal_cpt <- vulnerability_cpt

parent_labs <- c("Life History Strategy", "Maximum Length (mm)", "specialization_mp")

parents <- list(parent_priors[parent_priors$node == parent_labs[1], 1:2],
                parent_priors[parent_priors$node == parent_labs[2], 1:2],
                specialization_mp)

### indicator of observed state combinations
s1 <- node_vulnerability[node_vulnerability$node == parent_labs[1] & node_vulnerability[, 3] == 1, 2]
s2 <- specialization_mp[, 1]  # Need to propogate marginal probabilities of states
s3 <- node_vulnerability[node_vulnerability$node == parent_labs[2] & node_vulnerability[, 3] == 1, 2]


### indexing of observed states: because focal_cpt columns arent explicit, order of columns in Excel matters
observed_states <- ifelse(focal_cpt[, 1] %in% s1 & 
                            focal_cpt[, 2] %in% s2 & 
                            focal_cpt[, 3] %in% s3, 1, 0)  

### create marginal probability from two parent nodes
tmp <- lapply(parents, enumerate_states)
tmp <- data.frame(setorder(expand.grid(tmp)))

parent_probs <- tmp 
parent_probs[,] <- NA
for(k in 1:length(parent_probs)){
  parent_probs[, k] <- parents[[k]][tmp[,k],]$prob
}

## checks
# sum(apply(parent_probs,1,prod))                       # sums to 1, probability of each combination of states occurring
# length(apply(parent_probs,1,prod)) == nrow(focal_cpt) # Should be true

### marginal probability
child_mp <- t((focal_cpt[,-c(1:length(parents))])*observed_states)%*%apply(parent_probs,1,prod)
child_mp <- child_mp/sum(child_mp)             # Put on probability scale
child_mp <- data.frame(state = names(focal_cpt[,-c(1:length(parents))]),prob = child_mp) # clean up to mesh with child nodes
vulnerability_mp <- child_mp


############################## Extinction risk #################################
##### cpts of child nodes
es_extinction_cpt <- read.xlsx("bbn_inputs_pmt_cr.xlsx", sheet = "cpt_ES_extinction_risk")

##### taxon vulnerability
focal_cpt <- es_extinction_cpt

parents <- list(es_resiliency_l[[j]], es_redundancy_l[[j]], vulnerability_mp)

### create marginal probability from two parent nodes
tmp <- lapply(parents, enumerate_states)
tmp <- data.frame(setorder(expand.grid(tmp)))

parent_probs <- tmp 
parent_probs[,] <- NA
for(k in 1:length(parent_probs)){                   
  parent_probs[, k] <- parents[[k]][tmp[,k],]$prob
}

## checks
# sum(apply(parent_probs,1,prod))                       # sums to 1, probability of each combination of states occurring
# length(apply(parent_probs,1,prod)) == nrow(focal_cpt) # Should be true

### marginal probability
child_mp <- t((focal_cpt[,-c(1:length(parents))]))%*%apply(parent_probs,1,prod)
child_mp <- child_mp/sum(child_mp)                                                       # Put on probability scale
child_mp <- data.frame(state = names(focal_cpt[,-c(1:length(parents))]),prob = child_mp) # clean up to mesh with child nodes
es_extinction_mp <- child_mp
es_extinction_l[[j]] <- es_extinction_mp
names(es_extinction_l)[j] <- es_lab[j]
}                                                                                        # end outer loop over ecological settings


########################## Global at-risk of extinction ########################
#### create conditional probability table for ES extinction on the fly
g_extirpation_states <- c("Secure", "At_Risk")
g_extirpation_parents <- data.frame(matrix(nrow = 2, ncol = e_settings, 
                            rep(g_extirpation_states, times = e_settings)))
colnames(g_extirpation_parents) <- sapply(1:e_settings, function(x)paste("extirpation", x, sep = ""))

g_extinction_cpt <- setDT(expand.grid(g_extirpation_parents))  # state combinations for each population
g_extinction_cpt <- g_extinction_cpt[, lapply(.SD, factor, levels = g_extirpation_states)]
g_extinction_cpt <- data.frame(setorder(g_extinction_cpt))
g_extinction_cpt$Secure <- round((apply(g_extinction_cpt, 1, 
                                         function(x)sum((x == "Secure")))/e_settings), digits = 2)
g_extinction_cpt$At_Risk <- 1 - g_extinction_cpt$Secure

### focal cpt
focal_cpt <- g_extinction_cpt  # global extinction risk cpt

### indexing of observed states
observed_states <- rep(1, length = nrow(focal_cpt)) # incorporate all probabilities

### create marginal probability from two parent nodes
tmp <- lapply(es_extinction_l, enumerate_states)
tmp <- data.frame(setorder(expand.grid(tmp)))

parent_probs <- tmp 
parent_probs[,] <- NA
for(k in 1:length(parent_probs)){                                       
  parent_probs[, k]  <- es_extinction_l[[k]][tmp[,k],]$prob
}

## checks
#sum(apply(parent_probs,1,prod))                       # should sum to 1
#length(apply(parent_probs,1,prod)) == nrow(focal_cpt) # Should be true

### marginal probability
child_mp <- t((focal_cpt[,-c(1:length(es_extinction_l))]))%*%apply(parent_probs,1,prod)
child_mp <- child_mp/sum(child_mp)
child_mp <- data.frame(state = names(focal_cpt[,-c(1:length(es_extinction_l))]),prob = child_mp)
global_extinction <- child_mp


################# output results as Excel spreadsheet ##########################
### Add 
HUC_resiliency_DF <- data.frame("HUC" = rep(NA, times = sum(populations)), 
                                "Secure" = rep(NA, times = sum(populations)),
                                "At_risk" = rep(NA, times = sum(populations)))

### clean up population level resiliency
for(i in 1:length(pop_resiliency_l)){
  for(j in 1:populations[i]){
    HUC_resiliency_DF[pop_indexing[i, 3] + j,"HUC"] <- strsplit(names(pop_resiliency_l[[i]][j]), "_")[[1]][2]
    HUC_resiliency_DF[pop_indexing[i, 3] + j,"Secure"] <- data.frame(pop_resiliency_l[[i]][j])[1,2]
    HUC_resiliency_DF[pop_indexing[i, 3] + j,"At_risk"] <- data.frame(pop_resiliency_l[[i]][j])[2,2]
  }
}

pop_resiliency_df <- 
  data.frame("ecological_setting" = c(rep("obion", times = populations[1]), 
  rep("hatchie_wolf", times = populations[2]),
  rep("yazoo", times = populations[3]),
  rep("big_black", times = populations[4])), HUC_resiliency_DF)

### es resiliency
es_resiliency_df <- data.frame(attributes(es_resiliency_l),
                               as.data.frame(do.call("rbind", lapply(es_resiliency_l, '[[', 2))))
names(es_resiliency_df) <- c("ecological_setting", "p_secure", "p_at_risk")
rownames(es_resiliency_df) <- NULL

### es redundancy
es_redundancy_df <- data.frame(attributes(es_redundancy_l),
                               as.data.frame(do.call("rbind", lapply(es_redundancy_l, '[[', 2))))
names(es_redundancy_df) <- c("ecological_setting", "p_adequate", "p_inadequate")
rownames(es_redundancy_df) <- NULL

### es extirpation risk
es_extirpation_df <- data.frame(attributes(es_extinction_l),
                               as.data.frame(do.call("rbind", lapply(es_extinction_l, '[[', 2))))
names(es_extirpation_df) <- c("ecological_setting", "p_secure", "p_at_risk")
rownames(es_extirpation_df) <- NULL

#### write to excel spreadsheet for easy review
piebald_3Rs_ms <- createWorkbook("piebaldmadtom_results")

addWorksheet(piebald_3Rs_ms, "pop_resiliency")
addWorksheet(piebald_3Rs_ms, "es_resiliency")
addWorksheet(piebald_3Rs_ms, "es_redundancy")
addWorksheet(piebald_3Rs_ms, "es_extirpation_risk")
addWorksheet(piebald_3Rs_ms, "vulnerability")
addWorksheet(piebald_3Rs_ms, "global_extinction_risk")


writeData(piebald_3Rs_ms, sheet = "pop_resiliency", pop_resiliency_df)
writeData(piebald_3Rs_ms, sheet = "es_resiliency", es_resiliency_df)
writeData(piebald_3Rs_ms, sheet = "es_redundancy", es_redundancy_df)
writeData(piebald_3Rs_ms, sheet = "es_extirpation_risk", es_extirpation_df)
writeData(piebald_3Rs_ms, sheet = "vulnerability", vulnerability_mp)
writeData(piebald_3Rs_ms, sheet = "global_extinction_risk", global_extinction)

saveWorkbook(piebald_3Rs_ms, "piebaldmadtom_results_cr.xlsx", overwrite = T)


########################## print outputs #######################################
vulnerability_mp  # vulnerability
pop_resiliency_df # resiliency at HUC10
es_resiliency_df  # resiliency within Ecological Setting
es_redundancy_df  # redundancy within Ecological Setting
es_extirpation_df # extirpation risk by Ecological Setting
global_extinction # global extinction risk

