library(dplyr)
library(tidyr)

rm(list = ls())

### Control functions

# Individuals are stored in a data frame, containing age and alleles

# Function for calculating a vital rate
## CURRENTLY MISTAKE IS HERE
calc.vr = function(alls, rate, mal.rate) c(mal.rate, rate, rate)[apply(alls, 1, sum) + 1] 
# alls - alleles (coded as 0 or 1)
# rate - vital rate value for dominant phenotype (higher vital rate)
# mal.rate - vital rate for recessive phenotype (lower vital rate, `mal`)
# sum of alleles of 0 gives recessive, sum 1 or 2 gives dominant

# Function intakes data frame of individuals's alleles, outputs their vital rates
# (wrapper for use of calc.vr from above)
ins.vrs = function(alls, params) {
  suppressMessages(attach(params))
  df = alls %>%
    mutate(s1 = calc.vr(alls %>% select(s1.a, s1.b), s1, mal.s1), # see note for calc.vr
           s2 = calc.vr(alls %>% select(s2.a, s2.b), s2, mal.s2),
           phi = calc.vr(alls %>% select(phi.a, phi.b), phi, mal.phi))
  detach(params)
  return(df)
}
# alls is data frame of individuals and their alleles

# Function for initializing a population
init.popn = function(params) {
  # Returns population data frame (info below)

  suppressMessages(attach(params))
  # Attaches the following variables to namespace:
  
  # # Initial population size
  # n.init = params$n.init
  # # Frequency of each wild type allele
  # p.s1 = params$p.s1
  # p.s2 = params$p.s2
  # p.phi = params$p.phi
  # # Initial stage distribution (found from SSD of matrix model)
  # init.stage = params$init.stage
  
  # Data frame containing the population
  # At this stage, contains the following:
  #  i - individual identifier
  #  t - time step (population will be stored in long form)
  #  age - 0 or 1 for juvenile or adult
  #  s1.a, s1.b - two copies of allele for s1; 1 is wild type, 0 is mutant
  #  s2.a, s2.b, phi.a, phi.b - same coding as s1 locus
  
  popn = data.frame(i = 1:n.init,
                    t = 1,
                    age = sample(0:1, size = n.init, 
                                 replace = TRUE, 
                                 prob = c(init.stage1, init.stage2)),
                    s1.a = sample(rep(0:1, times = ceiling(n.init * c(1 - p.s1, p.s1))),
                                  size = n.init),
                    s1.b = sample(rep(0:1, times = ceiling(n.init * c(1 - p.s1, p.s1))),
                                  size = n.init),
                    s2.a = sample(rep(0:1, times = ceiling(n.init * c(1 - p.s2, p.s2))),
                                  size = n.init),
                    s2.b = sample(rep(0:1, times = ceiling(n.init * c(1 - p.s2, p.s2))),
                                  size = n.init),
                    phi.a = sample(rep(0:1, times = ceiling(n.init * c(1 - p.phi, p.phi))),
                                   size = n.init),
                    phi.b = sample(rep(0:1, times = ceiling(n.init * c(1 - p.phi, p.phi))),
                                   size = n.init))
  
  detach(params)
  
  return(popn)
}

sample.alleles = function(df, parents, rate.string) {
  
  alls = df[parents, grep(paste0(rate.string, '\\.'), names(df))] %>% unlist()
  
  return(alls[1:length(parents) + length(parents) * sample(0:1, length(parents), replace = TRUE)])

}
# A function to efficiently segregate alleles to offspring
# df is a data frame of previous generation (i.e., potential parents)
# parents is an array with rows to subset (i.e., parents for each offspring)
# rate.string is a string of 's1', 's2', 'phi' used for subsetting columns out of df

# Function to produce offspring
reproduce = function(popn, params) {

  # Use poisson draws to determine number of offspring for each parent
  # (rate of poisson is dependent on age; if age is 0, rate is 0 (no offspring),
  # if age is 1, then rate is phi)
  assign.parents = popn %>%
    ins.vrs(params = params) %>%
    mutate(off = rpois(nrow(.), lambda = ifelse(popn$age, phi, 0)))
  
  # Loop handles whether or not there are any offspring;
  # if there are no offspring, then return nothing
  if (sum(assign.parents$off)) {
    
    # Assign row indices for parents
    # Maternal parent for each offspring determined by `offspring` column
    moms = rep(1:nrow(assign.parents), times = assign.parents$off)
    # Paternal parent is randomly sampled from adults, with replacement
    # NOTE: this allows for chance of selfing, inversely proportional to population size
    dads = sample(which(assign.parents$age > 0), size = sum(assign.parents$off), replace = TRUE)
    
    # Generate offspring
    offspr = assign.parents %>%
      sample_n(size = sum(assign.parents$off)) %>%
      mutate_all(function(x) 0) %>%
      select(-c(s1, s2, phi, off)) %>%
      mutate(i = max(assign.parents$i) + 1:nrow(.),
             t = max(assign.parents$t) + 1,
             age = 0,
             s1.a = sample.alleles(assign.parents, moms, 's1'),
             s1.b = sample.alleles(assign.parents, dads, 's1'),
             s2.a = sample.alleles(assign.parents, moms, 's2'),
             s2.b = sample.alleles(assign.parents, dads, 's2'),
             phi.a = sample.alleles(assign.parents, moms, 'phi'),
             phi.b = sample.alleles(assign.parents, dads, 'phi'))
    # sample_n() and mutate_all() is used to neatly copy and clear all columns
    
    return(offspr)
  } else {
    
    return(NULL)
  
  }
  
  
}

# Takes in a population, assigns survvial (based on age), increments age and time
survival = function(popn) {
  popn %>%
    mutate(surv = rbinom(n = nrow(.), size = 1, prob = ifelse(age, s2, s1))) %>%
    filter(surv > 0) %>% select(-surv) %>%
    mutate(t = t + 1,
           age = 1)
}

### Structure of simulation

### First, initialize parameters and population

## Initialize parameters

params = list(n.init = 99, # initial population size
              s1 = 0.5, # survival of juveniles
              s2 = 0.6, # survival of adults
              phi = 1.0, # fecundity of adults
              p.s1 = 10/11, # proportion of wild-type s1 allele
              p.s2 = 10/11, # proportion of wild-type s2 allele
              p.phi = 10/11, # proportion of wild-type phi allele
              mu.mal = 0.95, # mean fitness for individual with ONE maladapted allele
              end.time = 10) # length of simulation
params = c(params, 
           # initial stage distribution
           init.stage = with(params, 
                             matrix(c(0, phi, s1, s2), byrow = TRUE, nrow = 2) %>%
                               eigen() %>%
                               (function(x) x$vectors[,1]) %>%
                               (function(x) x / sum(x))),
           # vital rates for recessive individuals
           mal.s1 = with(params, (mu.mal^2 - s2 * mu.mal) / phi),
           mal.s2 = with(params, mu.mal - s1*phi/mu.mal),
           mal.phi = with(params, (mu.mal^2 - s2 * mu.mal) / s1))

# These values were chosen to produce a mean deterministic growth rate (mu) of ~1.08
# However, there will be variance in these parameters;
# survival (for both age classes) is a binomial random variable
# recruits per adult is Poisson distributed
# (not sure yet how to calculate the variance around the growth rate...)

# Deterministic growth rate (mu)
with(params, 
     matrix(c(0, phi, s1, s2), byrow = TRUE, nrow = 2) %>%
       eigen() %>%
       (function(x) x$values[1]))

# Seed
set.seed(1009)

### Initialize population

### Detach is doing some funky stuff in here!!###
### Suppressing warnings to make this neater  ###
all.gens = init.popn(params = params)
# t = 1

# Create variable for previous generation
prev.gen = all.gens

for (t in 1:params$end.time) {
  
  ### Next, run through generations
  ### For generation in generations 1:t
  
  ## If there are individuals in the next generation,
  ## Produce offspring
  offspring = reproduce(prev.gen, params = params)
  
  
  ## Next, determine adult survival
  
  survivors = survival(prev.gen)
  # Age = 1, t = cur.gen
  
  all.gens = rbind(all.gens, offspring, survivors)
  prev.gen = rbind(offspring, survivors)
  
}

table(all.gens$t) %>% plot(type = 'l')
