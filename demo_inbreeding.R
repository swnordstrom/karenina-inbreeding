library(dplyr)
library(tidyr)

s1 = 0.5
s2 = 0.6
phi = 1.05

m = matrix(c(0, phi, s1, s2), byrow = TRUE, nrow = 2)
eigen(m)

init.stage = eigen(m)$vectors[,1] %>% (function(x) x/sum(x))

# s.s1 = phi / sqrt(s2^2 + 4*s1*phi)
# s.s2 = 1/2 + (1 + s2 / sqrt(s2^2 + 4*phi*s1))
# s.phi = s1 / sqrt(s2^2 + 4*s1*phi)
# s.s1; s.s2; s.phi

l.mal = 0.95

mal.s1 = (l.mal^2 - s2 * l.mal) / phi
mal.s2 = l.mal - s1*phi/l.mal
mal.phi = (l.mal^2 - s2 * l.mal) / s1
mal.s1; mal.s2; mal.phi

## simlate popn

n.init = 99;

set.seed(1009)

popn = data.frame(i = 1:n.init,
                  t = 1,
                  age = sample(0:1, size = n.init, replace = TRUE, prob = init.stage),
                  s1.a = sample(rep(0:1, times = c(9, 90))),
                  s1.b = sample(rep(0:1, times = c(9, 90))),
                  s2.a = sample(rep(0:1, times = c(9, 90))),
                  s2.b = sample(rep(0:1, times = c(9, 90))),
                  phi.a = sample(rep(0:1, times = c(9, 90))),
                  phi.b = sample(rep(0:1, times = c(9, 90))))

calc.vr = function(alls, rate, mal.rate) c(mal.rate, rate, rate)[apply(alls, 1, sum) + 1] 

ins.vrs = function(alls) {
	alls %>%
		mutate(s1 = calc.vr(alls %>% select(s1.a, s1.b), s1, mal.s1),
		   	   s2 = calc.vr(alls %>% select(s2.a, s2.b), s2, mal.s2),
		       phi = calc.vr(alls %>% select(phi.a, phi.b), phi, mal.phi))
}

## Hmm... okay matrix model is asexual but allele stuff is sexual so

popn = popn %>%
  ins.vrs() %>%
  mutate(off = rpois(nrow(.), lambda = ifelse(popn$age, phi, 0)))

prev.gen = popn

if (sum(prev.gen$off)) {
  offspr = prev.gen %>%
    sample_n(size = sum(prev.gen$off), replace = TRUE) %>%
    mutate_all(function(x) 0) %>%
    select(-c(s1, s2, phi, off)) %>%
    mutate(i = max(prev.gen$i) + 1:nrow(.),
           t = max(prev.gen$t) + 1,
           age = 0)
  
  moms = rep(1:nrow(prev.gen), times = prev.gen$off)# prev.gen %>% rep(1:nrow(.), times = off)
  dads = sample(1:nrow(prev.gen), size = sum(prev.gen$off), replace = TRUE)# with(prev.gen, sample(i, size = sum(off), replace = TRUE))
}

sample.alleles = function(df, pars, rate.string) {
  alls = df[pars, grep(paste0(rate.string, '\\.'), names(df))] %>% unlist()
  return(alls[1:length(pars) + length(pars) * sample(0:1, size = length(pars), replace = TRUE)])
}

offspr = offspr %>%
  mutate(s1.a = sample.alleles(prev.gen, moms, 's1'),
         s1.b = sample.alleles(prev.gen, dads, 's1'),
         s2.a = sample.alleles(prev.gen, moms, 's2'),
         s2.b = sample.alleles(prev.gen, dads, 's2'),
         phi.a = sample.alleles(prev.gen, moms, 'phi'),
         phi.b = sample.alleles(prev.gen, dads, 'phi')) %>%
  ins.vrs()

surviv = prev.gen %>% 
  mutate(surv = rbinom(n = nrow(.), size = 1, prob = ifelse(age, s2, s1))) %>%
  filter(surv > 0) %>%
  select(-surv, off) %>%
  mutate(t = t + 1,
         age = 1)

all.gen = rbind(all.gen, surviv, offspr)