# Demography and inbred, fragmented populations
# Scott Nordstrom
# October 25, 2019

# Project setup

This project idea arose from the *Echinacea* demography project (LTRE). A large part of that project was the observation that inbred populations will have lowered juvenile survival. We treated this as a most-likely outcome for populations reduced in size due to fragmentation.

This assumes that the primary demographic consequence of inbreeding is the fixation of alleles affecting one vital rate. But, the genome is complex, and different vital rates (probably) are controlled by multiple genes. Presumably, each of these genes are subject to possible inbreeding depression, meaning that there can be multiple demographic consequences of inbreeding. For example, in a plant there can be genes which control seedling germination, survival, and fecundity. Each one of these genes can be independtly subject to inbreeding depression.

This project wants to explore what happens in populations where multiple non-linked genes are subject to inbreeding depression. Here are some basic questions I am seeking to answer:

* Will a population with multiple, independent (uncorrelated) loci subject to inbreeding depression always reach the same genetic (and therefore demographic) outcomes after fragmentation? I.e., if one locus reaches fixation for a maladapted allele, will it always be the same locus?

This is the first, most fundamental question I am interested in answering. The setup here is to have one large population with genetic variation at three loci. These loci each affect different parts of the life cycle and are dominant for a normal phenotype (with a recessive phenotype with a lower vital rate). For the sake of simplicity, say there are two stage classes and three vital rates: survival of juveniles, survival of adults, and fecundity of adults (juveniles are not reproductively mature). Fragment the population and see what the genetic makeup of each subpopulation is over time. Outcomes to watch include:
  * Do populations reach fixation of any individual allele
  * If yes, which locus reaches fixation? Is it the same locus in each case?
  * Do populations reach fixation of only one locus?
  
I think the modal outcome would be most sensitive to he relative degree of fitness penalty for each vital rate. If one vital rate is the *least* unfit, then it should be favored over the other vital rates. If each vital rate decreases fitness by the same amount, then we may expect each outcome to be equally likely (on the other hand, although we can make the mean fitness penalty for each vital rate the same, it's difficult to do this for how sensitive *the variance* in the growth rate is to each vital rate). Also, whether the fitness penalty for each/any vital rate pulls the population below a mean fitness of 1 may also influence the outcomes of the populations. This is perhaps something to vary in our simulations.

Sources of stochasticity which may pull the population away from the modal outcome: 
  1) sampling variation in the initial distribution of alleles within a population
  2) sampling variation in mating and pairing of alleles
  3) stochasticity in demographic processes (i.e., process variance associated with each vital rate)
  
Seems like if we took out each of of these, we could determine the "mean" outcome pretty easily (question: is this the same as the modal outcome?). We can remove 1) by fixing the number of non-wild-type alleles in each population (and also assigning these to each subpopulation in equal proportion). 2) I am sure there is an easy fix for although I'm not quite sure what it is; this seems like it would just be through modeling change in allele frequencies as a determininstic process). 3) can be done by making this into a non-individual based model.

