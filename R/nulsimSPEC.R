#' Simple Neutral community (SN) saving certain precise timepoints
#'
#' @param tmax Arbitrary units of time the simulation should be run
#' @param b1 Birthrate
#' @param d1 Deathrate
#' @param k1 Carrying capacity
#' @param abun_original A vector of initial abundances for each species such as those created by generate_spat_abund
#' @param interval After how many evenst should the population state be saved to the output matrix
#' @param wantedtimes The specific timepoints that should be included in the final matrix
#' @author Timo van Eldijk
#'
#' @return A matrix denoting the abundances of all the species in the community over time, with certain precise timepoints saved (wantedtimes)
#' @export
#'
#' @examples
#'nulsimSPEC(10, 0.6, 0.1, 16000, generate_spat_abund(theta = 200,Ivec = rep(40,1),
#'Jvec = c(16000)), 20000000, c(1,3,5))
#'


nulsimSPEC = function (tmax, b1, d1, k1, abun_original, interval, wantedtimes) {
  wantedtimescounter=1
  print(Sys.time())
  nspec = length (abun_original) #Determine number of species in input community
  b=c(rep(b1,nspec))            #Vector of birth rates (this is so function can be rewritten to accept vector)
  d=c(rep(d1,nspec))            #Vector of deathrates

  k=c(rep(k1,nspec))            #Vector of carrying capacities (idem so vector could be input)
  specnrs = 1:nspec             #Number the species
  pop = data.frame()
  pop = specnrs
  pop = cbind(pop, abun_original) #Put the population into  a data frame

  t = 0 #Set time to zero

  keeper = data.frame() #empty dataframe to store data
  keep = c(t, pop[, 2]) #store initial abundances
  keeper = rbind (keeper, keep)#put them in dataframe
  saver = 0       #saver to 0


  totpopimpact = sum(pop[, 2])/k #Determine how much of the carrying capacity is occupied, never be lower than 0!
  probo = (max(0,(b*(1-totpopimpact))) + d) * pop[, 2] #Get vector of event probabilities (1 per species)
  timetoevent = stats::rexp(1, sum(probo)) #Sample time to first event
  t=t+timetoevent

  while (t < tmax) {
    #totpopimpact = sum(pop[, 2])/k #How much capacity occupied
    #probo = (max(0,(b*(1-totpopimpact))) + d) * pop[, 2] #Vector of event probabilites (1 per species)
    species = sample (c(1:nspec),size = 1,replace = T,prob = (c(probo)))#Sample which species has event

    event = sample(c(-1, 1),size = 1,replace = T,prob = c(d[species], max(0,(b[species]*(1-totpopimpact[species]))))) #Sample if the event is birth or death


    if (event < 2) {
      pop[species, 2] = pop[species, 2] + event #implement birth or death (done this way to easily allow for mutation)
    }



    totpopimpact = sum(pop[, 2])/k #sample the time to the next event
    probo = (max(0,(b*(1-totpopimpact))) + d) * pop[, 2]
    timetoevent = stats::rexp(1, (sum(probo)))
    tminone=t
    t = t + timetoevent
    saver = saver + 1

    if (t>wantedtimes[wantedtimescounter] && wantedtimescounter<=length(wantedtimes)){

      keep=c(tminone,pop[,2]) #store abundances
      keeper=rbind (keeper, keep)
      print (c(tminone,t,"WANTED"))
      wantedtimescounter=wantedtimescounter+1

    }


  }

  if (tmax!=wantedtimes[length(wantedtimes)]){
    keep = c(tminone, pop[, 2]) #store abundances for last timepoint before 100 (so what it was exactly at 100)
    keeper = rbind (keeper, keep)#put them in dataframe
  }

  return(t(keeper))
}
