#' Simulate a neutral community with a variable birthrate (VBN), saving certain precise timepoints
#'
#' @param tmax Arbitrary units of time the simulation should be run
#' @param b1 Normal birthrate
#' @param d1 Deathrate
#' @param k1 Carrying capacity
#' @param bneck Birthrate during the bottleneck
#' @param kneckstart Timepoint at which the bottleneck starts
#' @param kneckend  Timepoint at which the bottleneck ends
#' @param abun_original A vector of initial abundances for each species such as those created by generate_spat_abund
#' @param interval After how many evenst should the population state be saved to the output matrix
#' @param wantedtimes The specific timepoints that should be included in the final matrix
#'@author Timo van Eldijk
#' @return A matrix denoting the abundances of all the species in the community over time, with certain specific timepoints saved (wantedtimes)
#' @export
#'
#' @examples
#' simbvarSPEC( 10, 0.6, 0.1, 16000, 0.05,0, 5,
#' generate_spat_abund(theta = 200,Ivec = rep(40,1),Jvec = c(16000)),
#' 200,c(15,30,50,75,100))
#'

simbvarSPEC = function ( tmax, b1, d1, k1, bneck,kneckstart, kneckend, abun_original, interval, wantedtimes){
  wantedtimescounter=1
  print(Sys.time())
  nspec=length(abun_original) #Determine number of species in input community
  bn=c(rep(b1,nspec)) #vector of birthrates
  d=c(rep(d1,nspec)) #vector of deathrates

  k=c(rep(k1,nspec)) #Vector of normal carrying capacity
  bneck=c(rep(bneck,nspec)) #Vector of increased carrying capacity
  specnrs=1:nspec #Number species

  pop=data.frame()
  pop=specnrs
  pop=cbind(pop, abun_original) #create population data frame


  t=0 #set time to zero

  keeper=data.frame() #data frame for saving the data
  keep=c(t,pop[,2]) #store initial abundances
  keeper=rbind (keeper, keep)#put them in dataframe
  saver=0#saver to 0



  if (t > kneckstart && t < kneckend) {b=bneck} else {b=bn} # Check if time is at the inteval specified for increased K and change K accordingly

  totpopimpact=sum(pop[,2])/k #How much of carrying capacity is occupied
  probo=(max(0,(b*(1-totpopimpact)))+d)*pop[,2] #vector of event probabilites (1 per species)
  timetoevent=stats::rexp(1, sum(probo)) #Sample time to next event
  t=t+timetoevent #Add time to next event to time

  while (t<tmax){ #Loop through

    if (t > kneckstart && t < kneckend) {b=bneck} else {b=bn}# Check if time is at the inteval specified for increased K and change K accordingly

    #totpopimpact=sum(pop[,2])/k #How much carrying capacity occupied
    #probo=(max(0,(b*(1-totpopimpact)))+d)*pop[,2] #vector of event probabilities (1 per species)
    species=sample (c(1:nspec), size=1, replace=T, prob=(c(probo))) #Sample which species has an event

    event=sample(c(-1,1), size=1, replace=T, prob = c(d[species], max(0,(b*(1-totpopimpact))))) #Sample the type of event (birth or death)
    if (event<2){pop[species, 2]=pop[species, 2]+event} #Process this event on the species that has the event happen.



    totpopimpact=sum(pop[,2])/k#how much carrying capacity occupied
    probo=(max(0,(b*(1-totpopimpact)))+d)*pop[,2] #vector of event probabilites (1 per species)
    timetoevent=stats::rexp(1, (sum(probo))) #Sample time to next event
    tminone=t
    t=t+timetoevent
    saver=saver+1

    if (t>wantedtimes[wantedtimescounter] && wantedtimescounter<=length(wantedtimes)){

      keep=c(tminone,pop[,2]) #store abundances
      keeper=rbind (keeper, keep)
      #print (c(tminone,t,"WANTED"))
      wantedtimescounter=wantedtimescounter+1

    }


  }



  if (tmax!=wantedtimes[length(wantedtimes)]){
    keep=c(tminone,pop[,2]) #store abundances
    keeper=rbind (keeper, keep)#put them in dataframe
  }
  return(t(keeper))
}
