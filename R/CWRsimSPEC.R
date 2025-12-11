#' Community Wide evolutionary Rescue simulation (CWR), with certain precise timepoints saved
#'
#' @param tmax Arbitrary units of time the simulation should be run
#' @param b1 Birthrate of the mutants
#' @param d1 Deathrate of the mutants
#' @param b2 Birthrate of the residents
#' @param d2 Death rate of the residents
#' @param m12 Mutation rate of mutants to residents
#' @param m21 Mutation rate of residents to mutants
#' @param k1 Carrying capacity
#' @param abun_original A vector of initial abundances for each species such as those created by generate_spat_abund
#' @param interval After how many evenst should the population state be saved to the output matrix
#' @param wantedtimes The preciese timepoints that should be saved
#' @return A matrix denoting the abundances of all the species in the community over time, with certain precise timpoints saved
#' @export
#'
#' @examples
#' CWRsimSPEC(10,0.6 , 0.1, 0.05, 0.1, 0, 0.0005, 16000,
#' generate_spat_abund(theta = 200,Ivec = rep(40,1),Jvec = c(16000)),
#'  200,c(15,30,50,75,100))
#'


CWRsimSPEC = function (tmax, b1, d1, b2, d2, m12, m21, k1, abun_original, interval, wantedtimes){

  print(Sys.time())

  wantedtimescounter=1

  nspec=length(abun_original) #Get nr of species in start community
  abun_mutant=(rep(0, nspec)) #Set mutant abundances for each species to 0
  bo=c(rep(b2,nspec)) #Birth rate vector originals
  do=c(rep(d2,nspec)) #death rate vector originals
  bm=c(rep(b1,nspec)) #birth rate vector mutants
  dm=c(rep(d1,nspec)) #Death rate vector muntants

  om=c(rep(m12,nspec)) #Vector of mutation probabilities for original to mutant
  mo=c(rep(m21,nspec)) #Vector of mutation probabilities for mutant to originals

  k=c(rep(k1,(nspec))) #Vector of K's
  specnrs=1:nspec # Number the species
  pop=data.frame()
  pop=specnrs
  pop=cbind(pop, abun_original)#Make total community dataframe with orignals
  pop=cbind(pop, abun_mutant) #add the mutant

  t=0 #Set time to zero

  keeper=data.frame() #initialise data frame for storage
  keep=c(t,pop[,2], pop[,3]) #store initial abundances
  keeper=rbind (keeper, keep)#put them in dataframe
  saver=0 #saver to 0
  #Interval of steps after which saving occurs std 400


  totpopimpact=sum(sum(pop[,2]), sum(pop[,3]))/k #how much of carrying capacity is ocupied


  probo=(max(0,(bo*(1-totpopimpact)))+do+om)*pop[,2] #Vector of event probabilities for originals
  probm=(max(0,(bm*(1-totpopimpact)))+dm+mo)*pop[,3] #Vector of event probabilities for mutants
  timetoevent=stats::rexp(1, (sum(probo)+sum(probm))) #Sample time to next event
  t=t+timetoevent #Add time to next event to time



  while (t<tmax & sum(sum(pop[,2]), sum(pop[,3]))!=0){

    species=sample (c(1:(nspec*2)), size=1, replace=T, prob=(c(probo,probm))) #Sample which species has event occur and if this is a mutant or an original


    if (species<=nspec){ #If it is an original
      event=sample(c(-1,1,2), size=1, replace=T, prob = c(do[species], max(0, (bo[species]*(1-totpopimpact[species]))), om[species])) #is event a birth, death or mutation
      if (event<2){pop[species, 2]=pop[species, 2]+event}#implement birth or death
      if (event==2){pop[species, 2]=pop[species, 2]-1; pop[species, 3]=pop[species, 3]+1}#implement mutation
    }


    if (species>nspec){ #if it is a mutant
      event=sample(c(-1,1,2), size=1, replace=T, prob = c(dm[(species-nspec)],max(0, (bm[species-nspec]*(1-totpopimpact[species-nspec]))), mo[(species-nspec)]))#is event a birth, death or mutation
      if (event<2){pop[(species-nspec),3]=pop[(species-nspec), 3]+event} #implement birth or death
      if (event==2){pop[(species-nspec), 2]=pop[(species-nspec), 2]+1; pop[(species-nspec), 3]=pop[(species-nspec), 3]-1} #implement mutation
    }




    totpopimpact=sum(sum(pop[,2]), sum(pop[,3]))/k
    probo=(max(0,(bo*(1-totpopimpact)))+do+om)*pop[,2] #Vector of event probabilities for originals
    probm=(max(0,(bm*(1-totpopimpact)))+dm+mo)*pop[,3]#Vector of event probabilities for mutants

    timetoevent=stats::rexp(1, (sum(probo)+sum(probm)))#sample time to next event
    tminone=t
    t=t+timetoevent  #advance time by time to next event


    #Save certain dedicated times
    if (t>wantedtimes[wantedtimescounter] && wantedtimescounter<=length(wantedtimes)){

      keep=c(tminone,pop[,2], pop[,3]) #store abundances
      keeper=rbind (keeper, keep)
      #print (c(tminone,t,"WANTED"))
      wantedtimescounter=wantedtimescounter+1

    }


    saver=saver+1 #increase saver counter by 1
  }

  #SO AT THE END STORING FINAL STATE

  if (tmax!=wantedtimes[length(wantedtimes)]){
    keep=c(tminone,pop[,2], pop[,3]) #store abundances
    keeper=rbind (keeper, keep)#put them in dataframe
    #print (c(tminone, "finalpoint"))
  }

  return(t(keeper))
}

