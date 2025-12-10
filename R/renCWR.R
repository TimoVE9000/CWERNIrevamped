
#' Perform CWR simulation, Calculate Reny entropy (alfa=1, alfa=2), re-simulate neutral communities using SADISA parameters estimated on CWR simulation,
#' simulation neutral communities, calculate their Reny entropy (alfa=1, alfa=2) to generate a distribution.
#'
#'
#' @param tmax Arbitrary units of time the simulation should be run
#' @param b1 Birthrate of the mutants
#' @param d1 Deathrate of the mutants
#' @param b2 Birthrate of the residents
#' @param d2 Death rate of the residents
#' @param m12 Mutation rate of mutants to residents
#' @param m21 Mutation rate of residents to mutants
#' @param k1 Carrying capacity
#' @param interval After how many evenst should the population state be saved to the output matrix
#' @param orgsimmaxcount how many CWR simulations should be performed
#' @param totalresamp how many neutral communities should be simulated for each single CWR simulation
#' @author Timo van Eldijk
#' @return This function saves two .pdf plots for each CWR simulation performed, showing the Reny entropy of each CWR simulation and the distribution of Reny entorpy generated from the neutral re-simulations
#' @export
#'
#' @examples
#' renCWR(10, 0.6, 0.1, 16000,200, 0.05,0.1, 0, 0.0005, 1,1)
#'



renCWR=function(tmax, b1, d1, k1,interval, b2, d2, m12, m21, orgsimmaxcount, totalresamp){

  storallloglik=data.frame()
  storalltheta=data.frame()
  storalli=data.frame()
  orgsimcounter=1
  while (orgsimcounter<=orgsimmaxcount){
    print("start")
    #Generate initial abundances
    abun_original=generate_spat_abund(theta = 200,Ivec = rep(40,1),Jvec = c(16000))
    #plot(sort(abun_original, decreasing=T), log="y")


    #Do original CWR simulation
    x=CWRsim(tmax, b1, d1, b2, d2, m12, m21, k1, abun_original, interval) #simulate
    abundmut=x[(length(abun_original)+2):((length(abun_original)*2)+1),ncol(x)]
    abundres=x[2:((length(abun_original))+1),ncol(x)]
    abund=abundmut+abundres
    abund=abund [abund>0]
    # plot(sort (abund, decreasing=T),log= "y" , col="Red")
    #points (sort(abun_original, decreasing = T))



    #Do SADISA parameter estimation
    initpars=c(40, 2000)
    labelpars=c(1,1)
    idpars = c(1,1)
    pars=SADISA::SADISA_ML(abund,initpars, idpars,labelpars, model = c("pm", "dl"),
                           tol = c(1e-06, 1e-06, 1e-06), maxiter = 2000, optimmethod = "subplex")
    theta=pars$pars[1]
    I=pars$pars[2]
    originalloglik=pars$loglik
    originaltheta=pars$pars[1]
    originali=pars$pars[2]




    #Fit Renyi stuff

    #add the alpha =0 case

    ren1org=EntropyEstimation::Renyi.z(abund,1)
    ren2org=EntropyEstimation::Renyi.z(abund,2)


    #Simulate using estimated params
    storren1sim=numeric(0)
    storren2sim=numeric(0)
    currentresamp=1


    #pdf(paste(orgsimcounter,"rac.pdf"),
    #width=5,
    #height=4,
    #pointsize=12)
    #plot(sort (abund, decreasing=T),log='y' , col="Red")
    #points (sort(abun_original, decreasing = T), col="Blue")

    while (currentresamp<=totalresamp){
      simulation=generate_spat_abund(theta = theta,Ivec = rep(I,1),Jvec = c(16000))


      # points(sort (simulation, decreasing=T))

      #Fit renyi stuff to resimulation



      #plot(fits)
      ren1sim=EntropyEstimation::Renyi.z(simulation,1)
      ren2sim=EntropyEstimation::Renyi.z(simulation,2)


      storren1sim=c(storren1sim, ren1sim)
      storren2sim=c(storren2sim, ren2sim)
      print(currentresamp)
      currentresamp=currentresamp+1
    }

    #points(sort (abund, decreasing=T), col="Red")
    #dev.off()




    #Plot each renyi fit thing
    grDevices::pdf(paste(orgsimcounter,"_Ren1racs_rees_CWR.pdf", sep=""),
                   width=5,
                   height=4,
                   pointsize=12)
    graphics::plot(graphics::hist(storren1sim),xlim=c(0,5), main="Distribution of Ren1 re-estimates")
    graphics::abline ( v=ren1org , col="Red")
    grDevices::dev.off()

    #Plot each renyi fit thing
    grDevices::pdf(paste(orgsimcounter,"_Ren2racs_rees_CWR.pdf", sep=""),
                   width=5,
                   height=4,
                   pointsize=12)
    graphics::plot(graphics::hist(storren2sim),xlim=c(0,5), main="Distribution of Ren2 re-estimates")
    graphics::abline ( v=ren1org , col="Red")
    grDevices::dev.off()




    thisitloglik=c(originalloglik)
    thisittheta=c(originaltheta)
    thisiti=c(originali)

    storallloglik=rbind (storallloglik, thisitloglik)
    storalltheta=rbind (storalltheta, thisittheta)
    storalli=rbind (storalli, thisiti)

    print (paste ("Main",orgsimcounter))
    orgsimcounter=orgsimcounter+1
  }


  utils::write.table(storallloglik, file = "loglig_racs_rees_renCWR.csv", sep = ",", append = F, quote = TRUE,
                     col.names = F, row.names = FALSE)

  utils::write.table(storalltheta, file = "theta_racs_rees_renCWR.csv", sep = ",", append = F, quote = TRUE,
                     col.names = F, row.names = FALSE)

  utils::write.table(storalli, file = "i_loglik_racs_renCWR.csv", sep = ",", append = F, quote = TRUE,
                     col.names = F, row.names = FALSE)

}
