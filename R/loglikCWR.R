


#' Perform CWR simulation, estimate paramters using SADISA, re-simulate neutral communities to generate distribution of likelihoods
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
#' @param tolsadisa tolerances for SADISA estimation vector of 3 values, default c(1e-06, 1e-06, 1e-06)
#' @author Timo van Eldijk
#' @return This function saves a .pdf plot for each CWR simulation performed, showing the distribution of likelihoods
#' @export
#'
#' @examples
#' loglikCWR(5, 0.6, 0.1, 16000,200, 0.05,0.1, 0, 0.0005, 1,1,c(1e-1, 1e-1, 1e-1))
#'
#'
#'
#'
loglikCWR=function(tmax, b1, d1, k1,interval, b2, d2, m12, m21, orgsimmaxcount, totalresamp, tolsadisa=c(1e-06, 1e-06, 1e-06)){

  storallloglik=data.frame()
  storalltheta=data.frame()
  storalli=data.frame()
  orgsimcounter=1
  while (orgsimcounter<=orgsimmaxcount){
    print("start")

    #Generate original abundance and simulate CWR
    abun_original=generate_spat_abund(theta = 200,Ivec = rep(40,1),Jvec = c(16000))
    x=CWRsim(tmax, b1, d1, b2, d2, m12, m21, k1, abun_original, interval)
    abundmut=x[(length(abun_original)+2):((length(abun_original)*2)+1),ncol(x)]
    abundres=x[2:((length(abun_original))+1),ncol(x)]
    abund=abundmut+abundres
    abund=abund [abund>0]


    #Load sadisa library and do SADISA estimation
    initpars=c(40, 2000)
    labelpars=c(1,1)
    idpars = c(1,1)
    pars=SADISA::SADISA_ML(abund,initpars,idpars,labelpars, model = c("pm", "dl"),
                           tol = tolsadisa, maxiter = 2000, optimmethod = "subplex")
    #Store estimated parameters (on original CWR rac)
    originaltheta=pars$pars[1]
    originali=pars$pars[2]
    originalloglik=pars$loglik


    #Simulate neutral communities using the estimated parameters


    currentresamp=1
    saver=vector('numeric')
    savertheta=vector('numeric')
    saveri=vector('numeric')
    #pdf(paste(orgsimcounter,paste(orgsimcounter,"rac.pdf")),
    #width=5,
    #height=4,
    #pointsize=12)

    #plot(sort (abund, decreasing=T), log="y", col='Red')
    #points(sort(abun_original, decreasing=T), col="Blue")

    while  (currentresamp <=totalresamp){
      #Simulate using estimated params
      simulation=generate_spat_abund(theta = originaltheta,Ivec = rep(originali,1),Jvec = c(16000))
      #points(sort(simulation, decreasing=T))

      #Re-estimate parameters using SADISA
      initpars=c(20, 20)
      labelpars=c(1,1)
      idpars = c(1,1)
      comp=SADISA::SADISA_ML(simulation,initpars,idpars,labelpars, model = c("pm", "dl"),
                             tol = tolsadisa, maxiter = 2000, optimmethod = "subplex")
      saver=c(saver,comp$loglik)
      savertheta=c(savertheta,comp$pars[1])
      saveri=c(saveri,comp$pars[2])
      print(currentresamp)
      currentresamp=currentresamp+1
    }
    #dev.off()



    #Plot the distribution of logliklihoods
    grDevices::pdf(paste(orgsimcounter,"gsloglik_loglik_rees_CWR.pdf"),
                   width=5,
                   height=4,
                   pointsize=12)
    graphics::hist((saver), main="Loglikelihood of SADISA fit (pm, dl)",xlim=c(-100,-500))
    graphics::abline ( v= originalloglik, col="Red")
    grDevices::dev.off()



    #Save all the outcomes of all the SADISA estimates
    thisitloglik=c(originalloglik, saver)
    thisittheta=c(originaltheta,savertheta )
    thisiti=c(originali, saveri)
    storallloglik=rbind (storallloglik, thisitloglik)
    storalltheta=rbind (storalltheta, thisittheta)
    storalli=rbind (storalli, thisiti)

    print (paste ("Main",orgsimcounter))
    orgsimcounter=orgsimcounter+1
  }


  #Write the saved SADISA outcomes to .csv files
  utils::write.table(storallloglik, file = "loglig_loglik_rees_CWR.csv", sep = ",", append = F, quote = TRUE,
                     col.names = F, row.names = FALSE)
  utils::write.table(storalltheta, file = "theta_loglik_rees_CWR.csv", sep = ",", append = F, quote = TRUE,
                     col.names = F, row.names = FALSE)
  utils::write.table(storalli, file = "i_loglik_rees_CWR.csv", sep = ",", append = F, quote = TRUE,
                     col.names = F, row.names = FALSE)

}
