#' Multiple Simple Neutral community (SN) simulations, plotting the RAC at specific timepoints and also saving the output matrix for every simulation
#'
#'@param tmax Arbitrary units of time the simulation should be run
#' @param b1 Birthrate
#' @param d1 Deathrate
#' @param k1 Carrying capacity
#' @param abun_original A vector of initial abundances for each species such as those created by generate_spat_abund
#' @param nsim number of simulations that should be performed (integer)
#' @param interval After how many evenst should the population state be saved to the output matrix
#' @param wantedtimes Specific timepoints for which the RAC is to be plotted
#' @author Timo van Eldijk
#'
#' @return this function plots the RAC of the simulations for the times specified in Wantedtimes as .tiff files in the current working directory also writes a .csv file outputting the poplation matrix for every simulation performed
#' @export
#'
#' @examples
#' multinulsimSPEC(10, 0.6, 0.1, 1600,
#' generate_spat_abund(theta = 200,Ivec = rep(40,1),Jvec = c(16000)),
#' 2, 200, c(15,30,50,75,100))
#'
#'
#'

multinulsimSPEC=function (tmax, b1, d1, k1, abun_original, nsim ,interval,wantedtimes){

  counter=1
  while(counter<=nsim){
    print (counter)
    test=nulsimSPEC(tmax, b1, d1, k1, abun_original, interval, wantedtimes)
    utils::write.csv(test, file = paste(counter, "simdata_nulsim.csv", sep=""))
    counter=counter+1

  }


  ##### Plot racs of wantedtimes
  timepointcounter=1
  while(timepointcounter<=length(wantedtimes)){ #Cycle through the timepoints for which you want plots (whole integers!)

    wantedtime=wantedtimes[timepointcounter]
    # 2 Go throught the .csv files one by one to make RAC through time plots
    racstor=data.frame()
    counter=1
    while(counter<=nsim){
      MyData <- utils::read.csv(file=paste(counter, "simdata_nulsim.csv", sep=""), header=TRUE, sep=",")
      readdata=as.matrix(MyData)
      readdata=readdata[,2:ncol(readdata)]
      readdata=apply(readdata, 2, as.numeric)

      timepoint=readdata[1, readdata[1,]<=wantedtime]
      timepoint=max(timepoint)


      rac=readdata [2:nrow(readdata),readdata[1,]==timepoint]


      rac=sort((rac), decreasing = T)
      racstor=rbind(racstor, rac)

      #print (paste(counter, "rac", wantedtime))

      counter=counter+1

    }


    grDevices::tiff(filename=paste(wantedtime,"time_rac_nulsim.tiff", sep=""), width = 1000, height = 1000, units = "px", pointsize = 32)
    graphics::par(mfrow=c(1,1))
    graphics::plot(NA , xlab="Species Rank", ylab="Nr of individuals", ylim=c(1,10000), xlim=c(0,(length(abun_original))),log="y", pch = 20, cex = .8 )

    newracstor=sapply(racstor, function(x){stats::quantile(x,probs=c(0.05,0.25, 0.5, 0.75,0.95), type = 1)})
    vec1=newracstor[1,]
    vec2=newracstor[2,]
    vec3=newracstor[3,]
    vec4=newracstor[4,]
    vec5=newracstor[5,]
    graphics::lines (1:length(vec1[vec1>0]), vec1[vec1>0], cex = .6, col="grey", lwd=2, lty=2)
    graphics::lines (1:length(vec5[vec5>0]), vec5[vec5>0], cex = .6, col="grey", lwd=2, lty=2)
    graphics::lines (1:length(vec2[vec2>0]), vec2[vec2>0], cex = .6, col="grey", lwd=2, lty=3)
    graphics::lines (1:length(vec4[vec4>0]), vec4[vec4>0], cex = .6, col="grey", lwd=2, lty=3)
    graphics::lines (1:length(vec3[vec3>0]),vec3[vec3>0], pch = 20, cex = .8, col="black", lwd=2)
    graphics:: polygon(c((1:length(vec2[vec2>0])), rev(1:length(vec1[vec1>0]))), c(vec2[vec2>0], rev(vec1[vec1>0])),col=grDevices::gray(0.9))
    graphics::polygon(c((1:length(vec5[vec5>0])), rev(1:length(vec4[vec4>0]))), c(vec5[vec5>0], rev(vec4[vec4>0])),col=grDevices::gray(0.9))
    graphics:: polygon(c((1:length(vec3[vec3>0])), rev(1:length(vec2[vec2>0]))), c(vec3[vec3>0], rev(vec2[vec2>0])),col='lightblue')
    graphics:: polygon(c((1:length(vec4[vec4>0])), rev(1:length(vec3[vec3>0]))), c(vec4[vec4>0], rev(vec3[vec3>0])),col="lightblue")
    graphics:: points (1:ncol(racstor), newracstor[3,], pch = 20, cex = 0.8, col="black", lwd=2)
    graphics::points(1:ncol(racstor), sort(abun_original, decreasing=T), col="Red", pch = 20, cex = .8)


    grDevices::dev.off()

    timepointcounter=timepointcounter+1




  }

  #Trajectories
  grDevices::tiff(filename="Trajectory_totalcommsize_nulsim.tiff", width = 1000, height = 1000, units = "px", pointsize = 32)
  graphics::plot(1, type="n", xlab="time", ylab="Total community size", ylim=c(0,20000), xlim=c(0,tmax))

  counter=1 #cycle .csv files
  while(counter<=nsim){
    MyData <- utils::read.csv(file=paste(counter, "simdata_nulsim.csv", sep=""), header=TRUE, sep=",") #reading in .csv
    readdata=as.matrix(MyData)
    readdata=readdata[,2:ncol(readdata)]
    readdata=apply(readdata, 2, as.numeric)

    traj=cbind(readdata[1,], colSums(readdata [2:nrow(readdata), ])) #Calculte total community size trajectory over time

    graphics::points (traj[,1], traj[,2], pch = 20, cex = .8) #plot it all

    #print(paste(counter, "trajectory"))
    counter=counter+1

  }
  grDevices::dev.off()




}
