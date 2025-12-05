#' Simulate multiple neutral communities with a variable birthrate (VBN)
#'
#'
#'
#'@param tmax Arbitrary units of time the simulation should be run
#' @param b1 Normal birthrate
#' @param d1 Deathrate
#' @param k1 Carrying capacity
#' @param bneck Birthrate during the bottleneck
#' @param kneckstart Timepoint at which the bottleneck starts
#' @param kneckend  Timepoint at which the bottleneck ends
#' @param abun_original A vector of initial abundances for each species such as those created by generate_spat_abund
#' @param interval After how many evenst should the population state be saved to the output matrix
#' @param nsim number of simulations that should be performed
#' @param filenomen the name of the plot that is saved
#'
#'@author Timo van Eldijk
#' @return this function plots all the results of the simulations in a .tiff file in the current working directory
#' @export
#'
#'
#' @examples
#' multisimbvar( 10, 0.6, 0.1, 16000, 0.05,0, 5,
#' generate_spat_abund (theta = 200,Ivec = rep(40,1),Jvec = c(16000)),
#' 2, "bvar.tiff", 200)
#'
#'
multisimbvar=function (tmax, b1, d1, k1, bneck,kneckstart, kneckend, abun_original,nsim, filenomen, interval){
  #plot(1, type="n", xlab="rank", ylab="individuals", ylim=c(1,7000), xlim=c(0,((nrow(test)-1)/2)), log="y")
  #plot(1, type="n", xlab="time", ylab="individuals", ylim=c(0,20000), xlim=c(0,tmax))
  racstor=data.frame()
  trajstor=data.frame()
  count=0
  while (count<nsim){ #Preform nsim number of simulations
    test=simbvar(tmax, b1, d1, k1, bneck,kneckstart, kneckend, abun_original, interval) #simulate
    rac=sort(test [(2):(nrow(test)),length (test[1,])],decreasing = T) #calculate RAC
    racstor=rbind (racstor, rac) #Store RAC
    #points (rac)
    traj=cbind(test[1,], colSums(test [2:nrow(test), ])) #Calculate community size trajetory over time
    trajstor=rbind(trajstor, traj) #store commmunity size trajectory over time
    #points (traj[,1], traj[,2])
    count=count+1
  }

  #Plot the results
  grDevices::tiff(filename=filenomen, width = 1000, height = 1000, units = "px", pointsize = 21 )
  graphics::par(mfrow=c(2,2))


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
  graphics:: lines (1:length(vec4[vec4>0]), vec4[vec4>0], cex = .6, col="grey", lwd=2, lty=3)
  graphics::lines (1:length(vec3[vec3>0]),vec3[vec3>0], pch = 20, cex = .8, col="black", lwd=2)
  graphics::polygon(c((1:length(vec2[vec2>0])), rev(1:length(vec1[vec1>0]))), c(vec2[vec2>0], rev(vec1[vec1>0])),col=grDevices::gray(0.9))
  graphics::polygon(c((1:length(vec5[vec5>0])), rev(1:length(vec4[vec4>0]))), c(vec5[vec5>0], rev(vec4[vec4>0])),col=grDevices::gray(0.9))
  graphics::polygon(c((1:length(vec3[vec3>0])), rev(1:length(vec2[vec2>0]))), c(vec3[vec3>0], rev(vec2[vec2>0])),col='lightblue')
  graphics::polygon(c((1:length(vec4[vec4>0])), rev(1:length(vec3[vec3>0]))), c(vec4[vec4>0], rev(vec3[vec3>0])),col="lightblue")
  graphics::points (1:ncol(racstor), newracstor[3,], pch = 20, cex = 0.8, col="black", lwd=2)
  graphics::points(1:ncol(racstor), sort(abun_original, decreasing=T), col="Red", pch = 20, cex = .8)




  graphics::plot(1, type="n", xlab="time", ylab="Total community size", ylim=c(0,20000), xlim=c(0,tmax))
  graphics::points (trajstor[,1], trajstor[,2], pch = 20, cex = .8)
  grDevices::dev.off()
}
