#' Generate a neutral community
#' @description Uses sampling function for neutral community from Etienne et al. 2007 (Ecology Letters), used to initialize communities.
#' @param theta an integer that gives a value to Theta
#' @param Ivec  a vector of I values
#' @param Jvec  a vector of J values
#'
#' @return a vector of species abundances
#' @export
#'
#' @examples generate_spat_abund(theta = 200,Ivec = rep(40,1),Jvec = c(16000))

generate_spat_abund = function(theta,Ivec,Jvec)
{
  numsam = length(Jvec)
  J = sum(Jvec)
  locspecnum = matrix(0,nrow = numsam,ncol = J)
  mcspecnum = rep(0,J)
  abund = matrix(0,nrow = numsam,ncol = J)
  k = 0
  n = 0
  for(i in 1:numsam)
  {
    for(j in 1:Jvec[i])
    {
      bnd = Ivec[i]/(Ivec[i] + j - 1)
      if(stats::runif(1) > bnd)
      {
        locspecnum[i,j] = locspecnum[i,sample(1:j,1)]
      } else
      {
        k = k + 1
        if(stats::runif(1) <= theta/(theta + k - 1))
        {
          n = n + 1
          mcspecnum[k] = n
        } else
        {
          mcspecnum[k] = mcspecnum[sample(1:k,1)]
        }
        locspecnum[i,j] = mcspecnum[k]
      }
      abund[i,locspecnum[i,j]] = abund[i,locspecnum[i,j]] + 1
    }
  }
  zeros = 1
  numspec = J
  while(zeros == 1)
  {
    sumsites = sum(abund[1:numsam,numspec])
    if(sumsites == 0)
    {
      numspec = numspec - 1
    } else
    {
      zeros = 0
    }
  }
  abund = abund[1:numsam,1:numspec]
  return(abund)
}
