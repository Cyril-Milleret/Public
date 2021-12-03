#' Density and random generation for the categorical distribution of state transition with one alive state and 1 dead state
#' 
#'
#' The \code{dcatState1Dead} distribution is a NIMBLE custom distribution which can be used to model and simulate
#' individual state transition. It can be used in cases with one alive state and 1 dead states. 
#' If z_{i,t} = 1, individual i can be recruited (transition to state 2) with probability gamma_t, so z_{i,t+1} ~dcat(1- gamma_t, gamma_t, 0 , 0) where gamma_(t ) represent the probability of an unborn individual to be recruited.
#' If z_{i,t} = 2, individual i can survive with probability \phi and remain z_{i,t+1}=2 or die with probability 1- \phi and transition to z_{i,t+1}=3, the absorbing state. 
#' If gamma or phi are assumed to be spatially heterogeneous, a vector of probability should be provided for gammaSpatial or phiSpatial.
#' if gamma or phi are assumed to be spatially homogeneous, a scalar should be provided for gamma or phi.
#' 
#' 
#' @name dcatState1Dead 
#' 
#' @param x Scalar of the individual state z_{i,t+1}.
#' @param n Integer specifying the number of realizations to generate.  Only n = 1 is supported.
#' @param z Scalar of the initial individual state z_{i,t}
#' @param phi Scalar with probability phi to transition from z_{i,t} = 2 to z_{i,t+1} = 2.
#' @param phiSpatial Vector with probability phi_r to transition from z_{i,t} = 2 to z_{i,t+1} = 2 for each cell r of the habitat. The vector should be of the length of the number of habitat windows in \code{\link{habitatGrid}} .
#' @param habitatGrid Matrix of habitat window indices and only used if arguements hSpatial, wSpatial or phiSpatial are used.
#' Habitat window indices should match the order in \code{phiSpatial}, \code{hSpatial}, or \code{wSpatial}. 
#' @param log Logical argument, specifying whether to return the log-probability of the distribution.
#' @return 
#' \code{dcatState1Dead} gives the (log) probability density of \code{x}. 
#' \code{rcatState2Dead} gives a randomly generated individual states conditional on the initial state \code{z}.  
#' 
#' @author Cyril Milleret
#' 
#' 
#' @references
#' 
#' @examples
#' # Use the distribution in R
#' z <- 2
#' gamma <- 0.2
#' phi <- 0.7
#' 
#' lowerCoords <- matrix(c(0, 0, 1, 0, 0, 1, 1, 1), nrow = 4, byrow = TRUE)
#' upperCoords <- matrix(c(1, 1, 2, 1, 1, 2, 2, 2), nrow = 4, byrow = TRUE)  
#' logIntensities <- log(rep(1,4))
#' logSumIntensity <- log(sum(c(1:4))) 
#' habitatGrid <- matrix(c(1:4), nrow = 2, byrow = TRUE)
#' numGridRows <- nrow(habitatGrid)
#' numGridCols <- ncol(habitatGrid)
#' s <- rbernppAC(n=1, lowerCoords, upperCoords, logIntensities, logSumIntensity, 
#'                habitatGrid, numGridRows, numGridCols)
#' 
#' ## No spatial mortality 
#' zPlusOne <- rcatState1Dead( z = z
#'                             , gamma = gamma
#'                             , phi = phi
#'                             , s = s
#'                             , habitatGrid = habitatGrid)
#' zPlusOne
#' #
#' dcatState1Dead(  x = zPlusOne
#'                  , z = z
#'                  , gamma = gamma
#'                  , phi = phi
#'                  , s = s
#'                  , habitatGrid = habitatGrid)
#' 
#' ##  With spatial mortality
#' phiSpatial <- c(0.60, 0.70, 0.74, 0.65)
#' gammaSpatial <- c(0.4,0.5,0.1,0.3)
#' zPlusOne <- rcatState1Dead( z = z
#'                             , gammaSpatial = gammaSpatial
#'                             , phiSpatial = phiSpatial
#'                             , s = s
#'                             , habitatGrid = habitatGrid)
#' zPlusOne
#' dcatState1Dead(  x = zPlusOne
#'                  , z = z
#'                  , gammaSpatial = gammaSpatial
#'                  , phiSpatial = phiSpatial
#'                  , s = s
#'                  , habitatGrid = habitatGrid)
#' 
#' 
#' 
#' 
#' 
NULL
#' @rdname dcatState1Dead
#' @export



#### 1.Density function ####
dcatState1Dead  <- nimbleFunction(run = function( x = double(0)
                                                  , z = double(0)
                                                  , gamma = double(0, default=-999)
                                                  , gammaSpatial = double(1)
                                                  , phi = double(0, default=-999)
                                                  , phiSpatial = double(1)
                                                  , s = double(1)
                                                  , habitatGrid = double(2)
                                                  , log = integer(0, default = 0)){
  # Return type declaration
  returnType(double(0))
  
  if(z == 1){
    
    if(gamma == -999){
      sID <- habitatGrid[trunc(s[2])+1, trunc(s[1])+1]
      IndGamma <- gammaSpatial[sID]
      
    }else{
      ## EXTRACT LOCATION OF THE ID
      IndGamma <- gamma
    }
    
    logLikelihood <- dcat(x, prob = c(1 - IndGamma, IndGamma), log=1)
    if(log == 1){return(logLikelihood)}else{return(exp(logLikelihood))}
  }
  
  if(z == 2){
    ## EXTRACT LOCATION OF THE ID
    
    sID <- habitatGrid[trunc(s[2])+1, trunc(s[1])+1]
   
    #phi
    if(phi == -999){
      Indphi <- phiSpatial[sID]
    }else{
      Indphi <- phi
    }
    
    logLikelihood <- dcat(x, prob = c(0, Indphi, 1-Indphi), log=1)
    if(log == 1){return(logLikelihood)}else{return(exp(logLikelihood))}
  }
  
  if(z == 3){
    logLikelihood <- dcat(x, prob = c(0, 0, 1), log=1)
    if(log == 1){return(logLikelihood)}else{return(exp(logLikelihood))}
  }
  
  
})


NULL
#' @rdname rcatState1Dead
#' @export
#' 
#### 2.Sampling function ####
rcatState1Dead <- nimbleFunction(run = function( n = integer(0)
                                                 , z = double(0)
                                                 , gamma = double(0, default=-999)
                                                 , gammaSpatial = double(1)
                                                 , phi = double(0, default=-999)
                                                 , phiSpatial = double(1)
                                                 , s = double(1)
                                                 , habitatGrid = double(2)
){
  # Return type declaration
  returnType(double(0))
  
  if(z == 1){
    if(gamma == -999){
      sID <- habitatGrid[trunc(s[2])+1, trunc(s[1])+1]
      IndGamma <- gammaSpatial[sID]
    }else{
      ## EXTRACT LOCATION OF THE ID
      IndGamma <- gamma
    }
    state <- rcat(1, prob = c(1 - IndGamma, IndGamma, 0))
    return(state)
  }
  
  if(z == 2){
    ## EXTRACT LOCATION OF THE ID
    sID <- habitatGrid[trunc(s[2])+1, trunc(s[1])+1]
    #phi
    if(phi == -999){
      Indphi <- phiSpatial[sID]
    }else{
      Indphi <- phi
    }
    
    state <- rcat(1, prob = c(0, Indphi, 1-Indphi))
    return(state)
  }
  
  if(z == 3 ){
    state <- 3
    return(state)
  }
  
})

#### 3.Registration ####
registerDistributions(list(
  dcatState1Dead = list(
    BUGSdist = "dcatState1Dead(z, gamma       , gammaSpatial    , phi       , phiSpatial    , s, habitatGrid)",
    # question: how we can make the distribution work when "vec" is not in used and we dont provide habtiatGrid? We need to give habitatGrid=double(2)...
    Rdist = c( "dcatState1Dead(z, gamma       , gammaSpatial = s, phi       , phiSpatial = s, s, habitatGrid)",
               "dcatState1Dead(z, gamma       , gammaSpatial = s, phi = -999, phiSpatial    , s, habitatGrid)",
               "dcatState1Dead(z, gamma = -999, gammaSpatial    , phi       , phiSpatial = s, s, habitatGrid)",
               "dcatState1Dead(z, gamma = -999, gammaSpatial    , phi = -999, phiSpatial    , s, habitatGrid)"
    ),
    
    
    types = c( "value = double(0)", "z = double(0)","gamma = double(0)", "gammaSpatial = double(1)", "phi = double(0)", 
               "phiSpatial = double(1)", "s = double(1)", "habitatGrid = double(2)"
    ),
    discrete = TRUE,
    mixedSizes = TRUE,
    pqAvail = FALSE
  )))


