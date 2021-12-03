#' Density and random generation for the categorical distribution of state transition with one alive state and two dead state
#' 
#'
#' The \code{dcatState2Dead} distribution is a NIMBLE custom distribution which can be used to model and simulate
#' individual state transition. It can be used in cases with one alive state and two dead states. 
#' If z_{i,t} = 1, individual i can be recruited (transition to state 2) with probability gamma_t, so z_{i,t+1} ~dcat(1- gamma_t, gamma_t, 0 , 0) where gamma_(t ) represent the probability of an “unborn” individual to be recruited.
#' If z_{i,t} = 2, individual i can survive with probability \phi and remain z_{i,t+1}=2. If it does not survive, it can either die due to culling (or any other causes) and be recovered (transition to z_{i,t+1}=3) with probability hi, or die from other causes without being recovered (transition to z_{i,t+1} = 4) with probability w, so that z_{i,t+1} ~ dcat(0, \phi, h, w), where \phi = 1−h −w. 
#' All individuals in dead states (z_{i,t} = 3 or 4) transition to z_{i,t+1} = 4, the absorbing state, with probability 1.
#' If gamma, w, h or phi are assumed to be spatially heterogeneous, a vector of probability should be provided for gammaSpatial, wSpatial, hSpatial or phiSpatial.
#' if gamma, w, h or phi are assumed to be spatially homogeneous, a scalar should be provided for gamma, w, h or phi.
#' 
#' 
#' @name dcatState2Dead 
#' 
#' @param x Scalar of the individual state z_{i,t+1}.
#' @param n Integer specifying the number of realizations to generate.  Only n = 1 is supported.
#' @param z Scalar of the initial individual state z_{i,t}
#' @param h Scalar with probability h to transition from z_{i,t} = 2 to z_{i,t+1} = 3.
#' @param w Scalar with probability w to transition from z_{i,t} = 2 to z_{i,t+1} = 4.
#' @param phi Scalar with probability phi to transition from z_{i,t} = 2 to z_{i,t+1} = 2.
#' @param hSpatial Vector with probability h_r to transition from z_{i,t} = 2 to z_{i,t+1} = 3 for each cell r of the habitat. The vector should be of the length of the number of habitat windows in \code{\link{habitatGrid}} .
#' @param wSpatial Vector with probability w_r to transition from z_{i,t} = 2 to z_{i,t+1} = 4 for each cell r of the habitat. The vector should be of the length of the number of habitat windows in \code{\link{habitatGrid}} .
#' @param phiSpatial Vector with probability phi_r to transition from z_{i,t} = 2 to z_{i,t+1} = 2 for each cell r of the habitat. The vector should be of the length of the number of habitat windows in \code{\link{habitatGrid}} .
#' @param habitatGrid Matrix of habitat window indices and only used if arguements hSpatial, wSpatial or phiSpatial are used.
#' Habitat window indices should match the order in \code{phiSpatial}, \code{hSpatial}, or \code{wSpatial}. 
#' @param log Logical argument, specifying whether to return the log-probability of the distribution.
#' @return 
#' \code{dcatState2Dead} gives the (log) probability density of \code{x}. 
#' \code{rcatState2Dead} gives a randomly generated individual states conditional on the initial state \code{z}.  
#' 
#' @author Cyril Milleret
#' 
#' 
#' @references
#' 
#' @examples
#' # Use the distribution in R
#' 
#' z <- 2
#' gamma <- 0.2
#' h <- 0.4
#' w <- 0.1
#' phi <- 1-(h+w)
#' 
#' lowerCoords <- matrix(c(0, 0, 1, 0, 0, 1, 1, 1), nrow = 4, byrow = TRUE)
#' upperCoords <- matrix(c(1, 1, 2, 1, 1, 2, 2, 2), nrow = 4, byrow = TRUE)  
#' logIntensities <- log(rep(1,4))
#' logSumIntensity <- log(sum(c(1:4))) 
#' habitatGrid <- matrix(c(1:4), nrow = 2, byrow = TRUE)
#' numGridRows <- nrow(habitatGrid)
#' numGridCols <- ncol(habitatGrid)
#' s <- rbernppAC(n=1, lowerCoords, upperCoords, logIntensities, logSumIntensity, 
#'                  habitatGrid, numGridRows, numGridCols)
#' 
#' ## No spatial mortality 
#' zPlusOne <- rcatState2Dead( z = z
#'                               , gamma = gamma
#'                               , h = h
#'                               , w = w
#'                               , phi = phi
#'                               , s = s
#'                               , habitatGrid = habitatGrid)
#' #
#' dcatState2Dead(  x = zPlusOne
#'                  , z = z
#'                  , gamma = gamma
#'                  , h = h
#'                  , w = w
#'                  , phi = phi
#'                  , s = s
#'                  , habitatGrid = habitatGrid)
#' 
#' ##  With spatial mortality
#' hSpatial <- c(0.10, 0.20, 0.15, 0.30)
#' wSpatial <- c(0.13, 0.21, 0.12, 0.08)
#' phiSpatial <- 1-(hSpatial+wSpatial)
#' zPlusOne <- rcatState2Dead( z = z
#'                               , gammaSpatial = gammaSpatial
#'                               , hSpatial =  hSpatial
#'                               , wSpatial = wSpatial
#'                               , phiSpatial = phiSpatial
#'                               , s = s
#'                               , habitatGrid = habitatGrid)
#' 
#' dcatState2Dead(  x = zPlusOne
#'                  , z = z
#'                  , gammaSpatial = gammaSpatial
#'                  , hSpatial =  hSpatial
#'                  , wSpatial = wSpatial
#'                  , phiSpatial = phiSpatial
#'                  , s = s
#'                  , habitatGrid = habitatGrid)
#' 
#' 
#' 
#' 
NULL
#' @rdname dcatState2Dead
#' @export



#### 1.Density function ####
dcatState2Dead  <- nimbleFunction(run = function( x = double(0)
                                                  , z = double(0)
                                                  , gamma = double(0, default=-999)
                                                  , gammaSpatial = double(1)
                                                  , h = double(0, default=-999)
                                                  , hSpatial = double(1)
                                                  , w = double(0, default=-999)
                                                  , wSpatial = double(1)
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
    #h 
    if(h == -999){
      IndH <- hSpatial[sID]
        }else{
      IndH <- h
    }
    #w
    if(w == -999){
        IndW <- wSpatial[sID]
      }else{
        IndW <- w
    }
    #phi
    if(phi == -999){
        Indphi <- phiSpatial[sID]
      }else{
        Indphi <- phi
    }

    logLikelihood <- dcat(x, prob = c(0, Indphi, IndH, IndW), log=1)
    if(log == 1){return(logLikelihood)}else{return(exp(logLikelihood))}
  }
  
  if(z == 3 | z == 4 ){
    logLikelihood <- dcat(x, prob = c(0, 0, 0, 1), log=1)
    if(log == 1){return(logLikelihood)}else{return(exp(logLikelihood))}
  }
  
  
})


NULL
#' @rdname rcatState2Dead
#' @export
#' 
#### 2.Sampling function ####
rcatState2Dead <- nimbleFunction(run = function( n = integer(0)
                                                 , z = double(0)
                                                 , gamma = double(0, default=-999)
                                                 , gammaSpatial = double(1)
                                                 , h = double(0, default=-999)
                                                 , hSpatial = double(1)
                                                 , w = double(0, default=-999)
                                                 , wSpatial = double(1)
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
    state <- rcat(1, prob = c(1 - IndGamma, IndGamma))
    return(state)
  }
  
  if(z == 2){
    ## EXTRACT LOCATION OF THE ID
    sID <- habitatGrid[trunc(s[2])+1, trunc(s[1])+1]
    #h 
    if(h == -999){
      IndH <- hSpatial[sID]
      }else{
        IndH <- h
      }
    #w
    if(w == -999){
      IndW <- wSpatial[sID]
      }else{
        IndW <- w
      }
    #phi
    if(phi == -999){
      Indphi <- phiSpatial[sID]
      }else{
        Indphi <- phi
      }
    
    state <- rcat(1, prob = c(0, Indphi, IndH, IndW))
    return(state)
  }
  
  if(z == 3 | z == 4 ){
    state <- 4
    return(state)
  }
  
})

#### 3.Registration ####
registerDistributions(list(
  dcatState2Dead = list(
    BUGSdist = "dcatState2Dead(z, gamma       , gammaSpatial    , h       , hSpatial    , w       , wSpatial    , phi       , phiSpatial    , s, habitatGrid)",
    # question: how we can make the distribution work when "vec" is not in used and we dont provide habtiatGrid? We need to give habitatGrid=double(2)...
    Rdist = c( "dcatState2Dead(z, gamma       , gammaSpatial = s, h       , hSpatial = s, w       , wSpatial = s, phi       , phiSpatial = s, s, habitatGrid)",
               "dcatState2Dead(z, gamma       , gammaSpatial = s, h       , hSpatial = s, w = -999, wSpatial    , phi = -999, phiSpatial    , s, habitatGrid)",
               "dcatState2Dead(z, gamma       , gammaSpatial = s, h = -999, hSpatial    , w = -999, wSpatial    , phi = -999, phiSpatial    , s, habitatGrid)",
               "dcatState2Dead(z, gamma       , gammaSpatial = s, h = -999, hSpatial    , w       , wSpatial = s, phi = -999, phiSpatial    , s, habitatGrid)",
               "dcatState2Dead(z, gamma = -999, gammaSpatial    , h       , hSpatial = s, w       , wSpatial = s, phi       , phiSpatial = s, s, habitatGrid)",
               "dcatState2Dead(z, gamma = -999, gammaSpatial    , h       , hSpatial = s, w = -999, wSpatial    , phi = -999, phiSpatial    , s, habitatGrid)",
               "dcatState2Dead(z, gamma = -999, gammaSpatial    , h = -999, hSpatial    , w = -999, wSpatial    , phi = -999, phiSpatial    , s, habitatGrid)",
               "dcatState2Dead(z, gamma = -999, gammaSpatial    , h = -999, hSpatial    , w       , wSpatial = s, phi = -999, phiSpatial    , s, habitatGrid)"
               ),
               
               
    types = c( "value = double(0)", "z = double(0)","gamma = double(0)", "gammaSpatial = double(1)", "h = double(0)", "hSpatial = double(1)" , "w = double(0)", "wSpatial = double(1)", "phi = double(0)", "phiSpatial = double(1)", "s = double(1)", "habitatGrid = double(2)"
    ),
    discrete = TRUE,
    mixedSizes = TRUE,
    pqAvail = FALSE
  )))


