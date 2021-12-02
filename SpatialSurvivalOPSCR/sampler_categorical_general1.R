sampler_categorical_general1 <- nimbleFunction(
  name = 'sampler_categorical1',
  contains = sampler_BASE,
  setup = function(model, mvSaved, target, control) {
    ## node list generation
    targetAsScalar <- model$expandNodeNames(target, returnScalarComponents = TRUE)
    calcNodes  <- model$getDependencies(target)
    calcNodesNoSelf <- model$getDependencies(target, self = FALSE)
    ## numeric value generation
    k <- control[["numCategories"]]
    if(is.null(k)) stop("Must provide control argument numCategories")
    probs <- numeric(k)
    logProbs <- numeric(k)
    ## checks
    if(length(targetAsScalar) > 1)  stop('cannot use categorical sampler on more than one target node')
   #[CM] if(model$getDistribution(target) != 'dCatZwh') stop('can only use categorical sampler on node with dcat distribution')
  },
  run = function() {
    currentValue <- model[[target]]
    logProbs[currentValue] <<- getLogProb(model, calcNodes)
    for(i in 1:k) {
      if(i != currentValue) {
        model[[target]] <<- i
        logProbPrior <- calculate(model, target)
        if(logProbPrior == -Inf) {
          logProbs[i] <<- -Inf
        } else {
          if(is.nan(logProbPrior)) {
            logProbs[i] <<- -Inf
          } else {
            logProbs[i] <<- logProbPrior + calculate(model, calcNodesNoSelf)
            if(is.nan(logProbs[i])) logProbs[i] <<- -Inf
          }
        }
      }
    }
    logProbs <<- logProbs - max(logProbs)
    probs <<- exp(logProbs)
    newValue <- rcat(1, probs) #rcat normalizes the probs internally
    if(newValue != currentValue) {
      model[[target]] <<- newValue
      calculate(model, calcNodes)
      nimCopy(from = model, to = mvSaved, row = 1, nodes = calcNodes, logProb = TRUE)
    } else nimCopy(from = mvSaved, to = model, row = 1, nodes = calcNodes, logProb = TRUE)
  },
  methods = list(
    reset = function() { }
  )
)


