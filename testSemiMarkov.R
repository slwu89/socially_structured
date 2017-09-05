
lifespan = 12
waitTime = 2


semiMarkovOrganism <- function(N,waitTime,lifespan){
  
  states = vector(mode = "list",length = 100)
  times = vector(mode = "list",length = 100)
  duration = vector(mode = "list",length = 100)
  
  times[[1]] = 0
  states[[1]] = 1
  
  i=2
  while(states[[i-1]]!="D"){
    tau = rgamma(n = 1,shape = N,rate = (1/waitTime)*N)
    states[[i]] = states[[i-1]] + 1
    duration[[i-1]] = tau
    times[[i]] = times[[i-1]] + tau
    if(runif(1) < pexp(q=tau,rate = -log(1-(1/lifespan)))){
      states[[i]] = "D"
    }
    i=i+1
  }
  return(list(
    states=Filter(f = Negate(is.null),x = states),
    times=Filter(f = Negate(is.null),x = times),
    duration=Filter(f = Negate(is.null),x = duration)
  ))
}

xx = semiMarkovOrganism(N = 1,waitTime = waitTime,lifespan = lifespan)

yy= parallel::mclapply(X = 1:1e4,FUN = function(x,N,waitTime,lifespan){
  semiMarkovOrganism(N = N,waitTime = waitTime,lifespan = lifespan)
},N=1,waitTime=waitTime,lifespan=lifespan)

life = unlist(lapply(X = yy,FUN = function(x){tail(x$times,1)}))

mean(life)
median(life)
var(life)

dwell = unlist(lapply(X = yy,FUN = function(x){x$duration}))

mean(dwell)
median(dwell)
var(dwell)
