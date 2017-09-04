rm(list=ls());gc()


# follows notation from "Mt/G/inf Queues with Sinuoidal Arrival Rates (Eick, Massey, Whitt 1993)
# parameterized in days (hence psi=1 as denominator in sin curve)
sinRate <- function(t,alpha,lambda,off=0){
  if(abs(alpha)>1){stop("alpha parameter must be |alpha| <= 1")}
  lambda + lambda*alpha*sin((2*pi*t-(off/24))/1)
}

oneDay = seq(from=0,to=1,by=0.01)

lambdaMean = 1/3
sinCurve = sinRate(t=oneDay,alpha = 1,lambda = lambdaMean,off = 0)

plot({
  sinCurve
},type="l",xaxt="n",xlab = "Time of Day",ylab = "Sinusoidal Rate",main = "Periodic Transition Rate")
abline(h = lambdaMean,lty = 2,col="red")
abline(v = which.max(sinCurve),col = "purple")
text(x = which.max(sinCurve),y = lambdaMean+0.1,labels = paste0("max at ",oneDay[which.max(sinCurve)]*24),pos = 2,col="purple")
text(x = which.max(sinCurve),y = max(sinCurve),labels = paste0("max ",round(max(sinCurve),5)),pos = 1, col = "purple")
abline(v = which.min(sinCurve),col = "purple")
text(x = which.min(sinCurve),y = lambdaMean-0.1,labels = paste0("min at ",oneDay[which.min(sinCurve)]*24),pos = 2,col="purple")
text(x = which.min(sinCurve),y = min(sinCurve),labels = paste0("min ",round(min(sinCurve),5)),pos = 3, col = "purple")

inflexPts = rle(sign(diff(diff(sinCurve))))$lengths
for(i in 2:length(inflexPts)){
  pt = mean(rle(sign(diff(diff(sinCurve))))$lengths[i-1] + 1)
  abline(v = pt,col = "blue")
  text(x = pt,y = lambdaMean,pos = 1,col = "blue",labels = paste0("inflection point at ",oneDay[rle(sign(diff(diff(sinCurve))))$lengths[i-1]]))
}

axis(side = 1,at = which(oneDay%%0.1==0) ,labels = oneDay[oneDay%%0.1==0] * 24)


# gammaMosquito: a mosquito whose time-to-event (waiting time) distributions follow the Gamma(N,rate) distribution to reduce variance
gammaMosquito <- function(N=5,alpha=1,off=0){

  states = vector(mode = "list",length = 100)
  times = vector(mode = "list",length = 100)
  duration = vector(mode = "list",length = 100)

  times[[1]] = 0
  states[[1]] = "B"

  i = 2
  while(states[[i-1]]!="D"){

    print(paste0("iteration: ",i, " tNow: ",times[[i-1]]))
    switch(EXPR = states[[i-1]],
           B = {

             # diurnal forcing
             tNow = (times[[i-1]] %% 1)

             # duration of time until next event
             tDur = rgamma(n = 1,shape = N,rate = sinRate(t = tNow,alpha = alpha,lambda = Btime,off = off)*N)
             duration[[i-1]] = tDur

             # time of next event
             times[[i]] = times[[i-1]] + tDur

             # choose next event
             # pDie = pgamma(q = tDur,shape = N,rate = g*N)
             pDie = g*tDur
             if(runif(1) < pDie){
               states[[i]] = "D"
             } else {
               states[[i]] = "R"
             }

           },
           R = {

             # diurnal forcing
             tNow = (times[[i-1]] %% 1)

             # duration of time until next event
             tDur = rgamma(n = 1,shape = N,rate = sinRate(t = tNow,alpha = alpha,lambda = Rtime,off = off)*N)
             duration[[i-1]] = tDur

             # time of next event
             times[[i]] = times[[i-1]] + tDur

             # choose next event
             # pDie = pgamma(q = tDur,shape = N,rate = g*N)
             pDie = g*tDur
             if(runif(1) < pDie){
               states[[i]] = "D"
             } else {
               states[[i]] = "O"
             }

           },
           O = {

             # diurnal forcing
             tNow = (times[[i-1]] %% 1)

             # duration of time until next event
             tDur = rgamma(n = 1,shape = N,rate = sinRate(t = tNow,alpha = alpha,lambda = Otime,off = off)*N)
             duration[[i-1]] = tDur

             # time of next event
             times[[i]] = times[[i-1]] + tDur

             # choose next event
             # pDie = pgamma(q = tDur,shape = N,rate = g*N)
             pDie = g*tDur
             if(runif(1) < pDie){
               states[[i]] = "D"
             } else {
               states[[i]] = "B"
             }

           },
           {stop(paste0("unrecognized state",states[[i-1]]))}
    )

    i = i+1
  }
  return(list(
    states=unlist(Filter(Negate(is.null),states)),
    times=unlist(Filter(Negate(is.null),times)),
    duration=unlist(Filter(Negate(is.null),duration))
  ))
}


Btime = 1/1
Rtime = 1/1
Otime = 1/1
g = 1/12

mosyOut = gammaMosquito(N = 1,alpha = 1,off = 0)
plotMosyOut(mosyOut)

cohort = parallel::mclapply(X = 1:1e4,FUN = function(x){gammaMosquito(N=1,alpha=1)})



avgLife = vapply(X = cohort,FUN = function(x){
  tail(x$times,1)
},FUN.VALUE = numeric(1),USE.NAMES = FALSE)
mean(avgLife)
median(avgLife)




###############################################################################
# Old Hazard Forcing
###############################################################################








t=seq(from=0,to=24*60,by=0.01)

hourGrain = 120
xHours = which(t %% hourGrain == 0)
plot({
  (sin((2*pi*t)/1440) + 1)
},type="l",ylim=c(0,2),ylab="Relative Activity Levels",xlab="Time of Day",main="Mosquito Diurnal Forcing",xaxt="n")
abline(h = 1,lty = 2)
axis(side = 1,at = xHours,labels = paste0(t[xHours]/60,":00"))


# t: tNow
# o: offset (in hours)
hazFunc <- function(t,avgRate,offset=-2,scaling=0.25){
  (sin((2*pi*(t-offset*60))/1440)*scaling + 1) * avgRate
}
hazCurve = hazFunc(t = t,avgRate = 1/3,offset = 0)
plot(hazCurve,type="l",ylim=c(0,2),ylab="Relative Activity Levels",xlab="Time of Day",main="Mosquito Diurnal Forcing",xaxt="n")
abline(h = 1,lty = 2)
axis(side = 1,at = xHours,labels = paste0(t[xHours]/60,":00"))


Btime = 1/3
Rtime = 1/3
Otime = 1/3
g = 1/12

exponentialMosquito <- function(o=0,scaling=0.25){

  states = vector(mode = "list",length = 100)
  times = vector(mode = "list",length = 100)
  duration = vector(mode = "list",length = 100)

  times[[1]] = 0
  states[[1]] = "B"

  i = 2
  while(states[[i-1]]!="D"){

    print(paste0("iteration: ",i, " tNow: ",times[[i-1]]))
    switch(EXPR = states[[i-1]],
           B = {

             # diurnal forcing
             tNow = (times[[i-1]] %% 1)*24*60

             # duration of time until next event
             tDur = rexp(n = 1,rate = hazFunc(t = tNow,avgRate = Btime,offset = o,scaling = scaling))
             duration[[i-1]] = tDur

             # time of next event
             times[[i]] = times[[i-1]] + tDur

             # choose next event
             pDie = g*tDur
             if(runif(1) < pDie){
               states[[i]] = "D"
             } else {
               states[[i]] = "R"
             }

           },
           R = {

             # diurnal forcing
             tNow = (times[[i-1]] %% 1)*24*60

             # duration of time until next event
             tDur = rexp(n = 1,rate = hazFunc(t = tNow,avgRate = Rtime,offset = o,scaling = scaling))
             duration[[i-1]] = tDur

             # time of next event
             times[[i]] = times[[i-1]] + tDur

             # choose next event
             pDie = g*tDur
             if(runif(1) < pDie){
               states[[i]] = "D"
             } else {
               states[[i]] = "O"
             }

           },
           O = {

             # diurnal forcing
             tNow = (times[[i-1]] %% 1)*24*60

             # duration of time until next event
             tDur = rexp(n = 1,rate = hazFunc(t = tNow,avgRate = Otime,offset = o,scaling = scaling))
             duration[[i-1]] = tDur

             # time of next event
             times[[i]] = times[[i-1]] + tDur

             # choose next event
             pDie = g*tDur
             if(runif(1) < pDie){
               states[[i]] = "D"
             } else {
               states[[i]] = "B"
             }

           },
           {stop(paste0("unrecognized state",states[[i-1]]))}
    )

    i = i+1
  }
  return(list(
      states=unlist(Filter(Negate(is.null),states)),
      times=unlist(Filter(Negate(is.null),times)),
      duration=unlist(Filter(Negate(is.null),duration))
    ))
}


out = exponentialMosquito()
plotMosyOut(out,Btime,Rtime,Otime)


plotMosyOut = function(out,Btime,Rtime,Otime,scaling=0.25,o=0){

  tMax = ceiling(max(out$times))
  tSeq = seq(from=0,to=tMax*24*60,by=0.1)

  daySeq = vapply(X = 0:(tMax-1),FUN = function(x){x*24*60},FUN.VALUE = numeric(1),USE.NAMES = FALSE)

  # plot the diurnal forcing used
  hourGrain = 60*12
  xHours = which(tSeq %% hourGrain == 0)
  hazCurveB = hazFunc(t = tSeq,avgRate = Btime,offset = o,scaling = scaling)
  plot(x=tSeq,y=hazCurveB,type="l",ylim=c(0,2),ylab="Relative Activity Levels",xlab="Time of Day",main="Mosquito Diurnal Forcing",xaxt="n",col="red")
  abline(h = Btime,lty = 2,col="red")
  hazCurveR = hazFunc(t = tSeq,avgRate = Rtime,offset = o,scaling = scaling)
  lines(x=tSeq,y = hazCurveR,col="blue")
  abline(h = Rtime,lty =2,col="blue")
  hazCurveO = hazFunc(t = tSeq,avgRate = Otime,offset = o,scaling = scaling)
  lines(x=tSeq,y = hazCurveO,col="purple")
  abline(h=Otime,lty=2,col="purple")
  hourLabel = (tSeq[xHours]/60) %% 24
  dayLabel = paste0("Day",0:(tMax-1))
  axis(side = 1,at = tSeq[xHours],labels = paste0(hourLabel,":00"))
  text(x = daySeq,y=0.01,labels = dayLabel)

  # plot the mosquito time series
  pointsX = out$times*24*60
  pointsY = vector(mode = "numeric",length = length(pointsX))
  for(i in 1:length(pointsY)){
    pointsY[i] = switch(EXPR = out$states[i],
                          B = {
                              hazFunc(t = pointsX[i],avgRate = Btime,offset = o,scaling = scaling)
                            },
                          R = {
                              hazFunc(t = pointsX[i],avgRate = Rtime,offset = o,scaling = scaling)
                            },
                          O = {
                              hazFunc(t = pointsX[i],avgRate = Otime,offset = o,scaling = scaling)
                            },
                          D = {hazFunc(t = pointsX[i],avgRate = Otime,offset = o,scaling = scaling)}
                        )
  }

  # pointsY = hazFunc(t = pointsX,o = 0)
  pointsType = vapply(X = out$states,FUN = function(x){
      switch(x,
             B = 15,
             R = 16,
             O = 17,
             D = 18
       )
    },FUN.VALUE = numeric(1),USE.NAMES = FALSE)
  pointsCol = vapply(X = out$states,FUN = function(x){
    switch(x,
           B = "red",
           R = "blue",
           O = "purple",
           D = "darkgreen"
    )
  },FUN.VALUE = character(1),USE.NAMES = FALSE)
  points(pointsX,pointsY,pch=pointsType,col=pointsCol,cex=1.05)

}


cohort = parallel::mclapply(X = 1:1e4,FUN = function(x,o){exponentialMosquito(o=o)},o=0)



avgLife = vapply(X = cohort,FUN = function(x){
  tail(x$times,1)
},FUN.VALUE = numeric(1),USE.NAMES = FALSE)
avgLife = mean(avgLife); avgLife

ptsX = unlist(lapply(X = cohort,FUN = function(x){x$times}))
ptsY = hazFunc(t = ptsX*60*24,o = 0)
sum(ptsY>=1) / length(ptsY)
sum(ptsY<1) / length(ptsY)

plotMosyPopOut = function(cohort,o=0){

  tMax = vapply(X = cohort,FUN = function(x){
    tail(x$times,1)
  },FUN.VALUE = numeric(1),USE.NAMES = FALSE)
  tMax = ceiling(max(tMax))

  tSeq = seq(from=0,to=tMax*24*60,by=0.1)

  daySeq = vapply(X = 0:(tMax-1),FUN = function(x){x*24*60},FUN.VALUE = numeric(1),USE.NAMES = FALSE)

  # plot the diurnal forcing used
  hourGrain = 60*12
  xHours = which(tSeq %% hourGrain == 0)
  hazCurve = hazFunc(t = tSeq,o = o)
  plot(x=tSeq,y=hazCurve,type="l",ylim=c(0,2),ylab="Relative Activity Levels",xlab="Time of Day",main="Mosquito Diurnal Forcing",xaxt="n")
  abline(h = 1,lty = 2)
  hourLabel = (tSeq[xHours]/60) %% 24
  dayLabel = paste0("Day",0:(tMax-1))
  axis(side = 1,at = tSeq[xHours],labels = paste0(hourLabel,":00"))
  text(x = daySeq,y=0.01,labels = dayLabel)


  for(i in 1:length(cohort)){
    # plot the mosquito time series
    pointsX = cohort[[i]]$times*24*60
    pointsY = hazFunc(t = pointsX,o = 0)
    pointsType = vapply(X = cohort[[i]]$states,FUN = function(x){
      switch(x,
             B = 15,
             R = 16,
             O = 17,
             D = 18
      )
    },FUN.VALUE = numeric(1),USE.NAMES = FALSE)
    pointsCol = vapply(X = cohort[[i]]$states,FUN = function(x){
      switch(x,
             B = "red",
             R = "blue",
             O = "purple",
             D = "darkgreen"
      )
    },FUN.VALUE = character(1),USE.NAMES = FALSE)
    points(pointsX,pointsY,pch=pointsType,col=pointsCol,cex=1.05)
  }

}

plotMosyPopOut(cohort)





gammaMosquito <- function(N=5,o=0){

  states = vector(mode = "list",length = 100)
  times = vector(mode = "list",length = 100)
  duration = vector(mode = "list",length = 100)

  times[[1]] = 0
  states[[1]] = "B"

  i = 2
  while(states[[i-1]]!="D"){

    print(paste0("iteration: ",i, " tNow: ",times[[i-1]]))
    switch(EXPR = states[[i-1]],
           B = {

             # diurnal forcing
             tNow = (times[[i-1]] %% 1)*24*60
             force = hazFunc(t=tNow,o = o)

             # duration of time until next event
             tDur = rgamma(n = 1,shape = N,rate = Btime*force*N)
             duration[[i-1]] = tDur

             # time of next event
             times[[i]] = times[[i-1]] + tDur

             # choose next event
             pDie = g*tDur
             if(runif(1) < pDie){
               states[[i]] = "D"
             } else {
               states[[i]] = "R"
             }

           },
           R = {

             # diurnal forcing
             tNow = (times[[i-1]] %% 1)*24*60
             force = hazFunc(t=tNow,o = o)

             # duration of time until next event
             tDur = rgamma(n = 1,shape = N,rate = Rtime*force*N)
             duration[[i-1]] = tDur

             # time of next event
             times[[i]] = times[[i-1]] + tDur

             # choose next event
             pDie = g*tDur
             if(runif(1) < pDie){
               states[[i]] = "D"
             } else {
               states[[i]] = "O"
             }

           },
           O = {

             # diurnal forcing
             tNow = (times[[i-1]] %% 1)*24*60
             force = hazFunc(t=tNow,o = o)

             # duration of time until next event
             tDur = rgamma(n = 1,shape = N,rate = Otime*force*N)
             duration[[i-1]] = tDur

             # time of next event
             times[[i]] = times[[i-1]] + tDur

             # choose next event
             pDie = g*tDur
             if(runif(1) < pDie){
               states[[i]] = "D"
             } else {
               states[[i]] = "B"
             }

           },
           {stop(paste0("unrecognized state",states[[i-1]]))}
    )

    i = i+1
  }
  return(list(
    states=unlist(Filter(Negate(is.null),states)),
    times=unlist(Filter(Negate(is.null),times)),
    duration=unlist(Filter(Negate(is.null),duration))
  ))
}


cohort = parallel::mclapply(X = 1:100,FUN = function(x,N,o){gammaMosquito(N=N,o=o)},o=0,N=12)

plotMosyPopOut(cohort = cohort)



out = gammaMosquito(N = 64,o = 0)
plotMosyOut(out)
