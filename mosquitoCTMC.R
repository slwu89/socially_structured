rm(list=ls());gc()
# pmax(0,(b+cos(2*pi*(t-o*60)/24/60)))^p 


t=seq(from=0,to=24*60,by=0.01)

hourGrain = 120
xHours = which(t %% hourGrain == 0)
plot({
  (sin((2*pi*t)/1440)*0.5 + 1)
},type="l",ylim=c(0,2),ylab="Relative Activity Levels",xlab="Time of Day",main="Mosquito Diurnal Forcing",xaxt="n")
abline(h = 1,lty = 2)
axis(side = 1,at = xHours,labels = paste0(t[xHours]/60,":00"))


# t: tNow
# o: offset (in hours)
hazFunc <- function(t,o=-2){
  (sin((2*pi*(t-o*60))/1440)*0.75 + 1)
}
hazCurve = hazFunc(t = t,o = 0)
plot(hazCurve,type="l",ylim=c(0,2),ylab="Relative Activity Levels",xlab="Time of Day",main="Mosquito Diurnal Forcing",xaxt="n")
abline(h = 1,lty = 2)
axis(side = 1,at = xHours,labels = paste0(t[xHours]/60,":00"))

# hazFunc = function(t,o){
#   return(rep(1,length=length(t)))
# }

Btime = 1/1
Rtime = 1/1
Otime = 1/1
g = 1/10

exponentialMosquito <- function(o=0){
  
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
             tDur = rexp(n = 1,rate = Btime*force)
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
             tDur = rexp(n = 1,rate = Rtime*force)
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
             tDur = rexp(n = 1,rate = Otime*force)
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


mosyOut = exponentialMosquito()
plotMosyOut(mosyOut)


plotMosyOut = function(out,o=0){
  
  tMax = ceiling(max(out$times))
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
  
  # plot the mosquito time series
  pointsX = out$times*24*60
  pointsY = hazFunc(t = pointsX,o = 0)
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


cohort = parallel::mclapply(X = 1:100,FUN = function(x,o){exponentialMosquito(o=o)},o=0)



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
