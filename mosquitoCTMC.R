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
  (sin((2*pi*(t-o*60))/1440)*0.5 + 1)
}
hazCurve = hazFunc(t = t,o = 0)
plot(hazCurve,type="l",ylim=c(0,2),ylab="Relative Activity Levels",xlab="Time of Day",main="Mosquito Diurnal Forcing",xaxt="n")
abline(h = 1,lty = 2)
axis(side = 1,at = xHours,labels = paste0(t[xHours]/60,":00"))



Btime = 1/1
Rtime = 1/1
Otime = 1/1
g = 1/10

simpleMosquito <- function(o=0){
  
  states = vector(mode = "list",length = 100)
  times = vector(mode = "list",length = 100)
  duration = vector(mode = "list",length = 100)
  
  times[[1]] = 0
  states[[1]] = "B"
  duration[[1]] = rexp(n = 1,rate = Btime)
  
  i = 2
  while(states[[i-1]]!="D"){
    
    print(paste0("iteration: ",i, " tNow: ",times[[i-1]]))
    switch(EXPR = states[[i-1]],
           B = {
             
             # diurnal forcing
             tNow = (times[[i-1]] %% 1)*24*60
             force = hazFunc(t=tNow,o = o)
             
             # duration of time until next event
             tDur = rexp(n = 1,rate = Btime/force)
             duration[[i]] = tDur
             
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
             tDur = rexp(n = 1,rate = Rtime/force)
             duration[[i]] = tDur
             
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
             tDur = rexp(n = 1,rate = Otime/force)
             duration[[i]] = tDur
             
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


mosyOut = simpleMosquito()

out = mosyOut
plotMosyOut = function(out,o=0){
  
  tMax = ceiling(max(out$times))
  tSeq = seq(from=0,to=tMax*24*60,by=0.1)
  
  # plot the diurnal forcing used
  hourGrain = 120
  xHours = which(tSeq %% hourGrain == 0)
  hazCurve = hazFunc(t = tSeq,o = 0)
  plot(hazCurve,type="l",ylim=c(0,2),ylab="Relative Activity Levels",xlab="Time of Day",main="Mosquito Diurnal Forcing",xaxt="n")
  abline(h = 1,lty = 2)
  # fix so that it prints just hours, not cumulative hours
  # axis(side = 1,at = xHours,labels = paste0(tSeq[xHours]/60,":00"))
  
}


# plot({
#   ((sin((2*pi*t)/1440) + cos((2*pi*t)/1440)) ) + 1.5
# },type = "l",ylab = "Amplitude")
# abline(h = 1,lty=2)
# 
# plot({
#   ((sin((2*pi*t)/1440) + cos((2*pi*t)/1440)) * 0.5)
# },type = "l",ylab = "Amplitude")
# abline(h = 0,lty=2)
# 
# hazFunc = function(t){
#   sin((2*pi*t)/1440) + cos((2*pi*t)/1440)
# }
# 
# out = sin((2*pi*t)/1440) + cos((2*pi*t)/1440)
# integrate(f = hazFunc,lower = 0,upper = 1440,subdivisions = 1e3)