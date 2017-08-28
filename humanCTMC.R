dMin = 1:(60*24)

homeTime = 1/12
commute2Work = 1/1
commute2Home = 1/1
workTime = 1/10



simpleHuman <- function(maxSteps=1e3){
  states = vector(mode = "character",length = maxSteps+1)
  times = vector(mode = "numeric",length = maxSteps+1)
  duration = vector(mode = "numeric",length = maxSteps+1)
  
  times[1] = 0
  states[1] = "H"
  duration[1] = rexp(n = 1,rate = homeTime)
  
  for(i in 1:maxSteps){
    print(paste0("iteration: ",i))
    switch(EXPR = states[i],
            H = {
              tDur = rexp(n = 1,rate = homeTime)
              times[i+1] = times[i] + tDur
              duration[i+1] = tDur
              states[i+1] = "H2W"
            },
            H2W = {
              tDur = rexp(n = 1,rate = commute2Work)
              times[i+1] = times[i] + tDur
              duration[i+1] = tDur
              states[i+1] = "W"
            },
            W = {
              tDur = rexp(n = 1,rate = workTime)
              times[i+1] = times[i] + tDur
              duration[i+1] = tDur
              states[i+1] = "W2H"
            },
            W2H = {
              tDur = rexp(n = 1,rate = commute2Home)
              times[i+1] = times[i] + tDur
              duration[i+1] = tDur
              states[i+1] = "H"
            },
            {stop(paste0("unrecognized state",states[i]))}
           )
  }
  return(list(states=states,times=times,duration=duration))
}

out = simpleHuman(maxSteps = 50)

xlim = c(0,max(out$times)+1)
ylim = c(0,5)

plot(1,type="n",xaxt="n",yaxt="n",ylab="State Occupancy",xlab="Time (Days)",xlim=xlim,ylim=ylim)
ttMax = max(out$times)/24
axis(side = 1,at = (0:ttMax)*24,labels = 0:ttMax)
cols = c("blue","grey50","purple","grey20")
for(i in 1:length(out$states)){
  if(i==length(out$states)){
    switch(out$states[i],
           H = {
             segments(x0 = out$times[i],y0 = 1,x1 = out$times[i]+1,y1 = 1,col = cols[1],lwd = 3)
           },
           H2W = {
             segments(x0 = out$times[i],y0 = 2,x1 = out$times[i]+1,y1 = 2,col = cols[2],lwd = 3)
           },
           W = {
             segments(x0 = out$times[i],y0 = 3,x1 = out$times[i]+1,y1 = 3,col = cols[3],lwd = 3)
           },
           W2H = {
             segments(x0 = out$times[i],y0 = 4,x1 = out$times[i]+1,y1 = 4,col = cols[4],lwd = 3)
           }
    )
  } else {
    switch(out$states[i],
            H = {
              segments(x0 = out$times[i],y0 = 1,x1 = out$times[i+1],y1 = 1,col = cols[1],lwd = 3)
            },
            H2W = {
              segments(x0 = out$times[i],y0 = 2,x1 = out$times[i+1],y1 = 2,col = cols[2],lwd = 3)
            },
            W = {
              segments(x0 = out$times[i],y0 = 3,x1 = out$times[i+1],y1 = 3,col = cols[3],lwd = 3)
            },
            W2H = {
              segments(x0 = out$times[i],y0 = 4,x1 = out$times[i+1],y1 = 4,col = cols[4],lwd = 3)
            }
           )
  }
}

t=seq(from=0,to=24*60,by=0.01)
plot(sin((2*pi*t)/1440),type="l")
plot({
    sin((2*pi*t)/1440) + cos((2*pi*t)/1440)
  },type = "l")

hazFunc = function(t){
  sin((2*pi*t)/1440) + cos((2*pi*t)/1440)
}

out = sin((2*pi*t)/1440) + cos((2*pi*t)/1440)
integrate(f = hazFunc,lower = 0,upper = 1440,subdivisions = 1e3)




hazardFunction <- function(t,sinK,cosK,pow=1){
  sin((sinK*pi*t)/1440) + cos((cosK*pi*t)/1440)
}

tseq = 0:(60*24)

plotHazard <- function(sinK,cosK,pow=1){
  plot(x = tseq,y = hazardFunction(tseq,sinK,cosK,pow),type="l",xlab="Time of Day",ylab="Hazard")
}
plotHazard(sinK = 2,cosK = 2)

