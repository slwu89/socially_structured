dMin = 1:(60*24)

homeTime = 1/12
commute2Work = 1/1
commute2Home = 1/1
workTime = 1/10


simpleHuman <- function(maxSteps=1e3, N=1){
  states = vector(mode = "character",length = maxSteps+1)
  times = vector(mode = "numeric",length = maxSteps+1)
  duration = vector(mode = "numeric",length = maxSteps)

  times[1] = 0
  states[1] = "H"

  for(i in 1:maxSteps){
    print(paste0("iteration: ",i))
    switch(EXPR = states[i],
            H = {
              tDur = rgamma(n = 1,shape = N,rate = homeTime*N)
              times[i+1] = times[i] + tDur
              duration[i] = tDur
              states[i+1] = "H2W"
            },
            H2W = {
              tDur = rgamma(n = 1,shape = N,rate = commute2Work*N)
              times[i+1] = times[i] + tDur
              duration[i] = tDur
              states[i+1] = "W"
            },
            W = {
              tDur = rgamma(n = 1,shape = N,rate = workTime*N)
              times[i+1] = times[i] + tDur
              duration[i] = tDur
              states[i+1] = "W2H"
            },
            W2H = {
              tDur = rgamma(n = 1,shape = N,rate = commute2Home*N)
              times[i+1] = times[i] + tDur
              duration[i] = tDur
              states[i+1] = "H"
            },
            {stop(paste0("unrecognized state",states[i]))}
           )
  }
  return(list(states=states,times=times,duration=duration))
}


plotOneHuman <- function(out,lwd,cols){
  xlim = c(0,max(out$times)+1)
  ylim = c(0,0.4)

  plot(1,type="n",xaxt="n",yaxt="n",ylab="State Occupancy",xlab="Time (Days)",xlim=xlim,ylim=ylim,main="One Human (Assume Simulation begins at 20:00)\n Daily 20:00 in Red, Daily 08:00 in Blue")
  ttMax = max(out$times)/24
  axis(side = 1,at = (0:ttMax)*24,labels = 0:ttMax)
  axis(side = 2,at = c(0,0.125,0.25,0.375),labels = c("Home","H2W","Work","W2H"))
  abline(v =(0:ttMax)*24,col="red",lty=2 )
  abline(v = ((0+0.5):(ttMax+0.5))*24,col = "blue", lty=2)
  for(i in 1:length(out$states)){
    if(i==length(out$states)){
      switch(out$states[i],
             H = {
               segments(x0 = out$times[i],y0 = 0,x1 = out$times[i]+1,y1 = 0,col = cols[1],lwd = lwd,lend=2)
             },
             H2W = {
               segments(x0 = out$times[i],y0 = 0.125,x1 = out$times[i]+1,y1 = 0.125,col = cols[2],lwd = lwd,lend=2)
             },
             W = {
               segments(x0 = out$times[i],y0 = 0.25,x1 = out$times[i]+1,y1 = 0.25,col = cols[3],lwd = lwd,lend=2)
             },
             W2H = {
               segments(x0 = out$times[i],y0 = 0.375,x1 = out$times[i]+1,y1 = 0.375,col = cols[4],lwd = lwd,lend=2)
             }
      )
    } else {
      switch(out$states[i],
             H = {
               segments(x0 = out$times[i],y0 = 0,x1 = out$times[i+1],y1 = 0,col = cols[1],lwd = lwd,lend=2)
             },
             H2W = {
               segments(x0 = out$times[i],y0 = 0.125,x1 = out$times[i+1],y1 = 0.125,col = cols[2],lwd = lwd,lend=2)
             },
             W = {
               segments(x0 = out$times[i],y0 = 0.25,x1 = out$times[i+1],y1 = 0.25,col = cols[3],lwd = lwd,lend=2)
             },
             W2H = {
               segments(x0 = out$times[i],y0 = 0.375,x1 = out$times[i+1],y1 = 0.375,col = cols[4],lwd = lwd,lend=2)
             }
      )
    }
  }

}

cols = c("blue","grey50","purple","grey20")
out = simpleHuman(maxSteps = 25,N=24)
plotOneHuman(out = out,lwd = 15,cols = cols)

par(mfrow=c(1,2))
# exponential
set.seed(42)
out = simpleHuman(maxSteps = 50,N=1)
plotOneHuman(out = out,lwd = 15,cols = cols)

# gamma
out = simpleHuman(maxSteps = 50,N=12)
plotOneHuman(out = out,lwd = 15,cols = cols)
par(mfrow=c(1,1))
