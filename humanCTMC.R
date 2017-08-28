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

out = simpleHuman(maxSteps = 1000)



t=seq(from=0,to=24*60,by=0.01)
plot(sin((2*pi*t)/1440),type="l")
plot({sin((2*pi*t)/1440) + cos((2*pi*t)/1440)},type = "l")

hazFunc = function(t){
  sin((2*pi*t)/1440) + cos((2*pi*t)/1440)
}

out = sin((2*pi*t)/1440) + cos((2*pi*t)/1440)
integrate(f = hazFunc,lower = 0,upper = 1440,subdivisions = 1e3)
