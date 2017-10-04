
dMin = 1:(60*24)
d2Min = c(1:(60*48))

cosHaz = function(t, o=-2, b=.5, p=3){
  #NOTE :: offset o is in hours
  pmax(0,(b+cos(2*pi*(t-o*60)/24/60)))^p 
}

cumCosHaz = function(o=-2,b=0.5,p=2){
  #NOTE :: offset o is in hours
  #cs = c(0,cumsum(cosHaz(d2Min[-1],k,o,b,p)))
  cs = cumsum(cosHaz(d2Min,o,b,p))
  cs1 = 2*cs/max(cs)
  cs1
}

rCosHazTim = function(t, o=-2, b=1.1, p=2, cs2=NULL){
  t0 = ceiling((t-floor(t))*60*24)+1
  if(is.null(cs2)){
    cs2=cumCosHaz(o,b,p)
  }
  css2 = cs2[t0:c(t0+1440)]-cs2[t0]
  tm = min(which(css2>runif(1,0,1))) 
  # This occationally returns a number that is slightly larger than 1
  # which should never occur. 
  # print(c(check = diff(range(css2)))) 
  min(tm/60/24,1)
}

rCosHaz = function(t,P=1,o=-2,b=1,p=2,cs2=NULL){
  t+rCosHazTim(t,o,b,p,cs2)+rgeom(1,P)
}


# parameters
# t: hour of day / 24
# P: proportion of the total hazard (of the event) taking place in the next 24 hours
# o: offset/timing of maximum activity
# b: some sort of scaling parameter
# p: power
# N: number of humans to simulate
t = 6/24
P = 0.95
o = 12
b = 1.1
p = 2
N = 1e4

par(mfrow = c(2,2))
x <- c(dMin/60,dMin/60+24) # 2 days given in units of hours
y <- rep(cosHaz(dMin,o,b,p),2) # back to back cosHaz function
startx <- which.min((x-t*24)^2) # used for polygon 
endx <- startx + 24*60	 # used for polygon
plot(x,y, type = "l", ylab = "Relative Activity Levels", xaxt = "n", xlab = "Time of Day")
polygon(c(x[startx:endx],x[endx:startx]),c(y[startx:endx],rep(0,length(startx:endx))),col=rgb(0,0,0,.1),border=NA)
abline(v=x[startx],lty=3)
abline(v=x[endx],lty=3)
# lines(rep(o,2),c(0,cosHaz(dMin,o,b,p)[which(dMin/60==o)]),lty=2)
lines(rep(o,2),c(0,cosHaz(o*60,o,b,p)),lty=2)
if (o<12){
  text(o+.1,.1,paste("Max activity
                     occurs at",o),pos=4)
} else {
  text(o+.1,.1,paste("Max activity
                     occurs at",o),pos=2)
}

axis(1, 4*c(0:12), 4*c(0:12))

t0 = ceiling((t-floor(t))*60*24)+1
ixx = t0:(t0+1440)

cch = P*(cumCosHaz(o,b,p)[ixx] - cumCosHaz(o,b,p)[t0]) 
x = d2Min[ixx]/60
y = cch
plot(x,y, type = "l", ylab = "CDF of Activity (One Day)", xlab = "Time of Day", xaxt = "n", ylim = c(0,1)) 
polygon(c(x,rev(x)),c(y,rep(0,length(y))),col=rgb(0,0,0,.1),border=NA)
lines(c(0,max(d2Min[ixx]/60)),c(P,P),lty=2)
if (P>.5){
  text(min(d2Min[ixx]/60)+4,P,paste("Proportion of hazard
                                    in next 24 hours:",P),pos=1)
} else {
  text(min(d2Min[ixx]/60)+4,P,paste("Proportion of hazard
                                    in next 24 hours:",P),pos=3)
}
axis(1, 4*c(0:11), 4*c(0:11))

tt = replicate(N,rCosHazTim(1,o,b,p))
hist(tt, main = "Time of Day")

tt=replicate(N,rCosHaz(t,P,o,b,p))
hist(tt, main = paste("Time to Event, P=",P), xlab = "Time (Days)")