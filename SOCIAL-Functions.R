###############################################################################
#
#    _____   _ __  _    ___ _____       _   __
#   /__  /  (_) /_| |  / (_) ___/____  / | / /
#     / /  / / //_/ | / / /\__ \/ __ \/  |/ /
#    / /__/ / ,<  | |/ / /___/ / /_/ / /|  /
#   /____/_/_/|_| |___/_//____/\____/_/ |_/
#
#   ZikViSoN: Zika Virus Elimination through the use of Social Networks
#   Function Definitions
#   Héctor M. Sánchez C., Edgar E. Vallejo, Sean L. Wu
#   Model Adapted from https://www.ncbi.nlm.nih.gov/pubmed/24593919
#   August 26, 2017
#
###############################################################################

Neighbors <- function(loc){
	locx=floor((loc-1)/CityLength)+1
	locy=(loc-1)%%CityLength+1
	Nei=c(locx,locy)
	for (i in -1:1){
		for (j in -1:1){
			Nei=rbind(Nei,c(locx+i,locy+j))
		}
	}
	Nei=Nei[-1,]
	for (i in length(Nei[,1]):1){
		if (Nei[i,1]<1) Nei=Nei[-i,]
		else if (Nei[i,2]<1) Nei=Nei[-i,]
		else if (Nei[i,1]>CityWidth) Nei=Nei[-i,]
		else if (Nei[i,2]>CityLength) Nei=Nei[-i,]
	}
	NumNei=Nei[,2]+CityLength*(Nei[,1]-1)
	return(NumNei)
}


Neighbors2 <- function(loc){
	locx=floor((loc-1)/CityLength)+1
	locy=(loc-1)%%CityLength+1
	Nei=c(locx,locy)
	for (i in -2:2){
		for (j in -2:2){
			Nei=rbind(Nei,c(locx+i,locy+j))
		}
	}
	Nei=Nei[-1,]
	for (i in length(Nei[,1]):1){
		if (Nei[i,1]<1) Nei=Nei[-i,]
		else if (Nei[i,2]<1) Nei=Nei[-i,]
		else if (Nei[i,1]>CityWidth) Nei=Nei[-i,]
		else if (Nei[i,2]>CityLength) Nei=Nei[-i,]
	}
	NumNei=Nei[,2]+CityLength*(Nei[,1]-1)
	return(NumNei)
}

HousesleNAway <- function(loc,distance){
	if(distance==0){
    return(loc)
  } else {
	locx=floor((loc-1)/CityLength)+1
	locy=(loc-1)%%CityLength+1
	Nei=c(locx,locy)
	for (i in -distance:distance){
		for (j in -distance:distance){
			Nei=rbind(Nei,c(locx+i,locy+j))
		}
	}
	Nei=Nei[-1,]
	for (i in length(Nei[,1]):1){
		if (Nei[i,1]<1) Nei=Nei[-i,]
		else if (Nei[i,2]<1) Nei=Nei[-i,]
		else if (Nei[i,1]>CityWidth) Nei=Nei[-i,]
		else if (Nei[i,2]>CityLength) Nei=Nei[-i,]
	}
	NumNei=Nei[,2]+CityLength*(Nei[,1]-1)
	return(NumNei)
	}
}

HousesOnlyNAway<-function(loc,distance){
	if (distance==0) return(loc)
	else {
		return(setdiff(HousesleNAway(loc,distance),HousesleNAway(loc,distance-1)))
	}
}
