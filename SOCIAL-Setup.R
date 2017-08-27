###############################################################################
#
#    _____   _ __  _    ___ _____       _   __
#   /__  /  (_) /_| |  / (_) ___/____  / | / /
#     / /  / / //_/ | / / /\__ \/ __ \/  |/ /
#    / /__/ / ,<  | |/ / /___/ / /_/ / /|  /
#   /____/_/_/|_| |___/_//____/\____/_/ |_/
#
#   ZikViSoN: Zika Virus Elimination through the use of Social Networks
#   Simulation Initialization Functions
#   Héctor M. Sánchez C., Edgar E. Vallejo, Sean L. Wu
#   Model Adapted from https://www.ncbi.nlm.nih.gov/pubmed/24593919
#   August 26, 2017
#
###############################################################################

#' Generate Parameters for Simulation
#'
#' Parameter definitions below.
#'
#' @param IN governs the number of houses an individual moves to from their social group (including their home); in general, kappa = IN - 1
#' @param OUT governs the number of random houses an individual goes to
#' @param foray will indicate what percent of the population goes to one location outside their social group when IN=1; kappa = IN - 1  + (1 - foray)
#' @param spread governs what percent of the risk from infectious mosquitoes leaves the home
#' @param CityWidth total size of 'city' is CityWidth X CityLength
#' @param CityLength total size of 'city' is CityWidth X CityLength
#' @param AvgPerHouse average household size
#' @param PercentHomeMin minimum amount of time at home
#' @param PercentHomeMax maximum amount of time at home
#' @param
#' @param
#' @param
#' @param
#' @param
#'
#' @section Derived Parameters
#'
SOCIAL.Parameters <- function(
    IN = 5,
    OUT = 1,
    foray = 50,
    spread = 0.2,
    CityWidth = 20,
    CityLength = 20,
    AvgPerHouse = 6,
    PercentHomeMin = 0.4,
    PercentHomeMax = 0.75
  ){
    par = as.list(match.call())[-1]

    par$MinInSocial = IN # Home counts as an element of the Social
    MaxInSocial=IN
    MinOutOfSocial=0
    MaxOutOfSocial=OUT
    ### Sets the percent of people who leave their social group
    DistOutOfSocial=c(1-(0.01*Forray),0.01*Forray)
    SmallestSocial=max(MinInSocial,1)
    LargestSocial=IN
    if (MaxInSocial==1)LargestSocial=1


    return(par)
  }






  MinInSocial=IN # Home counts as an element of the Social
  MaxInSocial=IN
  MinOutOfSocial=0
  MaxOutOfSocial=OUT
  ### Sets the percent of people who leave their social group
  DistOutOfSocial=c(1-(0.01*Forray),0.01*Forray)
  SmallestSocial=max(MinInSocial,1)
  LargestSocial=IN
  if (MaxInSocial==1)LargestSocial=1
