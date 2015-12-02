## tutorial in the SYNCSA package by Valerio Pillar

# packages ----------------------------------------------------------------

library(SYNCSA)


# data --------------------------------------------------------------------

W <- read.table("data/W_14c_81spp.txt", nrows=81, header=TRUE) #read community data
B <- read.table("data/B_81spp_8t.txt", nrows=81, header=TRUE) #read trait data
E <- read.table("data/E_14c_1v.txt", nrows=14, header=TRUE) #read environmental data
DistPhyl <- read.table("data/81spp_PhylDissim.txt", nrows=81, header=TRUE) #read phylogenetic distance matrix

data <- organize.syncsa(W, B, DistPhyl, E)

set.seed(12345)
syncsa(comm = data$community,
       traits = data$traits,
       dist.spp = data$dist.spp,
       envir = data$environmental,
       method = "pearson", dist = "euclidean", scale = TRUE,
       scale.envir = TRUE, permutations = 999, na.rm = FALSE, notification = TRUE) 

#' # tcap
optimal(comm = data$community,
        traits = data$traits,
        envir = data$environmental,
        subset.min = 1, subset.max = 8, pattern = "tcap", dist = "euclidean",
        method = "pearson", scale = TRUE, scale.envir = TRUE, na.rm = FALSE,
        notification = TRUE, progressbar = TRUE)

#' maximum correlation is with traits `pi`, `he` and `ll`
#' So now we need to reduce the traits of the full matrix
#' 
good_traits <- colnames(data$traits) %in% c("pi", "he", "ll")
Btcap <- data$traits[,good_traits]

set.seed(12345)

syncsa(comm = data$community,
       traits = Btcap,
       dist.spp = data$dist.spp,
       envir = data$environmental,
       method = "pearson", dist = "euclidean", scale = TRUE,
       scale.envir = TRUE, permutations = 999, na.rm = FALSE, notification = TRUE)

#' # tdap
optimal(comm = data$community,
        traits = data$traits,
        envir = data$environmental,
        subset.min = 1, subset.max = 8, pattern = "tdap", dist = "euclidean", method = "pearson", scale = TRUE, scale.envir = TRUE, na.rm = FALSE, notification = TRUE, progressbar = TRUE)

#' # tcap and tdap
optimal(comm = data$community,
        traits = data$traits,
        envir = data$environmental,
        subset.min = 1, subset.max = 8, pattern = "tcap.tdap", dist = "euclidean", method = "pearson", scale = TRUE, scale.envir = TRUE, na.rm = FALSE, notification = TRUE, progressbar = TRUE)
