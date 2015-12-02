W <- read.table("W_14c_81spp.txt", nrows=81, header=TRUE) #read community data
B <- read.table("B_81spp_8t.txt", nrows=81, header=TRUE) #read trait data
E <- read.table("E_14c_1v.txt", nrows=14, header=TRUE) #read environmental data
DistPhyl <- read.table("81spp_PhylDissim.txt", nrows=81, header=TRUE) #read phylogenetic distance matrix

organize.syncsa(W, B, DistPhyl, E)
syncsa (W, B, DistPhyl, E, method = "pearson", dist = "euclidean", scale = TRUE, scale.envir = TRUE, permutations = 999, na.rm = FALSE, notification = TRUE) 
syncsa (W, B, DistPhyl, E, method = "pearson", dist = "euclidean", scale = TRUE, scale.envir = FALSE, permutations = 999, na.rm = FALSE, notification = TRUE) 

optimal(W, E, B, subset.min = 1, subset.max = 8, pattern = "tcap", dist = "euclidean", method = "pearson", scale = TRUE, scale.envir = TRUE, na.rm = FALSE, notification = TRUE, progressbar = TRUE)
optimal(W, E, B, subset.min = 1, subset.max = 8, pattern = "tdap", dist = "euclidean", method = "pearson", scale = TRUE, scale.envir = TRUE, na.rm = FALSE, notification = TRUE, progressbar = TRUE)
optimal(W, E, B, subset.min = 1, subset.max = 8, pattern = "tcap.tdap", dist = "euclidean", method = "pearson", scale = TRUE, scale.envir = TRUE, na.rm = FALSE, notification = TRUE, progressbar = TRUE)

