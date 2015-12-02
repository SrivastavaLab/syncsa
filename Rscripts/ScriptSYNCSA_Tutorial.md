# tutorial in the SYNCSA package by Valerio Pillar


```r
## 

# packages ----------------------------------------------------------------

suppressPackageStartupMessages(library(SYNCSA))
```

```
## Warning: package 'SYNCSA' was built under R version 3.1.3
```

```
## Warning: package 'vegan' was built under R version 3.1.3
```

```
## Warning: package 'permute' was built under R version 3.1.3
```

```
## Warning: package 'mice' was built under R version 3.1.3
```

```
## Warning: package 'Rcpp' was built under R version 3.1.3
```

```
## Warning: package 'FD' was built under R version 3.1.3
```

```
## Warning: package 'ade4' was built under R version 3.1.3
```

```
## Warning: package 'ape' was built under R version 3.1.3
```

```
## Warning: package 'geometry' was built under R version 3.1.3
```

```
## Warning: package 'magic' was built under R version 3.1.3
```

```
## Warning: package 'abind' was built under R version 3.1.3
```

```r
# data --------------------------------------------------------------------

W <- read.table("../data/W_14c_81spp.txt", nrows=81, header=TRUE) #read community data
B <- read.table("../data/B_81spp_8t.txt", nrows=81, header=TRUE) #read trait data
E <- read.table("../data/E_14c_1v.txt", nrows=14, header=TRUE) #read environmental data
DistPhyl <- read.table("../data/81spp_PhylDissim.txt", nrows=81, header=TRUE) #read phylogenetic distance matrix

data <- organize.syncsa(W, B, DistPhyl, E)

set.seed(12345)
syncsa(comm = data$community,
       traits = data$traits,
       dist.spp = data$dist.spp,
       envir = data$environmental,
       method = "pearson", dist = "euclidean", scale = TRUE,
       scale.envir = TRUE, permutations = 999, na.rm = FALSE, notification = TRUE) 
```

```
## $Notes
##             Correlation meanings                                                                                             
## note.roTE   "Trait-convergence assembly patterns (TCAP): roTE"                                                               
## note.roXE   "Both trait-convergence assembly patterns and trait-divergence assembly patterns: roXE"                          
## note.roXE.T "Trait-divergence assembly patterns (TDAP): roXE.T"                                                              
## note.roBF   "Phylogenetic signal at species level: roBF"                                                                     
## note.roPE   "Correlation of phylogenetically structured assembly patterns to ecological variables: roPE"                     
## note.roPT   "Correlation of phylogenetically structured assembly patterns to trait-convergence assembly patterns: roPT"      
## note.roPX.T "Correlation of phylogenetically structured assembly patterns to trait-divergence assembly patterns: roPX.T"     
## note.roTE.P "Removing phylogeny from trait-convergence assembly patterns: roTE.P"                                            
## note.roXE.P "Removing phylogeny from both trait-convergence assembly patterns and trait-divergence assembly patterns: roXE.P"
## 
## $Statistics
##        Obs         p    
## roTE   0.2525317   0.2  
## roXE   0.166629    0.307
## roPE   0.2259687   0.266
## roPT   0.5486798   0.025
## roPX.T 0.2323918   0.188
## roXE.T -0.09442922 0.734
## roTE.P 0.1578413   0.289
## roXE.P 0.04755375  0.439
## roBF   0.2303782   0.001
```

# tcap


```r
optimal(comm = data$community,
        traits = data$traits,
        envir = data$environmental,
        subset.min = 1, subset.max = 8, pattern = "tcap", dist = "euclidean",
        method = "pearson", scale = TRUE, scale.envir = TRUE, na.rm = FALSE,
        notification = TRUE, progressbar = FALSE)
```

```
##                      Subset           ro
## 48                 pi he ll  0.614129188
## 121             pi he la sh  0.602718520
## 50                 pi he sh  0.580111664
## 176          pi tx he la sh  0.579174114
## 108             pi ts he ll  0.577074755
## 98              pi tx he ll  0.576025247
## 110             pi ts he sh  0.569604123
## 119             pi he ll sh  0.568838399
## 11                    pi he  0.568751250
## 193          pi he ll la sh  0.568560310
## 184          pi ts he ll sh  0.564340700
## 12                    pi ll  0.562265183
## 186          pi ts he la sh  0.562142955
## 118             pi he ll la  0.560402811
## 229       pi tx he ll la sh  0.554225178
## 163          pi tx ts he ll  0.549636896
## 174          pi tx he ll sh  0.549552237
## 234       pi ts he ll la sh  0.548772352
## 100             pi tx he sh  0.548769100
## 220       pi tx ts he ll sh  0.548027394
## 222       pi tx ts he la sh  0.543313706
## 165          pi tx ts he sh  0.542605511
## 173          pi tx he ll la  0.540220658
## 247    pi tx ts he ll la sh  0.536812172
## 55                 pi la sh  0.519305492
## 53                 pi ll sh  0.516860362
## 183          pi ts he ll la  0.515223324
## 14                    pi sh  0.509044593
## 124             pi ll la sh  0.507191836
## 49                 pi he la  0.502685021
## 113             pi ts ll sh  0.500230596
## 39                 pi tx ll  0.498875804
## 219       pi tx ts he ll la  0.498716940
## 43                 pi ts he  0.498476630
## 44                 pi ts ll  0.494439971
## 38                 pi tx he  0.492068566
## 103             pi tx ll sh  0.489677926
## 46                 pi ts sh  0.485475302
## 179          pi tx ll la sh  0.485126822
## 105             pi tx la sh  0.483907627
## 189          pi ts ll la sh  0.480221530
## 52                 pi ll la  0.479456143
## 168          pi tx ts ll sh  0.477460183
## 115             pi ts la sh  0.471055297
## 41                 pi tx sh  0.467782478
## 99              pi tx he la  0.463800543
## 225       pi tx ts ll la sh  0.461524165
## 75                 ts he sh  0.450607899
## 96              pi tx ts sh  0.450157312
## 94              pi tx ts ll  0.450018391
## 93              pi tx ts he  0.447916650
## 170          pi tx ts la sh  0.442680013
## 102             pi tx ll la  0.441107129
## 149             ts he ll sh  0.439786312
## 109             pi ts he la  0.435321424
## 199          tx ts he ll sh  0.426983657
## 130             tx ts he sh  0.425162215
## 112             pi ts ll la  0.425143191
## 29                    he sh  0.416870047
## 151             ts he la sh  0.416197452
## 213          ts he ll la sh  0.412208205
## 84                 he ll sh  0.409541857
## 86                 he la sh  0.408326459
## 73                 ts he ll  0.407710950
## 164          pi tx ts he la  0.406279438
## 240       tx ts he ll la sh  0.401374804
## 139             tx he ll sh  0.400224554
## 65                 tx he sh  0.399155637
## 201          tx ts he la sh  0.397805114
## 158             he ll la sh  0.396516075
## 167          pi tx ts ll la  0.394479013
## 141             tx he la sh  0.391241794
## 128             tx ts he ll  0.390342434
## 1                        pi  0.388934912
## 208          tx he ll la sh  0.387103194
## 78                 ts ll sh  0.366038484
## 27                    he ll  0.356182831
## 63                 tx he ll  0.355757833
## 25                    ts sh  0.350997661
## 32                    ll sh  0.349042783
## 133             tx ts ll sh  0.347587115
## 13                    pi la  0.344133897
## 22                    ts he  0.342760883
## 4                        he  0.340164475
## 154             ts ll la sh  0.339061046
## 148             ts he ll la  0.335492106
## 7                        sh  0.334416505
## 68                 tx ll sh  0.331982302
## 10                    pi ts  0.330822620
## 89                 ll la sh  0.330188430
## 204          tx ts ll la sh  0.321149754
## 198          tx ts he ll la  0.319590115
## 34                    la sh  0.317371383
## 80                 ts la sh  0.316943425
## 61                 tx ts sh  0.316869485
## 144             tx ll la sh  0.312687933
## 23                    ts ll  0.312139878
## 83                 he ll la  0.311647902
## 20                    tx sh  0.304824913
## 138             tx he ll la  0.301159674
## 9                     pi tx  0.298354568
## 45                 pi ts la  0.292074541
## 58                 tx ts he  0.288576159
## 135             tx ts la sh  0.287350597
## 5                        ll  0.286364328
## 70                 tx la sh  0.286229788
## 17                    tx he  0.285161881
## 40                 pi tx la  0.284890566
## 37                 pi tx ts  0.269050403
## 59                 tx ts ll  0.269011704
## 236       pi ts he ll sh rn  0.267266240
## 253    pi ts he ll la sh rn  0.264389978
## 239       pi he ll la sh rn  0.257855082
## 195          pi he ll sh rn  0.257054519
## 249    pi tx ts he ll sh rn  0.254602083
## 18                    tx ll  0.254043890
## 255 pi tx ts he ll la sh rn  0.252531735
## 95              pi tx ts la  0.247812525
## 252    pi tx he ll la sh rn  0.245863953
## 231       pi tx he ll sh rn  0.244235738
## 74                 ts he la  0.242467866
## 77                 ts ll la  0.235523342
## 188          pi ts he sh rn  0.230376725
## 28                    he la  0.229302484
## 237       pi ts he la sh rn  0.229258433
## 196          pi he la sh rn  0.224288799
## 123             pi he sh rn  0.220298257
## 224       pi tx ts he sh rn  0.215712691
## 191          pi ts ll sh rn  0.215479276
## 250    pi tx ts he la sh rn  0.215215254
## 31                    ll la  0.213820203
## 238       pi ts ll la sh rn  0.211898277
## 129             tx ts he la  0.210380798
## 232       pi tx he la sh rn  0.209815347
## 126             pi ll sh rn  0.209297849
## 197          pi ll la sh rn  0.207958665
## 178          pi tx he sh rn  0.205397804
## 227       pi tx ts ll sh rn  0.201789655
## 132             tx ts ll la  0.201211406
## 251    pi tx ts ll la sh rn  0.198968993
## 64                 tx he la  0.198338172
## 181          pi tx ll sh rn  0.194681966
## 233       pi tx ll la sh rn  0.194406969
## 185          pi ts he ll rn  0.192255803
## 120             pi he ll rn  0.185848052
## 67                 tx ll la  0.181303411
## 235       pi ts he ll la rn  0.181091669
## 221       pi tx ts he ll rn  0.176038830
## 194          pi he ll la rn  0.175245960
## 175          pi tx he ll rn  0.168996756
## 117             pi ts sh rn  0.167762547
## 248    pi tx ts he ll la rn  0.166808183
## 192          pi ts la sh rn  0.166475802
## 127             pi la sh rn  0.163963156
## 57                 pi sh rn  0.161709460
## 230       pi tx he ll la rn  0.160649567
## 215          ts he ll sh rn  0.154027201
## 172          pi tx ts sh rn  0.152589898
## 228       pi tx ts la sh rn  0.151701773
## 182          pi tx la sh rn  0.148061243
## 246       ts he ll la sh rn  0.147542803
## 242       tx ts he ll sh rn  0.146259432
## 107             pi tx sh rn  0.145287895
## 254    tx ts he ll la sh rn  0.140281578
## 160             he ll sh rn  0.134233521
## 218          he ll la sh rn  0.130417125
## 210          tx he ll sh rn  0.129193794
## 111             pi ts he rn  0.126862380
## 114             pi ts ll rn  0.126000979
## 51                 pi he rn  0.125042635
## 245       tx he ll la sh rn  0.124984536
## 54                 pi ll rn  0.123883488
## 190          pi ts ll la rn  0.114148000
## 187          pi ts he la rn  0.113881122
## 122             pi he la rn  0.110054167
## 153             ts he sh rn  0.109737907
## 125             pi ll la rn  0.109690649
## 166          pi tx ts he rn  0.109633241
## 169          pi tx ts ll rn  0.108629441
## 3                        ts  0.108567961
## 101             pi tx he rn  0.106152761
## 104             pi tx ll rn  0.104500888
## 216          ts he la sh rn  0.103048478
## 203          tx ts he sh rn  0.100617917
## 156             ts ll sh rn  0.099657303
## 226       pi tx ts ll la rn  0.099180306
## 223       pi tx ts he la rn  0.098864478
## 243       tx ts he la sh rn  0.094401269
## 177          pi tx he la rn  0.094193404
## 217          ts ll la sh rn  0.093891741
## 180          pi tx ll la rn  0.093732802
## 206          tx ts ll sh rn  0.090489495
## 88                 he sh rn  0.090079442
## 161             he la sh rn  0.086875345
## 244       tx ts ll la sh rn  0.085281651
## 143             tx he sh rn  0.084687290
## 91                 ll sh rn  0.084268570
## 211          tx he la sh rn  0.080120941
## 162             ll la sh rn  0.079247621
## 146             tx ll sh rn  0.076457723
## 212          tx ll la sh rn  0.071741597
## 24                    ts la  0.071131740
## 150             ts he ll rn  0.047717537
## 15                    pi rn  0.043131056
## 47                 pi ts rn  0.042347051
## 82                 ts sh rn  0.042345953
## 200          tx ts he ll rn  0.040654503
## 157             ts la sh rn  0.038125883
## 214          ts he ll la rn  0.034840993
## 16                    tx ts  0.033520170
## 137             tx ts sh rn  0.032313866
## 116             pi ts la rn  0.030929355
## 241       tx ts he ll la rn  0.029303915
## 207          tx ts la sh rn  0.028380474
## 85                 he ll rn  0.027695116
## 36                    sh rn  0.027179417
## 56                 pi la rn  0.026349559
## 97              pi tx ts rn  0.025946396
## 140             tx he ll rn  0.025439677
## 92                 la sh rn  0.024189659
## 42                 pi tx rn  0.023694090
## 72                 tx sh rn  0.018612735
## 171          pi tx ts la rn  0.016920311
## 60                 tx ts la  0.016397303
## 147             tx la sh rn  0.015215390
## 159             he ll la rn  0.013575043
## 209          tx he ll la rn  0.011457961
## 106             pi tx la rn  0.010885404
## 6                        la -0.002032533
## 79                 ts ll rn -0.025688966
## 76                 ts he rn -0.028034469
## 131             tx ts he rn -0.033770306
## 134             tx ts ll rn -0.034038408
## 155             ts ll la rn -0.036600842
## 30                    he rn -0.040107777
## 33                    ll rn -0.041383944
## 205          tx ts ll la rn -0.042789210
## 66                 tx he rn -0.043114342
## 152             ts he la rn -0.044342902
## 69                 tx ll rn -0.047645255
## 202          tx ts he la rn -0.048044781
## 19                    tx la -0.051634438
## 90                 ll la rn -0.057729504
## 2                        tx -0.059336927
## 145             tx ll la rn -0.061322228
## 87                 he la rn -0.064254645
## 142             tx he la rn -0.064642219
## 26                    ts rn -0.127513581
## 62                 tx ts rn -0.131480755
## 81                 ts la rn -0.136261768
## 136             tx ts la rn -0.138283011
## 8                        rn -0.144640306
## 21                    tx rn -0.150065747
## 71                 tx la rn -0.161659437
## 35                    la rn -0.162508458
```

maximum correlation is with traits `pi`, `he` and `ll`
So now we need to reduce the traits of the full matrix



```r
good_traits <- colnames(data$traits) %in% c("pi", "he", "ll")
Btcap <- data$traits[,good_traits]

set.seed(12345)

syncsa(comm = data$community,
       traits = Btcap,
       dist.spp = data$dist.spp,
       envir = data$environmental,
       method = "pearson", dist = "euclidean", scale = TRUE,
       scale.envir = TRUE, permutations = 999, na.rm = FALSE, notification = TRUE)
```

```
## $Notes
##             Correlation meanings                                                                                             
## note.roTE   "Trait-convergence assembly patterns (TCAP): roTE"                                                               
## note.roXE   "Both trait-convergence assembly patterns and trait-divergence assembly patterns: roXE"                          
## note.roXE.T "Trait-divergence assembly patterns (TDAP): roXE.T"                                                              
## note.roBF   "Phylogenetic signal at species level: roBF"                                                                     
## note.roPE   "Correlation of phylogenetically structured assembly patterns to ecological variables: roPE"                     
## note.roPT   "Correlation of phylogenetically structured assembly patterns to trait-convergence assembly patterns: roPT"      
## note.roPX.T "Correlation of phylogenetically structured assembly patterns to trait-divergence assembly patterns: roPX.T"     
## note.roTE.P "Removing phylogeny from trait-convergence assembly patterns: roTE.P"                                            
## note.roXE.P "Removing phylogeny from both trait-convergence assembly patterns and trait-divergence assembly patterns: roXE.P"
## 
## $Statistics
##        Obs          p    
## roTE   0.6141292    0.007
## roXE   0.6284366    0.002
## roPE   0.2259687    0.266
## roPT   0.5257162    0.015
## roPX.T -0.001320636 0.586
## roXE.T 0.1820626    0.138
## roTE.P 0.5977546    0.007
## roXE.P 0.610931     0.005
## roBF   0.1861692    0.001
```

# tdap


```r
optimal(comm = data$community,
        traits = data$traits,
        envir = data$environmental,
        subset.min = 1, subset.max = 8, pattern = "tdap", dist = "euclidean", method = "pearson", scale = TRUE, scale.envir = TRUE, na.rm = FALSE, notification = TRUE, progressbar = FALSE)
```

```
##                      Subset            ro
## 1                        pi  6.465601e-01
## 57                 pi sh rn  5.757221e-01
## 123             pi he sh rn  5.358922e-01
## 14                    pi sh  5.001465e-01
## 126             pi ll sh rn  4.781566e-01
## 50                 pi he sh  4.573980e-01
## 195          pi he ll sh rn  4.404758e-01
## 127             pi la sh rn  4.226282e-01
## 196          pi he la sh rn  3.917902e-01
## 55                 pi la sh  3.447616e-01
## 197          pi ll la sh rn  3.375049e-01
## 121             pi he la sh  3.251723e-01
## 239       pi he ll la sh rn  3.169253e-01
## 91                 ll sh rn  3.167515e-01
## 107             pi tx sh rn  2.842948e-01
## 178          pi tx he sh rn  2.548702e-01
## 7                        sh  2.544834e-01
## 35                    la rn  2.418741e-01
## 11                    pi he  2.306914e-01
## 53                 pi ll sh  2.306546e-01
## 13                    pi la  2.229563e-01
## 41                 pi tx sh  2.221962e-01
## 12                    pi ll  2.208353e-01
## 162             ll la sh rn  2.201272e-01
## 181          pi tx ll sh rn  1.998940e-01
## 160             he ll sh rn  1.895780e-01
## 100             pi tx he sh  1.891194e-01
## 48                 pi he ll  1.820626e-01
## 231       pi tx he ll sh rn  1.811218e-01
## 119             pi he ll sh  1.752669e-01
## 36                    sh rn  1.693813e-01
## 71                 tx la rn  1.652083e-01
## 16                    tx ts  1.638733e-01
## 4                        he  1.513837e-01
## 218          he ll la sh rn  1.377480e-01
## 232       pi tx he la sh rn  1.372191e-01
## 182          pi tx la sh rn  1.280092e-01
## 9                     pi tx  1.220085e-01
## 136             tx ts la rn  1.201388e-01
## 21                    tx rn  1.195642e-01
## 233       pi tx ll la sh rn  1.178927e-01
## 124             pi ll la sh  1.151653e-01
## 49                 pi he la  1.126366e-01
## 60                 tx ts la  1.061156e-01
## 252    pi tx he ll la sh rn  1.055962e-01
## 3                        ts  9.956459e-02
## 117             pi ts sh rn  9.302374e-02
## 81                 ts la rn  9.198454e-02
## 2                        tx  8.733694e-02
## 62                 tx ts rn  8.623991e-02
## 188          pi ts he sh rn  8.281331e-02
## 193          pi he ll la sh  8.129066e-02
## 15                    pi rn  7.758430e-02
## 46                 pi ts sh  7.387032e-02
## 86                 he la sh  7.004100e-02
## 32                    ll sh  6.374576e-02
## 56                 pi la rn  6.110843e-02
## 34                    la sh  5.858587e-02
## 92                 la sh rn  5.765987e-02
## 105             pi tx la sh  5.238202e-02
## 90                 ll la rn  5.186724e-02
## 110             pi ts he sh  4.465651e-02
## 19                    tx la  4.368855e-02
## 146             tx ll sh rn  4.100044e-02
## 176          pi tx he la sh  3.882478e-02
## 191          pi ts ll sh rn  3.735606e-02
## 89                 ll la sh  3.618623e-02
## 212          tx ll la sh rn  3.147485e-02
## 26                    ts rn  2.134139e-02
## 5                        ll  2.124329e-02
## 236       pi ts he ll sh rn  2.082750e-02
## 29                    he sh  1.913319e-02
## 88                 he sh rn  1.397849e-02
## 6                        la  1.328522e-02
## 24                    ts la  1.158308e-02
## 103             pi tx ll sh  1.119434e-02
## 51                 pi he rn  5.478274e-03
## 237       pi ts he la sh rn  5.188775e-03
## 8                        rn  1.882440e-09
## 192          pi ts la sh rn -1.339592e-03
## 27                    he ll -2.043974e-03
## 238       pi ts ll la sh rn -6.593050e-03
## 161             he la sh rn -1.972021e-02
## 253    pi ts he ll la sh rn -2.238765e-02
## 40                 pi tx la -2.593300e-02
## 245       tx he ll la sh rn -2.964904e-02
## 210          tx he ll sh rn -2.975435e-02
## 174          pi tx he ll sh -3.212771e-02
## 106             pi tx la rn -3.378838e-02
## 84                 he ll sh -3.456379e-02
## 122             pi he la rn -3.705998e-02
## 158             he ll la sh -4.583761e-02
## 145             tx ll la rn -4.760375e-02
## 217          ts ll la sh rn -5.104893e-02
## 125             pi ll la rn -5.136046e-02
## 54                 pi ll rn -5.206279e-02
## 38                 pi tx he -5.257680e-02
## 42                 pi tx rn -5.378044e-02
## 179          pi tx ll la sh -6.725810e-02
## 224       pi tx ts he sh rn -7.096372e-02
## 156             ts ll sh rn -7.569799e-02
## 227       pi tx ts ll sh rn -7.885124e-02
## 251    pi tx ts ll la sh rn -8.085480e-02
## 202          tx ts he la rn -8.366131e-02
## 250    pi tx ts he la sh rn -8.366485e-02
## 172          pi tx ts sh rn -8.446719e-02
## 205          tx ts ll la rn -8.585178e-02
## 249    pi tx ts he ll sh rn -8.678054e-02
## 228       pi tx ts la sh rn -8.760766e-02
## 171          pi tx ts la rn -9.011574e-02
## 155             ts ll la rn -9.437523e-02
## 255 pi tx ts he ll la sh rn -9.442922e-02
## 68                 tx ll sh -9.482289e-02
## 244       tx ts ll la sh rn -9.601514e-02
## 120             pi he ll rn -9.618721e-02
## 186          pi ts he la sh -1.031478e-01
## 144             tx ll la sh -1.045120e-01
## 246       ts he ll la sh rn -1.052486e-01
## 116             pi ts la rn -1.097621e-01
## 115             pi ts la sh -1.130788e-01
## 177          pi tx he la rn -1.178355e-01
## 147             tx la sh rn -1.190691e-01
## 194          pi he ll la rn -1.193234e-01
## 229       pi tx he ll la sh -1.219724e-01
## 65                 tx he sh -1.274858e-01
## 142             tx he la rn -1.315328e-01
## 99              pi tx he la -1.318866e-01
## 58                 tx ts he -1.335859e-01
## 215          ts he ll sh rn -1.352172e-01
## 101             pi tx he rn -1.354397e-01
## 206          tx ts ll sh rn -1.370531e-01
## 254    tx ts he ll la sh rn -1.397830e-01
## 118             pi he ll la -1.430292e-01
## 157             ts la sh rn -1.431012e-01
## 113             pi ts ll sh -1.432745e-01
## 207          tx ts la sh rn -1.471536e-01
## 180          pi tx ll la rn -1.475708e-01
## 37                 pi tx ts -1.484121e-01
## 97              pi tx ts rn -1.498577e-01
## 72                 tx sh rn -1.507243e-01
## 17                    tx he -1.526112e-01
## 141             tx he la sh -1.539874e-01
## 211          tx he la sh rn -1.558413e-01
## 95              pi tx ts la -1.593098e-01
## 184          pi ts he ll sh -1.611581e-01
## 131             tx ts he rn -1.629436e-01
## 10                    pi ts -1.632946e-01
## 152             ts he la rn -1.655730e-01
## 28                    he la -1.661664e-01
## 139             tx he ll sh -1.676850e-01
## 129             tx ts he la -1.686676e-01
## 223       pi tx ts he la rn -1.695161e-01
## 31                    ll la -1.701657e-01
## 47                 pi ts rn -1.719584e-01
## 69                 tx ll rn -1.726112e-01
## 165          pi tx ts he sh -1.781380e-01
## 70                 tx la sh -1.782619e-01
## 45                 pi ts la -1.783914e-01
## 143             tx he sh rn -1.786516e-01
## 242       tx ts he ll sh rn -1.796210e-01
## 78                 ts ll sh -1.798941e-01
## 96              pi tx ts sh -1.805047e-01
## 243       tx ts he la sh rn -1.833154e-01
## 154             ts ll la sh -1.834936e-01
## 216          ts he la sh rn -1.837040e-01
## 208          tx he ll la sh -1.857747e-01
## 33                    ll rn -1.865642e-01
## 64                 tx he la -1.880119e-01
## 134             tx ts ll rn -1.888329e-01
## 75                 ts he sh -1.914041e-01
## 87                 he la rn -1.923078e-01
## 187          pi ts he la rn -1.933628e-01
## 20                    tx sh -1.985343e-01
## 22                    ts he -2.032856e-01
## 226       pi tx ts ll la rn -2.034113e-01
## 189          pi ts ll la sh -2.050702e-01
## 230       pi tx he ll la rn -2.054817e-01
## 82                 ts sh rn -2.058268e-01
## 241       tx ts he ll la rn -2.058797e-01
## 66                 tx he rn -2.087055e-01
## 204          tx ts ll la sh -2.094169e-01
## 190          pi ts ll la rn -2.095802e-01
## 133             tx ts ll sh -2.099553e-01
## 52                 pi ll la -2.116394e-01
## 104             pi tx ll rn -2.157672e-01
## 80                 ts la sh -2.168120e-01
## 151             ts he la sh -2.187454e-01
## 137             tx ts sh rn -2.196801e-01
## 135             tx ts la sh -2.242550e-01
## 166          pi tx ts he rn -2.261582e-01
## 74                 ts he la -2.290735e-01
## 159             he ll la rn -2.298849e-01
## 234       pi ts he ll la sh -2.311698e-01
## 25                    ts sh -2.314254e-01
## 209          tx he ll la rn -2.319249e-01
## 83                 he ll la -2.325886e-01
## 130             tx ts he sh -2.340933e-01
## 168          pi tx ts ll sh -2.352685e-01
## 61                 tx ts sh -2.382295e-01
## 149             ts he ll sh -2.383295e-01
## 214          ts he ll la rn -2.384871e-01
## 170          pi tx ts la sh -2.408493e-01
## 79                 ts ll rn -2.423306e-01
## 98              pi tx he ll -2.433455e-01
## 67                 tx ll la -2.452968e-01
## 203          tx ts he sh rn -2.495668e-01
## 213          ts he ll la sh -2.500962e-01
## 111             pi ts he rn -2.513377e-01
## 153             ts he sh rn -2.514684e-01
## 248    pi tx ts he ll la rn -2.548061e-01
## 175          pi tx he ll rn -2.561310e-01
## 43                 pi ts he -2.600694e-01
## 235       pi ts he ll la rn -2.617222e-01
## 201          tx ts he la sh -2.638972e-01
## 39                 pi tx ll -2.645847e-01
## 225       pi tx ts ll la sh -2.685241e-01
## 109             pi ts he la -2.696018e-01
## 222       pi tx ts he la sh -2.710639e-01
## 18                    tx ll -2.732937e-01
## 220       pi tx ts he ll sh -2.748287e-01
## 199          tx ts he ll sh -2.764773e-01
## 169          pi tx ts ll rn -2.773705e-01
## 132             tx ts ll la -2.787186e-01
## 240       tx ts he ll la sh -2.847524e-01
## 93              pi tx ts he -2.904740e-01
## 114             pi ts ll rn -3.013958e-01
## 77                 ts ll la -3.045120e-01
## 59                 tx ts ll -3.051487e-01
## 76                 ts he rn -3.053269e-01
## 30                    he rn -3.061558e-01
## 164          pi tx ts he la -3.089439e-01
## 200          tx ts he ll rn -3.153342e-01
## 221       pi tx ts he ll rn -3.168906e-01
## 63                 tx he ll -3.200793e-01
## 247    pi tx ts he ll la sh -3.236755e-01
## 23                    ts ll -3.332869e-01
## 185          pi ts he ll rn -3.335564e-01
## 138             tx he ll la -3.427281e-01
## 102             pi tx ll la -3.539527e-01
## 128             tx ts he ll -3.649077e-01
## 73                 ts he ll -3.696702e-01
## 140             tx he ll rn -3.707636e-01
## 173          pi tx he ll la -3.831132e-01
## 148             ts he ll la -3.858527e-01
## 150             ts he ll rn -3.916496e-01
## 198          tx ts he ll la -3.932905e-01
## 108             pi ts he ll -3.938768e-01
## 85                 he ll rn -3.944190e-01
## 94              pi tx ts ll -4.670619e-01
## 163          pi tx ts he ll -4.720595e-01
## 167          pi tx ts ll la -4.734093e-01
## 183          pi ts he ll la -4.792643e-01
## 44                 pi ts ll -4.868909e-01
## 112             pi ts ll la -4.947646e-01
## 219       pi tx ts he ll la -5.304715e-01
```

# tcap and tdap


```r
optimal(comm = data$community,
        traits = data$traits,
        envir = data$environmental,
        subset.min = 1, subset.max = 8, pattern = "tcap.tdap", dist = "euclidean", method = "pearson", scale = TRUE, scale.envir = TRUE, na.rm = FALSE, notification = TRUE, progressbar = FALSE)
```

```
##                      Subset            ro
## 50                 pi he sh  0.6664868850
## 121             pi he la sh  0.6561333087
## 1                        pi  0.6288751821
## 48                 pi he ll  0.6284365721
## 14                    pi sh  0.6221975733
## 11                    pi he  0.5989111195
## 55                 pi la sh  0.5893706058
## 119             pi he ll sh  0.5868183269
## 12                    pi ll  0.5863973380
## 100             pi tx he sh  0.5687585697
## 193          pi he ll la sh  0.5618834059
## 110             pi ts he sh  0.5497926631
## 53                 pi ll sh  0.5459830941
## 176          pi tx he la sh  0.5338591289
## 124             pi ll la sh  0.5148789221
## 174          pi tx he ll sh  0.5123100340
## 41                 pi tx sh  0.5060730771
## 49                 pi he la  0.4889704258
## 46                 pi ts sh  0.4854513253
## 184          pi ts he ll sh  0.4834354369
## 118             pi he ll la  0.4777932094
## 186          pi ts he la sh  0.4758474394
## 229       pi tx he ll la sh  0.4737855830
## 103             pi tx ll sh  0.4698184582
## 105             pi tx la sh  0.4587553898
## 113             pi ts ll sh  0.4377343527
## 179          pi tx ll la sh  0.4269708941
## 165          pi tx ts he sh  0.4264401194
## 234       pi ts he ll la sh  0.4226575204
## 38                 pi tx he  0.4096525935
## 86                 he la sh  0.4089594926
## 98              pi tx he ll  0.4069526723
## 29                    he sh  0.4061659172
## 52                 pi ll la  0.4054081384
## 115             pi ts la sh  0.4024269308
## 13                    pi la  0.3994465323
## 220       pi tx ts he ll sh  0.3916546434
## 84                 he ll sh  0.3875015427
## 189          pi ts ll la sh  0.3761246111
## 43                 pi ts he  0.3665899961
## 158             he ll la sh  0.3648089995
## 123             pi he sh rn  0.3645424116
## 195          pi he ll sh rn  0.3620040555
## 222       pi tx ts he la sh  0.3590564060
## 4                        he  0.3566177612
## 32                    ll sh  0.3540802481
## 108             pi ts he ll  0.3530413264
## 96              pi tx ts sh  0.3522947675
## 27                    he ll  0.3474864844
## 7                        sh  0.3463620987
## 168          pi tx ts ll sh  0.3445095117
## 75                 ts he sh  0.3425284208
## 239       pi he ll la sh rn  0.3400864167
## 39                 pi tx ll  0.3363208955
## 247    pi tx ts he ll la sh  0.3347305754
## 196          pi he la sh rn  0.3329334205
## 89                 ll la sh  0.3273670317
## 57                 pi sh rn  0.3269744802
## 126             pi ll sh rn  0.3236859714
## 139             tx he ll sh  0.3235822630
## 34                    la sh  0.3221911541
## 65                 tx he sh  0.3218540541
## 9                     pi tx  0.3198950875
## 99              pi tx he la  0.3173558949
## 149             ts he ll sh  0.3105369336
## 173          pi tx he ll la  0.2998211543
## 197          pi ll la sh rn  0.2982385507
## 170          pi tx ts la sh  0.2903744827
## 225       pi tx ts ll la sh  0.2902153217
## 141             tx he la sh  0.2902038444
## 231       pi tx he ll sh rn  0.2874939774
## 208          tx he ll la sh  0.2874131581
## 68                 tx ll sh  0.2872348761
## 5                        ll  0.2866380513
## 44                 pi ts ll  0.2822524305
## 109             pi ts he la  0.2758255230
## 78                 ts ll sh  0.2756108524
## 127             pi la sh rn  0.2752004668
## 10                    pi ts  0.2750415838
## 151             ts he la sh  0.2708371643
## 178          pi tx he sh rn  0.2694612443
## 252    pi tx he ll la sh rn  0.2663097069
## 25                    ts sh  0.2616955284
## 213          ts he ll la sh  0.2615693276
## 236       pi ts he ll sh rn  0.2563880095
## 183          pi ts he ll la  0.2558054797
## 130             tx ts he sh  0.2527792001
## 93              pi tx ts he  0.2469941777
## 144             tx ll la sh  0.2467028087
## 181          pi tx ll sh rn  0.2453789512
## 232       pi tx he la sh rn  0.2446727147
## 188          pi ts he sh rn  0.2440636350
## 199          tx ts he ll sh  0.2430414141
## 102             pi tx ll la  0.2328554582
## 154             ts ll la sh  0.2282772859
## 253    pi ts he ll la sh rn  0.2279114284
## 40                 pi tx la  0.2277840781
## 20                    tx sh  0.2250208355
## 233       pi tx ll la sh rn  0.2240130019
## 191          pi ts ll sh rn  0.2147469098
## 163          pi tx ts he ll  0.2138701297
## 107             pi tx sh rn  0.2122785022
## 237       pi ts he la sh rn  0.2103196817
## 133             tx ts ll sh  0.2058350103
## 112             pi ts ll la  0.2039312417
## 22                    ts he  0.2038828116
## 45                 pi ts la  0.2030517099
## 83                 he ll la  0.2019964004
## 240       tx ts he ll la sh  0.1966076706
## 80                 ts la sh  0.1964741035
## 164          pi tx ts he la  0.1913548972
## 70                 tx la sh  0.1895726579
## 249    pi tx ts he ll sh rn  0.1889145520
## 238       pi ts ll la sh rn  0.1886011890
## 201          tx ts he la sh  0.1885958675
## 117             pi ts sh rn  0.1873176068
## 182          pi tx la sh rn  0.1826186267
## 61                 tx ts sh  0.1706804031
## 37                 pi tx ts  0.1697020049
## 160             he ll sh rn  0.1696542357
## 224       pi tx ts he sh rn  0.1677924209
## 255 pi tx ts he ll la sh rn  0.1666289652
## 219       pi tx ts he ll la  0.1637296455
## 204          tx ts ll la sh  0.1628364648
## 120             pi he ll rn  0.1616682471
## 218          he ll la sh rn  0.1610375477
## 94              pi tx ts ll  0.1546338060
## 28                    he la  0.1534125454
## 192          pi ts la sh rn  0.1527751247
## 58                 tx ts he  0.1500208291
## 31                    ll la  0.1488132946
## 17                    tx he  0.1478338252
## 227       pi tx ts ll sh rn  0.1471036842
## 63                 tx he ll  0.1447123050
## 250    pi tx ts he la sh rn  0.1444970535
## 91                 ll sh rn  0.1394214266
## 95              pi tx ts la  0.1378491705
## 194          pi he ll la rn  0.1362502163
## 73                 ts he ll  0.1328443908
## 162             ll la sh rn  0.1322811801
## 3                        ts  0.1299180774
## 251    pi tx ts ll la sh rn  0.1293756715
## 135             tx ts la sh  0.1239723668
## 167          pi tx ts ll la  0.1230976744
## 51                 pi he rn  0.1207679740
## 210          tx he ll sh rn  0.1165113866
## 54                 pi ll rn  0.1149126180
## 16                    tx ts  0.1138911094
## 172          pi tx ts sh rn  0.1085978500
## 74                 ts he la  0.1078691507
## 245       tx he ll la sh rn  0.1074553653
## 23                    ts ll  0.1001533134
## 215          ts he ll sh rn  0.0966823666
## 148             ts he ll la  0.0953583366
## 125             pi ll la rn  0.0941328959
## 122             pi he la rn  0.0935504063
## 228       pi tx ts la sh rn  0.0912648432
## 88                 he sh rn  0.0911424668
## 246       ts he ll la sh rn  0.0878701012
## 129             tx ts he la  0.0857663365
## 146             tx ll sh rn  0.0842209106
## 175          pi tx he ll rn  0.0823130758
## 161             he la sh rn  0.0813350250
## 77                 ts ll la  0.0793714555
## 138             tx he ll la  0.0782823899
## 212          tx ll la sh rn  0.0781179032
## 18                    tx ll  0.0717601972
## 24                    ts la  0.0715147736
## 156             ts ll sh rn  0.0688453996
## 128             tx ts he ll  0.0688408520
## 15                    pi rn  0.0673220506
## 230       pi tx he ll la rn  0.0658132114
## 217          ts ll la sh rn  0.0637038625
## 60                 tx ts la  0.0630974124
## 101             pi tx he rn  0.0626954144
## 64                 tx he la  0.0574800327
## 242       tx ts he ll sh rn  0.0545949920
## 185          pi ts he ll rn  0.0533915560
## 198          tx ts he ll la  0.0504119490
## 254    tx ts he ll la sh rn  0.0498801101
## 36                    sh rn  0.0471094410
## 59                 tx ts ll  0.0467135461
## 235       pi ts he ll la rn  0.0444710399
## 177          pi tx he la rn  0.0442577477
## 143             tx he sh rn  0.0414048363
## 111             pi ts he rn  0.0409500755
## 56                 pi la rn  0.0408656634
## 104             pi tx ll rn  0.0391641540
## 132             tx ts ll la  0.0386559356
## 67                 tx ll la  0.0346158662
## 92                 la sh rn  0.0333725739
## 211          tx he la sh rn  0.0321938625
## 180          pi tx ll la rn  0.0296560578
## 187          pi ts he la rn  0.0287780467
## 153             ts he sh rn  0.0278406628
## 244       tx ts ll la sh rn  0.0266493724
## 206          tx ts ll sh rn  0.0264211436
## 216          ts he la sh rn  0.0223249274
## 114             pi ts ll rn  0.0153361961
## 190          pi ts ll la rn  0.0148896132
## 42                 pi tx rn  0.0101737841
## 166          pi tx ts he rn  0.0091529694
## 248    pi tx ts he ll la rn  0.0056674525
## 223       pi tx ts he la rn  0.0055721860
## 221       pi tx ts he ll rn  0.0048839990
## 6                        la -0.0001017273
## 106             pi tx la rn -0.0012311977
## 203          tx ts he sh rn -0.0029611093
## 243       tx ts he la sh rn -0.0039553346
## 2                        tx -0.0052540529
## 47                 pi ts rn -0.0068540963
## 116             pi ts la rn -0.0102054701
## 72                 tx sh rn -0.0112708268
## 82                 ts sh rn -0.0160400391
## 147             tx la sh rn -0.0169935747
## 157             ts la sh rn -0.0175983215
## 226       pi tx ts ll la rn -0.0185413754
## 19                    tx la -0.0233963114
## 171          pi tx ts la rn -0.0257652002
## 169          pi tx ts ll rn -0.0275831200
## 97              pi tx ts rn -0.0320549234
## 159             he ll la rn -0.0383703459
## 207          tx ts la sh rn -0.0401735380
## 85                 he ll rn -0.0418507517
## 90                 ll la rn -0.0465658829
## 137             tx ts sh rn -0.0468472732
## 33                    ll rn -0.0625638034
## 209          tx he ll la rn -0.0659226575
## 140             tx he ll rn -0.0690797801
## 155             ts ll la rn -0.0697028177
## 214          ts he ll la rn -0.0702238273
## 145             tx ll la rn -0.0728059422
## 241       tx ts he ll la rn -0.0767419132
## 205          tx ts ll la rn -0.0768185732
## 202          tx ts he la rn -0.0793058952
## 136             tx ts la rn -0.0820093466
## 69                 tx ll rn -0.0833933593
## 131             tx ts he rn -0.0893652541
## 150             ts he ll rn -0.0918787303
## 66                 tx he rn -0.0932187958
## 79                 ts ll rn -0.0943759581
## 200          tx ts he ll rn -0.0977447427
## 62                 tx ts rn -0.0983596681
## 152             ts he la rn -0.0997918342
## 142             tx he la rn -0.0998452479
## 134             tx ts ll rn -0.1008197026
## 81                 ts la rn -0.1020912630
## 30                    he rn -0.1021321678
## 87                 he la rn -0.1064919472
## 76                 ts he rn -0.1112232022
## 71                 tx la rn -0.1134701834
## 35                    la rn -0.1178244409
## 26                    ts rn -0.1186695612
## 21                    tx rn -0.1217994914
## 8                        rn -0.1446403056
```


---
title: "ScriptSYNCSA_Tutorial.R"
author: "User"
date: "Wed Dec 02 01:29:50 2015"
---
