

# Calculating dna p-dsitances and tree distances
# Progammer: Chun-I Chiu
install.packages(reshape2)
library(ape)
library(reshape2)

dist.calculation <-
  function(treeFileName,
           seqFileName,
           savingFileName) {
    # Reading the dna sequence & calculating p-distance
    seqs <- read.nexus.data(seqFileName) # Your input
    write.dna(seqs,
              file = "temp.fas",
              format = "fasta",
              append = FALSE)
    seqs <- read.dna("temp.fas", format = "fasta")
    pDist <-
      dist.dna(
        seqs,
        model = "raw",
        as.matrix = T,
        pairwise.deletion = T
      )
    taxonLength <- nrow(pDist)
    counter1 <- 1
    counter2 <- 1
    for (j in 1:taxonLength) {
      for (k in 1:taxonLength) {
        if (k <= counter2) {
          pDist[j, k] <- NA
        } else {
          pDist[j, k] <- pDist[j, k]
          counter1 <- counter1 + 1
        }
      }
      counter2 <- counter2 + 1
    }
    pDist.col <- melt(pDist)[seq(from = N + 1,
                                 to = N ^ 2,
                                 by = 1), ]
    pDist.col <- pDist.col[complete.cases(pDist.col), ]
    pDist.col <- pDist.col[order(pDist.col[, 1], pDist.col[, 2]), ]
    
    # Reading the tree & calculating distances
    tree <- read.nexus(treeFileName)
    MBDist <- cophenetic(tree)
    counter1 <- 1
    counter2 <- 1
    for (j in 1:taxonLength) {
      for (k in 1:taxonLength) {
        if (k <= counter2) {
          MBDist[j, k] <- NA
        } else {
          MBDist[j, k] <- MBDist[j, k]
          counter1 <- counter1 + 1
        }
      }
      counter2 <- counter2 + 1
    }
    MBDist.col <- melt(MBDist)[seq(from = N + 1,
                                   to = N ^ 2,
                                   by = 1), ]
    MBDist.col <- MBDist.col[complete.cases(MBDist.col), ]
    MBDist.col <- MBDist.col[order(MBDist.col[, 1], MBDist.col[, 2]), ]
    Table.length <- nrow(MBDist.col)
    Output.table <- cbind(MBDist.col, pDist.col)
    Output.table <- Output.table[,-(4:5)]
    colnames(Output.table) <-
      c("Species 1", "Species 2", "Tree distance", "DNA p-distance")
    write.csv(Output.table, savingFileName)
  }
#------------------------------------------------------------------------------------------
dist.calculation("C:/Users/Termite Chiu/Dropbox/Saturation_Plot/D2_tree.nexus", # Input tree
                 "C:/Users/Termite Chiu/Dropbox/Saturation_Plot/D2.nexus", # Input sequence
                 "C:/Users/Termite Chiu/Dropbox/Saturation_Plot/distTable.csv") # Output the csv file
