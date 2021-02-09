library(ape) #IMPORT LIBRARY ape; for phylogenetics and viz
setwd('/cloud/project') #SET WORKING DIRECTORY
myData <- read.tree('hackbio.nwk') #IMPORT DATA FROM WORKING DIRECTORY
par(mfrow=c(1,2)) # PARAMETER FOR VISUALIZATION
plot(myData) # INITIAL VIEW
ch <- chronos(myData) # calculate chronology
hc <- as.hclust(hc) # convert to hclust class
plot(hc)
labels = hc$labels # original fasta ID of each sequence
ct = cutree(hc, 10) # cut trees into 2, family cluster and non-family
n = length(labels) # number of input query sequences
dend = as.dendrogram(hc) # create dendogram from hclust
plot(dend) 
library(circlize) #import library for circular plot

circos.par(cell.padding = c(0, 0, 0, 0)) #parameters
circos.initialize("a", xlim = c(0, n)) # initialize with only one sector! important 
circos.track(ylim = c(0, 1), bg.border = NA, track.height = 0.3, 
             panel.fun = function(x, y) { 
               for(i in seq_len(n)) { #seqlen is an internal function
                 circos.text(i-0.5, 0, labels[i], adj = c(0, 0.5), 
                             facing = "clockwise", niceFacing = TRUE,
                             col = ct[labels[i]], cex = 0.5)
               }
             })

library(dendextend)

dend = color_branches(dend, k = 6, col = 1:6)
dend_height = attr(dend, "height")
circos.track(ylim = c(0, dend_height), bg.border = NA, 
             track.height = 0.4, panel.fun = function(x, y) {
               circos.dendrogram(dend)
             })
circos.dendrogram(dend, facing = "inside")

dev.off()
circos.clear()
