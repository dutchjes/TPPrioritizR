## Clustering of Nontarget peaks for Prioritization ####
## Jennifer Schollee, Eawag UChem, July 2016

### Required items ##
## 1) PCA.peaks - Feature matrix with rows for features, column with samples and populated by feature intensities
## 2) peaks - Data frame with all the characteristics of the features in PCA.peaks (format: peak id, mz, max int, and rt)

#######################################################################################################
####################Script from Barbara G and Aurea####################################################
###################For profile clustering##############################################################
#######################################################################################################

################# NORMALIZATION (three different strategies)
library(cluster)

pos <- PCA.peaks
dim(pos) 

## types      afterozone sandfilter afterbio beforebio
## sampleIDs      1          2          3         4
## reorder to     3          4          2         1

pos <- pos[,c(4,3,1,2)]
colnames(pos) <- c(1:4)
head(pos)


### correction for dilution factors
###   beforebio   x5
###   afterbio    x4
###   afterozone  x2
###   sandfilter  x2
###   GAK2        x2
###   GAK3        x2
###   movingbed   x2
###   fixedbed    x2

pos[,1] <- pos[,1]*5
pos[,2] <- pos[,2]*4
pos[,3] <- pos[,3]*2
pos[,4] <- pos[,4]*2

## normalization vectors

max_row_pos <- apply(pos, 1, max)

## do the division

normmax_pos  <- pos[,1:ncol(pos)]/max_row_pos

##########  Dissimilarity Matrix Calculation (daisy befehl, dist = distance matrix calculation)

dis_normmax_pos <- daisy(normmax_pos, metric = "euclidean", stand = FALSE)

# dot product matrix for AFP 
ip_matrix <- matrix(0, nrow=2400, ncol=2400)

for (i in 1:2400){
  for (j in 1:2400){
    v1 <- as.numeric(normmax_pos[i,1:4])
    v2 <- as.numeric(normmax_pos[j,1:4])
    vv <- v1 %*% v2
    ip_matrix[i,j] <- as.numeric(vv[1,1])
  }
}

cluster_ward_normmax_pos <- hclust(dis_normmax_pos, method="ward.D2")

### Cutting the three either at specific hight or number of clusters 
kp <- 30                                        # No. of clusters
method_pos <- cluster_ward_normmax_pos
method_name_pos <- "cluster_ward_normmax_pos"
y_value_pos <- normmax_pos

r.cutree_pos <- cutree(method_pos, k=kp)
tab_pos <- table(r.cutree_pos)

m <- as.matrix(normmax_pos)
my.dist <- function(c) dist(c, method="euclidean")
my.clust <- function(d) hclust(d, method="ward.D2")

png(file = "HeatMap.png")
heatmap(x=m, Colv=NA, distfun=my.dist, hclustfun=my.clust, scale='none',
        col=topo.colors(100))
dev.off()


##### Plotting the Clusters positive

pdf(file=paste(method_name_pos,"_",kp, ".pdf"),  paper="a4", pointsize=10, width=8, height=11)
par(mfrow=c(4,2))

plot(method_pos,  main = "Clusteranlyse", hang=0.1, labels=cutree(method_pos, k=kp), las=2)

plot(method_pos,  main = "Clusteranlyse", hang=0.1, labels=cutree(method_pos, k=kp), las=2)
rect.hclust(method_pos, k=kp, border="red")

for (i in 1:kp){
  plot(x=c(1,4), y=c(0,1), type="n", main=paste("normalized cluster",i),
       xlab="time", ylab="Intensity")
  
  for (j in 1:tab_pos[i]){
    list_in_i <- which(r.cutree_pos==i)
    lines(x=1:4, y=y_value_pos[list_in_i[j],1:4])
  }
  
  plot(x=c(1,4), y=c(0,40000000), type="n", main=paste("unnormalized cluster",i),
       xlab="time", ylab="Intensity")
  
  for (j in 1:tab_pos[i]){
    list_in_i <- which(r.cutree_pos==i)
    lines(x=1:4, y=pos[list_in_i[j],1:4])
  }
}

dev.off ()


############# Summary of clustering ######
summary <- data.frame(matrix(NA, nrow = kp, ncol = 1))
for(i in 1:kp){
  
  summary[i,1] <- length(which(r.cutree_pos==i))
  
}

clust.peaks <- list()
clust.id <- c()
for(i in 1:kp){
  
  clust <- as.data.frame(pos[which(r.cutree_pos==i),])
  clust.peaks[[i]] <- as.data.frame(peaks[which(r.cutree_pos==i),])
  clust.id <- c(clust.id, rep(i, nrow(clust.peaks[[i]])))
  
}

all.clust.peaks <- do.call(rbind, clust.peaks)
all.clust.peaks <- cbind(all.clust.peaks, clust.id)

summary(all.clust.peaks$clust.id) ## summary of the number of features per cluster
nfeatpclust <- as.data.frame(table(all.clust.peaks$clust.id)) ## number of features in each cluster
write.csv(nfeatpclust, file = "Featsperclust.csv")
all.clust.peaks$massdef <- all.clust.peaks$mz - round(all.clust.peaks$mz)


####### Colors assigned manually based on visual inspection of trend

## color coding
## grey - highest in 1
## red - highest in 2
## blue - highest in 3
## orange - highest in 4
## green - highest in 3&4
## pink - ND in 3&4
## yellow - present in 2,3&4
## purple - persistant
## black - other

#### sand filter clustering / trend assignment for max intensity adjusted
clust.col <- c("grey", "orange", "red", "grey", "black", "blue", "black", #7
               "green", "purple", "purple", "purple", "black", "green", "black", "blue", #15
               "red", "yellow", "black", "black", "black", "red", "pink", "pink", #23
               "red", "black", "pink", "grey", "black", "orange", "black")

#### GAK2 clustering / trend assignment for max intensity adjusted
clust.col <- c("grey", "red", "grey", "blue", "blue", "orange", "black",
               "red", "black", "purple", "black", "grey", "black", "pink", "black",
               "red", "blue", "yellow", "green", "pink", "black", "black", "yellow", 
               "black", "black", "black", "black", "pink", "red", "black")

#### GAK3 clustering / trend assignment for max intensity adjusted
clust.col <- c("grey", "red", "grey", "blue", "blue", "black", "orange",
               "purple", "yellow", "black", "black", "grey", "black", "red", "yellow",
               "black", "green", "black", "pink", "black", "pink", "blue", "red", 
               "orange", "yellow", "red", "black", "red", "red", "black")

#### movingbed clustering / trend assignment for max intensity adjusted
clust.col <- c("grey", "red", "grey", "blue", "black", "orange", "purple",
               "purple", "purple", "black", "grey", "grey", "black", "pink", "black",
               "black", "yellow", "green", "pink", "black", "pink", "black", "grey", 
               "pink", "green", "red", "black", "red", "blue", "red")

# #### fixedbed clustering / trend assignment for max intensity adjusted
clust.col <- c("grey", "red", "grey", "blue", "black", "orange", "purple",
               "purple", "purple", "black", "grey", "black", "red", "black", "black",
               "black", "pink", "green", "purple", "pink", "pink", "yellow", "black",
               "green", "green", "black", "red", "red", "black", "black")

table(clust.col) ## number of clusters with each trend

## Graphing of Clusters

a <- bwplot(all.clust.peaks[,"clust.id"] ~ all.clust.peaks[,"mz"], #par.settings = list(box.rectangle = list(fill=clust.col)),
       ylab = list(label = "Cluster ID", cex = 1.2), xlab = list(label = "m/z", cex = 1.2),
       key=list(space="top", columns = 2, cex = 1,
                rectangles = list(col = c("grey", "red", "yellow", "blue", "green", "orange",  "pink", "purple", "black")),
                text = list(c("Trend 1: Well removed during CAS treatment", 
                              "Trend 2: Biological TP, well removed during ozonation", 
                              "Trend 3: Biological TP, not removed during ozonation or post-treatment",
                              "Trend 4: Ozone TP, well removed during post-treatment", 
                              "Trend 5: Ozone TP, not removed during post-treatment", 
                              "Trend 6: TP formed during post-treatment",
                              "Trend 7: Persistent in CAS, removed during ozonation",  
                              "Trend 8:Persistent across all treatment steps", 
                              "Trend 9: Other")),
                title = "Assigned Trend", cex.title = 1
       ),
       par.settings = list(box.rectangle = list(fill=clust.col),
                           par.sub.text = list(font = 1, cex = 1.6, just = "left", x = grid::unit(20, "mm"), y = grid::unit(110, "mm"))),
       sub = paste("(", letters[1], ")", sep = ""), 
       scales = list(x = list(cex = 1.2), y = list(cex = 0.8))
       )

b <- bwplot(all.clust.peaks[,"clust.id"] ~ all.clust.peaks[,"rt"], #par.settings = list(box.rectangle = list(fill=clust.col)),
            ylab = list(label = "Cluster ID", cex = 1.2), xlab = list(label = "retention time (min)", cex = 1.2),
            # key=list(space="top", columns = 4, cex = 0.8,
            #          rectangles = list(col = c("white"), lwd = 0),
            #          text = list(c("Removed in biological treatment", 
            #                        "Biological TP, removed in ozone", 
            #                        "Biological TP, not removed in ozone",
            #                        "Ozone TP, removed in post-treatment", 
            #                        "Ozone TP, not removed in post-treatment", 
            #                        "TP formed in post-treatment",
            #                        "Persistant in bio, removed in ozone",  
            #                        "Persistant", 
            #                        "Other")),
            #          title = "Assigned Trend", cex.title = 0.8, col = "white"
            # ),
            par.settings = list(box.rectangle = list(fill=clust.col),
                                par.sub.text = list(font = 1, cex = 1.6, just = "left", x = grid::unit(20, "mm"), y = grid::unit(110, "mm"))),
            sub = paste("(", letters[2], ")", sep = ""), 
            scales = list(x = list(cex = 1.2), y = list(cex = 0.8))
)

c <- bwplot(all.clust.peaks[,"clust.id"] ~ log10(all.clust.peaks[,"into"]), #par.settings = list(box.rectangle = list(fill=clust.col)),
            ylab = list(label = "Cluster ID", cex = 1.2), xlab = list(label = "log10(intensity)", cex = 1.2),
            # key=list(space="top", columns = 4, cex = 0.8,
            #          rectangles = list(col = c("grey", "red", "yellow", "blue", "green", "orange",  "pink", "purple", "black")),
            #          text = list(c("Removed in biological treatment", 
            #                        "Biological TP, removed in ozone", 
            #                        "Biological TP, not removed in ozone",
            #                        "Ozone TP, removed in post-treatment", 
            #                        "Ozone TP, not removed in post-treatment", 
            #                        "TP formed in post-treatment",
            #                        "Persistant in bio, removed in ozone",  
            #                        "Persistant", 
            #                        "Other")),
            #          title = "Assigned Trend", cex.title = 0.8
            # ),
            par.settings = list(box.rectangle = list(fill=clust.col),
                                par.sub.text = list(font = 1, cex = 1.6, just = "left", x = grid::unit(20, "mm"), y = grid::unit(110, "mm"))),
            sub = paste("(", letters[3], ")", sep = ""), 
            scales = list(x = list(cex = 1.2), y = list(cex = 0.8))
            
)

d <- bwplot(all.clust.peaks[,"clust.id"] ~ all.clust.peaks[,"massdef"], #par.settings = list(box.rectangle = list(fill=clust.col)),
            ylab = list(label = "Cluster ID", cex = 1.2), xlab = list(label = "mass defect", cex = 1.2),
            # key=list(space="top", columns = 4, cex = 0.8,
            #          rectangles = list(col = c("grey", "red", "yellow", "blue", "green", "orange",  "pink", "purple", "black")),
            #          text = list(c("Removed in biological treatment", 
            #                        "Biological TP, removed in ozone", 
            #                        "Biological TP, not removed in ozone",
            #                        "Ozone TP, removed in post-treatment", 
            #                        "Ozone TP, not removed in post-treatment", 
            #                        "TP formed in post-treatment",
            #                        "Persistant in bio, removed in ozone",  
            #                        "Persistant", 
            #                        "Other")),
            #          title = "Assigned Trend", cex.title = 0.8
            # ),
            par.settings = list(box.rectangle = list(fill=clust.col),
                                par.sub.text = list(font = 1, cex = 1.6, just = "left", x = grid::unit(20, "mm"), y = grid::unit(110, "mm"))),
            sub = paste("(", letters[3], ")", sep = ""), 
            scales = list(x = list(cex = 1.2), y = list(cex = 0.8))
            
)

png("BWplots_Clusters.png", res = 100, height = 1600, width = 1200)
print(a, position = c(0,-0.19,1,1), split = c(1,1,1,3), more = TRUE)
print(b, position = c(0,-0.03,1,0.97), split = c(1,2,1,3), more = TRUE)
print(d, position = c(0,0,1,1), split = c(1,3,1,3), more = FALSE)
dev.off()


## Convert individual clusters to grouped trends

clust.trend <- clust.col
clust.trend <- replace(clust.trend, clust.col=="grey", 1)
clust.trend <- replace(clust.trend, clust.col=="red", 2)
clust.trend <- replace(clust.trend, clust.col=="yellow", 3)
clust.trend <- replace(clust.trend, clust.col=="blue", 4)
clust.trend <- replace(clust.trend, clust.col=="green", 5)
clust.trend <- replace(clust.trend, clust.col=="orange", 6)
clust.trend <- replace(clust.trend, clust.col=="pink", 7)
clust.trend <- replace(clust.trend, clust.col=="purple", 8)
clust.trend <- replace(clust.trend, clust.col=="black", 9)

trend.col <- c("grey", "red", "yellow", "blue", "green", "orange", "pink", "purple", "black")

all.clust.peaks$trend <- c()
for(i in 1:nrow(all.clust.peaks)){
  id <- all.clust.peaks$clust.id[i]
  all.clust.peaks$trend[i] <- clust.trend[id]
}

table(all.clust.peaks$trend) ## number of features in each trend
trend.sum <- as.data.frame(table(all.clust.peaks$trend))
trend.sum$color <- trend.col
trend.sum$descrip <- c("Well removed during CAS treatment", 
                       "Biological TP, well removed during ozonation", 
                       "Biological TP, not removed during ozonation or post-treatment",
                       "Ozone TP, well removed during post-treatment", 
                       "Ozone TP, not removed during post-treatment", 
                       "TP formed during post-treatment",
                       "Persistent in CAS, removed during ozonation",  
                       "Persistent across all treatment steps", 
                       "Other")

## Graphing of Trends

a <- bwplot(reorder(all.clust.peaks[,"trend"], - as.numeric(all.clust.peaks[,"trend"])) ~ all.clust.peaks[,"mz"],
            #par.settings = list(box.rectangle = list(fill=clust.col)),
            ylab = list(label = "Trend", cex = 1.2),
            xlab = list(label = "m/z", cex = 1.2),
            key=list(space="top", columns = 2, cex = 1,
                     rectangles = list(col = c("grey", "red", "yellow", "blue", "green", "orange",  "pink", "purple", "black")),
                     text = list(c("Trend 1: Well removed during CAS treatment", 
                                   "Trend 2: Biological TP, well removed during ozonation", 
                                   "Trend 3: Biological TP, not removed during ozonation or post-treatment",
                                   "Trend 4: Ozone TP, well removed during post-treatment", 
                                   "Trend 5: Ozone TP, not removed during post-treatment", 
                                   "Trend 6: TP formed during post-treatment",
                                   "Trend 7: Persistent in CAS, removed during ozonation",  
                                   "Trend 8: Persistent across all treatment steps", 
                                   "Trend 9: Other")),
                     title = "Assigned Trend", cex.title = 1
            ),
            par.settings = list(box.rectangle = list(fill=rev(trend.col)),
                                par.sub.text = list(font = 1, cex = 1.6, just = "left", x = grid::unit(20, "mm"), y = grid::unit(110, "mm"))),
            sub = paste("(", letters[1], ")", sep = ""), 
            scales = list(x = list(cex = 1.2), y = list(cex = 1.2) #, 
                                                        # labels = c("Removed in biological treatment", 
                                                        #                         "Biological TP, removed in ozone",
                                                        #                         "Biological TP, not removed in ozone",
                                                        #                         "Ozone TP, removed in post-treatment",
                                                        #                         "Ozone TP, not removed in post-treatment",
                                                        #                         "TP formed in post-treatment",
                                                        #                         "Persistant in bio, removed in ozone",
                                                        #                         "Persistant",
                                                        #                         "Other")
                                                       # )
                          )
)



b <- bwplot(reorder(all.clust.peaks[,"trend"], - as.numeric(all.clust.peaks[,"trend"])) ~ all.clust.peaks[,"rt"],
            #par.settings = list(box.rectangle = list(fill=clust.col)),
            ylab = list(label = "Trend", cex = 1.2), xlab = list(label = "retention time (min)", cex = 1.2),
            # key=list(space="top", columns = 4, cex = 0.8,
            #          rectangles = list(col = c("white"), lwd = 0),
            #          text = list(c("Removed in biological treatment", 
            #                        "Biological TP, removed in ozone", 
            #                        "Biological TP, not removed in ozone",
            #                        "Ozone TP, removed in post-treatment", 
            #                        "Ozone TP, not removed in post-treatment", 
            #                        "TP formed in post-treatment",
            #                        "Persistant in bio, removed in ozone",  
            #                        "Persistant", 
            #                        "Other")),
            #          title = "Assigned Trend", cex.title = 0.8, col = "white"
            # ),
            par.settings = list(box.rectangle = list(fill=rev(trend.col)),
                                par.sub.text = list(font = 1, cex = 1.6, just = "left", x = grid::unit(20, "mm"), y = grid::unit(110, "mm"))),
            sub = paste("(", letters[2], ")", sep = ""), 
            scales = list(x = list(cex = 1.2), y = list(cex = 1.2))
)

c <- bwplot(reorder(all.clust.peaks[,"trend"], - as.numeric(all.clust.peaks[,"trend"])) ~ log10(all.clust.peaks[,"into"]),
            #par.settings = list(box.rectangle = list(fill=clust.col)),
            ylab = list(label = "Trend", cex = 1.2), xlab = list(label = "log10(intensity)", cex = 1.2),
            # key=list(space="top", columns = 4, cex = 0.8,
            #          rectangles = list(col = c("grey", "red", "yellow", "blue", "green", "orange",  "pink", "purple", "black")),
            #          text = list(c("Removed in biological treatment", 
            #                        "Biological TP, removed in ozone", 
            #                        "Biological TP, not removed in ozone",
            #                        "Ozone TP, removed in post-treatment", 
            #                        "Ozone TP, not removed in post-treatment", 
            #                        "TP formed in post-treatment",
            #                        "Persistant in bio, removed in ozone",  
            #                        "Persistant", 
            #                        "Other")),
            #          title = "Assigned Trend", cex.title = 0.8
            # ),
            par.settings = list(box.rectangle = list(fill=rev(trend.col)),
                                par.sub.text = list(font = 1, cex = 1.6, just = "left", x = grid::unit(20, "mm"), y = grid::unit(110, "mm"))),
            sub = paste("(", letters[3], ")", sep = ""), 
            scales = list(x = list(cex = 1.2), y = list(cex = 1.2))
            
)

d <- bwplot(reorder(all.clust.peaks[,"trend"], - as.numeric(all.clust.peaks[,"trend"])) ~ all.clust.peaks[,"massdef"],
            #par.settings = list(box.rectangle = list(fill=clust.col)),
            ylab = list(label = "Trend", cex = 1.2), xlab = list(label = "mass defect", cex = 1.2),
            # key=list(space="top", columns = 4, cex = 0.8,
            #          rectangles = list(col = c("grey", "red", "yellow", "blue", "green", "orange",  "pink", "purple", "black")),
            #          text = list(c("Removed in biological treatment", 
            #                        "Biological TP, removed in ozone", 
            #                        "Biological TP, not removed in ozone",
            #                        "Ozone TP, removed in post-treatment", 
            #                        "Ozone TP, not removed in post-treatment", 
            #                        "TP formed in post-treatment",
            #                        "Persistant in bio, removed in ozone",  
            #                        "Persistant", 
            #                        "Other")),
            #          title = "Assigned Trend", cex.title = 0.8
            # ),
            par.settings = list(box.rectangle = list(fill=rev(trend.col)),
                                par.sub.text = list(font = 1, cex = 1.6, just = "left", x = grid::unit(20, "mm"), y = grid::unit(110, "mm"))),
            sub = paste("(", letters[3], ")", sep = ""), 
            scales = list(x = list(cex = 1.2), y = list(cex = 1.2))
            
)

png("BWplots_Trends.png", res = 300, height = 4800, width = 3600)
print(a, position = c(0,-0.19,1,1), split = c(1,1,1,3), more = TRUE)
print(b, position = c(0,-0.03,1,0.97), split = c(1,2,1,3), more = TRUE)
print(d, position = c(0,0,1,1), split = c(1,3,1,3), more = FALSE)
dev.off()

