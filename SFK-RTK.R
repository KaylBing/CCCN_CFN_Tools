# SFK-RTK links
# In Cytoscape, search for "receptor tyrosine kinase" and "SRC-family kinase" and create new network from the combined CFN/CCCN
# In each window:
rtks <- getAllNodes()
sfks <- getAllNodes()
PAG1 <- getSelectedNodes()
# Go back to CFN/CCCN
selectNodes(c(rtks, sfks), by ="id", preserve=F)
# createSubnetwork(nodes=c(rtks, sfks)) with deg 0, and deg 2, then re-select 
selectNodes(c(rtks, sfks), by ="id", preserve=F)
# Note interesting relationships with drugs:
rtksfk.df <- getTableColumns("node")
head(rtksfk.df[order(rtksfk.df$all.dasat, decreasing=T), 2:6])
dasat.changed <- head(rtksfk.df[order(rtksfk.df$all.dasat, decreasing=F), 2:6], 22) # all <-2.25
d.c <- rtksfk.df[which(abs(rtksfk.df$all.dasat)>=2.25), "shared name"]
c.c <- rtksfk.df[which(abs(rtksfk.df$all.criz)>=2.25), "shared name"]
e.c <- rtksfk.df[which(abs(rtksfk.df$all.erl)>=2.25), "shared name"]
selectNodes(d.c, by="id", preserve=F)
selectNodes(c.c, by="id", preserve=F)
selectNodes(e.c, by="id", preserve=F)
# Code to find PTMs and genes
selectNodes(getAllNodes()[getAllNodes() %in% rtklist], by="id", preserve=F)
selectNodes(getAllNodes()[getAllNodes() %in% sfklist], by="id", preserve=F)
# In PAG1 deg 2: are any focus pathways present? 
selectNodes(getAllNodes()[getAllNodes() %in% bioplanet[["Transmembrane transport of small molecules"]]], by="id", preserve=F)
# TWO: "SLC2A1"  "SLC12A2"
selectNodes(getAllNodes()[getAllNodes() %in% bioplanet[["Glycolysis and gluconeogenesis"]]], by="id", preserve=F)
# FIVE:  "MDH2"   "LDHB"   "SLC2A1" "ENO2"   "ENO1"  
# In SFK net deg 2:
selectNodes(getAllNodes()[getAllNodes() %in% bioplanet[["Transmembrane transport of small molecules"]]], by="id", preserve=F)
# "FLVCR1"  "TPR"     "NUP133"  "PRKAR1A" "SLC6A6"  "SLC3A2"  "SLC34A2" "SLC2A10" "SLC20A1" "TFRC"    "ATP8B1"  "ABCC3"  
# For new graphs:
selectNodes(c(rtks, sfks), by ="id", preserve=T)
createSubnetwork()  # NO direct edges but there are mutual friends from composite shortest paths From FYN to ABCC3
selectNodes(getAllNodes()[getAllNodes() %in% bioplanet[["Glycolysis and gluconeogenesis"]]], by="id", preserve=F)
# "GAPDH"   "ALDOA"   "ALDH1A3" "TPI1"    "DLD"     "DLAT"    "LDHA"    "PKM"     "ENO1"    "PGK1"    "PFKP"   
selectNodes(c(rtks, sfks), by ="id", preserve=T)
createSubnetwork()  # NO direct edges
# Include PAG1 for deg 0, 1, 2
selectNodes(c(rtks, sfks, PAG1), by ="id", preserve=F)
createSubnetwork()
# In deg 2 and 3: (@2: No connections)
selectNodes(getAllNodes()[getAllNodes() %in% bioplanet[["Transmembrane transport of small molecules"]]], by="id", preserve=F)
selectNodes(c(rtks, sfks, PAG1), by ="id", preserve=TRUE)
createSubnetwork()
toggleGraphicsDetails()
selectNodes(getAllNodes()[getAllNodes() %in% bioplanet[["Glycolysis and gluconeogenesis"]]], by="id", preserve=F)
selectNodes(c(rtks, sfks, PAG1), by ="id", preserve=TRUE)
createSubnetwork()
# Try between 
netgenes <- unique(sapply(c(rtks, sfks, PAG1),  function (x) unlist(strsplit(x, " ",  fixed=TRUE))[1]))
look.ts <- filter.edges.between( bioplanet[["Transmembrane transport of small molecules"]], netgenes, edge.file=gzalltgene.physical.cfn.merged)
look.gl <- filter.edges.between( bioplanet[["Glycolysis and gluconeogenesis"]], netgenes, edge.file=gzalltgene.physical.cfn.merged)
filter.edges.1("PAG1", edge.file=gzalltgene.physical.cfn.merged) # NA becasue co-clustered PTMs have no PPIs
graph.ts <-  graph.cfn.cccn.medians.check (look.ts, ld=FALSE, gz=TRUE, only.cfn=FALSE)
graph.gl <-  graph.cfn.cccn.medians.check (look.gl, ld=FALSE, gz=TRUE, only.cfn=FALSE)
# For figure
selectNodes(netgenes, by="id", preserve=F)
######
# RTK SFK Heat
rtksfk.tbl <- getTableColumns("node")
rtksfk.df <- rtksfk.tbl[,c(2,8,5,6,7)]
names(rtksfk.df)[1] <- "Gene.Name"
rtksf.df <- rtksfk.df[order(rtksfk.df$Gene.Name),]
rtksfk.m <- rtksfk.df[,c(2:5)]
rownames(rtksfk.m) <- rtksfk.df$Gene.Name
see <- graph.clust6nodnosort.l(rtksfk.m)
see <- graph.clust6d.la(rtksfk.m) #, "/Users/markgrimes/Dropbox/_Work/MG_Papers/KarenGuolin/Figs\ for\ paper ")

#_______________________________________________________________________________________________________________________________________
selectNodes(c(getAllNodes()[grep("CLTC", getAllNodes())], getAllNodes()[grep("EGFR", getAllNodes())]), by="id", preserve=F)