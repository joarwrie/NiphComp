####################
# Loading libraries
####################

# Libraries
library(ggplot2)
library(data.table)
library(cooccur)
#library(vegan)
#library(httr)
library(sf)
#library(igraph)
#library(CHGeoData)
#library(SwissRiverPlot)

##################
# Importing files
##################

entree=read.csv("tab_distri_ASVs.csv", header = T, sep="\t")
Dates=as.data.table(read.csv("Date_eDNA_sampling.txt", header=T, sep="\t"))
Blast=read.table("Blast_res_AmphipodDB.txt", header=F, sep="\t", dec=".")
colnames(Blast)=c("ASV", "ASV_length", "Ref", "Ref_length", "Ali_length", "gaps", "qcov", "pident")
Taxo=read.table("Taxo_AmphipodDB.txt", header=F, sep="\t")
colnames(Taxo)=c("Ref", "Species")
Ind=as.data.table(read.table("Data_citizen_science.txt", header=T, sep="\t"))
TotDuration=read.table("Duration.txt", header=T, sep="\t")
########!!!!!!!!!!!!!!!!!! The coordinates table is not provided for confidentiality reasons. Please contact the authors if you need them!!!!!!!!!!
Coord=read.table("Coord.txt", header=T, sep="\t")

#################################
# Creating metabarcoding dataset
#################################

# Transform into table
entree$Sample=row.names(entree)
entree=as.data.table(entree)
tableau=melt(entree, variable.name="ASV", value.name="Compte", id.vars="Sample")
tableau=tableau[tableau$Compte!=0,]
# Seperating sample info
tableau$Site=do.call(rbind, strsplit(as.character(tableau$Sample), "-"))[,1]
tableau$Replicat=do.call(rbind, strsplit(as.character(tableau$Sample), "_"))[,2]
tableau$Filter=do.call(rbind, strsplit(as.character(tableau$Sample), "_"))[,1]
tableau$Filter=do.call(rbind, strsplit(as.character(tableau$Filter), "-"))[,2]
tableau=tableau[tableau$Filter!="Undetermined",]

#################
# ASVs filtering
#################
# Calculation of the number of reads for each ASV or for each sample or both
tableau[,"TotASV":=sum(Compte), by=ASV]
tableau[,"TotEch":=sum(Compte), by=Site]
tableau[,"TotASVEch":=sum(Compte), by=.(Site,ASV)]
tableau[,"TotASVFilter":=sum(Compte), by=.(ASV, Site, Filter)]

# Calculation of the percentage of each ASV in each sample (filter unit)
tableau$PourcentASV=tableau$TotASVFilter/tableau$TotASV

# Identification of the maximum percentage of reads in a index control sample
neg=tableau[tableau$Site=="NA",]
maxIndex=max(neg$PourcentASV)

# Correction for index jump (maxIndex=0.0004)
tab_tagJump=tableau[tableau$PourcentASV>maxIndex,]

# Keep only ASVs that are present in at least 2 PCR replicates for each sample
tab_tagJump[,"RepCheck":=.N, by=.(Site, ASV)]
tab_final=tab_tagJump[tab_tagJump$RepCheck>1,]

# Production of Table SX
tab_tmp=as.data.table(table(unique(tab_final[,c(2,4)])$Site))
colnames(tab_tmp)=c("Site", "ASVs")
tab_final[,"TotReads":=sum(Compte), by=Site]
colnames(tab_final)[4]="SiteCode"
TableS2=merge(unique(tab_final[order(TotReads),c(4,13)]), tab_tmp, by.x="SiteCode", by.y="Site", all=T)
write.table(file="TableS2.txt", TableS2, col.names=T, row.names=F, sep="\t")

# Adjusting table to match barcoding sites
ASV=merge(tab_final[,1:6], Dates, by="SiteCode", all=T)

# Creation of a metabarcoding table without correction for index-jump and PCR errors
ASV_woCorr=merge(tableau[,1:6], Dates, by.x="Site", by.y="SiteCode", all=T)

######################
# Amphipod assignment
######################

# Add species names
Blast$Ref=do.call(rbind, strsplit(Blast$Ref, "\\|"))[,1]
Assign=merge(Blast, Taxo, by="Ref", all.x=T, all.y=F)
Assign=as.data.table(Assign)

# Selection of alignments with at least 99% query cover
Assign=Assign[Assign$qcov>=99,]
# Keeping only assignments with more than 90% identity
Assign=Assign[Assign$pident>=95,]
	# Then get the max identity for each ASV
Assign[,"Max":=max(pident), by=ASV]
	# And keep assignments within max and max-1
Assign=Assign[Assign$pident>=(Assign$Max - 1),]

# Listing assignments for all ASVs
Assign[,"Assignment":=paste(unique(.SD[,Species]),collapse="£"), by=ASV]
Assign=unique(Assign, by="ASV")

# Merging datasets
Metab=merge(ASV, Assign, by="ASV", all=F)
Metab_woCorr=merge(ASV_woCorr, Assign, by="ASV", all=F)

###############################
# Creation of barcoding tables
###############################

# Modifying table to get one line per sampling date, site and species
Ind[Ind$Duration==0,]$Duration=7
Ind[,"CountIndOneDate":=sum(nr_ind), by=.(Site, species, Day, Month, Year)]
Ind$CorrectedAbundanceOneDate=Ind$CountIndOneDate/Ind$Duration
Ind=unique(Ind, by=c("Site", "species", "Day", "Month", "Year"))

# Correcting for sampling effort (nb of days sampled)
Ind=merge(Ind, TotDuration, by="Site", all=T)
Ind[,"CountIndTotal":=sum(CountIndOneDate), by=.(Site, species)]
Ind$CorrectedAbundance=Ind$CountIndTotal/Ind$CountDayTotal

#########################################
# Figures Comparaison between approaches
#########################################

# Creation of a table for comparison
	# Table to compare with barcoding of all dates
AllDates=Ind[!(is.na(Ind$species)),]
AllDates=unique(AllDates, by=c("species", "Site"))

Metab[,"CountMetab":=sum(Compte), by=.(Assignment, SiteCode)]
Metab=unique(Metab[,c(19,7,20)])
colnames(Metab)=c("species", "Site", "CountMetab")

AllDates=merge(AllDates, Metab, by=c("species", "Site"), all=T)
AllDates[is.na(AllDates$CountMetab),]$CountMetab=0
AllDates$Type="AllDates"

	# Table to compare with barcoding at the same date
SameDate=merge(Ind, Dates, by=c("Site", "Day", "Month"), all=F)
SameDate=merge(SameDate, Metab, by=c("species", "Site"), all=T)
SameDate[is.na(SameDate$CountMetab),]$CountMetab=0
SameDate$Type="SameDate"

	# Combining both tables
colnames(AllDates)[14]="CountSite"
colnames(SameDate)[11]="CountSite"
tab_comp=rbind(AllDates[,c(1,2,14,15,16)], SameDate[,c(1,2,11,16,17)])
tab_comp[is.na(tab_comp$CountSite),]$CountSite=0
tab_comp=tab_comp[tab_comp$species!="Niphargus_sp." & !(is.na(tab_comp$species)),]

# Adding a column for point color
tab_comp$Color="White"
tab_comp[tab_comp$CountSite>0 & tab_comp$CountMetab>0,]$Color="Green"
tab_comp[tab_comp$CountSite>0 & tab_comp$CountMetab==0,]$Color="Red"
tab_comp[tab_comp$CountSite==0 & tab_comp$CountMetab>0,]$Color="Orange"

# Presence Absence plot for comparison of same dates
tab_graph=dcast(tab_comp[tab_comp$Type=="SameDate",], species~Site, fill="White", value.var="Color", fun.aggregate=paste)
	# Adding missing sites
tab_graph$T165="White"
tab_graph$T035="White"
tab_graph$T124="White"
tab_graph$T094="White"
tab_graph$T066="White"
tab_graph$T170="White"

tab_graph=melt(tab_graph, id.var="species", value.name="Color", variable.name="Site")
tab_graph$Presence="Detection"
tab_graph$Site=factor(tab_graph$Site, levels=rev(sort(unique(as.character(tab_graph$Site)))))

p1=ggplot(tab_graph, aes(x=Presence, y=Site, fill=Color)) +
  geom_point(color="black", shape=21, size=10) +
  facet_grid(Site~species, scales = "free_y", switch="y") +
  theme(panel.background=element_rect(fill="white", colour="black"), axis.line=element_line(colour="black")) +
  theme(strip.background.y=element_rect(fill="white", colour="black"), strip.background.x=element_rect(fill="white", colour="white"), strip.text.x=element_text(angle=90, size=8)) +
  theme(strip.text.y=element_text(size=10, face="bold")) +
  theme(axis.text=element_blank(), axis.title=element_blank(), axis.ticks=element_blank()) +
  theme(legend.position="none", panel.spacing=unit(0, "lines")) +
  scale_fill_manual(values=c("limegreen", "darkorange", "brown2", "white"))

pdf("Figure2A.pdf", width=3, height=9.5)
	p1
dev.off()

# Presence Absence plot for comparison of all dates
tab_graph=dcast(tab_comp[tab_comp$Type=="AllDates",], species~Site, fill="White", value.var="Color", fun.aggregate=paste)
	# Adding missing sites
tab_graph$T165="White"
tab_graph$T035="White"
tab_graph$T124="White"
tab_graph$T094="White"

tab_graph=reshape2::melt(tab_graph, id.var="species", value.name="Color", variable.name="Site")
tab_graph$Presence="Detection"
tab_graph$Site=factor(tab_graph$Site, levels=rev(sort(unique(as.character(tab_graph$Site)))))
Liste_sites=unique(tab_graph$Site)

p2=ggplot(tab_graph, aes(x=Presence, y=Site, fill=Color)) +
  geom_point(color="black", shape=21, size=10) +
  facet_grid(Site~species, scales = "free_y", switch="y") +
  theme(panel.background=element_rect(fill="white", colour="black"), axis.line=element_line(colour="black")) +
  theme(strip.background.y=element_rect(fill="white", colour="black"), strip.background.x=element_rect(fill="white", colour="white"), strip.text.x=element_text(angle=90, size=8)) +
  theme(strip.text.y=element_text(size=10, face="bold")) +
  theme(axis.text=element_blank(), axis.title=element_blank(), axis.ticks=element_blank()) +
  theme(legend.position="none", panel.spacing=unit(0, "lines")) +
  scale_fill_manual(values=c("limegreen", "brown2", "white"))

pdf("Figure2B.pdf", width=3, height=9.5)
	p2
dev.off()

# Barplot for abundance comparison
	# Comparison at same date
tab_graph=unique(tab_comp[tab_comp$Type=="SameDate",])
tab_graph$ModifCount=-tab_graph$CountSite
tab_graph[is.na(tab_graph$ModifCount),]$ModifCount=0
tab_graph=reshape2::melt(tab_graph[,c(1,2,4,7)], id.vars=c("species", "Site"), variable.name="Dataset", value.name="Value")
tab_graph$Dataset=factor(tab_graph$Dataset, levels=c("ModifCount", "CountMetab"))
tab_graph$Site=factor(tab_graph$Site, levels=rev(sort(as.character(Liste_sites))))
Approach=c("Citizen Science", "Metabarcoding")
names(Approach)=c("ModifCount", "CountMetab")
ExpandXlim=data.frame(Site=c("T005", "T005"), Value=c(0, -1.5), species=c(NA, NA), Dataset=c("CountMetab", "ModifCount"))
ExpandXlim$Dataset=factor(ExpandXlim$Dataset, levels=c("ModifCount", "CountMetab"))
ExpandXlim$Site=factor(ExpandXlim$Site, levels=rev(sort(as.character(Liste_sites))))

p3=ggplot(tab_graph, aes(x=Site, y=Value, fill=species)) +
  geom_point(data=ExpandXlim, size=0) +
  geom_bar(stat="identity", position="stack", color="black") +
  coord_flip() +
  facet_wrap(~Dataset, scales="free_x", labeller=labeller(Dataset=Approach)) +
  theme(panel.background=element_rect(fill="white", colour="black"), axis.line=element_line(colour="black"), panel.grid.major.x=element_line(color="gray70")) +
  theme(strip.background=element_rect(fill="white", colour="black"), strip.text.x=element_text(size=14, face="bold")) +
  theme(legend.position="top", legend.text=element_text(size=14), legend.title=element_blank()) +
  theme(axis.text=element_text(size=14), axis.title.y=element_blank(), axis.title.x=element_text(size=16, face="bold")) +
  ylab("\nAbundance") +
  scale_x_discrete(drop=F, labels=rev(sort(as.character(Liste_sites)))) +
  scale_fill_manual(values=c("darkred", "darkcyan", "peachpuff4", "royalblue4", "khaki", "cornflowerblue"))

pdf("FigureS1A.pdf", width=8, height=9.5)
	p3
dev.off()

	# Comparison at all dates
tab_graph=tab_comp[tab_comp$Type=="AllDates",]
tab_graph$ModifCount=-tab_graph$CountSite
tab_graph[is.na(tab_graph$ModifCount),]$ModifCount=0
tab_graph=reshape2::melt(tab_graph[,c(1,2,4,7)], id.vars=c("species", "Site"), variable.name="Dataset", value.name="Value")
tab_graph$Dataset=factor(tab_graph$Dataset, levels=c("ModifCount", "CountMetab"))
tab_graph$Site=factor(tab_graph$Site, levels=rev(sort(as.character(Liste_sites))))
Approach=c("Citizen Science", "Metabarcoding")
names(Approach)=c("ModifCount", "CountMetab")

p4=ggplot(tab_graph, aes(x=Site, y=Value, fill=species)) +
  geom_bar(stat="identity", position="stack", color="black") +
  coord_flip() +
  facet_wrap(~Dataset, scales="free_x", labeller=labeller(Dataset=Approach)) +
  theme(panel.background=element_rect(fill="white", colour="black"), axis.line=element_line(colour="black"), panel.grid.major.x=element_line(color="gray70")) +
  theme(strip.background=element_rect(fill="white", colour="black"), strip.text.x=element_text(size=14, face="bold")) +
  theme(legend.position="top", legend.text=element_text(size=14), legend.title=element_blank()) +
  theme(axis.text=element_text(size=14), axis.title.y=element_blank(), axis.title.x=element_text(size=16, face="bold")) +
  ylab("\nAbundance") +
  scale_x_discrete(drop=F, labels=rev(sort(as.character(Liste_sites)))) +
  scale_fill_manual(values=c("darkred", "darkcyan", "peachpuff4", "royalblue4", "khaki", "cornflowerblue"))

pdf("FigureS1B.pdf", width=8, height=9.5)
	p4
dev.off()

# Presence Absence plot for comparison of same date without correction for index-jump and PCR errors
AllDates=Ind[!(is.na(Ind$species)),]
AllDates=unique(AllDates, by=c("species", "Site"))

Metab_woCorr[,"CountMetab":=sum(Compte), by=.(Assignment, Site)]
Metab_woCorr=unique(Metab_woCorr[,c(19,7,20)])
colnames(Metab_woCorr)=c("species", "Site", "CountMetab")

AllDates=merge(AllDates, Metab_woCorr, by=c("species", "Site"), all=T)
AllDates[is.na(AllDates$CountMetab),]$CountMetab=0
AllDates$Type="AllDates"
colnames(AllDates)[14]="CountSite"

SameDate=merge(Ind, Dates, by=c("Site", "Day", "Month"), all=F)
SameDate=merge(SameDate, Metab_woCorr, by=c("species", "Site"), all=T)
SameDate[is.na(SameDate$CountMetab),]$CountMetab=0
SameDate$Type="SameDate"
colnames(SameDate)[11]="CountSite"

tab_comp=rbind(AllDates[,c(1,2,14,15,16)], SameDate[,c(1,2,11,16,17)])
tab_comp[is.na(tab_comp$CountSite),]$CountSite=0
tab_comp=tab_comp[tab_comp$species!="Niphargus_sp." & !(is.na(tab_comp$species)),]

# Adding a column for point color
tab_comp$Color="White"
tab_comp[tab_comp$CountSite>0 & tab_comp$CountMetab>0,]$Color="Green"
tab_comp[tab_comp$CountSite>0 & tab_comp$CountMetab==0,]$Color="Red"
tab_comp[tab_comp$CountSite==0 & tab_comp$CountMetab>0,]$Color="Orange"

# Graph for same date
tab_graph=dcast(tab_comp[tab_comp$Type=="SameDate",], species~Site, fill="White", value.var="Color", fun.aggregate=paste)
	# Adding missing sites
tab_graph$T035="White"
tab_graph$T124="White"
tab_graph$T094="White"
tab_graph$T066="White"
tab_graph$T170="White"

tab_graph=melt(tab_graph, id.var="species", value.name="Color", variable.name="Site")
tab_graph$Presence="Detection"
tab_graph$Site=factor(tab_graph$Site, levels=rev(sort(unique(as.character(tab_graph$Site)))))

p5=ggplot(tab_graph, aes(x=Presence, y=Site, fill=Color)) +
  geom_point(color="black", shape=21, size=10) +
  facet_grid(Site~species, scales = "free_y", switch="y") +
  theme(panel.background=element_rect(fill="white", colour="black"), axis.line=element_line(colour="black")) +
  theme(strip.background.y=element_rect(fill="white", colour="black"), strip.background.x=element_rect(fill="white", colour="white"), strip.text.x=element_text(angle=90, size=8)) +
  theme(strip.text.y=element_text(size=10, face="bold")) +
  theme(axis.text=element_blank(), axis.title=element_blank(), axis.ticks=element_blank()) +
  theme(legend.position="none", panel.spacing=unit(0, "lines")) +
  scale_fill_manual(values=c("limegreen", "darkorange", "brown2", "white"))

pdf("Figure3A.pdf", width=3, height=9.5)
	p5
dev.off()

# Graph for all dates
tab_graph=dcast(tab_comp[tab_comp$Type=="AllDates",], species~Site, fill="White", value.var="Color", fun.aggregate=paste)

	# Adding missing sites
tab_graph$T035="White"
tab_graph$T124="White"
tab_graph$T094="White"

tab_graph=reshape2::melt(tab_graph, id.var="species", value.name="Color", variable.name="Site")
tab_graph$Presence="Detection"
tab_graph$Site=factor(tab_graph$Site, levels=sort(unique(as.character(tab_graph$Site))))
Liste_sites=unique(tab_graph$Site)

p6=ggplot(tab_graph, aes(x=Presence, y=Site, fill=Color)) +
  geom_point(color="black", shape=21, size=10) +
  facet_grid(Site~species, scales = "free_y", switch="y") +
  theme(panel.background=element_rect(fill="white", colour="black"), axis.line=element_line(colour="black")) +
  theme(strip.background.y=element_rect(fill="white", colour="black"), strip.background.x=element_rect(fill="white", colour="white"), strip.text.x=element_text(angle=90, size=8)) +
  theme(strip.text.y=element_text(size=10, face="bold")) +
  theme(axis.text=element_blank(), axis.title=element_blank(), axis.ticks=element_blank()) +
  theme(legend.position="none", panel.spacing=unit(0, "lines")) +
  scale_fill_manual(values=c("limegreen", "darkorange", "brown2", "white"))

pdf("Figure3B.pdf", width=3, height=9.5)
	p6
dev.off()

###################################################
# Distribution of Niphargus species with barcoding
###################################################
# Distribution of species richness
tab_barplot=unique(Ind[Ind$Site%in%ASV$Site & Ind$species!="Niphargus_sp.",], by=c("species", "Site"))
tab_barplot$Site=factor(tab_barplot$Site, levels=Liste_sites)
tab_barplot=as.data.table(table(tab_barplot$Site))

p7=ggplot(data=tab_barplot, aes(x=N)) +
  geom_bar(fill="gray60", color="black") +
  coord_flip() +
  theme(panel.background=element_rect(fill="white", colour="black"), axis.line=element_line(colour="black"), panel.grid.major = element_line("gray90")) +
  theme(legend.position = "none") +
  theme(axis.text.x=element_text(size=14), axis.text.y=element_text(size=14), axis.title=element_text(size=18, face="bold")) +
  xlab("\nSpecies richness per site") + ylab("# sites\n") 

pdf("Figure4A.pdf")
  p7
dev.off()

# Map distribution per species (convex hulls)
########!!!!!!!!!!!!!!!!!! The coordinates table is not provided for confidentiality reasons. Please contact the authors you need them!!!!!!!!!!
tab_hull=merge(Ind, Coord, by="Site", all.x=T, all.y=F)
tab_hull=unique(tab_hull[tab_hull$CountIndTotal>0 & tab_hull$species!="Niphargus_sp.",], by=c("Site", "species"))
tab_hull=sf::st_as_sf(tab_hull, coords = c("Xcoord", "Ycoord"), crs=2056)
hull=st_sf(aggregate(tab_hull$geometry, list(tab_hull$species), function(x){st_convex_hull(st_union(x))}))

p8=ggplot() +
  geom_sf(data=hull, aes(fill=Group.1, colour=Group.1), alpha=0.5) +
  theme(panel.background=element_rect(fill="white", colour="black"), axis.line=element_line(colour="black")) +
  theme(axis.text.x=element_text(size=14), axis.text.y=element_text(size=14), axis.title=element_blank()) +
  theme(legend.title=element_blank()) +
  scale_fill_manual(values=c("darkred", "darkcyan", "peachpuff4", "royalblue4", "khaki", "cornflowerblue")) +
  scale_colour_manual(values=c("darkred", "darkcyan", "peachpuff4", "royalblue4", "khaki", "cornflowerblue"))

pdf("Figure4B.pdf")
  p8
dev.off()

########################
# Cooccurrence analysis
########################
# Creating the matrix for input
	# Barcoding
mat_temp=dcast(AllDates, species~Site, value.var="CountSite", fun.aggregate=sum, fill=0)
colnames(mat_temp)[1]="ASV"
mat_temp=mat_temp[mat_temp$ASV!="Niphargus_sp.",]
mat_temp$T165=0
mat_temp$T035=0
mat_temp$T124=0
mat_temp$T094=0
	# Selecting ASVs present at more than one site
ListeASV=unique(ASV[,c(3,7)])
ListeASV[,"NbSites":=.N, by=ASV]
ListeASV=unique(ListeASV[ListeASV$NbSites>1,]$ASV)
	# Combined matrix with both datasets
mat_cooccur=rbind(dcast(ASV[ASV%in%ListeASV,c(3,4,7)], ASV~Site, value.var="Compte", fun.aggregate=sum, fill=0), mat_temp[!(is.na(mat_temp$ASV)),])
ListeASV=mat_cooccur$ASV
mat_cooccur=as.matrix(mat_cooccur[,-1])
row.names(mat_cooccur)=ListeASV
rm(mat_temp)

# Transform into a presence-absence matrix
mat_cooccur[mat_cooccur>0]=1

# Calculation
Res_cooccur=cooccur(mat_cooccur, thresh=T, spp_names=T, true_rand_classifier=0.1, prob="hyper")

# Getting only significant results
Sign_cooccur=Res_cooccur$results[Res_cooccur$results$p_lt<=0.05 | Res_cooccur$results$p_gt<=0.05,]

  # Selecting only cooccurrences with an amphipod species
List_amphipods=c("Crangonyx_cf._subterraneus", "Niphargus_auerbachi", "Niphargus_fontanus", "Niphargus_puteanus", "Niphargus_thienemanni", "Niphargus_tonywhitteni")
Cooccur_amphipods=Sign_cooccur[Sign_cooccur$sp1_name%in%List_amphipods | Sign_cooccur$sp2_name%in%List_amphipods,]

tab_graph=Cooccur_amphipods[,8:11]
tab_graph$PN="Pos"
tab_graph[tab_graph$p_lt<=0.05,]$PN="Neg"
tab_graph$sp1_name=gsub("ASV_", "", tab_graph$sp1_name)
setorder(tab_graph, cols="sp2_name", "PN")
tab_graph$sp1_name=factor(tab_graph$sp1_name, levels=unique(tab_graph$sp1_name))

# Creation of a cooccurrence heatmap
p9=ggplot(data=tab_graph, aes(x=sp1_name, y=sp2_name, fill=PN))+
	geom_tile(color="black") +
	coord_polar() +
	theme(panel.background=element_blank(), axis.line=element_line(colour="black"), panel.grid.major=element_line(colour="gray60")) +
	theme(legend.position = "top", legend.text=element_text(size=20)) +
	theme(axis.text.x=element_text(size=15, angle=360/(2*pi)*rev(pi/2 + seq( pi/201, 2*pi-pi/201, len=201))), axis.text.y=element_blank()) +
	xlab("") + ylab("") +
	scale_fill_manual(values=c("darkorange", "purple"))
	
pdf("Figure5.pdf", width=15, height=15)
	p9
dev.off()

###################################################
# Correlation between amphipods and eDNA diversity
###################################################

# Preparing table
tab_diversity=unique(ASV, by=c("ASV", "Site"))
tab_diversity[, "NbASVs":=.N, by="Site"]
tab_diversity=unique(tab_diversity[,c(7,10)])
tab_graph=unique(Ind[Ind$species!="Niphargus_sp.",], by=c("Site", "species"))
tab_graph[,"Somme":=sum(CountIndTotal), by=Site]
tab_graph$CorrectedAbundanceTotal=tab_graph$Somme/tab_graph$CountDayTotal
tab_graph=merge(tab_diversity, unique(tab_graph[,c(1,16)]), by="Site", all=T)
tab_graph[is.na(tab_graph$CorrectedAbundanceTotal),]$CorrectedAbundanceTotal=0

Correlation_results=cor.test(log(tab_graph$NbASVs+1), log(tab_graph$CorrectedAbundanceTotal+1), method="pearson", alternative = "two.sided")

p10=ggplot(data=tab_graph, aes(x=log(CorrectedAbundanceTotal+1), y=log(NbASVs+1))) +
  geom_point(size=3) +
  geom_smooth(method="lm", formula="y~x") +
  theme(panel.background=element_rect(fill="white", colour="black"), axis.line=element_line(colour="black"), panel.grid.major = element_line("gray90")) +
  theme(axis.text.x=element_text(size=14), axis.text.y=element_text(size=14), axis.title=element_text(size=18, face="bold")) +
  xlab("\nlog(mean #amphipods per week + 1)") + ylab("log(#ASVs + 1)\n") +
  geom_text(aes(x=0.6, y=4.5), label=paste0("r = ", round(Correlation_results$estimate, 3)), size=6) +
  geom_text(aes(x=0.6, y=4.2), label=paste0("p = ", round(Correlation_results$p.value, 3)), size=6)

pdf("Figure6A.pdf", width=9, height=9)
	p10
dev.off()

tab_graph=unique(Ind[Ind$species!="Niphargus_sp.",], by=c("species", "Site"))
tab_graph[, "NbSpecies":=.N, by=Site]
tab_graph=unique(tab_graph[,c(1,15)])
tab_graph=merge(tab_diversity, tab_graph, by="Site", all.x=T, all.y=F)
tab_graph[is.na(tab_graph$NbSpecies),]$NbSpecies=0

Correlation_results=cor.test(log(tab_graph$NbASVs+1), tab_graph$NbSpecies, method="pearson", alternative = "two.sided")

p11=ggplot(data=tab_graph, aes(x=NbSpecies, y=log(NbASVs+1))) +
  geom_point(size=3) +
  geom_smooth(method="lm", formula="y~x") +
  theme(panel.background=element_rect(fill="white", colour="black"), axis.line=element_line(colour="black"), panel.grid.major = element_line("gray90")) +
  theme(axis.text.x=element_text(size=14), axis.text.y=element_text(size=14), axis.title=element_text(size=18, face="bold")) +
  xlab("\n# Species") + ylab("log(#ASVs)\n") +
  geom_text(aes(x=2.5, y=4.6), label=paste0("r = ", round(Correlation_results$estimate, 3)), size=6) +
  geom_text(aes(x=2.5, y=4.3), label=paste0("p = ", round(Correlation_results$p.value, 3)), size=6)

pdf("Figure6B.pdf", width=9, height=9)
  p11
dev.off()


****************************************









# END







#################### RETRIEVING AREA INFO ##########################
# Calculation of areas of distribution for each species
GetArea<-function(Sp){
  tableau=unique(Ind[Ind$Species==Sp,10:11])
  if(dim(tableau)[1]>2){
    tableau=sf::st_as_sf(tableau, coords = c("Xcoord", "Ycoord"), crs=2056)
    tableau=sf::st_geometry(tableau)
    area=st_area(st_convex_hull(st_union(tableau)))
    return(c(Sp, area))
  }
}

tab_area<-as.data.frame(do.call(rbind, (lapply(unique(Ind$Species), FUN=GetArea))))

##################### GETTING INFO ON AQUIFER TYPE ###################
# Retrieving shapefile
shp_file <- list.files("~/Hydrogeologische_LV95", pattern="shp", full.names=T) # retrieve filename for the shapefile 
aquifer <- sf::st_read(shp_file) # read shapefile. ATTENTION: To work it needs all the asociated files, not only the shapefile itself but also shx, dbf, etc
aquifer <- sf::st_transform(aquifer, crs=2056)
aquifer$aquifertyp <- as.character(aquifer$aquifertyp)

# Get aquifer for each point in a dataframe
  # Retrieving points coordinates
tab_map<-unique(Samples, by="SiteCode")
tab_map$Xcoord<-tab_map$Xcoord + 2000000
tab_map$Ycoord<-tab_map$Ycoord + 1000000
  # Comparing points position with the aquifer map
pts <- sf::st_as_sf(tab_map[,c(17,6,7)], coords = c("Xcoord", "Ycoord"), crs=2056) # transform sites dataframe to sf object; crs gives the coordonate system to be used
match_aquifer<-st_intersects(pts, aquifer, sparse=F)
colnames(match_aquifer)<-unique(aquifer$aquifertyp)
row.names(match_aquifer)<-tab_map$SiteCode
  # Melting data frame to get it long format
match_aquifer<-as.data.frame(match_aquifer)
match_aquifer$SiteCode<-row.names(match_aquifer)
match_aquifer<-reshape2::melt(match_aquifer, id.vars="SiteCode", value.name="match", variable.name="AquiferType")
match_aquifer<-match_aquifer[match_aquifer$match=="TRUE",]

################### GETTING INFO LAND USE TYPE #############################
# Retrieving pixel file : in this table, each lines gives the coordinates of the center of a 100 m2 square
landuse <- read.csv("ag-b-00.03-37-area-csv.csv", header=T, sep=";")

# Getting square coordinates for each point in the csv file
Toess<-landuse[landuse$E>2680000 & landuse$E<2720000 & landuse$N>1240000 & landuse$N<1270000,c(2,3,24)]

Squared<-function(vec){
  mat<-matrix(data=c(vec[1]-50, vec[1]-50, vec[1]+50, vec[1]+50, vec[1]-50, vec[2]-50, vec[2]+50, vec[2]+50, vec[2]-50, vec[2]-50, rep(vec[3], 5)), nrow=5)
  return(mat)
}

Toess<-do.call(rbind, apply(Toess, FUN=Squared, MARGIN=1, simplify=F))
tmpvec<-c()
for (i in 1:(nrow(Toess)/5)){
  tmpvec<-c(tmpvec, rep(i, 5))
}
Toess<-cbind(Toess, tmpvec)

# Transforming into a polygon file
Square_points<-st_as_sf(as.data.frame(Toess), coords=c("V1", "V2"))
Polygons<-st_sf(aggregate(Square_points$geometry, list(Square_points$tmpvec), function(x){st_cast(st_combine(x), "POLYGON")}))

# Removing internal lines between squares of the same category
# poly<-st_union(st_make_valid(poly), by_feature = T) # Too much memory used so step skipped

# Adding categories for landuse
Polygons<-cbind(Polygons, landuse=unique(Toess[,3:4])[,1])
Polygons$landuse<-as.character(Polygons$landuse)
Polygons<-sf::st_set_crs(Polygons, 2056)

# Get landuse for each point in a dataframe
match_landuse<-st_intersection(pts, Polygons, sparse=F)
# Removing double lines for points matched with two squares
match_landuse<-unique(as.data.frame(match_landuse)[,c(1,3)])
# !!!!!!!!!!!!!!!!!!! T143 is matched with both agriculture and forest 
match_landuse<-match_landuse[!(match_landuse$SiteCode=="T143" & match_landuse$landuse==2),]

##############################
# Diversity and Distribution
##############################
tab_map<-unique(Samples, by="SiteCode")

# Load map files
files <- list.files("~/SHAPEFILE_LV95_LN02", full.names=T) # what is in the donwloaded folder?
shp_file <- list.files("~/SHAPEFILE_LV95_LN02", pattern="shp", full.names=T) # retrieve filename for the shapefile 
border <- sf::st_read(shp_file) # read shapefile. ATTENTION: To work it needs all the asociated files, not only the shapefile itself but also shx, dbf, etc

# only select Swiss border
CH <- border[1,] # Only select Swiss parts. The initial file contains four objects, including the borders of Liechtenstein and enclaves of Germany (Buesingen) and Italy (Campione)

# create geometry object without attribute table
CH_geom <- sf::st_geometry(CH) # export the geometry part of the imported shapefile

# Small map of Switzerland
p1<-ggplot()+
  geom_sf(data=CH_geom) +
  theme(panel.background=element_rect(fill="white", colour="black"), axis.line=element_line(colour="black")) +
  theme(axis.text.x=element_text(size=14), axis.text.y=element_text(size=14), axis.title=element_blank()) +
  geom_point(data = tab_map, aes(x=Xcoord+2000000, y=Ycoord+1000000), color="indianred3")

pdf("Map_Switerland.pdf")
p1
dev.off()

# Map of Töss
CH <- CHGeoData:::CH
inn <- CHGeoData:::inn
rhone_VS <- CHGeoData:::rhone_VS
rhone_JU <- CHGeoData:::rhone_JU
rhine <- CHGeoData:::rhine
ri <- CHGeoData:::ri
lk <- CHGeoData:::lk

p2<-ggplot() +
  geom_polygon(data=CH, aes(x=x, y=y), fill="lightgrey") +
  geom_polygon(data=lk, aes(x=x, y=y), fill="blue") +
  geom_path(data=ri, aes(x=x, y=y, group=river_nr), color="blue") +
  scale_x_continuous(limits = c(681000, 717000), expand=c(0,0)) +
  scale_y_continuous(limits = c(240000, 272000), expand=c(0,0)) +
  geom_rect(mapping=aes(xmin=681000, xmax=717000, ymin=240000, ymax=272000), color="black", fill=NA, size=1.2) +
  theme(panel.background=element_rect(fill="white", colour="black"), axis.line=element_line(colour="black")) +
  theme(axis.text.x=element_text(size=14), axis.text.y=element_text(size=14), axis.title=element_blank()) +
  geom_point(data = tab_map[tab_map$Org_found=="yes",], aes(x=Xcoord, y=Ycoord), size=4, shape=21, fill="indianred3", color="black") +
  geom_point(data = tab_map[tab_map$Org_found=="no",], aes(x=Xcoord, y=Ycoord), size=4, shape=1, color="black") 

pdf("Map_Toss_Niph.pdf")  
p2
dev.off()

# Small map of Toess for aquifer types with sampled points
p1<-ggplot()+
  geom_sf(data=aquifer, aes(fill=aquifertyp)) +
  theme(panel.background=element_rect(fill="white", colour="black"), axis.line=element_line(colour="black")) +
  theme(axis.text.x=element_text(size=14), axis.text.y=element_text(size=14), axis.title=element_blank()) +
  scale_fill_manual(values=c("lightcoral", "lightgoldenrod1", "mediumpurple", "lightsteelblue1", "cornflowerblue", "darkseagreen2")) +
  geom_point(data = tab_map, aes(x=Xcoord+2000000, y=Ycoord+1000000), color="black") +
  scale_x_continuous(limits=c(2680000, 2720000)) +
  scale_y_continuous(limits=c(1240000, 1270000))

pdf("Map_Toess_aquifer.pdf")
  p1
dev.off()

# Small map of Toess for aquifer types with sampled points
p1bis<-ggplot()+
  geom_sf(data=Polygons, aes(fill=landuse), colour=NA) +
  theme(panel.background=element_rect(fill="white", colour="black"), axis.line=element_line(colour="black")) +
  theme(axis.text.x=element_text(size=14), axis.text.y=element_text(size=14), axis.title=element_blank()) +
  scale_fill_manual(values=c("indianred3", "lightgoldenrod1", "forestgreen", "cornflowerblue")) +
  geom_point(data = tab_map, aes(x=Xcoord, y=Ycoord), color="black") +
  scale_x_continuous(limits=c(2680000, 2720000)) +
  scale_y_continuous(limits=c(1240000, 1270000))

pdf("Map_Toess_landuse.pdf")
  p1bis
dev.off()


# Map of Töss
CH <- CHGeoData:::CH
inn <- CHGeoData:::inn
rhone_VS <- CHGeoData:::rhone_VS
rhone_JU <- CHGeoData:::rhone_JU
rhine <- CHGeoData:::rhine
ri <- CHGeoData:::ri
lk <- CHGeoData:::lk

p2<-ggplot() +
  geom_polygon(data=CH, aes(x=x, y=y), fill="lightgrey") +
  geom_polygon(data=lk, aes(x=x, y=y), fill="blue") +
  geom_path(data=ri, aes(x=x, y=y, group=river_nr), color="blue") +
  scale_x_continuous(limits = c(681000, 717000), expand=c(0,0)) +
  scale_y_continuous(limits = c(240000, 272000), expand=c(0,0)) +
  geom_rect(mapping=aes(xmin=681000, xmax=717000, ymin=240000, ymax=272000), color="black", fill=NA, size=1.2) +
  theme(panel.background=element_rect(fill="white", colour="black"), axis.line=element_line(colour="black")) +
  theme(axis.text.x=element_text(size=14), axis.text.y=element_text(size=14), axis.title=element_blank()) +
  geom_point(data = tab_map[tab_map$Org_found=="yes",], aes(x=Xcoord, y=Ycoord), size=4, shape=21, fill="indianred3", color="black") +
  geom_point(data = tab_map[tab_map$Org_found=="no",], aes(x=Xcoord, y=Ycoord), size=4, shape=1, color="black") 

pdf("Map_Toss_Niph.pdf")  
p2
dev.off()

# Map distribution per species (convex hulls)
tableau<-unique(Ind[Ind$Status=="Good" & Ind$Species!="Niphargus_sp." & Ind$Species!="Gammarus_fossarum_A" & Ind$Species!="Scutigerella_causeyae" & !(is.na(Ind$Species)),c(2,10,11)])
tableau<-sf::st_as_sf(tableau, coords = c("Xcoord", "Ycoord"), crs=2056)
hull<-st_sf(aggregate(tableau$geometry, list(tableau$Species), function(x){st_convex_hull(st_union(x))}))

pX<-ggplot() +
  geom_sf(data=hull, aes(fill=Group.1, colour=Group.1), alpha=0.5) +
  theme(panel.background=element_rect(fill="white", colour="black"), axis.line=element_line(colour="black")) +
  theme(axis.text.x=element_text(size=14), axis.text.y=element_text(size=14), axis.title=element_blank()) +
  scale_fill_manual(values=c("darkseagreen", "coral", "mediumorchid", "gray20", "limegreen", "hotpink", "goldenrod", "saddlebrown")) +
  scale_colour_manual(values=c("darkseagreen", "coral", "mediumorchid", "gray20", "limegreen", "hotpink", "goldenrod", "saddlebrown"))

pdf("Map_overlap_area.pdf")
  pX
dev.off()

# Barplot species distribution
tab_abund<-tab_distri[tab_distri$Status=="Good",]
tab_abund[,"Value":=sum(Mean), by=Species]
tab_abund<-unique(tab_abund, by="Species")[,c(3,17)]
tab_abund$Type<-"Abundance"

tab_sites<-tab_distri[tab_distri$Status=="Good",]
tab_sites[,"Value":=.N, by=Species]
tab_sites<-unique(tab_sites, by="Species")[,c(3,17)]
tab_sites$Type<-"#Sites"

tab_barplot<-rbind(tab_abund, tab_sites)
tab_barplot[tab_barplot$Type=="Abundance",]$Value <- -tab_barplot[tab_barplot$Type=="Abundance",]$Value
tab_barplot$Type<-factor(tab_barplot$Type, levels=c("Abundance", "#Sites"))
tab_barplot$Species<-factor(tab_barplot$Species, levels=rev(unique(tab_barplot[order(Value),]$Species)))
tab_barplot<-tab_barplot[tab_barplot$Species!="Scutigerella_causeyae" & tab_barplot$Species!="Niphargus_sp." & tab_barplot$Species!="Gammarus_fossarum_A" & !(is.na(tab_barplot$Species)),]

p3<-ggplot(data=tab_barplot, aes(x=Species, y=Value, fill=Type))+
  geom_bar(stat="identity") +
  facet_wrap(~Type, ncol=2, nrow=1, scales = "free_x") +
  coord_flip() +
  theme(panel.background=element_rect(fill="white", colour="black"), axis.line=element_line(colour="black"), legend.position = "none") +
  theme(axis.text.x=element_text(size=14), axis.text.y=element_text(size=14), axis.title=element_blank()) +
  theme(strip.background=element_rect(fill="white", colour="black"), strip.text=element_text(size=18, face="bold")) +
  scale_fill_manual(values=c("seagreen", "orange")) 

pdf("Barplot_distri_NiphSp_barcoding_correctedAbundances.pdf", width=15, height=8)
  p3
dev.off()

tab_barplot<-dcast(tab_barplot, Species~Type, fill=0, value.var = "Value")

p3bis<-ggplot(data=tab_barplot, aes(x=-Abundance, y=`#Sites`)) +
  geom_point(size=5) +
  theme(panel.background=element_rect(fill="white", colour="black"), axis.line=element_line(colour="black"), legend.position = "none") +
  theme(axis.text.x=element_text(size=14), axis.text.y=element_text(size=14), axis.title=element_text(size=18, face="bold")) +
  theme(strip.background=element_rect(fill="white", colour="black"), strip.text=element_text(size=18, face="bold"))

pdf("Corplot_distri_sites_NiphSp_barcoding_correctedAbundances.pdf")
  p3bis
dev.off()


# Barplot of species abundance against their distribution area
tab_abund<-tab_distri[tab_distri$Status=="Good",]
tab_abund[,"Value":=sum(Mean), by=Species]
tab_abund<-unique(tab_abund, by="Species")[,c(3,17)]
tab_abund$Type<-"Abundance"

colnames(tab_area)<-c("Species", "Value")
tab_area$Type<-"Area"

tab_barplot<-rbind(tab_abund, tab_area)
tab_barplot$Value<-as.numeric(tab_barplot$Value)
tab_barplot[tab_barplot$Type=="Abundance",]$Value <- -tab_barplot[tab_barplot$Type=="Abundance",]$Value
tab_barplot$Type<-factor(tab_barplot$Type, levels=c("Abundance", "Area"))
tab_barplot$Species<-factor(tab_barplot$Species, levels=rev(unique(tab_barplot[order(Value),]$Species)))
tab_barplot<-tab_barplot[tab_barplot$Species!="Scutigerella_causeyae" & tab_barplot$Species!="Niphargus_sp." & tab_barplot$Species!="Gammarus_fossarum_A" & !(is.na(tab_barplot$Species)) & tab_barplot$Species!="Unknown",]

p3<-ggplot(data=tab_barplot, aes(x=Species, y=Value, fill=Type))+
  geom_bar(stat="identity") +
  facet_wrap(~Type, ncol=2, nrow=1, scales = "free_x") +
  coord_flip() +
  theme(panel.background=element_rect(fill="white", colour="black"), axis.line=element_line(colour="black"), legend.position = "none") +
  theme(axis.text.x=element_text(size=14), axis.text.y=element_text(size=14), axis.title=element_blank()) +
  theme(strip.background=element_rect(fill="white", colour="black"), strip.text=element_text(size=18, face="bold")) +
  scale_fill_manual(values=c("seagreen", "orange")) 

pdf("Barplot_distriArea_NiphSp_barcoding_correctedAbundances.pdf", width=15, height=8)
  p3
dev.off()

tab_barplot<-dcast(tab_barplot, Species~Type, fill=0, value.var = "Value")

p3bis<-ggplot(data=tab_barplot, aes(x=log(-Abundance+1), y=Area)) +
  geom_point(size=5) +
  theme(panel.background=element_rect(fill="white", colour="black"), axis.line=element_line(colour="black"), legend.position = "none") +
  theme(axis.text.x=element_text(size=14), axis.text.y=element_text(size=14), axis.title=element_text(size=18, face="bold")) +
  theme(strip.background=element_rect(fill="white", colour="black"), strip.text=element_text(size=18, face="bold"))

pdf("Corplot_distri_area_NiphSp_barcoding_correctedAbundances.pdf")
  p3bis
dev.off()


# Barplot species richness and histogramme of abundance
tab_hist<-tab_distri
tab_hist[,"N":=sum(Mean), by=SiteCode]
tab_hist<-unique(tab_hist, by="SiteCode")
tab_hist<-tab_hist[tab_hist$N!=0,]

p4<-ggplot(data=tab_hist, aes(x=N)) +
  geom_histogram(binwidth= 1, fill="indianred3", color="black", boundary=0) +
  theme(panel.background=element_rect(fill="white", colour="black"), axis.line=element_line(colour="black"), panel.grid.major = element_line("gray90")) +
  theme(legend.position = "none") +
  theme(axis.text.x=element_text(size=14), axis.text.y=element_text(size=14), axis.title=element_text(size=18, face="bold")) +
  xlab("\nCorrected abundance per site") + ylab("# sites\n") + 
  scale_x_continuous(breaks=seq(0,30,5))

pdf("Histogram_corrected_abundance_perSite.pdf")
  p4
dev.off()

tab_barplot<-unique(Ind[Ind$Status=="Good" & Ind$Species!="Niphargus_sp." & Ind$Species!="Gammarus_fossarum_A" & Ind$Species!="Scutigerella_causeyae" & !(is.na(Ind$Species)),], by=c("Species", "SiteCode"))
tab_barplot<-as.data.table(table(tab_barplot$SiteCode))

p5<-ggplot(data=tab_barplot, aes(x=N)) +
  geom_bar(fill="indianred3", color="black") +
  theme(panel.background=element_rect(fill="white", colour="black"), axis.line=element_line(colour="black"), panel.grid.major = element_line("gray90")) +
  theme(legend.position = "none") +
  theme(axis.text.x=element_text(size=14), axis.text.y=element_text(size=14), axis.title=element_text(size=18, face="bold")) +
  xlab("\nSpecies richness per site") + ylab("# sites\n") 

pdf("Barplot_richness_perSite.pdf")
  p5
dev.off()

# Barplot richness compared to aquifer type
tab_barplot<-unique(Ind[Ind$Status=="Good" & Ind$Species!="Niphargus_sp." & Ind$Species!="Gammarus_fossarum_A" & Ind$Species!="Scutigerella_causeyae" & !(is.na(Ind$Species)),], by=c("Species", "SiteCode"))
tab_barplot<-as.data.table(table(tab_barplot$SiteCode))
colnames(tab_barplot)<-c("SiteCode", "Richness")
tab_barplot<-merge(tab_barplot, match_aquifer, by="SiteCode", all=T)
tab_barplot[is.na(tab_barplot$Richness),]$Richness<-0

p5<-ggplot(data=tab_barplot, aes(x=Richness, fill=AquiferType)) +
  geom_bar(position="stack", color="black") +
  theme(panel.background=element_rect(fill="white", colour="black"), axis.line=element_line(colour="black"), panel.grid.major = element_line("gray90")) +
  theme(legend.position = "none") +
  theme(axis.text.x=element_text(size=14), axis.text.y=element_text(size=14), axis.title=element_text(size=18, face="bold")) +
  xlab("\nSpecies richness per site") + ylab("# sites\n") +
  scale_fill_manual(values=c("cornflowerblue", "darkseagreen2"), labels=c("Unconsolidated", "Fissured"))

pdf("Barplot_richness_perSite_aquifer.pdf")
p5
dev.off()

# Barplot richness compared to landuse type
tab_barplot<-unique(Ind[Ind$Status=="Good" & Ind$Species!="Niphargus_sp." & Ind$Species!="Gammarus_fossarum_A" & Ind$Species!="Scutigerella_causeyae" & !(is.na(Ind$Species)),], by=c("Species", "SiteCode"))
tab_barplot<-as.data.table(table(tab_barplot$SiteCode))
colnames(tab_barplot)<-c("SiteCode", "Richness")
tab_barplot<-merge(tab_barplot, match_landuse, by="SiteCode", all=T)
tab_barplot[is.na(tab_barplot$Richness),]$Richness<-0

p5<-ggplot(data=tab_barplot, aes(x=Richness, fill=landuse)) +
  geom_bar(position="stack", color="black") +
  theme(panel.background=element_rect(fill="white", colour="black"), axis.line=element_line(colour="black"), panel.grid.major = element_line("gray90")) +
  theme(legend.position = "none") +
  theme(axis.text.x=element_text(size=14), axis.text.y=element_text(size=14), axis.title=element_text(size=18, face="bold")) +
  xlab("\nSpecies richness per site") + ylab("# sites\n") +
  scale_fill_manual(values=c("indianred3", "lightgoldenrod1", "forestgreen", "cornflowerblue"), labels=c("Urban", "Agriculture", "Forest", "Unproductive"))

pdf("Barplot_richness_perSite_landuse.pdf")
p5
dev.off()

# Comparing the overlap of distribution with the number of co-occurring sites

# Distribution over aquifer types
tab_barplot<-Ind[Ind$Status=="Good",]
tab_barplot<-merge(tab_barplot, match_aquifer, by="SiteCode", all.x=T, all.y=F)
tab_barplot[,"Abund_aquifer":=.N, by=.(Species, AquiferType)]
tab_barplot<-unique(tab_barplot, by=c("Species", "AquiferType"))
tab_barplot<-tab_barplot[tab_barplot$Species!="Scutigerella_causeyae" & tab_barplot$Species!="Niphargus_sp." & tab_barplot$Species!="Gammarus_fossarum_A" & !(is.na(tab_barplot$Species)),]
tab_barplot$Species<-factor(tab_barplot$Species, levels=c("Niphargus_tonywhitteni", "Niphargus_fontanus", "Niphargus_auerbachi", "Niphargus_puteanus", "Crangonyx_cf._subterraneus", "Niphargus_arolaensis", "Niphargus_thienemanni", "Niphargus_cf._thienemanni_2"))

pX<-ggplot(data=tab_barplot, aes(x=Species, y=Abund_landuse, fill=landuse)) +
  geom_bar(stat="identity", position="stack", color="black") +
  coord_flip() +
  theme(panel.background=element_rect(fill="white", colour="black"), axis.line=element_line(colour="black"), panel.grid.major = element_line("gray90")) +
  theme(legend.position = "top") +
  theme(axis.text.x=element_text(size=14), axis.text.y=element_text(size=14), axis.title=element_text(size=18, face="bold")) +
  scale_fill_manual(values=c("indianred3", "lightgoldenrod1", "forestgreen", "cornflowerblue"), labels=c("Urban", "Agriculture", "Forest", "Unproductive"))

pdf("Barplot_landuse_distribution.pdf", width=10, height=8)
  pX
dev.off()

# Boxplot of landuse types
tab_barplot<-tab_distri[tab_distri$Status!="Unsure",]
tab_barplot<-merge(tab_barplot, match_landuse, by="SiteCode", all.x=T, all.y=F)
tab_barplot[,"Somme":=sum(Mean), by=.(SiteCode, landuse)]
tab_barplot<-unique(tab_barplot, by=c("SiteCode", "landuse"))

pX<-ggplot(data=tab_barplot, aes(x=landuse, y=log(Mean + 1), fill=landuse)) +
  geom_boxplot(color="black") +
  theme(panel.background=element_rect(fill="white", colour="black"), axis.line=element_line(colour="black"), panel.grid.major = element_line("gray90")) +
  theme(legend.position = "top") +
  theme(axis.text.x=element_text(size=14), axis.text.y=element_text(size=14), axis.title=element_text(size=18, face="bold")) +
  scale_fill_manual(values=c("indianred3", "lightgoldenrod1", "forestgreen", "cornflowerblue"), labels=c("Urban", "Agriculture", "Forest", "Unproductive"))

pdf("Boxplot_distri_landuse.pdf")
  pX
dev.off()

# Distribution over landuse types
tab_barplot<-Ind[Ind$Status=="Good",]
tab_barplot<-merge(tab_barplot, match_landuse, by="SiteCode", all.x=T, all.y=F)
tab_barplot[,"Abund_landuse":=.N, by=.(Species, landuse)]
tab_barplot<-unique(tab_barplot, by=c("Species", "landuse"))
tab_barplot<-tab_barplot[tab_barplot$Species!="Scutigerella_causeyae" & tab_barplot$Species!="Niphargus_sp." & tab_barplot$Species!="Gammarus_fossarum_A" & !(is.na(tab_barplot$Species)),]
tab_barplot$Species<-factor(tab_barplot$Species, levels=c("Niphargus_tonywhitteni", "Niphargus_fontanus", "Niphargus_auerbachi", "Niphargus_puteanus", "Crangonyx_cf._subterraneus", "Niphargus_arolaensis", "Niphargus_thienemanni", "Niphargus_cf._thienemanni_2"))

pX<-ggplot(data=tab_barplot, aes(x=Species, y=Abund_landuse, fill=landuse)) +
  geom_bar(stat="identity", position="stack", color="black") +
  coord_flip() +
  theme(panel.background=element_rect(fill="white", colour="black"), axis.line=element_line(colour="black"), panel.grid.major = element_line("gray90")) +
  theme(legend.position = "top") +
  theme(axis.text.x=element_text(size=14), axis.text.y=element_text(size=14), axis.title=element_text(size=18, face="bold")) +
  scale_fill_manual(values=c("cornflowerblue", "darkseagreen2"), labels=c("Unconsolidated", "Fissured"))

pdf("Barplot_aquifer_distribution.pdf", width=10, height=8)
pX
dev.off()

# Boxplot of aquifer types
tab_barplot<-tab_distri[tab_distri$Status!="Unsure",]
tab_barplot<-merge(tab_barplot, match_aquifer, by="SiteCode", all.x=T, all.y=F)
tab_barplot[,"Somme":=sum(Mean), by=.(SiteCode, AquiferType)]
tab_barplot<-unique(tab_barplot, by=c("SiteCode", "AquiferType"))

pX<-ggplot(data=tab_barplot, aes(x=AquiferType, y=log(Mean + 1), fill=AquiferType)) +
  geom_boxplot(color="black") +
  theme(panel.background=element_rect(fill="white", colour="black"), axis.line=element_line(colour="black"), panel.grid.major = element_line("gray90")) +
  theme(legend.position = "top") +
  theme(axis.text.x=element_text(size=14), axis.text.y=element_text(size=14), axis.title=element_text(size=18, face="bold")) +
  scale_fill_manual(values=c("cornflowerblue", "darkseagreen2"), labels=c("Unconsolidated", "Fissured"))

pdf("Boxplot_distri_aquifer.pdf")
pX
dev.off()

#############################
# Correlation with diversity
#############################
# Presence-Absence contingency tables
Conting_barcod_PA<-as.matrix(Conting_barcod[,-1])
Conting_barcod_PA[Conting_barcod_PA>0]<-1
row.names(Conting_barcod_PA)<-Conting_barcod$Species

Conting_ASV_PA<-as.matrix(Conting_ASV[,-1])
Conting_ASV_PA[Conting_ASV_PA>0]<-1
row.names(Conting_ASV_PA)<-Conting_ASV$ASV

Conting_OTU_PA<-as.matrix(Conting_OTU[,-1])
Conting_OTU_PA[Conting_OTU_PA>0]<-1
row.names(Conting_OTU_PA)<-Conting_OTU$OTUs

# Correlation between amphipods and genetic diversity 
tab_graph<-as.data.frame(t(rbind(apply(Conting_ASV_PA, 2, sum), apply(Conting_barcod[,-1], 2, sum))))
colnames(tab_graph)<-c("ASVs", "Amphipods")
tab_graph$Code<-row.names(tab_graph)
tab_graph<-merge(tab_graph, Sites, by="Code", all=T)

Correlation_results<-cor.test(log(tab_graph$ASVs+1), log(tab_graph$Amphipods+1), method="pearson", alternative = "two.sided")

p6<-ggplot(data=tab_graph, aes(x=log(Amphipods+1), y=log(ASVs+1))) +
  geom_point(size=3) +
  geom_smooth(method="glm", formula="y~x") +
  theme(panel.background=element_rect(fill="white", colour="black"), axis.line=element_line(colour="black"), panel.grid.major = element_line("gray90")) +
  theme(axis.text.x=element_text(size=14), axis.text.y=element_text(size=14), axis.title=element_text(size=18, face="bold")) +
  xlab("\nlog(#Amphipods + 1)") + ylab("log(#ASVs + 1)\n") +
  geom_text(aes(x=0.5, y=8), label=paste0("r = ", round(Correlation_results$estimate, 3))) +
  geom_text(aes(x=0.5, y=7.8), label=paste0("p = ", round(Correlation_results$p.value, 3)))

p6bis<-ggplot(data=tab_graph, aes(x=log(Amphipods+1), y=log(ASVs+1))) +
  geom_point(aes(color=Land, shape=Aquifer), size=4) +
  theme(panel.background=element_rect(fill="white", colour="black"), axis.line=element_line(colour="black"), panel.grid.major = element_line("gray90")) +
  theme(axis.text.x=element_text(size=14), axis.text.y=element_text(size=14), axis.title=element_text(size=18, face="bold")) +
  theme(legend.position="top", legend.text=element_text(size=14), legend.title=element_text(size=18, face="bold")) +
  xlab("\nlog(#Amphipods + 1)") + ylab("log(#ASVs + 1)\n") +
  scale_color_manual(values=c("Purple4", "Darkorange")) +
  scale_shape_manual(values=c(16,2))

tab_graph<-as.data.frame(t(rbind(apply(Conting_ASV_PA, 2, sum), apply(Conting_barcod_PA[-c(1,9),], 2, sum))))  
colnames(tab_graph)<-c("ASVs", "Amphipods")
tab_graph$Code<-row.names(tab_graph)
tab_graph<-merge(tab_graph, Sites, by="Code", all=T)

Correlation_results<-cor.test(log(tab_graph$ASVs), tab_graph$Amphipods, method="pearson", alternative = "two.sided")

p7<-ggplot(data=tab_graph, aes(x=Amphipods, y=log(ASVs+1))) +
  geom_point(size=3) +
  geom_smooth(method="glm", formula="y~x") +
  theme(panel.background=element_rect(fill="white", colour="black"), axis.line=element_line(colour="black"), panel.grid.major = element_line("gray90")) +
  theme(axis.text.x=element_text(size=14), axis.text.y=element_text(size=14), axis.title=element_text(size=18, face="bold")) +
  xlab("\n# Species") + ylab("log(#ASVs)\n") +
  geom_text(aes(x=0.5, y=8), label=paste0("r = ", round(Correlation_results$estimate, 3))) +
  geom_text(aes(x=0.5, y=7.8), label=paste0("p = ", round(Correlation_results$p.value, 3)))

p7bis<-ggplot(data=tab_graph, aes(x=Amphipods, y=log(ASVs+1))) +
  geom_point(aes(color=Land, shape=Aquifer), size=4) +
  theme(panel.background=element_rect(fill="white", colour="black"), axis.line=element_line(colour="black"), panel.grid.major = element_line("gray90")) +
  theme(axis.text.x=element_text(size=14), axis.text.y=element_text(size=14), axis.title=element_text(size=18, face="bold")) +
  xlab("\n# Species") + ylab("log(#ASVs)\n") +
  theme(legend.position="top", legend.text=element_text(size=14), legend.title=element_text(size=18, face="bold")) +
  scale_color_manual(values=c("Purple4", "Darkorange")) +
  scale_shape_manual(values=c(16,2))

pdf("Correlation_amphipods_genetic_diversity_corrected.pdf", width=9, height=9)
  p6
  p6bis
  p7
  p7bis
dev.off()

# Correlation between amphipods and taxonomic diversity 
tab_graph<-as.data.frame(t(rbind(apply(Conting_OTU_PA, 2, sum), apply(Conting_barcod[,-1], 2, sum))))
colnames(tab_graph)<-c("OTUs", "Amphipods")
tab_graph$Code<-row.names(tab_graph)
tab_graph<-merge(tab_graph, Sites, by="Code", all=T)

Correlation_results<-cor.test(log(tab_graph$OTUs+1), log(tab_graph$Amphipods+1), method="pearson", alternative = "two.sided")

p6<-ggplot(data=tab_graph, aes(x=log(Amphipods+1), y=log(OTUs+1))) +
  geom_point(size=3) +
  geom_smooth(method="glm", formula="y~x") +
  theme(panel.background=element_rect(fill="white", colour="black"), axis.line=element_line(colour="black"), panel.grid.major = element_line("gray90")) +
  theme(axis.text.x=element_text(size=14), axis.text.y=element_text(size=14), axis.title=element_text(size=18, face="bold")) +
  xlab("\nlog(#Amphipods + 1)") + ylab("log(#OTUs + 1)\n") +
  geom_text(aes(x=0.5, y=8), label=paste0("r = ", round(Correlation_results$estimate, 3))) +
  geom_text(aes(x=0.5, y=7.8), label=paste0("p = ", round(Correlation_results$p.value, 3)))

p6bis<-ggplot(data=tab_graph, aes(x=log(Amphipods+1), y=log(OTUs+1))) +
  geom_point(aes(color=Land, shape=Aquifer), size=4) +
  theme(panel.background=element_rect(fill="white", colour="black"), axis.line=element_line(colour="black"), panel.grid.major = element_line("gray90")) +
  theme(axis.text.x=element_text(size=14), axis.text.y=element_text(size=14), axis.title=element_text(size=18, face="bold")) +
  theme(legend.position="top", legend.text=element_text(size=14), legend.title=element_text(size=18, face="bold")) +
  xlab("\nlog(#Amphipods + 1)") + ylab("log(#OTUs + 1)\n") +
  scale_color_manual(values=c("Purple4", "Darkorange")) +
  scale_shape_manual(values=c(16,2))

tab_graph<-as.data.frame(t(rbind(apply(Conting_OTU_PA, 2, sum), apply(Conting_barcod_PA[-c(1,9),], 2, sum))))  
colnames(tab_graph)<-c("OTUs", "Amphipods")
tab_graph$Code<-row.names(tab_graph)
tab_graph<-merge(tab_graph, Sites, by="Code", all=T)

Correlation_results<-cor.test(log(tab_graph$OTUs), tab_graph$Amphipods, method="pearson", alternative = "two.sided")

p7<-ggplot(data=tab_graph, aes(x=Amphipods, y=log(OTUs+1))) +
  geom_point(size=3) +
  geom_smooth(method="glm", formula="y~x") +
  theme(panel.background=element_rect(fill="white", colour="black"), axis.line=element_line(colour="black"), panel.grid.major = element_line("gray90")) +
  theme(axis.text.x=element_text(size=14), axis.text.y=element_text(size=14), axis.title=element_text(size=18, face="bold")) +
  xlab("\n# Species") + ylab("log(#OTUs)\n") +
  geom_text(aes(x=0.5, y=8), label=paste0("r = ", round(Correlation_results$estimate, 3))) +
  geom_text(aes(x=0.5, y=7.8), label=paste0("p = ", round(Correlation_results$p.value, 3)))

p7bis<-ggplot(data=tab_graph, aes(x=Amphipods, y=log(OTUs+1))) +
  geom_point(aes(color=Land, shape=Aquifer), size=4) +
  theme(panel.background=element_rect(fill="white", colour="black"), axis.line=element_line(colour="black"), panel.grid.major = element_line("gray90")) +
  theme(axis.text.x=element_text(size=14), axis.text.y=element_text(size=14), axis.title=element_text(size=18, face="bold")) +
  xlab("\n# Species") + ylab("log(#OTUs)\n") +
  theme(legend.position="top", legend.text=element_text(size=14), legend.title=element_text(size=18, face="bold")) +
  scale_color_manual(values=c("Purple4", "Darkorange")) +
  scale_shape_manual(values=c(16,2))

pdf("Correlation_amphipods_taxonomic_diversity_corrected.pdf", width=9, height=9)
  p6
  p6bis
  p7
  p7bis
dev.off()

# Correlation with sampling effort
tab_graph<-as.data.frame(t(rbind(apply(Conting_ASV_PA, 2, sum), apply(Conting_barcod[,-1], 2, sum))))
colnames(tab_graph)<-c("ASVs", "Amphipods")
tab_graph$Code<-row.names(tab_graph)
tab_graph<-merge(tab_graph, Sites, by="Code", all=T)

Correlation_results<-cor.test(log(tab_graph$ASVs+1), log(tab_graph$Amphipods+1), method="pearson", alternative = "two.sided")

p6<-ggplot(data=tab_graph, aes(x=log(Amphipods+1), y=log(ASVs+1))) +
  geom_point(aes(size=Effort)) +
  theme(panel.background=element_rect(fill="white", colour="black"), axis.line=element_line(colour="black"), panel.grid.major = element_line("gray90")) +
  theme(axis.text.x=element_text(size=14), axis.text.y=element_text(size=14), axis.title=element_text(size=18, face="bold")) +
  xlab("\nlog(#Amphipods + 1)") + ylab("log(#ASVs + 1)\n") +
  geom_text(aes(x=0.5, y=8), label=paste0("r = ", round(Correlation_results$estimate, 3))) +
  geom_text(aes(x=0.5, y=7.8), label=paste0("p = ", round(Correlation_results$p.value, 3)))

tab_graph<-as.data.frame(t(rbind(apply(Conting_ASV_PA, 2, sum), apply(Conting_barcod_PA[-c(1,9),], 2, sum))))  
colnames(tab_graph)<-c("ASVs", "Amphipods")
tab_graph$Code<-row.names(tab_graph)
tab_graph<-merge(tab_graph, Sites, by="Code", all=T)

Correlation_results<-cor.test(log(tab_graph$ASVs), tab_graph$Amphipods, method="pearson", alternative = "two.sided")

p7<-ggplot(data=tab_graph, aes(x=Amphipods, y=log(ASVs+1))) +
  geom_point(aes(size=Effort)) +
  theme(panel.background=element_rect(fill="white", colour="black"), axis.line=element_line(colour="black"), panel.grid.major = element_line("gray90")) +
  theme(axis.text.x=element_text(size=14), axis.text.y=element_text(size=14), axis.title=element_text(size=18, face="bold")) +
  xlab("\n# Species") + ylab("log(#ASVs)\n") +
  geom_text(aes(x=0.5, y=8), label=paste0("r = ", round(Correlation_results$estimate, 3))) +
  geom_text(aes(x=0.5, y=7.8), label=paste0("p = ", round(Correlation_results$p.value, 3)))

pdf("Correlation_with_sampling_effort.pdf")
  p6
  p7
dev.off()

#####################################################
# Comparison between taxonomic and genetic diversity
#####################################################
tab_graph<-as.data.frame(t(rbind(apply(Conting_OTU_PA, 2, sum), apply(Conting_ASV_PA, 2, sum))))  
colnames(tab_graph)<-c("OTUs", "ASVs")
tab_graph$Code<-row.names(tab_graph)
tab_graph<-merge(tab_graph, Sites, by="Code", all=T)

p8<-ggplot(data=tab_graph, aes(x=OTUs, y=ASVs)) +
  geom_point(aes(color=Land, shape=Aquifer), size=4) +
  theme(panel.background=element_rect(fill="white", colour="black"), axis.line=element_line(colour="black"), panel.grid.major = element_line("gray90")) +
  theme(axis.text.x=element_text(size=14), axis.text.y=element_text(size=14), axis.title=element_text(size=18, face="bold")) +
  theme(legend.position="top", legend.text=element_text(size=14), legend.title=element_text(size=18, face="bold")) +
  xlab("\n#OTUs") + ylab("#ASVs\n") +
  scale_color_manual(values=c("Purple4", "Darkorange")) +
  scale_shape_manual(values=c(16,2))

#########################
# Co-occurrence analyses
#########################
# Co-occurence between amphipod species presence and ASVs
  # Adding a total presence row on the barcoding contingency table (for looking at co-occurence between amphipods, whatever the species, and ASVs)
Conting_barcod_PA<-rbind(Conting_barcod_PA, apply(Conting_barcod_PA, 2, sum))
row.names(Conting_barcod_PA)[10]<-"Total"
Conting_barcod_PA[Conting_barcod_PA>0]<-1

# Calculation of co-occurences, personal method
CoOccur_PA<-CoOccurrencePA(eDNA_PA=Conting_ASV_PA, Amphipods_PA=Conting_barcod_PA)

# Heatmap showing the results
  # Sélection of only ASVs with a significant co-occurrence with at least one species
tab_pValues<-as.data.frame(CoOccur_PA$pValues[which(apply(CoOccur_PA$pValues, 1, min)<=0.05),])
tab_pValues$ASVs<-row.names(CoOccur_PA$pValues[which(apply(CoOccur_PA$pValues, 1, min)<=0.05),])
tab_pValues<-reshape2::melt(tab_pValues, value.name="pValue", variable.name="AmphiSP", id.vars="ASVs")

tab_graph<-as.data.frame(CoOccur_PA$match[which(apply(CoOccur_PA$pValues, 1, min)<=0.05),])
tab_graph$ASVs<-row.names(CoOccur_PA$match[which(apply(CoOccur_PA$pValues, 1, min)<=0.05),])
tab_graph<-reshape2::melt(tab_graph, value.name="Match", variable.name="AmphiSP", id.vars="ASVs")

tab_graph<-merge(tab_graph, tab_pValues, by=c("AmphiSP", "ASVs"))

p9<-ggplot(tab_graph, aes(x=AmphiSP, y=ASVs, fill=Match)) +
  geom_tile(color="black") +
  theme(panel.background=element_rect(fill="white", colour="black"), axis.line=element_line(colour="black")) +
  geom_text(aes(label=round(pValue,3))) +
  scale_fill_gradient(low="beige", high="darkred")

pdf("Heatmap_coOcur_PA_ASVs.pdf", width=10, height=30)
  p9
dev.off()

# Co-occurence between amphipod species presence and OTUs
# Calculation of co-occurences
CoOccur_PA<-CoOccurrencePA(eDNA_PA=Conting_OTU_PA, Amphipods_PA=Conting_barcod_PA)

# Heatmap showing the results
# Sélection of only ASVs with a significant co-occurrence with at least one species
tab_pValues<-as.data.frame(CoOccur_PA$pValues[which(apply(CoOccur_PA$pValues, 1, min)<=0.05),])
tab_pValues$OTUs<-row.names(CoOccur_PA$pValues[which(apply(CoOccur_PA$pValues, 1, min)<=0.05),])
tab_pValues<-reshape2::melt(tab_pValues, value.name="pValue", variable.name="AmphiSP", id.vars="OTUs")

tab_graph<-as.data.frame(CoOccur_PA$match[which(apply(CoOccur_PA$pValues, 1, min)<=0.05),])
tab_graph$OTUs<-row.names(CoOccur_PA$match[which(apply(CoOccur_PA$pValues, 1, min)<=0.05),])
tab_graph<-reshape2::melt(tab_graph, value.name="Match", variable.name="AmphiSP", id.vars="OTUs")

tab_graph<-merge(tab_graph, tab_pValues, by=c("AmphiSP", "OTUs"))

p10<-ggplot(tab_graph, aes(x=AmphiSP, y=OTUs, fill=Match)) +
  geom_tile(color="black") +
  theme(panel.background=element_rect(fill="white", colour="black"), axis.line=element_line(colour="black")) +
  geom_text(aes(label=round(pValue,3))) +
  scale_fill_gradient(low="beige", high="darkred")

pdf("Heatmap_coOcur_PA_OTUs.pdf", width=10, height=30)
  p10
dev.off()

##########################
# Test of package cooccur
##########################
# Creating the matrix for input
mat_cooccur<-rbind(Conting_barcod_PA, Conting_ASV_PA)
# Calculation
Res_cooccur<-cooccur(mat_cooccur, thresh=T, spp_names=T, true_rand_classifier=0.1, prob="comb")
# Getting only significant results
Sign_cooccur<-Res_cooccur$results[Res_cooccur$results$p_lt<0.05 | Res_cooccur$results$p_gt<0.05,]

# Creating network
  # Creating adjacency matrix
Pos_cooccur<-Sign_cooccur[Sign_cooccur$p_gt<=0.05,]
Neg_cooccur<-Sign_cooccur[Sign_cooccur$p_lt<=0.05,]

AdjMat_Pos<-matrix(0, nrow=length(unique(c(Pos_cooccur$sp1, Pos_cooccur$sp2))), ncol=length(unique(c(Pos_cooccur$sp1, Pos_cooccur$sp2))))
colnames(AdjMat_Pos)<-sort(unique(c(Pos_cooccur$sp1, Pos_cooccur$sp2)))
rownames(AdjMat_Pos)<-sort(unique(c(Pos_cooccur$sp1, Pos_cooccur$sp2)))

AdjMat_Neg<-matrix(0, nrow=length(unique(c(Neg_cooccur$sp1, Neg_cooccur$sp2))), ncol=length(unique(c(Neg_cooccur$sp1, Neg_cooccur$sp2))))
colnames(AdjMat_Neg)<-sort(unique(c(Neg_cooccur$sp1, Neg_cooccur$sp2)))
rownames(AdjMat_Neg)<-sort(unique(c(Neg_cooccur$sp1, Neg_cooccur$sp2)))

for(i in 1:nrow(Pos_cooccur)){
  AdjMat_Pos[as.character(Pos_cooccur$sp1[i]), as.character(Pos_cooccur$sp2[i])]<-Pos_cooccur$obs_cooccur[i]
}; rm(i)

for(i in 1:nrow(Neg_cooccur)){
  AdjMat_Neg[as.character(Neg_cooccur$sp1[i]), as.character(Neg_cooccur$sp2[i])]<-(Neg_cooccur$sp1_inc[i] + Neg_cooccur$sp2_inc[i] - Neg_cooccur$obs_cooccur[i])/(Neg_cooccur$sp1_inc[i] + Neg_cooccur$sp2_inc[i])
}; rm(i)

Mat.CO.Pos.binary = (Mat.CO.Pos > 0) + 0
Mat.CO.Neg.binary = (Mat.CO.Neg > 0) + 0
#####

# Correlation plot to find a detection threshold
tab_graph=tab_comp[tab_comp$Type=="AllDates",]

p3=ggplot(tab_graph, aes(x=log10(CountSite+1), y=log10(CountMetab+1), color=species)) +
  geom_point(size=4) +
  theme(panel.background=element_rect(fill="white", colour="black"), axis.line=element_line(colour="black")) +
  theme(legend.position="top", legend.text=element_text(size=14), legend.title=element_blank()) +
  theme(axis.text=element_text(size=14), axis.title.y=element_blank(), axis.title.x=element_text(size=16, face="bold")) +
  scale_color_manual(values=c("darkred", "darkcyan", "peachpuff4", "royalblue4", "khaki", "cornflowerblue"))

pdf("Corplot_AllDates.pdf")
	p3
dev.off()

