####################
# Loading libraries
####################

# Libraries
library(ggplot2)
library(data.table)
library(cooccur)
library(sf)

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

# END
