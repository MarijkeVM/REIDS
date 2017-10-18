# Processing the HTA CEL Files of the RASA paper
# Create folders
# rawData -> HTA_RASA -> hjay: put in CEL Files
# annotationData -> chipTypes -> hjay : put in CDF File

library(aroma.affymetrix)
verbose <- Arguments$getVerbose(-8, timestamp=TRUE)

chipType <- "hjay"
cdf <- AffymetrixCdfFile$byChipType(chipType)
print(cdf)


#Setting up CELset
cs <- AffymetrixCelSet$byName("HTA_RASA", cdf=cdf,verbose=TRUE)
samplenames=cs$names
setCdf(cs, cdf)

#Backgroundcorrection
bc <- RmaBackgroundCorrection(cs, tags="*,r1")
csBC <- process(bc,verbose=verbose)


#(Quantile) Normalization
qn <- QuantileNormalization(csBC, typesToUpdate="pm")
csN <- process(qn, verbose=verbose)


flattenCellIndices <- function(cells, ...) {
	# Flatten cell data
	cells <- unlist(cells) #ok
	
	# Do some tricks to clean up the names
	names(cells)<- gsub("[.](groups|indices)", "", names(cells))
	
	# Extract the vector of unit names
	unitNames <- gsub("[.].*", "", names(cells))
	
	# Extract the unit group names
	groupNames <- gsub(".*[.]", "", names(cells))
	
	# Merge data
	data.frame(unitNames=unitNames, groupNames=groupNames, cell=cells)
	
} 

cells1 <-getCellIndices(cdf, units=1:nbrOfUnits(cdf))
cells2 <-flattenCellIndices(cells1)

GeneID <- as.character(cells2$unitNames)
ExonID <- as.character(cells2$groupNames)
ProbeID<-as.character(cells2$cell)

probeintensities <- getData(csN, indices=cells2$cell, fields=c("intensities"))$intensities
probeintensities <- log2(probeintensities) 
colnames(probeintensities) <- cs$names

UniqueGeneID <- unique(GeneID)
UniqueExonID <- unique(ExonID)

HTAData_RASA=cbind(GeneID,ExonID,ProbeID,probeintensities)
HTAData_RASA=as.data.frame(HTAData_RASA)
head(HTAData_RASA)

save(HTAData_RASA,file="HTAData_RASA.RData")


# EXON LEVEL ESTIMATES
# NOTE: DATA NOT LOG2 TRANASFORMED
getCdf(csN)
plmEx <- ExonRmaPlm(csN, mergeGroups=FALSE)
print(plmEx)
fit(plmEx, verbose=verbose)

#### Getting the data out
cesEx <- getChipEffectSet(plmEx)
exFit <- extractDataFrame(cesEx, units=NULL, addNames=TRUE)
str(exFit)

HTA_ExonLevelSummarized_rma_RASA=exFit
colnames(HTA_ExonLevelSummarized_rma_RASA)[c(1,2)]=c("GeneID","ExonID")
head(HTA_ExonLevelSummarized_rma_RASA)

save(HTA_ExonLevelSummarized_rma_RASA,file="HTA_ExonLevelSummarized_rma_RASA.RData")

# GENELEVEL ESTIMATES
# NOTE: DATA NOT LOG2 TRANASFORMED
plmTr <- ExonRmaPlm(csN, mergeGroups=TRUE)
print(plmTr)
fit(plmTr, verbose=verbose)

#### Getting the data out
cesTr <- getChipEffectSet(plmTr)
TrFit <- extractDataFrame(cesTr, units=NULL, addNames=TRUE)
str(TrFit)

HTA_GeneLevelSummarized_rma_RASA=TrFit
colnames(HTA_GeneLevelSummarized_rma_RASA)[1]=c("GeneID")
head(HTA_GeneLevelSummarized_rma_RASA)

save(HTA_GeneLevelSummarized_rma_RASA,file="HTA_GeneLevelSummarized_rma_RASA.RData")