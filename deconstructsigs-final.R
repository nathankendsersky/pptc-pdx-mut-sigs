#### install/require packages ######################################################################

if (!requireNamespace("devtools", quietly = TRUE)) install.packages("devtools")
library(devtools)

if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
library(BiocManager)

if (!requireNamespace("Biostrings", quietly = TRUE)) BiocManager::install("Biostrings")
if (!requireNamespace("BSgenome.Hsapiens.UCSC.hg19", quietly = TRUE)) BiocManager::install("BSgenome.Hsapiens.UCSC.hg19")
if (!requireNamespace("deconstructSigs", quietly = TRUE)) devtools::install_github(repo = "raerose01/deconstructSigs",dependencies = F)
if (!requireNamespace("ggthemes", quietly = TRUE)) install.packages("ggthemes")
if (!requireNamespace("ggplot2", quietly = TRUE)) install.packages("ggplot2")
if (!requireNamespace("tidyr", quietly = TRUE)) install.packages("tidyr")
if (!requireNamespace("Hmisc", quietly = TRUE)) install.packages("Hmisc")
if (!requireNamespace("gdata", quietly = TRUE)) install.packages("gdata")
if (!requireNamespace("RCurl", quietly = TRUE)) install.packages("RCurl")

toLoad <- c("Biostrings","deconstructSigs","BSgenome.Hsapiens.UCSC.hg19","ggthemes","ggplot2","tidyr","Hmisc","gdata","RCurl")
lapply(toLoad,require,character.only=TRUE)

#### run mutsigs analysis - cosmic ######################################################################

# set directories for saving files, specify histology of interest
setwd("~")
home <- "~/Desktop/mutsig-demo/"
hist <- "all"

# load pptc.merge file; 240 PDXs
load(paste0(home,"2019-02-14-allpdx-clean-maf-240.rda"), verbose = T)
sig.df <- pptc.merge[,c("Tumor_Sample_Barcode","Chromosome","Start_position","Reference_Allele","Tumor_Seq_Allele2")]

# require Sample, chr, pos, ref, alt
names(sig.df) <- c("Sample", "chr", "pos", "ref", "alt")
unique(sig.df$Sample)

# list of samples
samplelist <- as.list(unique(sig.df$Sample))

# convert to deconstructSigs input - warning message for samples with <50 mutations
sigs.input <- mut.to.sigs.input(mut.ref = sig.df[1:5],
                                sample.id = "Sample", 
                                chr = "chr", 
                                pos = "pos", 
                                ref = "ref", 
                                alt = "alt",
                                bsg = BSgenome.Hsapiens.UCSC.hg19)

# remove samples with less than 50 muts (24 PDXs)
remove <- c("ALL-58", "OS-36-SJ", "OS-32", "Rh-30R", "COG-N-453x", 
            "COG-N-452x", "Rh-30", "COG-N-623x", "ALL-80", "NCH-CA-1", 
            "ALL-105", "ICB-1744MB", "ALL-82", "IC-22909PNET-rIII", "ALL-102", 
            "COG-N-618x", "ALL-46", "NCH-WT-7", "NCH-WT-4","IC-6634GBM", 
            "COG-N-603x", "ICb-9850PNET", "MLL-5", "COG-N-619x")
new.df <- subset(sig.df, !(sig.df$Sample %in% remove))
unique(new.df$Sample)
samplelist <- as.list(unique(new.df$Sample)) # write over sold samplelist

# convert to deconstructSigs input without samples <50 muts, no warning
sigs.input <- mut.to.sigs.input(mut.ref = new.df, 
                                sample.id = "Sample", 
                                chr = "chr", 
                                pos = "pos", 
                                ref = "ref", 
                                alt = "alt",
                                bsg = BSgenome.Hsapiens.UCSC.hg19)
head(t(sigs.input))

# determine the signatures contributing to each samples
for (each in samplelist) {
  
  test = whichSignatures(tumor.ref = sigs.input, 
                         signatures.ref = signatures.cosmic,
                         sample.id = each, 
                         contexts.needed = TRUE,
                         tri.counts.method = 'exome')
  names(test)
  
  weights <- as.data.frame(test$weights)
  weights$Model <- rownames(weights)
  weights <- weights[,c(28,1:27)]
  names(weights)
  weights$Unknown <- test$unknown
  write.table(weights, paste(home, "signatures/bysample/", Sys.Date(), "-", each, "-signature-weights.txt",
                             sep = ""), sep = "\t", col.names = T, row.names = T, quote = F)
  
  write.table(test$tumor, paste(home,  "signatures/bysample/", Sys.Date(), "-", each, "-tumor-trinucleotide-context.txt",
                                sep = ""), sep = "\t", col.names = T, row.names = T, quote = F)
  write.table(test$product, paste(home,  "signatures/bysample/", Sys.Date(), "-", each, "-product-trinucleotide-context.txt",
                                  sep = ""), sep = "\t", col.names = T, row.names = T, quote = F)
  
  pdf(paste(home,"figures/signatures/", Sys.Date(), "-", each, "-signatures.pdf", sep = ""), width = 11, height = 8.5)
  plotSignatures(test)
  makePie(test)
  dev.off()
  
}

# rbind all dataframes
allweights <- dir(path = paste(home, "signatures/bysample/", sep= ""),pattern = "weights.txt", recursive=F, full.names=T)
comb.weights <- do.call("rbind",lapply(allweights, FUN=function(files){read.table(files,header=TRUE, sep="\t")}))
comb.weights$Model <- rownames(comb.weights)
comb.weights <- comb.weights[,c(ncol(comb.weights),1:(ncol(comb.weights)-1))] # reorganize headers
head(comb.weights)
write.table(comb.weights, paste(home, Sys.Date(), "-allpdx-signature-weights.txt", sep = ""), 
            sep = "\t", col.names = T, row.names = F, quote = F)

sort(colnames(comb.weights))

# generate PDF histograms of each model, all signatures
for (each in colnames(comb.weights[,2:ncol(comb.weights)])){
  pdf(paste(home,"/figures/signatures/", Sys.Date(), "-", each, "-histograms.pdf", sep = ""), width = 11, height = 8.5)
  hist(comb.weights[,each], breaks =20, main = each)
  dev.off()
}


#### plot signatures by in each model by histology ####

home <- "~/Desktop/mutsig-demo/"

# import comb.weights, CHANGE name of file
sig.weights.file <- system(paste0("ls ",home,"*-allpdx-signature-weights.txt"),intern = T)
comb.weights <- read.delim(sig.weights.file,header=T,as.is=T,sep="\t")

# read in clinical data
laml.clin <- read.delim(paste0(home,"pptc-pdx-clinical-web.txt"),header=T,as.is=T,sep="\t")

# subset samples in PPTC with maf file
clin.pptc <- subset(laml.clin, DNA.Part.of.PPTC == "yes" & Have.maf == "yes")
setdiff(clin.pptc$Model, comb.weights$Model) # these are the samples with <50 mutations
col.to.keep <- c("Model","Histology","Histology.Detailed","Histology.Onco.New")
clin.pptc <- clin.pptc[,col.to.keep]

# merge weights and clinical data
weights.clin <- merge(clin.pptc, comb.weights)
table(weights.clin$Histology.Detailed)

# reshape the data to plot
headers <- names(weights.clin) # set gather to first Signature.X : last Signature.X
data.long <- gather(weights.clin, Signature, Measurement, Signature.28:Unknown, factor_key=TRUE)
data.long[data.long == 0] <- NA
data.complete <- data.long[complete.cases(data.long),]

# generate loop criteria (Detailed)
histos <- unique(data.complete$Histology.Onco.New)
head(histos)

# load leukemia/brain and solid tumors for specific order
leuk.brain <- read.table(text=getURL("https://raw.githubusercontent.com/nathankendsersky/pptc-pdx-mut-sigs/master/leukemia-brain-order.txt"),as.is=T)[[1]]
solid <- read.table(text=getURL("https://raw.githubusercontent.com/nathankendsersky/pptc-pdx-mut-sigs/master/solid-order.txt"),as.is=T)[[1]]


# select leukemia/brain OR solid to run analysis
histo <- c("leukemia","brain")
order <- leuk.brain
# histo <- c("solid")
# order <- solid

# subset data based on histology
data.complete.sub <- subset(data.complete, data.complete$Histology.Onco.New %in% histo)
data.complete.sub1 <- subset(data.complete.sub, data.complete.sub$Measurement >= .1) # subset samples that have >0.1 cosine similarity value
data.complete.sub1 <- data.complete.sub1[data.complete.sub1$Signature != "Unknown",] # remove unknown signatures
data.complete.sub1$Model <- factor(data.complete.sub1$Model, levels = (order))
data.complete.sub1[order(data.complete.sub1$Model), ]

# load signature colors and categories dataframe
fig.cat.color <- read.table(text=getURL("https://raw.githubusercontent.com/nathankendsersky/pptc-pdx-mut-sigs/master/signature-category-color.txt"),
                            sep="\t",header = T, comment.char = "")
fig.cat.color.sub <- subset(fig.cat.color, fig.cat.color$Signature %in% unique(data.complete.sub1$Signature))
plot.colors <- as.character(fig.cat.color.sub$Colors)

## plot graph
source(paste0(home,"theme.R"))
temp.plotname = paste0(home, "figures/", Sys.Date(), "-", "HistologyOnc-10cutoff-labels-",paste(histo,collapse=""),"-individ-proportions.pdf")
pdf(temp.plotname, width = 23, height = 8)
temp.plot <- ggplot() + geom_bar(data=data.complete.sub1,aes(y=data.complete.sub1$Measurement, x=data.complete.sub1$Model, fill=data.complete.sub1$Signature), 
                                 stat="identity", position = "fill") + 
  scale_fill_manual(values = plot.colors) +
  theme_Publication() + ylab("Proportion of Signature") + xlab("Model") +
  labs(fill='Signatures') + theme(axis.text.x = element_text(angle = 90, hjust = 1))
print(temp.plot)
dev.off()



#### stop ####
