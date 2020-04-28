library(ggplot2)

# import data summarized from BAM files ####

extractReadCoverage <- function(col){
  col <- sub("^[ACGT]:", "", col)
  col <- sub(":.*$", "", col)
  return(as.integer(col))
}

readBAMsummary <- function(file, poslookup = poslist, dir = folder){
  cc <- c("character", "integer", "numeric")[c(1, 2, 1, 1, 2, 2, 3, 1, 1, 1, 1)]
  mydata <- read.delim(file, colClasses = cc)
  
  # clean up file name to extract stuff from it
  file <- sub(folder, "", file, fixed = TRUE)
  # fill in sample name
  sample <- sub("^[[:upper:][:digit:]]*_[[:upper:][:digit:]]*-", "", file)
  sample <- sapply(strsplit(sample, "_"), function(x) paste(x[1], x[2], sep = "_"))
  mydata$Sample <- sample
  # fill in region name
  gene <- sub("^[[:upper:][:digit:]]*_[[:upper:][:digit:]]*-[[:digit:]]*_",
                "", file)
  gene <- sub("_[ACGT]*\\.BAM\\.summary\\.readcounts\\.txt$", "", gene)
  mydata$Gene <- gene
  region <- sub("-.*$", "", file)
  mydata$Region <- region
  
  # subset to desired positions
  positions <- poslookup[[region]]
  mydata <- mydata[mydata$POS %in% positions,]
  
  # extract read coverage
  mydata$A <- extractReadCoverage(mydata$Details.for.base.A)
  mydata$C <- extractReadCoverage(mydata$Details.for.base.C)
  mydata$G <- extractReadCoverage(mydata$Details.for.base.G)
  mydata$T <- extractReadCoverage(mydata$Details.for.base.T)
  
  # stack data frame
  stacked_coverage <- stack(mydata[,c("A", "C", "G", "T")])
  colnames(stacked_coverage) <- c("Coverage", "Base")
  mydataS <- cbind(mydata[rep(1:nrow(mydata), times = 4),
                          c("Gene", "Region", "Sample", "CHROM", "POS")],
                   stacked_coverage)
  
  return(mydataS)
}

# get all of the files
folder <- "results/counts4Lindsay/"
inputs <- list.files(path = folder, pattern = ".*\\.BAM\\.summary\\.readcounts\\.txt")
inputs <- paste0(folder, inputs)

# get a lookup for the region of interest
poslist <- list(BRAF584KF_BRAF584KR = 39639643:39639645,
                BRAF637VF_BRAF637VR = 39627782:39627784,
                EGFRF_EGFRR = 16869231:16869233,
                HRAS61QF_HRAS61QR = 141192549:141192551)

#test <- readBAMsummary(inputs[1])

# loop through files to make one data frame
allcoverage <- do.call(rbind, lapply(inputs, readBAMsummary))
rownames(allcoverage) <- NULL

# plots ####
ggplot(allcoverage,
       aes(x = POS, y = Coverage, fill = Base)) +
  geom_col(position = "fill", width = 1) +
  facet_grid(Sample ~ Region, scales = "free", space = "free") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1,size = 6)) +
  theme(axis.text.y.left = element_text(size = 5)) +
  theme(strip.text.y = element_text(angle = 0, size = 6)) +
  theme(strip.text.x = element_text(angle = 0, size = 6)) +
  labs(y = "Proportion coverage", x = "Position")
ggsave("anakk_codon_barplots3.pdf")
