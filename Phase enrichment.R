
# This script shows how to analyze the phase enrichment from a list of DRGs, as compared to circadian total mRNAs

results_TOT <- read.table("total_cycling.txt", header = T)

phases_TOT <- results_TOT[,c(1,5)]
names(phases_TOT) <- c("AGI","LAG")

# replace the phase 25.5 by 1.5 and 24 by 0 to avoid repeating the times phases 0 and 1.5
library(stringr)
phases_TOT$LAG <- str_replace_all(phases_TOT$LAG, '25.5', '1.5')
phases_TOT$LAG <- str_replace_all(phases_TOT$LAG, '24', '0')

ls.str(phases_TOT)
phases_TOT$LAG <- as.numeric(phases_TOT$LAG)

# count the number of transcript by phase in the reference
library(plyr)
library(dplyr)
phases_TOT_Freq <- count(phases_TOT, vars = "LAG")
names(phases_TOT_Freq) <- c("Phase","TOT")

# Count the number of transcript by phase in the uploaded list
# This requires to import the list of genes of interest prior this analysis
phases_upload <- merge.data.frame(upload,phases_TOT, by= "AGI")
phases_upload <- phases_upload[,c(1,2)]

phases_upload_Freq <- count(phases_upload, vars = "LAG")
names(phases_upload_Freq) <- c("Phase", "Upload")

# merge tables with the phase information
phases_summary <- Reduce(function(x,y) merge(x = x, y = y, by = "Phase"), 
       list(phases_TOT_Freq,phases_upload_Freq))

# calculate number of genes that do not have the phase (calculated line by line)
phases_summary$TOT.not <- (sum(phases_summary$TOT)-phases_summary$TOT)
phases_summary$Upload.not <- (sum(phases_summary$Upload)-phases_summary$Upload)

# calculate phase enrichment
phases_summary$TOT.Prop <- (phases_summary$TOT)/(sum(phases_summary$TOT))
phases_summary$Upload.Prop <- (phases_summary$Upload)/(sum(phases_summary$Upload))

phases_summary$Enrich <- (phases_summary$Upload.Prop)/(phases_summary$TOT.Prop)


# Chi-square test to calculate significant differences of enrichment
list.y.var <- unique(phases_summary$Phase)

Upload = list()
result = list()

for(i in list.y.var)
{
  sub = phases_summary[phases_summary$Phase == i, c(2,4,3,5)]
  names(sub) <- rep(c("N","non_N"),2)
  sub = rbind(sub[,1:2],sub[,3:4])
  sub = as.matrix(sub)
  khi2 <- chisq.test(sub)
  result[[i]] <- chisq.test(sub)
  Upload[[i]] <- khi2$p.value
}


# create a table with the results
Chi_square <- as.matrix(Chi_square_phases_up)
