#####################################
# Likelihood Ratio Test using DE-Seq2
#####################################

if (!requireNamespace("BiocManager", quietly = TRUE))
install.packages("BiocManager")

BiocManager::install("DESeq2")

library(DESeq2)


# Prepare the data 
#-----------------

# import the count table
data <- read.table("countDF_100_samples.txt", header = T)
row.names(data) <- data$Gene
data <- data[,-1]

# remove genes with total read number < 20 to increase the speed 
# Note: these lowly expressed genes were also removed for the Metacycle analysis and were not identified at differentially expressed in pairwise comparisons
data$Reads <- rowSums(data)
data <- data[data$Reads > 20,]
which(colnames(data)=="Reads")
data <- data[,-101]

# Note: one gene generated an error when performing the LRT. This gene has therefore been removed from the dataset
# remove the gene making error during LRT (see L78)
data <- data[-3699,]

# import the metadata file
metadata <- read.table("metadata.txt", header = T)

# We need to remove the time point T01 because it corresponds to ZT72 and there is data at 22C only 
# ZT72 was used for the time course at 22C but not to study the heat stress response
which(colnames(data)==c("TOT.T01.C1")) #73
which(colnames(data)==c("TOT.T01.C2")) #96
which(colnames(data)==c("TOT.T01.C3")) #92
which(colnames(data)==c("TR.T01.C1")) #98
which(colnames(data)==c("TR.T01.C2")) #54
which(colnames(data)==c("TR.T01.C3")) #93

data_AOV <- data[,-c(73,96,92,98,54,93)]

# Do the same in the metadata file
condition <- metadata[!(metadata$SampleName=="TOT.T01.C1"|
                          metadata$SampleName=="TOT.T01.C2"|
                          metadata$SampleName=="TOT.T01.C3"|
                          metadata$SampleName=="TR.T01.C1"|
                          metadata$SampleName=="TR.T01.C2"|
                          metadata$SampleName=="TR.T01.C3"),]

# We need to have the same order between columns in count data and target file, so you should check it before
row.names(condition) <- condition$SampleName
all(rownames(condition) %in% colnames(data_AOV))#If TRUE, row names of 'condition' and col names of 'data_AOV' match

# Transform the Time as a Factor
condition$Time <- factor(condition$Time)


# Build the model
#################
dds <- DESeqDataSetFromMatrix(countData = data_AOV,
                              colData = condition,
                              design= ~ Temperature + RNA + Time + Temperature:RNA + Temperature:Time)

# test the model to view it
t <- model.matrix(~ Temperature + RNA + Time + Temperature:RNA + Temperature:Time, condition)
colnames(t)


# run the model to test the effect of
#####################################
#Temperature:Time
dds_LRT_temp_time<- DESeq(dds, test="LRT", reduced = ~ Temperature + RNA + Time + Temperature:RNA)

#Temperature:RNA
dds_LRT_temp_rna <- DESeq(dds, test='LRT', reduced = ~ Temperature + RNA + Time + Temperature:Time)

#Temperature
dds_LRT_temp <- DESeq(dds, test='LRT', reduced = ~ RNA + Time + Tempeature:RNA + Temperature:Time)


## if get error "1 rows did not converge in beta, labelled in mcols(object)$fullBetaConv. Use larger maxit argument with nbinomLRT"
# Try to delete the problematic genes from raw dataset (they are "typically genes with very small counts and little power") see L26
which((mcols(dds_LRT_temp_time)$fullBetaConv)=="FALSE")
#or increase number of iterations (see maxit argument in the nbinomLRT() function)

#obtain contrasts names
resultsNames(dds_LRT_temp_time)
resultsNames(dds_LRT_temp_rna)

#Comparisons for single factors to generate a result table
##########################################################
# temp:time
#----------
Temp_Time <- as.data.frame(results(dds_LRT_temp_time, name="Time_T3_vs_T0", alpha = 0.05))#with that line, obtain p-value per gene of the impact of the factor(s) that have been excluded from the LRT model; if redo that line and change contrast name, obtain same p-value, but different log2 fold change
res_test1 <- as.data.frame(results(dds_LRT_temp_time, name="TemperatureH.TimeT9", alpha = 0.05))#demonstration: this should show same p-value that above

Temp_Time$AGI <- row.names(Temp_Time)

# export Temp_Time
write.table(Temp_Time, "Temp_Time.txt", sep = "\t", row.names = F)

# temp:RNA
#---------
Temp_RNA <- as.data.frame(results(dds_LRT_temp_rna, name="RNA_TR_vs_TOT", alpha = 0.05))#with that line, obtain p-value per gene of the impact of the factor(s) that have been excluded from the LRT model; if redo that line and change contrast name, obtain same p-value, but different log2 fold change
res_test2 <- as.data.frame(results(dds_LRT_temp_rna, name="TemperatureH.RNATR", alpha = 0.05))#demonstration: this should show same p-value that above

Temp_RNA$AGI <- row.names(Temp_RNA)

# export Temp_Time
write.table(Temp_RNA, "Temp_RNA.txt", sep = "\t", row.names = F)


