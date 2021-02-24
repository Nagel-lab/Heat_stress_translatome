
# ARS can not analyze unevenly sampled datasets, or evenly sampled datasets but with missing values, 
# or with replicate samples, or with non-integer sampling interval.
# In our datasets, we removed 2 samples after quality checking (one in Total mRNAs, one in TRAP mRNAs)

library("MetaCycle")
library("dplyr")
library("ggplot2")
library("cowplot")

# Metacycle total mRNAs
#######################
rlog_TOT <- read.csv("rlog_values_TOT_expressed.csv", header = T, sep= ";",na.strings = "NA")

# run Metacycle 
#--------------
# Here, we chose to not adjust the phase depending on the period 
# When adjusting the phase, the phase obtained for lhy for example is 2.4 while it is clearly at 0 when looking at the data
# We encourage anyone using this method or any other tool to detect oscillations to try the different options and to see what is more adapted to the dataset

TOT_Metacycle_2 <- meta2d(infile="txt",filestyle = "txt",
         minper=20, maxper=28, timepoints=rep(seq(0, 24, by=3), each=3),
         outputFile=FALSE, adjustPhase = "notAdjusted",
         combinePvalue = "fisher", inDF = TOT_expressed)

TOT_Metacycle2_results <- TOT_Metacycle_2$meta


# Metacycle TRAP mRNAs
#######################
rlog_TRAP <- read.csv("rlog_values_TRAP_expressed.csv", header = T, sep= ";",na.strings = "NA")


# run Metacycle without adjusting the phase in function of the period
#--------------------------------------------------------------------
TRAP_Metacycle <- meta2d(infile="txt",filestyle = "txt",
                          minper=20, maxper=28, timepoints=rep(seq(0, 24, by=3), each=3),
                          outputFile=FALSE, adjustPhase = "notAdjusted",
                          combinePvalue = "fisher", inDF = TRAP_expressed)

TRAP_Metacycle_results <- TRAP_Metacycle$meta



