setwd("C:/Users/Quinn Wong/Documents/SalkProjects")
library(ggplot2)
mydata <- read.table("Genes.txt")
genes <- unique(as.character(mydata$V1))
######
slice = 6
composite <- SpatialEnrichment(slice = 6, genes,reps = 20)
PlotBrain(composite@composite, Breaks = 9)
### Find the genes in that region
#regiongenes <- GetGenes(genes,tissueExp1,region.name="VZ")
boot <- Boot(tissueExp1 = composite@tissueExp1, tissueExp2 = composite@tissueExp2, random.matrix = composite@random.matrix)
write.table(x = boot,file = "boot.txt",quote=FALSE)

# Read a file and store into a variable
bootdata <- read.table("Mo_EGLU.boot.tsv", header=TRUE)
change <- bootdata$change

# Change parameter calculations
change_sd <- sd(change)
change_mean <- mean(change)

# Calculate the z-scores and p-values of the tissues
z <- (bootdata[, "change"] - change_mean) / change_sd
p <- pnorm(z)

bootdata$z_score <- z
bootdata$p_value <- p
bootdata[bootdata$change > change_mean, "new_p"] <- 1 - bootdata[bootdata$change > change_mean, "p_value"]
bootdata[bootdata$change > change_mean, "p_value"] <- bootdata[bootdata$change > change_mean, "new_p"]
bootdata$new_p <- NULL

write.table(bootdata, "Mo_EGLU.boot.tsv")


# Create a plot labelling all the tissues compared
bootdata$geneName <- row.names(bootdata)
p <- ggplot(bootdata, aes(reorder(geneName, change), change, fill=diff))+
  geom_bar(stat = "identity")
p + labs(x = "Tissue", y = "Change", title = "GPe_PaV")


## Create a density plot to pinpoint the tissues of interest
# Red for high difference in change, green for low
plot <- ggplot(bootdata, aes(change))+
  geom_density(fill = "grey")+
  geom_vline(size=2,xintercept = bootdata[bootdata$diff == "high","change"],color = "red")+
  geom_vline(size=2,xintercept = bootdata[bootdata$diff == "low","change"],color = "green")
plot + labs(title = "GPe_PaV")
