#################################################################
# @Course: Introduction to Systems Biology                      #
# @Rscript: inferringBiologicallyMeaningfulRelationships.R      #
# @Version: 1                                                   #
# @Author: Chris Plaisier                                       #
# @Sponsored by:                                                #
# Nitin Baliga, ISB                                             #
# Institute for Systems Biology                                 #
# 1441 North 34th Street                                        #
# Seattle, Washington  98103-8904                               #
# (216) 732-2139                                                #
# @Also Sponsored by:                                           #
# Luxembourg Systems Biology Grant                              #
# American Cancer Society Research Scholar Grant                #
#                                                               #
# If this program is used in your analysis please mention who   #
# built it. Thanks. :-)                                         #
#                                                               #
# Copyright (C) 2010 by Institute for Systems Biology,          #
# Seattle, Washington, USA.  All rights reserved.               #
#                                                               #
# This source code is distributed under the GNU Lesser          #
# General Public License, the text of which is available at:    #
#   http://www.gnu.org/copyleft/lesser.html                     #
#################################################################

### Load data for analysis from course web site ### 
# Load up image from earlier differential expresssion analysis
con = url('http://baliga.systemsbiology.net/events/sysbio/sites/baliga.systemsbiology.net.events.sysbio/files/uploads/differentialExpressionAnalysis.RData')
load(con)
close(con)

# read.csv can access web links just like any local system file
g1 = read.csv('http://baliga.systemsbiology.net/events/sysbio/sites/baliga.systemsbiology.net.events.sysbio/files/uploads/geneExp_shorter.csv', header=T, row.names=1)
dim(g1)
dimnames(g1)

### hsa-miR-26a PITA targets ###
# First load up the hsa-miR-26a PITA predicted target genes
mir26a = read.csv('http://baliga.systemsbiology.net/events/sysbio/sites/baliga.systemsbiology.net.events.sysbio/files/uploads/hsa_miR_26a_PITA.csv', header = T)

# Create data matrix for analysis
g2 = as.matrix(g1[sub('exp.','',rownames(g1)) %in% mir26a[,2],])
g3 = as.matrix(sapply(1:ncol(g2), function(i) {as.numeric(g2[,i])}))
dimnames(g3) = dimnames(g2)
g3 = rbind(g3, 'exp.hsa-miR-26a' = as.numeric(d1['exp.hsa-miR-26a',]))

# Calculate correlations between PITA predicted genes and hsa-miR-26a
c1.r = 1:(nrow(g3)-1)
c1.p = 1:(nrow(g3)-1)
for(i in 1:(nrow(g3)-1)) {
  c1 = cor.test(g3[i,], g3['exp.hsa-miR-26a',])
  c1.r[i] = c1$estimate
  c1.p[i] = c1$p.value
}

# Plot correlation coefficients
hist(c1.r, breaks = 15, main = 'Distribution of Correaltion Coefficients', xlab = 'Correlation Coefficient')

# Do testing correction
p.bonferroni = p.adjust(c1.p, method='bonferroni')
p.benjaminiHochberg = p.adjust(c1.p, method='BH')

# How many miRNA are considered significant via p-value only
print(paste('P-Value Only:  Uncorrected = ',sum(c1.p<=0.05),'; Bonferroni = ',sum(p.bonferroni<=0.05),'; Benjamini-Hochberg = ',sum(p.benjaminiHochberg<=0.05),sep=''))

# How many miRNAs are considered significant via both p-value and a negative correlation coefficient
print(paste('P-value and Rho:  Uncorrected = ',sum(c1.p<=0.05 & c1.r<=-0.15),'; Bonferroni = ',sum(p.bonferroni<=0.05 & c1.r<=-0.15),'; Benjamini-Hochberg = ',sum(p.benjaminiHochberg<=0.05 & c1.r<=-0.15 ),sep=''))

#The significantly negatively correlated miRNAs
sub('exp.','',rownames(g3)[which(p.benjaminiHochberg<=0.05 & c1.r<=-0.15)])

# Create index ordered by Benjamini-Hochberg corrected p-values to sort each vector
o1 = order(c1.r)

# Make a data.frame with the three columns
hsa_mir_26a_c1 = data.frame(rho = c1.r[o1], c.p = c1.p[o1], c.p.bonferroni = p.bonferroni[o1], c.p.benjaminiHochberg = p.benjaminiHochberg[o1])

# Add miRNA names as rownames
rownames(hsa_mir_26a_c1) = sub('exp.','',rownames(g3)[-480][o1])

# Take a look at the top results
head(hsa_mir_26a_c1)

# Plot top correlated gene
plot(as.numeric(g3['exp.ALS2CR2',]) ~ as.numeric(g3['exp.hsa-miR-26a',]), col = rgb(0, 0, 1, 0.5), pch = 20, xlab = 'miRNA Expresion', ylab = 'Gene Expression', main = 'ALS2CR2 vs. hsa-miR-26a')
# Make a trend line and plot it
lm1 = lm(as.numeric(g3['exp.ALS2CR2',]) ~ as.numeric(g3['exp.hsa-miR-26a',]))
abline(lm1, col = 'red', lty = 1, lwd = 1)

# Select out significantly correlated genes
corGenes = rownames((g3[-480,])[o1,])[which(p.benjaminiHochberg[o1]<=0.05 & c1.r[o1]<=-0.15)]

## Plot all significantly correlated genes
# Open a PDF device to output plots
pdf('genesNegativelyCorrelatedWith_hsa_miR_26a_gbm.pdf')
# Iterate through all correlated genes
for(cg1 in corGenes) {
  plot(as.numeric(g3[cg1,]) ~ as.numeric(g3['exp.hsa-miR-26a',]), col = rgb(0, 0, 1, 0.5), pch = 20, xlab = 'miRNA Expresion', ylab = 'Gene Expression', main = paste(sub('exp.','',cg1),' vs. hsa-miR-26a:\n R = ',round(hsa_mir_26a_c1[sub('exp.','',cg1),1],2),', P-Value = ',signif(hsa_mir_26a_c1[sub('exp.','',cg1),4],2),sep=''))
  # Make a trend line and plot it
  lm1 = lm(as.numeric(g3[cg1,]) ~ as.numeric(g3['exp.hsa-miR-26a',]))
  abline(lm1, col = 'red', lty = 1, lwd = 1)
}
# Close PDF device
dev.off()

# Write out results file
write.csv(hsa_mir_26a_c1, file = 'hsa_mir_26a_correlated_target_genes_PITA.csv')


### Correlation of CNV with downstream gene expression ###
# Create data matrix for analysis
g2 = as.matrix(g1[sub('exp.','',rownames(g1)) %in% mir26a[,2],])
g3 = as.matrix(sapply(1:ncol(g2), function(i) {as.numeric(g2[,i])}))
dimnames(g3) = dimnames(g2)
g3 = rbind(g3, 'cnv.hsa-miR-26a' = as.numeric(d1['cnv.hsa-miR-26a',]))

# Calculate correlations between PITA predicted genes and hsa-miR-26a
c1.r = 1:(nrow(g3)-1)
c1.p = 1:(nrow(g3)-1)
for(i in 1:(nrow(g3)-1)) {
  c1 = cor.test(g3[i,], g3['cnv.hsa-miR-26a',])
  c1.r[i] = c1$estimate
  c1.p[i] = c1$p.value
}

# Plot correlation coefficients
hist(c1.r, breaks = 15, main = 'Distribution of Correaltion Coefficients', xlab = 'Correlation Coefficient')

# Do testing correction
p.bonferroni = p.adjust(c1.p, method='bonferroni')
p.benjaminiHochberg = p.adjust(c1.p, method='BH')

# How many miRNA are considered significant via p-value only
print(paste('P-Value Only:  Uncorrected = ',sum(c1.p<=0.05),'; Bonferroni = ',sum(p.bonferroni<=0.05),'; Benjamini-Hochberg = ',sum(p.benjaminiHochberg<=0.05),sep=''))

# How many miRNAs are considered significant via both p-value and a negative correlation coefficient
print(paste('P-value and Rho:  Uncorrected = ',sum(c1.p<=0.05 & c1.r<=-0.15),'; Bonferroni = ',sum(p.bonferroni<=0.05 & c1.r<=-0.15),'; Benjamini-Hochberg = ',sum(p.benjaminiHochberg<=0.05 & c1.r<=-0.15 ),sep=''))

#The significantly negatively correlated miRNAs
sub('exp.','',rownames(g3)[which(p.benjaminiHochberg<=0.05 & c1.r<=-0.15)])

# Create index ordered by Benjamini-Hochberg corrected p-values to sort each vector
o1 = order(c1.r)

# Make a data.frame with the three columns
hsa_mir_26a_cnv1 = data.frame(rho = c1.r[o1], c.p = c1.p[o1], c.p.bonferroni = p.bonferroni[o1], c.p.benjaminiHochberg = p.benjaminiHochberg[o1])

# Add miRNA names as rownames
rownames(hsa_mir_26a_cnv1) = sub('exp.','',rownames(g3)[-480][o1])

# Take a look at the top results
head(hsa_mir_26a_cnv1)

# Plot top correlated gene
plot(as.numeric(g3['exp.ALS2CR2',]) ~ as.numeric(g3['cnv.hsa-miR-26a',]), col = rgb(0, 0, 1, 0.5), pch = 20, xlab = 'miRNA CNV', ylab = 'Gene Expression', main = 'ALS2CR2 vs. hsa-miR-26a CNV')
# Make a trend line and plot it
lm1 = lm(as.numeric(g3['exp.ALS2CR2',]) ~ as.numeric(g3['cnv.hsa-miR-26a',]))
abline(lm1, col = 'red', lty = 1, lwd = 1)

# Select out significantly correlated genes
corGenes = rownames((g3[-480,])[o1,])[which(p.benjaminiHochberg[o1]<=0.05 & c1.r[o1]<=-0.15)]

## Plot all significantly correlated genes
# Open a PDF device to output plots
pdf('genesNegativelyCorrelatedWith_CNV_hsa_miR_26a_gbm.pdf')
# Iterate through all correlated genes
for(cg1 in corGenes) {
  plot(as.numeric(g3[cg1,]) ~ as.numeric(g3['cnv.hsa-miR-26a',]), col = rgb(0, 0, 1, 0.5), pch = 20, xlab = 'miRNA CNV', ylab = 'Gene Expression', main = paste(sub('exp.','',cg1),' vs. hsa-miR-26a CNV:\n R = ',round(hsa_mir_26a_cnv1[sub('exp.','',cg1),1],2),', P-Value = ',signif(hsa_mir_26a_cnv1[sub('exp.','',cg1),4],2),sep=''))
  # Make a trend line and plot it
  lm1 = lm(as.numeric(g3[cg1,]) ~ as.numeric(g3['cnv.hsa-miR-26a',]))
  abline(lm1, col = 'red', lty = 1, lwd = 1)
}
# Close PDF device
dev.off()

# Write out results file
write.csv(hsa_mir_26a_cnv1, file = 'CNV_hsa_mir_26a_correlated_target_genes_PITA.csv')


### Causality of association ###
# Create data matrix for analysis
g2 = as.matrix(g1[sub('exp.','',rownames(g1)) %in% mir26a[,2],])
g3 = as.matrix(sapply(1:ncol(g2), function(i) {as.numeric(g2[,i])}))
dimnames(g3) = dimnames(g2)
g3 = rbind(g3, 'exp.hsa-miR-26a' = as.numeric(d1['exp.hsa-miR-26a',]), 'cnv.hsa-miR-26a' = as.numeric(d1['cnv.hsa-miR-26a',]))

## Prepare to make a 2x2 plot
pdf('causality_hsa_miR_26a_CNV_EXP_ALS2CR2_EXP.pdf', width = 11, height = 8.5)
par(mfrow=c(2,2))

## Plot (1,1) - CNV vs. miRNA expression
# Calculate correlation between miRNA expression and miRNA copy number
c1 = cor.test(g3['cnv.hsa-miR-26a',], g3['exp.hsa-miR-26a',])
# Plot correlated miRNA expression vs. copy number variation
plot(g3['exp.hsa-miR-26a',] ~ g3['cnv.hsa-miR-26a',], col = rgb(0, 0, 1, 0.5), pch = 20, xlab = 'Copy Number', ylab = 'miRNA Expression', main = paste('miRNA Expression vs. miRNA Copy Number:\n R = ',round(c1$estimate,2),', P-Value = ',signif(c1$p.value,2),sep=''))
# Make a trend line and plot it
lm1 = lm(g3['exp.hsa-miR-26a',] ~ g3['cnv.hsa-miR-26a',])
abline(lm1, col = 'red', lty = 1, lwd = 1)

## Plot (1,2) - CNV vs. ALS2CR2 Expression
# Calculate correlation between miRNA target gene ALS2CR2 expression and miRNA copy number
c1 = cor.test(g3['cnv.hsa-miR-26a',], g3['exp.ALS2CR2',])
# Plot correlated miRNA target gene ALS2CR2 expression vs. copy number variation
plot(g3['exp.ALS2CR2',] ~ g3['cnv.hsa-miR-26a',], col = rgb(0, 0, 1, 0.5), pch = 20, xlab = 'Copy Number', ylab = 'Gene Expression', main = paste('ALS2CR2 Expression vs. miRNA Copy Number:\n R = ',round(c1$estimate,2),', P-Value = ',signif(c1$p.value,2),sep=''))
# Make a trend line and plot it
lm1 = lm(g3['exp.ALS2CR2',] ~ g3['cnv.hsa-miR-26a',])
abline(lm1, col = 'red', lty = 1, lwd = 1)

## Plot (2,1) - miRNA expression vs. CNV
# Calculate correlation between miRNA target gene ALS2CR2 expression and miRNA expression
c1 = cor.test(g3['exp.hsa-miR-26a',], g3['exp.ALS2CR2',])
# Plot correlated miRNA target gene ALS2CR2 expression vs. miRNA expression
plot(g3['exp.ALS2CR2',] ~ g3['exp.hsa-miR-26a',], col = rgb(0, 0, 1, 0.5), pch = 20, xlab = 'miRNA Expression', ylab = ' Gene Expression', main = paste('ALS2CR2 Expression vs. miRNA Expression:\n R = ',round(c1$estimate,2),', P-Value = ',signif(c1$p.value,2),sep=''))
# Make a trend line and plot it
lm1 = lm(g3['exp.ALS2CR2',] ~ g3['exp.hsa-miR-26a',])
abline(lm1, col = 'red', lty = 1, lwd = 1)

## Plot (2,2) - miRNA expression vs. CNV
# Get rid of NA's
g4 = t(na.omit(t(g3)))
# Calcualte regression model for miRNA target gene ALS2CR2 expression vs. miRNA expression
r1 = lm(g4['exp.ALS2CR2',] ~ g4['exp.hsa-miR-26a',])
# Calculate correlation between residual variation after conditioning miRNA target gene ALS2CR2 expression
# on miRNA expression against miRNA copy number
c1 = cor.test(g4['cnv.hsa-miR-26a',], r1$residuals)
# Plot correlated miRNA target gene expression vs. copy number variation
plot(r1$residuals ~ g4['cnv.hsa-miR-26a',], col = rgb(0, 0, 1, 0.5), pch = 20, xlab = 'Copy Number', ylab = 'Residual', main = paste('Residual vs. miRNA Copy Number:\n R = ',round(c1$estimate,2),', P-Value = ',signif(c1$p.value,2),sep=''))
# Make a trend line and plot it
lm1 = lm(r1$residual ~ g4['exp.hsa-miR-26a',])
abline(lm1, col = 'red', lty = 1, lwd = 1)

# Close PDF device
dev.off()

### Save environment as .RData file ###
save.image('inferringBiologicallyMeaningfulRelationships.RData')

#######################
### End of Tutorial ###
#######################
