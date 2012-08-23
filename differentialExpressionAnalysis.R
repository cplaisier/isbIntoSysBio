#################################################################
# @Course: Introduction to Systems Biology                      #
# @Rscript: differentialExpressionAnalysis.R                    #
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
# read.csv can access web links just like any local system file
d1 = read.csv('http://baliga.systemsbiology.net/events/sysbio/sites/baliga.systemsbiology.net.events.sysbio/files/uploads/cnvData_miRNAExp_0.csv', header=T, row.names=1)
dim(d1)
dimnames(d1)

### Divy up data into consituent vectors / matrices ###
# Case or control status (1 = case, 0 = control)
case_control = d1[1,]
# miRNA expression levels
mirna = d1[325:nrow(d1),]

### Take a look at the data ###
mu1 = 1:nrow(mirna)
sd1 = 1:nrow(mirna)
for(i in 1:nrow(mirna)) {
  mu1[i] = mean(as.numeric(mirna[i,]), na.rm = T)
  sd1[i] = sd(as.numeric(mirna[i,]), na.rm = T)
}
plot(mu1, sd1, col = rgb(0,0,1,0.5), pch = 20, ylab = 'Standard Deviation', xlab = 'Mean')

# Or

plot(apply(mirna,1,mean,na.rm=T), apply(mirna,1,sd,na.rm=T), col = rgb(0,0,1,0.5), pch = 20, ylab = 'Standard Deviation', xlab = 'Mean')


### Calculate Fold-Change ###
# http://bioinfostore.unc.edu/tcga/miRNA_Normailzation_Level_2_3_description.doc -> They use log2 normalization which would affect the fold changes
fc = rep(NA,nrow(mirna))
for(i in 1:nrow(mirna)) {
  fc[i] = median(2^as.numeric(mirna[i,which(case_control==1)]),na.rm=T)/median(2^as.numeric(mirna[i,which(case_control==0)]),na.rm=T)
}

# Or

fc.2 = sapply(1:nrow(mirna), function(i) { return(median(2^as.numeric(mirna[i,which(case_control==1)]),na.rm=T)/median(2^as.numeric(mirna[i,which(case_control==0)]),na.rm=T)) } )


### Calculate Student's T-test ###
t1.t = rep(NA, nrow(mirna))
t1.p = rep(NA, nrow(mirna))
for(i in 1:nrow(mirna)) {
  t1 = t.test(mirna[i,which(case_control==1)],mirna[i,which(case_control==0)])
  t1.t[i] = t1$statistic
  t1.p[i] = t1$p.value
}

### Adjusting p-values for multiple testing (multiple testing correction) ###
# Do Bonferroni multiple testing correction (FWER)
p.bonferroni = p.adjust(t1.p, method='bonferroni')

# Do Benjamini-Hochberg multiple testing correction (FDR)
p.benjaminiHochberg = p.adjust(t1.p, method='BH')

# How many miRNA are considered significant via p-value only
print(paste('P-Value Only:  Uncorrected = ',sum(t1.p<=0.05),'; Bonferroni = ',sum(p.bonferroni<=0.05),'; Benjamini-Hochberg = ',sum(p.benjaminiHochberg<=0.05),sep=''))

# How many miRNAs are considered significant via both p-value and fold-change
print(paste('P-value and Fold-Change:  Uncorrected = ',sum(t1.p<=0.05 & (fc<=0.5 | fc>=2)),'; Bonferroni = ',sum(p.bonferroni<=0.05 & (fc<=0.5 | fc>=2)),'; Benjamini-Hochberg = ',sum(p.benjaminiHochberg<=0.05 & (fc<=0.5 | fc>=2)),sep=''))

#The significant miRNAs
sub('exp.','',rownames(mirna)[which(p.benjaminiHochberg<=0.05 & (fc<=0.5 | fc>=2))])


### Combine into one data.frame and save it out to a file ###
# Create index ordered by Benjamini-Hochberg corrected p-values to sort each vector
o1 = order(p.benjaminiHochberg)

# Make a data.frame with the three columns
td1 = data.frame(fold.change = fc[o1], t.stats = t1.t[o1], t.p = t1.p[o1], t.p.bonferroni = p.bonferroni[o1], t.p.benjaminiHochberg = p.benjaminiHochberg[o1])

# Add miRNA names as rownames
rownames(td1) = sub('exp.','',rownames(mirna)[o1])

# Write out results file
write.csv(td1,file='t.test.fc.mirna.csv')


### Volcano Plot ###
# Open a PDF output device to store the volcano plot
pdf('volcanoPlotTCGAmiRNAs.pdf')

# Make a volcano plot
plot(log(fc,2), -log(t1.p, 10) , ylab = '-log10(p-value)', xlab = 'log2(Fold Change)', axes = F, col = rgb(0, 0, 1, 0.25), pch = 20, main = "TCGA miRNA Differential Expresion", xlim = c(-6, 5.5), ylim=c(0, 120))

# Get some plotting information for later
p1 = par()

# Add the axes
axis(1)
axis(2)

## Label significant miRNAs on the plot
# Don't make a new plot just write over top of the current plot
par(new=T)

# Choose the significant miRNAs
included = c(intersect(rownames(td1)[which(td1[, 't.p']<=(0.05/534))], rownames(td1)[which(td1[, 'fold.change']>=2)]), intersect(rownames(td1)[which(td1[, 't.p']<=(0.05/534))], rownames(td1)[which(td1[, 'fold.change']<=0.5)]))

# Plot the red highlighting circles
plot(log(td1[included, 'fold.change'], 2), -log(td1[included, 't.p'], 10), ylab = '-log10(p-value)', xlab = 'log2(Fold Change)', axes = F, col = rgb(1, 0, 0, 1), pch = 1, main = "TCGA miRNA Differential Expresion", cex = 1.5, xlim = c(-6, 5.5), ylim = c(0, 120))

# Add labels as text
text((log(td1[included, 'fold.change'], 2)), ((-log(td1[included, 't.p'], 10))+-3), included, cex = 0.4)

# Close PDF output device, closes PDF file so it can be opened
dev.off()


### Integrating miRNA expression with CNVs ###
# Subset the original data matrix to get create a matrix of the Copy Number Variation data for each miRNA
# This was calcualted by tabulating the CNV levels for each miRNA across the genome
cnv = d1[2:324,]

# Run the analysis for hsa-miR-10b
c1 = cor.test(as.numeric(cnv['cnv.hsa-miR-10b',]), as.numeric(mirna['exp.hsa-miR-10b',]), na.rm = T)

# Plot hsa-miR-10b expression vs. Copy Number levels
plot(as.numeric(mirna['exp.hsa-miR-10b',]) ~ as.numeric(cnv['cnv.hsa-miR-10b',]), col = rgb(0,0,1,0.5), pch = 20, xlab = 'Copy Number', ylab = 'hsa-miR-10b Expression', main = 'hsa-miR-10b:\n Expression vs. Copy Number')

# Add a trend line to the plot
lm1 = lm(as.numeric(mirna['exp.hsa-miR-10b',]) ~ as.numeric(cnv['cnv.hsa-miR-10b',]))
abline(lm1,col='red',lty=1,lwd=1)

# Create a matrix to strore the output
m1 = matrix(nrow=nrow(cnv),ncol=2,dimnames=list(rownames(cnv),c('cor.r','cor.p')))

# Run the analysis
for(i in rownames(cnv)) {
  # Try function catches errors and doesn't stop execution, issue because of missing data
  c1 = try(cor.test(as.numeric(cnv[i,]), as.numeric(mirna[sub('cnv','exp',i),]), na.rm = T), silent = T)
  # If there is no error in analysis then add values to matrix, otherwise keep them as NA
  m1[i, 'cor.r'] = ifelse(class(c1)=='try-error', 'NA', c1$estimate)
  m1[i, 'cor.p'] = ifelse(class(c1)=='try-error', 'NA', c1$p.value)
}

# Adjust p-values and get rid of NA's using na.omit
m2 = na.omit(cbind(m1,p.adjust(as.numeric(m1[,2]), method = 'bonferroni'), p.adjust(as.numeric(m1[,2]), method = 'BH')))

# Create index ordered by correlation coefficient to sort the entire matrix
o1 = order(m2[,1], decreasing = T)

# Write out results file
write.csv(m2[o1,], file = 'corTestCnvExp_miRNA_gbm.csv')

# Select a subset to plot
top15 = sub('cnv.', '', rownames(head(m2[order(as.numeric(m2[,1]), decreasing = T), ], n = 15)))

## Plot top 15 correlations
# Open a PDF device to output plots
pdf('corTestCnvExp_miRNA_gbm.pdf')
# Iterate through all the top 15 miRNAs
for(mi1 in top15) {
  # Plot correlated miRNA expression vs. copy number variation
  plot(as.numeric(mirna[paste('exp.', mi1, sep = ''),]) ~ as.numeric(cnv[paste('cnv.', mi1, sep = ''),]), col = rgb(0, 0, 1, 0.5), pch = 20, xlab = 'Copy Number', ylab = 'Expression', main = paste(mi1, '\n Expression vs. Copy Number'), sep = '')
  # Make a trend line and plot it
  lm1 = lm(as.numeric(mirna[paste('exp.', mi1, sep = ''),]) ~ as.numeric(cnv[paste('cnv.', mi1, sep = ''),]))
  abline(lm1, col = 'red', lty = 1, lwd = 1)
}
# Close PDF device
dev.off()


### Save environment as .RData file ###
save.image('differentialExpressionAnalysis.RData')

#######################
### End of Tutorial ###
#######################
