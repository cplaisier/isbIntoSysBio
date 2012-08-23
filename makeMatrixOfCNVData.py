import os, math

# Load CNV data
print 'Reading in CNV data...'
ls1 = os.listdir('MSKCC__HG-CGH-244A/Level_3')
cnvs = {}
for file1 in ls1:
    if file1.split('.')[-1]=='txt':
        inFile = open('MSKCC__HG-CGH-244A/Level_3/'+file1, 'r')
        inFile.readline()
        for line in inFile.readlines():
            splitUp = line.strip().split('\t')
            tmp = splitUp[0].split('-')
            patient = tmp[0]+'-'+tmp[1]+'-'+tmp[2]
            if not patient in cnvs:
                cnvs[patient] = {}
            if not splitUp[1] in cnvs[patient]:
                cnvs[patient][splitUp[1]] = []
            if not splitUp[5]=='NA':
                cnvs[patient][splitUp[1]].append([int(splitUp[2]), int(splitUp[3]), float(splitUp[5])])
        inFile.close()

# Load miRNA expression data
print 'Reading in miRNA expression data...'
ls1 = os.listdir('../Expression-miRNA/UNC__H-miRNA_8x15K/Level_3')
mirnas = {}
mirnasKeys = []
for file1 in ls1:
    if file1.split('.')[-1]=='txt':
        inFile = open('../Expression-miRNA/UNC__H-miRNA_8x15K/Level_3/'+file1, 'r')
        inFile.readline()
        lines = [line.strip().split('\t') for line in inFile.readlines()]
        tmp = lines[0][0].split('-')
        patient = tmp[0]+'-'+tmp[1]+'-'+tmp[2]
        mirnaNames = [line[1] for line in lines]
        expression = [line[2] for line in lines]
        mirnas[patient] = dict(zip(mirnaNames,expression))
        mirnasKeys = list(set(mirnasKeys+mirnaNames))
        inFile.close()
     
# Load gene expression data
print 'Reading in gene expression [1] data...'
ls1 = os.listdir('../Expression-Genes/UNC__AgilentG4502A_07_1/Level_3')
genesKeys = []
genes = {}
for file1 in ls1:
    if file1.split('.')[-1]=='txt':
        inFile = open('../Expression-Genes/UNC__AgilentG4502A_07_1/Level_3/'+file1, 'r')
        inFile.readline()
        lines = [line.strip().split('\t') for line in inFile.readlines()]
        tmp = lines[0][0].split('-')
        patient = tmp[0]+'-'+tmp[1]+'-'+tmp[2]
        geneNames = [line[1] for line in lines]
        expression = [line[2] for line in lines]
        genes[patient] = dict(zip(geneNames,expression))
        genesKeys = list(set(genesKeys+geneNames))
        inFile.close()
print 'Reading in gene expression [2] data...'
ls1 = os.listdir('../Expression-Genes/UNC__AgilentG4502A_07_2/Level_3')
for file1 in ls1:
    if file1.split('.')[-1]=='txt':
        inFile = open('../Expression-Genes/UNC__AgilentG4502A_07_2/Level_3/'+file1, 'r')
        inFile.readline()
        lines = [line.strip().split('\t') for line in inFile.readlines()]
        tmp = lines[0][0].split('-')
        patient = tmp[0]+'-'+tmp[1]+'-'+tmp[2]
        geneNames = [line[1] for line in lines]
        expression = [line[2] for line in lines]
        genes[patient] = dict(zip(geneNames,expression))
        genesKeys = list(set(genesKeys+geneNames))
        inFile.close()

# Read in genomic coordinates for the miRNAs
genomic_miRNA = {'cdk4':{'chr':'12', 'start':56428270, 'end':56432431}}
inFile = open('miRNACoordinates.csv','r')
inFile.readline() # ,fold.change,t.stats,t.p,chr,start,end
for line in inFile.readlines():
    splitUp = line.strip().split(',')
    genomic_miRNA[splitUp[0]] = {'chr':splitUp[4], 'start':int(splitUp[5]), 'end':int(splitUp[6])}
inFile.close()

# Write out results for:
# hsa-miR-26a-2 = chr12:56,504,659-56,504,742
# CDK4 = chr12:56,428,270-56,432,431
outFile = open('cnvData_miRNAExp_geneExp.csv','w')
#outFile.write(',cnv.cdk4,cnv.hsa-miR-26a-2,exp.hsa-miR-26a,'+','.join(['exp.'+i for i in genesKeys])+'\n')
patients = list(set(cnvs.keys()+mirnas.keys()+genes.keys()))
outFile.write(','+','.join(patients))
writeMe = []
cdk4 = []
hsa_mir_26a = []
#genomic_miRNA = {'hsa-miR-10b':{'chr':'2', 'start':176723276, 'end':176723386}, 'cdk4':{'chr':'12', 'start':56428270, 'end':56432431}, 'hsa-miR-26a':{'chr':'12', 'start':56499977, 'end':56527014}}

for miRNA in genomic_miRNA:
    values = []
    for patient in patients:
        value = 'NA'
        if patient in cnvs and genomic_miRNA[miRNA]['chr'] in cnvs[patient]:
            for cnv in cnvs[patient][genomic_miRNA[miRNA]['chr']]:
                if ((cnv[0]<=genomic_miRNA[miRNA]['start'] and genomic_miRNA[miRNA]['end']<=cnv[1]) or genomic_miRNA[miRNA]['start']<=cnv[0]<=genomic_miRNA[miRNA]['end'] or genomic_miRNA[miRNA]['start']<=cnv[1]<=genomic_miRNA[miRNA]['end']):
                    if value=='NA' or math.fabs(cnv[2]) > math.fabs(value):
                        value = cnv[2]
        values.append(str(value))
    if not sum([1 for i in values if not i=='NA'])==0:
        outFile.write('\ncnv.'+miRNA+','+','.join(values))
for mirna in mirnasKeys: # ['hsa-miR-26a']:
    writeMe = []
    for patient in patients:
        if patient in mirnas:
            writeMe.append(mirnas[patient][mirna])
        else:
            writeMe.append('NA')
    outFile.write('\nexp.'+mirna+','+','.join(writeMe))
for gene in genesKeys: # ['CDK4','CTDSP2']:
    writeMe = []
    for patient in patients:
        if patient in genes:
            writeMe.append(genes[patient][gene])
        else:
            writeMe.append('NA')
    outFile.write('\nexp.'+gene+','+','.join(writeMe))

outFile.close()

