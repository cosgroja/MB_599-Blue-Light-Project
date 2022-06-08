import sys, re, os
import math

import seaborn as sns

import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

#################
## SUBROUTINES ##
#################

def filterdata (data):
    n = 0
    LCData = []
    for number in data:
        if number == 'N/A':
            n +=1
        else:
            LCData.append(float(number))
    return LCData, n

def sumNormalization (data):
    newdata = []
    length = len(data)
    number = len(data[0])
    first = True
    for n in range(number):
        total = 0
        for l in range(length):
            print data[l]
            total += data[l][n]
        for l in range(length):
            if first:
                newdata.append([data[l][n]/total])
            else:
                newdata[l].append(data[l][n]/total)
        first = False
    return newdata

def logTransormation (data):
    newdata = []
    length = len(data)
    number = len(data[0])
    for l in range(length):
        newdata.append([])
        for n in range(number):
            #print data[l][n]
            newdata[l].append(math.log10(data[l][n]))
    return newdata

def rangeScaling (dataD, dataB):
    newdataD = []
    newdataB = []
    length = len(dataD)
    numberD = len(dataD[0])
    numberB = len(dataB[0])
    for l in range(length):
        print dataD[l]
        newdataD.append([])
        newdataB.append([])
        mindata = min(dataD[l]+dataB[l])
        print mindata
        maxMmin = max(dataD[l]+dataB[l]) - min(dataD[l]+dataB[l])
        for n in range(numberD):
            newdataD[l].append((dataD[l][n]-mindata)/maxMmin)
        for n in range(numberB):
            newdataB[l].append((dataB[l][n]-mindata)/maxMmin)
    return newdataD, newdataB

def createDensityData (data):
    newdata = []
    length = len(data)
    number = len(data[0])
    for l in range(length):
        for n in range(number):
            newdata.append(data[l][n])
    return newdata
##########
## MAIN ##
##########

#python scripts/00_DataNormalization.py LC10.txt



usage = "Usage: " + sys.argv[0] + "<LC Original Data.txt>"
if len(sys.argv) != 2:
    print usage
    sys.exit()


###########
## INPUT ##
###########

#1. original data from LCMS
#OriginalData.txt
#Current MS Compounds	HMDB ID	KEGG ID	
#10B-1	10B-2	10B-3	10B-4	10B-5	10B-6	10B-7	10B-8	10D-1	10D-2	10D-3	10D-4	10D-5	10D-6	10D-7	10D-8
#14B-1	14B-2	14B-4	14B-5	14B-6	14B-7	14B-8	14D-1	14D-2	14D-3	14D-4	14D-5	14D-6	14D-7	14D-8

LCMS = sys.argv[1]
data = {}
HMDBID = []
dataD = []
dataB = []
day = re.search('LC(\d+)_OriginalData',LCMS).group(1)
with open(LCMS,'r') as F:
    for line in F:
        if not line.startswith('Current'):
            compoundName, HMDB, KEGG = line.strip().split('\t')[0:3]
            term = line.strip().split('\t')[3:]
            if HMDB != 'N/A':
                if day == '10':
                    B,nb = filterdata(term[0:8])
                    D,nd = filterdata(term[8:])
                    if nb<4 and nd<4:
                        data[HMDB] = {}
                        if len(B) < 8:
                            for x in range(8-len(B)):
                                B.append(min(B))
                        if len(D) < 8:
                            for x in range(8-len(D)):
                                D.append(min(D))
                        data[HMDB]['B'] = B
                        data[HMDB]['D'] = D
                        HMDBID.append(HMDB)
                        dataB.append(B)
                        dataD.append(D)
                elif day == '14':
                    B,nb = filterdata(term[0:7])
                    D,nd = filterdata(term[7:])
                    if nb<4 and nd<4:
                        data[HMDB] = {}
                        if len(B) < 7:
                            for x in range(7-len(B)):
                                B.append(min(B))
                        if len(D) < 8:
                            for x in range(8-len(D)):
                                D.append(min(D))
                        data[HMDB]['B'] = B
                        data[HMDB]['D'] = D
                        HMDBID.append(HMDB)
                        dataB.append(B)
                        dataD.append(D)

# data[HMDB][B] = [0,1,2,3,4,5,6,7]
# HMBDID = [ID,ID, .....]
# dataD = [[HMDB1 data],[],[],.....]
# dataB = [[HMDB1 data],[],[],.....]

############
## OutPUT ##
############

#print dataB
#print dataD

#normalize by sum
dataD_normalize1 = sumNormalization (dataD)
dataB_normalize1 = sumNormalization (dataB)

#print dataD_normalize1

#Log transormation
dataD_normalize2 = logTransormation (dataD_normalize1)
dataB_normalize2 = logTransormation (dataB_normalize1)

#Range scaling
dataD_normalize3, dataB_normalize3 = rangeScaling (dataD_normalize2, dataB_normalize2)

#Data for density plot
dataD_density = createDensityData (dataD)
dataD_density_AfterNormalize = createDensityData (dataD_normalize3)
dataB_density = createDensityData (dataB)
dataB_density_AfterNormalize = createDensityData (dataB_normalize3)

#Ploting Density plot
pdfBase = '00_LC' + day + 'DensityPlot'
pdfName = pdfBase + '.pdf'
with PdfPages(pdfName) as pdf:
    pageName1 = 'LC ' + day + ' Density Plot before Normalization'
    # page 1 plot density plot
    fig, ax = plt.subplots(1, 1)
    fig.suptitle(pageName1, fontsize = 10)
    sns.distplot(dataD_density, color = "red", hist = False, kde = True, kde_kws = {'linewidth': 3}, label = 'D')
    sns.distplot(dataB_density, color = "blue", hist = False, kde = True, kde_kws = {'linewidth': 3}, label = 'B')
    plt.title('Density Plot before Normalization')
    plt.xlabel('Compound relative concentration')
    plt.ylabel('Density')
    plt.legend()
    pdf.savefig()
    plt.close()
    
    pageName2 = 'LC ' + day + 'Density Plot after Normalization'
    # page 2 plot density plot
    fig, ax = plt.subplots(1, 1)
    fig.suptitle(pageName2, fontsize = 10)
    sns.distplot(dataD_density_AfterNormalize, color = "red", hist = True, kde = True, kde_kws = {'linewidth': 3}, label = 'D')
    sns.distplot(dataB_density_AfterNormalize, color = "blue", hist = True, kde = True, kde_kws = {'linewidth': 3}, label = 'B')
    plt.title('Density Plot after Normalization')
    plt.xlabel('Compound relative concentration')
    plt.ylabel('Density')
    plt.legend()
    pdf.savefig()
    plt.close()

#Ploting Concentration map
pdfBase = '00_LC' + day + 'ConcentrationMap'
pdfName = pdfBase + '.pdf'
with PdfPages(pdfName) as pdf:
    pageName1 = 'LC ' + day + 'D Concentration Map before Normalization'
    # page 1 plot density plot D
    fig, ax = plt.subplots(1, 1)
    fig.suptitle(pageName1, fontsize = 10)
    plt.boxplot(dataD,vert = False, labels = HMDBID)
    plt.title('Density Plot with D before Normalization')
    plt.yticks([])  # Disable xticks.
    plt.ylabel('HMDB ID')
    pdf.savefig()
    plt.close()

    pageName2 = 'LC ' + day + 'B Concentration Map before Normalization'
    # page 2 plot density plot B
    fig, ax = plt.subplots(1, 1)
    fig.suptitle(pageName2, fontsize = 10)
    plt.boxplot(dataB,vert = False, labels = HMDBID)
    plt.title('Density Plot with B before Normalization')
    plt.yticks([])  # Disable xticks.
    plt.ylabel('HMDB ID')
    pdf.savefig()
    plt.close()
    
    pageName3 = 'LC ' + day + 'D Concentration Map after Normalization'
    # page 3 plot density plot D
    fig, ax = plt.subplots(1, 1)
    fig.suptitle(pageName3, fontsize = 10)
    plt.boxplot(dataD_normalize3,vert = False, labels = HMDBID)
    plt.title('Density Plot with D after Normalization')
    plt.yticks([])  # Disable xticks.
    plt.ylabel('HMDB ID')
    pdf.savefig()
    plt.close()
    
    pageName4 = 'LC ' + day + 'B Concentration Map after Normalization'
    # page 4 plot density plot D
    fig, ax = plt.subplots(1, 1)
    fig.suptitle(pageName4, fontsize = 10)
    plt.boxplot(dataB_normalize3,vert = False, labels = HMDBID)
    plt.title('Density Plot with B after Normalization')
    plt.yticks([])  # Disable xticks.
    plt.ylabel('HMDB ID')
    pdf.savefig()
    plt.close()
    
outFileD ='00_LC_' + day + 'DData.txt'
with open(outFileD , 'w+') as outputFile:
    outputFile.write('HMDBID\t%s' % (day+'D'))
    for l in range(len(dataD_normalize3)):
        outputFile.write('\n')
        outputFile.write(HMDBID[l])
        for x in dataD_normalize3[l]:
            outputFile.write('\t%f' % (x))
outputFile.close()

outFileB ='00_LC_' + day + 'BData.txt'
with open(outFileB , 'w+') as outputFile:
    outputFile.write('HMDBID\t%s' % (day+'B'))
    for l in range(len(dataB_normalize3)):
        outputFile.write('\n')
        outputFile.write(HMDBID[l])
        for x in dataB_normalize3[l]:
            outputFile.write('\t%f' % (x))
outputFile.close()