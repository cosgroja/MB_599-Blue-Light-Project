import sys, re, os
import matplotlib.pyplot as plt
import numpy as np
import math
from matplotlib.colors import LinearSegmentedColormap
from scipy import stats

###############
# SUBROUTINES #
###############    

def readCompNameFile (compNameFile):
    compNameStore = {}
    with open(compNameFile, 'r') as F:
        for line in F:
            KEGG, HMDB, normalName = line.strip().split("\t")
            compNameStore[normalName] = HMDB
    return compNameStore
       
def readLCGCDataFile (LCGCDataFile):
    LCGCDataStore = {}
    with open(LCGCDataFile, 'r') as F:
        for line in F:
            if not line.startswith ('HMDBID') and not line.startswith('B/D'):
                HMDBID = line.strip().split("\t")[0]
                print HMDBID
                BData = line.strip().split("\t")[1:9]
                DData = line.strip().split("\t")[9:]
                LCGCDataStore[HMDBID] = {}
                LCGCDataStore[HMDBID]['treatment'] = [float(x) for x in BData]
                LCGCDataStore[HMDBID]['control'] = [float(x) for x in DData]
    return LCGCDataStore

def plotFigure(figureName, pdfName, xValue, yValue, cols, sizes, cm, label):
    sizes, xValue, yValue, cols, label = zip(*sorted(zip(sizes, xValue, yValue, cols, label)))
    sizes = list(sizes)
    sizes.reverse()
    xValue = list(xValue)
    xValue.reverse()
    yValue = list(yValue)
    yValue.reverse()
    cols = list(cols)
    cols.reverse()
    label = list(label)
    label.reverse()
    plt.scatter(xValue,yValue,c=cols,s=sizes,cmap=cm,edgecolor='none',alpha=0.8)
    plt.xlabel("Rate of Hits/Total Metabolite %")
    plt.ylabel("-log 10 (Geometric Mean)")
    plt.colorbar(label='-log10(FDR)')
    for a in range(len(label)):
        if cols[a] >= 0.0:
            plt.annotate(label[a], (xValue[a], yValue[a]), fontsize=10)
    plt.title(figureName)
    plt.savefig(pdfName)
    plt.close()
########
# Main #
########

#python scripts/06_PathwayVisual_10d.py 06_d10_pathway_results.txt 06_d10_name_map.txt 06_d10Data.txt

usage = "Usage: " + sys.argv[0] + " <pathway result. txt> <normal company Name.txt> <LCGC data .txt>"
if len(sys.argv) != 4:
    print usage
    sys.exit()

# create same gene list, load diff files, and store the reads
pathwayResultFile = sys.argv[1]
compNameFile = sys.argv[2]
compNameStore = readCompNameFile (compNameFile)
LCGCDataFile = sys.argv[3]
LCGCDataStore = readLCGCDataFile (LCGCDataFile)
print LCGCDataStore

colors = ['m','black']
cm = LinearSegmentedColormap.from_list('dah', colors, N=100)

xValue = []
yValue = []
sizes = []
cols = []
labels = []

with open(pathwayResultFile , 'r') as pathwayFile:
    for line in pathwayFile:
        if not line.startswith ('\tTotal'):
            pathwayName, total, expected, hits, pValue, log10P, holm, FDR, impact = line.strip().split('\t')[0:9]
            terms = line.strip().split('\t')[9:]
            
            sizes.append(int(hits)*50)
            per = (float(hits)/float(total))*100
            xValue.append(float(per))
            P_Store = []
            for term in terms:
                HMDBID = compNameStore[term]
                DData = LCGCDataStore[HMDBID]['control']
                BData = LCGCDataStore[HMDBID]['treatment']
                #print (stats.ttest_ind(DData, BData, equal_var=False))
                p_Value = stats.ttest_ind(DData, BData, equal_var=False)[1]
                #print p_Value
                P_Store.append(p_Value)
            #product
            multSum = 1
            for p in P_Store:
                multSum = multSum * p
            #print multSum
            geometricMean = multSum**(1.0/(len(P_Store)))
            #print len(P_Store)
            #print geometricMean
            y_Value = -math.log10(geometricMean)
            yValue.append(float(y_Value))
            log10FDR = -math.log10(float(FDR))
            if float(log10FDR) >= 3.0:
                cols.append(float(3.0))
            else:
                cols.append(float(log10FDR))
            if ' metabolism' in pathwayName:
                pathwayName = pathwayName.replace(" metabolism", "")
            labels.append(pathwayName)
pathwayFile.close()

figureName = '10 days Blue Light Metabolomic Pathway Analysis'

pdfName = '06_10dMetabolomicPathwayResult.pdf'
#pdfName = '06_10dMetabolomicPathwayResult.pdf'

plotFigure(figureName, pdfName, xValue, yValue, cols, sizes, cm, labels)