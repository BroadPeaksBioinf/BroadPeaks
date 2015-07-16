__author__ = 'dima'

import sys
import numpy
import pysam
import logging
import os.path
import time
import scipy
import scipy.stats



def inputAndIndex(bamPath):
    # Check if there is an BAM file at entered dir
    if os.path.isfile(bamPath) == False:
        logging.error("BAM file does not exist")

    # Check if there is an index file, create one if there isn't
    if os.path.isfile(bamPath+".bai") == False:
        pysam.index(bamPath)
        logging.info('No index was found, new index was generated')

    # Take chromosome data from BAM index (ref.seq. name, ref.seq. length, number of mapped reads and number of unmapped reads)
    chromosomesInfo = []
    for chr in pysam.idxstats(bamPath):
        chromosomesInfo.append(chr.split("\t"))
    # Last line is unmapped reads, we don't need them
    chromosomesInfo.pop()
    #print(chromosomesInfo)

    bamfile = pysam.AlignmentFile(bamPath, "rb")

    return (bamfile, chromosomesInfo)

def countUniqueReads(bamfile, chromosomesInfo):

    totalUniqueReadsCount = 0
    i = 0
    for chromosome in chromosomesInfo:
        beginningOfThePreviousRead = 0
        currentChromosomeName = chromosome[0]
        currentChromosomeSize = int(chromosome[1])
        allReadsInChromosome = bamfile.fetch(currentChromosomeName)

        i = 0

        for read in allReadsInChromosome:
                readStr = str(read)
                beginningOfTheRead = ([int(s) for s in readStr.split() if s.isdigit()][2])
                if beginningOfTheRead != beginningOfThePreviousRead:
                    beginningOfThePreviousRead = beginningOfTheRead
                    totalUniqueReadsCount = totalUniqueReadsCount +1




    #print("Unique reads counted")

    return(totalUniqueReadsCount)

# Makes a simple list of windows, where each window is a list [WindowStart, ReadsInWindow].
# Chromosomes are separated by [-1,-1] window
# Sequences of ineligible windows loger than gap+1 are not stored
def makeWindowsList(bamfile, chromosomesInfo,l0, windowSize, gap):
    i = 0
    windowList = []
    for chromosome in chromosomesInfo:
        beginningOfThePreviousRead = 0
        currentChromosomeName = chromosome[0]
        currentChromosomeSize = int(chromosome[1])
        #print([currentChromosomeName, currentChromosomeSize, len(windowList)])
        allReadsInChromosome = bamfile.fetch(currentChromosomeName)

        gapCount = 0
        i = 0
        windowReadsCount = 0
        for read in allReadsInChromosome:
                readStr = str(read)
                beginningOfTheRead = ([int(s) for s in readStr.split() if s.isdigit()][2])
                if beginningOfTheRead != beginningOfThePreviousRead:
                    beginningOfThePreviousRead = beginningOfTheRead

                    while True:
                        if (i <= beginningOfTheRead) and (beginningOfTheRead < i + windowSize):
                            windowReadsCount = windowReadsCount + 1
                            break
                        elif ( beginningOfTheRead < i):
                            break
                        else:
                            if windowReadsCount< l0:
                                gapCount = gapCount + 1
                            windowList.append([i, windowReadsCount])
                            # If we have a g+1 sized gap, go and delete last g windows
                            if gapCount > gap:
                                while gapCount>0:
                                    windowList.pop()
                                    gapCount = gapCount-1
                            i = i + windowSize
                            windowReadsCount = 0

        # Next chromosome marker just in case
        windowList.append([-1,-1])
        print([currentChromosomeName, currentChromosomeSize, len(windowList)], "READY")
    windowList.append([1,1])
    return(windowList)

def makeIslandsList(windowList, lambdaa, windowSize,l0, chromosomesInfo,islandScoreThreshold):

    chromosomeCounter = 0
    currentChromosomeName = chromosomesInfo[chromosomeCounter][0]
    islandsList = []
    islandScore = 0
    windowStart = windowList[0][0] - windowSize
    islandStart = windowList[0][0]

    for i, window in enumerate(windowList):
        windowStartNew = window[0]

        #New chromosome check
        if window[0] == -1:
            print (currentChromosomeName + " done")

            windowStart = windowList[i+1][0] - windowSize

            chromosomeCounter = chromosomeCounter +1
            if chromosomeCounter<len(chromosomesInfo):
                currentChromosomeName = chromosomesInfo[chromosomeCounter][0]
            #print ("start " + currentChromosomeName)
        else:

            if windowStartNew != windowStart + windowSize:
                # A bug here: loads of 0-score islands are generated
                if islandScore>=islandScoreThreshold:
                    islandsList.append([currentChromosomeName,islandStart, windowStart+ windowSize, int(islandScore)])
                islandScore = 0
                islandStart = window[0]
                windowStart = windowStartNew
            else:
                # Check eligibility
                if window[1]>=l0:
                    #sometimes 0 and therefore inf in -log  is generated
                    temp = scipy.stats.poisson.pmf(window[1],lambdaa)
                    if temp == 0:
                        windowScore = 10
                    else:
                        windowScore = -numpy.log(temp)
                else:
                    windowScore = 0
                islandScore = islandScore + windowScore

                windowStart = windowStartNew

    return(islandsList)


startTime = time.time()

# Input arguments: path to BAM file
#bamPath = sys.argv[1]
#windowSize = sys.argv[2]
#gap = sys.argv[3]
#gap = int(float(gap)/float(windowSize))
#bamPath = "/home/dima/BAMfiles/Bernstein_H1_hESC_CTCF.bam"
bamPath = "/home/dima/BAMfiles/h3k4me3_rep1.bam"
windowSize = 200
p0 = 0.05
gap = 1
islandScoreThreshold = 100

# Log file
logging.basicConfig(filename=(os.path.dirname(bamPath) + '/SICER_log.log'),level=logging.DEBUG)

bamfile, chromosomesInfo = inputAndIndex(bamPath)
print("Counting unique reads")
totalUniqueReadsCount = countUniqueReads(bamfile, chromosomesInfo)

# Effective genome length
L = 0.77 * sum(int(row[1]) for row in chromosomesInfo)
# Lambda for poisson dist
lambdaa = float(windowSize) * float(totalUniqueReadsCount)/ float(L)
# Minimum #reads in a window for eligibility
# Formula (1), finding l0
l0 = scipy.stats.poisson.ppf(1-p0,lambdaa)
#print(L, totalUniqueReadsCount, lambdaa, l0)


print("Finished counting reads, now making window list")
windowList = makeWindowsList(bamfile, chromosomesInfo, l0, windowSize, gap)
print("Finished window list, now making island list")
islandList = makeIslandsList(windowList, lambdaa, windowSize,l0, chromosomesInfo,islandScoreThreshold)

#print(len(windowList), sys.getsizeof(windowList)/1024)
#print(len(islandList), sys.getsizeof(islandList)/1024)

f = open(bamPath[:-4] +'_peaks.bed', 'wb')
for island in islandList:
    islandString = str(island[0]) + "\t" + str(island[1]) + "\t" + str(island[2]) + "\n"
    f.write(islandString)
#print(islandList[11104])

bamfile.close()

print("Finished. Elapsed time, minutes: " + str((time.time() - startTime)/60))