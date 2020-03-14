#!/usr/bin/env python
"""
Created on Mon Aug  5 13:25:02 2019

@author: huyn
"""
import json
import argparse
import matplotlib.pyplot as plt # for plot
from collections import Counter
import numpy as np
from file_handle import traverseAll,parsing
import os
def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("--naiveInput","-n",help="naive directory (either time or result)",default ="result_naive/E_Coli/" )
    parser.add_argument("--approxInput","-a",help="aprrox directory (either time or result)",default = "result_approx/E_Coli/")
    parser.add_argument("--ILPInput","-i",help="ILP directory (either time or result)",default ="result_approx_ILP/E_Coli/" )
    parser.add_argument("--LPInput","-l",help="LP directory (either time or result)",default = "result_approx_LP/E_Coli/")
    parser.add_argument("--geneBlock","-b",help="gene_block_names_and_genes.txt file",default = "./E_Coli/gene_block_names_and_genes.txt")
#    parser.add_argument("--isTime","-t",help="making graph for time (Y or N)",default = "N") 
#    parser.add_argument("--isReconstruction","-r",help="making graph for reconstruction (Y or N)",default = "N")  
    parser.add_argument("-o", "--output", default='analysis',
                help="graph file")                                
    return parser.parse_args()
    
###############################################################################
## given directory, compute a dictionary to store 
###############################################################################
def computeDeletion(string1,string2):
    s1 = set(string1)
    s2 = set(string2)
    try:
        s1.remove("|")
    except:
        pass
    try:
        s2.remove("|")
    except:
        pass    
    return len(s1.symmetric_difference(s2))
    
    

def computeSplit(string1,string2):
    string1 = string1.split("|")
    string2 = string2.split("|")
    return abs(len(string1)-len(string2))
    
    
def computeDuplication(string1,string2):
    string1,string2 = string1.split("|"),string2.split("|")
    duplicated1,duplicated2= set(),set()
    for block in string1:
        c = Counter(block)
        for item in c:
            if c[item]>=2:
                duplicated1.add(item)
    for block in string2:
        c = Counter(block)
        for item in c:
            if c[item]>=2:
                duplicated2.add(item)  
#    print (duplicated1,duplicated2)
    return len(duplicated1.symmetric_difference(duplicated2))
    
def computePairWiseDistance(dictionary,geneSet):
    pairWiseDistance = [0,0,0]
    distance = [0,0,0]
    myList   = list(dictionary.keys())
#    print (dictionary)
    for i in range(len(myList)-1):
        currentGeneBlock = dictionary[myList[i]]
        numBlock = len(currentGeneBlock.split("|"))    
        genes = set(currentGeneBlock) 
        try:
            genes.remove("|")
        except:
            pass
        deletion = len(geneSet)-len(genes)
        split    = abs(numBlock-1)
        distance[0] += deletion
        distance[2] += split
        for j in range(i+1,len(myList)):
            nextGeneBlock =dictionary[myList[j]]
            deletion     = computeDeletion(currentGeneBlock,nextGeneBlock)
            dup          = computeDuplication(currentGeneBlock,nextGeneBlock)
            split        = computeSplit(currentGeneBlock,nextGeneBlock)
            pairWiseDistance[0]+=deletion
            pairWiseDistance[1]+= dup
            pairWiseDistance[2]+=split
    return pairWiseDistance,distance   
    
def computeDistanceVSILP(dictionary1,dictionary2):
    pairWiseDistance = [0,0,0]
    distance = [0,0,0]
    myList   = list(dictionary.keys())
#    print (dictionary)
    for i in range(len(myList)-1):
        currentGeneBlock1 = dictionary1[myList[i]]
        currentGeneBlock2 = dictionary2[myList[i]]
        deletion     = computeDeletion(currentGeneBlock1,currentGeneBlock2)
        split        = computeSplit(currentGeneBlock1,currentGeneBlock2)
        distance[0] += deletion
        distance[2] += split

    return pairWiseDistance,distance      
###############################################################################
## parse data from the 3 input file
###############################################################################
def parseDirectory(directory,geneData):
    timeFile = directory + "time.txt"
    new_result = directory + "new_result"
    distanceDictionary = {"pairWiseDistance":{},"referenceDistance":{}}
    for file in traverseAll(new_result):
        operonName = file.split("/")[-1]
        mapping,genomes = parsing(file)
        pairWiseDistance,distance   = computePairWiseDistance(genomes,set(geneData[operonName]))
        distanceDictionary["pairWiseDistance"][operonName] = pairWiseDistance
        distanceDictionary["referenceDistance"][operonName] = distance
    with open(timeFile,"r") as infile:
        timeDictionary = json.load(infile)
    return distanceDictionary,timeDictionary
def parseData(naiveInput,approxInput,geneBlock):
    dictionary = {"time":{},"pairWiseDistance":{},"referenceDistance":{}}
    geneData = {}
    with open(geneBlock,"r") as textFile:
        for line in textFile.readlines():
            line = line.strip().split()
            geneData[line[0]]=line[1:]
    naiveDistanceDictionary,naiveTimeDictionary = parseDirectory(naiveInput,geneData)
    approxDistanceDictionary,approxTimeDictionary = parseDirectory(approxInput,geneData)
    for operon in geneData:
        dictionary["pairWiseDistance"][operon] = {"naive":naiveDistanceDictionary["pairWiseDistance"][operon],
                              "approx":approxDistanceDictionary["pairWiseDistance"][operon],
                              "genes":geneData[operon]      }
        dictionary["referenceDistance"][operon] = {"naive":naiveDistanceDictionary["referenceDistance"][operon],
                              "approx":approxDistanceDictionary["referenceDistance"][operon],
                              "genes":geneData[operon]      }
        dictionary["time"][operon] = {"naive":naiveTimeDictionary[operon],
                              "approx":approxTimeDictionary[operon],
                              "genes":geneData[operon]      }
    return dictionary
    
###############################################################################
## Explore data
###############################################################################
# given x,y data, draw 1 figure and save it
def drawOne(x,operonName,dictionary,name,index,indices,isTime,directory,labelTop,labelBot):
    plt.figure()
    plt.title(name.split()[0])
    if index!=None: 
#        print (dictionary)
        naiveData  = [dictionary[operon]["naive"][index] for operon in operonName]
        approxData = [dictionary[operon]["approx"][index] for operon in operonName]
    else:
        naiveData  = np.log10([dictionary[operon]["naive"] for operon in operonName])
        approxData = np.log10([dictionary[operon]["approx"] for operon in operonName])    
    differences =[]
    for i in range(len(naiveData)):
        differences.append(naiveData[i]-approxData[i])
    # change xticks
    geneIndices = sorted(indices.keys())
    flipD =  {}
    for gene in geneIndices:
        flipD[indices[gene]] = gene
    plt.xticks(x,["{} ({})".format(operonName[i][:3],flipD[i]) if i in flipD else operonName[i][:3]  for i in range(len(operonName))],rotation='vertical',fontsize = 10)
    plt.axvline(x=indices[geneIndices[0]],color='k',label = "Number of genes in the operon")
    for i in geneIndices[1:]:
        plt.axvline(x=indices[i],color='k')
    plt.scatter(x, naiveData, s= 20,c="r",label=labelTop)
    plt.scatter(x, approxData, s=20,c="b",label=labelBot)
        
#    plt.plot(x, differences, 'g-o',label="Differences")
#    print (differences)
#        plt.plot(operonName, geneData, 'g-o',label="Rep3")
    plt.xlabel('Operon Name')
    plt.ylabel(name)
    # Place a legend to the right of the plot
    plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left', borderaxespad=0.)
    # resize the plot
    
    plt.xticks(range(len(dictionary))) # add loads of ticks
    plt.grid()
    
    plt.gca().margins(x=0)
    plt.gcf().canvas.draw()
    tl = plt.gca().get_xticklabels()
    maxsize = max([t.get_window_extent().width for t in tl])
    m = 0.2 # inch margin
    s = maxsize/plt.gcf().dpi*len(dictionary)+2*(m+1)
    margin = m/plt.gcf().get_size_inches()[0]
    
    plt.gcf().subplots_adjust(left=margin, right=1.-margin)
    plt.gcf().set_size_inches(s, plt.gcf().get_size_inches()[1])        
    
    # Pad margins so that markers don't get clipped by the axes
    plt.margins(0.1)
    # Tweak spacing to prevent clipping of tick-labels
    plt.subplots_adjust(bottom=0.1)
    # save somewhere we want lol
    if isTime=="Y":
        plt.savefig(directory+"/{}Plot".format(name),bbox_inches='tight',quality=100,dpi=500)
    else:
        plt.savefig(directory+"/{}Plot".format(name),bbox_inches='tight',quality=100,dpi=500)
        
        
# given x,y data, draw 1 figure and save it
def drawDifference(x,operonName,dictionary,name,index,indices,isTime,directory,label):
    plt.figure()
    plt.title(name.split()[0]+"Differences")
    if index!=None: 
#        print (dictionary)
        naiveData  = [dictionary[operon]["naive"][index] for operon in operonName]
        approxData = [dictionary[operon]["approx"][index] for operon in operonName]
    else:
        naiveData  = np.log10([dictionary[operon]["naive"] for operon in operonName])
        approxData = np.log10([dictionary[operon]["approx"] for operon in operonName])    
    differences =[]
    info = []
    for i in range(len(naiveData)):
        differences.append(naiveData[i]-approxData[i])
        info.append([operonName[i],naiveData[i]-approxData[i]])
#    print (info)
    print ()
    # change xticks
    geneIndices = sorted(indices.keys())
    flipD =  {}
    for gene in geneIndices:
        flipD[indices[gene]] = gene
    plt.xticks(x,["{} ({})".format(operonName[i][:3],flipD[i]) if i in flipD else operonName[i][:3]  for i in range(len(operonName))],rotation='vertical',fontsize = 10)
    plt.axvline(x=indices[geneIndices[0]],color='k',label = "Number of genes in the operon")
    for i in geneIndices[1:]:
        plt.axvline(x=indices[i],color='k')
    plt.scatter(x, differences, s= 20,c="black",label=label)
        
#    plt.plot(x, differences, 'g-o',label="Differences")
#    print (differences)
#        plt.plot(operonName, geneData, 'g-o',label="Rep3")
    plt.xlabel('Operon Name')
    plt.ylabel(name)
    # Place a legend to the right of the plot
    plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left', borderaxespad=0.)
    # resize the plot
    
    plt.xticks(range(len(dictionary))) # add loads of ticks
    plt.grid()
    
    plt.gca().margins(x=0)
    plt.gcf().canvas.draw()
    tl = plt.gca().get_xticklabels()
    maxsize = max([t.get_window_extent().width for t in tl])
    m = 0.2 # inch margin
    s = maxsize/plt.gcf().dpi*len(dictionary)+2*(m+1)
    margin = m/plt.gcf().get_size_inches()[0]
    
    plt.gcf().subplots_adjust(left=margin, right=1.-margin)
    plt.gcf().set_size_inches(s, plt.gcf().get_size_inches()[1])        
    
    # Pad margins so that markers don't get clipped by the axes
    plt.margins(0.1)
    # Tweak spacing to prevent clipping of tick-labels
    plt.subplots_adjust(bottom=0.1)
    # save somewhere we want lol
    if isTime=="Y":
        plt.savefig(directory+"/{}DifferencesPlot".format(name),bbox_inches='tight',quality=100,dpi=500)
    else:
        plt.savefig(directory+"/{}DifferencesPlot".format(name),bbox_inches='tight',quality=100,dpi=500)
    return differences
    
#
def drawDifferenceSum(x,operonName,dictionary,name,index,indices,isTime,directory,differences,label):
    plt.figure()
    plt.title(name.split()[0]+"SumDifferences")
    # change xticks
    geneIndices = sorted(indices.keys())
    flipD =  {}
    for gene in geneIndices:
        flipD[indices[gene]] = gene
    plt.xticks(x,["{} ({})".format(operonName[i][:3],flipD[i]) if i in flipD else operonName[i][:3]  for i in range(len(operonName))],rotation='vertical',fontsize = 10)
    plt.axvline(x=indices[geneIndices[0]],color='k',label = "Number of genes in the operon")
    for i in geneIndices[1:]:
        plt.axvline(x=indices[i],color='k')
    plt.scatter(x, differences, s= 20,c="black",label=label)
        
#    plt.plot(x, differences, 'g-o',label="Differences")
#    print (differences)
#        plt.plot(operonName, geneData, 'g-o',label="Rep3")
    plt.xlabel('Operon Name')
    plt.ylabel(name)
    # Place a legend to the right of the plot
    plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left', borderaxespad=0.)
    # resize the plot
    
    plt.xticks(range(len(dictionary))) # add loads of ticks
    plt.grid()
    
    plt.gca().margins(x=0)
    plt.gcf().canvas.draw()
    tl = plt.gca().get_xticklabels()
    maxsize = max([t.get_window_extent().width for t in tl])
    m = 0.2 # inch margin
    s = maxsize/plt.gcf().dpi*len(dictionary)+2*(m+1)
    margin = m/plt.gcf().get_size_inches()[0]
    
    plt.gcf().subplots_adjust(left=margin, right=1.-margin)
    plt.gcf().set_size_inches(s, plt.gcf().get_size_inches()[1])        
    
    # Pad margins so that markers don't get clipped by the axes
    plt.margins(0.1)
    # Tweak spacing to prevent clipping of tick-labels
    plt.subplots_adjust(bottom=0.1)
    # save somewhere we want lol
    if isTime=="Y":
        plt.savefig(directory+"/{}SumDifferencesPlot".format(name),bbox_inches='tight',quality=100,dpi=500)
    else:
        plt.savefig(directory+"/{}SumDifferencesPlot".format(name),bbox_inches='tight',quality=100,dpi=500)

    
# given dictionary, draw out figures 
def drawAll(dictionary,isTime,directory,field,label,labelTop,labelBot):
    x = [i for i in range(len(dictionary))]
    operonName = sorted([key for key in dictionary],key = lambda operon: len(dictionary[operon]["genes"]))
    indices    = {}
    for i in range(len(operonName)):
        operon = operonName[i]
        indices[len(dictionary[operon]["genes"])]=i
    output = [[0,0,0] for operon in operonName]
    if not isTime:
        graphName  = ["DeletionEvents","DuplicationEvents","SplitEvents"]
        differences = [0]*len(operonName)
        # plot with color, and legend
        # plot the deletion cost
        for index in range(3):
            drawOne(x,operonName,dictionary,field+graphName[index],index,indices,isTime,directory,labelTop,labelBot)
            dif = drawDifference(x,operonName,dictionary,field+graphName[index],index,indices,isTime,directory,label)
#            print (320,dif)
            
            for i in range(len(operonName)):
                differences[i]+=dif[i]
                output[i][index] += dif[i]
        drawDifferenceSum(x,operonName,dictionary,"reference",index,indices,isTime,directory,differences,label)
    else:
        drawOne(x,operonName,dictionary,"Time(log10)",None,indices,isTime,directory,labelTop,labelBot)
        output = drawDifference(x,operonName,dictionary,"Time(log10)",None,indices,isTime,directory,label)
    return operonName,output
    
###############################################################################
## Main program
###############################################################################
if __name__ == "__main__":
    args        = parse_args()
    naiveInput  = args.naiveInput
    approxInput = args.approxInput
    ILPInput    = args.ILPInput
    LPInput     = args.LPInput
    output_file = args.output
    geneBlock   = args.geneBlock

    # parse object for naive and approx
    dictionary1  = parseData(naiveInput,approxInput,geneBlock)
    # draw graph
    analysis = output_file+"NaiveApproxReference"
    try:
        os.mkdir(analysis)
    except:
        print ("directory analysis is already created")
    print ("Time")
    drawAll(dictionary1["time"],True,analysis,None,"Differences(Naive - Approx)","Naive","Approx")
#    drawAll(dictionary["pairWiseDistance"],False,analysis,"pairWise")
    print ("Events")
    drawAll(dictionary1["referenceDistance"],False,analysis,"reference","Differences(Naive - Approx)","Naive","Approx")
    
    # compare ILP and naive
    dictionary2  = parseData(naiveInput,ILPInput,geneBlock)
    # draw graph
    analysis = output_file+"NaiveILPReference"
    try:
        os.mkdir(analysis)
    except:
        print ("directory analysis is already created")
    print ("Time")
    operonName,timeDifferencesNaiveILP = drawAll(dictionary2["time"],True,analysis,None,"Differences(Naive - ILP)","Naive","ILP")
#    drawAll(dictionary["pairWiseDistance"],False,analysis,"pairWise")
    print ("Events")
    operonName,differencesNaiveILP = drawAll(dictionary2["referenceDistance"],False,analysis,"reference","Differences(Naive - ILP)","Naive","ILP")
    
    # compare ILP and approx
    dictionary3  = parseData(approxInput,ILPInput,geneBlock)
    # draw graph
    analysis = output_file+"ApproxILPReference"
    try:
        os.mkdir(analysis)
    except:
        print ("directory analysis is already created")
    print ("Time")
    operonName,timeDifferencesApproxILP = drawAll(dictionary3["time"],True,analysis,None,"Differences(Approx - ILP)","Approx","ILP")
#    drawAll(dictionary["pairWiseDistance"],False,analysis,"pairWise")
    print ("Events")
    operonName,differencesApproxILP = drawAll(dictionary3["referenceDistance"],False,analysis,"reference","Differences(Approx - ILP)","Approx","ILP")
    
    print (differencesNaiveILP)
    print ("**************************")
    print (differencesApproxILP)
    
    dictionary4 = {}
    for index,operon in enumerate(operonName):
        dictionary4[operon] = {"naive":differencesNaiveILP[index],"approx":differencesApproxILP[index]}
        dictionary4[operon]["genes"] = dictionary1["referenceDistance"][operon]["genes"]
    analysis = output_file+ "NaiveILPvsApproxILP"
    try:
        os.mkdir(analysis)
    except:
        print ("directory analysis is already created")    
    print ("Events")
    drawAll(dictionary4,False,analysis,"reference","Differences(NaiveILP - ApproxILP)","NaiveILP","ApproxILP")