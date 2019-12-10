# -*- coding: utf-8 -*-
"""
Created on Mon Sep  9 14:03:35 2019

@author: huyn
"""

from file_handle import traverseAll,parsing
import json
naive_result = "new_result"
approx_result = "new_result_approx"
def writeFile(inDirectory, outfile):
    dictionary = {}
    for file in traverseAll(inDirectory):
        operonName = file.split("/")[-1]
        mapping,genomes = parsing(file)
        distance = computeDistance(genomes)
        dictionary[operonName] = distance
    with open(outfile,"w") as outfile:
        json.dump(dictionary,outfile)  
def computeDistance(dictionary):
    distance = [0,0,0]
    myList   = list(dictionary.keys())
    for i in range(len(myList)-1):
        currentGeneBlock = dictionary[myList[i]]
        for j in range(i+1,len(myList)):
            nextGeneBlock =dictionary[myList[j]]
            deletion     = computeDeletion(currentGeneBlock,nextGeneBlock)
            split        = computeSplit(currentGeneBlock,nextGeneBlock)
            distance[0]+=deletion
            distance[2]+=split
    return distance
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
# run the main
writeFile(naive_result,"naive_result.txt")
writeFile(approx_result,"approx_result.txt")