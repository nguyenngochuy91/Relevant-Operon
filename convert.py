#!/usr/bin/env python
''' Author  : Huy Nguyen, and David C.Ream
    Program : Given the operon directory, for each operons file, get the info about the gene in each genomes, map it
              into alphabet letter , get the gap, map gap to '|' and write to ouputfile.
    Start   : 05/04/2016
    End     : 05/05/2016
'''

import os
import argparse
import time
import uuid
import json
import pulp # interger linear programming 
# traverse and get the file
def traverseAll(path):
    res=[]
    for root,dirs,files in os.walk(path):
        for f in files:
            res.append(root+f)
    return res

class readable_dir(argparse.Action):
    def __call__(self,parser, namespace, values, option_string=None):
        prospective_dir=values
        if not os.path.isdir(prospective_dir):
           try:
               os.mkdir(prospective_dir)
           except OSError:
               print (argparse.ArgumentTypeError("readable_dir:{0} is not a readable dir".format(prospective_dir)))
        if os.access(prospective_dir, os.R_OK):
            setattr(namespace,self.dest,prospective_dir)
        else:
            raise argparse.ArgumentTypeError("readable_dir:{0} is not a readable dir".format(prospective_dir))

def get_arguments():
    parser = argparse.ArgumentParser()
    parser.add_argument("--OperonDataDirectory","-i",action=readable_dir,help="This directory should contain files with gene name, start, stop, strand direction information for each genome.")
    parser.add_argument("--splitDistance","-d", type=int ,default=500,help="Splitting distance")
    parser.add_argument("--OutputDirectory","-o", help="Output of this program will be stored in the path supplied here. It will make a new directory if path given is valid or it will raise an error")
    parser.add_argument("--approx","-a", help="Using approx method (Y,N)",default = "N")      
    args = parser.parse_args()
    return args


def chk_output_directory_path(OutputDirectory,sessionID):
    if not os.path.exists(OutputDirectory + "_" + str(sessionID)):
        try:
           #os.mkdir(OutputDirectory + "_" + str(sessionID))
           return True
        except OSError:
           print ("Unable to create directory:", OutputDirectory)
           sys.exit()

# convert the file into dictionary with useful info    
def toDict(file):
    infile = open(file,'r')
    map_code=''
    mapping =''
    dic_map ={}
    main_dic={}
    # create 3 main key
    line = infile.readline()
    map_code = line
    mapping = line.split('\t')[:-1] # 
    for item in mapping:
        item_split = item.split(',')
        key = item_split[0]
        value = item_split[1]
        dic_map[key]=value # {'astA': 'a'}
    for line in infile.readlines():
        genome = line.split(':')[0] # (line.split(':')= ['NC_002696', '(astA,634700,635744,1)\t(astD,635730,637149,1)\t(astB,637145,638426,1)\t(astC,638435,639614,1)\t']
        main_dic[genome]={}
        # 3 distance sub dictionary
        main_dic[genome]['+1']={}
        main_dic[genome]['-1']={}
        main_dic[genome]['1']={}
        genes_string = line.split(':')[1]
        # to deal with each genes, they are in tuple, where first is the name of the gene, follow by the position, and the strand it is on
        # should consider 2 type of strand (so i can put a gap correctly
        genes_string = genes_string.split('\t')[:-1] # ['(astA,634700,635744,1)', '(astD,635730,637149,1)', '(astB,637145,638426,1)', '(astC,638435,639614,1)']
        genes_string = list(set(genes_string))
        for item in genes_string:
            info= item.split(',') #['dppA', '402362', '400796', '+1']
            position=(int(info[1]),int(info[2]))
            position=(min(position),max(position))
            main_dic[genome][info[3]][position]=dic_map[info[0]]
    return (main_dic,map_code)

# from dic, create string
def toString(dic,map_code,f,timeDictionary):
    mapping = map_code.split('\t')[:-1] 
 
    reference = ""
    for item in mapping:
        item_split = item.split(',')
        key = item_split[0]
        value = item_split[1]
        reference+=value # {'astA': 'a'}
#    print (dic)
    S = set(reference)
    wholestring=''
    wholestring+=map_code
    newWhole   = map_code 
    accumulate = 0
    ILPString = ""
    timeILP = 0
    LPString = ""
    timeLP  = 0
    for genome in dic:
        whole = genome + ':'
        string= genome + ':' # the string to be written
        LPString += genome +":"
        ILPString += genome +":"
        substring = ''
        flag = False # check if it has a gene block
        for key in dic[genome]:
            if len(dic[genome][key])==0:
                continue
            else:
                myList=[]
                for position in dic[genome][key]:
                    myList.append(position)
                myList.sort()
                substring += dic[genome][key][myList[0]]
                for index in range(len(myList)-1):
                    dif = abs(myList[index+1][0] - myList[index][1]) 
                    if dif >500:
                        substring += '|'
                    else:
                        flag = True   
                    substring += dic[genome][key][myList[index+1]]
                substring += '|' # changing strand
        if flag:
            string += substring[:-1] # only add if there is a gene block
#            print (130,string)
            blocks = set(substring[:-1].split("|"))
            C= set()
            for block in blocks:
                C.add(block)
#            print ("C",C)
            start = time.time()
            # solve using approx
            block = approxSolve(S,C)
            stop = time.time()
            accumulate+= (stop-start)/1000
            if check(block):
                whole+="|".join(block)
            else:
                whole+=""
                
            start = time.time()
            ILPBlock = ILPSolve(blocks)
            stop = time.time()
            timeILP+= (stop-start)/1000
            if check(ILPBlock):
                ILPString += "|".join(ILPBlock)
            
            start = time.time()
            LPBlock = LPSolve(blocks)
            stop = time.time()
            timeLP+= (stop-start)/1000
            if check(LPBlock):
                LPString += "|".join(LPBlock)
#            print ("S: {}, C: {} \napproxBlock: {} \nILPBlock: {} \nLPBlock: {}\n".format(S,C,block,ILPBlock,LPBlock))

        string += '\n'
        wholestring += string
        whole+='\n'
        newWhole+=whole
        LPString+= "\n"
        ILPString += "\n"
    timeDictionary["approx"][f]=accumulate
    timeDictionary["ILP"][f]=timeILP
    timeDictionary["LP"][f]=timeLP
#    print (timeDictionary)
    return wholestring,newWhole,ILPString,LPString
def check(block):
    for b in block:
        if len(b)>=2:
            return True
    return False
# solve using approx
def approxSolve(S,C):
    res = []
    cover = set()
    while cover!=S:
        maxSize = 0
        for item in C:
            subSet = set(item)-cover
            if len(subSet)>maxSize:
                maxSize = len(subSet)
                toAdd = item
            else:
                continue
                
        if maxSize == 0 or not C:
            break
        else:
            res.append(toAdd)
            cover.update(set(toAdd))
            C.remove(toAdd)
    return res
    
# solve using integer programing
def ILPSolve(C):
    S = set()
    for block in C:
        S.update(set(block))
    problem = pulp.LpProblem("Relevant Gene Block Integer", pulp.LpMinimize)
    # create a binary variable to state that a block is being used
    x = pulp.LpVariable.dicts("block",C,lowBound= 0, upBound= 1, cat = pulp.LpInteger)
    problem += sum([x[block] for block in C])
    # ensure that each gene in S is cover
    for gene in S:
        problem += sum(x[block] for block in C if gene in block) >=1, gene
    problem.solve()
    return [block for block in C if x[block].value() == 1.0]
# solve using relaxation
def LPSolve(C):
    S = set()
    for block in C:
        S.update(set(block))
    problem = pulp.LpProblem("Relevant Gene Block Integer", pulp.LpMinimize)
    # create a binary variable to state that a block is being used
    x = pulp.LpVariable.dicts("block",C,lowBound= 0, upBound= 1, cat = pulp.LpContinuous)
    problem += sum([x[block] for block in C])
    # ensure that each gene in S is cover
    for gene in S:
        problem += sum(x[block] for block in C if gene in block) >=1, gene
    problem.solve()
    maxCount = 0
    for gene in S:
        count = 0
        for block in C:
            if gene in block:
                count +=1
        maxCount = max(maxCount,count)
    return [block for block in C if x[block].value() >= 1/maxCount]    
if __name__ == "__main__":

    start = time.time()
    args = get_arguments()
    sessionID = uuid.uuid1()
    approx                     = args.approx
    outputsession = args.OutputDirectory
    res = traverseAll(args.OperonDataDirectory)
    # make a dictionary to compare 
    timeDictionary= {"approx":{},"ILP":{},"LP":{}}
    parentDir = outputsession.split("/")
    species = parentDir[-3]
    myDir = "/".join(parentDir[:-3])
    ILPDir =  "/".join(parentDir[:-3])+'_ILP/'+species
    LPDir = "/".join(parentDir[:-3])+'_LP/'+species
    try:
        os.mkdir("/".join(parentDir[:-3])+'_ILP/')
        os.mkdir(ILPDir)
        os.mkdir(ILPDir+"/new_result/")
    except:
        print ("Error 286")
    try:
        os.mkdir("/".join(parentDir[:-3])+'_LP/')
        os.mkdir(LPDir)
        os.mkdir(LPDir+"/new_result/")
    except:
        print ("Error 291")    
    for r in res:
        root,f = os.path.split(r)
#        print (root)
        result= toDict(r)
#        print (result)
        # compute both
        wholestring,whole,ILPString,LPString = toString(result[0],result[1],f,timeDictionary)
#        print (timeDictionary)
        if approx == "N":
            try:
                os.mkdir(outputsession)
            except:
                pass
            outfile = open(outputsession+'/'+f,'w')
            outfile.write(wholestring)
            outfile.close()
        else:  
            try:
                os.mkdir(outputsession)
            except:
                pass
    
            outfile = open(outputsession+'/'+f,'w')
            outfile.write(whole)
            outfile.close()
            
            outfile = open(ILPDir+"/new_result/"+f,'w')
            outfile.write(ILPString)
            outfile.close()
            outfile = open(LPDir+"/new_result/"+f,'w')
            outfile.write(LPString)
            outfile.close()
    if approx == "Y":
        root = "/".join(root.split("/")[:-1])
        with open(root+"/time.txt","w") as outfile:
            json.dump(timeDictionary["approx"],outfile)

        with open(myDir+'_ILP/'+species+"/time.txt","w") as outfile:
            json.dump(timeDictionary["ILP"],outfile)

        with open(LPDir+"/time.txt","w") as outfile:
            json.dump(timeDictionary["LP"],outfile)
