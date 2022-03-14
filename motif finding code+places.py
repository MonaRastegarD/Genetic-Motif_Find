# -*- coding: utf-8 -*-
"""
Created on Sun Feb  6 01:36:26 2022

@author: tahsin - rastegar
"""

import random 

def randomizedmotifsearch(lines): #Finds random l_mer
    global sTable
  
    sTable=[]
    motifs = []
  
    for n in lines:
       i = random.randint(0, 365)
       lmer=n[i:i+35]
       motifs.append(lmer)     

       sTable.append({
            "seq": lmer,
            "place": i,                      
        })
   
    for i in range(len(sTable)):
        print("p",i+1,"  ",sTable[i]["seq"] )   
    
    return  sTable
#**************************************
def crossover(lines,sTable): #random cut point crossover
    
    Table=[]
       
    for n in range(25):
       k = random.randint(2, 33) 
       
       i = random.randint(1, len(sTable)-1)  #selecting parent
       j = random.randint(1, len(sTable)-1) 
       
       p1 = sTable[i]["seq"]
       p2 = sTable[j]["seq"]
       
       #interchanging the genes
       o1=p1[0:k+1]+p2[k+1:k+35]    #creating offspring
       o2=p2[0:k+1]+p1[k+1:k+35]
                    
       Table.append({
            "seq": o1,
            "place": sTable[i]["place"],                       
        })
       
       Table.append({
            "seq": o2,
            "place": sTable[j]["place"],                       
        })
 
    print('Offsprings '.center(60, '-'))
    for i in range(len(Table)):
        print("O",i+1,"  ",  Table[i]["seq"])
   
    sTable = Table    
    
    return  Table

#************************************************************
# calculation of FS=Fitness score and TFS=Total fitness score
# Counts the Number of matching characters in a pair of string

def Score(sTable,lines):
   
    l=[]
    seqs= [sTable[i]["seq"] ]
    FSlist=[]
    TFSlist=[]
    
    for line in lines:
       
       sj=[]
       for j in range(364):    
           sj.append(line[j:j+35])    
  #  print("len sj",len(sj)) 
  #  print("sj",sj)
        
   #FS calculation with match scores -> match=1   mismatch=0  
    FS=[]
    TFS=[]
    for seq in seqs: 
       for k in range(len(sj)):     
          s1= seq
          s2= sj[k]
          a = sum(1 if i == k else 0 for i, k in zip(s1, s2)) / 35
          FSlist.append(a)        
          FS=max(FSlist)               
          index=FSlist.index(FS)         
          indexmax = seq 
           #l.append([indexmax,index])
   
    #TFS calculation with match scores
    TFSlist.append(sum(FSlist))
    TFS=sum(TFSlist)
    l.append([indexmax,index,TFS])
    
    """print("TFS= ",TFS)
    print("FS values: ",FSlist)
    for d in range(len(FSlist)):
       print("%6.4f "%FSlist[d])"""
            
    print(" result ".center(60,"-"))
    #print("FS = %6.4f  "%FS,"  TFS= %6.4f"%TFS)
    print("Potential motif : ",indexmax ,"\tscore= %6.4f"%TFS," \n")   
    return l

#**********************************************
def motiffind(bestmotif):
    import difflib
    motif=bestmotif[i]["seq"]
    match_strings=[]
    #places=[]
    table=[]
    for line in lines:        
       sj=[]
       for j in range(364):    
            sj.append(line[j:j+35])  
         
       #n is the maximum number of close matches to return in each line   
       match_strings=difflib.get_close_matches(motif,sj,n=1,cutoff=0.6)[0]     
       t=lines.index(line)
       #print("%d"%t,match_strings)  
       index = line.find(match_strings)   
       table.append({
                    "seq": match_strings,
                    "line" : t,
                    "place": index
                    })   
    
    #print(table) 
    for t in range(len(lines)):
      
      #return positions starting from 1
        print(" line: ",table[t]["line"]+1,"\t\tplace:",table[t]["place"]+1)   
    return           
    

#**********************************************
if __name__ == "__main__":
    
    motif_lenght=35    
    bestmotif=[]
    with open('motif_seq_main.txt') as f:
        lines = f.readlines()
    f.close()
    iteration_number=3
    
    for i in range(iteration_number):
        # print("")
         print(f'iteration {i+1} '.center(60, '-'))
         print("")
         seqs = randomizedmotifsearch(lines)
         l=Score(seqs,lines)
         seqs = crossover(lines,sTable)
         l=l+Score(seqs,lines)  
         
    m=max(map(lambda x: x[2], l)) #finds max score to select motif
         #l=[indexmax,index,TFS]
    for seq in l:
        if seq[2] == m:
            bestmotif.append({
                      "seq"  : seq[0],
                      "place": seq[1], 
                      "score": seq[2]
                       })
    
    print(" final result ".center(60,"-"))            
    for i in range(len(bestmotif)):
         print("best motif is:", bestmotif[i]["seq"],
               "best score: %6.4f "%bestmotif[i]["score"], sep = '\t') 
    print("")          
    print(" places of best motif ".center(60,"-"),"\n")
    
    motiffind(bestmotif)
       
    

