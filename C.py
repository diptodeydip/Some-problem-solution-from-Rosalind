#Functions
def inputMatrix(row,col):
    mat = []
    for i in range(row):
        mat.append([float(j) for j in input().split()])
        
    return mat

def findProbability(temp,mat):
    store = {'A':0,'C':1,'G':2,'T':3}
    score = float(1)
    for i in range(len(temp)):
        score *= mat[store[temp[i]]][i]
    return score
    
def findScore(motifSet,k,t):
    score = 0
    for i in range(k):
        flag = {'A':0,'C':0,'G':0,'T':0}
        mx = 0
        for j in range(t):
            flag[motifSet[j][i]]+=1
            if flag[motifSet[j][i]]>=mx:
                mx = flag[motifSet[j][i]]
        score += t - mx
    return score

def createProfile(motifSet,k):
    profile = [[],[],[],[]]
    for i in range(k):
        flag = {'A':1.0,'C':1.0,'G':1.0,'T':1.0} # pseudocounts
        for j in range(len(motifSet)):
            flag[motifSet[j][i]]+=1.0
        profile[0].append(flag['A']/(len(motifSet)+4.0))
        profile[1].append(flag['C']/(len(motifSet)+4.0))
        profile[2].append(flag['G']/(len(motifSet)+4.0))
        profile[3].append(flag['T']/(len(motifSet)+4.0))
    return profile.copy()

def findProfileMostProbableKmer(dna,k,profile):
    mx = 0.0
    kmer = "none"

    for i in range(len(dna)-k+1):
        temp = dna[i:i+k]
        score = findProbability(temp,profile)
        if score >= mx :
            mx = score
            kmer = temp
        #elif score == mx:
            #if kmer == "none":
                #kmer = temp
    return kmer

def getMotifSet(data,k,profile):
    tempMotif = []
    for j in range(len(data)):
        pmpKmer = findProfileMostProbableKmer(data[j],k,profile)
        tempMotif.append(pmpKmer)
    return tempMotif.copy()
            

#ba2f
def ba2f():
    #RandomizedMotifSearch
    from random import randint

    k,t = input("Enter k , t: ").split()
    k = int(k)
    t = int(t)

    data = [j for j in input("Enter {0} strings: ".format(t)).split()]

    finalBestMotif = []

    i=1000
    #print(data)

    while i:
        i-=1
        bestMotif = []
        for x in data:
            rand = randint(0, len(data[0])-k)
            bestMotif.append(x[rand:rand+k])

        if i==999:
            finalBestMotif = list(bestMotif) # initialize finalBestMotif one time
        #print(bestMotif)    

        while True:
            profile = createProfile(bestMotif,k)
            tempMotifs = getMotifSet(data,k,profile)

            #print(bestMotif)  
            #print(tempMotifs)

            if findScore(tempMotifs,k,t) < findScore(bestMotif,k,t):
                bestMotif = list(tempMotifs)
                #print(tempMotif)
            else:
                break

        if findScore(bestMotif,k,t) < findScore(finalBestMotif,k,t):
            finalBestMotif = list(bestMotif)
        #print(i)

    for x in finalBestMotif:
        print(x)   
        

#ba2g
def ba2g():
    #GibbsSampler
    from random import randint

    k,t,n = input("Enter k , t , n: ").split()
    k = int(k)
    t = int(t)
    n = int(n)

    data = [j for j in input("Enter {0} strings: ".format(t)).split()]


    bestMotif = []
    repeat = 50

    while repeat:

        repeat-= 1
        i=n
        motifs = []
        for x in data:
            rand = randint(0, len(data[0])-k)
            motifs.append(x[rand:rand+k])

        if repeat==49:
            bestMotif = list(motifs)

        while i:
            i-=1
            ind = randint(0,t-1)
            #print(ind)
            tempMotifs = list(motifs)
            tempMotifs.pop(ind)

            profile = createProfile(tempMotifs,k)
            motifs[ind] = findProfileMostProbableKmer(data[ind],k,profile)

            if findScore(motifs,k,t) < findScore(bestMotif,k,t):
                bestMotif = list(motifs)

    for x in bestMotif:
        print(x)   
        

#ba3a
def ba3a():
    k= input("Enter k: ")
    k = int(k)
    data = input("Enter string: ")

    for i in range(len(data)-k+1):
        temp = data[i:i+k]
        print(temp)    

#ba3b
def ba3b():
    data = [j for j in input("Enter  strings: ").split()]

    dna = data[0]

    for i in range(1,len(data)):
        dna+= data[i][-1:]
    print(dna)
    
#ba3c
def ba3c():
    data = [j for j in input("Enter  strings: ").split()]

    for x in data:
        for y in data:
            if x[1:]==y[:-1]:
                print("{0} -> {1}".format(x,y))    

#ba3d
def ba3d():
    k= input("Enter k: ")
    k = int(k)
    data = input("Enter  strings: ")

    kmers = []

    for i in range(len(data)-k+1):
        kmers.append(data[i:i+k])

    dicto = {}

    for x in kmers:
        if x[:-1] not in dicto:
            ls = []
            ls.append(x[1:])
            dicto[x[:-1]] = ls
        else:
            dicto[x[:-1]].append(x[1:])

    for i in dicto:
        temp = dicto[i]
        print("{0} -> {1}".format(i,temp[0]), end ="")
        for j in range(1,len(temp)):
            print(",{0}".format(temp[j]), end ="")
        print()
    
    
#ba3e
def ba3e():
    kmers = [j for j in input("Enter  strings: ").split()]

    dicto = {}

    for x in kmers:
        if x[:-1] not in dicto:
            ls = []
            ls.append(x[1:])
            dicto[x[:-1]] = ls
        else:
            dicto[x[:-1]].append(x[1:])

    for i in dicto:
        temp = dicto[i]
        print("{0} -> {1}".format(i,temp[0]), end ="")
        for j in range(1,len(temp)):
            print(",{0}".format(temp[j]), end ="")
        print()
    
    