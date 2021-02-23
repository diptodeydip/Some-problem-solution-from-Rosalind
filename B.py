#Functions
def occurrences(string, pattern):
    ls  = []
    indx = 0
    while True:
        indx = string.find(pattern, indx)
        if indx != -1:
            ls.append(indx)
            indx+=1
        else:
            return ls
        
def findClumps(string, pattern,l,t):
    ls = occurrences(string,pattern)
    queue = []
    for x in ls:
        try:
            while x-queue[0]+1>l:  #eta correct while x-queue[0]+len(pattern)>l: len(pattern) add korte hoise karon pattern
                                    # jekhane paoa gese sei index diye difference ber kortesi
                queue.pop(0)
        except:
            pass
        queue.append(x)    
        if len(queue)>=t:
               return True
        
    return False

def hammingDis(p,q):
    length = len(p)
    hm = 0
    for i in range(length):
        if p[i]!=q[i]:
            hm+=1
    return hm

def genPattern(k,ls,txt):
    if k == 0:
        ls.append(txt)
        return
    arr = ['A','C','G','T']
    for x in arr:
        genPattern(k-1,ls,txt+x)
        
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
        flag = {'A':0,'C':0,'G':0,'T':0} #flag = {'A':0.0,'C':0.0,'G':0.0,'T':0.0}
        for j in range(len(motifSet)):
            flag[motifSet[j][i]]+=1
        profile[0].append(flag['A']/len(motifSet))
        profile[1].append(flag['C']/len(motifSet))
        profile[2].append(flag['G']/len(motifSet))
        profile[3].append(flag['T']/len(motifSet))
    return profile

def findProfileMostProbableKmer(dna,k,profile):
    mx = 0.0
    kmer = "none"

    for i in range(len(dna)-k+1):
        temp = dna[i:i+k]
        score = findProbability(temp,profile)
        if score > mx :
            mx = score
            kmer = temp
        elif score == mx:
            if kmer == "none":
                kmer = temp
    return kmer
                
                


#ba1f
def ba1f():
    dna = input("Enter sequence: ")
    g = c = 0
    ls = [0]

    for x in dna:
        if x == 'G':
            g+=1
        elif x == 'C':
            c+=1
        ls.append(g-c)
    mn = min(ls)

    for i in range(len(ls)):
        if ls[i]==mn:
            print(i)
#ba1e
def ba1e():
    data = input("Enter a string: ")
    dataLen = len(data)
    k,l,t = input("Enter k, l , t : ").split()
    k = int(k)
    l = int(l)
    t = int(t)

    store = {}


    for i in range(dataLen-k+1):
        tempPattern = data[i:i+k]
        if tempPattern not in store:
            store[tempPattern] = True
            x = findClumps(data,tempPattern,l,t)
            if x:
                print(tempPattern)

#ba2b
def ba2b():
    k = input("Enter k : ")
    k = int(k)
    data = [j for j in input("Enter strings: ").split()]

    # while True:
    #     x = input("Enter string or press enter to execute : ")
    #     if x=='':
    #         break
    #     else:
    #         data.append(x)

    ls = []
    genPattern(k,ls,"")

    mn = 100000000
    medianS = "none"

    for pattern in ls:
        totalDis= 0
        for string in data:
            tempMn = 100000000
            for i in range(len(string)-k+1):
                x = string[i:i+k]
                hD = hammingDis(pattern,x)
                if hD <= tempMn:
                    tempMn = hD
            totalDis+=tempMn
        if totalDis <= mn:
            mn = totalDis
            medianS = pattern

    print(medianS)
                
#ba2c
def ba2c():
    dna = input("Enter sequence: ")
    k = input("Enter k : ")
    k = int(k)

    profile = inputMatrix(4,k)


    kmer = findProfileMostProbableKmer(dna,k,profile)

    print(kmer)

#ba2d
def ba2d():
    k,t = input("Enter k , t: ").split()
    k = int(k)
    t = int(t)

    data = [j for j in input("Enter {0} strings: ".format(t)).split()]

    bestMotif = []
    for x in data:
        bestMotif.append(x[0:k])

    for i in range(len(data[0])-k+1):
        kmer = data[0][i:i+k]
        tempMotif = []
        tempMotif.append(kmer)
        for j in range(1,t,1):
            profile = createProfile(tempMotif,k)
            pmpKmer = findProfileMostProbableKmer(data[j],k,profile)
            tempMotif.append(pmpKmer)
        if findScore(tempMotif,k,t) < findScore(bestMotif,k,t):
            bestMotif = tempMotif

    for x in bestMotif:
        print(x)

    