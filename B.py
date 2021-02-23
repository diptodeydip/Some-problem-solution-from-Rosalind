#Functions


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

    