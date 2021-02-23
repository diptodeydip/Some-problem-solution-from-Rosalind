#Functions
def getMass(st,symbolValue):
    mass = 0
    for x in st:
        mass+=symbolValue[x]
    return mass

def getSymbol(st):
    dict = {"G":57,"A":71,"S":87,"P":97,"V":99,"T":101,"C":103,"I":113,"L":113,"N":114,"D":115,"K":128,"Q":128,"E":129,"M":131,
           "H":137,"F":147,"R":156,"Y":163,"W":186}
    for  x in dict:
        if dict[x]==st:
            return x
    return "none"

def score(peptide , spectrum,symbolValue):
    cycloPep = peptide + peptide

    score = 0
    if 0 in spectrum:
        spectrum.remove(0)
        score+=1
    if getMass(peptide,symbolValue) in spectrum:
        score+=1
        spectrum.remove(getMass(peptide,symbolValue))
    
    for i in range(1,len(peptide)):
        for j in range(len(peptide)):
            temp = cycloPep[j:j+i]
            #print(getMass(temp))
            if getMass(temp,symbolValue) in spectrum:
                score+=1
                spectrum.remove(getMass(temp,symbolValue))
    return score




def expand(aminoAcids, leaderBoard):
    newLeaderBoard = []
    for x in leaderBoard:
        for acid in aminoAcids:
            newLeaderBoard.append(x+acid)
    return newLeaderBoard


    
def cut(leaderBoard ,spectrum, N,symbolValue):
    dict = {}
    for item in leaderBoard:
        dict[item] = score(item,spectrum.copy(),symbolValue)
    
    dict = sorted(dict.items(), key=lambda x: x[1], reverse=True) # sort descending order (reverse=true) #returns list of tuple
    leaderBoard.clear()
    i = 1
    #print(dict)
    nthELement = 0
    for x in dict: # x is a tuple from dict and dict is now a list
        if i <= N:
            leaderBoard.append(x[0])
            nthELement = x[1]
        elif i>N:
            if x[1]==nthELement:
                leaderBoard.append(x[0])
            else:
                break
        i+=1
    return leaderBoard




    
def leaderBoardCyclopeptide(spectrum, N,aminoAcids,symbolValue):
    leaderBoard = [""]
    leaderPeptide = ""
    parentMass = 0
    for acid in spectrum:
        if acid > parentMass:
            parentMass = acid
    
    flagLeaderBoard = [""]

    while leaderBoard:
        leaderBoard = expand(aminoAcids,leaderBoard)
        flagLeaderBoard = leaderBoard.copy()
   
        for peptide in leaderBoard:
            if getMass(peptide,symbolValue) == parentMass:
                if score(peptide,spectrum.copy(),symbolValue) > score(leaderPeptide , spectrum.copy(),symbolValue):
                    leaderPeptide = peptide
            elif getMass(peptide,symbolValue) > parentMass:
                flagLeaderBoard.remove(peptide)
        
        leaderBoard = flagLeaderBoard.copy()
        leaderBoard = cut(leaderBoard,spectrum.copy() , N,symbolValue)

    return leaderPeptide

def conv(spectrum):
    
    cnv = []

    for item in spectrum:
        for x in spectrum:
            if x > item:
                cnv.append(x-item)
    return cnv
    

#ba4f
def ba4f():
    peptide = input("Enter peptide: ")

    spectrum = [int(j) for j in input("Enter masses: ").split()]

    sc = score(peptide,spectrum.copy())
    print(sc)

    
#ba4h
def ba4h():
    spectrum = [int(j) for j in input("Enter masses: ").split()]

    cnv = conv(spectrum.copy())
    dict = {}

    for item in cnv:
        if item not in dict:
            dict[item] = 1
        else:
            dict[item]+=1

    dict = sorted(dict.items(), key=lambda x: x[1], reverse=True)  #sort dictionary reverse true for descending order
    # now dict is list of tuple

    for x in dict: # x is a tuple from dict
        for i in range(0,x[1]):
            print("{0} ".format(x[0]),end="")
    
#ba4g
def ba4g():
    N = input("Enter N: ")
    N = int(N)
    spectrum = [int(j) for j in input("Enter masses: ").split()]

    aminoAcids = ["G","A","S","P","V","T","C","I","N","D","Q","E","M",
               "H","F","R","Y","W"]
    symbolValue = {"G":57,"A":71,"S":87,"P":97,"V":99,"T":101,"C":103,"I":113,"L":113,"N":114,"D":115,"K":128,"Q":128,"E":129,"M":131,
               "H":137,"F":147,"R":156,"Y":163,"W":186}

    leaderPeptide = leaderBoardCyclopeptide(spectrum,N,aminoAcids,symbolValue)
    for x in leaderPeptide:
        print("{0}-".format(getMass(x,symbolValue)),end="") # have to remove last '-'

    
#ba4i
def ba4i():
    #convLeaderBoardCycloPeptide
    M = input("Enter M: ")
    M = int(M)
    N = input("Enter N: ")
    N = int(N)
    spectrum = [int(j) for j in input("Enter masses: ").split()]

    symbolValue = {"G":57,"A":71,"S":87,"P":97,"V":99,"T":101,"C":103,"I":113,"L":113,"N":114,"D":115,"K":128,"Q":128,"E":129,"M":131,
               "H":137,"F":147,"R":156,"Y":163,"W":186}

    #################################
    cnv = conv(spectrum.copy())
    dict = {}
    for item in cnv:
        if item>=57 and item<=200:
            if item not in dict:
                dict[item] = 1
            else:
                dict[item]+=1

    dict = sorted(dict.items(), key=lambda x: x[1], reverse=True)  #sort dictionary reverse true for descending order
    #################################
    i = 1
    aminoAcids = []
    mthELement = 0

    reserveSym = 97 # small a
    checkSym = {}
    for x in dict: # x is a tuple from dict and dict is now a list
        sym = getSymbol(x[0])
        if i <= M:
            if sym == "none":
                aminoAcids.append(chr(reserveSym))
                symbolValue[chr(reserveSym)] = x[0]
                reserveSym+=1
            else:
                aminoAcids.append(sym)

            mthELement = x[1]
        elif i>N:
            if x[1]==mthELement:
                if sym == "none":
                    aminoAcids.append(chr(reserveSym))
                    symbolValue[chr(reserveSym)] = x[0]
                    reserveSym+=1
                else:
                    aminoAcids.append(sym)
            else:
                break
        i+=1

    ####################################

    leaderPeptide = leaderBoardCyclopeptide(spectrum,N,aminoAcids,symbolValue)
    for x in leaderPeptide:
        print("{0}-".format(getMass(x,symbolValue)),end="") # have to remove last '-'
    
#ba9g
def ba9g():
    string = input("Enter string: ")

    suffixList = []


    for i in range(len(string)):
        temp = string[i:len(string)]
        tuppple = (i,temp)
        suffixList.append(tuppple)

    #print(suffixList)

    suffixList.sort(key = lambda x: x[1])

    #print(suffixList)
    #for x in suffixList:
        #print(" {0},".format(x[0]),end="") # have to remove last , and first space

    with open('output.txt', 'w') as f:
        for x in suffixList:
            f.write("%d, " % x[0]) # have to remove last , and space
    
#ba9a
def ba9a():
    #creating trie
    strings = [j for j in input("Enter strings: ").split()]

    pNode = 0
    cNode = 1
    dict = {}


    for item in strings:

        found = False
        pNode = 0

        for i in range(len(item)):

            if item[0:i+1] in dict:
                found = True
                pNode = dict[item[0:i+1]]   # updating parentNode number
                pass

            else:
                print("{0}->{1}:{2}".format(pNode,cNode,item[i]))
                dict[item[0:i+1]] = cNode
                pNode=cNode
                cNode+=1

            
    
#ba9b
def ba9b():
    #trie matching
    text = input("Enter text: ")

    strings = [j for j in input("Enter strings: ").split()]

    pNode = 0
    cNode = 1
    dict = {}


    for item in strings:
        dict[item] = True

    for i in range(len(text)):
        for j in range(i+1,len(text)+1):
            temp = text[i:j]
            if temp in dict:
                print("{0} ".format(i),end="")
