#Functions
def getAminoAcid(st):
    if st =="UCG" or st =="UCA" or st =="UCC" or st =="UCU":
        return "S"
    elif st=="UAC" or st=="UAU":
        return "Y"
    elif st=="UGG":
        return "W"
    elif st=="UGC" or st=="UGU":
        return "C"
    elif st=="UUG" or st=="UUA":
        return "L"
    elif st=="UUC" or st=="UUU":
        return "F"
    elif st=="GUG" or st=="GUA" or st=="GUC" or st=="GUU":
        return "V"
    elif st=="GGG" or st=="GGA" or st=="GGC" or st=="GGU":
        return "G"
    elif st=="GCG" or st=="GCA" or st=="GCC" or st=="GCU":
        return "A"
    elif st=="GAA" or st=="GAG":
        return "E"
    elif st=="GAU" or st=="GAC":
        return "D"
    elif st=="CCG" or st=="CCA" or st=="CCC" or st=="CCU":
        return "P"
    elif st=="CGG" or st=="CGA" or st=="CGC" or st=="CGU":
        return "R"
    elif st=="CUG" or st=="CUA" or st=="CUC" or st=="CUU":
        return "L"
    elif st=="CAU" or st=="CAC":
        return "H"
    elif st=="CAA" or st=="CAG":
        return "Q"
    elif  st=="AUA" or st=="AUC" or st=="AUU":
        return "I"
    elif st=="AUG":
        return "M"
    elif st=="ACG" or st=="ACA" or st=="ACC" or st=="ACU":
        return "T"
    elif st=="AGG" or st=="AGA":
        return "R"
    elif st=="AGU" or st=="AGC":
        return "S"
    elif st=="AAG" or st=="AAA":
        return "K"
    elif st=="AAU" or st=="AAC":
        return "N"
    else:
        return False
    
def dnaRevCompli(st):
    st = st[::-1]
    length = len(st)
    for i in range(length):
        if st[i]=="A" or st[i]=='a':
            st = st[:i]+'T'+st[i+1:]
        elif st[i]=='T' or st[i]=='t':
            st = st[:i]+'A'+st[i+1:]
        elif st[i]=='G' or st[i]=='g':
            st = st[:i]+'C'+st[i+1:]
        else:
            st = st[:i]+'G'+st[i+1:]
    return st

def transcribe(st):
    length = len(st)
    for i in range(length):

        if st[i]=='T' or st[i]=='t':
            st = st[:i]+'U'+st[i+1:]

    return st

def getMass(st):
    dict = {"G":57,"A":71,"S":87,"P":97,"V":99,"T":101,"C":103,"I":113,"L":113,"N":114,"D":115,"K":128,"Q":128,"E":129,"M":131,
           "H":137,"F":147,"R":156,"Y":163,"W":186}
    mass = 0
    for x in st:
        mass+=dict[x]
    return mass


#ba4a
def ba4a():
    rna = input("Enter RNA: ")

    for i in range(0,len(rna)-2,3):
        st = rna[i:i+3]
        x = getAminoAcid(st)
        if x == False:
            break
        print(x, end ="")
    
#ba4b
def ba4b():
    dna = input("Enter DNA: ")
    revDna = dnaRevCompli(dna)
    peptide = input("Enter peptide: ")
    rna = transcribe(dna)
    revRna = transcribe(revDna)

    rnaLen = len(peptide)*3

    for j in range(0,len(rna)-rnaLen+1):
        segment1 = rna[j:j+rnaLen]
        segment2 = revRna[j:j+rnaLen]

        x=""
        y=""

        for i in range(0,len(segment1)-2,3):
            st = segment1[i:i+3]
            t = getAminoAcid(st)
            if t!=False:
                x += t

        for i in range(0,len(segment2)-2,3):
            st = segment2[i:i+3]
            t = getAminoAcid(st)
            if t!=False:
                y += t

        if x == peptide:
            print(dna[j:j+rnaLen])
        if y == peptide:
            print(dnaRevCompli(revDna[j:j+rnaLen]))

    
#ba4c
def ba4c():
    peptide = input("Enter peptide: ")

    cycloPep = peptide + peptide

    print(0)
    for i in range(1,len(peptide)):
        for j in range(len(peptide)):
            temp = cycloPep[j:j+i]
            print(getMass(temp))

    print(getMass(peptide))
    
#ba4d
def ba4d():
    peptideMass = input("Enter peptideMass: ")
    peptideMass = int(peptideMass)

    count = 0
    dp = {}

    dict = {"G":57,"A":71,"S":87,"P":97,"V":99,"T":101,"C":103,"I":113,"N":114,"D":115,"Q":128,"E":129,"M":131,
               "H":137,"F":147,"R":156,"Y":163,"W":186} #I and L , K and Q are identical. so take them as one

    def recursion(currentMass):
        #global count
        if currentMass == peptideMass:
            dp[currentMass] = 1
            return 1
        elif currentMass > peptideMass:
            return 0
        elif currentMass in dp:
            return dp[currentMass]

        for x in dict:
            ls.append(x)
            y = recursion(dict[x]+currentMass)
            ls.remove(x)
            if currentMass in dp:
                dp[currentMass] += y
            else:
                dp[currentMass] = y

        return dp[currentMass]

    x = recursion(0)

    print(x)
    #print(count)
    
#ba3f
def ba3f():
    def parse_data(data, dict):
        for line in data:
            temp = line.split(" -> ")
            dict[temp[0]] = temp[1].split(",")
        return dict

    edges = []

    #while True:
        #x = input("Enter edges or press enter to execute : ")
        #if x=='':
            #break
       # else:
            #edges.append(x)

    f = open("data.txt", "r")
    for x in f:
        edges.append(x.rstrip('\n'))

    dict = {}
    dict = parse_data(edges,dict)

    # data input complete

    flag = {}
    stack = ['0']
    path = []

    #print(dict)


    def heirholzer(flag,stack,path): # uses dfs
        while stack:
            #print(stack)
            #print(stack[-1])
            boolean = False
            lastElem = stack[-1]
            for x in dict[lastElem]:
                if (lastElem,x) not in flag:
                    flag[(lastElem,x)] = 1
                    stack.append(x)
                    boolean = True
                    break

            if boolean == False:
                path.append(stack[-1])
                stack.pop(-1)
        return

    heirholzer(flag,stack,path)

    #print(flag)

    path.reverse()
    #print(path[0], end="")
    #for i in range(1,len(path)):
        #print(" -> {0}".format(path[i]), end="")

    with open('output.txt', 'w') as f:
        for item in path:
            f.write("%s->" % item) # have to remove last '->' from output file
    
    
#ba3h
def ba3h():
    k = input("Enter length of Kmer: ")
    k = int(k)
    data = [j for j in input("Enter  strings: ").split()]

    start = ""
    for item in data:
        flag = True
        for x in data:
            if item[0:k-1]==x[1:k]:
                flag = False
        if flag == True:
            start = item
            #print(start)
            break

    string = start
    while True:
        flag = True
        for item in data:
            if start[1:k]==item[0:k-1]:
                string+=item[-1]
                start = item
                flag = False
        if flag == True:
            break
    print(string)
    
#ba3k
def ba3k():
    data = [j for j in input("Enter  strings: ").split()]
    k = len(data[0])

    outDeg = {}
    inDegCount = {}
    inDeg = {}
    nodes = {}
    nodes = set()

    for x in data:
        if x[0:k-1] not in outDeg:
            outDeg[x[0:k-1]] =[x[1:k]]
        else:
            outDeg[x[0:k-1]].append(x[1:k])

        if x[1:k] not in inDeg:
            inDegCount[x[1:k]] = 1
            inDeg[x[1:k]] = [x[0:k-1]]
        else:
            inDegCount[x[1:k]]+=1
            inDeg[x[1:k]].append( x[0:k-1])
        nodes.add(x[0:k-1])
        nodes.add(x[1:k])

    #print(outDeg)
    #print(inDeg)
    #print(inDegCount)


    def dfs(st,currentNode):
        if currentNode in outDeg:
            if len(outDeg[currentNode])>1 or len(inDeg.get(currentNode,[]))>1:
                print(st)
                for x in outDeg[currentNode]:
                    if inDegCount[x]>0:
                        inDegCount[x]-=1
                        dfs(currentNode+x[-1],x)


            elif inDegCount[outDeg[currentNode][0]]>0:
                inDegCount[outDeg[currentNode][0]]-=1
                dfs(st+outDeg[currentNode][0][-1],outDeg[currentNode][0])

        else:
            print(st)
        return


    for x in nodes:
        if x not in inDeg:
            dfs(x,x)