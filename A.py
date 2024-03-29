#Functions

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


def hammingDis(p,q):
    length = len(p)
    hm = 0
    for i in range(length):
        if p[i]!=q[i]:
            hm+=1
    return hm

def subStrPrint(pattern,dna,d):
    length = len(dna)
    pattLen = len(pattern)
    for i in range(length-pattLen+1):
        temp = dna[i:i+pattLen]
        tempD = len([[i,j] for i,j in zip(pattern, temp) if i!=j])# or use hammingDis() method
        if tempD<=d:
            print(i)
            
def occurrences(pattern,dna,d):
    length = len(dna)
    pattLen = len(pattern)
    cnt = 0
    for i in range(length-pattLen+1):
        temp = dna[i:i+pattLen]
        tempD = len([[i,j] for i,j in zip(pattern, temp) if i!=j])# or use hammingDis() method
        if tempD<=d:
            cnt+=1
            if(pattern == "ATGT" or pattern == "ACAT"):
                print(i)
    return cnt

def genPattern(k,ls,txt):
    if k == 0:
        ls.append(txt)
        return
    arr = ['A','C','G','T']
    for x in arr:
        genPattern(k-1,ls,txt+x)
    

#ba1c
def ba1c():
    dna = input("Enter your sequence: ")
    revCompDna = dnaRevCompli(dna)
    print(revCompDna)


#ba1g
def ba1g():
    dna1 = input("Enter sequence one: ")
    dna2 = input("Enter sequence two: ")
    hm = hammingDis(dna1,dna2)
    print(hm)

#ba1h
def ba1h():
    pattern = input("Enter pattern: ")
    dna = input("Enter sequence: ")
    d = int(input("Enter d: "))
    subStrPrint(pattern,dna,d)

#ba1i
def ba1i():
    dna = input("Enter sequence: ")
    k,d = input("Enter k and d: ").split()
    k= int(k)
    d= int(d)

    mx = 0
    store = []
    patterns = ["none"]

    genPattern(k,store,"")

    for i in store:
            cnt = occurrences(i,dna,d)
            if cnt > mx:
                mx = cnt
                patterns.clear()
                patterns.append(i)
            elif cnt == mx:
                patterns.append(i)
        
    print("\n")
    for x in patterns:
        print(x)

#ba1j
def ba1j():
    dna = input("Enter sequence: ")
    k,d = input("Enter k and d: ").split()
    k= int(k)
    d= int(d)

    mx = 0
    store = []
    patterns = ["none"]
    flag = {}

    genPattern(k,store,"")

    for i in store:
        if i not in flag:
            cnt = occurrences(i,dna,d) + occurrences(dnaRevCompli(i),dna,d)
            flag[i]= True
            flag[dnaRevCompli(i)]= True
            if cnt > mx:
                mx = cnt
                patterns.clear()
                patterns.append((i,dnaRevCompli(i)))
            elif cnt == mx:
                patterns.append((i,dnaRevCompli(i)))

    print(mx)
    print("\n")
    for x,y in patterns:
        print(x+" "+y)

#ba1l
def ba1l():
    dna = input("Enter sequence: ")
    store = {'A':0,'C':1,'G':2,'T':3}
    dna = dna[::-1]
    indx = 0
    for i in range(len(dna)):
        indx+=store[dna[i]]*pow(2,2*i)
        
    print(indx)
    
#ba1m
def ba1m():
    indx =  int(input("Enter index: "))
    k =  int(input("Enter k: "))
    store = {0:'A',1:'C',2:'G',3:'T'}
    txt = ''

    for i in range(k-1,-1,-1):
        div = int(indx/pow(2,2*i))
        indx = indx%pow(2,2*i)
        txt+=store[div]

    print(txt)
