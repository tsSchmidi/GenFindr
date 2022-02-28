def fasta(msg):
    output=""
    #Instruction for user on what to input.
    print(msg)
    while True:
        #Looping through an input enables user to paste multiple rows of text.
        #Otherwise x=paste only registers the first row.
        temp=input().upper()
        if len(temp)==0:
            break
        #Split by spaces because numbers can be on the same row but with spacing.
        for i in temp.split(" "):
            #Discard input fragments that contain numbers etc.
            if i.isalpha():
                output+=i
    return(output)

def frames(dna,circular):
    if circular:
        #Circular DNA will be translated twice to ensure a full readout is possible from even the end.
        #Duplicates have to be removed later.
        output=2*dna,dna[1:]+dna,dna[2:]+dna
    else:
        output=dna,dna[1:],dna[2:]
    return(output)

def translate(frame):
    codons={"TTT":"F","TTC":"F","TTA":"L","TTG":"L","CTT":"L","CTC":"L","CTA":"L","CTG":"L","ATT":"I","ATC":"I","ATA":"I","ATG":"M","GTT":"V","GTC":"V","GTA":"V","GTG":"V","TCT":"S","TCC":"S","TCA":"S","TCG":"S","AGT":"S","AGC":"S","CCT":"P","CCC":"P","CCA":"P","CCG":"P","ACT":"T","ACC":"T","ACA":"T","ACG":"T","GCT":"A","GCC":"A","GCA":"A","GCG":"A","TAT":"Y","TAC":"Y","CAT":"H","CAC":"H","CAA":"Q","CAG":"Q","AAT":"N","AAC":"N","AAA":"K","AAG":"K","GAT":"D","GAC":"D","GAA":"E","GAG":"E","TGT":"C","TGC":"C","TGG":"W","CGT":"R","CGC":"R","CGA":"R","CGG":"R","AGA":"R","AGG":"R","GGT":"G","GGC":"G","GGA":"G","GGG":"G","TAA":"*","TAG":"*","TGA":"*"}
    output=""
    for pos in range(len(frame)//3):
        output+=codons[frame[pos*3:3+pos*3]]
    return(output)

def complementary(sequence):
    return("".join(reversed(sequence)).replace("A","t").replace("T","a").replace("G","c").replace("C","g").upper())

def find_all(template,find):
    n=0
    output=[]
    while True:
        n=template.find(find,n)
        if n==-1: break
        output+=n,
        n+=1
    return(output)

def compute(protein,dna,circular):
    #Create linear DNA reading frames
    f1,f2,f3=frames(dna,circular)
    #The entire DNA sequence is converted to codon outcomes (len/3)
    f1=translate(f1)
    f2=translate(f2)
    f3=translate(f3)
    #Repeat for reverse strand
    r1,r2,r3=frames(complementary(dna),circular)
    r1=translate(r1)
    r2=translate(r2)
    r3=translate(r3)

    f1hits=[]
    f2hits=[]
    f3hits=[]
    r1hits=[]
    r2hits=[]
    r3hits=[]
#For each frame, for each protein 10mer, return a list of positions (/3) where 10mer is found
#f1:[[pos1_all][pos2_all]...[posn_all]]. All is usually 1.
    for pos in range(len(protein)-9):
        f1hits+=find_all(f1,protein[pos:pos+10]),
        f2hits+=find_all(f2,protein[pos:pos+10]),
        f3hits+=find_all(f3,protein[pos:pos+10]),
        r1hits+=find_all(r1,protein[pos:pos+10]),
        r2hits+=find_all(r2,protein[pos:pos+10]),
        r3hits+=find_all(r3,protein[pos:pos+10]),
#Create empty lists for each key, a key being the common predicted start position of the gene
    result={}
    for ori in list(set([(n-i)*3 for i in range(len(f1hits)) for n in f1hits[i]])):
        if ori>len(dna):
            result[ori-len(dna)]=[]
        else:
            result[ori]=[]
    for ori in list(set([(n-i)*3+1 for i in range(len(f2hits)) for n in f2hits[i]])):
        if ori>len(dna):
            result[ori-len(dna)]=[]
        else:
            result[ori]=[]
    for ori in list(set([(n-i)*3+2 for i in range(len(f3hits)) for n in f3hits[i]])):
        if ori>len(dna):
            result[ori-len(dna)]=[]
        else:
            result[ori]=[]
    for ori in list(set([-(n-i)*3 for i in range(len(r1hits)) for n in r1hits[i]])):
        if -ori>len(dna):
            result[ori+len(dna)]=[]
        else:
            result[ori]=[]
    for ori in list(set([-(n-i)*3-1 for i in range(len(r2hits)) for n in r2hits[i]])):
        if -ori>len(dna):
            result[ori+len(dna)]=[]
        else:
            result[ori]=[]
    for ori in list(set([-(n-i)*3-2 for i in range(len(r3hits)) for n in r3hits[i]])):
        if -ori>len(dna):
            result[ori+len(dna)]=[]
        else:
            result[ori]=[]

#Convert 10mer hits (/3) into DNA position, correcting for frame and compiling into dictionary
    for i in range(len(f1hits)):
        for n in f1hits[i]:
            temp=n*3
            if temp>len(dna):
                temp+=-len(dna)
            result[temp-3*i]+=temp,
    for i in range(len(f2hits)):
        for n in f2hits[i]:
            temp=n*3+1
            if temp>len(dna):
                temp+=-len(dna)
            result[temp-3*i]+=temp,
    for i in range(len(f3hits)):
        for n in f3hits[i]:
            temp=n*3+2
            if temp>len(dna):
                temp+=-len(dna)
            result[temp-3*i]+=temp,
#Here, reverse strand positions converted into negative numbers to keep distinction
    for i in range(len(r1hits)):
        for n in r1hits[i]:
            temp=-n*3
            if -temp>len(dna):
                temp+=len(dna)
            result[temp+3*i]+=temp,
    for i in range(len(r2hits)):
        for n in r2hits[i]:
            temp=-n*3-1
            if -temp>len(dna):
                temp+=len(dna)
            result[temp+3*i]+=temp,
    for i in range(len(r3hits)):
        for n in r3hits[i]:
            temp=-n*3-2
            if -temp>len(dna):
                temp+=len(dna)
            result[temp+3*i]+=temp,

    print("Protein length: "+str(len(protein))+" codons (including stop codon)\n")
    #Loop through each key, i.e. possible gene or shifted fragment
    for key in result.keys():
        #Remove duplicates from reading circular DNA twice
        #result[key] is a list of integers
        result[key]=list(set(result[key]))
        result[key].sort()
    #Convert lists of consecutive numbers into start and end of each fragment
        start=[result[key][0]]
        end=[]
    #i is codon position, d[k][i] is DNA position. First position is skipped.
        for i in range(1,len(result[key])):
            if not result[key][i]==result[key][i-1]+3:
                end+=result[key][i-1],
                start+=result[key][i],
        end+=result[key][-1],
    #Now result[key] is a list of (fewer) pairs of integers
        result[key]=[[start[i],end[i]] for i in range(len(start))]
    #Position indexing translated from machine to human here. Also DNA -> codon.
        for fragment in result[key]:
            if key>0:
                print("Codons "+str(int(1+(fragment[0]-key)/3))+"-"+str(int(1+9+(fragment[1]-key)/3))+" in bases "+str(int(1+fragment[0]))+"-"+str(int(1+29+fragment[1]))+" (frame f"+str(int(fragment[0]%3+1))+") "+str(int(((fragment[1]-fragment[0])/3+10)/len(protein)*100))+"%")
            else:
                print("Codons "+str(int(1-(fragment[1]-key)/3))+"-"+str(int(1+9-(fragment[0]-key)/3))+" in bases "+str(int(len(dna)+fragment[1]))+"-"+str(int(len(dna)-29+fragment[0]))+" (frame r"+str(int(((-fragment[1]+1)-1)%3+1))+") "+str(int(((fragment[1]-fragment[0])/3+10)/len(protein)*100))+"%")

#Test if protein of interest is/was encoded on a DNA sequence
protein=fasta("Paste protein sequence and press enter twice")+"*"
print(protein)
dna=fasta("Paste DNA sequence and press enter twice")
circular=input("Is the DNA circular? Y/N\n>>> ").lower()=="y"
compute(protein,dna,circular)

while True:
    x=input("\nType 'protein', 'dna' or 'circular' to change parameters, 'go' to compute or 'quit' to quit.\n>>> ")
    if x.lower()=="quit" or x.lower()=="q":
        break
    elif x.lower()=="protein" or x.lower()=="p":
        protein=fasta("Paste protein sequence and press enter twice")+"*"
    elif x.lower()=="dna" or x.lower()=="d":
        dna=fasta("Paste DNA sequence and press enter twice")
    elif x.lower()=="circular" or x.lower()=="c":
        circular=input("Is the DNA circular? Y/N\n>>> ").lower()=="y"
    elif x.lower()=="go" or x.lower()=="g":
        compute(protein,dna,circular)
