#Generer un tableau contenant toutes les sequences du fichier sequences.fasta
def generateSequenceArray():

    filename = "sequences.fasta"
    seqFile = open(filename,"r")
    sequences = [] #la liste de toutes les sequences
    sequence = "" #la sequence courante
    
    for line in seqFile:
        
        if(line[0] == ">"): #cette ligne indique le numero de la sequence
            continue
        elif(line.isspace()): #la sequence courante a ete completement extraite du fichier, l'ajouter au tableau des sequences
            sequences.append(sequence)
            sequence = ""
        else:
            if(line[-1] == "\n"):
                line = line[:-1] #enlever le caracter \n 
            sequence = sequence + line
            
    sequences.append(sequence) #ajouter la derniere sequence
            
    seqFile.close()     
    return sequences
    
    
         
               
