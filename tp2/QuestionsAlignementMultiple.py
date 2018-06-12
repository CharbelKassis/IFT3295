import numpy as np
from Sequences import generateSequenceArray
from AlignementGlobal import initMatrix,alignementGlobal,indel
from MatrixScore import getCharacters

#Question 1
seqs = generateSequenceArray() #le tableau contenant les sequences
scores = initMatrix(len(seqs),len(seqs),0).astype(int) #matrice des scores

for i in range(len(scores)):
    for j in range(len(scores[i])):
        if(i!=j):
            scores[i][j],seq1,seq2 = alignementGlobal(seqs[i],seqs[j])

print("Question 1:")
print()
print(scores)
print()

#Question 2

sommes = scores.sum(axis=1) #vecteur contenant la somme de chaque ligne
index = np.argmax(sommes) #chercher l'index de la sequence central

print("Question 2:")
print()
print("seq"+str(index+1)+":")
print(seqs[index])
print()

#Question 3

print("Quesiton 3:")
print()

for j in range(len(seqs)):

    if(j!=index):

        #alignement global entre s* et s_j
        score,seq1,seq2 = alignementGlobal(seqs[index],seqs[j])

        #stocker en memoire l'ancienne s*
        old_s = seqs[index]

        compare = "Comparer S* avec S"+str(j+1)

        print("*"*len(compare))
        print(compare)
        print("*"*len(compare))
        print()
        print("seq"+str(index+1)+":")
        print(seq1)
        print()
        print("seq"+str(j+1)+":")
        print(seq2)
        print()

        #mettre a jour s* et s_j
        new_s = seqs[index] = seq1
        seqs[j] = seq2

        indel_pos = [] #tableau contenant la position ou un nouveau indel ete insere dans s*

        if(len(old_s) != len(new_s)): #si il y a eu indel: taille de vieux s* differente de nouveau s*
            i2=0
            for i in range(len(new_s)):
                if(new_s[i] == old_s[i2]):
                    i2 = i2+1
                else:
                    indel_pos.append(i)

        maj = "Mis a jour des "+str(max(2,j+1))+" sequences precedentes:"
        print("*"*len(maj))
        print(maj)
        print("*"*len(maj))

        print("seq"+str(index+1)+": ")
        print(seqs[index])
        
        #mettre a jour les anciennes sequences en prenant compte les nouveaux indels de s*
        for i in range(0,j):
            if(i!=index):
                old_seq_i = seqs[i] #ancienne seq_i
                new_seq_i = ""      #nouvelle seq_i
                #pas de indel ajouter a s* donc les sequences precedentes resteront les memes
                if(len(indel_pos) == 0):
                    new_seq_i = old_seq_i
                
                #construire la nouvelle sequence en ajoutant les indels
                else:
                    start = 0
                    for k in range(len(indel_pos)):
                        end = indel_pos[k]
                        new_seq_i = new_seq_i + old_seq_i[start:end] + indel
                        start = end
                    new_seq_i = new_seq_i + old_seq_i[start:]
                
                #mettre a jour
                seqs[i] = new_seq_i
                
                print("seq"+str(i+1)+": ")
                print(seqs[i])
                
        print("seq"+str(j+1)+": ")
        print(seqs[j])
        print()
        
#Question 4

print("Question 4:")
print()

sp_array = [] #Matrice qui contiendra le nombre de missmatch de chaque colonne
nb_seqs = len(seqs) #Le nombre de sequences
seq_length = len(seqs[0]) #la taille d'une sequence

#pour chaque position de caractere j, verifier pour chaque sequence i le nombre d'erreurs qu'il y a avec les autres sequences k
for j in range(seq_length):
    missmatch = 0
    for i in range(nb_seqs):
        for k in range(i,nb_seqs):
            if(seqs[i][j] != seqs[k][j]):
                missmatch = missmatch+1
    sp_array.append(missmatch)

sp = np.array(sp_array).sum()
    
print("tableau: ")
print(sp_array)
print()
print("sp: ")
print(sp)
print()

#Question 5

print("Question 5:")
print()

characters = getCharacters()


frequences = initMatrix(len(characters)+1,seq_length+1," ").astype(object)

for i in range(1,len(frequences)):
    frequences[i][0] = characters[i-1]
for j in range(1,len(frequences[0])):
    frequences[0][j] = "c"+str(j) 

#pour chaque position j de caractere, calculer la frequence de la lettre k pour chaque sequence i
for j in range(seq_length):
    for k in range(len(characters)):
        frequence = 0.0
        for i in range(nb_seqs):
            if(characters[k] == seqs[i][j]):
                frequence = frequence+1
        frequence = frequence/nb_seqs
        frequences[k+1][j+1] = repr(frequence)
        
np.savetxt("frequences.txt",frequences,delimiter=" ",fmt='%s')
print(frequences)
print()
print("Table des frequences enregistree dans le fichier sequences.txt")
            
        
        
        
        
            
            
            
    

    
        
        

        
                
                
            
            
