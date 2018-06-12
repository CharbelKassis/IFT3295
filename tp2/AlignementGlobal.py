import numpy as np
from MatrixScore import getScore
from Sequences import generateSequenceArray

a = 10 # le cout d'une ouverture
b = 1 # le cout d'un gap
indel = "-" # le caractere de l'indel
moinsInfini = -999999 # la valeur qui sera toujours la plus petite

deplacements = [[-1,-1],[-1,0],[0,-1]] # les vecteures de deplacements, utiles pour le backtracking
VERTICAL = 1
HORIZONTAL = 2          

#fonction qui applique l'alignement global sur seq1 et seq2
#retourne un triplet: (Score,seq1,seq2)
def alignementGlobal(seq1,seq2):
    
    E = initMatrix(len(seq1)+1,len(seq2)+1,0).astype(int) # Matrice des indels verticals
    F = initMatrix(len(seq1)+1,len(seq2)+1,0).astype(int) # Matrice des indels horizontals
    G = initMatrix(len(seq1)+1,len(seq2)+1,0).astype(int) # Matrice des matchs et mismatchs
    V = initMatrix(len(seq1)+1,len(seq2)+1,0).astype(int) # Matrice des scores
    fleches = initMatrix(len(seq1)+1,len(seq2)+1,[0,0]).astype(int)# Matrice des fleches
    
    _casDeBase(E,F,G,V,fleches,seq1,seq2) # initialiser les matrices
    _recurrence(E,F,G,V,fleches,seq1,seq2) # remplir les tables de programmation dynamique

    # trouver la solution en utilisant les fleches, la matrice de scores et les sequences originales
    return _backtracking(V,fleches,seq1,seq2)
   

#fonction servant a generer une matrice numpy avec dimension size1 x size2 remplis de "value" 
def initMatrix(size1,size2,value):

    matrix = np.array([[value for j in range(size2)] for i in range(size1)])
    return matrix

#fonction servant a initialiser les matrices
def _casDeBase(E,F,G,V,fleches,seq1,seq2):
	
    # E(0,j) = V(0,j) = -(a+jb)
    # F(0,j) = G(0,j) = -infini
    for j in range(1,len(seq2)+1):
        e = E[0][j] = V[0][j] = -(a+j*b)
        f = F[0][j] = moinsInfini
        g = G[0][j] = moinsInfini
        V[0][j] = max(e,f,g)
        fleches[0][j] = deplacements[HORIZONTAL]

    # F(i,0) = V(i,0) = -(a+ib)
    # E(i,0) = G(i,0) = -infini
    for i in range(1,len(seq1)+1):
        f = F[i][0] = V[i][0] = -(a+i*b)
        e = E[i][0] = moinsInfini
        g = G[i][0] = moinsInfini
        V[i][0] = max(e,f,g)
        fleches[i][0] = deplacements[VERTICAL]

#fonction qui applique la recurrence : remplir la table de programmation dynamique V ainsi que les fleches
def _recurrence(E,F,G,V,fleches,seq1,seq2):
    
    for i in range(1,len(V)):
        for j in range(1,len(V[i])):
            
            # E(i,j) = max( E(i-1,j) - b , V(i-1,j) - a - b )
            e = E[i][j] = max (E[i-1][j] - b , V[i-1][j] - a - b)
        
            # F(i,j) = max( F(i,j-1) - b , V(i,j-1) - a - b )
            f = F[i][j] = max (F[i][j-1] - b , V[i][j-1] - a - b)
        
            # G(i,j) = V(i-1,j-1) + delta(v_i , w_j)
            g = G[i][j] = V[i-1][j-1] + getScore(seq1[i-1],seq2[j-1])
        
            # V(i,j) = max(G(i,j) , E(i,j) , F(i,j))
            allScores = np.array([g,e,f]).astype(int)
            V[i][j] = np.amax(allScores)
            maxIndex = np.argmax(allScores)
            fleches[i][j] = deplacements[maxIndex]
            
def _backtracking(V,fleches,seq1,seq2):
    
    i = len(V)-1
    j = len(V[0])-1
    score = V[i][j]
    newSeq1 = ""
    newSeq2 = ""
    
    while(not(i==0 and j==0)):

        #stocker en memoire les caracteres courrants et la fleche
        char1 = seq1[i-1]
        char2 = seq2[j-1]
        fleche = fleches[i][j]

        #on suppose que le nouveau caractere va etre un indel
        newChar1 = indel
        newChar2 = indel

        #si on se deplace en haut, alors on garde le caractere provenant de seq1
        if(fleche[0] == -1):
            newChar1 = char1

        #si on se deplace a gauche, alors on garde le caractere provenant de seq2
        if(fleche[1] == -1):
            newChar2 = char2

        #ajouter les caracteres a la reponse
        newSeq1 = newSeq1+newChar1
        newSeq2 = newSeq2+newChar2
        
        #on change les indexes des fleches
        i=i+fleche[0]
        j=j+fleche[1]

    #remettre les sequences dans l'ordre
    newSeq1 = newSeq1[::-1]
    newSeq2 = newSeq2[::-1]

    return (score,newSeq1,newSeq2)
