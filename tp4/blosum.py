indel = "-"
scoreMatrix = False
characterIndex = {}
characterList = []
     
#generer la matrice des score a partir du fichier
def _generateScoreMatrix():

    filename = "BLOSUM62.txt"
    blosum62 = open(filename,"r")
    lines = blosum62.readlines()
    global scoreMatrix
    scoreMatrix = []
 
    #remplir la matrice score
    for i in range (1,len(lines)):
        values = lines[i].split()
        del values[0] #enlever la lettre
 
        for j in range (0,len(values)):
            values[j] = int(values[j])
         
        scoreMatrix.append(values)
 
    #remplir le dictionaire des indexes de caracteres. Utile pour trouver le score dans la matrice des scores
    global characters
    characters = lines[0]
    index = 0
    for char in characters:
        global characterIndex
        global characterList
        if char.isalpha():
            characterIndex[char] = index
            characterList.append(char)
            index = index+1
        elif char == "*":
            characterIndex[indel] = index
            characterList.append(indel)
            index = index+1
 
    blosum62.close()
 
#retourne le score entre le character a et b
def getScore(a,b):
     
    #generer la matrice si ce n'est pas deja fait
    global scoreMatrix
    if not scoreMatrix:
        _generateScoreMatrix()
         
    #regarder la matrice des scores apres avoir extrait l'index des deux caracteres
    index_a = characterIndex[a]
    index_b = characterIndex[b]
 
    return scoreMatrix[index_a][index_b]
 
#retourne les caracteres
def getCharacters():
 
    if characterList:
 
        return characterList
 
    else:
 
        return None
