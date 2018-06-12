from fasta import Fasta
from ete3 import Tree
from blosum  import getScore
from question3 import createTrees,createBP,rfMatrix
import numpy as np

nw_array = ["((PCDHA1_Humain, OR2J3_Humain), ((PCDHA1_Rat, PCDHA1_Souris), PCDHA1_Bonobo));","(((PCDHA1_Humain, PCDHA1_Rat), (PCDHA1_Souris, PCDHA1_Bonobo)), OR2J3_Humain);","(((OR2J3_Humain, PCDHA1_Rat), (PCDHA1_Souris, PCDHA1_Bonobo)), PCDHA1_Humain);","((((PCDHA1_Rat,PCDHA1_Souris),OR2J3_Humain),PCDHA1_Humain),PCDHA1_Bonobo);","((((PCDHA1_Humain,PCDHA1_Souris),PCDHA1_Rat),OR2J3_Humain),PCDHA1_Bonobo);"]


#question 2.1
def poids(seq1,seq2):

    score = 0
    
    for i in range(len(seq1)):

        score += getScore(seq1[i],seq2[i])

    return score

#initialisation de la matrice de distance
def initMatriceDist(seq):

    d = [[0 for j in range(len(seq))] for i in range(len(seq))]

    for i in range(len(d)):

        for j in range(len(d[i])):

            p = poids(seq[i],seq[j])
            qi = poids(seq[i],seq[i])
            qj = poids(seq[j],seq[j])

            d[i][j] = 1-(p/max(qi,qj))

    return d

#calcul la matrice des poids Q
def QMatrix(d):

    q = [[0 for j in range(len(d))] for i in range(len(d))]

    for i in range(len(q)):

        for j in range(len(q[i])):

            if(i != j):

                ri = somme_r(i,j,d)
                rj =somme_r(j,i,d)
                q[i][j] = d[i][j] - (ri+rj)

    return q

#calcul la somme r pour la fonction QMatrix
def somme_r(i,j,d):

    diviseur = len(d)-2
    somme = 0
            
    for k in range(len(d)):
        
         if(k != i and k != j):

             somme += d[i][k]

    return somme/diviseur

#fusioner deux colonnes dans la matrice de distance, la colonne fusionee sera a la position j=0
def fusion(dist,a,b):
    
    new_dist = [[0 for j in range(len(dist)-1)] for i in range(len(dist)-1)]
    u = [0]
        
    #generer le vecteur u
    for c in range(len(dist)):
    
        if(c!=a and c!=b):

            u.append((dist[a][c] + dist[b][c] - dist[a][b] ) / 2)
    
    #mettre u dans la premiere ligne
    new_dist[0] = u
    
    #mettre u dans la premiere colonne
    for i in range(len(new_dist)):
        
        new_dist[i][0] = u[i]
        
    i2=1
    j2=1
    
    for i in range(len(dist)):

        if(i == a or i == b):

            continue
        
        for j in range(len(dist[i])):

            if(j == a or j == b):

                continue
           
            new_dist[i2][j2] = dist[i][j]
            j2 += 1
            
        i2+=1
        j2=1

    return new_dist
    
#retourne l'index de la case avec la plus petite valeur dans une matrice
def getMinIndex(arr):
    
    a=0;b=0
    min = arr[0][0]
        
    for i in range(len(arr)):
        
        for j in range(len(arr[i])):
            
            if(arr[i][j] < min):

                min = arr[i][j]
                a = i
                b = j
                
    return (a,b)


#combiner les noeuds ou feuilles a,b pour en former qu'un seul
def combineSeq(tab,a,b):

    newTab = []

    newTab.append("("+tab[a]+", "+tab[b]+")") #le nouveau noeud sera de forme (a, b)

    for i in range(len(tab)): #ajouter le reste

        if(i != a and i != b):

            newTab.append(tab[i])

    return newTab


if(__name__ == "__main__"):
    
    fasta = Fasta("proteines.fa")
                   
    #les indel on ete enleves
    fasta.filtrer()

    #Quesiton 4
    print("Question 4")
    print("afficher les sequences filtrees")
    print()

    seqs = [] #ce tableau sera la liste des sequences
    seqTree = [] #ce tableau sera utilise pour generer une expression newick
    
    for nom_seq in fasta:

        seqs.append(fasta[nom_seq]) #ajoute la sequence dans le tableau seq
        seqTree.append(nom_seq) #ajoute le nom de la sequence dans le tableau de l'arbre des sequences
        print(nom_seq) #le nom de la sequence
        print(fasta[nom_seq]) #la sequence

    print()
    print("_"*50)

    #Question 2.1
    print()
    print("Question 2.1")
    print("Afficher la matrice de distance")
    dist = initMatriceDist(seqs)

    print("Matrice distance initiale")
    print()
    for i in range(len(dist)):
        print(dist[i])

    print()
    print("_"*50)
    print()
    print("Question 2.2")

    print(seqTree)
    
    #imprimer chaque nouvelle matrice Q et chaque nouvelle matrice de score
    i=0
    while(len(dist) > 2): #tant qu'il reste plus que 2 feuilles

        print()
        print("Matrice q" + str(i))
        q = QMatrix(dist)
        
        for j in range(len(q)): #imprimer la nouvelle matrice Q
            print(q[j])
            
        print()
        print("Matrice de score " + str(i))
        a,b = getMinIndex(q)
        seqTree = combineSeq(seqTree,a,b) #mettre a jour le tableau qui contiendra les expressions newick
        dist = fusion(dist,a,b) #fusioner les valeurs des distances

        for j in range(len(dist)): #imprimer la nouvelle matrice de distances
            print(dist[j])
        
        i+=1

    '''fusioner les enfants de la racine, puis aller extraire le seul element restant du tableau
    lui ajouter ";" a la fin pour obtenir une expression newick valide'''
    expr_nw = combineSeq(seqTree,0,1)[0]+";"

    print()
    print("Voici l'expression newick obtenue: ")
    print(expr_nw)
    
    tree = Tree(expr_nw)
    print(tree)
    print()
    print("_"*50)

    print()
    print("Question 2.3")
    print("Afficher la nouvelle matrice des score rf")
    
    nw_array.append(expr_nw) #ajouter cet arbre a ceux utilises dans les questions precedentes
    arbres = createTrees(nw_array) #cree les arbres
    bp_tab = createBP(arbres) #chercher leur bipartion
    rf_score = rfMatrix(bp_tab) #generer la matrice des scores rf

    for i in range(len(rf_score)):
        print(rf_score[i])


    
    



        


        
    

        
        
    

    
    
    

        
    

    
    



        
