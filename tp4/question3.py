from ete3 import Tree
import numpy as np

#generer un tableau d'arbre ETE a partir d'un tableau d'expresasion nw
def createTrees(nw_array):

    trees = []
    
    for i in range(len(nw_array)):

        trees.append(Tree(nw_array[i]))

    return trees


#retourne un ensemble de toutes les feuilles decendantes d'un arbre ETE
def getDecendants(arbre):

    decendants = []

    for node in arbre.iter_descendants("postorder"):

        if(node.is_leaf()): #pour eviter de retourner un noeud interne, qui a un nom vide

            decendants.append(node.name)

    return set(decendants)

#retourne toute les bipartions non triviales d'un arbre non enracine. On suppose que le noeud parent a toujours 2 enfants
def getBipartitions(arbre):

    #une bipartion = [decendant de i , (decendant de l'arbre - decendant de i)]
    allDecendants = getDecendants(arbre)
    partitions = [] #garder en memoire les partitions qui ont deja etaient utilisee
    bipartitions = [] #toutes les bipartitions

    for node in arbre.iter_descendants("postorder"):

        if(node.is_leaf()):#ne rien faire si c'est une feuille

            continue
  
        nodeDecendants = getDecendants(node) #toutes les feuilles decendantes du noeud
        reste = allDecendants-nodeDecendants #toutes les feuilles de l'arbre sauf les decendants du noeud
        
        if(not len(reste) == 1): #s'assurer que l'autre partition n'est pas triviale

            if(not(nodeDecendants in partitions)): #s'assurer qu'on a pas deja utilisee la meme partition

                bipartitions.append([nodeDecendants,reste])
                partitions.append(nodeDecendants)
                partitions.append(reste)

    return bipartitions

#retourne le score RF entre bipartition1 et bipartition2
def getScoreRF(bp1,bp2):

    score = len(bp1) + len(bp2)
    
    for partition1 in bp1:

        for partition2 in bp2:

            if(partition1 == partition2):

                score = score-2

    return score
                       
if(__name__ == "__main__"):

    nw = ["((PCDHA1_Humain, OR2J3_Humain), ((PCDHA1_Rat, PCDHA1_Souris), PCDHA1_Bonobo));","(((PCDHA1_Humain, PCDHA1_Rat), (PCDHA1_Souris, PCDHA1_Bonobo)), OR2J3_Humain);","(((OR2J3_Humain, PCDHA1_Rat), (PCDHA1_Souris, PCDHA1_Bonobo)), PCDHA1_Humain);","((((PCDHA1_Rat,PCDHA1_Souris),OR2J3_Humain),PCDHA1_Humain),PCDHA1_Bonobo);","((((PCDHA1_Humain,PCDHA1_Souris),PCDHA1_Rat),OR2J3_Humain),PCDHA1_Bonobo);"]
    
    arbres = createTrees(nw) #creer des arbre pour chaque nw
    bp = [] #tableau des bipartition
    
    #creer les bipartitions pour chaque arbre
    for arbre in arbres:
        
        bp.append(getBipartitions(arbre)) 

    score_rf = [[0 for j in range(len(arbre))] for i in range(len(arbre))]

    #remplir la matrice score RF de tous les arbres
    for i in range(len(score_rf)):
        
        for j in range(i,len(score_rf[i])):

            if(i != j):
                
                score_rf[i][j] = getScoreRF(bp[i],bp[j])

        print(score_rf[i])

    




    














