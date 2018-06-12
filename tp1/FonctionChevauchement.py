import numpy as np

#Code ecrit en python36
#Par Charbel Kassis et Boumediene Boukharouba
#TP1

match = 4
missMatch = -4
indel = -8

# index=0 => provient de (miss)match => se deplacer diagonalement
# index=1 => provient de indelHaut => se deplacer vers le haut
# index=2 => provient de indelGauche => se deplacer vers la gauche
#les trois fleches possibles
fleches = [ [-1,-1] , [-1,0] , [0,-1] ] 


####################################################################################################################################
# Fonction qui retourne un quintuplet: (score, allignement(les deux sequences) la longeur du chevauchement) et si mot1 est suffixe #
####################################################################################################################################
# mot1: le premier mot, sur l'axe horizontal
# mot2: le deuxieme mot, sur l'axe vertical

def CalculeChevauchement(mot1, mot2):
	
        scores = np.zeros((len(mot1)+1, len(mot2)+1)).astype(int) # le score de tous les [i,j]
        traces = np.zeros((len(mot1)+1, len(mot2)+1,2)).astype(int) # contiendra une seule fleche a chaque [i,j]
        maxScorePosition = [0,0] # la position [x,y] du prochain score maximal
	
        for i in range(1, len(mot1) + 1):
                for j in range(1, len(mot2) + 1):
		
                        #values est un tableau qui contiendra les trois valeurs possibles que pourra avoir la case [i,j]
                        values = getValues(scores,mot1,mot2,i,j)
		
                        #le score de la case [i,j] est le maximum des trois valeurs calculer plus haut
                        scores[i][j] = np.amax(values)
			
                        #trouver l'index de la valeur maximum, l'index indiquera quelle fleche utiliser
                        maxIndex = np.argmax(values)
                        traces[i][j] = fleches[maxIndex]
			
                        #mettre a jour la position du score maximum
                        x = maxScorePosition[0]
                        y = maxScorePosition[1]
			
                        if(scores[i][j] > scores[x][y]):
                                maxScorePosition = [i,j]
				
	#trouver les sequence 1 et sequence 2
        length,sequence1,sequence2,isSuffixe,suffixeLength = getSequences(maxScorePosition,traces,mot1,mot2,scores)  
	
        #retourner le score, sequence1, sequence2, la longeur du chevauchement et si le mot1 est suffixe
        score = scores[maxScorePosition[0],maxScorePosition[1]]
        return (score,sequence1,sequence2,length,isSuffixe,suffixeLength)

#################################
# fonction qui retourne des gap #
#################################
# nbGap: le nombre de gap        
def gap(nbGap):
        return "-"*nbGap

######################################################################################################################
# fonction qui calcule les trois valeurs possibles que pourra contenir la case [i,j] et les retourne dans un tableau #
######################################################################################################################
# score: le tableau qui contient les scores
# mot1: le premier mot											
# mot2: le deuxieme mot 
# i: la position i courante
# j: la position j courante
def getValues(scores,mot1,mot2,i,j):

	matchingscore = scores[i-1][j-1] + getMatchingValue(mot1,mot2,i,j) 
	indelHaut	= scores[i-1][j] + indel
	indelGauche	= scores[i][j-1] + indel
	return np.array([matchingscore,indelHaut,indelGauche]).astype(int)
	
################################################################################################################################
# fonction qui verifie s'il y a un match entre 2 caracteres, si oui retourner la valeur de match, sinon la valeur de missmatch #
################################################################################################################################
# mot1: le premier mot
# mot2: le deuxieme mot
# i: la position i courante
# j: la position j courante
def getMatchingValue(mot1,mot2,i,j):

	if(mot1[j-1] == mot2[i-1]):
	
		return match
		
	return missMatch

############################
# envoie le chevauchement  #
############################
# maxScorePosition: la position dans le tableau "traces" correspondant au score maximum
# traces: le tableau contenant les fleches qui serviront a determine les deux sequences
# mot1: le premier mot
# mot2: le deuxieme mot
def getSequences(maxScorePosition,traces,mot1,mot2,scores):
	
        sequence1 = ""
        sequence2 = ""
	
        #le point de depart
        x = maxScorePosition[0]
        y = maxScorePosition[1]

        #taille du suffixe
        suffixeLength = 0

        #tant que la fleche n'est pas egale a [0,0], on continue de chercher les sequences, retourne les deux sequences
        while(traces[x][y][0] != 0 or traces[x][y][1] != 0):
		
                #stocker en memoire l'ancienne position. l'ancienne position indiquera la position du caractere a ecrire
                oldx = x
                oldy = y
		
                #se deplacer en utilisant la fleche
                fleche = traces[x][y]
                x = x+fleche[0]
                y = y+fleche[1]

                #par defaut, on mettra un gap
                nextChar1 = '-'
                nextChar2 = '-'
		
                #on est aller a gauche, donc la sequence2 prendra le caractere a la position oldx du mot2
                if(x == oldx-1):
                        nextChar2 = mot2[x]
                #on est aller en haut, donc la sequence1 prend le caractere a la position oldy du mot1
                if(y == oldy-1):
                        nextChar1 = mot1[y]

                sequence1 = sequence1 + nextChar1
                sequence2 = sequence2 + nextChar2
                suffixeLength = suffixeLength + 1
                
                length = len(sequence1)

                maxCol = max(scores[:,len(mot1)])
                maxLigne = max(scores[len(mot2),:])
        
                if maxCol >= maxLigne:
                        isSuffixe = True
                else:
                        isSuffixe = False

        #renverser les sequences pour retrouver l'ordre
        sequence1 = sequence1[::-1]
        sequence2 = sequence2[::-1]

        #ajouter le reste du mot et les gaps seulement si mot1 est suffixe. Si le mot1 n'est pas suffixe, alors seulement la partie
        #commune au mot1 et mot2 sera envoyee, car ce cas ne nous interesse pas.
        if(isSuffixe):
                debutmot1 = mot1[0:y]
                debut_et_sequence1 = debutmot1+sequence1

                finmot2 = mot2[maxScorePosition[0]:]
                sequence2_et_fin = sequence2+finmot2

                debut_et_sequence1 = debut_et_sequence1 + gap(len(finmot2))
                sequence2_et_fin = gap(len(debutmot1)) + sequence2_et_fin

                sequence1 = debut_et_sequence1
                sequence2 = sequence2_et_fin

        return (length,sequence1,sequence2,isSuffixe,suffixeLength)



                
                        
                        
		   

