from FonctionChevauchement import CalculeChevauchement
from methodes_utiles import writeEdge
import numpy as np
import time

#Code ecrit en python36
#Par Charbel Kassis et Boumediene Boukharouba
#TP1

########################################################################################################################################################
# retourne des matrices de taille (nb reads x nb reads) avec le score de chevauchement du mot i et j, et leurs sequences et genere un fichier EdgeList #
########################################################################################################################################################
# input_file: le fichier de format FASTQ
# output_file: le fichier de sortie de format EdgeList

def calculTabReads(input_file,output_file):

        input_array = input_file.readlines()
        
        #la matrice contenant tous les reads (a l'index 1) ainsi que leur index (a l'index 0)
        reads = []
        #la matrice qui servira a stocker le score entre le read i et le read j
        scoreMatrix = []
        #la matrice qui servira a stocker le sequencement entre le read i et le read j ainsi que la taille du suffixe (prefixe)
        #indice 0: sequence i
        #indice 1: sequence j
        #indice 2: taille du suffixe
        seqMatrix = []
        #la matrice qui indiquera si pour la paire [i,j], la sequence i est un suffixe du mot i, utile pour eviter de calculer [j,i]
        suffixeMatrix = []
        
        for i in range(0,len(input_array),4):
                #chercher l'index du read qui se trouve tout de suite apres "@READS_"
                read_index = int(float(input_array[i][7:-1]))

                #ajouter l'index ainsi que la sequence dans le tableau reads
                reads.append([read_index,input_array[i+1][:-1]])

        length = len(reads)

        for i in range(length):
                seqRow = []
                for j in range(length):
                        seqRow.append([""]*2+[0])
                seqMatrix.append(seqRow)
                scoreMatrix.append([0] * length)
                suffixeMatrix.append([False] * length)

        print("Generating EdgeList file "+output_file.name)
        for i in (range(length)):
                for j in (range(length)):

		        #ne pas comparer un read avec lui-meme
                        if( i != j ):

                                #calculer seulement si [j][i] n'est pas un suffixe
                                if(not suffixeMatrix[j][i]):
                                        timerStart = time.time()
                                        score,seq1,seq2,size,isSuffixe,suffixeLength = CalculeChevauchement(reads[i][1],reads[j][1])
                                        timerEnd = time.time()
                                        suffixeMatrix[i][j] = isSuffixe
                                        print("\nComparing "+str(reads[i][0])+" to "+str(reads[j][0]))
                                        print("\t(time elapsed: " + "%.2f"%(timerEnd-timerStart) + " seconds)")
                                        print("\tscore: "+str(score))
				
                                        if(score <= 80):
                                                scoreMatrix[i][j] = 0
                                        else:
                                                if(isSuffixe):
                                                        scoreMatrix[i][j] = score
                                                        seqMatrix[i][j][0] = seq1
                                                        seqMatrix[i][j][1] = seq2
                                                        seqMatrix[i][j][2] = suffixeLength
                                                        writeEdge(i,j,scoreMatrix,reads,output_file)
                                                else:
                                                        scoreMatrix[j][i] = score
                                                        seqMatrix[j][i][0] = seq1
                                                        seqMatrix[j][i][1] = seq2
                                                        seqMatrix[j][i][2] = suffixeLength
                                                        
                                                        writeEdge(j,i,scoreMatrix,reads,output_file)
                                                        suffixeMatrix[j][i] = True

                                else:
                                        print("\nSkipping seq"+str(reads[i][0])+" "+ "seq"+str(reads[j][0]))
        print("\nFinished generating EdgeList file "+output_file.name)

        return scoreMatrix,seqMatrix
        

                                
			



