#Code ecrit en python36
#Par Charbel Kassis et Boumediene Boukharouba
#TP1

def complementReverse(read):
        
        complements = {"A":"T", "C":"G", "T":"A", "G":"C"}
        complement = ""
        
        for i in range(len(read)):

            complement = complement+complements[read[i]]

        return complement[::-1]

def reductionTransitive(scoreMatrix):

        length = len(scoreMatrix)
        
        for i in range(length):
                for j in range(length):
                        if(scoreMatrix[j][i] > 0):
                                for k in range(length):
                                        if(scoreMatrix[i][k] > 0):
                                                scoreMatrix[j][k] = 0

########################################################################
# fonction servant a produire une entree dans un fichier de type FASTQ #
########################################################################
# output_file: le fichier .fq dans lequel on va ajoute une nouvelle entree
# index: le numero de la sequence
# sequence: la sequence
def writeFASTQ(output_file,index,sequence):
        output_file.write("@READS_"+str(index)+'\n')
        output_file.write(sequence)
        output_file.write("+\n")
        output_file.write("2"*len(sequence)+'\n')


#####################################################################################
# fonction servant a ecrire dans le fichier output_file: seq(i) seq(j) scores[i][j] #
#####################################################################################
# i: sequence i
# j: sequence j
# readMatrix: la matrice des scores entre le read i et le read j
# reads: la listes des reads ainsi que leur index
# output_file: le fichier de sortie de format EdgeList
def writeEdge(i,j,readMatrix,reads,output_file):

        #Si le tableau des reads n'est pas donne, alors on suppose que les reads sont dans l'ordre
        if(reads == None):
                i_value = i+1
                j_value = j+1
        else:
                i_value = reads[i][0]
                j_value = reads[j][0]

        output = "seq"+str(i_value)+" "+ "seq"+str(j_value) + " "+str(readMatrix[i][j])
        print("\n\tWriting "+output)
        output_file.write(output+"\n")
        output_file.flush()
                                                
                                                





        
    
